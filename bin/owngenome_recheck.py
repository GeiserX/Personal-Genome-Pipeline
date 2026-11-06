#!/usr/bin/env python3
"""owngenome_recheck.py: say which of a finished sample's results a newer ClinVar or catalog would change.

  owngenome_recheck.py --sample S [--genome-dir DIR] [--offline FILE] [--owngenome-repo OWNER/NAME]
                 [--rerun [--refresh-clinvar]]
  owngenome_recheck.py --affected KEY[,KEY...]      the steps a change of these manifest keys redoes
  owngenome_recheck.py --launch-dir DIR             the last bare `nextflow run` there, with -resume

It reads the sample's run_manifest.tsv (bin/write_manifest.sh): the databases
and pinned data the run used. For each one it asks the upstream source for its
newest release, through the lookups of scripts/ci/freshness.py (the monthly CI
check), plus one PGS Catalog lookup per score file. It also reads the newest
release of this pipeline and the versions.env of that tag, to say whether a
newer pin is already in a release. It sends release names only, never anything
about the sample, and prints every host it asks. --offline FILE reads the
answers from a JSON file instead (an air-gapped server, the tests).

Each item gets one class:
  rerun-ready          what is installed here is newer than what the run used:
                       a rerun picks it up
  refresh-available    ClinVar only: NCBI has a newer file than the one installed
  in-pipeline-release  a newer pin is in a pipeline release: update the pipeline
  upstream-only        upstream has a newer release that no pipeline release pins
  no-newer-release     the lookup answered and nothing is newer
  unknown              a pin with no upstream lookup that the run's manifest
                       does not name: nothing was compared
  error                the lookup failed (never read as up to date); exit 1
The steps an item feeds come from the reads column of assets/steps.tsv, then
every step downstream of them through its upstream column, then the reports.

The result goes to <sample>/recheck/state.json (rewritten whole) and one row per
item to <sample>/recheck/log.tsv (only ever appended). Nothing schedules itself:
docs/recheck.md has a crontab line, a systemd timer and a launchd job to copy.

--rerun replays the last run-all.sh command, which run-all.sh writes into
logs/run_status.tsv ('meta command'), with the same switches. Nextflow's
-resume then redoes only the tasks whose inputs changed. It refuses while a run
is going on. Afterwards it compares the steps the run redid with the steps this
check named, and lists the named script steps (04b, 14, 29, 03a-03d), which
skip when their output exists, as not redone. --refresh-clinvar runs
setup.sh --refresh clinvar first, and refuses while a run of any sample of
GENOME_DIR is going on, since they share the ClinVar files. It holds
GENOME_DIR/clinvar/.refresh.lock while it refreshes: run-all.sh and step 06
wait on it. Standard library only.
"""
import argparse
import fcntl
import gzip
import importlib.util
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import urllib.parse
from datetime import datetime, timedelta, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STEPS_TSV = os.path.join(ROOT, "assets", "steps.tsv")
sys.path.insert(0, HERE)
import run_state  # noqa: E402  the steps catalogue and the run's state

SCHEMA_VERSION = 1
DEFAULT_REPO = "GeiserX/OwnGenome"
PGS_API = "https://www.pgscatalog.org/rest/score"
RAW = "https://raw.githubusercontent.com"
REPORT_STEP = "24"  # the HTML report; the text report is written after it
# ClinVar's ##fileDate is the release; NCBI's Last-Modified comes a day or two later.
CLINVAR_SLACK_DAYS = 2
CLASSES = ("rerun-ready", "refresh-available", "in-pipeline-release", "upstream-only", "no-newer-release", "unknown",
           "error")
LOG_HEADER = "utc\tevent\titem\tclass\tdetail\n"
# Held exclusive by --refresh-clinvar, shared by scripts/run-all.sh and scripts/06-clinvar-screen.sh.
CLINVAR_LOCK = "clinvar/.refresh.lock"
CLINVAR_FILES = ("clinvar/clinvar_pathogenic_chr.vcf.gz", "clinvar/clinvar_chr.vcf.gz", "clinvar/clinvar.vcf.gz")
# The pinned items: (manifest key, label, freshness lookup, how to read its upstream).
PINNED = (
    ("VEP_CACHE_RELEASE", "Ensembl / VEP cache", "ensembl", r"Ensembl release (\d+)"),
    ("PCGR_DATA_BUNDLE", "PCGR / CPSR data bundle", "pcgr", r"bundle (\d+)"),
    ("PYPGX_BUNDLE_VERSION", "pypgx and its bundle", "pypgx", None),
    ("PHARMCAT_IMAGE", "PharmCAT", "pharmcat", r"^(\d+(?:\.\d+)+)"),
    ("hla_database", "IPD-IMGT/HLA", "imgt", r"^(\d+(?:\.\d+)+)"),
)


# The opt-in steps' pins: the setup.sh option that installs their data and image (the main setup leaves
# them out). VNtyper's reference data is inside its image, which the run pulls, so it needs no option.
OPT_IN_SETUP = {"MEGANE_IMAGE": "--megane-kmers", "MOCHA_PANEL_VERSION": "--mocha-data", "MOCHA_IMAGE": "--mocha-data",
                "LOCITYPER_PANEL_VERSION": "--locityper-data", "LOCITYPER_IMAGE": "--locityper-data",
                "VNTYPER_IMAGE": None}


class Failed(Exception):
    """A lookup with no usable answer."""


def now_utc():
    return datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_freshness():
    path = os.path.join(ROOT, "scripts", "ci", "freshness.py")
    spec = importlib.util.spec_from_file_location("freshness", path)
    mod = importlib.util.module_from_spec(spec)
    keep, sys.dont_write_bytecode = sys.dont_write_bytecode, True
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.dont_write_bytecode = keep
    if os.environ.get("OWNGENOME_RECHECK_RETRIES"):
        mod.RETRIES = int(os.environ["OWNGENOME_RECHECK_RETRIES"])
    return mod


def ver(text):
    return tuple(int(x) for x in re.findall(r"\d+", str(text or "")))


def date8(text):
    """The first date in text as YYYYMMDD ('20261005', '2026-10-05'), or None."""
    m = re.search(r"(\d{4})-?(\d{2})-?(\d{2})", str(text or ""))
    if not m:
        return None
    try:
        return datetime.strptime("".join(m.groups()), "%Y%m%d").strftime("%Y%m%d")
    except ValueError:  # an impossible date (20260230) reads as no date, so the item is an error, not a crash
        return None


class Remembered:
    """The answers of one check, kept, so the check again after --rerun asks nobody a second time."""

    def __init__(self, inner):
        self.inner, self.hosts, self._seen = inner, inner.hosts, {}

    def _ask(self, name, *args):
        k = (name,) + args
        if k not in self._seen:
            try:
                self._seen[k] = (True, getattr(self.inner, name)(*args))
            except Failed as e:
                self._seen[k] = (False, e)
        ok, v = self._seen[k]
        if not ok:
            raise v
        return v

    def clinvar(self):
        return self._ask("clinvar")

    def pinned(self, key, lookup, pattern):
        return self._ask("pinned", key, lookup, pattern)

    def pgs(self, pgs_id):
        return self._ask("pgs", pgs_id)

    def release(self):
        return self._ask("release")


# --- the files ------------------------------------------------------------------------------

def read_manifest(path):
    rows = {}
    with open(path, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and not c[0].startswith("#"):
                rows[(c[0], c[1])] = "\t".join(c[2:])
    return rows


def header_value(path, key):
    """The value of '#key=' or '##key=' in the leading comment lines of a (gzipped) text file."""
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", errors="replace") as f:
            for n, line in enumerate(f):
                if not line.startswith("#") or n > 500:
                    break
                m = re.match(r"^#+" + re.escape(key) + r"=(.*)$", line.rstrip("\n"))
                if m:
                    return m.group(1).strip()
    except (OSError, EOFError, ValueError):
        return None
    return None


def installed_clinvar(genome_dir):
    for f in CLINVAR_FILES:
        p = os.path.join(genome_dir, f)
        if os.path.isfile(p):
            d = header_value(p, "fileDate")
            if d:
                return date8(d)
    return None


def pgs_field(value, field):
    m = re.search(r"(?:^|;)" + re.escape(field) + r"=([^;]*)", value or "")
    return m.group(1) if m else None


def read_status_meta(sample_dir):
    st = run_state.read_run_status(os.path.join(sample_dir, "logs", "run_status.tsv"))
    return (st or {}).get("meta", {})


# --- what the steps read -------------------------------------------------------------------

def reads_of(step):
    return [x for x in (step.get("reads") or "").split(",") if x]


def key_matches(key, read):
    """read is a reads entry: KEY, a prefix PREFIX* (pgs:*), either with @self."""
    read = read[:-5] if read.endswith("@self") else read
    return key == read or (read.endswith("*") and key.startswith(read[:-1]))


def affected_steps(steps, keys):
    """The steps that read one of keys, every step downstream of them, then the HTML report. A step
    that reads a key only for a table of its own (KEY@self) is marked, and nothing below it."""
    direct = {s["step"]: all(r.endswith("@self") for r in reads_of(s) if any(key_matches(k, r) for k in keys))
              for s in steps if any(key_matches(k, r) for k in keys for r in reads_of(s))}
    if not direct:
        return []
    hit = {n for n, own in direct.items() if not own}
    grew = True
    while grew:
        grew = False
        for s in steps:
            if s["step"] not in hit and set(s["upstream"]) & hit:
                hit.add(s["step"])
                grew = True
    hit |= set(direct) | {REPORT_STEP}
    return [s["step"] for s in steps if s["step"] in hit]


# --- the upstream answers ------------------------------------------------------------------

class Online:
    """Answers through freshness.py, printing each host the first time it is asked."""

    def __init__(self, fr, pins, owngenome_repo, pgs_api=PGS_API, raw=RAW, out=None):
        self.fr, self.pins, self.repo, self.pgs_api, self.raw = fr, pins, owngenome_repo, pgs_api, raw
        self.hosts, self.out = [], out or sys.stderr
        orig = fr.request

        def request(url, *a, **k):
            host = urllib.parse.urlsplit(url).netloc
            if host not in self.hosts:
                self.hosts.append(host)
                print(f"recheck: asking {host} (release names only, nothing about the sample)", file=self.out)
            return orig(url, *a, **k)
        fr.request = request

    def _row(self, name, *args):
        rep = self.fr.Report()
        try:
            getattr(self.fr, name)(rep, *args)
        except self.fr.LOOKUP_ERRORS as e:
            raise Failed(str(e))
        if not rep.db_rows:
            raise Failed(f"freshness.{name} wrote no row")
        return rep.db_rows[-1]

    def clinvar(self):
        row = self._row("clinvar", os.path.join(ROOT, "scripts", "setup.sh"))
        m = re.search(r"dated (\d{4}-\d{2}-\d{2})", row[2])
        if not m:
            raise Failed(f"no date in the ClinVar answer: {row[2]}")
        return date8(m.group(1))

    def pinned(self, key, lookup, pattern):
        args = (self.pins,)
        row = self._row(lookup, *args)
        if lookup == "pypgx":  # the newest pypgx with an image and a bundle tag of its version
            m = re.match(r"can move to (\S+)", row[3])
            return m.group(1) if m else self.pins.get("PYPGX_BUNDLE_VERSION")
        m = re.search(pattern, row[2])
        if not m:
            raise Failed(f"cannot read a release from the {lookup} answer: {row[2]}")
        return m.group(1)

    def pgs(self, pgs_id):
        try:
            body = self.fr.request(f"{self.pgs_api}/{pgs_id}")[0]
            if not isinstance(body, dict) or body.get("id") != pgs_id:
                raise Failed(f"the PGS Catalog has no score {pgs_id}")
            url = body["ftp_harmonized_scoring_files"]["GRCh38"]["positions"]
            _, h = self.fr.request(url, method="HEAD", want="head")
            lm = self.fr.need(h.get("Last-Modified"), f"Last-Modified of {url}")
        except self.fr.LOOKUP_ERRORS as e:
            raise Failed(str(e))
        import email.utils
        try:  # a malformed date is a failed lookup of this score, never a crash of the whole check
            return email.utils.parsedate_to_datetime(lm).strftime("%Y%m%d")
        except (TypeError, ValueError, AttributeError, OverflowError):
            raise Failed(f"cannot read the Last-Modified of {url}: {lm!r}")

    def release(self):
        try:
            tag = self.fr.need(self.fr.gh_api(f"/repos/{self.repo}/releases/latest").get("tag_name"),
                               f"latest release of {self.repo}")
            text = self.fr.request(f"{self.raw}/{self.repo}/{tag}/versions.env", want="text")[0]
            common = self.fr.request(f"{self.raw}/{self.repo}/{tag}/scripts/lib/common.sh", want="text")[0]
        except self.fr.LOOKUP_ERRORS as e:
            raise Failed(str(e))
        pins = dict(re.findall(r'^([A-Z][A-Z0-9_]*)="([^"]*)"', text, re.M))
        m = re.search(r"^HLA_DB_RELEASE=\$\{HLA_DB_RELEASE:-(\d+(?:\.\d+)+)\}", common, re.M)
        if m:
            pins["hla_database"] = m.group(1)
        return tag, pins


class Offline:
    """Answers from a JSON file: {item: {"upstream": value} or {"error": text}, "release": {"tag", "pins"}}."""

    def __init__(self, path):
        with open(path) as f:
            self.data = json.load(f)
        self.hosts = []

    def _get(self, key):
        e = self.data.get(key)
        if not isinstance(e, dict):
            raise Failed(f"{key} is not in the offline answers")
        if e.get("error"):
            raise Failed(str(e["error"]))
        return e

    def clinvar(self):
        return date8(self._get("clinvar_file_date").get("upstream"))

    def pinned(self, key, lookup, pattern):
        return self._get(key).get("upstream")

    def pgs(self, pgs_id):
        return date8(self._get(f"pgs:{pgs_id}").get("upstream"))

    def release(self):
        e = self._get("release")
        tag, pins = e.get("tag"), e.get("pins") or {}
        if not isinstance(tag, str) or not tag or not isinstance(pins, dict):
            raise Failed("the offline answers' release needs a tag and a pins object")
        return tag, dict(pins)


# --- the check -----------------------------------------------------------------------------

def image_tag(ref):
    ref = (ref or "").split("@", 1)[0]
    return ref.rsplit(":", 1)[1] if ":" in ref.rsplit("/", 1)[-1] else ref


def pin_value(key, value):
    """A manifest or versions.env value as the release it names."""
    if value is None:
        return None
    if key.endswith("_IMAGE"):
        return image_tag(value)
    if key == "hla_database":
        m = re.search(r"(\d+\.\d+\.\d+)", value)
        return m.group(1) if m else None
    return value


def classify(item, release_tag, release_pins, release_error, clinvar=False):
    if item.get("error"):
        return "error"
    ran, inst, up = item.get("ran_with"), item.get("installed"), item.get("upstream")
    if inst and ran and ver(inst) > ver(ran):
        return "rerun-ready"
    # An image re-pinned at the same tag (a new digest): the tags compare equal, the references do not.
    ran_ref, inst_ref = item.get("ran_with_ref"), item.get("installed_ref")
    if inst and ran and ver(inst) == ver(ran) and ran_ref and inst_ref and ran_ref != inst_ref:
        return "rerun-ready"
    # ClinVar's date is in the manifest whenever step 06 had a file, so a run without one cannot be compared
    # with the file here. A pin with no run value (HLA data not installed then) compares this checkout's.
    if clinvar and not ran:
        item.setdefault("error", "the manifest does not say which ClinVar file the run used")
        return "error"
    base = inst or ran
    if not base:
        item.setdefault("error", "the manifest does not say which release the run used")
        return "error"
    if not ver(up):
        item.setdefault("error", "the answer names no release")
        return "error"
    if clinvar:
        newer = datetime.strptime(up, "%Y%m%d") - datetime.strptime(base, "%Y%m%d") > timedelta(days=CLINVAR_SLACK_DAYS)
        return "refresh-available" if newer else "no-newer-release"
    if ver(up) <= ver(base):
        return "no-newer-release"
    if item["kind"] == "file":
        return "upstream-only"
    if release_error:
        item["error"] = f"the pipeline's newest release could not be read: {release_error}"
        return "error"
    rel = pin_value(item["item"], release_pins.get(item["item"]))
    if rel and ver(rel) > ver(base):
        item["in_owngenome_release"] = {"tag": release_tag, "pin": rel}
        return "in-pipeline-release"
    return "upstream-only"


def steps_not_run(sample, sample_dir):
    """The steps this run did not run: the ones run-all.sh recorded as skipped (no BAM, missing
    data), and the ones it wrote no row for, which is how it leaves out an opt-in step not in TOOLS
    (run_state's "not planned"). Empty when logs/run_status.tsv records no step at all (a run of an
    older version), so nothing is left out on a guess."""
    try:
        status = run_state.read_run_status(os.path.join(sample_dir, "logs", "run_status.tsv"))
        if not (status or {}).get("steps"):
            return set()
        return {s["step"] for s in run_state.collect(sample, sample_dir)["steps"]
                if s["state"] in ("skipped", "not planned")}
    except (OSError, ValueError, KeyError):
        return set()


def action(item, sample):
    c, key = item["class"], item["item"]
    if c == "rerun-ready":
        if item["kind"] == "pin" and key in OPT_IN_SETUP:
            opt = OPT_IN_SETUP[key]
            how = (f"run ./scripts/setup.sh {opt} $GENOME_DIR (the main setup leaves this step's data out)"
                   if opt else "the run pulls the image")
            return f"this checkout pins {item['installed']}: {how}, then rerun: bin/owngenome_recheck.py --sample {sample} --rerun"
        if item["kind"] == "pin":
            return (f"this checkout pins {item['installed']}: run ./scripts/setup.sh (it installs the data), "
                    f"then rerun: bin/owngenome_recheck.py --sample {sample} --rerun")
        return f"rerun: bin/owngenome_recheck.py --sample {sample} --rerun"
    if c == "refresh-available":
        return (f"refresh it with ./scripts/setup.sh --refresh clinvar $GENOME_DIR, then rerun "
                f"(bin/owngenome_recheck.py --sample {sample} --rerun --refresh-clinvar does both)")
    if c == "in-pipeline-release":
        r = item["in_owngenome_release"]
        return f"update the pipeline to {r['tag']} (it pins {r['pin']}), run setup.sh, then rerun"
    if c == "upstream-only":
        if key.startswith("pgs:"):
            return f"download {key[4:]} again from the PGS Catalog (docs/25-prs.md), then rerun"
        return "no pipeline release pins it yet; nothing to do until one does"
    if c == "error":
        return "the lookup failed; check again later"
    if c == "unknown":
        return ("the run's manifest does not name this pin (a run of an older version): nothing to compare; "
                "the next run records it")
    return ""


def check(sample, genome_dir, answers, steps, this_pins, this_hla):
    sample_dir = os.path.join(genome_dir, sample)
    man = read_manifest(os.path.join(sample_dir, "run_manifest.tsv"))
    version = {k: v for (sec, k), v in man.items() if sec == "version"}
    items = []
    try:
        tag, rel_pins = answers.release()
        rel_err = None
    except Failed as e:
        tag, rel_pins, rel_err = None, {}, str(e)

    def add(item, fn, *args, clinvar=False):
        try:
            item["upstream"] = fn(*args)
        except Failed as e:
            item["error"] = str(e)
        item["class"] = classify(item, tag, rel_pins, rel_err, clinvar=clinvar)
        item["steps"] = affected_steps(steps, [item["item"]])
        item["action"] = action(item, sample)
        items.append(item)

    add({"item": "clinvar_file_date", "label": "ClinVar (step 06's file)", "kind": "file",
         "ran_with": date8(man.get(("data", "clinvar_file_date"))), "installed": installed_clinvar(genome_dir)},
        answers.clinvar, clinvar=True)
    for key, label, lookup, pattern in PINNED:
        ran = pin_value(key, man.get(("data", key)) if key == "hla_database" else version.get(key))
        inst = this_hla if key == "hla_database" else pin_value(key, this_pins.get(key))
        item = {"item": key, "label": label, "kind": "pin", "ran_with": ran, "installed": inst}
        if key.endswith("_IMAGE"):
            item["ran_with_ref"], item["installed_ref"] = version.get(key) or None, this_pins.get(key) or None
        if key == "hla_database" and not this_hla:
            item["error"] = "could not read this checkout's HLA_DB_RELEASE (scripts/lib/common.sh)"
        add(item, answers.pinned, key, lookup, pattern)
    # The other versions.env pins a step reads (steps 50 to 53): no freshness source exists for them, so
    # nothing is looked up; a pin here that differs from the run's (an image by its full reference, so a
    # new digest at the same tag counts) makes its step rerun-ready. A pin the run's manifest lacks is
    # unknown, never no-newer-release. A pin read only by steps this run did not run (an opt-in step not
    # in TOOLS, or one skipped for missing data) is left out: a rerun would not run that step either.
    known = {k for k, *_ in PINNED}
    skipped = steps_not_run(sample, sample_dir)
    for key in dict.fromkeys(r[:-5] if r.endswith("@self") else r for s in steps for r in reads_of(s)):
        if key in known or key not in this_pins:
            continue
        readers = [s["step"] for s in steps if any(key_matches(key, r) for r in reads_of(s))]
        if readers and all(n in skipped for n in readers):
            continue
        ran_ref, inst_ref = version.get(key) or None, this_pins.get(key) or None
        ran, inst = pin_value(key, ran_ref), pin_value(key, inst_ref)
        cls = "unknown" if not (ran_ref and inst_ref) else "rerun-ready" if ran_ref != inst_ref else "no-newer-release"
        item = {"item": key, "label": f"{key} (no upstream lookup)", "kind": "pin", "ran_with": ran,
                "installed": inst, "class": cls}
        if key.endswith("_IMAGE"):
            item["ran_with_ref"], item["installed_ref"] = ran_ref, inst_ref
        item["steps"] = affected_steps(steps, [key])
        item["action"] = action(item, sample)
        items.append(item)
    for (sec, key), value in sorted(man.items()):
        if sec != "data" or not key.startswith("pgs:"):
            continue
        pid = key[4:]
        now = header_value(os.path.join(genome_dir, "prs_scores", f"{pid}.txt.gz"), "HmPOS_date")
        add({"item": key, "label": f"PGS Catalog score {pid}", "kind": "file",
             "ran_with": date8(pgs_field(value, "HmPOS_date")), "installed": date8(now)}, answers.pgs, pid)
    act = [i["item"] for i in items if i["class"] in ("rerun-ready", "refresh-available")]
    return {
        "schema_version": SCHEMA_VERSION, "sample": sample, "checked_utc": now_utc(),
        "manifest_written_utc": man.get(("run", "written_utc")),
        "offline": isinstance(getattr(answers, "inner", answers), Offline), "hosts": list(answers.hosts),
        "pipeline_release": {"tag": tag} if not rel_err else {"error": rel_err},
        "items": items,
        "affected_steps": affected_steps(steps, act),
        "rerun_ready_steps": affected_steps(steps, [i["item"] for i in items if i["class"] == "rerun-ready"]),
    }


# --- the files it writes -------------------------------------------------------------------

def write_json(path, obj):
    d = os.path.dirname(path)
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".state.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, indent=1)
            f.write("\n")
        mask = os.umask(0)
        os.umask(mask)
        os.chmod(tmp, 0o666 & ~mask)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def append_log(path, rows):
    """One write per row, opened for append only: the log never shrinks."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    new = not os.path.exists(path)
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o666)
    try:
        if new:
            os.write(fd, LOG_HEADER.encode())
        for r in rows:
            os.write(fd, ("\t".join(str(x if x is not None else "").replace("\t", " ").replace("\n", " ")
                                    for x in r) + "\n").encode())
    finally:
        os.close(fd)


def print_check(state, out=None):
    w = (out or sys.stdout).write
    w(f"Database re-check of {state['sample']} ({state['checked_utc']}):\n")
    for i in state["items"]:
        detail = i.get("error") or f"run used {i.get('ran_with') or 'unknown'}, here {i.get('installed') or 'unknown'}, upstream {i.get('upstream') or 'unknown'}"
        if not i.get("error") and i.get("ran_with_ref") and i.get("installed_ref") and i["ran_with_ref"] != i["installed_ref"] \
                and i.get("ran_with") == i.get("installed"):
            detail += f"; same tag, new image: run used {i['ran_with_ref']}, here {i['installed_ref']}"
        w(f"  {i['label']}: {i['class']} ({detail})\n")
        if i["class"] not in ("no-newer-release",) and i.get("action"):
            w(f"    {i['action']}\n")
        if i["class"] in ("rerun-ready", "refresh-available") and i["steps"]:
            w(f"    steps it changes: {', '.join(i['steps'])} and the text report\n")
    if state["affected_steps"]:
        w(f"A rerun would redo steps {', '.join(state['affected_steps'])} and the text report.\n")
    else:
        w("Nothing here is newer than what the run used.\n")


# --- the rerun -----------------------------------------------------------------------------

def redone_steps(rs):
    """The steps of which a task ran in this run (COMPLETED, not taken from the cache)."""
    return [s["step"] for s in rs["steps"] if any(c.get("done", 0) > 0 for c in s["processes"].values())]


def compare_rerun(predicted, rs, steps):
    by = {s["step"]: s for s in steps}
    redone = redone_steps(rs)
    script_steps = [n for n in predicted if not by.get(n, {}).get("processes")]
    return {
        "predicted": predicted, "redone": redone,
        "predicted_not_redone": [n for n in predicted if n not in redone and n not in script_steps],
        "redone_not_predicted": [n for n in redone if n not in predicted],
        "not_redone_script_steps": script_steps,
    }


def rerun(a, genome_dir, sample, steps, state, log):
    sample_dir = os.path.join(genome_dir, sample)
    os.makedirs(os.path.join(sample_dir, "recheck"), exist_ok=True)
    # One --rerun per sample at a time: held from the running check through the end of the replay,
    # so two of them cannot both pass the check and run over each other's run_status.tsv.
    with open(os.path.join(sample_dir, "recheck", "rerun.lock"), "a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            print(f"ERROR: another --rerun of {sample} is going on: not rerunning.", file=sys.stderr)
            return 3
        return _rerun_locked(sample, sample_dir, steps, state, log)


def _rerun_locked(sample, sample_dir, steps, state, log):
    rs = run_state.collect(sample, sample_dir)
    if rs["run"]["state"] == "running":
        print(f"ERROR: a run of {sample} is going on (logs/run_status.tsv has no finished_utc): not rerunning. "
              "When it has ended, run this again.", file=sys.stderr)
        return 3
    meta = read_status_meta(sample_dir)
    cmd, names = meta.get("command"), meta.get("switches", "")
    if not cmd:
        print("ERROR: logs/run_status.tsv names no run-all.sh command (a run of an older version, or a bare "
              "Nextflow run): run the same command again by hand.", file=sys.stderr)
        return 2
    env = {k: v for k, v in os.environ.items() if k not in names.split()}
    started = now_utc()
    predicted = state["rerun_ready_steps"]  # after --refresh-clinvar the new file is installed, so rerun-ready
    print(f"Rerunning: {cmd}")
    append_log(log, [(started, "rerun-start", "", "", cmd)])
    rc = subprocess.run(["bash", "-c", cmd], env=env).returncode
    after = run_state.collect(sample, sample_dir)
    cmp = compare_rerun(predicted, after, steps)
    state["rerun"] = {"started_utc": started, "finished_utc": now_utc(), "rc": rc, "command": cmd,
                      "predicted_vs_redone": cmp}
    append_log(log, [(state["rerun"]["finished_utc"], "rerun-end", "", f"rc {rc}",
                      "predicted " + ",".join(cmp["predicted"]) + "; redone " + ",".join(cmp["redone"]))])
    print(f"Predicted steps: {', '.join(cmp['predicted']) or 'none'}; redone: {', '.join(cmp['redone']) or 'none'}")
    for k, what in (("predicted_not_redone", "named but not redone"), ("redone_not_predicted", "redone but not named"),
                    ("not_redone_script_steps", "script steps not redone (they skip when their output exists; "
                                                "delete it and run the script)")):
        if cmp[k]:
            print(f"  {what}: {', '.join(cmp[k])}")
    return rc


def launch_dir_command(d):
    """The last command in DIR/.nextflow/history, with -resume."""
    p = os.path.join(d, ".nextflow", "history")
    try:
        lines = [x for x in open(p, errors="replace").read().splitlines() if x.strip()]
    except OSError:
        return None
    if not lines:
        return None
    cmd = lines[-1].split("\t")[-1].strip()
    try:
        words = shlex.split(cmd)
    except ValueError:
        words = cmd.split()
    return cmd if "-resume" in words else cmd + " -resume"


def running_samples(genome_dir):
    """The samples of genome_dir (a folder with logs/run_status.tsv) whose run is going on."""
    try:
        names = sorted(os.listdir(genome_dir))
    except OSError:
        return []
    return [n for n in names if os.path.isfile(os.path.join(genome_dir, n, "logs", "run_status.tsv"))
            and run_state.collect(n, os.path.join(genome_dir, n))["run"]["state"] == "running"]


def refresh_clinvar(genome_dir):
    """setup.sh --refresh clinvar under the ClinVar lock of GENOME_DIR (CLINVAR_LOCK). The ClinVar files
    are shared by every sample, so a run of any of them blocks the refresh: run-all.sh and step 06 hold
    the lock shared while they read the files, and wait while this holds it. A run with no lock (a bare
    nextflow run, an older checkout) is caught by its logs/run_status.tsv. The lock is let go before
    --rerun starts run-all.sh, which takes it shared."""
    os.makedirs(os.path.join(genome_dir, "clinvar"), exist_ok=True)
    with open(os.path.join(genome_dir, CLINVAR_LOCK), "a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            print(f"ERROR: a run or step 06 holds {os.path.join(genome_dir, CLINVAR_LOCK)}: not refreshing the "
                  "ClinVar files it reads.", file=sys.stderr)
            return 3
        busy = running_samples(genome_dir)
        if busy:
            print(f"ERROR: a run of {', '.join(busy)} is going on: not refreshing the ClinVar files it reads.",
                  file=sys.stderr)
            return 3
        r = subprocess.run([os.path.join(ROOT, "scripts", "setup.sh"), "--refresh", "clinvar", genome_dir]).returncode
        if r != 0:
            print("ERROR: setup.sh --refresh clinvar failed: see above.", file=sys.stderr)
        return r


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample")
    ap.add_argument("--genome-dir", default=os.environ.get("GENOME_DIR"), help="default: $GENOME_DIR")
    ap.add_argument("--offline", help="read the upstream answers from this JSON file")
    ap.add_argument("--owngenome-repo", default=DEFAULT_REPO, help=f"the pipeline's repository (default {DEFAULT_REPO})")
    ap.add_argument("--pgs-api", default=PGS_API, help=argparse.SUPPRESS)
    ap.add_argument("--raw", default=os.environ.get("OWNGENOME_RECHECK_RAW", RAW), help=argparse.SUPPRESS)
    ap.add_argument("--rerun", action="store_true", help="replay the last run-all.sh command")
    ap.add_argument("--refresh-clinvar", action="store_true", help="with --rerun: setup.sh --refresh clinvar first")
    ap.add_argument("--affected", help="print the steps a change of these manifest keys redoes")
    ap.add_argument("--launch-dir", help="print the last nextflow run of this launch directory with -resume")
    ap.add_argument("--steps", default=STEPS_TSV, help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    steps = run_state.read_steps(a.steps)
    if a.affected is not None:
        print(" ".join(affected_steps(steps, [k for k in a.affected.split(",") if k])))
        return 0
    if a.launch_dir:
        cmd = launch_dir_command(a.launch_dir)
        if not cmd:
            print(f"ERROR: no .nextflow/history in {a.launch_dir}", file=sys.stderr)
            return 2
        print("The last run there, with -resume (not run; check it, then run it from that directory):")
        print(cmd)
        return 0
    if not a.sample or not a.genome_dir:
        ap.error("give --sample and --genome-dir (or GENOME_DIR)")
    if a.refresh_clinvar and not a.rerun:
        ap.error("--refresh-clinvar goes with --rerun")
    sample_dir = os.path.join(a.genome_dir, a.sample)
    if not os.path.isfile(os.path.join(sample_dir, "run_manifest.tsv")):
        print(f"ERROR: no {sample_dir}/run_manifest.tsv: run the pipeline on {a.sample} first.", file=sys.stderr)
        return 2
    if a.refresh_clinvar:
        r = refresh_clinvar(a.genome_dir)
        if r != 0:
            return r
    fr = load_freshness()
    this_pins, _ = fr.read_pins(os.path.join(ROOT, "versions.env"))
    try:
        this_hla = fr.hla_db_release()
    except fr.LookupFailed:
        this_hla = None
    if a.offline:
        answers = Offline(a.offline)
    else:
        man = read_manifest(os.path.join(sample_dir, "run_manifest.tsv"))
        pins = {k: v for (sec, k), v in man.items() if sec == "version"}
        answers = Online(fr, pins, a.owngenome_repo, a.pgs_api, a.raw)
    answers = Remembered(answers)
    state = check(a.sample, a.genome_dir, answers, steps, this_pins, this_hla)
    d = os.path.join(sample_dir, "recheck")
    log = os.path.join(d, "log.tsv")
    try:
        with open(os.path.join(d, "state.json")) as f:
            old = json.load(f)
        if old.get("rerun"):
            state["rerun"] = old["rerun"]  # the last rerun stays until another one replaces it
    except (OSError, ValueError):
        pass
    append_log(log, [(state["checked_utc"], "check", i["item"], i["class"],
                      i.get("error") or f"{i.get('ran_with')} {i.get('installed')} {i.get('upstream')}")
                     for i in state["items"]])
    write_json(os.path.join(d, "state.json"), state)
    print_check(state)
    print(f"Written: {os.path.join(d, 'state.json')}")
    errors = [i for i in state["items"] if i["class"] == "error"]
    if a.rerun:
        before = state.get("rerun")
        rc = rerun(a, a.genome_dir, a.sample, steps, state, log)
        if state.get("rerun") is not before:
            # The replay wrote a new run_manifest.tsv: check again against it, with the same answers,
            # so state.json no longer calls ready what the rerun has just used.
            fresh = check(a.sample, a.genome_dir, answers, steps, this_pins, this_hla)
            fresh["rerun"] = state["rerun"]
            state = fresh
            append_log(log, [(state["checked_utc"], "check", i["item"], i["class"],
                              i.get("error") or f"{i.get('ran_with')} {i.get('installed')} {i.get('upstream')}")
                             for i in state["items"]])
        write_json(os.path.join(d, "state.json"), state)
        if rc:
            return rc
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
