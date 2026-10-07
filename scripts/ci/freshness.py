#!/usr/bin/env python3
"""freshness.py: say which pinned databases, bundles and images have fallen behind.

Renovate can bump an image tag, but it cannot see the Ensembl cache release,
the PCGR bundle date, the pypgx-bundle tag, the Nextflow line or a ClinVar file
that has gone a month without a refresh. This script reads the pins from
versions.env (and the ClinVar URL from scripts/setup.sh), asks each upstream,
and prints a Markdown report. The monthly workflow (.github/workflows/
freshness.yml) puts that report into one issue labelled `freshness`.

Every lookup fails closed: an HTTP error, a timeout, an empty answer or an
answer without the expected field is an ERROR row, never "current", and the
script then exits 1. Being behind is not an error; it is what the report is for.

Usage:
  scripts/ci/freshness.py [--out FILE] [--sections a,b] [--versions F] [--setup F]
      sections: databases, clinvar (the ClinVar row alone), images, digests
  scripts/ci/freshness.py --update-issue BODY_FILE --repo OWNER/NAME
  scripts/ci/freshness.py --self-test

Standard library only. GITHUB_TOKEN, when set, is sent to api.github.com (and
is required by --update-issue).
"""

import argparse
import datetime as dt
import email.utils
import http.client
import http.server
import json
import os
import re
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
UA = "Personal-Genome-Pipeline freshness check"
CLINVAR_MAX_AGE_DAYS = 35
ISSUE_LABEL = "freshness"
ISSUE_MARKER = "<!-- freshness-report: edited by .github/workflows/freshness.yml -->"
ISSUE_TITLE = "Freshness report: databases, bundles, images and links"
# pypgx 0.27.0 pulled pandas 3.0 and broke every gene (versions.env hold note),
# so it never counts as a version the held pin could move to.
PYPGX_BROKEN = {"0.27.0"}
SECTIONS = ("databases", "clinvar", "images", "digests")  # clinvar: the ClinVar row alone

TIMEOUT = 30
RETRIES = 3


class LookupFailed(Exception):
    """A lookup that gave no usable answer."""


# An answer in an unexpected shape (a list where a dict was expected, a missing
# key) is a failed lookup too: it becomes an error row for that lookup instead
# of a traceback that loses the whole report.
LOOKUP_ERRORS = (LookupFailed, KeyError, IndexError, TypeError, ValueError, AttributeError)


# --------------------------------------------------------------------------
# HTTP

def request(url, method="GET", headers=None, data=None, want="json", retries=None, timeout=None):
    """Return (parsed body, response headers). Fails closed on anything odd.

    want: "json" (a non-empty JSON document), "text" (non-empty text) or
    "head" (no body; the headers are the answer). 4xx answers other than 429
    fail at once; 5xx, 429, timeouts and connection errors are retried.
    """
    retries = RETRIES if retries is None else retries
    timeout = TIMEOUT if timeout is None else timeout
    hdrs = {"User-Agent": UA}
    hdrs.update(headers or {})
    last = None
    for attempt in range(retries + 1):
        if attempt:
            time.sleep(5 * attempt)
        req = urllib.request.Request(url, method=method, headers=hdrs, data=data)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                raw = resp.read() if want != "head" else b""
                resp_headers = resp.headers
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}"
            if 400 <= e.code < 500 and e.code != 429:
                break
            continue
        except (urllib.error.URLError, TimeoutError, OSError, http.client.HTTPException) as e:
            last = f"no answer ({getattr(e, 'reason', e)})"
            continue
        if want == "head":
            return None, resp_headers
        text = raw.decode("utf-8", "replace")
        if not text.strip():
            last = "empty answer"
            continue
        if want == "text":
            return text, resp_headers
        try:
            return json.loads(text), resp_headers
        except ValueError:
            last = "answer is not JSON"
            continue
    raise LookupFailed(f"{url}: {last}")


def gh_api(path, **kw):
    headers = {"Accept": "application/vnd.github+json"}
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    base = os.environ.get("FRESHNESS_GITHUB_API", "https://api.github.com")
    return request(base + path, headers=headers, **kw)[0]


def need(value, what):
    """An empty value from a lookup is an error, never 'up to date'."""
    if value in (None, "", [], {}):
        raise LookupFailed(f"{what}: empty")
    return value


# --------------------------------------------------------------------------
# Pins

def read_pins(path):
    """Return ({NAME: value}, {NAME: note}) from versions.env.

    The note is the trailing comment on the line, or a '# hold:' comment block
    directly above it.
    """
    pins, notes, block = {}, {}, []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("#"):
                block.append(s.lstrip("#").strip())
                continue
            m = re.match(r'^([A-Z][A-Z0-9_]*)="([^"]*)"\s*(?:#\s*(.*))?$', s)
            if m:
                name, value, trailing = m.groups()
                pins[name] = value
                note = trailing or ""
                if not note and block and block[0].startswith("hold:"):
                    note = " ".join(block)
                if note:
                    notes[name] = first_sentence(note)
            block = []
    return pins, notes


def first_sentence(text):
    m = re.match(r"(.+?[.;])(\s|$)", text)
    return (m.group(1) if m else text).rstrip(".;")


def pin(pins, name):
    if not pins.get(name):
        raise LookupFailed(f"{name} is not set in versions.env")
    return pins[name]


# --------------------------------------------------------------------------
# Versions. Never `sort -V` over raw tag lists: it ranks conda build strings
# and suffixed tags (3.14-slim, 1.10.0-gpu, 3.15.0rc1) as versions. A tag is
# compared only with tags of the same shape as the pin.

BIOCONDA = re.compile(r"^(?P<v>\d+(?:\.\d+)*)--(?P<hash>[A-Za-z0-9]+)_(?P<build>\d+)$")


def tag_key(tag, like):
    """Sort key for `tag` when it has the same shape as the pinned tag `like`;
    None when it does not."""
    pm = BIOCONDA.match(like)
    if pm:
        m = BIOCONDA.match(tag)
        if not m:
            return None
        return (tuple(int(x) for x in m.group("v").split(".")), int(m.group("build")))
    if not re.search(r"\d", like) or re.sub(r"\d+", "N", tag) != re.sub(r"\d+", "N", like):
        return None
    return (tuple(int(x) for x in re.findall(r"\d+", tag)), 0)


def newer_tags(pin_tag, tags):
    """Return (newest tag overall, newest tag on the pin's major), each newer
    than the pin or None."""
    pk = tag_key(pin_tag, pin_tag)
    if pk is None:
        raise LookupFailed(f"cannot read a version from the pinned tag {pin_tag}")
    keyed = [(k, t) for t in tags for k in [tag_key(t, pin_tag)] if k is not None and k > pk]
    if not keyed:
        return None, None
    newest = max(keyed)[1]
    same = [(k, t) for k, t in keyed if k[0][:1] == pk[0][:1]]
    return newest, (max(same)[1] if same else None)


def change_kind(pin_tag, new_tag):
    a, b = tag_key(pin_tag, pin_tag), tag_key(new_tag, pin_tag)
    if a[0][:1] != b[0][:1]:
        return "major"
    if a[0][:2] != b[0][:2]:
        return "minor"
    if a[0] != b[0]:
        return "patch"
    return "build"


def ver(text):
    return tuple(int(x) for x in re.findall(r"\d+", text))


# --------------------------------------------------------------------------
# Container registries (the v2 API every registry speaks)

def split_image(image):
    """'quay.io/biocontainers/x:1--h_0' -> ('quay.io', 'biocontainers/x', '1--h_0', None);
    'python:3.11.17@sha256:ab' -> ('docker.io', 'library/python', '3.11.17', 'sha256:ab'): a
    tag pinned to a digest keeps its tag; 'a/b@sha256:ab' (digest only) has tag None."""
    digest = None
    if "@" in image:
        image, digest = image.split("@", 1)
    if ":" in image.split("/")[-1]:
        image, _, tag = image.rpartition(":")
    else:
        tag = None if digest else "latest"
    parts = image.split("/")
    if len(parts) > 1 and ("." in parts[0] or ":" in parts[0]):
        host, repo = parts[0], "/".join(parts[1:])
    else:
        host, repo = "docker.io", image
    if host == "docker.io" and "/" not in repo:
        repo = "library/" + repo
    return host, repo, tag, digest


def registry(host):
    return "registry-1.docker.io" if host == "docker.io" else host


def registry_headers(host, repo):
    if host == "docker.io":
        tok = request("https://auth.docker.io/token?service=registry.docker.io&scope="
                   f"repository:{repo}:pull")[0].get("token")
    elif host == "ghcr.io":
        tok = request(f"https://ghcr.io/token?scope=repository:{repo}:pull")[0].get("token")
    else:
        return {}
    return {"Authorization": "Bearer " + need(tok, f"registry token for {host}/{repo}")}


def list_tags(host, repo):
    headers = registry_headers(host, repo)
    url = f"https://{registry(host)}/v2/{repo}/tags/list?n=1000"
    origin = urllib.parse.urlsplit(url)[:2]
    tags = []
    for _ in range(200):
        body, h = request(url, headers=headers)
        tags += body.get("tags") or []
        m = re.search(r'<([^>]+)>;\s*rel="next"', h.get("Link", "") or "")
        if not m:
            break
        url = urllib.parse.urljoin(url, m.group(1))
        # The registry token goes only to the registry that issued it.
        if urllib.parse.urlsplit(url)[:2] != origin:
            raise LookupFailed(f"tag list of {host}/{repo} points its next page at another host: {url}")
    return need(tags, f"tag list of {host}/{repo}")


def manifest_digest(host, repo, ref):
    headers = registry_headers(host, repo)
    headers["Accept"] = ", ".join([
        "application/vnd.oci.image.index.v1+json",
        "application/vnd.docker.distribution.manifest.list.v2+json",
        "application/vnd.oci.image.manifest.v1+json",
        "application/vnd.docker.distribution.manifest.v2+json",
    ])
    _, h = request(f"https://{registry(host)}/v2/{repo}/manifests/{ref}", method="HEAD",
                headers=headers, want="head")
    return need(h.get("Docker-Content-Digest"), f"digest of {host}/{repo}:{ref}")


# --------------------------------------------------------------------------
# Report rows

def esc(text):
    return str(text).replace("|", "\\|").replace("\n", " ")


class Report:
    def __init__(self):
        self.db_rows = []        # (item, pinned, upstream, state, todo)
        self.image_rows = []     # (var, pinned, newest, kind, note)
        self.current = []        # var names
        self.rebuilds = []       # var names: only a new build of the pinned version
        self.line_rows = []      # (var, pinned, newest): tags that float on a release line
        self.digest_rows = []    # (var, pinned, publisher, state)
        self.errors = []         # (what, message)

    def error(self, what, exc):
        msg = str(exc) if isinstance(exc, LookupFailed) else f"unexpected answer ({type(exc).__name__}: {exc})"
        self.errors.append((what, msg))

    def db(self, item, pinned, upstream, state, todo=""):
        self.db_rows.append((item, pinned, upstream, state, todo))


def run_lookup(report, what, fn, *args):
    try:
        fn(report, *args)
    except LOOKUP_ERRORS as e:
        report.error(what, e)
        report.db(what, "", "", "ERROR", "lookup failed, see Errors")


def clinvar(report, setup_path):
    text = open(setup_path).read()
    m = re.search(r'CLINVAR_URL="(https?://[^"]+)"', text)
    if not m:
        raise LookupFailed(f"no CLINVAR_URL in {setup_path}")
    url = m.group(1)
    _, h = request(url, method="HEAD", want="head")
    lm = need(h.get("Last-Modified"), f"Last-Modified of {url}")
    when = email.utils.parsedate_to_datetime(lm)
    age = (dt.datetime.now(dt.timezone.utc) - when).days
    state = "current" if age <= CLINVAR_MAX_AGE_DAYS else "behind"
    report.db("ClinVar (steps 6, 23)", "downloaded once by setup.sh",
              f"{url.rsplit('/', 1)[-1]} dated {when:%Y-%m-%d} ({age} days ago)", state,
              f"NCBI refreshes it monthly. Refresh a local copy older than {CLINVAR_MAX_AGE_DAYS} days "
              "(its `##fileDate`): delete `$GENOME_DIR/clinvar/` and rerun `./scripts/setup.sh`."
              + ("" if state == "current" else
                 f" Upstream itself is older than {CLINVAR_MAX_AGE_DAYS} days: check that the URL still is the live file."))


def vep_cache_template():
    text = open(os.path.join(ROOT, "scripts", "lib", "common.sh")).read()
    m = re.search(r"vep_cache_url\(\)\s*\{\s*printf '([^']+)'", text)
    if not m:
        raise LookupFailed("no vep_cache_url template in scripts/lib/common.sh")
    return m.group(1)


def ensembl(report, pins):
    cache = pin(pins, "VEP_CACHE_RELEASE")
    image_tag = split_image(pin(pins, "VEP_IMAGE"))[2]
    releases = need(request("https://rest.ensembl.org/info/data/?content-type=application/json")[0]
                    .get("releases"), "Ensembl releases")
    current = max(int(r) for r in releases)
    url = vep_cache_template().replace("%s", cache)
    request(url, method="HEAD", want="head")
    state = "current" if current <= int(cache) else "behind"
    report.db("Ensembl / VEP cache (step 13)", f"VEP_CACHE_RELEASE={cache}, VEP_IMAGE {image_tag}",
              f"Ensembl release {current}", state,
              "" if state == "current" else
              f"Move VEP_IMAGE to release_{current}.x and VEP_CACHE_RELEASE to {current} in one change "
              "(a major: Renovate waits for approval); rerun steps 13, 30, 23 and 31.")


def pcgr(report, pins):
    tag = split_image(pin(pins, "PCGR_IMAGE"))[2]
    bundle, vep = pin(pins, "PCGR_DATA_BUNDLE"), pin(pins, "PCGR_VEP_CACHE_RELEASE")
    latest = need(gh_api("/repos/sigven/pcgr/releases/latest").get("tag_name"), "PCGR latest release")
    version = latest.lstrip("v")
    src = request(f"https://raw.githubusercontent.com/sigven/pcgr/{latest}/pcgr/pcgr_vars.py", want="text")[0]
    db = re.search(r"^DB_VERSION\s*=\s*'(\d+)'", src, re.M)
    vv = re.search(r"^VEP_VERSION\s*=\s*'(\d+)'", src, re.M)
    if not db or not vv:
        raise LookupFailed(f"no DB_VERSION or VEP_VERSION in pcgr_vars.py at {latest}")
    where = [name for name, host, repo in (("Docker Hub sigven/pcgr", "docker.io", "sigven/pcgr"),
                                            ("ghcr.io/sigven/pcgr", "ghcr.io", "sigven/pcgr"))
             if version in list_tags(host, repo)]
    behind = ver(version) > ver(tag) or db.group(1) != bundle or vv.group(1) != vep
    report.db("PCGR / CPSR (step 17)", f"{tag}, bundle {bundle}, VEP {vep}",
              f"{version}, bundle {db.group(1)}, VEP {vv.group(1)}; image on "
              + (", ".join(where) if where else "no registry yet"),
              "behind" if behind else "current",
              "" if not behind else
              "PCGR_IMAGE, PCGR_DATA_BUNDLE and PCGR_VEP_CACHE_RELEASE move together; rerun step 17."
              + ("" if where else " Wait for a published image."))


def pypgx(report, pins):
    held = split_image(pin(pins, "PYPGX_IMAGE"))[2].split("--")[0]
    bundle = pin(pins, "PYPGX_BUNDLE_VERSION")
    bundle_tags = [t["name"] for t in need(gh_api("/repos/sbslee/pypgx-bundle/tags?per_page=100"),
                                           "pypgx-bundle tags")]
    bundle_versions = {t.lstrip("v") for t in bundle_tags if re.fullmatch(r"v?\d+\.\d+\.\d+", t)}
    if bundle not in bundle_versions:
        raise LookupFailed(f"PYPGX_BUNDLE_VERSION {bundle} is not a pypgx-bundle tag")
    pypi = need(request("https://pypi.org/pypi/pypgx/json")[0].get("info", {}).get("version"), "pypgx on PyPI")
    images = need({m.group("v") for t in list_tags("quay.io", "biocontainers/pypgx")
                   for m in [BIOCONDA.match(t)] if m}, "pypgx biocontainer versions")
    newest_bundle = max(bundle_versions, key=ver)
    movable = sorted((v for v in images & bundle_versions
                      if ver(v) > ver(held) and v not in PYPGX_BROKEN), key=ver)
    report.db("pypgx and pypgx-bundle (step 32)", f"pypgx {held} (held), bundle {bundle}",
              f"PyPI {pypi}, newest biocontainer {max(images, key=ver)}, newest bundle tag {newest_bundle}",
              f"can move to {movable[-1]}" if movable else "held, cannot move yet",
              (f"pypgx {movable[-1]} has an image and a bundle tag: move PYPGX_IMAGE and "
               "PYPGX_BUNDLE_VERSION together and rerun step 32." if movable else
               "The hold waits for a pypgx release past " + ", ".join(sorted(PYPGX_BROKEN))
               + " with a biocontainer and a pypgx-bundle tag of the same version."))


def pharmcat(report, pins):
    tag = split_image(pin(pins, "PHARMCAT_IMAGE"))[2]
    latest = need(gh_api("/repos/PharmGKB/PharmCAT/releases/latest").get("tag_name"),
                  "PharmCAT latest release").lstrip("v")
    has_image = latest in list_tags("docker.io", "pgkb/pharmcat")
    behind = ver(latest) > ver(tag)
    report.db("PharmCAT (steps 7, 27)", tag,
              f"{latest}" + ("" if has_image else " (no pgkb/pharmcat image yet)"),
              "behind" if behind else "current",
              "" if not behind else "Revalidate steps 7 and 27 on a known sample before moving PHARMCAT_IMAGE.")


def imgt(report, pins):
    text = request("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/release_version.txt", want="text")[0]
    v = re.search(r"#\s*version:\s*IPD-IMGT/HLA\s+(\S+)", text)
    d = re.search(r"#\s*date:\s*(\S+)", text)
    if not v:
        raise LookupFailed("no version line in IMGTHLA release_version.txt")
    report.db("IPD-IMGT/HLA (step 8)", "whatever was current when the T1K index was built",
              f"{v.group(1)}" + (f" ({d.group(1)})" if d else ""), "report only",
              "To type against this release, delete `$GENOME_DIR/t1k_idx/` and rerun step 8; "
              "the release is in the header of `t1k_idx/hlaidx/hla.dat`.")


def gnomad(report, versions_path):
    m = re.search(r"gnomad_v(\d+(?:\.\d+)*)_constraint", open(versions_path).read())
    if not m:
        raise LookupFailed("no gnomad_v<version>_constraint hold note in versions.env")
    pinned = m.group(1)
    api = "https://storage.googleapis.com/storage/v1/b/gcp-public-data--gnomad/o"
    prefixes = need(request(f"{api}?prefix=release/&delimiter=/")[0].get("prefixes"), "gnomAD releases")
    found = []
    for p in prefixes:
        rel = p.split("/")[1]
        if re.fullmatch(r"\d+(\.\d+)*", rel) and ver(rel) > ver(pinned):
            items = request(f"{api}?prefix=release/{rel}/constraint/&maxResults=1")[0].get("items")
            if items:
                found.append(rel)
    newest = max(found, key=ver) if found else None
    report.db("gnomAD constraint (steps 23, 31)", f"v{pinned} (held)",
              f"newest release with constraint data: {newest}" if newest else f"no release past v{pinned}",
              "report only",
              "Held: the 4.1.2 table has a different layout; read it and adapt steps 23 and 31 before moving."
              if newest else "")


def nextflow(report, pins):
    pinned = pin(pins, "NEXTFLOW_VERSION")
    rels = need(gh_api("/repos/nextflow-io/nextflow/releases?per_page=100"), "Nextflow releases")
    stable = [r["tag_name"].lstrip("v") for r in rels
              if not r.get("prerelease") and re.fullmatch(r"v?\d+\.\d+\.\d+", r["tag_name"])]
    need(stable, "stable Nextflow releases")
    line = ".".join(pinned.split(".")[:2])
    in_line = max((s for s in stable if s.startswith(line + ".")), key=ver, default=None)
    newer_lines = {}
    for s in stable:
        ln = ".".join(s.split(".")[:2])
        if ver(ln) > ver(line) and ver(s) > ver(newer_lines.get(ln, "0")):
            newer_lines[ln] = s
    newest_line = max(newer_lines.values(), key=ver) if newer_lines else None
    behind = in_line and ver(in_line) > ver(pinned)
    report.db("Nextflow", pinned,
              f"{line} line: {in_line or 'none'}; newest stable line: {newest_line or line}",
              "behind" if behind else ("newer line" if newest_line else "current"),
              ("Move NEXTFLOW_VERSION within the line. " if behind else "")
              + (f"Moving to {newest_line} is a new line: run the stub and e2e jobs on it first." if newest_line else ""))


def cyrius(report, pins):
    pinned = pin(pins, "CYRIUS_VERSION")
    latest = need(request("https://pypi.org/pypi/cyrius/json")[0].get("info", {}).get("version"), "Cyrius on PyPI")
    behind = ver(latest) > ver(pinned)
    report.db("Cyrius (step 21, legacy)", pinned, f"PyPI {latest}", "behind" if behind else "current",
              "" if not behind else "Check scripts/cyrius-constraints.txt still resolves before moving CYRIUS_VERSION.")


def images(report, pins, notes):
    for var in sorted(v for v in pins if v.endswith("_IMAGE")):
        image = pins[var]
        host, repo, tag, digest = split_image(image)
        if digest and not tag:
            continue  # digest only: the digests section reports it
        try:
            tags = list_tags(host, repo)
            if tag not in tags:
                raise LookupFailed(f"{image}: the pinned tag is not in the registry's tag list")
            newest, same_major = newer_tags(tag, tags)
        except LOOKUP_ERRORS as e:
            report.error(var, e)
            continue
        if newest is None:
            report.current.append(var)
            continue
        # python:3.11 follows 3.11.x by itself; 3.14 is a new line, a choice
        # rather than a missed update, so it is reported apart.
        if re.fullmatch(r"\d+\.\d+", tag) and any(re.fullmatch(re.escape(tag) + r"\.\d+", t) for t in tags):
            report.line_rows.append((var, tag, newest))
            continue
        if change_kind(tag, newest) == "build":
            report.rebuilds.append(var)
            continue
        if same_major and same_major != newest:
            report.image_rows.append((var, tag, same_major, change_kind(tag, same_major), notes.get(var, "")))
        report.image_rows.append((var, tag, newest, change_kind(tag, newest), notes.get(var, "")))


def digests(report, pins):
    """A tag pinned to a digest (name:tag@sha256:...) is compared with what its
    own tag points at now; a digest-only pin with the publisher's `latest`."""
    for var in sorted(v for v in pins if v.endswith("_IMAGE")):
        host, repo, tag, digest = split_image(pins[var])
        if not digest:
            continue
        try:
            now = manifest_digest(host, repo, tag or "latest")
            versioned = [] if tag else [t for t in list_tags(host, repo) if re.search(r"\d", t)]
        except LOOKUP_ERRORS as e:
            report.error(var, e)
            continue
        if tag:
            state = (f"`{tag}` still points at the pinned digest" if now == digest else
                     f"`{tag}` was rebuilt: move the digest after the image test passes on it")
        else:
            state = "same as publisher's latest" if now == digest else "publisher's latest moved"
        if versioned:
            state += "; versioned tags now exist: " + ", ".join(sorted(versioned, key=ver)[-3:])
        pinned = f"{repo}:{tag}@{digest[:19]}" if tag else f"{repo}@{digest[:19]}"
        report.digest_rows.append((var, pinned, f"{tag or 'latest'}: {now[:19]}", state))


def build(pins, notes, versions_path, setup_path, sections):
    report = Report()
    if "databases" in sections or "clinvar" in sections:
        run_lookup(report, "ClinVar", clinvar, setup_path)
    if "databases" in sections:
        for what, fn in (("Ensembl / VEP cache", ensembl), ("PCGR / CPSR", pcgr), ("pypgx", pypgx),
                         ("PharmCAT", pharmcat), ("IPD-IMGT/HLA", imgt), ("Nextflow", nextflow),
                         ("Cyrius", cyrius)):
            run_lookup(report, what, fn, pins)
        run_lookup(report, "gnomAD constraint", gnomad, versions_path)
    if "images" in sections:
        images(report, pins, notes)
    if "digests" in sections:
        digests(report, pins)
    return report


def render(report, sections):
    out = []
    if "databases" in sections or "clinvar" in sections:
        out += ["### Databases, bundles and tools Renovate cannot see", "",
                "| Item | Pinned | Upstream | State | What to do |", "|---|---|---|---|---|"]
        out += ["| " + " | ".join(esc(c) for c in row) + " |" for row in report.db_rows]
        out.append("")
    if "images" in sections:
        majors = [r for r in report.image_rows if r[3] == "major"]
        rest = [r for r in report.image_rows if r[3] != "major"]
        out += ["### Images behind on a major version", "",
                "Renovate never opens a major update on its own: it waits for approval on its "
                "dependency dashboard, or renovate.json disables the image. Each needs its "
                "revalidation run.", ""]
        if majors:
            out += ["| Variable | Pinned | Newest | Note |", "|---|---|---|---|"]
            out += [f"| {v} | `{p}` | `{n}` | {esc(note)} |" for v, p, n, _, note in majors]
        else:
            out.append("None.")
        out += ["", "### Images behind on a minor, patch or build", "",
                "Renovate opens most of these as pull requests. The images renovate.json disables, "
                "holds or sends to the dependency dashboard (PCGR, Python) wait there. They are "
                "listed so a stalled or silent Renovate shows.", ""]
        if rest:
            out += ["| Variable | Pinned | Newest | Change | Note |", "|---|---|---|---|---|"]
            out += [f"| {v} | `{p}` | `{n}` | {k} | {esc(note)} |" for v, p, n, k, note in rest]
        else:
            out.append("None.")
        out += [""]
        if report.line_rows:
            out += ["Floating on a release line (the tag follows its own patch releases; a newer line is "
                    "a deliberate move, so it is not counted as behind): "
                    + ", ".join(f"{v} `{p}` (newest line `{n}`)" for v, p, n in report.line_rows) + ".", ""]
        if report.rebuilds:
            out += ["Only a new build of the pinned version (not listed above): "
                    + ", ".join(report.rebuilds) + ".", ""]
        out += [f"Current: {len(report.current)} images ({', '.join(report.current) or 'none'}).", ""]
    if "digests" in sections:
        out += ["### Images pinned by digest", "",
                "| Variable | Pinned | What the tag points at now (`latest` for a digest-only pin) | State |",
                "|---|---|---|---|"]
        out += [f"| {v} | `{p}` | `{d}` | {esc(s)} |" for v, p, d, s in report.digest_rows]
        out.append("")
    if report.errors:
        out += ["### Errors", "", "A lookup that fails is an error, never \"current\". This run exits 1.", "",
                "| Lookup | Error |", "|---|---|"]
        out += [f"| {esc(w)} | {esc(m)} |" for w, m in report.errors]
        out.append("")
    return "\n".join(out)


# --------------------------------------------------------------------------
# The one issue

def update_issue(repo, body_path, api=None):
    """Edit the one open issue that carries ISSUE_MARKER; reopen the newest
    closed one when none is open; create it when none exists. Issues with the
    label but without the marker (another workflow's) are never touched."""
    if api:
        os.environ["FRESHNESS_GITHUB_API"] = api
    if not os.environ.get("GITHUB_TOKEN"):
        raise LookupFailed("--update-issue needs GITHUB_TOKEN")
    body = open(body_path).read()
    if ISSUE_MARKER not in body:
        body = ISSUE_MARKER + "\n" + body
    if len(body) > 65000:
        body = body[:64000] + "\n\n(Cut at 64,000 characters; the full report is in the run log.)\n"
    try:
        gh_api(f"/repos/{repo}/labels/{ISSUE_LABEL}", retries=1)
    except LookupFailed as e:
        if "HTTP 404" not in str(e):
            raise
        gh_api(f"/repos/{repo}/labels", method="POST", data=json.dumps({
            "name": ISSUE_LABEL, "color": "c5def5",
            "description": "Monthly report of pins and databases behind upstream"}).encode())
    found, page = [], 1
    while True:
        batch = gh_api(f"/repos/{repo}/issues?labels={ISSUE_LABEL}&state=all&per_page=100&page={page}")
        if not isinstance(batch, list):
            raise LookupFailed("issue list is not a list")
        found += [i for i in batch if "pull_request" not in i and ISSUE_MARKER in (i.get("body") or "")]
        if len(batch) < 100:
            break
        page += 1
    open_ = sorted((i for i in found if i["state"] == "open"), key=lambda i: i["number"])
    closed = sorted((i for i in found if i["state"] != "open"), key=lambda i: i["number"])
    payload = {"title": ISSUE_TITLE, "body": body}
    if open_:
        target, action = open_[0], "edited"
        if len(open_) > 1:
            print(f"WARNING: {len(open_)} open freshness issues; editing #{target['number']} only.")
    elif closed:
        target, action = closed[-1], "reopened and edited"
        payload["state"] = "open"
    else:
        target, action = None, "created"
    data = json.dumps(dict(payload, labels=[ISSUE_LABEL]) if target is None else payload).encode()
    if target is None:
        # No retry: a create whose answer was lost may have landed, and a
        # second POST would open a second issue. A failed create fails the run.
        res = gh_api(f"/repos/{repo}/issues", method="POST", data=data, retries=0)
    else:
        res = gh_api(f"/repos/{repo}/issues/{target['number']}", method="PATCH", data=data)
    print(f"Issue {action}: {need(res.get('html_url'), 'issue URL')}")
    return action, res["number"]


# --------------------------------------------------------------------------
# Self-test

class FakeGitHub(http.server.BaseHTTPRequestHandler):
    """Just enough of the issues API for update_issue, kept in memory."""
    issues = []
    labels = set()
    lose_create_answer = False

    def log_message(self, *a):
        pass

    def reply(self, code, obj):
        raw = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def body(self):
        return json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")

    def do_GET(self):
        path = urllib.parse.urlparse(self.path)
        if path.path == "/empty":
            self.send_response(200)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        if path.path == "/truncated":
            # Promises 100 bytes and sends 10, so the read ends in IncompleteRead.
            self.send_response(200)
            self.send_header("Content-Length", "100")
            self.end_headers()
            self.wfile.write(b'{"a": 1, "')
            return
        if path.path.startswith("/repos/o/r/labels/"):
            name = path.path.rsplit("/", 1)[1]
            return self.reply(200, {"name": name}) if name in self.labels else self.reply(404, {})
        if path.path == "/repos/o/r/issues":
            return self.reply(200, [i for i in self.issues if ISSUE_LABEL in i["labels"]])
        self.reply(404, {})

    def do_POST(self):
        b = self.body()
        if self.path == "/repos/o/r/labels":
            self.labels.add(b["name"])
            return self.reply(201, b)
        if self.path == "/repos/o/r/issues":
            n = max([i["number"] for i in self.issues] + [0]) + 1
            self.issues.append({"number": n, "state": "open", "title": b["title"], "body": b["body"],
                                "labels": b.get("labels", []), "html_url": f"https://x/{n}"})
            if self.lose_create_answer:
                # The issue is created, but the answer says the server failed.
                return self.reply(502, {})
            return self.reply(201, self.issues[-1])
        self.reply(404, {})

    def do_PATCH(self):
        b = self.body()
        n = int(self.path.rsplit("/", 1)[1])
        for i in self.issues:
            if i["number"] == n:
                i.update(b)
                return self.reply(200, i)
        self.reply(404, {})


def self_test():
    fails = []

    def check(ok, what):
        print(("ok   " if ok else "FAIL ") + what)
        if not ok:
            fails.append(what)

    # 1. Version order: the sort -V traps.
    bio = ["1.3.6--h43da1c4_0", "1.3.10--h5ca1c30_0", "1.3.9--hdeadbe_2", "1.3.6--habc123_1", "latest"]
    check(newer_tags("1.3.6--h43da1c4_0", bio) == ("1.3.10--h5ca1c30_0", "1.3.10--h5ca1c30_0"),
          "biocontainers: 1.3.10 ranks above 1.3.9, whatever the build hash")
    check(change_kind("1.3.6--h43da1c4_0", "1.3.6--habc123_1") == "build", "biocontainers: a new build is a build")
    py = ["3.11", "3.9", "3.14", "3.14-slim", "3.15.0rc1", "3.11.9", "3.15-rc", "windowsservercore", "3.14.0"]
    check(newer_tags("3.11", py) == ("3.14", "3.14"), "python 3.11: only plain N.N tags count (3.14, not rc or -slim)")
    rl = Report()
    _orig = globals()["list_tags"]
    globals()["list_tags"] = lambda host, repo: py if repo == "library/python" else ["1.20", "1.24"]
    try:
        images(rl, {"PYTHON_IMAGE": "python:3.11", "SAMTOOLS_IMAGE": "staphb/samtools:1.20"}, {})
    finally:
        globals()["list_tags"] = _orig
    check(rl.line_rows == [("PYTHON_IMAGE", "3.11", "3.14")] and [r[0] for r in rl.image_rows] == ["SAMTOOLS_IMAGE"],
          "python:3.11 is a release-line pin, not behind; samtools 1.20 -> 1.24 is behind")
    vep = ["release_116.0", "release_116.2", "release_117.0", "latest", "116.2"]
    check(newer_tags("release_116.0", vep) == ("release_117.0", "release_116.2"), "VEP: major and same-major found")
    check(change_kind("release_116.0", "release_117.0") == "major", "VEP 116 -> 117 is a major")
    check(newer_tags("1.21", ["1.21", "1.20", "1.21-2"]) == (None, None), "current pin has nothing newer")
    try:
        newer_tags("latest", ["latest"])
        check(False, "a pin with no version is an error")
    except LookupFailed:
        check(True, "a pin with no version is an error")
    check(split_image("python:3.11")[:3] == ("docker.io", "library/python", "3.11"), "Docker Hub library image")
    check(split_image("quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0")[:3]
          == ("quay.io", "biocontainers/t1k", "1.0.9--h5ca1c30_0"), "quay image")
    check(split_image("a/b@sha256:00")[2:] == (None, "sha256:00"), "digest-only pin has no tag")
    check(split_image("python:3.11.17@sha256:ab") == ("docker.io", "library/python", "3.11.17", "sha256:ab"),
          "a tag pinned to a digest keeps its tag")
    rt = Report()
    _orig = globals()["list_tags"]
    globals()["list_tags"] = lambda host, repo: ["3.11.17", "3.11.18", "3.14.0"]
    try:
        images(rt, {"PYTHON_IMAGE": "python:3.11.17@sha256:ab"}, {})
    finally:
        globals()["list_tags"] = _orig
    check([(r[0], r[1], r[2]) for r in rt.image_rows] == [("PYTHON_IMAGE", "3.11.17", "3.14.0")],
          "a tag pinned to a digest is still checked for newer tags")
    asked = []

    def fake_digest(host, repo, ref):
        asked.append(ref)
        return {"latest": "sha256:" + "1" * 64, "3.11.17": "sha256:" + "2" * 64}.get(ref, "sha256:" + "3" * 64)
    _orig_d, _orig_t = globals()["manifest_digest"], globals()["list_tags"]
    globals()["manifest_digest"], globals()["list_tags"] = fake_digest, (lambda host, repo: ["latest"])
    try:
        rd = Report()
        digests(rd, {"PYTHON_IMAGE": "python:3.11.17@sha256:" + "2" * 64,
                     "OLD_IMAGE": "a/b@sha256:" + "1" * 64})
        rm = Report()
        digests(rm, {"PYTHON_IMAGE": "python:3.11.17@sha256:" + "4" * 64})
    finally:
        globals()["manifest_digest"], globals()["list_tags"] = _orig_d, _orig_t
    states = {r[0]: r[3] for r in rd.digest_rows}
    check("3.11.17" in asked and "still points at the pinned digest" in states.get("PYTHON_IMAGE", ""),
          "a tag pinned to a digest is compared with its own tag, not with latest")
    check(states.get("OLD_IMAGE") == "same as publisher's latest", "a digest-only pin is compared with latest")
    check(rm.digest_rows and "was rebuilt" in rm.digest_rows[0][3], "a rebuilt tag is reported")
    sent = []

    def fake_request(url, headers=None, **kw):
        sent.append(url)
        return {"tags": ["1.0"]}, {"Link": '<https://elsewhere.example/v2/x/tags/list?last=1.0>; rel="next"'}
    _orig = globals()["request"]
    globals()["request"] = fake_request
    try:
        list_tags("quay.io", "biocontainers/x")
        crossed = "followed"
    except LookupFailed:
        crossed = "refused"
    finally:
        globals()["request"] = _orig
    check(crossed == "refused" and not any("elsewhere.example" in u for u in sent),
          "a tag-list next page on another host is refused, and nothing is sent there")
    rs = Report()
    run_lookup(rs, "odd shape", lambda report: [].get("x"))
    check(rs.errors and rs.errors[0][0] == "odd shape" and "unexpected answer" in rs.errors[0][1]
          and rs.db_rows[0][3] == "ERROR", "an answer in an unexpected shape is an error row, not a crash")

    # 2. Fake server: an empty 200 must fail, and the issue logic.
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), FakeGitHub)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_port}"
    try:
        request(base + "/empty", retries=0)
        check(False, "an empty 200 answer fails closed")
    except LookupFailed:
        check(True, "an empty 200 answer fails closed")
    try:
        request(base + "/truncated", retries=0)
        check(False, "a truncated answer fails closed")
    except LookupFailed:
        check(True, "a truncated answer fails closed")
    except Exception as e:  # noqa: BLE001 - a raw protocol error is the bug this control catches
        check(False, f"a truncated answer fails closed (raised {type(e).__name__} instead)")
    try:
        request("https://192.0.2.1/", retries=0, timeout=5)
        check(False, "no answer (timeout) fails closed")
    except LookupFailed:
        check(True, "no answer (timeout) fails closed")

    saved = os.environ.get("GITHUB_TOKEN"), os.environ.get("FRESHNESS_GITHUB_API")
    os.environ["GITHUB_TOKEN"] = "fake"
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as fh:
        fh.write("report body\n")
    FakeGitHub.issues = [{"number": 7, "state": "open", "title": "other", "body": "image test failed",
                          "labels": [ISSUE_LABEL], "html_url": "https://x/7"}]
    a1 = update_issue("o/r", fh.name, api=base)
    a2 = update_issue("o/r", fh.name, api=base)
    mine = [i for i in FakeGitHub.issues if ISSUE_MARKER in i["body"]]
    check(a1[0] == "created" and a2 == ("edited", a1[1]) and len(mine) == 1,
          "first run creates the issue, second run edits the same one")
    check(FakeGitHub.issues[0]["body"] == "image test failed", "a freshness-labelled issue without the marker is untouched")
    check(ISSUE_LABEL in FakeGitHub.labels, "a missing label is created")
    mine[0]["state"] = "closed"
    a3 = update_issue("o/r", fh.name, api=base)
    check(a3 == ("reopened and edited", a1[1]) and mine[0]["state"] == "open",
          "a closed report is reopened, not duplicated")
    FakeGitHub.issues, FakeGitHub.lose_create_answer = [], True
    try:
        update_issue("o/r", fh.name, api=base)
        lost = "no error"
    except LookupFailed as e:
        lost = str(e)
    FakeGitHub.lose_create_answer = False
    created = [i for i in FakeGitHub.issues if ISSUE_MARKER in i["body"]]
    check("HTTP 502" in lost and len(created) == 1,
          f"a create whose answer is lost fails the run and is not sent twice ({len(created)} issue(s))")
    srv.shutdown()
    os.unlink(fh.name)
    for k, v in zip(("GITHUB_TOKEN", "FRESHNESS_GITHUB_API"), saved):
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v

    # 3. Controls on the real registries: a bogus image name, a pinned tag that
    #    does not exist and a bogus ClinVar URL must make the run exit non-zero,
    #    and the real image beside them must not be an error.
    with tempfile.TemporaryDirectory() as d:
        ve, se, out = os.path.join(d, "versions.env"), os.path.join(d, "setup.sh"), os.path.join(d, "out.md")
        with open(ve, "w") as fh:
            fh.write('SAMTOOLS_IMAGE="staphb/samtools:1.20"\n'
                     'BOGUS_IMAGE="quay.io/biocontainers/pgp-freshness-no-such-image:1.0--h0_0"\n'
                     'BADTAG_IMAGE="staphb/samtools:0.0.0-pgp-no-such-tag"\n')
        with open(se, "w") as fh:
            fh.write('CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/pgp-no-such-file.vcf.gz"\n')
        run = subprocess.run([sys.executable, __file__, "--versions", ve, "--setup", se,
                              "--sections", "clinvar,images", "--out", out],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        rc = run.returncode
        text = open(out).read() if os.path.exists(out) else ""
        print("\n".join("     | " + line for line in run.stdout.splitlines()))
        check(rc == 1, "the run with a bogus image and a bogus URL exits 1")
        errors = text.split("### Errors", 1)[-1] if "### Errors" in text else ""
        check("BOGUS_IMAGE" in errors, "the bogus image is an error")
        check("BADTAG_IMAGE" in errors and "pinned tag is not in" in errors, "a missing pinned tag is an error")
        check("ClinVar" in errors and "HTTP 404" in errors, "the bogus ClinVar URL is an error")
        check("SAMTOOLS_IMAGE" not in errors and "SAMTOOLS_IMAGE" in text, "the real image is looked up without error")

    print("SELF-TEST " + ("OK" if not fails else f"FAILED: {len(fails)} check(s)"))
    return 0 if not fails else 1


# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--versions", default=os.path.join(ROOT, "versions.env"))
    ap.add_argument("--setup", default=os.path.join(ROOT, "scripts", "setup.sh"))
    ap.add_argument("--sections", default=",".join(SECTIONS))
    ap.add_argument("--out")
    ap.add_argument("--update-issue", metavar="BODY_FILE")
    ap.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY"))
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)

    if a.self_test:
        return self_test()
    if a.update_issue:
        if not a.repo:
            ap.error("--update-issue needs --repo or GITHUB_REPOSITORY")
        try:
            update_issue(a.repo, a.update_issue)
        except LookupFailed as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1
        return 0

    sections = [s for s in a.sections.split(",") if s]
    bad = [s for s in sections if s not in SECTIONS]
    if bad or not sections:
        ap.error(f"unknown section(s) {bad}; choose from {', '.join(SECTIONS)}")
    pins, notes = read_pins(a.versions)
    report = build(pins, notes, a.versions, a.setup, sections)
    text = render(report, sections)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(text + "\n")
    print(text)
    if report.errors:
        print(f"ERROR: {len(report.errors)} lookup(s) failed.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
