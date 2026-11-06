#!/usr/bin/env python3
"""Check the 1.0 output contract, assets/output_contract.tsv (docs/output-contract.md).

Two checks:

  --static [--repo DIR]
      The command surface and the Nextflow file list against the code, with no
      run: every param of nextflow_schema.json, every --tools name of main.nf
      (known_tools) and every run-all.sh switch (its SWITCHES list) has a row,
      and every row still names one; every file row of a process names a file
      that process still writes and publishes; and every file a module
      publishes has a row. A stable row that no longer matches is a break. A
      new param, tool, switch or published file with no row fails too, so the
      contract stays the full list.

  --outputs DIR --leg bash|nextflow --sample S [--observe]
      The files a run wrote. DIR is the folder that holds <sample>/ (GENOME_DIR
      for the bash steps, --outdir for Nextflow). For every stable row whose e2e
      column names this leg: the file must be there under its contracted name,
      and its required columns (tsv, csv), INFO/FORMAT fields (vcf) or top-level
      keys (json) must be in it. Experimental rows are listed, never enforced.
      --observe also prints every file under DIR/<sample> with what it has, the
      lines to read when the contract is written again.

Exit 0 when nothing is broken, 1 when something is, 2 on a usage error. Only
the standard library: it runs on the CI host and in self-tests.
"""
import argparse
import fnmatch
import glob
import gzip
import json
import os
import re
import sys

COLUMNS = ["kind", "name", "step", "process", "legs", "format", "stability", "e2e", "required"]
KINDS = {"file", "param", "tool", "switch"}
STABILITY = {"stable", "experimental", "internal"}
LEGS = {"bash", "nextflow", "both"}
FORMATS = {"tsv", "csv", "vcf", "json", "text", "html", "bam", "cram", "index", "bed", "dir",
           "zip", "svg", "png", "md", "yml", "fasta", "bin"}


def read_contract(path):
    """The rows of the contract as dicts; raises ValueError on a malformed file."""
    rows, header = [], None
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cells = line.split("\t")
            if header is None:
                if cells != COLUMNS:
                    raise ValueError(f"{path}:{n}: header must be {' '.join(COLUMNS)}, got {' '.join(cells)}")
                header = cells
                continue
            if len(cells) != len(COLUMNS):
                raise ValueError(f"{path}:{n}: {len(cells)} cells, want {len(COLUMNS)}")
            r = dict(zip(COLUMNS, cells))
            r["line"] = n
            bad = []
            if r["kind"] not in KINDS:
                bad.append(f"kind {r['kind']}")
            if r["stability"] not in STABILITY:
                bad.append(f"stability {r['stability']}")
            if r["kind"] == "file":
                if r["legs"] not in LEGS:
                    bad.append(f"legs {r['legs']}")
                if r["format"] not in FORMATS:
                    bad.append(f"format {r['format']}")
                if r["e2e"] not in ("-", "bash", "nextflow", "both"):
                    bad.append(f"e2e {r['e2e']}")
                if r["legs"] == "bash" and r["process"] != "-":
                    bad.append("a bash-only row names a process")
            if bad:
                raise ValueError(f"{path}:{n}: bad " + ", ".join(bad))
            rows.append(r)
    if header is None:
        raise ValueError(f"{path}: no header line")
    seen = {}
    for r in rows:
        key = (r["kind"], r["name"], r["legs"] if r["kind"] == "file" else "")
        if key in seen:
            raise ValueError(f"{path}:{r['line']}: {r['kind']} {r['name']} is listed twice (line {seen[key]})")
        seen[key] = r["line"]
    return rows


# ---------------------------------------------------------------- the code side

def braces(pattern):
    """Expand one level of {a,b} alternatives, as Nextflow globs do."""
    m = re.search(r"\{([^{}]*)\}", pattern)
    if not m:
        return [pattern]
    out = []
    for alt in m.group(1).split(","):
        out.extend(braces(pattern[:m.start()] + alt + pattern[m.end():]))
    return out


def concrete(pattern, sample="SAMPLE"):
    """A Nextflow path with ${meta.id} and the family name filled in."""
    p = pattern.replace("${meta.id}", sample).replace("{sample}", sample)
    p = re.sub(r"\$\{params\.family_name[^}]*\}", "FAMILY", p).replace("{family}", "FAMILY")
    p = p.replace("${params.outdir}/", "").replace("${params.outdir}", "")
    return p.strip("/")


def module_processes(repo):
    """{process: {'publish': [(dir, [patterns] or None, has_saveas)], 'outputs': [glob]}} of modules/local."""
    procs = {}
    for f in sorted(glob.glob(os.path.join(repo, "modules", "local", "*", "main.nf"))):
        src = open(f).read()
        starts = [(m.start(), m.group(1)) for m in re.finditer(r"^process\s+(\w+)\s*\{", src, re.M)]
        for i, (pos, name) in enumerate(starts):
            body = src[pos:starts[i + 1][0] if i + 1 < len(starts) else len(src)]
            head = re.split(r"^\s*(?:script|shell|exec|stub):", body, maxsplit=1, flags=re.M)[0]
            pubs = []
            lines = head.split("\n")
            for j, line in enumerate(lines):
                if not re.match(r"\s*publishDir\b", line):
                    continue
                stmt = line
                k = j
                while stmt.rstrip().endswith(",") or stmt.count("{") > stmt.count("}"):
                    k += 1
                    if k >= len(lines):
                        break
                    stmt += "\n" + lines[k]
                d = re.search(r'publishDir\s+\{?\s*"([^"]+)"', stmt)
                pat = re.search(r"pattern:\s*[\"']([^\"']+)[\"']", stmt)
                pubs.append((concrete(d.group(1)) if d else "",
                             braces(pat.group(1)) if pat else None,
                             "saveAs" in stmt))
            out = re.search(r"^\s*output:\s*$(.*?)(?:^\s*(?:when|topic):|\Z)", head, re.M | re.S)
            outputs = []
            if out:
                for a, b in re.findall(r"path\s*\(?\s*\"([^\"]+)\"|path\s*\(?\s*'([^']+)'", out.group(1)):
                    g = a or b
                    if g != "versions.yml":
                        outputs.append(g)
            procs[name] = {"publish": pubs, "outputs": outputs, "file": os.path.relpath(f, repo)}
    return procs


def published(proc):
    """(contracted-path glob, basename-only) pairs a process publishes; basename-only when a saveAs renames."""
    items = []
    for d, pats, saveas in proc["publish"]:
        for g in proc["outputs"]:
            c = concrete(g)
            if pats is not None and not any(fnmatch.fnmatch(c.replace("*", "X"), p) or fnmatch.fnmatch(c, p)
                                            for p in pats):
                continue
            items.append((g, (d + "/" + c) if d else c, saveas))
    return items


def schema_params(repo):
    s = json.load(open(os.path.join(repo, "nextflow_schema.json")))
    defs = s.get("definitions") or s.get("$defs") or {}
    return {p for group in defs.values() for p in group.get("properties", {})}


def known_tools(repo):
    src = open(os.path.join(repo, "main.nf")).read()
    m = re.search(r"def known_tools = \[(.*?)\]", src, re.S)
    return set(re.findall(r"'([a-z0-9_]+)'", m.group(1))) if m else set()


def run_all_switches(repo):
    src = open(os.path.join(repo, "scripts", "run-all.sh")).read()
    m = re.search(r"^SWITCHES=\((.*?)\)", src, re.S | re.M)
    return set(m.group(1).split()) if m else set()


def check_static(repo, rows):
    problems = []
    have = {"param": schema_params(repo), "tool": known_tools(repo), "switch": run_all_switches(repo)}
    where = {"param": "nextflow_schema.json", "tool": "main.nf known_tools", "switch": "scripts/run-all.sh SWITCHES"}
    for kind in ("param", "tool", "switch"):
        if not have[kind]:
            problems.append(f"found no {kind}s in {where[kind]}: the reader no longer matches the code")
            continue
        listed = {r["name"]: r for r in rows if r["kind"] == kind}
        for name, r in sorted(listed.items()):
            if name not in have[kind]:
                what = "BREAKING: the stable" if r["stability"] == "stable" else "the"
                problems.append(f"{what} {kind} {name} (contract line {r['line']}) is gone from {where[kind]}")
        for name in sorted(have[kind] - set(listed)):
            problems.append(f"{kind} {name} of {where[kind]} has no row in the contract: add one")

    procs = module_processes(repo)
    if not procs:
        problems.append("found no processes under modules/local")
    file_rows = [r for r in rows if r["kind"] == "file" and r["process"] != "-"]
    for r in file_rows:
        if r["process"] not in procs:
            what = "BREAKING: " if r["stability"] == "stable" else ""
            problems.append(f"{what}{r['name']} (contract line {r['line']}): process {r['process']} is not in modules/local")
            continue
        p = procs[r["process"]]
        want = concrete(r["name"])
        if r["stability"] == "internal":
            ok = any(fnmatch.fnmatch(os.path.basename(want), os.path.basename(concrete(g))) for g in p["outputs"])
        else:
            ok = any(fnmatch.fnmatch(want, path) if not saveas
                     else fnmatch.fnmatch(os.path.basename(want), os.path.basename(path))
                     for _, path, saveas in published(p))
        if not ok:
            what = "BREAKING: the stable file" if r["stability"] == "stable" else "the file"
            problems.append(f"{what} {r['name']} (contract line {r['line']}) is not what {r['process']} "
                            f"publishes ({p['file']}): renamed, moved or no longer written")
    for name, p in sorted(procs.items()):
        rows_p = [concrete(r["name"]) for r in file_rows if r["process"] == name]
        if not p["publish"]:
            continue
        for g in p["outputs"]:
            c = os.path.basename(concrete(g))
            if not any(fnmatch.fnmatch(os.path.basename(w), c) for w in rows_p):
                problems.append(f"{name} ({p['file']}) writes {g} and has no contract row for it: "
                                "add one (internal when it is not published)")
    return problems


# ---------------------------------------------------------------- the run side

def opener(path):
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def fields(path, fmt):
    """What a file has: its columns, its INFO/FORMAT ids or its top-level keys."""
    if fmt in ("tsv", "csv"):
        sep = "\t" if fmt == "tsv" else ","
        with opener(path) as f:
            for i, line in enumerate(f):
                if i > 200:
                    break
                line = line.rstrip("\n")
                if line.startswith("##") or sep not in line:
                    continue
                return [c.strip() for c in line.lstrip("#").split(sep)]
        return []
    if fmt == "vcf":
        out = []
        with opener(path) as f:
            for line in f:
                if not line.startswith("##"):
                    break
                m = re.match(r"##(INFO|FORMAT|FILTER)=<ID=([^,>]+)", line)
                if m:
                    out.append(f"{m.group(1)}/{m.group(2)}")
        return out
    if fmt == "json":
        try:
            data = json.load(opener(path))
        except ValueError:
            return None
        return list(data.keys()) if isinstance(data, dict) else []
    return []


def guess_format(path):
    b = os.path.basename(path)
    for ext, fmt in ((".vcf.gz", "vcf"), (".vcf", "vcf"), (".tsv.gz", "tsv"), (".tsv", "tsv"), (".csv", "csv"),
                     (".json", "json"), (".txt", "tsv")):
        if b.endswith(ext):
            return fmt
    return "-"


def check_outputs(rows, outdir, leg, sample, observe):
    problems, checked = [], 0
    for r in rows:
        if r["kind"] != "file" or r["stability"] != "stable" or r["e2e"] not in (leg, "both"):
            continue
        if r["legs"] not in (leg, "both"):
            problems.append(f"contract line {r['line']}: e2e {r['e2e']} on a {r['legs']} row")
            continue
        checked += 1
        name = r["name"].replace("{sample}", sample)
        hits = sorted(glob.glob(os.path.join(outdir, name)))
        if not hits:
            problems.append(f"BREAKING: the stable file {name} (contract line {r['line']}, step {r['step']}) "
                            f"is missing from {outdir}: renamed, moved or not written")
            continue
        req = [] if r["required"] == "-" else r["required"].split(",")
        if r["format"] == "dir":
            if not all(os.path.isdir(h) for h in hits):
                problems.append(f"{name} (contract line {r['line']}) is not a folder")
            continue
        for h in hits:
            if os.path.isdir(h):
                problems.append(f"{os.path.relpath(h, outdir)} (contract line {r['line']}) is a folder, want a file")
                continue
            if not req:
                continue
            have = fields(h, r["format"])
            if have is None:
                problems.append(f"{os.path.relpath(h, outdir)} (contract line {r['line']}) is not valid {r['format']}")
                continue
            for col in req:
                if col not in have:
                    problems.append(f"BREAKING: {os.path.relpath(h, outdir)} (contract line {r['line']}) "
                                    f"has no {col} (it has: {','.join(have) or 'nothing'})")
    if observe:
        base = os.path.join(outdir, sample)
        print(f"# observed under {base} (leg {leg}): path, format, fields")
        for root, dirs, files in os.walk(base):
            dirs[:] = sorted(d for d in dirs if d not in ("work", "nf-work", ".nextflow"))
            depth = os.path.relpath(root, base).count(os.sep)
            for fn in sorted(files):
                p = os.path.join(root, fn)
                rel = os.path.relpath(p, outdir).replace(sample, "{sample}")
                fmt = guess_format(p)
                try:
                    have = fields(p, fmt) if fmt != "-" and os.path.getsize(p) < 2 << 30 else []
                except (OSError, EOFError, UnicodeDecodeError):
                    have = []
                print(f"OBSERVED\t{leg}\t{rel}\t{fmt}\t{','.join(have or [])}"[:4000])
            if depth >= 3:
                dirs[:] = []
    return checked, problems


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--contract", help="default: <repo>/assets/output_contract.tsv")
    ap.add_argument("--repo", default=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    ap.add_argument("--static", action="store_true")
    ap.add_argument("--outputs")
    ap.add_argument("--leg", choices=["bash", "nextflow"])
    ap.add_argument("--sample")
    ap.add_argument("--observe", action="store_true")
    a = ap.parse_args(argv)
    if not a.static and not a.outputs:
        ap.error("give --static or --outputs")
    if a.outputs and not (a.leg and a.sample):
        ap.error("--outputs needs --leg and --sample")
    try:
        rows = read_contract(a.contract or os.path.join(a.repo, "assets", "output_contract.tsv"))
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}")
        return 1
    problems = []
    if a.static:
        p = check_static(a.repo, rows)
        print(f"static: {len(rows)} contract rows, {len(p)} problems")
        problems += p
    if a.outputs:
        n, p = check_outputs(rows, a.outputs, a.leg, a.sample, a.observe)
        print(f"outputs ({a.leg}): {n} stable files checked under {a.outputs}, {len(p)} problems")
        if n == 0:
            p.append(f"no stable row has e2e {a.leg}: the check checked nothing")
        problems += p
    for p in problems:
        print(f"FAIL {p}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
