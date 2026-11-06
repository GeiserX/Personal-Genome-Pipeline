#!/usr/bin/env python3
"""assets/output_contract.tsv and scripts/ci/check-output-contract.py.

  1. the static check on this checkout passes: every param, --tools name,
     run-all.sh switch and published Nextflow file has its row, and every row
     still matches the code;
  2. static controls, each on a copy of the code: a renamed published file, a
     removed param, a removed switch, a new --tools name with no row, a new
     published file with no row and a schema default that is not
     nextflow.config's are each refused, and the stable ones are named
     BREAKING; so are, on a copy of the contract, a file row in the folder of
     a publishDir whose saveAs drops it and a nextflow row whose step script
     writes the same path (also on a line of seven shell variables);
     --pedigree is a declared param; a saveAs closure is read up to its own
     brace and a script path counts only as a whole path;
  3. the run check on a planted results folder made from the contract's own
     e2e rows passes for both legs; then a dropped column, a dropped VCF
     field, a dropped JSON key and a renamed file are each refused, naming the
     file and what is missing; an experimental row is never enforced;
  4. a contract with a row listed twice or a bad stability is refused.

Run: python3 tests/test_output_contract.py
"""
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHECK = os.path.join(REPO, "scripts", "ci", "check-output-contract.py")
CONTRACT = os.path.join(REPO, "assets", "output_contract.tsv")
SAMPLE = "S1"
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {detail}"))
    if not ok:
        FAILS += 1


def run(*args):
    p = subprocess.run([sys.executable, CHECK, *args], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def rows():
    out, header = [], None
    for line in open(CONTRACT):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        cells = line.split("\t")
        if header is None:
            header = cells
            continue
        out.append(dict(zip(header, cells)))
    return out


def code_copy(tmp):
    """The files the static check reads, copied, so a control can change them."""
    dst = os.path.join(tmp, "repo")
    shutil.copytree(os.path.join(REPO, "modules"), os.path.join(dst, "modules"))
    os.makedirs(os.path.join(dst, "scripts"))
    for f in ("main.nf", "nextflow_schema.json", "nextflow.config"):
        shutil.copy(os.path.join(REPO, f), os.path.join(dst, f))
    for f in os.listdir(os.path.join(REPO, "scripts")):
        if f.endswith(".sh"):
            shutil.copy(os.path.join(REPO, "scripts", f), os.path.join(dst, "scripts", f))
    return dst


def contract_copy(tmp, name, change):
    """The contract with one row changed: change(cells) edits the cells of the file row called name."""
    path = os.path.join(tmp, "contract.tsv")
    hit = False
    with open(CONTRACT) as f, open(path, "w") as g:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[0] == "file" and c[1] == name:
                change(c)
                line = "\t".join(c) + "\n"
                hit = True
            g.write(line)
    assert hit, f"no file row {name}"
    return path


def edit(path, old, new, count=1):
    src = open(path).read()
    assert old in src, f"{old!r} not in {path}"
    open(path, "w").write(src.replace(old, new, count))


def static_controls():
    rc, out = run("--static")
    check("static check passes on this checkout", rc == 0 and "0 problems" in out, out[-2000:])

    cases = [
        ("a renamed published file (ROH's summary)",
         lambda d: edit(os.path.join(d, "modules", "local", "roh", "main.nf"),
                        'path("*_roh_summary.txt")', 'path("*_roh_segments.txt")'),
         [r"BREAKING: the stable file \{sample\}/roh/\{sample\}_roh_summary\.txt", r"writes \*_roh_segments\.txt and has no contract row"]),
        ("a removed param (--intervals)",
         lambda d: _drop_param(d, "intervals"),
         [r"BREAKING: the stable param intervals"]),
        ("a removed run-all.sh switch (SKIP_TRIM)",
         lambda d: edit(os.path.join(d, "scripts", "run-all.sh"), " SKIP_TRIM ", " "),
         [r"BREAKING: the stable switch SKIP_TRIM"]),
        ("a new --tools name with no row",
         lambda d: edit(os.path.join(d, "main.nf"), "'pharmcat', 'cpic',", "'pharmcat', 'newtool', 'cpic',"),
         [r"tool newtool of main.nf known_tools has no row"]),
        ("a new published file with no row",
         lambda d: edit(os.path.join(d, "modules", "local", "roh", "main.nf"),
                        'path("*_roh.txt"),', 'path("*_roh.txt"),\n    tuple val(meta), path("*_roh_extra.tsv"), emit: extra\n'),
         [r"ROH .* writes \*_roh_extra\.tsv and has no contract row"]),
        ("a schema default that is not nextflow.config's (traits left out of --tools)",
         lambda d: edit(os.path.join(d, "nextflow_schema.json"), "mito_variants,traits,html_report", "mito_variants,html_report"),
         [r"param tools: nextflow_schema\.json's default .* is not nextflow\.config's"]),
    ]
    for desc, change, want in cases:
        tmp = tempfile.mkdtemp()
        try:
            d = code_copy(tmp)
            change(d)
            rc, out = run("--static", "--repo", d, "--contract", CONTRACT)
            missing = [w for w in want if not re.search(w, out)]
            check(f"static control refused: {desc}", rc == 1 and not missing, f"rc={rc}, missing {missing}: {out[-1500:]}")
        finally:
            shutil.rmtree(tmp)


    def to_aligned(c):
        c[1] = "{sample}/aligned/{sample}_read_count.tsv"

    def to_nextflow(c):
        c[4] = "nextflow"

    contract_cases = [
        ("the read count in aligned/, whose publishDir drops it (it is published to qc/)",
         "{sample}/qc/{sample}_read_count.tsv", to_aligned,
         [r"BREAKING: the stable file \{sample\}/aligned/\{sample\}_read_count\.tsv .* is not what ALIGN_MINIMAP2 publishes"]),
        ("legs nextflow on a file step 49 writes under the same name",
         "{sample}/chip/{sample}_chip.tsv", to_nextflow,
         [r"\{sample\}/chip/\{sample\}_chip\.tsv .* has legs nextflow, but scripts/49-chip\.sh writes it"]),
        ("legs nextflow on a file step 50 writes on a line of seven shell variables",
         "{sample}/mei/{sample}_megane.log", to_nextflow,
         [r"\{sample\}/mei/\{sample\}_megane\.log .* has legs nextflow, but scripts/50-mei\.sh writes it"]),
    ]
    for desc, name, change, want in contract_cases:
        tmp = tempfile.mkdtemp()
        try:
            rc, out = run("--static", "--contract", contract_copy(tmp, name, change))
            missing = [w for w in want if not re.search(w, out)]
            check(f"static control refused: {desc}", rc == 1 and not missing, f"rc={rc}, missing {missing}: {out[-1500:]}")
        finally:
            shutil.rmtree(tmp)

    config = open(os.path.join(REPO, "nextflow.config")).read()
    schema = json.load(open(os.path.join(REPO, "nextflow_schema.json")))
    in_schema = any("pedigree" in g.get("properties", {}) for g in schema["definitions"].values())
    check("--pedigree is declared in nextflow.config and nextflow_schema.json",
          re.search(r"^\s*pedigree\s*=\s*null\b", config, re.M) is not None and in_schema)


def _drop_param(d, name):
    path = os.path.join(d, "nextflow_schema.json")
    s = json.load(open(path))
    for group in s["definitions"].values():
        group.get("properties", {}).pop(name, None)
    json.dump(s, open(path, "w"))


def plant(root, leg):
    """A results folder with every stable e2e file of a leg, each carrying its required fields."""
    for r in rows():
        if r["kind"] != "file" or r["stability"] != "stable" or r["e2e"] not in (leg, "both"):
            continue
        path = os.path.join(root, r["name"].replace("{sample}", SAMPLE).replace("*", "x"))
        req = [] if r["required"] == "-" else r["required"].split(",")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if r["format"] == "dir":
            os.makedirs(path, exist_ok=True)
            continue
        if r["format"] == "vcf":
            text = "##fileformat=VCFv4.2\n" + "".join(
                f"##{f.split('/')[0]}=<ID={f.split('/')[1]},Number=1,Type=String,Description=\"x\">\n" for f in req) + \
                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n"
        elif r["format"] == "json":
            text = json.dumps({k: None for k in req})
        elif r["format"] in ("tsv", "csv") and req:
            sep = "\t" if r["format"] == "tsv" else ","
            text = "# a comment line\n" + sep.join(req) + "\n" + sep.join("1" for _ in req) + "\n"
        else:
            text = "x\n"
        if path.endswith(".gz"):
            with gzip.open(path, "wt") as f:
                f.write(text)
        else:
            open(path, "w").write(text)


def run_controls():
    legs = {"bash": None, "nextflow": None}
    tmp = tempfile.mkdtemp()
    try:
        for leg in legs:
            root = os.path.join(tmp, leg)
            plant(root, leg)
            rc, out = run("--outputs", root, "--leg", leg, "--sample", SAMPLE)
            n = re.search(r"(\d+) stable files checked", out)
            check(f"run check passes on a planted {leg} folder", rc == 0 and n and int(n.group(1)) >= 5, out[-1500:])
            legs[leg] = root

        nf = legs["nextflow"]
        hits = os.path.join(nf, SAMPLE, "clinvar", f"{SAMPLE}_clinvar_hits.tsv")
        lines = open(hits).read().split("\n")
        cols = lines[1].split("\t")
        i = cols.index("clnsig")
        open(hits, "w").write("\n".join("\t".join(c for j, c in enumerate(l.split("\t")) if j != i) for l in lines))
        rc, out = run("--outputs", nf, "--leg", "nextflow", "--sample", SAMPLE)
        check("a dropped column (clnsig of the ClinVar hits) is refused and named",
              rc == 1 and re.search(r"BREAKING: S1/clinvar/S1_clinvar_hits\.tsv .* has no clnsig", out), out[-1500:])

        roh = os.path.join(nf, SAMPLE, "roh", f"{SAMPLE}_roh_summary.txt")
        os.rename(roh, roh.replace("_roh_summary.txt", "_roh_segments.txt"))
        rc, out = run("--outputs", nf, "--leg", "nextflow", "--sample", SAMPLE)
        check("a renamed file (ROH summary) is refused as missing",
              rc == 1 and re.search(r"BREAKING: the stable file S1/roh/S1_roh_summary\.txt .* is missing", out), out[-1500:])

        ann = os.path.join(nf, SAMPLE, "vep", f"{SAMPLE}_annotated.vcf.gz")
        with gzip.open(ann, "rt") as f:
            text = f.read()
        with gzip.open(ann, "wt") as f:
            f.write(re.sub(r"##INFO=<ID=REVEL,[^\n]*\n", "", text))
        rc, out = run("--outputs", nf, "--leg", "nextflow", "--sample", SAMPLE)
        check("a dropped VCF field (INFO/REVEL) is refused and named",
              rc == 1 and re.search(r"S1_annotated\.vcf\.gz .* has no INFO/REVEL", out), out[-1500:])

        summ = os.path.join(legs["bash"], SAMPLE, "summary.json")
        data = json.load(open(summ))
        data.pop("not_assessed")
        json.dump(data, open(summ, "w"))
        rc, out = run("--outputs", legs["bash"], "--leg", "bash", "--sample", SAMPLE)
        check("a dropped JSON key (not_assessed of summary.json) is refused and named",
              rc == 1 and re.search(r"summary\.json .* has no not_assessed", out), out[-1500:])

        # An experimental row with e2e set is still not enforced.
        exp = os.path.join(tmp, "contract.tsv")
        with open(CONTRACT) as f, open(exp, "w") as g:
            for line in f:
                if line.startswith("file\t{sample}/bloodtype/{sample}_bloodtype.tsv\t"):
                    c = line.rstrip("\n").split("\t")
                    c[7] = "both"
                    line = "\t".join(c) + "\n"
                g.write(line)
        fresh = os.path.join(tmp, "fresh")
        plant(fresh, "bash")
        rc, out = run("--outputs", fresh, "--leg", "bash", "--sample", SAMPLE, "--contract", exp)
        check("an experimental file is listed but not enforced", rc == 0, out[-1500:])
        rc, out = run("--outputs", os.path.join(tmp, "empty"), "--leg", "bash", "--sample", SAMPLE)
        check("an empty results folder fails", rc == 1 and "is missing" in out, out[-500:])
    finally:
        shutil.rmtree(tmp)


def reader_controls():
    """The two readers the static check is built on, on inputs the code does not have yet."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("contract_check", CHECK)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    stmt = "publishDir \"x\", saveAs: { fn -> fn.endsWith('.log') ? null : fn }, enabled: { true }"
    body = m.closure_body(stmt)
    check("a saveAs closure ends at its own brace (an option after it is not part of it)",
          body is not None and "enabled" not in body and m.save_as(body, "a.log") == "drop", repr(body))
    tmp = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(tmp, "scripts"))
        open(os.path.join(tmp, "scripts", "11-roh.sh"), "w").write(
            'S="${GENOME_DIR}/${SAMPLE}"\nbgzip -c x > "${S}/roh/${SAMPLE}_roh.txt.gz"\n')
        short = m.bash_writes(tmp, "11", "{sample}/roh/{sample}_roh.txt", {})
        full = m.bash_writes(tmp, "11", "{sample}/roh/{sample}_roh.txt.gz", {})
        check("a step script writing x.txt.gz does not count as writing x.txt",
              short is None and full == "scripts/11-roh.sh", f"{short!r} {full!r}")
    finally:
        shutil.rmtree(tmp)


def malformed():
    tmp = tempfile.mkdtemp()
    try:
        lines = open(CONTRACT).read().split("\n")
        first = next(i for i, l in enumerate(lines) if l.startswith("file\t"))
        dup = os.path.join(tmp, "dup.tsv")
        open(dup, "w").write("\n".join(lines[:first + 1] + [lines[first]] + lines[first + 1:]))
        rc, out = run("--static", "--contract", dup)
        check("a row listed twice is refused", rc == 1 and "listed twice" in out, out[-500:])
        bad = os.path.join(tmp, "bad.tsv")
        open(bad, "w").write("\n".join(lines[:first] + [lines[first].replace("\tstable\t", "\tsolid\t")] + lines[first + 1:]))
        rc, out = run("--static", "--contract", bad)
        check("an unknown stability is refused", rc == 1 and "bad stability solid" in out, out[-500:])
    finally:
        shutil.rmtree(tmp)


static_controls()
reader_controls()
run_controls()
malformed()
print(f"{FAILS} failed")
sys.exit(1 if FAILS else 0)
