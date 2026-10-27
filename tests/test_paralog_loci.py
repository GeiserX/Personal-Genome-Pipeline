#!/usr/bin/env python3
"""assets/paralog_loci.tsv against assets/parascopy_models_v1.7.txt.

Step 35 runs `parascopy cn-using` once per row of the catalog, with the model
file the row names. The model list is list.txt of Parascopy's pinned data
(models_GRCh38_1KGP_v1.7.tar.gz, Zenodo record 15019940, archive sha256
2e9be08b69936f13a49d52c7764f2414beac27a99b95e0c48b2e2ccafbb8c846; the five
populations' lists are the same file).
  1. the model list: 167 names, each once, each <name>.gz, and the sha256 of
     the file as committed;
  2. the catalog: its columns in order, and every row valid: a model on the
     list, a GRCh38 region, pscn_reliable yes or no, a known tag, a table
     name used once, a source; the SMN1 row writes the table `smn` the reports
     read;
  3. controls: a catalog naming a model not on the list (GBA1.gz, where the
     model file is GBA.gz), a bad region, a repeated table name and a catalog
     without the SMN1 row are each refused;
  4. the report reader (bin/collect_summary.py): a not-assessed SMN table
     from step 35 gives a reading that names the reason, and a called one
     still reads as before.

Run: python3 tests/test_paralog_loci.py
"""
import hashlib
import os
import re
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODELS = os.path.join(REPO, "assets", "parascopy_models_v1.7.txt")
CATALOG = os.path.join(REPO, "assets", "paralog_loci.tsv")
MODELS_SHA256 = "4ef9d7a4f3c9749281348509de8f6cf695b896073836944f5305a665e0a179a2"
COLUMNS = ["model", "genes", "region", "reading", "pscn_reliable", "tag", "cross_check", "table", "source"]
TAGS = {"carrier", "SF", "PGx", "trait"}
CONTIGS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {detail}"))
    if not ok:
        FAILS += 1


def read_models(path):
    with open(path) as f:
        return [line.rstrip("\n") for line in f if line.strip()]


def problems(lines, models):
    """Every problem of a catalog (its lines, comments included), as text."""
    out = []
    rows = [line.rstrip("\n").split("\t") for line in lines if line.strip() and not line.startswith("#")]
    if not rows or rows[0] != COLUMNS:
        return [f"header is {rows[0] if rows else None}, want {COLUMNS}"]
    tables = set()
    for n, r in enumerate(rows[1:], 2):
        if len(r) != len(COLUMNS):
            out.append(f"row {n}: {len(r)} columns")
            continue
        c = dict(zip(COLUMNS, r))
        if c["model"] not in models:
            out.append(f"row {n}: model {c['model']} is not in the pinned model list")
        m = re.fullmatch(r"(chr[0-9XY]+):(\d+)-(\d+)", c["region"])
        if not m or m.group(1) not in CONTIGS or int(m.group(2)) < 1 or int(m.group(2)) > int(m.group(3)):
            out.append(f"row {n}: region {c['region']} is not chrom:start-end on GRCh38")
        if c["pscn_reliable"] not in ("yes", "no"):
            out.append(f"row {n}: pscn_reliable {c['pscn_reliable']}")
        if c["tag"] not in TAGS:
            out.append(f"row {n}: tag {c['tag']}")
        # model, region and table reach a shell loop and file names: the
        # patterns of scripts/35-paralogs.sh (a table is never "." or "..")
        for k, pat in (("model", r"[A-Za-z0-9._-]+\.gz"), ("region", r"chr[0-9XY]+:[0-9]+-[0-9]+"),
                       ("table", r"[A-Za-z0-9][A-Za-z0-9_-]*")):
            if not re.fullmatch(pat, c[k]):
                out.append(f"row {n}: {k} {c[k]!r} is not what the step accepts")
        if c["table"] in tables or c["table"] in ("hba", "paraphase"):
            out.append(f"row {n}: table {c['table']} is used twice")
        tables.add(c["table"])
        for k in ("genes", "reading", "cross_check", "source"):
            if not c[k].strip():
                out.append(f"row {n}: empty {k}")
    if not any(r[:1] == ["SMN1.gz"] and len(r) == len(COLUMNS) and r[7] == "smn" for r in rows[1:]):
        out.append("no SMN1.gz row writing the table smn, which the reports read")
    return out


def main():
    # 1. the model list
    models = read_models(MODELS)
    with open(MODELS, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    check("the model list is the committed list.txt (sha256)", sha == MODELS_SHA256, sha)
    check("the model list has 167 names", len(models) == 167, len(models))
    check("each model is named once", len(set(models)) == len(models))
    check("each name is <model>.gz", all(re.fullmatch(r"[A-Za-z0-9._-]+\.gz", m) for m in models),
          [m for m in models if not m.endswith(".gz")])
    check("SMN1.gz and GBA.gz are on the list", {"SMN1.gz", "GBA.gz"} <= set(models))
    models = set(models)

    # 2. the catalog
    with open(CATALOG) as f:
        lines = f.readlines()
    found = problems(lines, models)
    check("assets/paralog_loci.tsv: every row is valid", not found, "; ".join(found))

    # 3. controls: each planted fault is refused
    head = "\t".join(COLUMNS) + "\n"
    smn = "SMN1.gz\tSMN1,SMN2\tchr5:70951946-70951946\tx\tyes\tcarrier\tparaphase\tsmn\ts\n"
    check("control: the SMN1 row alone is valid", not problems([head, smn], models), problems([head, smn], models))
    planted = {
        "a model not on the list (GBA1.gz)": "GBA1.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba\ts\n",
        "a region on no GRCh38 contig": "GBA.gz\tGBA1\tchr1_alt:1-2\tx\tno\tcarrier\tnone\tgba\ts\n",
        "a region that ends before it starts": "GBA.gz\tGBA1\tchr1:20-10\tx\tno\tcarrier\tnone\tgba\ts\n",
        "a repeated table name": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tsmn\ts\n",
        "a table name with a space": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tg a\ts\n",
        "a table name the step refuses (gba:2)": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba:2\ts\n",
        "a table name that is a path (..)": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\t..\ts\n",
        "an unknown tag": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcancer\tnone\tgba\ts\n",
    }
    for desc, row in planted.items():
        check(f"control: {desc} is refused", bool(problems([head, smn, row], models)))
    check("control: a catalog without the SMN1 row is refused",
          bool(problems([head, "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba\ts\n"], models)))

    # 4. the report reader on step 35's tables
    sys.path.insert(0, os.path.join(REPO, "bin"))
    import collect_summary  # noqa: E402
    head11 = ("chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\t"
              "homologous_regions\n")
    with tempfile.TemporaryDirectory() as d:
        os.makedirs(os.path.join(d, "paralogs"))
        t = os.path.join(d, "paralogs", "S_smn_copy_number.tsv")
        with open(t, "w") as f:
            f.write(head11 + ".\t.\t.\tSMN1\tnot assessed: no read over chr5:70951946-70951946\t.\t.\t.\t.\t.\t.\n")
        _, data = collect_summary.sec_parascopy(d, "S")
        check("report: a not-assessed SMN table reads as not assessed, with the reason",
              data["reading"] == "SMN1 and SMN2 not assessed: no read over chr5:70951946-70951946: no carrier reading"
              and data["smn1"]["cn"] is None, data.get("reading"))
        with open(t, "w") as f:
            f.write(head11 + "chr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tPASS\t1,3\t35,30\t"
                    "chr5:70049000-70077000:-\n")
        _, data = collect_summary.sec_parascopy(d, "S")
        check("report: a called SMN table still reads one SMN1 copy as a carrier",
              data["smn1"]["cn"] == 1 and "carrier" in data["reading"] and "not assessed" not in data["reading"],
              data.get("reading"))

    if FAILS:
        print(f"{FAILS} check(s) failed")
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
