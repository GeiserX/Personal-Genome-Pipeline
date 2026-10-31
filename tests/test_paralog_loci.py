#!/usr/bin/env python3
"""assets/paralog_loci.tsv against the models of assets/parascopy_models_v1.7.txt.

Step 35 runs `parascopy cn-using` once per row of the catalog, with the model
file the row names. The models are the genes of
assets/parascopy_models_v1.7.txt (every region of each model of Parascopy's
pinned data, models_GRCh38_1KGP_v1.7.tar.gz, Zenodo record 15019940), each as
<gene>.gz, the file name list.txt gives it (the e2e case paralog-catalog-1-gba
compares the two names lists for the five populations).
  1. the model list: 167 genes, each with exactly one main region, each name
     <gene>.gz;
  2. the catalog: its columns in order, and every row valid: a model on the
     list, a GRCh38 region, pscn_reliable yes or no, a known tag, a table
     name used once, a source; the SMN1 row writes the table `smn` the reports
     read;
  3. controls: a catalog naming a model not on the list (GBA1.gz, where the
     model file is GBA.gz), a bad region, a repeated table name and a catalog
     without the SMN1 row are each refused;
  4. the report readers: the not-assessed row step 35 writes (its own lines,
     cut from scripts/35-paralogs.sh and run) gives a reading that names the
     reason in bin/collect_summary.py, and in bin/paraphase_report.py and
     bin/couple_report.py (the Paraphase side-by-side and the couple screen),
     and a called table still reads as before. Control: the old row with "."
     for its place is not found over SMN1 exon 7 there.

Run: python3 tests/test_paralog_loci.py
"""
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODELS = os.path.join(REPO, "assets", "parascopy_models_v1.7.txt")
CATALOG = os.path.join(REPO, "assets", "paralog_loci.tsv")
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
    """One <gene>.gz per main line of the region file."""
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip() and not line.startswith("#")]
    return [r[0] + ".gz" for r in rows if len(r) == 3 and r[2] == "main"]


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
        for k, pat in (("model", r"[A-Za-z0-9._-]+\.gz"), ("region", r"chr[0-9XY]+:[1-9][0-9]*-[1-9][0-9]*"),
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
    with open(MODELS) as f:
        genes = {line.split("\t")[0] for line in f if line.strip() and not line.startswith("#")}
    check("the model list has 167 names", len(models) == 167, len(models))
    check("each model has exactly one main region", len(set(models)) == len(models) == len(genes),
          (len(set(models)), len(models), len(genes)))
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
        # the lines of step 35 that write a not-assessed table, run as they are
        with open(os.path.join(REPO, "scripts", "35-paralogs.sh")) as f:
            src = f.read()
        m = re.search(r"\n(  chrom=\$\{region%%:\*\}.*?\n  \} > \"\$out\"\n)", src, re.S)
        check("step 35 holds the lines that write a not-assessed table", m is not None)
        reason = "no read over chr5:70951946-70951946"
        snippet = ("header() { printf 'chrom\\tstart\\tend\\tlocus\\tagCN_filter\\tagCN\\tagCN_qual\\t"
                   "psCN_filter\\tpsCN\\tpsCN_qual\\thomologous_regions\\n'; }\n" + (m.group(1) if m else "false\n"))
        subprocess.run(["bash", "-c", snippet], check=False, env=dict(
            os.environ, region="chr5:70951946-70951946", model="SMN1.gz", reason=reason, out=t))
        import paraphase_report  # noqa: E402
        import couple_report  # noqa: E402
        _, data = collect_summary.sec_parascopy(d, "S")
        check("report: a not-assessed SMN table reads as not assessed, with the reason",
              data["reading"] == f"SMN1 and SMN2 not assessed: {reason}: no carrier reading"
              and data["smn1"]["cn"] is None, data.get("reading"))
        p = paraphase_report.parascopy_smn1(t)
        check("Paraphase side-by-side: a not-assessed SMN table gives no SMN1 copy number, with the reason",
              p["cn"] is None and reason in p["text"], p)
        found, txt = couple_report.smn_finding(t, {"SMN1"})
        check("couple screen: a not-assessed SMN table is a gap that names the reason, not a carrier",
              found == [] and reason in txt and "not resolved" in txt, txt)
        with open(t, "w") as f:
            f.write(head11 + f".\t.\t.\tSMN1\tnot assessed: {reason}\t.\t.\t.\t.\t.\t.\n")
        check("control: the old not-assessed row, with . for its place, loses the reason there",
              reason not in paraphase_report.parascopy_smn1(t)["text"])
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
