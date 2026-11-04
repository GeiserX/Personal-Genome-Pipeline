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
     model file is GBA.gz; CYP21A2.gz, which the data has no model for), a bad
     region, a region outside its model's main region, a repeated table name
     and a catalog without the SMN1 row are each refused;
  5. one owner per locus: no catalog model's regions (main or copy) overlap
     the targets or controls of a bin/cn_loci.tsv row the depth engine owns
     (owner cn-depth), except RHCE, whose model covers RHD as a cross-check
     of the RHD row (the exemption is needed: without it RHCE is refused);
     every catalog gene assets/callability_genes.tsv lists has the owner
     parascopy there (a planted '.' owner is refused);
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
TAGS = {"carrier", "SF", "disease", "PGx", "trait"}
CROSS_CHECKS = {"paraphase", "cn-depth", "none"}
CN_LOCI = os.path.join(REPO, "bin", "cn_loci.tsv")
CALLABILITY = os.path.join(REPO, "assets", "callability_genes.tsv")
# The one catalog model allowed to overlap a cn-depth row of bin/cn_loci.tsv:
# RHCE's model covers most of RHD, and the catalog reads it as a cross-check.
CN_DEPTH_EXEMPT = {("RHCE.gz", "RHD")}
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


def span(text):
    """'chr1:10-20' -> ('chr1', 10, 20), 1-based and inclusive; None if it is not one."""
    m = re.fullmatch(r"(chr[0-9XY]+):(\d+)-(\d+)", text)
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None


def model_regions(path):
    """<gene>.gz -> [(chrom, start, end, kind)] of every region of the model, 1-based."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                g, r, kind = line.rstrip("\n").split("\t")
                m = re.fullmatch(r"([^:]+):(\d+)-(\d+)", r)
                out.setdefault(g + ".gz", []).append((m.group(1), int(m.group(2)), int(m.group(3)), kind))
    return out


def table_rows(path):
    with open(path) as f:
        lines = [x.rstrip("\n").split("\t") for x in f if x.strip() and not x.startswith("#")]
    return [dict(zip(lines[0], r)) for r in lines[1:]]


def cn_depth_parts(path):
    """[(locus, chrom, start, end)] of every target and control of the cn_loci rows
    the depth engine owns, 1-based inclusive (the file is BED: 0-based, end exclusive)."""
    out = []
    for r in table_rows(path):
        if r["owner"] != "cn-depth":
            continue
        for part in (r["targets"] + ";" + r["controls"]).split(";"):
            c = part.split(":")
            if len(c) >= 2 and "-" in c[1]:
                s, e = c[1].split("-")
                out.append((r["locus"], r["chrom"], int(s) + 1, int(e)))
    return out


def owner_problems(catalog_models, regions, cn_parts, exempt):
    """Every overlap of a catalog model's region with a cn-depth row's target or
    control that the exemptions do not name."""
    out = []
    for model in catalog_models:
        for c, s, e, kind in regions.get(model, []):
            for locus, cc, cs, ce in cn_parts:
                if c == cc and s <= ce and cs <= e and (model, locus) not in exempt:
                    out.append(f"{model} ({kind} {c}:{s}-{e}) overlaps cn_loci row {locus} ({cc}:{cs}-{ce})")
    return out


def callability_problems(catalog_genes, rows):
    """Every catalog gene the callability catalog lists without the owner parascopy."""
    return [f"{r['gene']} has owner {r['owner']}" for r in rows
            if r["gene"] in catalog_genes and r["owner"] != "parascopy"]


def problems(lines, models, regions=None):
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
        elif regions is not None and c["model"] in regions:
            r = span(c["region"])
            if not any(k == "main" and ch == r[0] and s <= r[1] and r[2] <= e for ch, s, e, k in regions[c["model"]]):
                out.append(f"row {n}: region {c['region']} is not inside the main region of {c['model']}")
        if c["cross_check"] not in CROSS_CHECKS:
            out.append(f"row {n}: cross_check {c['cross_check']}")
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
        for k in ("genes", "reading", "source"):
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
    regions = model_regions(MODELS)

    # 2. the catalog
    with open(CATALOG) as f:
        lines = f.readlines()
    found = problems(lines, models, regions)
    check("assets/paralog_loci.tsv: every row is valid", not found, "; ".join(found))
    rows = table_rows(CATALOG)
    check("the catalog has the SMN1 row and 36 more", len(rows) == 37, len(rows))

    # 3. controls: each planted fault is refused
    head = "\t".join(COLUMNS) + "\n"
    smn = "SMN1.gz\tSMN1,SMN2\tchr5:70951946-70951946\tx\tyes\tcarrier\tparaphase\tsmn\ts\n"
    check("control: the SMN1 row alone is valid", not problems([head, smn], models), problems([head, smn], models))
    planted = {
        "a model not on the list (GBA1.gz)": "GBA1.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba\ts\n",
        "a model the data does not have (CYP21A2.gz)":
            "CYP21A2.gz\tCYP21A2\tchr6:32038265-32041670\tx\tno\tcarrier\tnone\tcyp21a2\ts\n",
        "a region outside its model's main region": "GBA.gz\tGBA1\tchr5:70951946-70951946\tx\tno\tcarrier\tnone\tgba\ts\n",
        "an unknown cross_check": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tgauchian\tgba\ts\n",
        "a region on no GRCh38 contig": "GBA.gz\tGBA1\tchr1_alt:1-2\tx\tno\tcarrier\tnone\tgba\ts\n",
        "a region that ends before it starts": "GBA.gz\tGBA1\tchr1:20-10\tx\tno\tcarrier\tnone\tgba\ts\n",
        "a repeated table name": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tsmn\ts\n",
        "a table name with a space": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tg a\ts\n",
        "a table name the step refuses (gba:2)": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba:2\ts\n",
        "a table name that is a path (..)": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\t..\ts\n",
        "an unknown tag": "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcancer\tnone\tgba\ts\n",
    }
    for desc, row in planted.items():
        check(f"control: {desc} is refused", bool(problems([head, smn, row], models, regions)))
    check("control: a catalog without the SMN1 row is refused",
          bool(problems([head, "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba\ts\n"], models)))

    # 5. one owner per locus
    cn_parts = cn_depth_parts(CN_LOCI)
    check("bin/cn_loci.tsv has cn-depth rows to compare with (HBA, RHD among them)",
          {"HBA", "RHD"} <= {x[0] for x in cn_parts}, sorted({x[0] for x in cn_parts}))
    catalog_models = [r["model"] for r in rows]
    found = owner_problems(catalog_models, regions, cn_parts, CN_DEPTH_EXEMPT)
    check("no catalog model overlaps a cn-depth row of bin/cn_loci.tsv but the RHCE cross-check", not found,
          "; ".join(found[:5]))
    check("control: without the exemption RHCE's model is refused (it covers RHD)",
          any("RHCE.gz" in x and "RHD" in x for x in owner_problems(catalog_models, regions, cn_parts, set())))
    planted_regions = dict(regions, **{"HBAX.gz": [("chr16", 172000, 178000, "main")]})
    check("control: a model over HBA1 and HBA2 is refused",
          bool(owner_problems(["HBAX.gz"], planted_regions, cn_parts, CN_DEPTH_EXEMPT)))
    rhce = next((r for r in rows if r["model"] == "RHCE.gz"), {})
    check("the RHCE row says it cross-checks the depth engine, and does not list RHD",
          rhce.get("cross_check") == "cn-depth" and "RHD" not in rhce.get("genes", "").split(","), rhce)
    catalog_genes = {g for r in rows for g in r["genes"].split(",")}
    cal = table_rows(CALLABILITY)
    found = callability_problems(catalog_genes, cal)
    check("every catalog gene assets/callability_genes.tsv lists has the owner parascopy", not found, "; ".join(found))
    check("the callability catalog lists catalog genes (GBA1, C4A, AMY1C, SMN1 among them)",
          {"GBA1", "C4A", "AMY1C", "SMN1"} <= {r["gene"] for r in cal if r["gene"] in catalog_genes})
    planted_cal = [dict(r, owner=".") if r["gene"] == "C4A" else r for r in cal]
    check("control: a catalog gene without the owner parascopy is refused",
          bool(callability_problems(catalog_genes, planted_cal)))

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
