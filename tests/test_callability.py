#!/usr/bin/env python3
"""The gene callability of step 16b and the ALT contig refusal.

assets/callability_genes.tsv lists the genes short reads often cannot place
(scripts/ci/build-callability.sh writes it; the e2e case
tests/e2e/callability-1-catalog.sh rebuilds it and compares). Step 16b counts,
per gene, the bases with 10 or more reads of MAPQ 20 or more
(<sample>_callability.bed.gz, mosdepth's thresholds table), and
bin/collect_summary.py turns that into the report's section.
  1. the catalog: its header, one row per gene, known reasons, an HGNC id per
     row, a position on chr1-22, X or Y unless ALT-only, and owners as the
     build script's rule gives them (assets/paralog_loci.tsv,
     assets/pgx_genes.tsv, bin/cn_loci.tsv);
  2. a gene whose reads all had MAPQ 0 (no base at 10x of MAPQ >= 20) is not
     read; a unique-sequence gene with every base at 10x is read; an ALT-only
     gene is never read, whatever the table holds. Control: a reader that
     takes an ALT-only gene as read (the planted bug) fails this check;
  3. the count: the report's list of genes not read and with no owner has as
     many genes as a query of the catalog gives (rows with owner '.', minus
     the ones the table reads), in the summary and in both reports. Control:
     a reader that drops the ALT-only genes from that list fails it;
  4. a table from another catalog (a gene missing) makes the section
     unreadable, not a reading;
  5. bin/copy_number.py check-bam and call --bam-header refuse a header with
     an ALT or an HLA contig (exit 4, a REFUSED line) and pass a header
     without one (the positive control).

Run: python3 tests/test_callability.py
"""
import csv
import gzip
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402
import render_report  # noqa: E402

CATALOG = os.path.join(REPO, "assets", "callability_genes.tsv")
COLUMNS = ["gene", "hgnc_id", "chrom", "start", "end", "reasons", "owner"]
REASONS = {"CMRG", "segdup", "false_duplication", "collapsed", "ALT-only"}
CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {str(detail)[:400]}"))
    if not ok:
        FAILS += 1


def rows_of(path):
    with open(path) as f:
        return list(csv.DictReader((x for x in f if not x.startswith("#")), delimiter="\t", quoting=csv.QUOTE_NONE))


def table_of(path):
    with open(path) as f:
        return list(csv.DictReader((x for x in f if x.strip() and not x.startswith("#")), delimiter="\t",
                                   quoting=csv.QUOTE_NONE))


def write_thresholds(path, catalog, zero=(), drop=()):
    """mosdepth's thresholds table over the catalog's genes: every base at
    10x, except the genes in zero (none) and the ones in drop (no row)."""
    with gzip.open(path, "wt") as f:
        f.write("#chrom\tstart\tend\tregion\t10X\n")
        for r in catalog:
            if r["start"] == "." or r["gene"] in drop:
                continue
            n = int(r["end"]) - int(r["start"])
            f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['gene']}\t{0 if r['gene'] in zero else n}\n")


# --- 1. the catalog --------------------------------------------------------------------
with open(CATALOG) as f:
    header = next(x for x in f if not x.startswith("#")).rstrip("\n").split("\t")
check("1. the catalog's columns", header == COLUMNS, header)
cat = rows_of(CATALOG)
genes = [r["gene"] for r in cat]
check("1. one row per gene", len(genes) == len(set(genes)) and len(genes) > 1000, len(genes))
check("1. every row has an HGNC id", all(r["hgnc_id"].startswith("HGNC:") for r in cat))
check("1. every reason is known", all(set(r["reasons"].split(",")) <= REASONS for r in cat))
bad = [r["gene"] for r in cat if ("ALT-only" in r["reasons"]) != (r["start"] == ".")
       or (r["start"] == "." and not r["chrom"].endswith("_alt"))
       or (r["start"] != "." and (r["chrom"] not in CHROMS or not int(r["start"]) < int(r["end"])))]
check("1. a position on chr1-22, X or Y, or an ALT scaffold with no position for an ALT-only gene", not bad, bad)
check("1. 273 CMRG genes", sum("CMRG" in r["reasons"] for r in cat) == 273)
by = {r["gene"]: r for r in cat}
check("1. GSTT1 is ALT-only and pypgx's; HLA-DRB3 is ALT-only with no owner",
      by.get("GSTT1", {}).get("reasons") == "ALT-only" and by["GSTT1"]["owner"] == "pypgx"
      and by.get("HLA-DRB3", {}).get("reasons") == "ALT-only" and by["HLA-DRB3"]["owner"] == ".")
check("1. KCNE1, CBS and CRYAA are in a false duplication of GRCh38",
      all("false_duplication" in by.get(g, {}).get("reasons", "") for g in ("KCNE1", "CBS", "CRYAA")))
parascopy = {g for r in table_of(os.path.join(REPO, "assets", "paralog_loci.tsv")) for g in r["genes"].split(",")}
pypgx = {r["gene"] for r in table_of(os.path.join(REPO, "assets", "pgx_genes.tsv")) if r["mode"] == "bam"}
wrong = [r["gene"] for r in cat if (r["gene"] in parascopy and r["owner"] != "parascopy")
         or (r["gene"] not in parascopy and r["gene"] in pypgx and r["owner"] != "pypgx")
         or (r["owner"] not in (".", "parascopy", "pypgx") and not r["owner"].startswith("cn_loci:"))]
check("1. owners: parascopy for assets/paralog_loci.tsv's genes, then pypgx for the bam genes", not wrong, wrong)
check("1. HBA1 and HBA2 belong to the depth engine's HBA row",
      by.get("HBA1", {}).get("owner") == "cn_loci:HBA" and by.get("HBA2", {}).get("owner") == "cn_loci:HBA")


# --- 2 to 4. the reader ----------------------------------------------------------------
def alt_read_bug(catalog, table):
    """The planted bug: an ALT-only gene taken as read."""
    d = collect_summary.callability(catalog, table)
    d["not_read"] = [g for g in d["not_read"] if g["status"] != "not on this reference"]
    d["read"] += len(d["not_on_reference"])
    return d


def readings(d):
    return {g["gene"]: g["status"] for g in d["not_read"] + d["not_read_owned"]}


def mapq0_and_alt(reader, tsv):
    """Check 2 against reader: the MAPQ 0 gene and the ALT-only genes not
    read, the unique gene read. [] when it holds, else what is wrong."""
    d = reader(cat, tsv)
    got = readings(d)
    out = []
    if got.get(MAPQ0) != "not read":
        out.append(f"{MAPQ0}: {got.get(MAPQ0)}")
    if UNIQUE in got:
        out.append(f"{UNIQUE}: {got[UNIQUE]}")
    for g in ("GSTT1", "HLA-DRB3"):
        if got.get(g) != "not on this reference":
            out.append(f"{g}: {got.get(g, 'read')}")
    return out


def count_matches(reader, tsv):
    """Check 3: the unowned genes not read are the catalog query's."""
    d = reader(cat, tsv)
    query = [r["gene"] for r in cat if r["owner"] == "." and (r["start"] == "." or r["gene"] == MAPQ0)]
    return [g["gene"] for g in d["not_read"]] == query, d


MAPQ0 = "KCNE1"   # all its reads MAPQ 0, as the false duplication of GRCh38 makes them
UNIQUE = "CBS"    # every base at 10x
with tempfile.TemporaryDirectory() as work:
    tsv = os.path.join(work, "S_callability.bed.gz")
    write_thresholds(tsv, cat, zero={MAPQ0})
    check("2. control: a reader that takes an ALT-only gene as read fails check 2",
          bool(mapq0_and_alt(alt_read_bug, tsv)), "the planted bug passed")
    problems = mapq0_and_alt(collect_summary.callability, tsv)
    check("2. a gene of MAPQ 0 reads is not read, a unique one is read, an ALT-only one is never read",
          not problems, problems)
    check("3. control: the planted bug's list differs from the catalog query",
          not count_matches(alt_read_bug, tsv)[0])
    ok, d = count_matches(collect_summary.callability, tsv)
    check("3. the genes not read with no owner are the catalog query's", ok, d["not_read"][:5])
    check("3. owned genes not read are listed apart (GSTT1 by pypgx)",
          any(g["gene"] == "GSTT1" and g["owner"] == "pypgx" for g in d["not_read_owned"]), d["not_read_owned"])

    # the section and both reports, from a sample folder (the Nextflow layout)
    sd = os.path.join(work, "S")
    os.makedirs(os.path.join(sd, "coverage"))
    os.replace(tsv, os.path.join(sd, "coverage", "S_callability.bed.gz"))
    s = collect_summary.collect("S", sd)
    sec = s["sections"]["callability"]
    check("3. the section reads the Nextflow folder (coverage/)", sec["state"] == "ok", sec)
    n = len(d["not_read"])
    txt, html = render_report.text_report(s), render_report.html_report(s)
    line = f"Not read, and no other step calls them ({n}): "
    check("3. the text report lists them with the count", line in txt and f"{MAPQ0} (0%)" in txt
          and "HLA-DRB3 (not on this reference)" in txt, txt[:600])
    check("3. the HTML report lists them with the count", render_report.E(line) in html
          and "Genes This Genome Could Not Read Reliably" in html)

    # 4. a table written from another catalog
    os.makedirs(os.path.join(sd, "mosdepth"))
    write_thresholds(os.path.join(sd, "mosdepth", "S_callability.bed.gz"), cat, drop={UNIQUE})
    sec = collect_summary.collect("S", sd)["sections"]["callability"]
    check("4. a table missing a catalog gene is unreadable", sec["state"] == "unreadable"
          and UNIQUE in (sec["note"] or ""), sec)


# --- 5. the ALT contig refusal ---------------------------------------------------------
def copy_number(*args, stdin=None):
    p = subprocess.run([sys.executable, os.path.join(REPO, "bin", "copy_number.py"), *args],
                       capture_output=True, text=True, input=stdin)
    return p.returncode, p.stderr


with tempfile.TemporaryDirectory() as work:
    noalt = os.path.join(work, "noalt.txt")
    with open(noalt, "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr16\tLN:90338345\n@SQ\tSN:chrM\tLN:16569\n")
    alt = os.path.join(work, "alt.txt")
    with open(alt, "w") as f:
        f.write("@SQ\tSN:chr16\tLN:90338345\n@SQ\tSN:chr22_KI270879v1_alt\tLN:304135\n@RG\tID:x\tSM:y\n")
    rc, err = copy_number("check-bam", "--header", noalt)
    check("5. positive control: a header without ALT contigs passes", rc == 0 and not err, (rc, err))
    rc, err = copy_number("check-bam", "--header", alt)
    check("5. a header with an ALT contig is refused (exit 4, a REFUSED line naming it)",
          rc == 4 and err.startswith("REFUSED:") and "chr22_KI270879v1_alt" in err, (rc, err))
    rc, err = copy_number("check-bam", "--header", "-", stdin="chr1\nHLA-A*01:01:01:01\n")
    check("5. contig names on stdin with an HLA contig are refused", rc == 4 and "HLA-A" in err, (rc, err))
    bed = os.path.join(work, "bedcov.tsv")
    open(bed, "w").close()
    rc, err = copy_number("call", "--locus", "HBA", "--bedcov", bed, "--bedcov-mapq1", bed, "--bam-header", alt)
    check("5. call --bam-header refuses before it reads any depth", rc == 4 and err.startswith("REFUSED:"),
          (rc, err))
    rc, err = copy_number("call", "--locus", "HBA", "--bedcov", bed, "--bedcov-mapq1", bed, "--bam-header", noalt)
    check("5. call --bam-header without ALT contigs goes on (two empty tables: not assessed)", rc == 0, (rc, err))

print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
sys.exit(1 if FAILS else 0)
