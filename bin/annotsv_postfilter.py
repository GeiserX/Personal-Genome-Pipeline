#!/usr/bin/env python3
"""The checks a class 4 or 5 structural variant needs before anyone reads it.

AnnotSV gives a class to every SV of its input, and three kinds of call reach
class 4 or 5 without being a real rare deletion:

  - a call only one caller made;
  - a "deletion" whose two breakpoints sit on exon edges of one transcript:
    the reads of a processed transcript (a retrocopy, or cDNA) aligned across
    an intron;
  - an SV that many people carry.

Step 05 and the Nextflow ANNOTSV run this script on their class 4 and 5 table
(<sample>_sv_class4_5.tsv: AnnotSV's header and one "full" row per SV) and
write <sample>_sv_class4_5_checked.tsv, one row per SV with:

  Caller             the caller of the input record (INFO SVCALLER, which
                     step 05 writes: manta, delly, cnvpytor or sv_vcf)
  Exon_edge          for a deletion: yes when its start lies within 5 bp of
                     the end of one exon and its end within 5 bp of the
                     start of a later exon of the same transcript, else no;
                     "." for other types; "not checked" without the genes
                     file. Exon_edge_detail names the transcript and exons.
  Population_AF      AnnotSV's benign-set frequency for its type
                     (B_loss_AFmax for a deletion, B_gain_AFmax for a
                     duplication, B_ins_AFmax, B_inv_AFmax), empty when no
                     benign record covers it; Population_source names the
                     record. Common is yes at 1% or more, AnnotSV's own
                     -benignAF default.

--genes is AnnotSV's transcript file, Annotations_Human/Genes/GRCh38/
genes.RefSeq.sorted.bed of the annotation bundle (chrom, txStart, txEnd,
strand, gene, transcript, cdsStart, cdsEnd, exonStarts, exonEnds; starts
0-based, ends 1-based, as UCSC writes them).

How many callers saw each SV is a fact of step 22's consensus, which runs
after step 05 (and beside ANNOTSV in Nextflow). The reports join it when they
are written (bin/collect_summary.py, consensus_support below); --consensus
adds the same two columns here:

  Support            SURVIVOR's SUPP of the consensus record with the same
                     type whose two breakpoints lie within 1,000 bp of the
                     SV's (step 22's own distance), else 1: a call not in the
                     consensus was made by one caller
  Support_callers    the callers of that record (SUPP_VEC, named by the
                     consensus VCF's sample columns, which step 22 names
                     after the callers), else the SV's own caller

Usage: annotsv_postfilter.py --class45 CLASS45.tsv [--genes GENES.bed]
       [--consensus CONSENSUS.vcf.gz] --out CHECKED.tsv
"""
import argparse
import csv
import gzip
import os
import sys

EDGE_BP = 5
COMMON_AF = 0.01
CONSENSUS_BP = 1000
AF_COLUMNS = {"DEL": "B_loss", "DUP": "B_gain", "INS": "B_ins", "INV": "B_inv"}
KEEP = ["AnnotSV_ID", "ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type", "Gene_name", "ACMG_class",
        "AnnotSV_ranking_score"]
ADDED = ["Caller", "Exon_edge", "Exon_edge_detail", "Population_AF", "Population_source", "Common"]
SUPPORT = ["Support", "Support_callers"]


def open_text(path):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def chrom_key(c):
    """AnnotSV writes chromosomes without 'chr'; step 22 writes them with it."""
    return c[3:] if c.lower().startswith("chr") else c


def sv_type(row):
    t = (row.get("SV_type") or "").upper()
    return "DUP" if t.startswith("DUP") else t


def info_value(info, key):
    for kv in (info or "").split(";"):
        k, _, v = kv.partition("=")
        if k == key:
            return v
    return ""


def read_transcripts(path, chroms):
    """{chrom: [(start0, end1, transcript, gene, strand, [(exon_start0, exon_end1), ...])]}
    of AnnotSV's genes file, for the chromosomes in chroms only."""
    out = {}
    with open_text(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 10 or line.startswith("#") or chrom_key(c[0]) not in chroms:
                continue
            try:
                starts = [int(x) for x in c[8].split(",") if x]
                ends = [int(x) for x in c[9].split(",") if x]
                out.setdefault(chrom_key(c[0]), []).append(
                    (int(c[1]), int(c[2]), c[5], c[4], c[3], sorted(zip(starts, ends))))
            except ValueError:
                continue
    return out


def exon_edge(chrom, start, end, transcripts):
    """(flag, detail) for a deletion from start to end (AnnotSV's SV_start and
    SV_end: the base before the deleted bases, and the last deleted base).
    A deletion that removes an intron exactly starts on the last base of an
    exon and ends on the base before the next exon (that exon's 0-based start)."""
    for ts, te, tx, gene, strand, exons in transcripts.get(chrom_key(chrom), []):
        if te < start - EDGE_BP or ts > end + EDGE_BP:
            continue
        n = len(exons)
        for i, (_, e_end) in enumerate(exons):
            if abs(start - e_end) > EDGE_BP:
                continue
            for j in range(i + 1, n):
                if abs(end - exons[j][0]) <= EDGE_BP:
                    a, b = (i + 1, j + 1) if strand != "-" else (n - j, n - i)
                    return "yes", f"{tx} ({gene}): exon {a} and exon {b} edges"
    return "no", ""


def population(row):
    """(AF or None, source) of AnnotSV's benign set for the SV's type."""
    col = AF_COLUMNS.get(sv_type(row))
    if not col:
        return None, ""
    raw = (row.get(col + "_AFmax") or "").strip()
    try:
        af = float(raw)
    except ValueError:
        return None, ""
    return af, (row.get(col + "_source") or "").strip()


def check_rows(rows, transcripts=None):
    """The class 4 and 5 rows with the added columns. transcripts None: the
    genes file was not given, so the exon edges are not checked."""
    out = []
    for r in rows:
        o = {k: r.get(k, "") for k in KEEP}
        o["Caller"] = info_value(r.get("INFO"), "SVCALLER") or "."
        if sv_type(r) != "DEL":
            o["Exon_edge"], o["Exon_edge_detail"] = ".", ""
        elif transcripts is None:
            o["Exon_edge"], o["Exon_edge_detail"] = "not checked", ""
        else:
            try:
                o["Exon_edge"], o["Exon_edge_detail"] = exon_edge(r["SV_chrom"], int(r["SV_start"]), int(r["SV_end"]),
                                                                  transcripts)
            except (KeyError, ValueError):
                o["Exon_edge"], o["Exon_edge_detail"] = "not checked", ""
        af, src = population(r)
        o["Population_AF"] = "" if af is None else f"{af:g}"
        o["Population_source"] = src
        o["Common"] = "yes" if af is not None and af >= COMMON_AF else "no"
        out.append(o)
    return out


def read_consensus(path):
    """(caller names, records) of step 22's consensus VCF; a record is
    (chrom, pos, end, type, supp, supp_vec)."""
    names, recs = [], []
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names = c[9:]
                continue
            if len(c) < 8:
                continue
            t = info_value(c[7], "SVTYPE").upper()
            try:
                pos = int(c[1])
                end = int(info_value(c[7], "END") or pos)
                supp = int(info_value(c[7], "SUPP") or 0)
            except ValueError:
                continue
            recs.append((chrom_key(c[0]), pos, end, "DUP" if t.startswith("DUP") else t, supp,
                         info_value(c[7], "SUPP_VEC")))
    return names, recs


def consensus_support(rows, path):
    """Support and Support_callers for each checked row from the consensus
    at path; None for both when there is no consensus file."""
    if not path or not os.path.isfile(path):
        return [{"Support": None, "Support_callers": None} for _ in rows]
    names, recs = read_consensus(path)
    by = {}
    for rec in recs:
        by.setdefault((rec[0], rec[3]), []).append(rec)
    out = []
    for r in rows:
        best = None
        try:
            start, end = int(r["SV_start"]), int(r["SV_end"])
        except (KeyError, ValueError):
            start = end = None
        if start is not None:
            for rec in by.get((chrom_key(r.get("SV_chrom", "")), sv_type(r)), []):
                d = abs(rec[1] - start) + abs(rec[2] - end)
                if abs(rec[1] - start) <= CONSENSUS_BP and abs(rec[2] - end) <= CONSENSUS_BP and (best is None or d < best[0]):
                    best = (d, rec)
        if best:
            vec = best[1][5]
            callers = [names[i] if i < len(names) else f"caller {i + 1}" for i, x in enumerate(vec) if x == "1"]
            out.append({"Support": best[1][4], "Support_callers": ",".join(callers)})
        else:
            out.append({"Support": 1, "Support_callers": r.get("Caller") or "."})
    return out


def read_tsv(path):
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--class45", required=True, help="the class 4 and 5 table (AnnotSV header, full rows)")
    ap.add_argument("--genes", help="AnnotSV's genes.RefSeq.sorted.bed")
    ap.add_argument("--consensus", help="step 22's consensus VCF")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    rows = read_tsv(a.class45) if os.path.getsize(a.class45) else []
    if rows and not {"SV_chrom", "SV_start", "SV_end", "SV_type", "ACMG_class"} <= set(rows[0]):
        sys.exit(f"ERROR: {a.class45} lacks one of SV_chrom, SV_start, SV_end, SV_type, ACMG_class")
    transcripts = None
    if a.genes:
        if not os.path.isfile(a.genes):
            sys.exit(f"ERROR: no genes file {a.genes}")
        transcripts = read_transcripts(a.genes, {chrom_key(r.get("SV_chrom", "")) for r in rows})
    checked = check_rows(rows, transcripts)
    cols = KEEP + ADDED
    if a.consensus:
        for r, s in zip(checked, consensus_support(checked, a.consensus)):
            r.update({k: "" if v is None else v for k, v in s.items()})
        cols += SUPPORT
    tmp = a.out + ".part"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(checked)
    os.replace(tmp, a.out)
    edge = sum(r["Exon_edge"] == "yes" for r in checked)
    common = sum(r["Common"] == "yes" for r in checked)
    print(f"Class 4-5 SVs checked: {len(checked)}; deletions on exon edges: "
          f"{edge if transcripts is not None else 'not checked'}; common (benign-set AF >= 1%): {common}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
