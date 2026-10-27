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

Some class 4 and 5 calls are common or unreliable for a known reason, and
the table moves them to class 3. ACMG_class is the class after these checks;
AnnotSV_class keeps AnnotSV's own. Rule names each check that fired and
Rule_detail says why:

  Depth_support      for a deletion or duplication: yes or no from the
                     record's own read depth, duphold's DHFFC (a deletion
                     needs < 0.7) or DHBFC (a duplication needs > 1.3), step
                     15's thresholds; else, for a Delly call over 1,000 bp,
                     its RDCN (a deletion needs 0 or 1 copies, a duplication
                     3 or more, or 2 or more on chrX and chrY, where the sex
                     is not known here); a CNVpytor call is a read-depth call,
                     so yes. A Delly call step 05 added (FORMAT GT only)
                     over 1,000 bp is yes: cnvpytor_to_vcf.py set aside the
                     ones whose RDCN does not show the change. "not checked"
                     when the record has no such value (step 15 did not run,
                     a Delly call of 1,000 bp or less); "." for other types.
                     "no" moves the call to class 3 (rule no_depth_support).
  --rules            assets/sv_rule_regions.tsv (scripts/ci/build-sv-rule-
                     regions.sh writes it), one region per row:
    tandem_repeat    an SV that lies inside a coding tandem repeat (each
                     breakpoint within 10 bp of it) is a change in the
                     repeat's length, which short-read SV callers do not size
                     and AnnotSV can read as a frameshift: class 3, and the
                     length is a question for a repeat genotyper;
    ig_tr_locus      an SV more than half inside IGH, IGK, IGL, TRA, TRB or
                     TRG: blood cells cut DNA out of these loci when they
                     rearrange their receptor genes (V(D)J recombination), so
                     a blood sample carries such deletions without any
                     inherited change: class 3;
    mhc              an SV more than half inside the MHC, whose haplotypes
                     differ by tens of kilobases between people; gnomAD-SV
                     marks its SVs low confidence (FILTER IGH_MHC_OVERLAP,
                     over 50% inside IGH or the MHC): class 3;
    gene_alias       a pseudogene whose HGNC alias or previous symbol is a
                     disease gene's symbol: when Gene_name lists that symbol
                     and the SV reaches no transcript of the disease gene in
                     the genes file (so --genes is needed), Gene_name names
                     the pseudogene instead, and the class goes to 3 when no
                     other gene is left.

Usage: annotsv_postfilter.py --class45 CLASS45.tsv [--genes GENES.bed]
       [--rules RULES.tsv] [--consensus CONSENSUS.vcf.gz] --out CHECKED.tsv
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
ADDED = ["Caller", "Exon_edge", "Exon_edge_detail", "Population_AF", "Population_source", "Common", "AnnotSV_class",
         "Depth_support", "Rule", "Rule_detail"]
DHFFC_MAX = 0.7
DHBFC_MIN = 1.3
DELLY_DEPTH_BP = 1000
REPEAT_SLACK_BP = 10
LOCUS_FRACTION = 0.5
RULES = {"tandem_repeat", "ig_tr_locus", "mhc", "gene_alias"}
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


def sample_format(row):
    """{FORMAT key: value} of the first sample column, the one after FORMAT
    in AnnotSV's header; {} when the row has none."""
    keys = list(row)
    if "FORMAT" not in keys:
        return {}
    i = keys.index("FORMAT") + 1
    if i >= len(keys) or keys[i] == "Annotation_mode":
        return {}
    return dict(zip((row.get("FORMAT") or "").split(":"), (row.get(keys[i]) or "").split(":")))


def number(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if x == x else None


def sv_span(row):
    """(start0, end1) of the SV: AnnotSV's SV_start is the base before a
    deletion, so the deleted bases are SV_start..SV_end-1 0-based. An
    insertion gets one base."""
    start, end = int(row["SV_start"]), int(row["SV_end"])
    return start, max(end, start + 1)


def depth_support(row, caller):
    """(yes|no|not checked|., detail) for a deletion or duplication."""
    t = sv_type(row)
    if t not in ("DEL", "DUP"):
        return ".", ""
    if caller == "cnvpytor":
        return "yes", "a read-depth call (CNVpytor)"
    fmt = sample_format(row)
    key, ok = ("DHFFC", lambda v: v < DHFFC_MAX) if t == "DEL" else ("DHBFC", lambda v: v > DHBFC_MIN)
    v = number(fmt.get(key))
    if v is not None:
        return ("yes" if ok(v) else "no"), f"duphold {key} {v:g}"
    rdcn = number(fmt.get("RDCN"))
    try:
        length = abs(int(row.get("SV_length") or 0)) or abs(int(row["SV_end"]) - int(row["SV_start"]))
    except (KeyError, ValueError):
        length = 0
    if rdcn is not None and rdcn >= 0 and length > DELLY_DEPTH_BP:
        sex_chrom = chrom_key(row.get("SV_chrom", "")).upper() in ("X", "Y")
        yes = rdcn <= 1 if t == "DEL" else rdcn >= (2 if sex_chrom else 3)
        return ("yes" if yes else "no"), f"Delly RDCN {rdcn:g}"
    if caller == "delly" and rdcn is None and length > DELLY_DEPTH_BP:
        # Step 05 adds a Delly call through cnvpytor_to_vcf.py, which keeps
        # FORMAT GT only, after it has set aside every call over 1,000 bp
        # whose RDCN does not show the change.
        return "yes", "RDCN checked when step 05 added the call (bin/cnvpytor_to_vcf.py)"
    return "not checked", "no duphold or Delly depth value in the record"


def read_rules(path, chroms):
    """{chrom: [(start0, end1, rule, name)]} of the rule table, for the
    chromosomes in chroms only."""
    out = {}
    with open_text(path) as f:
        for n, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6 or c[3] not in RULES:
                sys.exit(f"ERROR: {path} line {n}: want chrom, start, end, rule ({', '.join(sorted(RULES))}), "
                         f"action, name, source")
            if chrom_key(c[0]) in chroms:
                out.setdefault(chrom_key(c[0]), []).append((int(c[1]), int(c[2]), c[3], c[5]))
    return out


def region_rules(row, rules, transcripts):
    """([(rule, detail)], new Gene_name or None, downgrade) of the rule table
    for one row."""
    try:
        start, end = sv_span(row)
    except (KeyError, ValueError):
        return [], None, False
    chrom = chrom_key(row.get("SV_chrom", ""))
    genes = [g for g in (row.get("Gene_name") or "").split(";") if g]
    hits, seen, remap = [], set(), {}
    for rs, re_, rule, name in rules.get(chrom, []):
        if rs >= end + REPEAT_SLACK_BP or re_ <= start - REPEAT_SLACK_BP:
            continue
        where = f"{row.get('SV_chrom')}:{rs + 1}-{re_}"
        if rule == "tandem_repeat" and rule not in seen and start >= rs - REPEAT_SLACK_BP and end <= re_ + REPEAT_SLACK_BP:
            seen.add(rule)
            hits.append((rule, f"inside the coding tandem repeat {name} at {where} (TRExplorer): a change in the "
                               f"repeat's length, which SV callers do not size; read it with a repeat genotyper"))
        elif rule in ("ig_tr_locus", "mhc") and rule not in seen:
            frac = (min(end, re_) - max(start, rs)) / (end - start)
            if frac > LOCUS_FRACTION:
                seen.add(rule)
                if rule == "mhc":
                    hits.append((rule, f"{frac:.0%} inside the MHC ({where}), whose haplotypes differ widely between "
                                       f"people; gnomAD-SV marks such SVs low confidence (IGH_MHC_OVERLAP)"))
                else:
                    hits.append((rule, f"{frac:.0%} inside the {name} locus ({where}), where blood cells cut DNA out "
                                       f"when they rearrange their receptor genes (V(D)J recombination)"))
        elif rule == "gene_alias" and transcripts is not None and min(end, re_) > max(start, rs):
            alias, pseudo, pseudo_id, gene_id = (name.split("|") + ["", "", "", ""])[:4]
            if alias not in genes or alias in remap:
                continue
            real = any(g == alias and te > start and ts < end and not (ts >= rs and te <= re_)
                       for ts, te, _, g, _, _ in transcripts.get(chrom, []))
            if not real:
                remap[alias] = pseudo
                hits.append((rule, f"{alias} here is the pseudogene {pseudo} ({pseudo_id}), which carries {alias} as "
                                   f"an alias; the gene {alias} ({gene_id}) lies elsewhere"))
    new_genes = None
    downgrade = any(rule != "gene_alias" for rule, _ in hits)
    if remap:
        new_genes = ";".join(dict.fromkeys(remap.get(g, g) for g in genes))
        if all(g in remap for g in genes):
            downgrade = True
        else:
            hits = [(rule, d + ("; the class may rest on the other genes, so it stays" if rule == "gene_alias" else ""))
                    for rule, d in hits]
    return hits, new_genes, downgrade


def check_rows(rows, transcripts=None, rules=None):
    """The class 4 and 5 rows with the added columns. transcripts None: the
    genes file was not given, so the exon edges and gene aliases are not
    checked. rules: read_rules() of the rule table, or None."""
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
        o["AnnotSV_class"] = r.get("ACMG_class", "")
        o["Depth_support"], depth_detail = depth_support(r, o["Caller"])
        hits, new_genes, downgrade = region_rules(r, rules or {}, transcripts)
        if o["Depth_support"] == "no":
            hits.insert(0, ("no_depth_support", f"read depth does not show the change ({depth_detail})"))
            downgrade = True
        if new_genes is not None:
            o["Gene_name"] = new_genes
        cls = number(o["AnnotSV_class"])
        if downgrade and cls is not None and cls > 3:
            o["ACMG_class"] = "3"
        o["Rule"] = ";".join(h[0] for h in hits)
        o["Rule_detail"] = " | ".join(h[1] for h in hits)
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
    ap.add_argument("--rules", help="assets/sv_rule_regions.tsv")
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
    rules = None
    if a.rules:
        if not os.path.isfile(a.rules):
            sys.exit(f"ERROR: no rule table {a.rules}")
        rules = read_rules(a.rules, {chrom_key(r.get("SV_chrom", "")) for r in rows})
    checked = check_rows(rows, transcripts, rules)
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
          f"{edge if transcripts is not None else 'not checked'}; common (benign-set AF >= 1%): {common}; "
          f"moved to class 3: {sum(r['ACMG_class'] != r['AnnotSV_class'] for r in checked)}"
          f"{'' if rules is not None else ' (no rule table)'}; depth support no: "
          f"{sum(r['Depth_support'] == 'no' for r in checked)}, not checked: "
          f"{sum(r['Depth_support'] == 'not checked' for r in checked)}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
