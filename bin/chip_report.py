#!/usr/bin/env python3
"""chip_report.py: clonal haematopoiesis (CHIP) calls of step 49, with the detection limit at the sample's depth.

  chip_report.py prepare --gff ENSEMBL.gff3.gz --genes assets/chip_genes.tsv --fai REF.fai \\
      --out-gff chip.gff3 --out-bed chip.bed
  chip_report.py report --genes assets/chip_genes.tsv --gff chip.gff3 --bed chip.bed \\
      --calls CALLS.tsv --depth DEPTH.tsv --sample S --out S_chip.tsv --depth-out S_chip_depth.tsv
  chip_report.py check-table assets/chip_genes.tsv

prepare picks, for each gene of the table, its MANE Select transcript in
Ensembl's GFF3 (the transcript tagged MANE_Select) and writes the gene, that
transcript and its exons, CDS and UTRs with UCSC contig names (1 -> chr1,
MT -> chrM), the file bcftools csq reads. It also writes the coding exons, 10
bases wider on each side so the splice sites are inside, as a BED of the
contigs the reference has: Mutect2 calls inside it, and the depth is counted
over it. A table gene with no MANE Select transcript stops the run.

report reads two files bcftools writes (scripts/49-chip.sh and the Nextflow
CHIP_ANNOTATE process run the same commands):
  CALLS.tsv  one line per Mutect2 record after FilterMutectCalls, split to one
             ALT per line and annotated by bcftools csq: CHROM POS REF ALT
             FILTER BCSQ AD AF DP
  DEPTH.tsv  bcftools mpileup over the BED, reads of mapping quality 20 and
             bases of quality 10 and up, duplicates left out (Mutect2's own
             defaults): CHROM POS DP. A base with no line has depth 0.

Each call that changes a protein of a table gene on its MANE transcript gets
one row with its class and a status:
  chip              PASS, in the table's classes, VAF under GERMLINE_VAF and
                    at least MIN_ALT_READS reads with the change
  likely_germline   VAF GERMLINE_VAF or more: an inherited variant sits near
                    50% or 100%, and a tumour-only call cannot tell the two
                    apart
  filtered          FilterMutectCalls did not pass it (the filter names are in
                    the row)
  low_support       fewer than MIN_ALT_READS reads with the change
  other_change      a protein change outside the table's classes
Synonymous, intronic and UTR changes are left out.

The depth table has one row per gene: its coding bases on the reference, the
mean depth over them and the detection limit at that depth, the smallest VAF
a clone needs to be found with probability DETECT_POWER, taken as at least
MIN_ALT_READS reads with the change out of the mean depth (binomial). Below it
a clone is more likely missed than found, so a negative result never rules
one out.

Standard library only (it runs in the plain python image).
"""
import argparse
import gzip
import math
import os
import re
import sys

TABLE_COLUMNS = ["gene", "truncating", "missense", "note", "pmids"]
OUT_COLUMNS = ["gene", "chrom", "pos", "ref", "alt", "consequence", "protein_change", "exon", "vaf",
               "alt_reads", "depth", "filter", "status", "class", "note"]
DEPTH_COLUMNS = ["gene", "transcript", "chrom", "coding_bases", "mean_depth", "detection_limit_vaf", "status"]
GERMLINE_VAF = 0.35   # at or above: read as a possible inherited variant, not counted
MIN_ALT_READS = 3     # reads with the change a call needs, and the detection limit assumes
DETECT_POWER = 0.8    # the detection limit: the VAF found this often at the mean depth
CHIP_MIN_VAF = 0.02   # CHIP is defined from a VAF of 2% (Steensma et al. 2015, PMID 25931582)
PAD = 10              # bases each side of a coding exon: the splice donor and acceptor sites
TRUNCATING = {"frameshift", "stop_gained", "splice_donor", "splice_acceptor"}
MISSENSE = {"missense"}
PROTEIN = TRUNCATING | MISSENSE | {"inframe_deletion", "inframe_insertion", "inframe_altering", "start_lost",
                                   "stop_lost", "feature_elongation"}
AA = set("ACDEFGHIKLMNPQRSTVWY")


def die(msg):
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# ─── the gene table ─────────────────────────────────────────────────────────
def read_table(path):
    rows = []
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if not ln.startswith("#") and ln.strip()]
    if not lines or lines[0].split("\t") != TABLE_COLUMNS:
        die(f"{path}: the header must be {' '.join(TABLE_COLUMNS)} (tab separated)")
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) != len(TABLE_COLUMNS):
            die(f"{path}: '{ln[:60]}' has {len(f)} columns, not {len(TABLE_COLUMNS)}")
        rows.append(dict(zip(TABLE_COLUMNS, f)))
    return rows


def parse_range(tok):
    m = re.fullmatch(r"(\d+)-(\d+)", tok)
    if not m or int(m.group(1)) > int(m.group(2)) or int(m.group(1)) < 1:
        return None
    return int(m.group(1)), int(m.group(2))


def parse_truncating(spec):
    """'all' -> True, 'none' -> False, 'exons:A-B' -> (A, B); None when malformed."""
    if spec in ("all", "none"):
        return spec == "all"
    return parse_range(spec[6:]) if spec.startswith("exons:") else None


def parse_missense(spec):
    """A list of rules: ('all',), ('codon', ref, n), ('change', ref, n, alt), ('range', a, b),
    ('exons', a, b); [] for none; None when a token is malformed."""
    if spec == "none":
        return []
    rules = []
    for tok in spec.split(","):
        if tok == "all":
            rules.append(("all",))
        elif tok.startswith("exons:"):
            r = parse_range(tok[6:])
            if not r:
                return None
            rules.append(("exons",) + r)
        elif re.fullmatch(r"\d+-\d+", tok):
            r = parse_range(tok)
            if not r:
                return None
            rules.append(("range",) + r)
        else:
            m = re.fullmatch(r"([A-Z])(\d+)([A-Z]?)", tok)
            if not m or m.group(1) not in AA or (m.group(3) and m.group(3) not in AA):
                return None
            if m.group(3):
                rules.append(("change", m.group(1), int(m.group(2)), m.group(3)))
            else:
                rules.append(("codon", m.group(1), int(m.group(2))))
    return rules


def check_table(path):
    """The problems of the table, one line each ([] when it is fine)."""
    problems, seen = [], set()
    for r in read_table(path):
        g = r["gene"]
        if not re.fullmatch(r"[A-Z0-9-]+", g):
            problems.append(f"{g}: not a gene symbol")
        if g in seen:
            problems.append(f"{g}: listed twice")
        seen.add(g)
        if parse_truncating(r["truncating"]) is None:
            problems.append(f"{g}: truncating '{r['truncating']}' is not all, none or exons:A-B")
        if parse_missense(r["missense"]) is None:
            problems.append(f"{g}: missense '{r['missense']}' has a token that is not none, all, R882, V617F, "
                            "102-292 or exons:A-B")
        if r["truncating"] == "none" and r["missense"] == "none":
            problems.append(f"{g}: counts no class at all")
        if not r["note"].strip():
            problems.append(f"{g}: no note")
        if not re.fullmatch(r"\d{6,9}(;\d{6,9})*", r["pmids"]):
            problems.append(f"{g}: pmids '{r['pmids']}' is not a ;-separated list of PubMed ids")
    if not seen:
        problems.append("the table lists no gene")
    return problems


# ─── prepare: the MANE transcripts and the BED ──────────────────────────────
def attrs(col9):
    out = {}
    for kv in col9.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def ucsc(seqid):
    return "chrM" if seqid in ("MT", "M", "chrM") else (seqid if seqid.startswith("chr") else "chr" + seqid)


def prepare(gff, genes, fai, out_gff, out_bed):
    want = {r["gene"] for r in read_table(genes)}
    gene_line, gene_of = {}, {}      # gene symbol -> its line; gene id -> symbol
    tx_line, children = {}, {}       # gene symbol -> its MANE transcript line; transcript id -> lines
    with opener(gff) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            a = attrs(f[8])
            if f[2] == "gene" and a.get("Name") in want:
                gene_line[a["Name"]] = f
                gene_of[a.get("ID", "")] = a["Name"]
            elif "Parent" in a and a["Parent"] in gene_of and "MANE_Select" in a.get("tag", "").split(","):
                tx_line[gene_of[a["Parent"]]] = f
                children.setdefault(a.get("ID", ""), [])
            elif a.get("Parent") in children:
                children[a["Parent"]].append(f)
    missing = sorted(want - set(tx_line))
    if missing:
        die(f"{gff}: no MANE Select transcript for {', '.join(missing)} (an Ensembl GFF3 of release 108 or later "
            "tags them MANE_Select)")
    contigs = set()
    with open(fai) as fh:
        for ln in fh:
            contigs.add(ln.split("\t", 1)[0])
    bed = []
    with open(out_gff, "w") as out:
        out.write("##gff-version 3\n")
        for g in sorted(tx_line, key=lambda s: (ucsc(gene_line[s][0]), int(gene_line[s][3]))):
            tx = tx_line[g]
            kids = sorted(children[attrs(tx[8])["ID"]], key=lambda f: int(f[3]))
            for f in [gene_line[g], tx] + kids:
                out.write("\t".join([ucsc(f[0])] + f[1:]) + "\n")
            out.write("###\n")
            chrom = ucsc(tx[0])
            if chrom in contigs:
                for f in kids:
                    if f[2] == "CDS":
                        bed.append((chrom, max(0, int(f[3]) - 1 - PAD), int(f[4]) + PAD, g))
    with open(out_bed, "w") as out:
        for b in sorted(bed, key=lambda b: (b[0], b[1])):
            out.write("%s\t%d\t%d\t%s\n" % b)


# ─── report ─────────────────────────────────────────────────────────────────
def read_gff(path):
    """gene -> {'tx', 'chrom', 'strand', 'exons': [(start, end, rank)]}, from prepare's GFF3."""
    out, tx_gene = {}, {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            a = attrs(f[8])
            if f[2] == "gene":
                gid = a.get("ID", "")
                tx_gene[gid] = a.get("Name", "")
            elif a.get("Parent", "") in tx_gene and "transcript_id" in a:
                g = tx_gene[a["Parent"]]
                out[g] = {"tx": a["transcript_id"], "chrom": f[0], "strand": f[6], "exons": []}
                tx_gene["transcript:" + a["transcript_id"]] = g
            elif f[2] == "exon" and a.get("Parent", "").startswith("transcript:"):
                g = tx_gene.get(a["Parent"])
                if g:
                    out[g]["exons"].append((int(f[3]), int(f[4]), int(a.get("rank", 0))))
    return out


def exon_of(gene, pos):
    """The rank of the exon holding pos, or of the nearest one (a splice site sits in the intron)."""
    best = None
    for s, e, rank in gene["exons"]:
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if best is None or d < best[0]:
            best = (d, rank)
    return best[1] if best else 0


def parse_bcsq(field, tx):
    """The consequences, protein change, codon and residues of the transcript tx, from bcftools csq's
    BCSQ (run with --local-csq, so each record stands alone): 'missense|DNMT3A|ENST00000321117|
    protein_coding|-|882R>882H|25234373C>T' -> (['missense'], 'R882H', 882, 'R', 'H'). A frameshift
    reads 'R882fs', a change of several residues keeps the first one.
    A consequence bcftools marks as downstream of a stop (a leading '*') is no protein change."""
    for entry in (field or "").split(","):
        p = entry.split("|")
        if len(p) < 3 or p[2].split(".")[0] != tx or entry.startswith(("@", "*")):
            continue
        csq = p[0].split("&")
        aa = p[5] if len(p) > 5 else ""
        m = re.match(r"(\d+)([A-Z*]*)(?:\.\.\d+[A-Z*]*)?>(\d+)?([A-Z*]*)", aa)
        if m:
            ref, alt = m.group(2)[:1], m.group(4) if len(m.group(2)) == 1 else ""
            label = f"{ref}{m.group(1)}fs" if "frameshift" in csq else f"{ref}{m.group(1)}{alt}"
            return csq, label, int(m.group(1)), ref, alt
        return csq, aa, 0, "", ""
    return [], "", 0, "", ""


def in_class(row, csq, codon, ref_aa, alt_aa, exon, nexons):
    """The class a change belongs to ('' when none) and the reason."""
    trunc = parse_truncating(row["truncating"])
    hit = TRUNCATING.intersection(csq)
    if hit and trunc:
        if trunc is True or trunc[0] <= exon <= trunc[1]:
            return "truncating", f"{sorted(hit)[0]}{'' if trunc is True else f' in exon {exon}'}"
        return "", f"truncating, but in exon {exon} of {nexons} (the table counts exons {trunc[0]}-{trunc[1]})"
    if hit:
        return "", "truncating changes of this gene are not counted"
    if "missense" in csq:
        for rule in parse_missense(row["missense"]):
            kind = rule[0]
            if kind == "all":
                return "missense", "missense"
            if kind == "range" and rule[1] <= codon <= rule[2]:
                return "missense", f"missense in codons {rule[1]}-{rule[2]}"
            if kind == "exons" and rule[1] <= exon <= rule[2]:
                return "missense", f"missense in exon {exon}"
            if kind == "codon" and codon == rule[2] and ref_aa == rule[1]:
                return "missense", f"hotspot {rule[1]}{rule[2]}"
            if kind == "change" and codon == rule[2] and ref_aa == rule[1] and alt_aa == rule[3]:
                return "missense", f"{rule[1]}{rule[2]}{rule[3]}"
        return "", "a missense change outside the table's hotspots"
    return "", "a protein change outside the table's classes"


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def classify(calls, table, genes):
    """The output rows of the calls (CHROM POS REF ALT FILTER BCSQ AD AF DP), in input order."""
    rows_by_gene = {r["gene"]: r for r in table}
    out = []
    for c in calls:
        chrom, pos, ref, alt, flt, bcsq, ad, af, dp = (c + [""] * 9)[:9]
        for g, info in genes.items():
            if g not in rows_by_gene or info["chrom"] != chrom:
                continue
            csq, change, codon, ref_aa, alt_aa = parse_bcsq(bcsq, info["tx"])
            if not PROTEIN.intersection(csq):
                continue
            exon = exon_of(info, int(pos))
            cls, why = in_class(rows_by_gene[g], csq, codon, ref_aa, alt_aa, exon, len(info["exons"]))
            ads = [num(x) for x in ad.split(",")] if ad and ad != "." else []
            alt_reads = int(ads[1]) if len(ads) > 1 and ads[1] is not None else 0
            depth = num(dp)
            if depth is None and ads and None not in ads:
                depth = sum(ads)
            vaf = num(af.split(",")[0]) if af and af != "." else None
            if vaf is None and depth:
                vaf = alt_reads / depth
            if vaf is not None and vaf >= GERMLINE_VAF:
                status, note = "likely_germline", f"VAF {vaf:.2f}: near an inherited 50% or 100%, not counted"
            elif flt not in ("PASS", "."):
                status, note = "filtered", f"FilterMutectCalls: {flt}"
            elif alt_reads < MIN_ALT_READS:
                status, note = "low_support", f"{alt_reads} reads with the change, under {MIN_ALT_READS}"
            elif cls:
                status, note = "chip", why
            else:
                status, note = "other_change", why
            out.append({"gene": g, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                        "consequence": "&".join(csq), "protein_change": change, "exon": str(exon),
                        "vaf": "" if vaf is None else f"{vaf:.3f}", "alt_reads": str(alt_reads),
                        "depth": "" if depth is None else str(int(depth)), "filter": flt, "status": status,
                        "class": cls, "note": note})
            break
    return out


def p_at_least(k, n, f):
    """P(X >= k) for X ~ Binomial(n, f)."""
    return 1.0 - sum(math.comb(n, i) * f ** i * (1 - f) ** (n - i) for i in range(k))


def detection_limit(depth, k=MIN_ALT_READS, power=DETECT_POWER):
    """The smallest VAF found with probability `power` at `depth` (k reads needed); None below k reads."""
    n = int(depth)
    if n < k:
        return None
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if p_at_least(k, n, mid) >= power:
            hi = mid
        else:
            lo = mid
    return hi


def depth_rows(table, genes, bed_path, depth_path):
    """One row per table gene: coding bases on the reference, mean depth and detection limit."""
    regions = {}
    with open(bed_path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 4:
                regions.setdefault(f[3], []).append((f[0], int(f[1]), int(f[2])))
    where = {}
    for g, rs in regions.items():
        for chrom, s, e in rs:
            for p in range(s + 1, e + 1):
                where[(chrom, p)] = g
    total = {g: 0 for g in regions}
    with open(depth_path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 3:
                g = where.get((f[0], int(f[1])))
                d = num(f[2])
                if g and d is not None:
                    total[g] += d
                    where[(f[0], int(f[1]))] = None   # a base counted twice (overlapping records) counts once
    out = []
    for r in table:
        g = r["gene"]
        info = genes.get(g, {})
        bases = len({(c, p) for c, s, e in regions.get(g, []) for p in range(s + 1, e + 1)})
        if not bases:
            out.append({"gene": g, "transcript": info.get("tx", ""), "chrom": info.get("chrom", ""),
                        "coding_bases": "0", "mean_depth": "", "detection_limit_vaf": "",
                        "status": "not_on_reference"})
            continue
        mean = total.get(g, 0) / bases
        lim = detection_limit(mean)
        out.append({"gene": g, "transcript": info.get("tx", ""), "chrom": info.get("chrom", ""),
                    "coding_bases": str(bases), "mean_depth": f"{mean:.1f}",
                    "detection_limit_vaf": "" if lim is None else f"{lim:.3f}",
                    "status": "covered" if lim is not None else "too_shallow"})
    return out


def summary_line(depth):
    """The detection-limit sentence for the depth rows."""
    covered = sorted(float(r["mean_depth"]) for r in depth if r["status"] == "covered")
    if not covered:
        return ("No CHIP gene has the depth to find a clone (at least %d reads with the change): "
                "no result either way." % MIN_ALT_READS)
    med = covered[len(covered) // 2] if len(covered) % 2 else (covered[len(covered) // 2 - 1] + covered[len(covered) // 2]) / 2
    lim = detection_limit(med)
    return (f"Detection limit at this depth (median {med:.0f}x over the CHIP genes' coding exons): a clone with "
            f"a VAF of about {lim * 100:.0f}% or more is found {int(DETECT_POWER * 10)} times in 10. CHIP starts at "
            f"a VAF of {CHIP_MIN_VAF * 100:.0f}%, so smaller clones are likely missed: no call here does not rule "
            "one out. Not a diagnosis.")


def write_tsv(path, columns, rows):
    with open(path + ".tmp", "w") as out:
        out.write("\t".join(columns) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in columns) + "\n")
    os.replace(path + ".tmp", path)


def read_lines(path):
    with open(path) as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    for a in ("--gff", "--genes", "--fai", "--out-gff", "--out-bed"):
        p.add_argument(a, required=True)
    r = sub.add_parser("report")
    for a in ("--genes", "--gff", "--bed", "--calls", "--depth", "--sample", "--out", "--depth-out"):
        r.add_argument(a, required=True)
    c = sub.add_parser("check-table")
    c.add_argument("table")
    args = ap.parse_args(argv)
    if args.cmd == "check-table":
        problems = check_table(args.table)
        for pr in problems:
            sys.stderr.write(f"{args.table}: {pr}\n")
        return 1 if problems else 0
    problems = check_table(args.genes)
    if problems:
        die(f"{args.genes}: " + "; ".join(problems))
    if args.cmd == "prepare":
        prepare(args.gff, args.genes, args.fai, args.out_gff, args.out_bed)
        return 0
    table = read_table(args.genes)
    genes = read_gff(args.gff)
    rows = classify(read_lines(args.calls), table, genes)
    depth = depth_rows(table, genes, args.bed, args.depth)
    write_tsv(args.out, OUT_COLUMNS, rows)
    write_tsv(args.depth_out, DEPTH_COLUMNS, depth)
    n = sum(1 for x in rows if x["status"] == "chip")
    print(f"{args.sample}: {n} CHIP call{'s' if n != 1 else ''} in {len(table)} genes "
          f"({len(rows) - n} other protein change{'s' if len(rows) - n != 1 else ''} listed, not counted)")
    print(summary_line(depth))
    return 0


if __name__ == "__main__":
    sys.exit(main())
