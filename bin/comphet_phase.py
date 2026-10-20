#!/usr/bin/env python3
"""comphet_phase.py: are the two variants of a compound-het candidate pair on
the same copy of the gene (cis) or on different copies (trans)?

slivar compound-hets lists every pair of heterozygous candidates in a gene.
From one unphased genome a pair in cis (both variants on one copy, the other
copy intact) looks the same as a pair in trans (one variant on each copy, no
intact copy left). The reads can tell when one read, or one read pair, covers
both positions: in cis a read shows both alt alleles or neither, in trans it
shows one of them.

For each pair this counts the reads and read pairs that show an allele at both
positions: CIS_READS show alt and alt or ref and ref, TRANS_READS show alt and
ref. A read pair counts once, with what its two mates show together (a pair
whose mates disagree at a position is not counted). The verdict is trans or
cis when at least MIN_READS reads or pairs show the two positions and at least
MIN_SHARE of them agree; otherwise unknown, with the reason in PHASE_NOTE.
Two positions further apart than a read pair spans are never covered by one
fragment and stay unknown: read-backed phasing over the whole gene (WhatsHap
on long reads, or a parent's genome, step 38) can phase those.

Reads that count: mapped, primary, not duplicates, not QC-failed, mapping
quality at least MIN_MAPQ. The allele a read shows is bin/clip_check.py's rule
(SNV, MNV, or an indel that shares its first base); a complex allele is not
phased.

  comphet_phase.py --vcf COMPHET.vcf.gz --bam BAM --out PAIRS.tsv \\
      [--variants COMPHET.tsv --variants-out OUT.tsv]
  comphet_phase.py --vcf COMPHET.vcf.gz --skip REASON --out PAIRS.tsv [...]

COMPHET.vcf.gz is slivar compound-hets' VCF: its INFO/slivar_comphet entries
(sample/gene/pair_id/chrom/pos/ref/alt) name the pairs. PAIRS.tsv has one row
per pair:

  GENE CHROM POS_1 REF_1 ALT_1 POS_2 REF_2 ALT_2 DISTANCE PHASE CIS_READS TRANS_READS PHASE_NOTE

--variants takes step 31's per-variant compound-het table (columns CHROM, POS,
REF, ALT read by name) and writes it to --variants-out with a PHASE column
after GT: trans when the variant is in trans with at least one partner, else
unknown when one of its pairs is unknown, else cis (every partner is on the
same copy). --skip writes every pair as unknown with REASON (no BAM).
Needs pysam for --bam (PYPGX_IMAGE); the counting is plain Python.
"""
import argparse
import gzip
import os
import sys

from clip_check import allele, complex_allele, walk

MIN_MAPQ = 20
MIN_READS = 2
MIN_SHARE = 0.8
# A read pair is one fragment: its two mates are joined by name. The unit
# test turns this off to show that a trans pair seen only across mates then
# reads unknown.
JOIN_MATES = True

COLUMNS = ["GENE", "CHROM", "POS_1", "REF_1", "ALT_1", "POS_2", "REF_2", "ALT_2", "DISTANCE",
           "PHASE", "CIS_READS", "TRANS_READS", "PHASE_NOTE"]


def open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_pairs(path):
    """The pairs of slivar's VCF: a sorted list of (gene, var1, var2), each var
    a (chrom, pos, ref, alt) tuple. The two records of a pair share its id;
    each entry names the record's partner, so both variants come from it."""
    groups = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            own = (f[0], int(f[1]), f[3], f[4])
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in f[7].split(";"))
            for entry in filter(None, info.get("slivar_comphet", "").split(",")):
                p = entry.split("/")
                if len(p) < 7:
                    raise ValueError(f"unreadable slivar_comphet entry {entry!r}")
                gene, pid = "/".join(p[1:-5]), p[-5]
                other = (p[-4], int(p[-3]), p[-2], p[-1])
                groups.setdefault((gene, pid), set()).update((own, other))
    pairs = set()
    for (gene, pid), vs in groups.items():
        vs = sorted(vs, key=var_key)
        if len(vs) != 2:
            print(f"WARNING: pair {pid} in {gene} names {len(vs)} variants, not 2; skipped", file=sys.stderr)
            continue
        pairs.add((gene, vs[0], vs[1]))
    return sorted(pairs, key=lambda p: (p[0], var_key(p[1]), var_key(p[2])))


def var_key(v):
    return (v[0], v[1], v[2], v[3])


def counts_read(r):
    return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail
                or not r.cigartuples or r.mapping_quality < MIN_MAPQ)


def site_alleles(reads, var):
    """{fragment: 'alt' or 'ref'} for the reads at var. A fragment is a read
    pair (or one read with JOIN_MATES off); one whose reads disagree is left
    out."""
    out, bad = {}, set()
    _, pos, ref, alt = var
    for r in reads:
        if not counts_read(r):
            continue
        bases, ins = walk(r.reference_start, r.cigartuples, r.query_sequence)
        a = allele(bases, ins, pos - 1, ref, alt)
        if a is None:
            continue
        key = r.query_name if JOIN_MATES else (r.query_name, r.is_read1)
        if key in out and out[key] != a:
            bad.add(key)
        out[key] = a
    for k in bad:
        del out[k]
    return out


def phase_counts(a1, a2):
    """(cis, trans): fragments that show both positions, from two site_alleles maps."""
    cis = trans = 0
    for k, x in a1.items():
        y = a2.get(k)
        if y is None:
            continue
        if x == y:
            cis += 1
        else:
            trans += 1
    return cis, trans


def verdict(cis, trans):
    """(PHASE, PHASE_NOTE) from the two counts."""
    n = cis + trans
    if n == 0:
        return "unknown", "no read or read pair covers both positions"
    if n < MIN_READS:
        return "unknown", f"only {n} read or read pair covers both positions"
    if trans >= MIN_SHARE * n:
        return "trans", f"{trans} of {n} reads or read pairs show the two variants on different copies"
    if cis >= MIN_SHARE * n:
        return "cis", f"{cis} of {n} reads or read pairs show the two variants on the same copy"
    return "unknown", f"the reads disagree: {cis} cis, {trans} trans"


class Reads:
    """The fragments at each variant, read once from the BAM."""

    def __init__(self, bam):
        self.bam = bam
        self.refs = set(bam.references)
        self.cache = {}

    def contig(self, chrom):
        for c in (chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom):
            if c in self.refs:
                return c
        return None

    def at(self, var):
        if var not in self.cache:
            chrom, pos, ref, _ = var
            c = self.contig(chrom)
            self.cache[var] = site_alleles(self.bam.fetch(c, pos - 1, pos - 1 + max(len(ref), 1) + 1), var)
        return self.cache[var]


def phase(pairs, reads, skip):
    rows = []
    for gene, v1, v2 in pairs:
        if skip:
            ph, c, t, note = "unknown", ".", ".", skip
        elif reads.contig(v1[0]) is None:
            ph, c, t, note = "unknown", ".", ".", f"the BAM has no contig {v1[0]}"
        elif v1[0] != v2[0]:
            ph, c, t, note = "unknown", ".", ".", "the two variants are on different contigs"
        elif complex_allele(v1[2], v1[3]) or complex_allele(v2[2], v2[3]):
            ph, c, t, note = "unknown", ".", ".", "a complex allele: not phased from the reads"
        else:
            c, t = phase_counts(reads.at(v1), reads.at(v2))
            ph, note = verdict(c, t)
        dist = str(abs(v2[1] - v1[1])) if v1[0] == v2[0] else "."
        rows.append([gene, v1[0], str(v1[1]), v1[2], v1[3], str(v2[1]), v2[2], v2[3], dist, ph, str(c), str(t), note])
    return rows


def variant_phase(rows):
    """{(chrom, pos, ref, alt): PHASE} over each variant's pairs."""
    seen = {}
    for r in rows:
        for v in ((r[1], r[2], r[3], r[4]), (r[1], r[5], r[6], r[7])):
            seen.setdefault(v, set()).add(r[9])
    order = ("trans", "unknown", "cis")
    return {v: next(x for x in order if x in s) for v, s in seen.items()}


def add_phase_column(src, dst, rows):
    best = variant_phase(rows)
    tmp = dst + ".tmp"
    with open(src) as fh, open(tmp, "w") as out:
        head = fh.readline().rstrip("\n").split("\t")
        for k in ("CHROM", "POS", "REF", "ALT"):
            if k not in head:
                raise ValueError(f"{src} has no {k} column")
        at = head.index("GT") + 1 if "GT" in head else len(head)
        col = {k: head.index(k) for k in ("CHROM", "POS", "REF", "ALT")}
        out.write("\t".join(head[:at] + ["PHASE"] + head[at:]) + "\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = tuple(f[col[k]] for k in ("CHROM", "POS", "REF", "ALT"))
            f.insert(at, best.get(v, "."))
            out.write("\t".join(f) + "\n")
    os.replace(tmp, dst)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf", required=True, help="slivar compound-hets' VCF")
    ap.add_argument("--out", required=True, help="the pairs table")
    ap.add_argument("--bam")
    ap.add_argument("--skip", help="write every pair as unknown, with this reason")
    ap.add_argument("--variants", help="step 31's per-variant compound-het table")
    ap.add_argument("--variants-out")
    a = ap.parse_args(argv)
    if bool(a.bam) == bool(a.skip):
        ap.error("give --bam or --skip REASON")
    if bool(a.variants) != bool(a.variants_out):
        ap.error("--variants and --variants-out go together")
    try:
        pairs = read_pairs(a.vcf)
        reads = None
        if a.bam:
            import pysam
            reads = Reads(pysam.AlignmentFile(a.bam, "rb"))
        rows = phase(pairs, reads, a.skip)
        if a.variants:
            add_phase_column(a.variants, a.variants_out, rows)
    except (OSError, ValueError, EOFError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    os.replace(tmp, a.out)
    n = {k: sum(1 for r in rows if r[9] == k) for k in ("trans", "cis", "unknown")}
    print(f"Phase of {len(rows)} compound-het pairs from the reads: "
          f"{n['trans']} trans, {n['cis']} cis, {n['unknown']} unknown"
          + (f" ({a.skip})" if a.skip else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
