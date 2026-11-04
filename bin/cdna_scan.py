#!/usr/bin/env python3
"""cdna_scan.py: where in the genome do reads jump from one exon to the next?

Genomic DNA has introns, so a read from it that reaches an exon's end runs on
into the intron. Spliced cDNA, an expression plasmid or a processed
retrocopy has none: a read across one of its exon junctions aligns as two
pieces, the first clipped at an exon's end and the rest, a supplementary
alignment (SA tag), starting at a later exon's start. A spliced aligner
writes the same read as one alignment with an N gap from exon end to exon
start. bin/clip_check.py looks for that at each ClinVar hit; this script
counts it over every exon in the genome, so contamination shows even where
it carries no ClinVar allele.

  cdna_scan.py --bam BAM --exons EXONS.bed --out OUT.tsv
  cdna_scan.py --skip REASON --out OUT.tsv

It reads the reads over each exon (EXONS.bed: chrom, 0-based start, end, the
GENCODE basic exons setup.sh installs). A read counts for a junction when one
of its pieces has a clip of at least MIN_CLIP bases within EDGE_SLACK bases
of an exon end, another piece on the same contig and strand, later in the
read (overlapping it by at most EDGE_SLACK read bases), has such a clip at a
later exon start, at least MIN_INTRON and at most MAX_INTRON bases on, or when
an N gap of its CIGAR runs from an exon end to an exon start. Supplementary
records count too, so a read whose primary alignment lies outside every exon
(on a vector contig, say) still shows its junctions. A junction counts
when at least MIN_READS reads show it. Junctions less than LOCUS_GAP bases
apart form one locus (one gene, usually), and a locus is listed when it has
at least MIN_JUNCTIONS junctions: one junction alone can be a chimeric read
pair or a deletion, two or more in one place are a spliced copy of the gene.
A listed locus is cDNA or plasmid contamination, or a processed retrocopy
the person carries (ten Berk de Boer 2023, Front Genet, PMID 37323659); the
reads tell which, not this count.

OUT.tsv starts with one comment line, "# cdna_scan: ...", saying what was
read or why nothing was, then one row per listed locus:

  chrom start end junctions reads junction_list

start is the first junction's exon end and end the last junction's exon
start (1-based, inclusive); reads counts the locus's distinct reads (a read
across two junctions counts once); junction_list is "end-start:reads" per
junction, comma-separated. --skip writes the comment line with REASON and no rows.
Needs pysam (run it in PYPGX_IMAGE); the junction logic is plain Python.
"""
import argparse
import bisect
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clip_check import EDGE_SLACK, MIN_CLIP, SKIP, parse_cigar, segment  # noqa: E402

MIN_READS = 2
MIN_JUNCTIONS = 2
MIN_INTRON = 30
MAX_INTRON = 1_000_000
LOCUS_GAP = 20_000
COLUMNS = ["chrom", "start", "end", "junctions", "reads", "junction_list"]


class Edges:
    """Exon starts and ends per contig, and the exon intervals merged."""

    def __init__(self, path):
        starts, ends, spans = {}, {}, {}
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or not f[1].isdigit() or not f[2].isdigit():
                    continue
                a, b = int(f[1]), int(f[2])
                starts.setdefault(f[0], set()).add(a)
                ends.setdefault(f[0], set()).add(b)
                spans.setdefault(f[0], []).append((a, b))
        if not starts:
            raise ValueError(f"no exon intervals in {path}")
        self.starts = {c: sorted(v) for c, v in starts.items()}
        self.ends = {c: sorted(v) for c, v in ends.items()}
        self.regions = {}
        for c, v in spans.items():
            merged = []
            for a, b in sorted(v):
                if merged and a <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            self.regions[c] = merged

    @staticmethod
    def nearest(edges, x):
        """The edge within EDGE_SLACK of x, nearest first; None when there is none."""
        if not edges:
            return None
        i = bisect.bisect_left(edges, x - EDGE_SLACK)
        best = None
        while i < len(edges) and edges[i] <= x + EDGE_SLACK:
            if best is None or abs(edges[i] - x) < abs(best - x):
                best = edges[i]
            i += 1
        return best


def gap_ok(donor, acceptor):
    return MIN_INTRON <= acceptor - donor <= MAX_INTRON


def piece(ref_start, cigar, reverse):
    """(ref_start, ref_end, left_clip, right_clip, reverse, query_start,
    query_end) of one alignment. The query interval is in the reference's
    orientation, as the CIGAR is, so two pieces on the same strand compare
    directly: query_start is the left clip, query_end adds the read bases the
    alignment holds (M, I, = and X)."""
    s = segment(ref_start, cigar)
    q = s[2] + sum(n for op, n in cigar if op in (0, 1, 7, 8))
    return s + (reverse, s[2], q)


def sa_pieces(sa, chrom):
    """The pieces of the SA tag 'chr,pos,strand,CIGAR,mapq,nm;...' on chrom."""
    out = []
    for part in (sa or "").split(";"):
        f = part.split(",")
        if len(f) >= 4 and f[0] == chrom and f[1].isdigit():
            out.append(piece(int(f[1]) - 1, parse_cigar(f[3]), f[2] == "-"))
    return out


def junctions(chrom, ref_start, cigar, sa, edges, reverse=False):
    """The (exon end, exon start) junctions one read shows on chrom: from its
    pieces (the alignment and its SA alignments) and from its N gaps. Two
    pieces pair only on the same strand and in read order: the one clipped at
    the exon end comes first in the read, the other starts where it stops
    (EDGE_SLACK bases of overlap allowed, the bases an aligner can put on
    either side of a junction)."""
    out = set()
    ends, starts = edges.ends.get(chrom), edges.starts.get(chrom)
    if not ends or not starts:
        return out
    r = ref_start
    for op, n in cigar:
        if op == SKIP:
            d, a = Edges.nearest(ends, r), Edges.nearest(starts, r + n)
            if d is not None and a is not None and gap_ok(d, a):
                out.add((d, a))
        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X move along the reference
            r += n
    pieces = [piece(ref_start, cigar, reverse)] + sa_pieces(sa, chrom)
    for p in pieces:
        if p[3] < MIN_CLIP:
            continue
        d = Edges.nearest(ends, p[1])
        if d is None:
            continue
        for q in pieces:
            if q is p or q[2] < MIN_CLIP or q[4] != p[4]:
                continue
            if not (p[5] < q[5] and p[6] <= q[5] + EDGE_SLACK):
                continue
            a = Edges.nearest(starts, q[0])
            if a is not None and gap_ok(d, a):
                out.add((d, a))
    return out


def scan(reads_by_contig, edges):
    """{(chrom, exon end, exon start): set of read keys} over pysam-like reads.
    reads_by_contig yields (chrom, reads); secondary, duplicate, QC-failed
    and unmapped alignments are skipped. Primary and supplementary records
    both count, each with its SA tag, so a read whose primary lies outside
    every exon region still shows its junctions; a read counts once per
    junction (the key is its name and mate), however many of its records
    show it."""
    found = {}
    for chrom, reads in reads_by_contig:
        for r in reads:
            if r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail or not r.cigartuples:
                continue
            sa = r.get_tag("SA") if r.has_tag("SA") else ""
            if not sa and not any(op == SKIP for op, _ in r.cigartuples):
                continue
            for d, a in junctions(chrom, r.reference_start, r.cigartuples, sa, edges, r.is_reverse):
                found.setdefault((chrom, d, a), set()).add((r.query_name, r.is_read1))
    return found


def loci(found):
    """The listed loci: rows of COLUMNS, in genome order of their contig as
    first seen and then position. A locus's reads are its distinct reads: a
    read across two of its junctions counts once there, and on each junction."""
    kept = sorted((c, d, a) for (c, d, a), k in found.items() if len(k) >= MIN_READS)
    groups, cur = [], None
    for c, d, a in kept:
        if cur and cur[0] == c and d <= cur[2] + LOCUS_GAP:
            cur[1].append((d, a))
            cur[2] = max(cur[2], a)
        else:
            cur = [c, [(d, a)], a]
            groups.append(cur)
    out = []
    for c, js, _ in groups:
        if len(js) < MIN_JUNCTIONS:
            continue
        reads = set().union(*(found[(c, d, a)] for d, a in js))
        out.append([c, str(min(d for d, _ in js)), str(max(a for _, a in js) + 1), str(len(js)),
                    str(len(reads)), ",".join(f"{d}-{a + 1}:{len(found[(c, d, a)])}" for d, a in js)])
    return out


def bam_reads(bam, edges):
    """(chrom, reads) over the merged exon regions of each contig the BAM has."""
    for chrom in bam.references:
        for a, b in edges.regions.get(chrom, []):
            yield chrom, bam.fetch(chrom, a, b)


def write(path, comment, rows):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(f"# cdna_scan: {comment}\n")
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--bam")
    ap.add_argument("--exons")
    ap.add_argument("--skip", help="write no rows, with this reason")
    a = ap.parse_args(argv)
    if a.skip:
        write(a.out, f"not run: {a.skip}", [])
        print(f"cDNA scan: not run ({a.skip})")
        return 0
    if not (a.bam and a.exons):
        ap.error("give --bam and --exons, or --skip REASON")
    try:
        import pysam  # only here: --skip runs without it
        edges = Edges(a.exons)
        bam = pysam.AlignmentFile(a.bam)
        shared = [c for c in bam.references if c in edges.regions]
        if not shared:
            raise ValueError(f"the BAM and {a.exons} have no contig name in common")
        rows = loci(scan(bam_reads(bam, edges), edges))
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    n_regions = sum(len(edges.regions[c]) for c in shared)
    write(a.out, f"read {n_regions} exon regions on {len(shared)} contigs; {len(rows)} loci listed", rows)
    print(f"cDNA scan: {len(rows)} loci with reads across {MIN_JUNCTIONS} or more exon junctions "
          f"({n_regions} exon regions read)")
    for r in rows:
        print(f"  SPLICED READS {r[0]}:{r[1]}-{r[2]}: {r[3]} exon junctions, {r[4]} reads "
              "(cDNA or plasmid contamination, or a processed retrocopy: check the reads)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
