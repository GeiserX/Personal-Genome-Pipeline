#!/usr/bin/env python3
"""clip_check.py: do the alt reads of a ClinVar or CPSR hit look like spliced cDNA?

Spliced cDNA or expression-vector DNA in a sample (a known lab contaminant)
aligns to a gene's exons. Its reads that cross an exon junction are clipped
exactly at the exon edge, with the rest of the read on the next exon as a
supplementary alignment (SA tag). A variant the construct carries is then
called, and the screen lists it as a ClinVar hit. Genomic reads at the same
site run on into the intron and are not clipped there.

For each hit this counts, among the reads that carry the alt allele and among
those that carry the ref allele, the reads with a clip of at least MIN_CLIP
bases (soft or hard) whose clip point lies within EDGE_SLACK bases of an
exon edge on the same side (a clip on the right end at an exon's end, on the
left end at an exon's start), in the read's own alignment or in one of its SA
alignments. A hit is flagged when at least MIN_ALT_READS reads carry the alt
allele, at least ALT_SHARE_MIN of them are clipped that way, and at most
REF_SHARE_MAX of the ref reads are.

  clip_check.py --hits HITS.tsv --bam BAM --exons EXONS.bed --out OUT.tsv
  clip_check.py --cpsr CLASSIFICATION.tsv.gz --bam BAM --exons EXONS.bed --out OUT.tsv
  clip_check.py --hits HITS.tsv --skip REASON --out OUT.tsv

HITS.tsv is step 06's hits table (columns chrom, pos, ref, alt are read by
name). --cpsr reads step 17's CPSR classification table instead and checks
its P/LP rows: the rows whose final class (CLASSIFICATION) or CPSR's own
class (CPSR_CLASSIFICATION) is Pathogenic or Likely pathogenic, at the
position of their VAR_ID (chrom_pos_ref_alt) or GENOMIC_CHANGE
(chrom:g.posREF>ALT). CPSR writes contig names without "chr"; a name the BAM
lacks is also looked up with "chr" added or removed, and OUT.tsv keeps
CPSR's name. EXONS.bed holds exon intervals, 0-based and half-open (chrom, start,
end), such as the GENCODE basic exons setup.sh installs. --skip writes every
hit as not checked, with REASON. OUT.tsv has one row per hit:

  chrom pos ref alt alt_reads alt_clipped ref_reads ref_clipped cdna_check cdna_note

cdna_check is flagged, clear or not_checked; cdna_note says why. Needs pysam
(run it in PYPGX_IMAGE); the counting itself is plain Python.
"""
import argparse
import bisect
import csv
import gzip
import os
import re
import sys

MIN_CLIP = 10
EDGE_SLACK = 3
MIN_ALT_READS = 2
ALT_SHARE_MIN = 0.5
REF_SHARE_MAX = 0.2

FLAG_NOTE = "possible cDNA/transgene contamination: check the reads"
COLUMNS = ["chrom", "pos", "ref", "alt", "alt_reads", "alt_clipped", "ref_reads", "ref_clipped",
           "cdna_check", "cdna_note"]

# CIGAR operation codes, as pysam and the BAM format number them.
M, INS, DEL, SKIP, SOFT, HARD, PAD, EQ, DIFF = range(9)
CIGAR_OPS = {"M": M, "I": INS, "D": DEL, "N": SKIP, "S": SOFT, "H": HARD, "P": PAD, "=": EQ, "X": DIFF}


def parse_cigar(text):
    """'75M25S' -> [(0, 75), (4, 25)]."""
    out, n = [], ""
    for ch in text:
        if ch.isdigit():
            n += ch
        else:
            out.append((CIGAR_OPS[ch], int(n)))
            n = ""
    return out


def segment(ref_start, cigar):
    """(ref_start, ref_end, left_clip, right_clip) of one alignment: 0-based,
    end exclusive, clips counting soft and hard clipped bases."""
    ref_len = sum(n for op, n in cigar if op in (M, DEL, SKIP, EQ, DIFF))
    left = right = 0
    for op, n in cigar:
        if op not in (SOFT, HARD):
            break
        left += n
    for op, n in reversed(cigar):
        if op not in (SOFT, HARD):
            break
        right += n
    return ref_start, ref_start + ref_len, left, right


def sa_segments(sa):
    """The SA tag 'chr,pos,strand,CIGAR,mapq,nm;...' as (chrom, segment) pairs."""
    out = []
    for part in (sa or "").split(";"):
        f = part.split(",")
        if len(f) >= 4 and f[1].isdigit():
            out.append((f[0], segment(int(f[1]) - 1, parse_cigar(f[3]))))
    return out


class Exons:
    """Exon edges per contig: starts (left edges) and ends (right edges)."""

    def __init__(self, path):
        starts, ends = {}, {}
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or not f[1].isdigit() or not f[2].isdigit():
                    continue
                starts.setdefault(f[0], set()).add(int(f[1]))
                ends.setdefault(f[0], set()).add(int(f[2]))
        if not starts:
            raise ValueError(f"no exon intervals in {path}")
        self.starts = {c: sorted(v) for c, v in starts.items()}
        self.ends = {c: sorted(v) for c, v in ends.items()}

    @staticmethod
    def _near(edges, x):
        if not edges:
            return False
        i = bisect.bisect_left(edges, x - EDGE_SLACK)
        return i < len(edges) and edges[i] <= x + EDGE_SLACK

    def edge_clipped(self, chrom, seg):
        """True when the alignment has a clip of MIN_CLIP or more at an exon edge
        on the same side."""
        start, end, left, right = seg
        return ((left >= MIN_CLIP and self._near(self.starts.get(chrom), start))
                or (right >= MIN_CLIP and self._near(self.ends.get(chrom), end)))

    def read_clipped(self, chrom, seg, sa):
        return self.edge_clipped(chrom, seg) or any(self.edge_clipped(c, s) for c, s in sa_segments(sa))


def walk(ref_start, cigar, seq):
    """The read base at each reference position (None where the read has a
    deletion) and the bases inserted after each reference position."""
    bases, ins = {}, {}
    r, q = ref_start, 0
    for op, n in cigar:
        if op in (M, EQ, DIFF):
            for i in range(n):
                bases[r + i] = seq[q + i]
            r += n
            q += n
        elif op == INS:
            ins[r - 1] = ins.get(r - 1, "") + seq[q:q + n]
            q += n
        elif op in (DEL, SKIP):
            for i in range(n):
                bases[r + i] = None
            r += n
        elif op == SOFT:
            q += n
    return bases, ins


def allele(bases, ins, pos0, ref, alt):
    """'alt', 'ref' or None (the read does not show either, or does not span
    the site) for the VCF record at 0-based pos0."""
    if len(ref) == len(alt):
        got = [bases.get(pos0 + i) for i in range(len(ref))]
        if None in got:
            return None
        got = "".join(got).upper()
        return "alt" if got == alt.upper() else "ref" if got == ref.upper() else None
    anchor = min(len(ref), len(alt))
    if ref[:anchor].upper() != alt[:anchor].upper():
        return None
    last = pos0 + anchor - 1  # the last anchor base
    if any(bases.get(pos0 + i) is None for i in range(anchor)):
        return None
    if len(ref) > len(alt):
        gap = range(last + 1, pos0 + len(ref))
        after = pos0 + len(ref)
        if last in ins or after not in bases or bases[after] is None:
            return None
        if all(p in bases and bases[p] is None for p in gap):
            return "alt"
        if all(bases.get(p) is not None for p in gap):
            return "ref"
        return None
    if last + 1 not in bases or bases[last + 1] is None:
        return None
    got = ins.get(last, "").upper()
    if got == alt[anchor:].upper():
        return "alt"
    return "ref" if got == "" else None


def complex_allele(ref, alt):
    n = min(len(ref), len(alt))
    return len(ref) != len(alt) and ref[:n].upper() != alt[:n].upper() or not ref.isalpha() or not alt.isalpha()


def verdict(alt_n, alt_c, ref_n, ref_c):
    """(cdna_check, cdna_note) from the four counts."""
    if alt_n < MIN_ALT_READS:
        return "not_checked", f"fewer than {MIN_ALT_READS} reads carry the alt allele"
    if alt_c >= ALT_SHARE_MIN * alt_n and (ref_n == 0 or ref_c <= REF_SHARE_MAX * ref_n):
        return "flagged", FLAG_NOTE
    return "clear", ""


def count_reads(reads, exons, chrom, pos0, ref, alt):
    """[alt_reads, alt_clipped, ref_reads, ref_clipped] over pysam-like reads.
    A read counts once (its name and mate number), as its first alignment
    seen; secondary, duplicate, QC-failed and unmapped alignments are skipped."""
    seen = set()
    n = {"alt": [0, 0], "ref": [0, 0]}
    for r in reads:
        if r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail or not r.cigartuples:
            continue
        key = (r.query_name, r.is_read1, r.is_read2)
        if key in seen:
            continue
        seq = r.query_sequence or ""
        a = allele(*walk(r.reference_start, r.cigartuples, seq), pos0, ref, alt) if seq else None
        if a is None:
            continue
        seen.add(key)
        n[a][0] += 1
        sa = r.get_tag("SA") if r.has_tag("SA") else ""
        if exons.read_clipped(chrom, segment(r.reference_start, r.cigartuples), sa):
            n[a][1] += 1
    return n["alt"] + n["ref"]


def read_hits(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        missing = [k for k in ("chrom", "pos", "ref", "alt") if k not in r]
        if missing:
            raise ValueError(f"{path} has no column {missing[0]}")
    return rows


def norm_class(v):
    """A CPSR or ClinVar class in one spelling: 'Likely_Pathogenic' and 'LP'
    -> 'likely pathogenic'; '' for an empty or NA cell."""
    v = (v or "").replace("_", " ").strip().lower()
    return {"p": "pathogenic", "lp": "likely pathogenic", "b": "benign", "lb": "likely benign",
            "na": "", ".": ""}.get(v, v)


PLP = ("pathogenic", "likely pathogenic")


def cpsr_variant(row):
    """(chrom, pos, ref, alt) of a CPSR row from VAR_ID or GENOMIC_CHANGE;
    None when neither can be read."""
    f = (row.get("VAR_ID") or "").split("_")
    if len(f) >= 4 and f[-3].isdigit():
        return "_".join(f[:-3]), int(f[-3]), f[-2], f[-1]
    m = re.match(r"^(\S+):g\.(\d+)([A-Za-z]+)>([A-Za-z]+)$", row.get("GENOMIC_CHANGE") or "")
    if m:
        return m.group(1), int(m.group(2)), m.group(3), m.group(4)
    return None


def read_cpsr(path):
    """The P/LP rows of a CPSR classification table as hits (chrom, pos, ref, alt)."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    for r in rows:
        classes = [norm_class(r.get(c)) for c in ("CLASSIFICATION", "FINAL_CLASSIFICATION", "CPSR_CLASSIFICATION")]
        v = cpsr_variant(r)
        if v and any(c in PLP for c in classes):
            out.append(dict(zip(("chrom", "pos", "ref", "alt"), (v[0], str(v[1]), v[2], v[3]))))
    return out


def bam_contig(chrom, references):
    """The BAM's name for chrom: itself, or with "chr" added or removed."""
    for c in (chrom, "chr" + chrom, chrom[3:] if chrom.startswith("chr") else None):
        if c and c in references:
            return c
    return None


def check(hits, bam_path, exons_path, skip):
    """One output row per hit."""
    out = []
    bam = exons = None
    if not skip:
        import pysam  # only here: --skip runs without it
        exons = Exons(exons_path)
        bam = pysam.AlignmentFile(bam_path)
    for h in hits:
        chrom, pos, ref, alt = h["chrom"], int(h["pos"]), h["ref"], h["alt"]
        row = [chrom, str(pos), ref, alt, "", "", "", ""]
        if skip:
            out.append(row + ["not_checked", skip])
            continue
        if complex_allele(ref, alt):
            out.append(row + ["not_checked", "not a SNV, MNV or simple indel"])
            continue
        contig = bam_contig(chrom, bam.references)
        if contig is None:
            out.append(row + ["not_checked", f"the BAM has no contig {chrom}"])
            continue
        pos0 = pos - 1
        counts = count_reads(bam.fetch(contig, pos0, pos0 + max(len(ref), 1) + 1), exons, contig, pos0, ref, alt)
        out.append([chrom, str(pos), ref, alt] + [str(x) for x in counts] + list(verdict(*counts)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--hits", help="step 06's hits table")
    src.add_argument("--cpsr", help="step 17's CPSR classification table (its P/LP rows are checked)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--bam")
    ap.add_argument("--exons")
    ap.add_argument("--skip", help="write every hit as not checked, with this reason")
    a = ap.parse_args(argv)
    if not a.skip and not (a.bam and a.exons):
        ap.error("give --bam and --exons, or --skip REASON")
    try:
        hits = read_cpsr(a.cpsr) if a.cpsr else read_hits(a.hits)
        rows = check(hits, a.bam, a.exons, a.skip)
    except (OSError, ValueError, EOFError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    os.replace(tmp, a.out)
    flagged = [r for r in rows if r[8] == "flagged"]
    clear = sum(1 for r in rows if r[8] == "clear")
    print(f"cDNA read check: {len(flagged)} flagged, {clear} clear, {len(rows) - len(flagged) - clear} not checked")
    for r in flagged:
        print(f"  FLAGGED {r[0]}:{r[1]} {r[2]}>{r[3]}: {r[5]} of {r[4]} alt reads and {r[7]} of {r[6]} ref reads "
              f"clipped at an exon edge ({FLAG_NOTE})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
