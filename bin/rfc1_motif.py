#!/usr/bin/env python3
"""rfc1_motif.py: which repeat motif do the reads at RFC1 carry?

ExpansionHunter sizes RFC1 with the degenerate unit AARRG, and Stranger's
catalog marks it full_mutation from a low count, so a common benign
expansion is flagged as readily as a pathogenic one. What matters is the
motif. The reference allele is (AAAAG)11. Expanded AAAAG and AAAGG alleles
are common and benign; biallelic AAGGG expansions cause CANVAS (Cortese et
al., Nat Genet 2019;51:649-658); ACAGG expansions are pathogenic too
(Scriba et al., Brain 2020;143:2904-2910; Beecroft et al., Brain
2020;143:2673-2680). Short reads cannot size a long expansion, but they
can say which motif it is made of.

For every read at the locus this takes the read's bases over the repeat
and FLANK bases either side, the soft-clipped bases next to them (a read
whose repeat part did not align is clipped there), and the whole sequence
of two kinds of read the aligner did not place in the repeat: an unmapped
read placed at its mate's position within FETCH_PAD bases of the repeat,
and the mate of a read there when that mate is mapped somewhere else
(mate rescue). It finds the
longest run of one pentamer repeated in tandem in those bases, and when
the run has at least MIN_UNITS units it counts the read for that motif.
A motif is named by its known form whatever its rotation and strand
(AAGGG, GGGAA and CCCTT are one motif). Duplicate, secondary,
supplementary and QC-failed alignments are not counted.

The verdict, from the reads that carry a run:
  pathogenic_motif  at least MIN_READS reads and PATHOGENIC_SHARE of them
                    carry AAGGG or ACAGG: confirm with repeat-primed PCR
  benign_motif      otherwise, when a benign motif (AAAAG, AAAGG) is the
                    most common and carries BENIGN_SHARE of the reads: the
                    flag is not supported
  uncertain_motif   otherwise, with MIN_READS reads or more: the dominant
                    motif is not a known benign one
  not_resolved      fewer than MIN_READS reads carry a run
  not_checked       --skip REASON (no BAM, for one)

  rfc1_motif.py --bam BAM --vcf STRANGER_VCF --out OUT.tsv
  rfc1_motif.py --vcf STRANGER_VCF --skip REASON --out OUT.tsv

The locus is the RFC1 record of the VCF (REPID=RFC1): its POS to its END.
OUT.tsv has a header and one row:

  locus chrom start end repeat_reads dominant_motif dominant_class
  dominant_share motif_shares verdict note

motif_shares lists every motif seen as MOTIF=share, most common first.
--no-softclip and --no-mates turn the two rescues off (for the tests).
Needs pysam (run it in PYPGX_IMAGE); the counting itself is plain Python.
"""
import argparse
import os
import sys

FLANK = 50
FETCH_PAD = 1000
MIN_UNITS = 4
MIN_READS = 3
PATHOGENIC_SHARE = 0.2
BENIGN_SHARE = 0.5

KNOWN = {"AAAAG": "benign", "AAAGG": "benign", "AAGGG": "pathogenic", "ACAGG": "pathogenic"}
COLUMNS = ["locus", "chrom", "start", "end", "repeat_reads", "dominant_motif", "dominant_class",
           "dominant_share", "motif_shares", "verdict", "note"]

# CIGAR operation codes, as pysam and the BAM format number them.
M, INS, DEL, SKIP, SOFT, HARD, PAD, EQ, DIFF = range(9)
CIGAR_OPS = {"M": M, "I": INS, "D": DEL, "N": SKIP, "S": SOFT, "H": HARD, "P": PAD, "=": EQ, "X": DIFF}
COMP = str.maketrans("ACGTN", "TGCAN")


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


def revcomp(s):
    return s.translate(COMP)[::-1]


def forms(m):
    """Every rotation of a motif on both strands."""
    out = set()
    for x in (m, revcomp(m)):
        out.update(x[i:] + x[:i] for i in range(len(x)))
    return out


FORM_TO_KNOWN = {f: k for k in KNOWN for f in forms(k)}


def motif_name(m):
    """The known name of a pentamer, or its smallest rotation on either strand."""
    return FORM_TO_KNOWN.get(m) or min(forms(m))


def longest_run(seq):
    """(motif, units) of the longest tandem run of one pentamer in seq, or
    (None, 0). A run of a homopolymer (period 1) is not a pentamer run."""
    seq = seq.upper()
    best, best_len, i, n = None, 0, 0, len(seq)
    while i + 5 < n:
        j = i
        while j + 5 < n and seq[j] == seq[j + 5] and seq[j] != "N":
            j += 1
        length = j - i + 5
        unit = seq[i:i + 5]
        if j > i and length > best_len and len(set(unit)) > 1:
            best, best_len = unit, length
        i = j + 1 if j > i else i + 1
    return (motif_name(best), best_len // 5) if best else (None, 0)


def repeat_bases(read, lo, hi, softclip=True):
    """The read's bases aligned within [lo, hi) (0-based), with the soft clip
    next to them when softclip is set; '' when none are."""
    seq = read.query_sequence or ""
    cig = read.cigartuples or []
    r, q, first, last = read.reference_start, 0, None, None
    for op, n in cig:
        if op in (M, EQ, DIFF):
            a, b = max(r, lo), min(r + n, hi)
            if a < b:
                qa, qb = q + (a - r), q + (b - r)
                first = qa if first is None else first
                last = qb
            r += n
            q += n
        elif op in (DEL, SKIP):
            r += n
        elif op in (INS, SOFT):
            if op == INS and first is not None and lo <= r < hi:
                last = q + n
            q += n
    if first is None:
        return ""
    if softclip and False:  # RED CONTROL: soft-clip rescue off; reverted in the next commit
        if cig and cig[0][0] == SOFT and first == cig[0][1]:
            first = 0
        if cig and cig[-1][0] == SOFT and last == len(seq) - cig[-1][1]:
            last = len(seq)
    return seq[first:last]


def usable(read):
    return not (read.is_duplicate or read.is_secondary or read.is_supplementary or read.is_qcfail)


def read_segments(reads, lo, hi, mate_of=None, softclip=True, fetched=None):
    """(key, bases) for every read fetched around the locus, once per read:
    the bases over [lo, hi) for a mapped read, the whole sequence for an
    unmapped read placed in the fetched window [fetched) and for a mate
    mapped outside it (mate_of(read) returns it or None; None turns the
    mate rescue off)."""
    flo, fhi = fetched or (lo, hi)
    seen = set()
    for read in reads:
        if not usable(read):
            continue
        key = (read.query_name, read.is_read1)
        if key in seen:
            continue
        seen.add(key)
        if read.is_unmapped:
            yield key, read.query_sequence or ""
            continue
        yield key, repeat_bases(read, lo, hi, softclip)
        if mate_of is None or not read.is_paired or read.mate_is_unmapped:
            continue
        mkey = (read.query_name, not read.is_read1)
        if mkey in seen:
            continue
        if read.next_reference_id == read.reference_id and flo <= read.next_reference_start < fhi:
            continue  # the mate is in the fetched window: the fetch returns it
        mate = mate_of(read)
        if mate is not None and usable(mate):
            seen.add(mkey)
            yield mkey, mate.query_sequence or ""


def count_motifs(segments):
    counts = {}
    for _, bases in segments:
        motif, units = longest_run(bases)
        if motif and units >= MIN_UNITS:
            counts[motif] = counts.get(motif, 0) + 1
    return counts


def verdict(counts):
    """(dominant_motif, dominant_class, dominant_share, shares, verdict, note)."""
    total = sum(counts.values())
    order = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    shares = ",".join(f"{m}={c / total:.2f}" for m, c in order) if total else "."
    if total < MIN_READS:
        return (".", ".", ".", shares, "not_resolved",
                f"{total} read(s) carry {MIN_UNITS} or more repeat units at RFC1, fewer than {MIN_READS}: "
                "the motif is not resolved; confirm the flag with repeat-primed PCR if symptoms fit")
    dom, n = order[0]
    cls = KNOWN.get(dom, "unknown")
    share = n / total
    path = {m: c for m, c in counts.items() if KNOWN.get(m) == "pathogenic"}
    npath = sum(path.values())
    if npath >= MIN_READS and npath / total >= PATHOGENIC_SHARE:
        names = "/".join(sorted(path, key=lambda m: -path[m]))
        return (dom, cls, f"{share:.2f}", shares, "pathogenic_motif",
                f"{names} in {npath} of {total} repeat reads ({npath / total:.0%}): confirm with repeat-primed PCR "
                "(CANVAS needs a pathogenic expansion on both copies)")
    if cls == "benign" and share >= BENIGN_SHARE:
        return (dom, cls, f"{share:.2f}", shares, "benign_motif",
                f"motif {dom} (benign) in {n} of {total} repeat reads ({share:.0%}): the flag is not supported")
    return (dom, cls, f"{share:.2f}", shares, "uncertain_motif",
            f"most common motif {dom} ({cls}) in {n} of {total} repeat reads ({share:.0%}): "
            "not a known benign expansion; confirm with repeat-primed PCR if symptoms fit")


def rfc1_locus(vcf):
    """(chrom, start, end) of the REPID=RFC1 record: 0-based start, end exclusive."""
    import gzip
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = dict(x.split("=", 1) if "=" in x else (x, "") for x in f[7].split(";"))
            if info.get("REPID") == "RFC1" or info.get("VARID") == "RFC1":
                start = int(f[1]) - 1
                end = int(info["END"]) if info.get("END", "").isdigit() else start + len(f[3])
                return f[0], start, end
    return None


def check(bam_path, vcf, skip=None, softclip=True, mates=True):
    locus = rfc1_locus(vcf)
    if locus is None:
        return ["RFC1", ".", ".", ".", "0", ".", ".", ".", ".", "not_checked", "no RFC1 record in the VCF"]
    chrom, start, end = locus
    head = ["RFC1", chrom, str(start + 1), str(end)]
    if skip:
        return head + ["0", ".", ".", ".", ".", "not_checked", skip]
    import pysam
    lo, hi = start - FLANK, end + FLANK
    with pysam.AlignmentFile(bam_path) as bam:
        if chrom not in bam.references:
            return head + ["0", ".", ".", ".", ".", "not_checked", f"the BAM has no contig {chrom}"]
        other = pysam.AlignmentFile(bam_path)

        def mate_of(read):
            try:
                return other.mate(read)
            except ValueError:
                return None
        flo, fhi = max(0, start - FETCH_PAD), end + FETCH_PAD
        try:
            counts = count_motifs(read_segments(bam.fetch(chrom, flo, fhi), lo, hi,
                                                mate_of if mates else None, softclip, (flo, fhi)))
        finally:
            other.close()
    total = sum(counts.values())
    return head + [str(total)] + list(verdict(counts))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf", required=True, help="the Stranger (or ExpansionHunter) VCF with the RFC1 record")
    ap.add_argument("--out", required=True)
    ap.add_argument("--bam")
    ap.add_argument("--skip", help="write the locus as not checked, with this reason")
    ap.add_argument("--no-softclip", action="store_true", help="do not add soft-clipped bases")
    ap.add_argument("--no-mates", action="store_true", help="do not fetch mates mapped elsewhere")
    a = ap.parse_args(argv)
    if not a.skip and not a.bam:
        ap.error("give --bam, or --skip REASON")
    try:
        row = check(a.bam, a.vcf, a.skip, not a.no_softclip, not a.no_mates)
    except (OSError, ValueError, EOFError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n" + "\t".join(row) + "\n")
    os.replace(tmp, a.out)
    print(f"RFC1 motif count: {row[9]} ({row[10]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
