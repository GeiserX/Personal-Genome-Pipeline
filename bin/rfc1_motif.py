#!/usr/bin/env python3
"""rfc1_motif.py: which repeat motif do the reads at RFC1 (or another motif locus) carry?

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
(mate rescue). It finds the longest run of one pentamer repeated in
tandem in those bases, and when the run is longer than the reference
allele (more than its 11 units, and at least MIN_UNITS) it counts the read
for that motif. A read across an unexpanded allele is not counted, so the
reads of a reference allele cannot outvote those of an expanded one.
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

  rfc1_motif.py --bam BAM --vcf STRANGER_VCF --out OUT.tsv [--catalog CATALOG.json]
  rfc1_motif.py --vcf STRANGER_VCF --skip REASON --out OUT.tsv
  rfc1_motif.py --bam BAM --vcf STRANGER_VCF --catalog CATALOG.json --flagged --out OUT.tsv

The locus is the RFC1 record of the VCF (REPID=RFC1): its POS to its END.
With --catalog (the STRchive catalog bin/strchive_convert.py writes, which
ExpansionHunter and Stranger ran on) the motifs come from the locus's
BenignMotifs and PathogenicMotifs, the repeat unit is their length, and the
locus is every repeat of the catalog entry, start of the first to end of
the last (at RFC1 ExpansionHunter's RFC1 record is then the AAGGG part
alone). Without it RFC1 is read with the motifs above. --flagged reads
every other locus of the catalog that has motifs and that Stranger flags
(pre_mutation or full_mutation), one row each, the FAME loci (TTTCA
pathogenic, TTTTA benign; Ishiura et al., Nat Genet 2018;50:581-590) and
SCA31 and SCA37 among them, and writes no file when none is flagged.
OUT.tsv has a header and one row per locus:

  locus chrom start end repeat_reads dominant_motif dominant_class
  dominant_share motif_shares verdict note

motif_shares lists every motif seen as MOTIF=share, most common first.
--no-softclip and --no-mates turn the two rescues off (for the tests).
Needs pysam (run it in PYPGX_IMAGE); the counting itself is plain Python.
"""
import argparse
import json
import os
import re
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


def motif_name(m, known=None):
    """The known name of a motif, or its smallest rotation on either strand."""
    table = FORM_TO_KNOWN if known is None else {f: k for k in known for f in forms(k)}
    return table.get(m) or min(forms(m))


def longest_run(seq, period=5, known=None):
    """(motif, units) of the longest tandem run of one unit of period bases
    in seq, or (None, 0). A homopolymer run is not a run of the unit."""
    seq = seq.upper()
    best, best_len, i, n = None, 0, 0, len(seq)
    while i + period < n:
        j = i
        while j + period < n and seq[j] == seq[j + period] and seq[j] != "N":
            j += 1
        length = j - i + period
        unit = seq[i:i + period]
        if j > i and length > best_len and len(set(unit)) > 1:
            best, best_len = unit, length
        i = j + 1 if j > i else i + 1
    return (motif_name(best, known), best_len // period) if best else (None, 0)


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
    if softclip:
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


def min_units(start, end, period=5):
    """Units a run needs to count: more than the reference allele of [start, end)."""
    return max(MIN_UNITS, (end - start) // period + 1)


def count_motifs(segments, units_min=MIN_UNITS, period=5, known=None):
    counts = {}
    for _, bases in segments:
        motif, units = longest_run(bases, period, known)
        if motif and units >= units_min:
            counts[motif] = counts.get(motif, 0) + 1
    return counts


# RFC1 as read without a catalog: the motifs above, CANVAS, recessive.
RFC1 = {"locus": "RFC1", "known": KNOWN, "period": 5, "disease": "CANVAS", "recessive": True, "span": None}


def verdict(counts, units_min=MIN_UNITS, spec=None):
    """(dominant_motif, dominant_class, dominant_share, shares, verdict, note)."""
    spec = spec or RFC1
    known, locus = spec["known"], spec["locus"]
    total = sum(counts.values())
    order = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    shares = ",".join(f"{m}={c / total:.2f}" for m, c in order) if total else "."
    if total < MIN_READS:
        return (".", ".", ".", shares, "not_resolved",
                f"{total} read(s) carry {units_min} or more repeat units at {locus}, fewer than {MIN_READS}: "
                "the motif is not resolved; confirm the flag with repeat-primed PCR if symptoms fit")
    dom, n = order[0]
    cls = known.get(dom, "unknown")
    share = n / total
    path = {m: c for m, c in counts.items() if known.get(m) == "pathogenic"}
    npath = sum(path.values())
    if npath >= MIN_READS and npath / total >= PATHOGENIC_SHARE:
        names = "/".join(sorted(path, key=lambda m: -path[m]))
        both = f" ({spec['disease']} needs a pathogenic expansion on both copies)" if spec["recessive"] else ""
        return (dom, cls, f"{share:.2f}", shares, "pathogenic_motif",
                f"{names} in {npath} of {total} repeat reads ({npath / total:.0%}): confirm with repeat-primed PCR"
                + both)
    if cls == "benign" and share >= BENIGN_SHARE:
        return (dom, cls, f"{share:.2f}", shares, "benign_motif",
                f"motif {dom} (benign) in {n} of {total} repeat reads ({share:.0%}): the flag is not supported")
    return (dom, cls, f"{share:.2f}", shares, "uncertain_motif",
            f"most common motif {dom} ({cls}) in {n} of {total} repeat reads ({share:.0%}): "
            "not a known benign expansion; confirm with repeat-primed PCR if symptoms fit")


def vcf_records(vcf):
    """REPID -> (chrom, start, end, STR_STATUS) of every record: 0-based start, end exclusive."""
    import gzip
    out = {}
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = dict(x.split("=", 1) if "=" in x else (x, "") for x in f[7].split(";"))
            start = int(f[1]) - 1
            end = int(info["END"]) if info.get("END", "").isdigit() else start + len(f[3])
            for key in (info.get("REPID"), info.get("VARID")):
                if key and key not in out:
                    out[key] = (f[0], start, end, info.get("STR_STATUS", ""))
    return out


def rfc1_locus(vcf):
    """(chrom, start, end) of the REPID=RFC1 record: 0-based start, end exclusive."""
    r = vcf_records(vcf).get("RFC1")
    return r[:3] if r else None


REGION = re.compile(r"^(chr[^:]+|[^:]+):(\d+)-(\d+)$")


def catalog_specs(path):
    """LocusId -> how to read it, for every catalog entry with motifs."""
    with open(path) as f:
        cat = json.load(f)
    out = {}
    for e in cat:
        if not e.get("PathogenicMotifs"):
            continue
        known = {m: "benign" for m in e.get("BenignMotifs", [])}
        known.update({m: "pathogenic" for m in e["PathogenicMotifs"]})
        regions = [REGION.match(r).groups() for r in
                   (e["ReferenceRegion"] if isinstance(e["ReferenceRegion"], list) else [e["ReferenceRegion"]])]
        span = (regions[0][0], min(int(r[1]) for r in regions), max(int(r[2]) for r in regions))
        out[e["LocusId"]] = {"locus": e["LocusId"], "known": known, "period": len(e["PathogenicMotifs"][0]),
                             "disease": e.get("Disease", e["LocusId"]),
                             "recessive": e.get("InheritanceMode", "").split("/") == ["AR"], "span": span}
    return out


def flagged(status):
    return any(x in ("pre_mutation", "full_mutation") for x in status.split(","))


def check(bam_path, vcf, skip=None, softclip=True, mates=True, spec=None, records=None):
    spec = spec or RFC1
    name = spec["locus"]
    rec = (records if records is not None else vcf_records(vcf)).get(name)
    if rec is None:
        return [name, ".", ".", ".", "0", ".", ".", ".", ".", "not_checked", f"no {name} record in the VCF"]
    chrom, start, end = spec["span"] or rec[:3]
    head = [name, chrom, str(start + 1), str(end)]
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
        units = min_units(start, end, spec["period"])
        try:
            counts = count_motifs(read_segments(bam.fetch(chrom, flo, fhi), lo, hi,
                                                mate_of if mates else None, softclip, (flo, fhi)),
                                  units, spec["period"], spec["known"])
        finally:
            other.close()
    total = sum(counts.values())
    return head + [str(total)] + list(verdict(counts, units, spec))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf", required=True, help="the Stranger (or ExpansionHunter) VCF with the RFC1 record")
    ap.add_argument("--out", required=True)
    ap.add_argument("--bam")
    ap.add_argument("--skip", help="write the locus as not checked, with this reason")
    ap.add_argument("--no-softclip", action="store_true", help="do not add soft-clipped bases")
    ap.add_argument("--no-mates", action="store_true", help="do not fetch mates mapped elsewhere")
    ap.add_argument("--catalog", help="the STRchive catalog (bin/strchive_convert.py): motifs and locus from it")
    ap.add_argument("--flagged", action="store_true",
                    help="every flagged locus of --catalog with motifs other than RFC1, one row each")
    a = ap.parse_args(argv)
    if not a.skip and not a.bam:
        ap.error("give --bam, or --skip REASON")
    if a.flagged and not a.catalog:
        ap.error("--flagged needs --catalog")
    try:
        specs = catalog_specs(a.catalog) if a.catalog else {}
        records = vcf_records(a.vcf)
        if a.flagged:
            todo = [sp for name, sp in specs.items()
                    if name != "RFC1" and name in records and flagged(records[name][3])]
        else:
            todo = [specs.get("RFC1", RFC1)]
        rows = [check(a.bam, a.vcf, a.skip, not a.no_softclip, not a.no_mates, sp, records) for sp in todo]
    except (OSError, ValueError, EOFError, KeyError, AttributeError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    if not rows:
        if os.path.exists(a.out):
            os.remove(a.out)
        print("No flagged locus with motifs in the catalog: no motif count")
        return 0
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    os.replace(tmp, a.out)
    for row in rows:
        print(f"{row[0]} motif count: {row[9]} ({row[10]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
