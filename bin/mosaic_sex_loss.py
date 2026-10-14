#!/usr/bin/env python3
"""Mosaic loss of chromosome Y (men) or of one X (women) in the cells the DNA
came from, estimated from the depth step 16b's mosdepth wrote and from the
allele fractions of the VCF. Step 33 (and the Nextflow SAMPLE_QC process) runs
it after the verdict and adds its rows to the sample QC table.

Usage:
  mosaic_sex_loss.py --table S_sample_qc.tsv [--regions S.regions.bed.gz] [--vcf S.vcf.gz]

The sex is the one somalier infers from the reads (inferred_sex in the table),
else the declared one; with neither nothing is estimated.

Depth. mosdepth's 500 bp windows are grouped into 10 kb blocks, each block is
the median of its windows with reads, and the ratio is the median block of the
sex chromosome over the median block of chr1-22. Only the single-copy part of
the sex chromosome counts: GIAB's GRCh38 stratifications v3.3, chrY or chrX
outside the pseudoautosomal regions, the X-transposed region and the
ampliconic regions, without the centromere, pieces of 200 kb or more. That is
the normalised ratio; a whole-chromosome mean would count chrY's gaps, its
repeats and the copies of chrX it holds. A man with one Y in every cell reads
0.5, a woman with two X 1.0. With a fraction f of the cells missing that copy
the ratio drops by f/k of its expected value, k the normal number of copies
(1 for Y, 2 for X): f = k * (1 - ratio / expected).

Allele fractions (women only). At heterozygous chrX SNVs, a loss of one X in a
fraction f of the cells moves the allele fraction from 0.5 to 0.5 +- d, with
d = f / (2 (2 - f)), so f = 4 d / (1 + 2 d). Per site, u = (alt - ref)^2 / n
has mean phi + 4 (n - 1) d^2, phi the spread of the library at that depth.
phi is read from the heterozygous SNVs of chr1-22, so a noisy library does not
read as a shift: d^2 = (mean u on chrX - mean u on chr1-22) / (4 (mean n - 1)).

Uncertainty: the standard error of each median (1.2533 * 1.4826 * MAD / sqrt
of the number of blocks) and of each mean u, as a 95% range. The floor is
fixed: at 30x the depth of a sex chromosome against the autosomes varies by a
few percent between libraries (GC content, mappability), which is 5 to 10% of
the estimate, so a loss counts as detected only when the low end of its range
is above 10%.
"""
import argparse
import bisect
import gzip
import math
import os
import statistics
import sys
import zlib
from array import array

FLOOR = 0.10          # below this, nothing can be told at 30x
BLOCK = 10_000        # bp per block of mosdepth windows
MIN_WINDOWS = 10      # windows with reads a block needs (a block holds 20 at --by 500)
MIN_BLOCKS = 10       # blocks the sex chromosome and the autosomes each need
MIN_DP = 10           # reads a heterozygous site needs
MIN_X_HETS = 100      # heterozygous chrX sites the allele-fraction estimate needs
MIN_AUTO_HETS = 1000  # heterozygous autosomal sites it needs for the library's spread
Z95 = 1.96

# GRCh38, 0-based half-open. GIAB genome stratifications v3.3, GRCh38@all/XY:
# chrY_nonPAR (or chrX_nonPAR) minus _XTR and _ampliconic, minus the
# centromere (chrY:10,316,944-10,544,039, chrX:58,100,000-61,000,000), pieces of
# 200 kb or more.
SINGLE_COPY = {
    "Y": [(2781479, 3050044), (6748713, 7574481), (11037032, 11592902), (11674123, 13984473),
          (14058733, 15875093), (15905214, 16159795), (16425965, 17456266), (18870334, 20207793),
          (20257793, 21305953), (26415435, 26673214)],
    "X": [(2781479, 37099262), (37285837, 48343310), (48434604, 49157254), (49528394, 50228964),
          (50278964, 52032464), (52223402, 52489200), (52520537, 55437684), (55547739, 58100000),
          (61000000, 62412542), (62462542, 63185573), (63193371, 71674267), (71835832, 72721314),
          (73105236, 89140830), (93428068, 101563740), (101648990, 102180805), (102519463, 103940531),
          (104117650, 114281198), (114331198, 115738949), (115838949, 116557779), (116595566, 120879381),
          (120929381, 139870917), (139871081, 144425606), (144475606, 155701383)],
}
STARTS = {c: [s for s, _ in iv] for c, iv in SINGLE_COPY.items()}
AUTOSOMES = {str(i) for i in range(1, 23)}
# sex -> (chromosome, its expected depth ratio, its normal number of copies)
EXPECTED = {"male": ("Y", 0.5, 1), "female": ("X", 1.0, 2)}

KEYS = ["mosaic_sex", "mosaic_sex_from", "mosaic_chrom", "mosaic_depth_ratio", "mosaic_depth_ratio_expected",
        "mosaic_depth_blocks", "mosaic_loss_depth", "mosaic_loss_depth_low", "mosaic_loss_depth_high",
        "mosaic_x_het_sites", "mosaic_loss_af", "mosaic_loss_af_low", "mosaic_loss_af_high",
        "mosaic_floor", "mosaic_call", "mosaic_reason"]

READ_ERRORS = (OSError, ValueError, EOFError, zlib.error, UnicodeDecodeError)


def name(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def single_copy(chrom, start, end):
    """True when [start, end) lies inside one single-copy piece of chrom (X or Y)."""
    i = bisect.bisect_right(STARTS[chrom], start) - 1
    return i >= 0 and end <= SINGLE_COPY[chrom][i][1]


def open_text(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def median_se(xs):
    m = statistics.median(xs)
    mad = statistics.median(abs(x - m) for x in xs)
    return m, 1.2533 * 1.4826 * mad / math.sqrt(len(xs))


def depth_blocks(path, sex_chrom):
    """Block medians of mosdepth's windows with reads: (on chr1-22, on the
    single-copy part of sex_chrom)."""
    auto, sex = [], []
    key, vals = None, []

    def flush():
        if key is not None and len(vals) >= MIN_WINDOWS:
            (sex if key[0] else auto).append(statistics.median(vals))

    with open_text(path) as f:
        for line in f:
            c, s, e, d = line.split("\t", 4)[:4]
            n = name(c)
            if n in AUTOSOMES:
                is_sex = False
            elif n == sex_chrom and single_copy(n, int(s), int(e)):
                is_sex = True
            else:
                continue
            dp = float(d)
            if dp <= 0:   # a gap, or outside what a sliced or targeted BAM holds
                continue
            k = (is_sex, c, int(s) // BLOCK)
            if k != key:
                flush()
                key, vals = k, []
            vals.append(dp)
    flush()
    return auto, sex


def depth_estimate(path, sex):
    """(ratio, blocks on the sex chromosome, loss, low, high), or a reason string."""
    chrom, expected, copies = EXPECTED[sex]
    auto, sx = depth_blocks(path, chrom)
    if len(auto) < MIN_BLOCKS or len(sx) < MIN_BLOCKS:
        return (f"too few 10 kb blocks with reads: {len(sx)} on chr{chrom}'s single-copy part and "
                f"{len(auto)} on chr1-22 (each needs {MIN_BLOCKS}; a targeted or sliced BAM may have none)")
    ma, sea = median_se(auto)
    ms, ses = median_se(sx)
    r = ms / ma
    se_r = r * math.sqrt((ses / ms) ** 2 + (sea / ma) ** 2) if ms > 0 else 0.0
    f = copies * (1 - r / expected)
    se_f = copies * se_r / expected
    return r, len(sx), f, f - Z95 * se_f, f + Z95 * se_f


def het_sites(path):
    """(n, u) arrays of the PASS heterozygous SNVs of chr1-22 and of the
    single-copy part of chrX, from the first sample's GT and AD."""
    an, au, xn, xu = array("I"), array("d"), array("I"), array("d")
    with open_text(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t", 10)
            if len(p) < 10:
                continue
            n = name(p[0])
            if n in AUTOSOMES:
                to_x = False
            elif n == "X":
                pos = int(p[1])
                if not single_copy("X", pos - 1, pos):
                    continue
                to_x = True
            else:
                continue
            if p[6] not in ("PASS", ".") or len(p[3]) != 1 or len(p[4]) != 1 or p[4] in ".*":
                continue
            keys = p[8].split(":")
            vals = p[9].rstrip("\n").split(":")
            if "GT" not in keys or "AD" not in keys or len(vals) < len(keys):
                continue
            if vals[keys.index("GT")].replace("|", "/") not in ("0/1", "1/0"):
                continue
            ad = vals[keys.index("AD")].split(",")
            try:
                ref, alt = int(ad[0]), int(ad[1])
            except (ValueError, IndexError):
                continue
            dp = ref + alt
            if dp < MIN_DP:
                continue
            (xn if to_x else an).append(dp)
            (xu if to_x else au).append((alt - ref) ** 2 / dp)
    return an, au, xn, xu


def af_estimate(path):
    """(chrX sites, loss, low, high), or a reason string."""
    an, au, xn, xu = het_sites(path)
    if len(an) < MIN_AUTO_HETS:
        return f"{len(an)} heterozygous SNVs on chr1-22 with {MIN_DP} or more reads (needs {MIN_AUTO_HETS})"
    # Sites far deeper than the genome are collapsed copies; their allele
    # fractions are not one locus'.
    cap = 2 * statistics.median(an)
    a = [u for n, u in zip(an, au) if n <= cap]
    x = [(n, u) for n, u in zip(xn, xu) if n <= cap]
    if len(x) < MIN_X_HETS:
        return f"{len(x)} heterozygous SNVs on chrX's single-copy part (needs {MIN_X_HETS})"
    xs = [u for _, u in x]
    scale = 4 * (statistics.fmean(n for n, _ in x) - 1)
    d2 = (statistics.fmean(xs) - statistics.fmean(a)) / scale
    se = math.sqrt(statistics.variance(xs) / len(xs) + statistics.variance(a) / len(a)) / scale

    def loss(v):
        d = math.sqrt(max(v, 0.0))
        return 4 * d / (1 + 2 * d)
    return len(x), loss(d2), loss(d2 - Z95 * se), loss(d2 + Z95 * se)


def fmt(v):
    return f"{v:.3f}"


def estimate(table, regions=None, vcf=None):
    """The mosaic rows for one sample, as ordered (key, value) pairs."""
    q = dict(table)
    sex, source = q.get("inferred_sex"), "reads"
    if sex not in EXPECTED:
        sex, source = q.get("declared_sex"), "declared"
    row = dict.fromkeys(KEYS, "")
    row["mosaic_floor"] = fmt(FLOOR)
    if sex not in EXPECTED:
        row.update(mosaic_call="not_estimated",
                   mosaic_reason="no sex: none declared and somalier could not tell it from the reads")
        return list(row.items())
    chrom, expected, _ = EXPECTED[sex]
    row.update(mosaic_sex=sex, mosaic_sex_from=source, mosaic_chrom=f"chr{chrom}",
               mosaic_depth_ratio_expected=f"{expected:.1f}")
    why = []
    dep_lo = dep_hi = af_lo = None
    if not regions:
        why.append("no mosdepth windows (step 16b's <sample>.regions.bed.gz): depth not read")
    else:
        try:
            dep = depth_estimate(regions, sex)
        except READ_ERRORS as e:
            dep = f"could not read {os.path.basename(regions)}: {type(e).__name__}: {e}"
        if isinstance(dep, str):
            why.append(f"depth not estimated: {dep}")
        else:
            r, blocks, f, dep_lo, dep_hi = dep
            row.update(mosaic_depth_ratio=fmt(r), mosaic_depth_blocks=str(blocks), mosaic_loss_depth=fmt(f),
                       mosaic_loss_depth_low=fmt(dep_lo), mosaic_loss_depth_high=fmt(dep_hi))
    if sex == "female":
        if not vcf:
            why.append("no VCF: chrX allele fractions not read")
        else:
            try:
                af = af_estimate(vcf)
            except READ_ERRORS as e:
                af = f"could not read {os.path.basename(vcf)}: {type(e).__name__}: {e}"
            if isinstance(af, str):
                why.append(f"allele fractions not estimated: {af}")
            else:
                n, f, af_lo, hi = af
                row.update(mosaic_x_het_sites=str(n), mosaic_loss_af=fmt(f), mosaic_loss_af_low=fmt(af_lo),
                           mosaic_loss_af_high=fmt(hi))
    what = "the Y" if sex == "male" else "one X"
    if dep_lo is None and af_lo is None:
        call = "not_estimated"
    elif dep_hi is not None and dep_hi < -FLOOR:
        # An extra copy also shifts chrX allele fractions (1/3 and 2/3 in 47,XXX).
        call = "gain"
        why.insert(0, f"more chr{chrom} depth than {'one copy' if sex == 'male' else 'two copies'}: "
                      "an extra copy (47,XYY or 47,XXX) or a gain, not a loss")
    elif dep_lo is not None and dep_lo > FLOOR:
        call = "detected"
        why.insert(0, f"the low end of the 95% range is above {FLOOR:.0%}: a share of the cells lost {what}")
    elif af_lo is not None and af_lo > FLOOR:
        call = "detected"
        why.insert(0, f"chrX allele fractions are shifted beyond {FLOOR:.0%} but its depth is not: in a share of "
                      "the cells one X was lost and the other copied (copy-neutral), or the depth cannot tell")
    else:
        call = "not_detected"
        why.insert(0, f"no loss of {what} above {FLOOR:.0%}; below that nothing can be told at 30x")
    row.update(mosaic_call=call, mosaic_reason="; ".join(why))
    return list(row.items())


def read_table(path):
    with open(path) as f:
        lines = [l.rstrip("\n").split("\t", 1) for l in f if l.strip()]
    if not lines or lines[0] != ["key", "value"]:
        raise ValueError(f"{path} is not a key/value table")
    return [(r[0], r[1] if len(r) > 1 else "") for r in lines[1:]]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", required=True, help="the step 33 table; the rows are added to it")
    ap.add_argument("--regions", help="mosdepth's <sample>.regions.bed.gz (step 16b, --by 500)")
    ap.add_argument("--vcf", help="the sample's VCF (GT and AD), for the chrX allele fractions of a woman")
    a = ap.parse_args(argv)
    try:
        table = read_table(a.table)
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    # A rerun replaces the rows of the last one.
    kept = [(k, v) for k, v in table if not k.startswith("mosaic_")]
    rows = estimate(kept, a.regions, a.vcf)
    tmp = a.table + ".tmp"
    with open(tmp, "w") as f:
        f.write("key\tvalue\n")
        for k, v in kept + rows:
            f.write(f"{k}\t{v}\n")
    os.replace(tmp, a.table)
    for k, v in rows:
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
