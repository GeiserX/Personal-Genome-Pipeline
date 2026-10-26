#!/usr/bin/env python3
"""lpa_kiv2.py: the LPA KIV-2 repeat count from short-read depth (step 46).

Lp(a) is a mostly inherited cardiovascular risk factor, and the number of
kringle IV type 2 (KIV-2) repeats in LPA sets the size of apo(a): fewer
repeats, a smaller isoform and, on average, a higher Lp(a) level. The KIV-2
repeat is a 5.5 kb unit with two exons; a Danish population study counted 6
to 99 on both chromosomes together (Kamstrup et al. 2009, PMID 19509380).
GRCh38 holds six units, and their exons are near-identical, so an aligner
that reports every read places each KIV-2 read on one of the six (with
mapping quality 0). The depth over those exons, against the depth over
single-copy LPA exons, counts the units the sample carries: the approach of
the published short-read estimators (Zekavat et al. 2018, PMID 29973585;
Behera et al. 2024, PMID 39449055, which also gives the GRCh38 region and the
six reference units).

  lpa_kiv2.py bed
      print the regions (GRCh38, 0-based BED, name in column 4) to stdout
  lpa_kiv2.py estimate --bedcov BEDCOV.tsv --bedcov-mapq1 BEDCOV_Q1.tsv --out LPA.tsv
      read `samtools bedcov -j` over those regions (chrom, start, end, name,
      summed depth; all MAPQ, duplicates left out) and `samtools bedcov -j -Q 1`
      over the same regions, and write LPA.tsv (metric<TAB>value). Exit 0
      whatever the status; 2 on unreadable input. Two empty tables mean the
      BAM has no chr6: status 'not assessed'.

The regions are exons of the MANE Select transcript ENST00000316300.10,
numbered from its first exon:
  kiv2_e03 to kiv2_e16  the six KIV-2 units of GRCh38 (exons 4 to 15), plus
                        exon 3 (the last exon of KIV-1) and exon 16 (the first
                        of KIV-3), which are near-identical to KIV-2 exons and
                        take some of its reads. Seven 160 bp and seven 182 bp
                        exons, 2,394 bp: six KIV-2 units plus one single-copy
                        pair per chromosome.
  flank_e19 to flank_e33  single-copy kringle exons downstream (the second
                        exon of KIV-4 to KV), under 90% identical to any KIV-2
                        exon. Exons 17 and 18 (about 95% identical) are in
                        neither group, nor are exons 1 and 2 or the protease
                        exons, which resemble plasminogen (PLG).

With D the mean depth over the flank exons (two copies) and M the mean depth
over the KIV-2 block, the block holds N + 2 copies of a 342 bp pair spread
over seven reference pairs, so M = (N + 2) x D / 14 and
    N = 14 x M / D - 2      (the total over both chromosomes)
The 95% interval comes from the spread of the per-exon depths in each group:
relative standard error sqrt(cv_kiv2^2 / 14 + cv_flank^2 / 15), applied to
N + 2 on the log scale. It covers sampling noise only, not a systematic bias
of an aligner or a library.

status is 'estimated', 'too little depth' or 'not assessed'; the message says
why. Research estimate, not a diagnostic test: the actionable test is a blood
Lp(a) measurement. Standard library only: it runs in the Python image and
under tests/test_lpa_kiv2.py.
"""
import argparse
import math
import sys

# GRCh38, 0-based, from Ensembl's ENST00000316300.10 (checked 2026-10-25).
REGIONS = [
    ("chr6", 160646213, 160646395, "kiv2_e03"),
    ("chr6", 160644846, 160645006, "kiv2_e04"),
    ("chr6", 160640666, 160640848, "kiv2_e05"),
    ("chr6", 160639299, 160639459, "kiv2_e06"),
    ("chr6", 160635122, 160635304, "kiv2_e07"),
    ("chr6", 160633752, 160633912, "kiv2_e08"),
    ("chr6", 160629573, 160629755, "kiv2_e09"),
    ("chr6", 160628206, 160628366, "kiv2_e10"),
    ("chr6", 160624029, 160624211, "kiv2_e11"),
    ("chr6", 160622662, 160622822, "kiv2_e12"),
    ("chr6", 160618483, 160618665, "kiv2_e13"),
    ("chr6", 160617116, 160617276, "kiv2_e14"),
    ("chr6", 160612931, 160613113, "kiv2_e15"),
    ("chr6", 160611561, 160611721, "kiv2_e16"),
    ("chr6", 160600916, 160601098, "flank_e19"),
    ("chr6", 160599499, 160599659, "flank_e20"),
    ("chr6", 160595353, 160595535, "flank_e21"),
    ("chr6", 160593957, 160594117, "flank_e22"),
    ("chr6", 160590943, 160591101, "flank_e23"),
    ("chr6", 160589552, 160589712, "flank_e24"),
    ("chr6", 160586448, 160586630, "flank_e25"),
    ("chr6", 160585045, 160585205, "flank_e26"),
    ("chr6", 160578522, 160578704, "flank_e27"),
    ("chr6", 160577135, 160577295, "flank_e28"),
    ("chr6", 160557389, 160557571, "flank_e29"),
    ("chr6", 160556024, 160556184, "flank_e30"),
    ("chr6", 160548477, 160548659, "flank_e31"),
    ("chr6", 160547788, 160547937, "flank_e32"),
    ("chr6", 160545439, 160545533, "flank_e33"),
]
# Copies of the 342 bp exon pair in the KIV-2 block of GRCh38: six KIV-2
# units and the KIV-1/KIV-3 pair. One single-copy pair per chromosome is
# taken off the total.
BLOCK_PAIRS = 7
SINGLE_COPY_PAIRS = 2
# Below this mean flank depth the ratio is noise.
MIN_FLANK_DEPTH = 10.0
# Part of the KIV-2 block depth from reads of MAPQ 0 under which the BAM is
# taken to leave out the reads that fit several KIV-2 copies (novoalign does
# by default): the count would come out low.
MIN_MAPQ0_SHARE = 0.2
# About 70% of the 1000 Genomes samples had a total of 30 to 45 KIV-2 copies
# (Behera et al. 2024, DRAGEN's caller; mean 37.5 in Europeans, 35.7 in
# Africans, 44.6 in East Asians).
LOW, HIGH = 30, 45
Z95 = 1.96

READINGS = {
    "low": (f"fewer KIV-2 repeats than most people (the whole 95% interval is under {LOW}; about 70% of "
            f"the 1000 Genomes samples had {LOW} to {HIGH}): a low count goes with small apo(a) isoforms "
            "and a higher Lp(a) level"),
    "typical": (f"within or overlapping the range most people have ({LOW} to {HIGH} for about 70% of the "
                "1000 Genomes samples)"),
    "high": (f"more KIV-2 repeats than most people (the whole 95% interval is over {HIGH}): a high count "
             "goes with large apo(a) isoforms and, on average, a lower Lp(a) level"),
}
ADVICE = ("A research estimate from read depth, not a diagnostic test. Lp(a) is measured in blood; a 2022 "
          "European Atherosclerosis Society consensus advises measuring it at least once in every adult.")


def bed_lines():
    return ["\t".join(map(str, r)) for r in REGIONS]


def read_bedcov(path):
    """{name: mean depth} from `samtools bedcov` output over REGIONS, which
    must all be there and nothing else; {} for an empty table (no chr6)."""
    means, seen = {}, set()
    with open(path, errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                raise ValueError(f"{path}: expected chrom, start, end, name, depth sum; got {line.strip()!r}")
            chrom, start, end, name = c[0], int(c[1]), int(c[2]), c[3]
            seen.add((chrom, start, end, name))
            means[name] = float(c[4]) / (end - start)
    if seen and seen != set(REGIONS):
        raise ValueError(f"{path}: its regions differ from the ones `lpa_kiv2.py bed` prints")
    return means


def _cv(values):
    """Coefficient of variation (sample standard deviation over the mean)."""
    m = sum(values) / len(values)
    if m <= 0 or len(values) < 2:
        return 0.0
    var = sum((v - m) ** 2 for v in values) / (len(values) - 1)
    return math.sqrt(var) / m


def assess(m, q1):
    """The result dict for the mean depths of read_bedcov: m over every read,
    q1 over reads of MAPQ 1 or more."""
    out = {"status": "not assessed", "kiv2_copies": ".", "ci95_low": ".", "ci95_high": ".", "class": ".",
           "kiv2_depth": ".", "flank_depth": ".", "mapq0_share": ".", "message": ""}
    if not m:
        out["message"] = "not assessed: the BAM has no chr6, so LPA has no reads"
        return out
    kiv = [r for r in REGIONS if r[3].startswith("kiv2_")]
    fl = [r for r in REGIONS if r[3].startswith("flank_")]
    klen = sum(e - s for _, s, e, _ in kiv)
    flen = sum(e - s for _, s, e, _ in fl)
    kdepth = sum(m[n] * (e - s) for _, s, e, n in kiv) / klen
    fdepth = sum(m[n] * (e - s) for _, s, e, n in fl) / flen
    kq1 = sum(q1.get(n, 0.0) * (e - s) for _, s, e, n in kiv) / klen
    share = 1 - kq1 / kdepth if kdepth > 0 else 0.0
    out.update(kiv2_depth=f"{kdepth:.2f}", flank_depth=f"{fdepth:.2f}", mapq0_share=f"{share:.3f}")
    if fdepth < MIN_FLANK_DEPTH:
        out["status"] = "too little depth"
        out["message"] = (f"too little depth: {fdepth:.1f}x over the single-copy LPA exons, under the "
                          f"{MIN_FLANK_DEPTH:g}x the count needs")
        return out
    if share < MIN_MAPQ0_SHARE:
        out["message"] = (f"not assessed: reads of mapping quality 0 make up {100 * share:.1f}% of the depth over "
                          "the KIV-2 exons, so the aligner left out the reads that fit several KIV-2 copies (as "
                          "novoalign does by default) and the count would come out low; realign with minimap2 "
                          "(step 02) or BWA-MEM")
        return out
    total = BLOCK_PAIRS * 2 * kdepth / fdepth          # N + 2
    est = total - SINGLE_COPY_PAIRS
    if est <= 0:
        out["message"] = (f"not assessed: the KIV-2 exons have {kdepth:.1f}x against {fdepth:.1f}x over the "
                          "single-copy exons, too little for any repeat count (a deletion or a mapping problem)")
        return out
    rel = math.sqrt(_cv([m[n] for *_, n in kiv]) ** 2 / len(kiv) + _cv([m[n] for *_, n in fl]) ** 2 / len(fl))
    lo = max(0.0, total * math.exp(-Z95 * rel) - SINGLE_COPY_PAIRS)
    hi = total * math.exp(Z95 * rel) - SINGLE_COPY_PAIRS
    cls = "low" if hi < LOW else "high" if lo > HIGH else "typical"
    out.update(status="estimated", kiv2_copies=f"{est:.1f}", ci95_low=f"{lo:.1f}", ci95_high=f"{hi:.1f}",
               **{"class": cls},
               message=f"about {est:.0f} KIV-2 repeats in total (95% interval {lo:.1f} to {hi:.1f}): {READINGS[cls]}")
    return out


ORDER = ("status", "kiv2_copies", "ci95_low", "ci95_high", "class", "kiv2_depth", "flank_depth",
         "mapq0_share", "message", "advice")


def estimate(args):
    try:
        m, q1 = read_bedcov(args.bedcov), read_bedcov(args.bedcov_mapq1)
        if bool(m) != bool(q1):
            raise ValueError(f"{args.bedcov_mapq1 if m else args.bedcov} is empty but the other table is not: "
                             "both must cover the regions, or both be empty (no chr6)")
        res = assess(m, q1)
    except (OSError, ValueError, KeyError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res["advice"] = ADVICE
    try:
        with open(args.out, "w") as f:
            f.write("metric\tvalue\n")
            for k in ORDER:
                f.write(f"{k}\t{res[k]}\n")
    except OSError as e:
        print(f"ERROR: cannot write {args.out}: {e}", file=sys.stderr)
        return 2
    print(f"LPA KIV-2: {res['message']} (KIV-2 exons {res['kiv2_depth']}x, single-copy exons "
          f"{res['flank_depth']}x, MAPQ 0 share {res['mapq0_share']})")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="lpa_kiv2.py", description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bed", help="print the regions")
    c = sub.add_parser("estimate", help="count the KIV-2 repeats from samtools bedcov over the regions")
    c.add_argument("--bedcov", required=True)
    c.add_argument("--bedcov-mapq1", required=True, help="samtools bedcov -j -Q 1 over the same regions")
    c.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "bed":
        print("\n".join(bed_lines()))
        return 0
    return estimate(args)


if __name__ == "__main__":
    sys.exit(main())
