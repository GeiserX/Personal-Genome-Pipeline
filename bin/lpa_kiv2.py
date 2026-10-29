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

The regions and the method are the LPA_KIV2 row of bin/cn_loci.tsv, run by
bin/copy_number.py (step 46 and the LPA processes call that engine); this
file lays out the step's table and keeps its own command line.

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
import sys

import copy_number as cn

LOCUS = "LPA_KIV2"
_ROW = []
ADVICE = ("A research estimate from read depth, not a diagnostic test. Lp(a) is measured in blood; a 2022 "
          "European Atherosclerosis Society consensus advises measuring it at least once in every adult.")


def row():
    """The LPA_KIV2 row of the catalog, read once."""
    if not _ROW:
        _ROW.append(cn.load_catalog()[LOCUS])
    return _ROW[0]


def __getattr__(name):
    # The regions and thresholds come from the catalog row, read on first use.
    r = row()
    if name == "REGIONS":   # in catalog order: the KIV-2 block, then the single-copy exons
        return [(r["chrom"], s, e, n) for n, s, e, _ in r["_targets"] + r["_controls"]]
    vals = {"MIN_FLANK_DEPTH": r["_min_control"], "MIN_MAPQ0_SHARE": r["_mapq0_min"],
            "LOW": r["_call"][1][0], "HIGH": r["_call"][1][1], "READINGS": r["_readings"]}
    if name in vals:
        return vals[name]
    raise AttributeError(name)


def bed_lines():
    return ["\t".join(map(str, r)) for r in cn.regions(row())]


def read_bedcov(path):
    """{name: mean depth} from `samtools bedcov` output over the regions,
    which must all be there and nothing else; {} for an empty table (no chr6)."""
    return cn.read_bedcov(path, row())


def table(b):
    """([(metric, value)], summary line) of the engine's LPA_KIV2 block."""
    r = row()
    out = {"status": "not assessed", "kiv2_copies": ".", "ci95_low": ".", "ci95_high": ".", "class": ".",
           "kiv2_depth": ".", "flank_depth": ".", "mapq0_share": ".", "message": "", "advice": ADVICE}
    reason = b["reason"]
    if reason == "no_data":
        out["message"] = "not assessed: the BAM has no chr6, so LPA has no reads"
    else:
        kdepth, fdepth, share = b["target_depth"], b["control_depth"], b["mapq0_share"]
        out.update(kiv2_depth=f"{kdepth:.2f}", flank_depth=f"{fdepth:.2f}", mapq0_share=f"{share:.3f}")
        if reason == "low_control_depth":
            out["status"] = "too little depth"
            out["message"] = (f"too little depth: {fdepth:.1f}x over the single-copy LPA exons, under the "
                              f"{r['_min_control']:g}x the count needs")
        elif reason == "low_mapq0_share":
            out["message"] = (f"not assessed: reads of mapping quality 0 make up {100 * share:.1f}% of the depth "
                              "over the KIV-2 exons, so the aligner left out the reads that fit several KIV-2 "
                              "copies (as novoalign does by default) and the count would come out low; realign "
                              "with minimap2 (step 02) or BWA-MEM")
        elif reason == "below_floor":
            out["message"] = (f"not assessed: the KIV-2 exons have {kdepth:.1f}x against {fdepth:.1f}x over the "
                              "single-copy exons, too little for any repeat count (a deletion or a mapping problem)")
        else:
            est, (lo, hi) = b["estimate"], b["ci95"]
            out.update(status="estimated", kiv2_copies=f"{est:.1f}", ci95_low=f"{lo:.1f}", ci95_high=f"{hi:.1f}",
                       **{"class": b["class"]},
                       message=(f"about {est:.0f} KIV-2 repeats in total (95% interval {lo:.1f} to {hi:.1f}): "
                                f"{b['reading']}"))
    line = (f"LPA KIV-2: {out['message']} (KIV-2 exons {out['kiv2_depth']}x, single-copy exons "
            f"{out['flank_depth']}x, MAPQ 0 share {out['mapq0_share']})")
    return [(k, out[k]) for k in ORDER], line


def assess(m, q1):
    """The result dict for the mean depths of read_bedcov: m over every read,
    q1 over reads of MAPQ 1 or more."""
    out = dict(table(cn.assess(row(), m, q1))[0])
    del out["advice"]
    return out


ORDER = ("status", "kiv2_copies", "ci95_low", "ci95_high", "class", "kiv2_depth", "flank_depth",
         "mapq0_share", "message", "advice")


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
    return cn.main(["call", "--locus", LOCUS, "--bedcov", args.bedcov, "--bedcov-mapq1", args.bedcov_mapq1,
                    "--table", args.out])


if __name__ == "__main__":
    sys.exit(main())
