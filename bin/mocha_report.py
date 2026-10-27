#!/usr/bin/env python3
"""mocha_report.py: the mosaic chromosomal alterations of step 51 from MoChA's tables.

A mosaic chromosomal alteration (mCA) is a large gain, loss or copy-neutral
loss of heterozygosity (CN-LOH) carried by a share of the blood cells, not by
every cell. In a population study of 151,202 UK Biobank participants they
were more common with age and went with a higher risk of a later blood cancer
(Loh et al. 2018, PMID 29995854). MoChA (freeseek/mocha) calls them from a
phased VCF: in the cells that carry the event one parental haplotype is lost or
gained, so the allele balance at every heterozygous site leans the same way
along the phased haplotype, and the depth moves with a gain or a loss.

  mocha_report.py --sample S --calls CALLS.tsv --stats STATS.tsv --sites SITES.tsv \\
      --phasing beagle|input --out-calls S_mocha.tsv --out-summary S_mocha_summary.tsv

CALLS.tsv and STATS.tsv are MoChA's --calls and --stats tables. SITES.tsv has
one row per chromosome MoChA read: chrom, first and last position, sites,
heterozygous sites and phased heterozygous sites (no header). Exit 0 whatever
the result; 2 on unreadable input.

S_mocha.tsv has one row per MoChA call, `reported` yes or no and the reason.
A call is reported when all of these hold:
  - it is on chr1 to chr22 (mosaic loss of Y and of X are step 33's);
  - MoChA's own filter of the UK Biobank callset keeps it (README of
    freeseek/mocha, "Filter callset"): not a germline copy-number polymorphism
    (type CNP_*, which MoChA gives only to a call over a region of its --cnp
    list, the panel's common germline CNVs of setup.sh --mocha-data); a BAF deviation under 0.1, or heterozygous sites closer than
    200 kb (N50), or a phase-concordance LOD over 10, so a constitutional
    event is left out; a BAF deviation standard error, or a BAF-and-phase LOD
    over 10; and not a likely germline duplication (relative coverage 2.1 or
    more with a large BAF deviation, with the size limits of that filter);
  - MoChA estimated its cell fraction and it is at or above the floor.

The floor is the larger of FLOOR_MIN (0.05, our choice: below it the
pipeline reports nothing, whatever the depth) and the read-sampling floor at
the sample's depth: the cell fraction at which a 10 Mb gain or loss would move
the allele balance by three standard errors of read sampling alone,
6 / sqrt(D x H), with D MoChA's median depth and H the heterozygous sites per
10 Mb over the chromosomes read. A loss or gain of one copy in a fraction c of
the cells moves the allele balance by about c / 4, and the standard error of
the mean balance over H sites read D times is 0.5 / sqrt(D x H). CN-LOH moves
it by c / 2, so it shows at half that floor. Phasing errors and bias in the
reads are not in it: it is a lower bound, never a promise that an event above
it is found.

The sample is not assessed when MoChA's sample filter of the same README
fails: a call rate under 0.97, or a phased BAF autocorrelation over 0.03.

S_mocha_summary.tsv is metric<TAB>value: status (assessed or not assessed),
message, depth, het_sites, phased_het_sites, het_per_10mb, analysed_mb,
chromosomes, call_rate, baf_auto, floor_read, floor_min, floor, calls_mocha,
calls_reported, phasing.

Standard library only: it runs in the Python image and under
tests/test_mocha_report.py.
"""
import argparse
import csv
import math
import sys

FLOOR_MIN = 0.05
FLOOR_Z = 3.0
FLOOR_MB = 10.0
CALL_RATE_MIN = 0.97
BAF_AUTO_MAX = 0.03
AUTOSOMES = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}
TYPES = {"Loss": "loss", "Gain": "gain", "CN-LOH": "CN-LOH", "Undetermined": "undetermined"}
CALL_COLUMNS = ("chrom", "start", "end", "length_mb", "arms", "type", "cell_fraction", "bdev", "rel_cov",
                "n_hets", "lod_baf_phase", "lod_lrr_baf", "reported", "reason")
NEED_CALLS = ("chrom", "length", "p_arm", "q_arm", "n50_hets", "bdev", "bdev_se", "rel_cov", "lod_baf_phase",
              "lod_lrr_baf", "lod_baf_conc", "n_hets", "type", "cf")
NEED_STATS = ("sample_id", "call_rate", "cov_median", "baf_auto", "n_hets")


class Unreadable(Exception):
    pass


def num(x):
    """MoChA's number, or nan for 'nan', '-nan', '.' or an empty cell."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return math.nan
    return v


def read_table(path, need):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = reader.fieldnames or []
        rows = list(reader)
    missing = [k for k in need if k not in header]
    if missing:
        raise Unreadable(f"{path} lacks the column{'s' if len(missing) > 1 else ''} {', '.join(missing)}")
    return rows


def read_sites(path):
    out = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for n, line in enumerate(f, 1):
            c = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            if len(c) != 6 or not all(x.isdigit() for x in c[1:]):
                raise Unreadable(f"{path} line {n}: want chrom, first, last, sites, hets, phased hets")
            out.append((c[0], int(c[1]), int(c[2]), int(c[3]), int(c[4]), int(c[5])))
    return out


def floor_read(depth, het_per_10mb):
    """The read-sampling floor: 6 / sqrt(D x H), the cell fraction of a one-copy
    gain or loss at FLOOR_Z standard errors over FLOOR_MB; nan without data."""
    if not (depth > 0 and het_per_10mb > 0):
        return math.nan
    return 4 * FLOOR_Z * 0.5 / math.sqrt(depth * het_per_10mb)


def arms(chrom, p_arm, q_arm):
    """'20p', '20q' or '20p and 20q', with the telomere it reaches."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    on = [a for a, v in (("p", p_arm), ("q", q_arm)) if v in ("Y", "C", "T")]
    tel = [a for a, v in (("p", p_arm), ("q", q_arm)) if v == "T"]
    if not on:
        return f"{c} (interstitial)"
    text = " and ".join(f"{c}{a}" for a in on)
    if tel:
        text += ", to the " + " and ".join(f"{a} telomere" for a in tel)
    return text


def mocha_filter(r):
    """The reason MoChA's UK Biobank filter drops a call, or None when it keeps it."""
    if r["type"].startswith("CNP"):
        return "germline copy-number polymorphism (MoChA type " + r["type"] + ")"
    bdev, n50, conc = num(r["bdev"]), num(r["n50_hets"]), num(r["lod_baf_conc"])
    if not (bdev < 0.1 or n50 < 2e5 or (not math.isnan(conc) and conc > 10.0)):
        return "likely constitutional (BAF deviation 0.1 or more with sparse heterozygous sites)"
    if math.isnan(num(r["bdev_se"])) and not num(r["lod_baf_phase"]) > 10.0:
        return "weak: no BAF deviation estimate and a BAF-and-phase LOD of 10 or less"
    rel, length = num(r["rel_cov"]), num(r["length"])
    if not (rel < 2.1 or bdev < 0.05 or (length > 5e5 and bdev < 0.1 and rel < 2.5)
            or (length > 5e6 and bdev < 0.15)):
        return "likely germline duplication (relative coverage 2.1 or more)"
    return None


def fmt(x, digits=4):
    return "." if isinstance(x, float) and math.isnan(x) else f"{x:.{digits}f}"


def report(sample, calls, stats, sites, phasing):
    """(call rows, summary dict) from MoChA's rows, the sites table and how the VCF was phased."""
    st = next((r for r in stats if r["sample_id"] == sample), None)
    if st is None and len(stats) == 1:
        st = stats[0]
    if st is None:
        raise Unreadable(f"MoChA's stats table has no row for {sample}")
    auto = [s for s in sites if s[0] in AUTOSOMES]
    span = sum(max(0, last - first + 1) for _, first, last, *_ in auto)
    hets = sum(s[4] for s in auto)
    phased = sum(s[5] for s in auto)
    depth = num(st["cov_median"])
    h10 = hets * FLOOR_MB * 1e6 / span if span else math.nan
    fr = floor_read(depth, h10)
    floor = max(FLOOR_MIN, fr) if not math.isnan(fr) else math.nan
    call_rate, baf_auto = num(st["call_rate"]), num(st["baf_auto"])
    summary = {
        "status": "assessed", "message": ".",
        "depth": fmt(depth, 1), "het_sites": str(hets), "phased_het_sites": str(phased),
        "het_per_10mb": fmt(h10, 0), "analysed_mb": fmt(span / 1e6, 1),
        "chromosomes": ",".join(s[0] for s in auto) or ".",
        "call_rate": fmt(call_rate), "baf_auto": fmt(baf_auto),
        "floor_read": fmt(fr, 3), "floor_min": fmt(FLOOR_MIN, 2), "floor": fmt(floor, 3),
        "calls_mocha": str(len(calls)), "calls_reported": "0", "phasing": phasing,
    }
    why = None
    if not auto or not hets:
        why = "no heterozygous site on chr1 to chr22 reached MoChA"
    elif not phased:
        why = "no heterozygous site is phased, so MoChA cannot follow a haplotype"
    elif math.isnan(floor):
        why = "MoChA gave no median depth, so no floor can be set"
    elif call_rate < CALL_RATE_MIN:
        why = f"call rate {call_rate:.4f} is under {CALL_RATE_MIN} (MoChA's sample filter)"
    elif baf_auto > BAF_AUTO_MAX:
        why = f"phased BAF autocorrelation {baf_auto:.4f} is over {BAF_AUTO_MAX} (MoChA's sample filter): noisy data"
    if why:
        summary.update(status="not assessed", message=why)
    rows = []
    for r in calls:
        if r.get("sample_id") not in (None, sample) and len(stats) > 1:
            continue
        cf = num(r["cf"])
        reason = None
        if r["chrom"] not in AUTOSOMES:
            reason = "not an autosome (mosaic loss of Y and X: step 33)"
        if reason is None:
            reason = mocha_filter(r)
        if reason is None and math.isnan(cf):
            reason = "cell fraction not estimated"
        if reason is None and not math.isnan(floor) and cf < floor:
            reason = f"cell fraction {cf:.3f} under the floor {floor:.3f}"
        if reason is None and why:
            reason = "sample not assessed"
        beg = next((v for k, v in r.items() if k.startswith("beg_")), ".")
        end = next((v for k, v in r.items() if k.startswith("end_")), ".")
        rows.append({
            "chrom": r["chrom"], "start": beg, "end": end, "length_mb": fmt(num(r["length"]) / 1e6, 2),
            "arms": arms(r["chrom"], r["p_arm"], r["q_arm"]), "type": TYPES.get(r["type"], r["type"]),
            "cell_fraction": fmt(cf, 3), "bdev": fmt(num(r["bdev"])), "rel_cov": fmt(num(r["rel_cov"])),
            "n_hets": r["n_hets"], "lod_baf_phase": r["lod_baf_phase"], "lod_lrr_baf": r["lod_lrr_baf"],
            "reported": "no" if reason else "yes", "reason": reason or ".",
        })
    summary["calls_reported"] = str(sum(1 for x in rows if x["reported"] == "yes"))
    return rows, summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--phasing", required=True, choices=("beagle", "input"))
    ap.add_argument("--out-calls", required=True)
    ap.add_argument("--out-summary", required=True)
    a = ap.parse_args(argv)
    try:
        rows, summary = report(a.sample, read_table(a.calls, NEED_CALLS), read_table(a.stats, NEED_STATS),
                               read_sites(a.sites), a.phasing)
    except (OSError, Unreadable, csv.Error) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    with open(a.out_calls, "w", encoding="utf-8") as f:
        f.write("\t".join(CALL_COLUMNS) + "\n")
        for x in rows:
            f.write("\t".join(str(x[k]) for k in CALL_COLUMNS) + "\n")
    with open(a.out_summary, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        for k, v in summary.items():
            f.write(f"{k}\t{v}\n")
    print(f"{a.sample}: {summary['status']}; {summary['calls_reported']} of {summary['calls_mocha']} MoChA calls "
          f"reported (floor {summary['floor']}, depth {summary['depth']}x)"
          + ("" if summary["message"] == "." else f": {summary['message']}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
