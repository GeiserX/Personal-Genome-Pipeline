#!/usr/bin/env python3
"""Step 52's table: VNtyper's MUC1 VNTR result in plain words.

ADTKD-MUC1 is caused by frameshifts inside the MUC1 coding VNTR, a repeat of
60 bp units, most often one extra C in a run of seven (dupC; Kirby et al.,
Nat Genet 2013, PMID 23396133). Short reads cannot be placed inside the
repeat, so standard callers do not see it. VNtyper (Saei et al., iScience
2023, PMID 37456840) genotypes it with k-mers (Kestrel); this script reads
the run folder VNtyper writes and says positive, possible, negative or not
assessed, with the share of known positives VNtyper found at the depth this
sample has over the VNTR.

    vntyper_report.py --run-dir DIR --out TABLE [--exit-code N --log FILE]

DIR is VNtyper's --output-dir: kestrel/kestrel_result.tsv and
coverage/coverage_summary.tsv. --exit-code is VNtyper's own exit status; a
non-zero one, or a missing file, gives "not assessed"; --log, VNtyper's
console output, tells a BAM with no read over MUC1 (VNtyper 2.0.47 stops on
it with "FASTQ layout 'empty'") from other errors. TABLE gets metric<TAB>
value rows. Standard library only: it runs in VNTYPER_IMAGE's Python, in the
bash step and in the VNTYPER process alike.
"""
import argparse
import csv
import os
import sys

CALLS = ("High_Precision*", "High_Precision", "Low_Precision")

# What VNtyper 2.0.47 says about sensitivity at depth (docs/user-guide/
# configuration.md, "Coverage thresholds"): confirmed positive GRCh38 exomes
# downsampled and re-run, 42 samples at nine fractions each. Rows: lowest
# mean depth over the VNTR, detected runs, runs.
SENSITIVITY = [(100.0, 160, 161), (50.0, 95, 108), (0.0, 48, 109)]
SENSITIVITY_SOURCE = ("VNtyper 2.0.47's downsampling of 42 confirmed positive GRCh38 exomes "
                      "(repeated runs, not independent tests)")

FIELDS = ("status", "confidence", "variant", "flag", "alt_depth", "region_depth", "depth_score",
          "vntr_mean_depth", "vntr_uncovered_pct", "coverage_qc", "sensitivity", "message", "advice")

ADVICE_CALL = ("A research call, not a diagnosis. VNtyper's authors ask for confirmation by another method: "
               "SNaPshot for dupC, long-read sequencing for other variants. ADTKD-MUC1 is dominant: "
               "discuss it with a nephrologist or a clinical geneticist.")
ADVICE_NONE = ("A negative here lowers the chance of the frameshifts VNtyper looks for; it does not rule "
               "ADTKD-MUC1 out. A family history of kidney failure with gout or a bland urine sediment "
               "still deserves a clinical test (SNaPshot for dupC).")


def rows(path):
    """The data rows of a VNtyper TSV (its ## banner lines skipped)."""
    with open(path, newline="") as f:
        lines = [x for x in f if not x.startswith("#") and x.strip()]
    return list(csv.DictReader(lines, delimiter="\t"))


def number(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else x  # NaN is no number


def sensitivity_line(mean):
    if mean is None:
        return "Sensitivity not known: VNtyper measured no depth over the VNTR."
    for low, hit, runs in SENSITIVITY:
        if mean >= low:
            band = (f"at {low:.0f}x or more" if low == 100 else
                    f"between {low:.0f}x and 100x" if low else "below 50x")
            return (f"Mean depth over the MUC1 VNTR: {mean:.0f}x. At that depth ({band}) VNtyper found the "
                    f"known variant in {hit} of {runs} runs ({100 * hit / runs:.0f}%), in {SENSITIVITY_SOURCE}.")
    raise AssertionError("unreachable")


NO_READS = "FASTQ layout 'empty'"
REGION = "chr1:155,184,000-155,194,000"


def read(run_dir, exit_code, log_text=""):
    out = dict.fromkeys(FIELDS, ".")
    kpath = os.path.join(run_dir, "kestrel", "kestrel_result.tsv")
    cpath = os.path.join(run_dir, "coverage", "coverage_summary.tsv")
    if exit_code:
        why = (f"the BAM has no read in the MUC1 region ({REGION}), so VNtyper had nothing to read"
               if NO_READS in log_text else f"VNtyper stopped with exit code {exit_code}: see its log")
        out.update(status="not assessed", message=why[0].upper() + why[1:] + ".", advice=ADVICE_NONE)
        return out
    missing = [os.path.relpath(p, run_dir) for p in (kpath, cpath) if not os.path.isfile(p)]
    if missing:
        out.update(status="not assessed", message=f"VNtyper wrote no {' and no '.join(missing)}.",
                   advice=ADVICE_NONE)
        return out
    cov = (rows(cpath) or [{}])[0]
    mean, uncovered = number(cov.get("mean")), number(cov.get("percent_uncovered"))
    qc = (cov.get("coverage_qc") or "").strip() or "NOT_EVALUATED"
    k = rows(kpath)
    best = k[0] if k else {}
    conf = (best.get("Confidence") or "").strip() or "Negative"
    flag = (best.get("Flag") or "").strip()
    flagged = bool(flag) and flag != "Not flagged"
    covered = qc in ("PASS", "REDUCED")
    out.update(confidence=conf, coverage_qc=qc, sensitivity=sensitivity_line(mean),
               vntr_mean_depth=f"{mean:.2f}" if mean is not None else ".",
               vntr_uncovered_pct=f"{uncovered:.2f}" if uncovered is not None else ".")
    if conf in CALLS:
        name = (best.get("Nomenclature") or "").strip()
        out.update(variant=name or (best.get("Variant") or ".").strip() or ".", flag=flag or ".",
                   alt_depth=best.get("Estimated_Depth_AlternateVariant") or ".",
                   region_depth=best.get("Estimated_Depth_Variant_ActiveRegion") or ".",
                   depth_score=best.get("Depth_Score") or ".", advice=ADVICE_CALL)
        what = f"a MUC1 VNTR frameshift ({out['variant']}), Kestrel confidence {conf}"
        if conf.startswith("High_Precision") and not flagged and (uncovered is None or uncovered <= 50):
            out.update(status="positive", message=f"VNtyper calls {what}.")
        else:
            why = ([f"flagged: {flag}"] if flagged else []) \
                + (["low precision"] if conf == "Low_Precision" else []) \
                + ([f"{uncovered:.0f}% of the VNTR has no read"] if uncovered is not None and uncovered > 50 else [])
            out.update(status="possible", message=f"VNtyper calls {what}, with caveats ({'; '.join(why)}).")
        return out
    if covered:
        out.update(status="negative", message="VNtyper calls no MUC1 VNTR frameshift.", advice=ADVICE_NONE)
    else:
        depth = f"{mean:.0f}x" if mean is not None else "no measured depth"
        out.update(status="not assessed", advice=ADVICE_NONE,
                   message=(f"VNtyper calls no frameshift, but its coverage check is {qc} ({depth} over the "
                            "VNTR, 50x is its floor on GRCh38), so a frameshift cannot be ruled out."))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run-dir", required=True, help="VNtyper's --output-dir")
    ap.add_argument("--out", required=True, help="the table to write")
    ap.add_argument("--exit-code", type=int, default=0, help="VNtyper's exit status")
    ap.add_argument("--log", help="VNtyper's console output")
    a = ap.parse_args(argv)
    try:
        log_text = open(a.log, errors="replace").read() if a.log and os.path.isfile(a.log) else ""
        t = read(a.run_dir, a.exit_code, log_text)
    except (OSError, csv.Error, UnicodeDecodeError) as e:
        print(f"ERROR: cannot read VNtyper's results in {a.run_dir}: {e}", file=sys.stderr)
        return 2
    tmp = a.out + ".tmp"
    with open(tmp, "w") as f:
        f.write("metric\tvalue\n")
        for k in FIELDS:
            f.write(f"{k}\t{str(t[k]).replace(chr(9), ' ').replace(chr(10), ' ')}\n")
    os.replace(tmp, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
