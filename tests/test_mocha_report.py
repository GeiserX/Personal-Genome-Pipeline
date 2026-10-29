#!/usr/bin/env python3
"""bin/mocha_report.py: step 51's mosaic chromosomal alterations from MoChA's tables.

Each case writes MoChA's --stats and --calls tables and the per-chromosome
site counts with chosen numbers, runs the script, and checks what it writes:

  1. the floor: 6 / sqrt(D x H) at depth D and H heterozygous sites per 10 Mb
     of the chromosomes read, computed here by hand; at 30x with 3,000 sites
     per 10 Mb it is under 0.05, so the floor is 0.05; at 5x with 300 sites it
     is the read-sampling floor (0.155)
  2. a CN-LOH of 20 Mb at a cell fraction of 0.15 is reported, typed CN-LOH,
     with its arm; the same call at 0.04 is under the floor and not reported,
     with the reason; at 5x with 300 sites per 10 Mb a 0.15 call is under the
     floor of that depth
  3. MoChA's UK Biobank filter: a germline CNP, a constitutional event (BAF
     deviation 0.2 with sparse hets), a likely germline duplication (relative
     coverage 3, BAF deviation 0.16) and a call on chrX are each kept out, with their reason; a
     call with no cell fraction is not reported
  4. the sample filter: a call rate under 0.97 or a BAF autocorrelation over
     0.03 is 'not assessed', and then no call is reported; no phased het, or a
     call rate or BAF autocorrelation MoChA left as nan, is 'not assessed'
  5. the CLI: a stats table without cov_median, or a short sites row, exits 2
  6. collect_summary.py's mosaic section and render_report.py's lines read the
     tables: the reported call only; with no reported call the not-assessed
     list names autosomal mosaic alterations unless step 51 assessed the sample;
     a missing call table, one with fewer reported calls than the summary, or
     one without MoChA's columns is unreadable
  7. versions.env pins the MoChA image by digest and the plugin by commit,
     those of the e2e case (tests/e2e/mocha-1.sh)

Run: python3 tests/test_mocha_report.py
"""
import math
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "mocha_report.py")
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary as cs  # noqa: E402
import mocha_report as mr  # noqa: E402
import render_report as rr  # noqa: E402

FAILS = []
STATS_HEAD = ("sample_id\tcomputed_gender\tcall_rate\tcov_median\tcov_sd\tcov_auto\tbaf_corr\tbaf_conc\tbaf_auto"
              "\tn_sites\tn_hets\n")
CALLS_HEAD = ("sample_id\tcomputed_gender\tchrom\tbeg_GRCh38\tend_GRCh38\tlength\tp_arm\tq_arm\tn_sites\tn_hets"
              "\tn50_hets\tbdev\tbdev_se\trel_cov\trel_cov_se\tlod_lrr_baf\tlod_baf_phase\tn_flips\tbaf_conc"
              "\tlod_baf_conc\ttype\tcf\n")


def check(name, got, want):
    if got == want:
        print(f"PASS {name}")
    else:
        print(f"FAIL {name}: got {got!r}, want {want!r}")
        FAILS.append(name)


def call(chrom="chr20", beg=60000, end=20000000, p="T", q="N", n50=3000, bdev=0.075, bdev_se=0.001,
         rel=2.0, lod=900.0, typ="CN-LOH", cf=0.15):
    return (f"S\tU\t{chrom}\t{beg}\t{end}\t{end - beg}\t{p}\t{q}\t20000\t6600\t{n50}\t{bdev}\t{bdev_se}\t{rel}"
            f"\t0.01\tnan\t{lod}\t0\tnan\tnan\t{typ}\t{cf}\n")


def run(calls, depth=30.0, hets_per_10mb=3000, call_rate=0.999, baf_auto=0.001, phased=True, stats=None,
        sites=None):
    """Write the three inputs in a fresh folder, run the script; (exit, rows, summary, folder)."""
    d = tempfile.mkdtemp()
    span = 60_000_000
    h = int(hets_per_10mb * span / 1e7)
    with open(f"{d}/stats.tsv", "w") as f:
        f.write(stats or STATS_HEAD + f"S\tU\t{call_rate}\t{depth}\t3\t0.01\t0.001\t0.5\t{baf_auto}\t100000\t{h}\n")
    with open(f"{d}/calls.tsv", "w") as f:
        f.write(CALLS_HEAD + "".join(calls))
    with open(f"{d}/sites.tsv", "w") as f:
        f.write(sites or f"chr20\t1\t{span}\t{3 * h}\t{h}\t{h if phased else 0}\n")
    os.makedirs(f"{d}/S/mosaic")
    rc = subprocess.run([sys.executable, SCRIPT, "--sample", "S", "--calls", f"{d}/calls.tsv", "--stats",
                         f"{d}/stats.tsv", "--sites", f"{d}/sites.tsv", "--phasing", "input",
                         "--out-calls", f"{d}/S/mosaic/S_mocha.tsv", "--out-summary", f"{d}/S/mosaic/S_mocha_summary.tsv"],
                        capture_output=True, text=True).returncode
    if rc:
        return rc, [], {}, d
    rows = [dict(zip(mr.CALL_COLUMNS, ln.rstrip("\n").split("\t"))) for ln in open(f"{d}/S/mosaic/S_mocha.tsv")][1:]
    summ = dict(ln.rstrip("\n").split("\t") for ln in open(f"{d}/S/mosaic/S_mocha_summary.tsv"))
    return rc, rows, summ, d


# 1. the floor
check("1 30x, 3,000 hets per 10 Mb: the read-sampling floor by hand", round(mr.floor_read(30, 3000), 4),
      round(6 / math.sqrt(30 * 3000), 4))
_, _, s, _ = run([call()])
check("1 at 30x the floor is the minimum, 0.05", (s["floor"], s["floor_min"], s["floor_read"]), ("0.050", "0.05", "0.020"))
_, _, s, _ = run([call()], depth=5, hets_per_10mb=300)
check("1 at 5x with 300 hets per 10 Mb the floor is the read-sampling one", s["floor"], f"{6 / math.sqrt(5 * 300):.3f}")
check("1 the het density comes from the sites table", s["het_per_10mb"], "300")

# 2. reported, under the floor
rc, rows, s, _ = run([call()])
check("2 a 0.15 CN-LOH is reported", [(r["type"], r["cell_fraction"], r["reported"]) for r in rows],
      [("CN-LOH", "0.150", "yes")])
check("2 with its arm", rows[0]["arms"], "20p, to the p telomere")
check("2 the summary counts it", (s["status"], s["calls_reported"], s["calls_mocha"]), ("assessed", "1", "1"))
_, rows, s, _ = run([call(bdev=0.02, cf=0.04)])
check("2 a 0.04 CN-LOH is under the floor", (rows[0]["reported"], rows[0]["reason"]),
      ("no", "cell fraction 0.040 under the floor 0.050"))
_, rows, _, _ = run([call()], depth=5, hets_per_10mb=300)
check("2 at 5x with 300 hets the same 0.15 call is under that depth's floor", rows[0]["reported"], "no")

# 3. MoChA's filter
_, rows, s, _ = run([call(typ="CNP_Loss"), call(bdev=0.2, n50=500000, cf=0.4), call(rel=3.0, bdev=0.16, typ="Gain",
                     cf=0.6), call(chrom="chrX"), call(cf="nan"), call(beg=30000000, end=60000000, p="N", q="T",
                     typ="Loss", bdev=0.03, cf=0.12)])
check("3 a CNP is left out", rows[0]["reason"].startswith("germline copy-number polymorphism"), True)
check("3 a constitutional event is left out", rows[1]["reason"].startswith("likely constitutional"), True)
check("3 a likely germline duplication is left out", rows[2]["reason"].startswith("likely germline duplication"), True)
check("3 chrX is step 33's", rows[3]["reason"].startswith("not an autosome"), True)
check("3 no cell fraction: not reported", rows[4]["reason"], "cell fraction not estimated")
check("3 the loss on 20q is the one reported", [(r["type"], r["arms"]) for r in rows if r["reported"] == "yes"],
      [("loss", "20q, to the q telomere")])
check("3 calls_reported 1 of 6", (s["calls_reported"], s["calls_mocha"]), ("1", "6"))

# 4. the sample filter
for kw, why in (({"call_rate": 0.95}, "call rate"), ({"baf_auto": 0.05}, "phased BAF autocorrelation"),
                ({"phased": False}, "no heterozygous site is phased"),
                ({"call_rate": "nan"}, "MoChA gave no call rate or phased BAF autocorrelation"),
                ({"baf_auto": "nan"}, "MoChA gave no call rate or phased BAF autocorrelation")):
    _, rows, s, _ = run([call()], **kw)
    check(f"4 {why}: not assessed", (s["status"], s["message"].startswith(why.split(':')[0]) or why in s["message"]),
          ("not assessed", True))
    check(f"4 {why}: no call reported", [r["reported"] for r in rows], ["no"])

# 5. the CLI
rc, *_ = run([call()], stats="sample_id\tcall_rate\nS\t0.99\n")
check("5 a stats table without cov_median exits 2", rc, 2)
rc, *_ = run([call()], sites="chr20\t1\t100\n")
check("5 a short sites row exits 2", rc, 2)

# 6. both reports
_, _, _, d = run([call(), call(chrom="chrX")])
path, data = cs.sec_mocha(f"{d}/S", "S")
check("6 the section lists the reported call only", [(x["type"], x["chrom"], x["cell_fraction"]) for x in data["calls"]],
      [("CN-LOH", "chr20", "0.150")])
check("6 and counts the other", data["not_reported"], 1)
lines = rr.mosaic_mca_lines(data)
check("6 the text names the call", lines[1], "CN-LOH on 20p, to the p telomere, chr20:60000-20000000 (19.94 Mb): "
      "about 0.150 of the cells")
check("6 the first line gives the floor", lines[0].startswith("Reported at a cell fraction of 0.050 or more at 30.0x"), True)
summ = cs.collect("S", f"{d}/S")
check("6 assessed: no not-assessed line", cs.MCA_NOT_ASSESSED in summ["not_assessed"], False)
os.remove(f"{d}/S/mosaic/S_mocha_summary.tsv")
summ = cs.collect("S", f"{d}/S")
check("6 without step 51: the not-assessed line", cs.MCA_NOT_ASSESSED in summ["not_assessed"], True)
_, _, _, d = run([call()], baf_auto=0.05)
data = cs.sec_mocha(f"{d}/S", "S")[1]
check("6 not assessed: the text says why", rr.mosaic_mca_lines(data)[0].startswith("Not assessed: phased BAF"), True)
check("6 and the not-assessed line stays", cs.MCA_NOT_ASSESSED in cs.collect("S", f"{d}/S")["not_assessed"], True)
check("6 no reported call: the text says none is reported, not that none exists",
      rr.mosaic_mca_lines(cs.sec_mocha(f"{run([call(cf=0.04)])[3]}/S", "S")[1])[1].startswith("No mosaic alteration reported"), True)
check("6 the pipeline-wide not-assessed list keeps large mosaic alterations out of its chrM-only line",
      [x for x in cs.NOT_ASSESSED if x.startswith("Mosaic")], ["Mosaic and low-fraction small variants (only chrM heteroplasmy is reported, and with step 49 the CHIP genes)"])


def unreadable(d):
    try:
        cs.sec_mocha(f"{d}/S", "S")
    except cs.Unreadable as e:
        return str(e)
    return None


_, _, _, d = run([call()])
os.remove(f"{d}/S/mosaic/S_mocha.tsv")
check("6 a summary without its call table is unreadable", bool(unreadable(d)), True)
_, _, _, d = run([call()])
with open(f"{d}/S/mosaic/S_mocha.tsv") as f:
    head = f.readline()
with open(f"{d}/S/mosaic/S_mocha.tsv", "w") as f:
    f.write(head)
check("6 a call table with fewer reported calls than the summary is unreadable", bool(unreadable(d)), True)
with open(f"{d}/S/mosaic/S_mocha.tsv", "w") as f:
    f.write("chrom\tstart\n")
check("6 a header-only call table without MoChA's columns is unreadable", bool(unreadable(d)), True)

# 7. the pins
env = open(os.path.join(REPO, "versions.env")).read()
img = re.search(r'^MOCHA_IMAGE="([^"]+)"', env, re.M)
name, _, digest = (img.group(1) if img else "").partition("@")
check("7 MOCHA_IMAGE is the authors' image, tag 1.22-20250819", (name.split("/")[0], name.rsplit(":", 1)[-1]),
      ("us.gcr.io", "1.22-20250819"))
check("7 and pinned by its digest", digest, "sha256:c9398385cfe3a968918ddf470e2769b7ee42cbc0bb8f564b4744c460343550a5")
com = re.search(r'^MOCHA_COMMIT="([^"]+)"', env, re.M)
check("7 MOCHA_COMMIT is freeseek/mocha 95686b7b", com and com.group(1), "95686b7b65f53a490513be76bb120e5fc20a8bcf")
ver = re.search(r'^MOCHA_VERSION="([^"]+)"', env, re.M)
check("7 MOCHA_VERSION is the define of mocha.c at that commit", ver and ver.group(1), "2025-08-19")

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
