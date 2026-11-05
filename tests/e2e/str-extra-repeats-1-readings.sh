#!/usr/bin/env bash
# The counts of the repeats STRchive leaves out, read into results, with the
# real tools: ExpansionHunter 5.0.0 (step 09) sizes the CFTR TG/poly-T tract
# and the UGT1A1 TA repeat of assets/str_extra_loci.json on made-up reads,
# and the reports (step 24 and generate-report.sh, bin/collect_summary.py and
# bin/render_report.py) read them:
#   SXA  CFTR TG12-5T on one copy and TG11-7T on the other, UGT1A1 7 and 8 TA
#        units, no R117H copy (a made-up step 42 traits row), PharmCAT *1/*28
#        (a made-up step 27 table)
#        -> the two CFTR records named CFTR_TG and CFTR_T, TG11/TG12 5T/7T with
#           the 5T line, UGT1A1 *1/*28 that agrees with PharmCAT, and the
#           5-HTTLPR not-tested line
#   SXB  TG11-7T on both copies, UGT1A1 7 units on both (the red control)
#        -> TG11-7T/TG11-7T, no 5T line, UGT1A1 *1/*1
# The reference is made up: chr7 and chr2 are N except a window around each
# tract, the GRCh38 tract at its GRCh38 place (CFTR TG11 T7, UGT1A1 A(TA)6TAA).
# The catalog is the asset itself, copied under GENOME_DIR (setup.sh's
# merge with STRchive is the genotype-rule-rows-1-str-extra case).
# What this cannot show: reads of a real person at these tracts, which TG
# count sits on the 5T copy when both counts differ (ExpansionHunter writes
# the two tracts as two records), or R117H phase (17.6 kb from the tract).
# Writes two sample folders and a reference folder under GENOME_DIR, all
# removed at the end; no other case reads them.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
RD="${G}/reference-str-readings"
rm -rf "$RD"
mkdir -p "$RD"
cp "${REPO}/assets/str_extra_loci.json" "${G}/reference-str-readings.json"
for s in SXA SXB; do rm -rf "${G:?}/${s}"; mkdir -p "${G}/${s}/aligned"; done
python3 - "$RD" "$G" <<'PY'
import random, sys
rd, g = sys.argv[1], sys.argv[2]
rng = random.Random(11)
W = 1500  # random bases each side of a tract
# chrom: (contig length, 0-based tract start, left bases, reference tract units, right bases)
LOCI = {"chr7": (117551000, 117548606, "TTTTGA", [("TG", 11), ("T", 7)], "AACAGG"),
        "chr2": (233763000, 233760234, "TTGCCA", [("TA", 7)], "AGTAGG")}
win = {}
for c, (length, start, left, units, right) in LOCI.items():
    lf = "".join(rng.choice("ACGT") for _ in range(W - len(left))) + left
    rf = right + "".join(rng.choice("ACGT") for _ in range(W - len(right)))
    win[c] = (start - W, lf, rf)
    tract = "".join(m * n for m, n in units)
    with open(f"{rd}/ref.fa", "a") as f:
        f.write(f">{c}\n")
        seq_start = start - W
        body = "N" * seq_start + lf + tract + rf
        body += "N" * (length - len(body))
        for i in range(0, len(body), 60):
            f.write(body[i:i + 60] + "\n")

def reads(c, haps, tag):
    """Read pairs from each haplotype's tract; the CIGAR puts the tract's length change at its end."""
    w0, lf, rf = win[c]
    ref_tract = "".join(m * n for m, n in LOCI[c][3])
    out = []
    for h, units in enumerate(haps):
        tract = "".join(m * n for m, n in units)
        seq = lf + tract + rf
        t0, t1 = len(lf), len(lf) + len(tract)
        d = len(tract) - len(ref_tract)

        def place(a, n=150):
            """(1-based POS, CIGAR) of hap bases [a, a+n), or None when the read starts in the tract."""
            if a + n <= t0:
                return w0 + a + 1, f"{n}M"
            if a >= t1:
                return w0 + a - d + 1, f"{n}M"
            if a >= t0:
                return None
            if a + n <= t1:
                return None
            pre = t0 - a + min(len(tract), len(ref_tract))
            post = n - (t0 - a) - len(tract)
            mid = f"{d}I" if d > 0 else (f"{-d}D" if d < 0 else "")
            cig = f"{pre}M{mid}{post}M" if mid else f"{n}M"
            return w0 + a + 1, cig
        k = 0
        for a in range(0, len(seq) - 400, 4):
            f_len = 300 + (a % 7) * 10
            b = a + f_len - 150
            p1, p2 = place(a), place(b)
            if not p1 or not p2:
                continue
            k += 1
            name = f"{tag}{c}h{h}r{k}"
            r1, r2 = seq[a:a + 150], seq[b:b + 150]
            tlen = p2[0] + 150 - p1[0]
            out.append(f"{name}\t99\t{c}\t{p1[0]}\t60\t{p1[1]}\t=\t{p2[0]}\t{tlen}\t{r1}\t{'I' * 150}")
            out.append(f"{name}\t147\t{c}\t{p2[0]}\t60\t{p2[1]}\t=\t{p1[0]}\t{-tlen}\t{r2}\t{'I' * 150}")
    return out

SAMPLES = {"SXA": ([[("TG", 12), ("T", 5)], [("TG", 11), ("T", 7)]], [[("TA", 7)], [("TA", 8)]]),
           "SXB": ([[("TG", 11), ("T", 7)]] * 2, [[("TA", 7)]] * 2)}
for s, (cftr, ugt) in SAMPLES.items():
    with open(f"{g}/{s}/aligned/{s}.sam", "w") as f:
        f.write("@HD\tVN:1.6\tSO:unsorted\n")
        for c, (length, *_r) in LOCI.items():
            f.write(f"@SQ\tSN:{c}\tLN:{length}\n")
        f.write(f"@RG\tID:{s}\tSM:{s}\n")
        for line in reads("chr7", cftr, s) + reads("chr2", ugt, s):
            f.write(line + f"\tRG:Z:{s}\n")
PY
check "the made-up reference was written" test -s "${RD}/ref.fa"
sam faidx reference-str-readings/ref.fa
for s in SXA SXB; do
  in_genome "$SAMTOOLS_IMAGE" bash -c "samtools sort -o '${s}/aligned/${s}_sorted.bam' '${s}/aligned/${s}.sam' && samtools index '${s}/aligned/${s}_sorted.bam'"
  check "${s}: the BAM is indexed" test -s "${G}/${s}/aligned/${s}_sorted.bam.bai"
done

# --- step 09, then the made-up step 42 and step 27 tables of SXA -------------------------------
export REF_FASTA="${RD}/ref.fa"
for s in SXA SXB; do
  EH_CATALOG="${G}/reference-str-readings.json" run_step 09-expansion-hunter.sh "$s" female
  check_step_exit 09-expansion-hunter.sh
done
unset REF_FASTA
python3 - "$G" "$REPO" <<'PY'
import os, sys
g, repo = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(repo, "bin"))
import traits_report
os.makedirs(f"{g}/SXA/traits", exist_ok=True)
os.makedirs(f"{g}/SXA/cpic", exist_ok=True)
with open(f"{g}/SXA/traits/SXA_traits.tsv", "w") as f:
    f.write("\t".join(traits_report.OUT_COLUMNS) + "\n")
    f.write("rs78655421\tCFTR\tchr7\t117530975\tCFTR R117H\tG/G\tcalled\tNo R117H copy.\tno\t-\t32404922\n")
with open(f"{g}/SXA/cpic/SXA_phenotypes.tsv", "w") as f:
    f.write("Gene\tDiplotype\tPhenotype\tStatus\nUGT1A1\t*1/*28\tIntermediate Metabolizer\tnon-normal\n")
PY

# --- the reports ------------------------------------------------------------------------------
for s in SXA SXB; do
  run_step 24-html-report.sh "$s"
  check_step_exit 24-html-report.sh
  run_step generate-report.sh "$s"
  check_step_exit generate-report.sh
done
A=$(cat "${G}/SXA/SXA_report.txt" 2>/dev/null)
B=$(cat "${G}/SXB/SXB_report.txt" 2>/dev/null)
AH=$(sed 's/&#x27;/'"'"'/g; s/&amp;/\&/g' "${G}/SXA/SXA_report.html" 2>/dev/null)
grep -A3 'CFTR intron 9\|UGT1A1 promoter\|5-HTTLPR' <<< "$A" || true
check "SXA: the CFTR records are named CFTR_TG and CFTR_T" has '^    CFTR_TG +11/12' "$A"
check "SXA: and CFTR_T" has '^    CFTR_T +(5/7|7/5)' "$A"
check "SXA: CFTR TG11/TG12 5T/7T" has 'CFTR intron 9 \(legacy IVS8\) TG and poly-T tract: TG11/TG12 5T/7T$' "$A"
check "SXA: the 5T line, TG12-5T/TG11-7T among the pairings, no R117H copy" \
  has '^    5T: 1 copy with the 5T poly-T tract, TG11 or TG12 on the 5T copy.*TG12-5T/TG11-7T.*No R117H copy' "$A"
check "SXA: UGT1A1 *1/*28 (TA6/TA7), the same as PharmCAT" \
  has 'UGT1A1 promoter TA repeat: \*1/\*28$' "$A"
check "SXA: and the agreement with PharmCAT" has 'PharmCAT \(step 7\) called \*1/\*28: the same number' "$A"
check "SXA: the 5-HTTLPR line" has 'SLC6A4 5-HTTLPR promoter length polymorphism: not tested' "$A"
check "SXA: the HTML report has the 5T line" has '5T: 1 copy with the 5T poly-T tract' "$AH"
check "SXB (red control): TG11-7T/TG11-7T" has 'TG and poly-T tract: TG11-7T/TG11-7T$' "$B"
check "SXB (red control): no 5T line" lacks '5T: [0-9] cop' "$B"
check "SXB: UGT1A1 *1/*1, not compared without a PharmCAT call" has 'UGT1A1 promoter TA repeat: \*1/\*1$' "$B"

in_genome "$SAMTOOLS_IMAGE" rm -rf SXA SXB reference-str-readings reference-str-readings.json
finish
