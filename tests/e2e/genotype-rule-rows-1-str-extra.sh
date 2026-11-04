#!/usr/bin/env bash
# The repeats STRchive leaves out (assets/str_extra_loci.json) with the real tools:
#   1. STR_EXTRA_LOCI=yes setup.sh --strchive-data appends them to STRchive's
#      catalog (82 + 2 entries, the asset's sha256 in catalog.stamp); the same
#      folder without the opt-in converts again and is back to 82 entries.
#   2. ExpansionHunter 5.0.0 (step 09) on made-up reads, against a made-up
#      reference whose chr7 and chr2 are N except a window around each tract,
#      the GRCh38 tract at its GRCh38 place (CFTR TG11 T7, UGT1A1 A(TA)6TAA):
#        EXA  CFTR TG12-T5 on one copy and TG11-T7 on the other, UGT1A1 7 and
#             8 TA units (TA6/TA7, *1/*28)
#             -> CFTR_TG 11/12, CFTR_T 5/7, UGT1A1 7/8
#        EXB  TG11-T7 on both copies, UGT1A1 7 units on both
#             -> CFTR_TG 11/11, CFTR_T 7/7 (no 5T: the red control), UGT1A1 7/7
#      and step 09 prints both loci as not graded, with the asset's reason.
# What this cannot show: reads of a real person at these tracts, or a reading
# of the counts (the 5T line and *28 are not written by any step yet).
# Writes two sample folders and a reference folder under GENOME_DIR, all removed
# at the end; no other case reads them.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
EXTRA="${REPO}/assets/str_extra_loci.json"
SETUP="${CASE_TMP}/setup"
rm -rf "$SETUP"
mkdir -p "$SETUP"
CAT="${SETUP}/reference/strchive-${STRCHIVE_VERSION}/catalog.json"
count() { python3 -c 'import json, sys; print(len(json.load(open(sys.argv[1]))))' "$1" 2>/dev/null || echo 0; }

# --- 1. setup, with and without the opt-in ------------------------------------------------
STR_EXTRA_LOCI=yes "${REPO}/scripts/setup.sh" --strchive-data "$SETUP" > "${CASE_TMP}/setup.log" 2>&1
check_eq "STR_EXTRA_LOCI=yes setup.sh --strchive-data exits 0" "$?" 0
cat "${CASE_TMP}/setup.log"
check_eq "the catalog has STRchive's 82 entries and the 2 extra loci" "$(count "$CAT")" 84
check "catalog.stamp names the extra catalog's sha256" has "extra=$(sha256sum "$EXTRA" | cut -d' ' -f1)" \
  "$(cat "${CAT%/*}/catalog.stamp" 2>/dev/null)"
python3 - "$CAT" "${G}/reference-str-extra.json" <<'PY'
import json, sys
json.dump([e for e in json.load(open(sys.argv[1])) if e.get("Source") != "STRchive"], open(sys.argv[2], "w"), indent=1)
PY
check_eq "the extra loci of the merged catalog" \
  "$(python3 -c 'import json, sys; print(",".join(e["LocusId"] for e in json.load(open(sys.argv[1]))))' "${G}/reference-str-extra.json")" \
  "CFTR,UGT1A1"
"${REPO}/scripts/setup.sh" --strchive-data "$SETUP" > "${CASE_TMP}/setup-off.log" 2>&1
check_eq "control: setup.sh again without the opt-in exits 0" "$?" 0
check_eq "control: without the opt-in the catalog is back to 82 entries" "$(count "$CAT")" 82
check "control: catalog.stamp says extra=off" has 'extra=off$' "$(cat "${CAT%/*}/catalog.stamp" 2>/dev/null)"

# --- 2. made-up reference and reads ---------------------------------------------------------
RD="${G}/reference-str-extra"
rm -rf "$RD"
mkdir -p "$RD"
for s in EXA EXB; do rm -rf "${G:?}/${s}"; mkdir -p "${G}/${s}/aligned"; done
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

SAMPLES = {"EXA": ([[("TG", 12), ("T", 5)], [("TG", 11), ("T", 7)]], [[("TA", 7)], [("TA", 8)]]),
           "EXB": ([[("TG", 11), ("T", 7)]] * 2, [[("TA", 7)]] * 2)}
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
sam faidx reference-str-extra/ref.fa
check_eq "ref.fa.fai: chr7 and chr2" "$(cut -f1 "${RD}/ref.fa.fai" 2>/dev/null | tr '\n' ' ')" "chr7 chr2 "
check_eq "the CFTR tract on the made-up chr7 is GRCh38's TG11 T7" \
  "$(sam faidx reference-str-extra/ref.fa chr7:117548601-117548641 2>/dev/null | tail -n +2 | tr -d '\n')" \
  "TTTTGATGTGTGTGTGTGTGTGTGTGTGTTTTTTTAACAGG"
for s in EXA EXB; do
  in_genome "$SAMTOOLS_IMAGE" bash -c "samtools sort -o '${s}/aligned/${s}_sorted.bam' '${s}/aligned/${s}.sam' && samtools index '${s}/aligned/${s}_sorted.bam'"
  check "${s}: the BAM is indexed" test -s "${G}/${s}/aligned/${s}_sorted.bam.bai"
done

# --- 3. step 09 on the extra loci -----------------------------------------------------------
repcn() {  # repcn <sample> <VARID>: REPCN of that record
  python3 - "${G}/$1/expansion_hunter/$1_eh.vcf" "$2" <<'PY' 2>/dev/null
import sys
for line in open(sys.argv[1]):
    if line.startswith("#"):
        continue
    c = line.rstrip("\n").split("\t")
    info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
    if info.get("VARID") == sys.argv[2]:
        print(dict(zip(c[8].split(":"), c[9].split(":"))).get("REPCN", "."))
PY
}
export REF_FASTA="${RD}/ref.fa"
for s in EXA EXB; do
  EH_CATALOG="${G}/reference-str-extra.json" run_step 09-expansion-hunter.sh "$s" female
  check_step_exit 09-expansion-hunter.sh
  cp "$STEP_LOG" "${CASE_TMP}/${s}.log"
  check "${s}: step 09 prints CFTR as not graded, with the asset's reason" \
    has 'Not graded.*CFTR: the CFTR intron 9 TG and poly-T tract' "$(cat "$STEP_LOG")"
  check "${s}: step 09 prints UGT1A1 as not graded, with the asset's reason" \
    has 'Not graded.*UGT1A1: the UGT1A1 promoter TA repeat' "$(cat "$STEP_LOG")"
  check_eq "${s}: one ExpansionHunter record per repeat (CFTR_TG, CFTR_T, UGT1A1)" \
    "$(grep -vc '^#' "${G}/${s}/expansion_hunter/${s}_eh.vcf" 2>/dev/null || echo 0)" 3
done
unset REF_FASTA
echo "counts: EXA $(repcn EXA CFTR_TG) $(repcn EXA CFTR_T) $(repcn EXA UGT1A1); EXB $(repcn EXB CFTR_TG) $(repcn EXB CFTR_T) $(repcn EXB UGT1A1)"
check_eq "EXA: CFTR TG counts" "$(repcn EXA CFTR_TG)" "11/12"
check_eq "EXA: CFTR poly-T counts, one 5T copy" "$(repcn EXA CFTR_T)" "5/7"
check_eq "EXA: UGT1A1 TA units 7/8 (TA6/TA7, *1/*28)" "$(repcn EXA UGT1A1)" "7/8"
check_eq "EXB: CFTR TG counts" "$(repcn EXB CFTR_TG)" "11/11"
check_eq "EXB (red control): CFTR poly-T 7/7, no 5T copy" "$(repcn EXB CFTR_T)" "7/7"
check_eq "EXB: UGT1A1 TA units 7/7 (*1/*1)" "$(repcn EXB UGT1A1)" "7/7"

in_genome "$SAMTOOLS_IMAGE" rm -rf EXA EXB reference-str-extra reference-str-extra.json
finish
