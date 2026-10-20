#!/usr/bin/env bash
# PharmCAT's MT-RNR1 call from step 20's chrM calls, and no logback noise.
#
# PharmCAT 3.4.0 reads MT-RNR1 only as an outside call. Step 07 runs
# bin/pgx_mito_merge.py on mito/<sample>_chrM_filtered.vcf.gz when it exists
# and adds the MT-RNR1 line to the outside calls. Four samples, each HG002's
# variants-only VCF (case 21) with an invented Mutect2-style chrM call set,
# in a GENOME_DIR of their own (the fixture reference hard-linked into it),
# so no later case sees them:
#   MTHOM  m.1555A>G at AF 0.99, PASS: PharmCAT reports MT-RNR1 m.1555A>G as
#          an outside call, with an aminoglycoside recommendation for it;
#   MTREF  PASS chrM records, none at an MT-RNR1 allele: MT-RNR1 Reference;
#   MTHET  m.1555A>G at AF 0.40: nothing passed, the log says heteroplasmic,
#          MT-RNR1 stays No Result;
#   MTNONE no chrM call set (the control, the state before the fix for every
#          sample): the log says so and MT-RNR1 stays No Result.
# The jar ran in /pharmcat, which the calling user cannot write, and logback
# printed a "pharmcat.log (Permission denied)" stack trace on every run. Step
# 07 now runs it in /tmp: no step log here has the line. Control, in the same
# run: the jar on MTHOM's preprocessed VCF as before (calling user, the
# image's working directory) prints it.
. "$(dirname "$0")/lib.sh"

G="${E2E_WORK}/genome-mtrnr1"
REF="reference/GRCh38_no_alt_analysis_set.fasta"
SRC="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
check "case 21 wrote HG002's VCF" test -s "$SRC"
rm -rf "$G"
mkdir -p "${G}/reference"
for f in "$REF" "${REF}.fai"; do
  ln -f "${GENOME_DIR}/${f}" "${G}/${f}" 2>/dev/null || cp "${GENOME_DIR}/${f}" "${G}/${f}"
done

HEADER='##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=weak_evidence,Description="Mutation does not meet likelihood threshold">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##contig=<ID=chrM,length=16569>'
row() { printf 'chrM\t%s\t.\t%s\t%s\t.\t%s\t.\tGT:AF:DP\t0/1:%s:3000\n' "$@"; }
BACKGROUND="$(row 263 A G PASS 0.998)
$(row 750 A G PASS 0.999)
$(row 3010 G A weak_evidence 0.012)"

# sample P [CHRM_RECORDS]: the sample's vcf/ (HG002's VCF) and, with records,
# its mito/ chrM call set as step 20 leaves it.
sample() {
  local p=$1
  mkdir -p "${G}/${p}/vcf"
  cp "$SRC" "${G}/${p}/vcf/${p}.vcf.gz"
  cp "${SRC}.tbi" "${G}/${p}/vcf/${p}.vcf.gz.tbi"
  [ $# -gt 1 ] || return 0
  mkdir -p "${G}/${p}/mito"
  printf '%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n%s\n' "$HEADER" "$p" "$2" \
    > "${G}/${p}/mito/${p}_chrM.vcf"
  GENOME_DIR="$G" bcf view -Oz -o "${p}/mito/${p}_chrM_filtered.vcf.gz" "${p}/mito/${p}_chrM.vcf"
  rm -f "${G}/${p}/mito/${p}_chrM.vcf"
  check_eq "${p}: the chrM call set holds its records" \
    "$(GENOME_DIR="$G" vcf_count "${p}/mito/${p}_chrM_filtered.vcf.gz")" "$(grep -c . <<< "$2")"
}

# mt REPORT: PharmCAT's callSource and diplotype label for MT-RNR1, then the
# drugs with a recommendation for that label, one per line.
mt() {
  python3 - "$1" <<'PY'
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception as e:
    print(f"unreadable report.json: {e}")
    sys.exit(0)
g = (d.get("genes") or {}).get("MT-RNR1") or {}
dips = g.get("sourceDiplotypes") or g.get("recommendationDiplotypes") or [{}]
label = dips[0].get("label") or (dips[0].get("allele1") or {}).get("name") or "?"
print(f"{g.get('callSource', 'absent')}\t{label}")
for src, by_drug in (d.get("drugs") or {}).items():
    for name, rep in (by_drug or {}).items():
        for gl in rep.get("guidelines") or []:
            for ann in gl.get("annotations") or []:
                labels = {dp.get("label") for gt in ann.get("genotypes") or []
                          for dp in gt.get("diplotypes") or [] if dp.get("gene") == "MT-RNR1"}
                if label in labels and " ".join(str(ann.get("drugRecommendation") or "").split()):
                    print(f"drug\t{rep.get('name') or name}\t{src}")
PY
}
denied() { grep -c 'pharmcat\.log (Permission denied)' <<< "$1" || true; }

sample MTHOM "${BACKGROUND}
$(row 1555 A G PASS 0.990)"
sample MTREF "$BACKGROUND"
sample MTHET "${BACKGROUND}
$(row 1555 A G PASS 0.400)"
sample MTNONE

declare -A LOGS=()
for p in MTHOM MTREF MTHET MTNONE; do
  GENOME_DIR="$G" run_step 07-pharmacogenomics.sh "$p"
  check_step_exit "07-pharmacogenomics.sh (${p})"
  LOGS[$p]=$(cat "$STEP_LOG")
  check_eq "${p}: no 'pharmcat.log (Permission denied)' line in the step log" "$(denied "${LOGS[$p]}")" 0
  check "${p}: PharmCAT wrote report.json" test -s "${G}/${p}/vcf/${p}.report.json"
  check_eq "${p}: the merged outside-call file is deleted" \
    "$(find "${G}/${p}/vcf" -name '*.pharmcat_outside_calls.tsv' | wc -l | tr -d ' ')" 0
done

R=$(mt "${G}/MTHOM/vcf/MTHOM.report.json")
echo "MTHOM MT-RNR1: ${R}"
check "MTHOM: step 07 read the chrM calls" has '^chrM calls \(step 20\): ' "${LOGS[MTHOM]}"
check "MTHOM: the merge passed m.1555A>G, homoplasmic" \
  has '^pgx_mito_merge: MT-RNR1 m\.1555A>G passed to PharmCAT as an outside call \(homoplasmic, chrM:1555 A>G AF 0\.990' "${LOGS[MTHOM]}"
check "MTHOM: PharmCAT got the outside calls" has '^Outside calls \(step 36, and MT-RNR1 from step 20\): ' "${LOGS[MTHOM]}"
check_eq "MTHOM: PharmCAT reports MT-RNR1 m.1555A>G as an outside call" "$(head -n 1 <<< "$R")" $'OUTSIDE\tm.1555A>G'
check_ge "MTHOM: drugs with a recommendation for MT-RNR1 m.1555A>G" "$(grep -c '^drug' <<< "$R" || true)" 1
check "MTHOM: amikacin is one of them" has $'^drug\tamikacin\t' "$R"

R=$(mt "${G}/MTREF/vcf/MTREF.report.json")
echo "MTREF MT-RNR1: ${R}"
check "MTREF: the merge passed Reference" has '^pgx_mito_merge: MT-RNR1 Reference passed to PharmCAT' "${LOGS[MTREF]}"
check_eq "MTREF: PharmCAT reports MT-RNR1 Reference as an outside call" "$(head -n 1 <<< "$R")" $'OUTSIDE\tReference'

R=$(mt "${G}/MTHET/vcf/MTHET.report.json")
echo "MTHET MT-RNR1: ${R}"
check "MTHET: the merge holds m.1555A>G back as heteroplasmic" \
  has '^pgx_mito_merge: held back: m\.1555A>G \(chrM:1555 A>G AF 0\.400\) is heteroplasmic' "${LOGS[MTHET]}"
check "MTHET: and says MT-RNR1 stays No Result" has '^pgx_mito_merge: no MT-RNR1 call passed to PharmCAT: it stays No Result' "${LOGS[MTHET]}"
check "MTHET: PharmCAT got no outside call" has '^No outside calls ' "${LOGS[MTHET]}"
check "MTHET: PharmCAT has no MT-RNR1 call" has '^NONE' "$(head -n 1 <<< "$R")"

R=$(mt "${G}/MTNONE/vcf/MTNONE.report.json")
echo "MTNONE MT-RNR1 (control): ${R}"
check "MTNONE: the step says there are no chrM calls" \
  has '^No chrM calls \(step 20 not run, no .*MTNONE_chrM_filtered\.vcf\.gz\): MT-RNR1 stays No Result\.' "${LOGS[MTNONE]}"
check "MTNONE: no merge ran" lacks '^pgx_mito_merge:' "${LOGS[MTNONE]}"
check "MTNONE: PharmCAT has no MT-RNR1 call" has '^NONE' "$(head -n 1 <<< "$R")"

# --- control: the jar as step 07 ran it before (image working directory) --------
PRE="${G}/MTHOM/vcf/MTHOM.preprocessed.vcf.bgz"
check "MTHOM's preprocessed VCF is there for the control" test -s "$PRE"
CTL="${CASE_TMP}/control"
rm -rf "$CTL"
mkdir -p "$CTL"
# The header rewrite of step 07, then the jar with no -w.
HEADER_FIX='/^##/ { gsub(/\\"/, "\047"); gsub(/\\/, "/") } { print }'
# shellcheck disable=SC2016  # $1 and $2 belong to the inner sh
docker run --rm --network none --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "${G}/MTHOM/vcf:/data:ro" -v "${CTL}:/ctl" "$PHARMCAT_IMAGE" \
  sh -c 'gzip -dc "$1" | awk "$2" > /ctl/in.vcf && java -jar /pharmcat/pharmcat.jar -vcf /ctl/in.vcf -o /ctl -bf control -matcher' \
  sh /data/MTHOM.preprocessed.vcf.bgz "$HEADER_FIX" > "${CTL}/log.txt" 2>&1
CTL_RC=$?
CTL_LOG=$(cat "${CTL}/log.txt")
grep -m 3 'pharmcat\.log' <<< "$CTL_LOG"
check_eq "control: the jar runs" "$CTL_RC" 0
check_ge "control: run as before, logback prints 'pharmcat.log (Permission denied)'" "$(denied "$CTL_LOG")" 1

finish
