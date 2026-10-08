#!/usr/bin/env bash
# Step 12 after step 20: the haplogroup comes from the Mutect2 chrM calls, and
# haplocheck's contamination status reaches the text report. Case 33 ran
# step 12 before step 20 existed, from the DeepVariant chrM records; the two
# inputs must give the same top-level haplogroup.
. "$(dirname "$0")/lib.sh"

M="${GENOME_DIR}/${SAMPLE}/mito"
check "step 20's Mutect2 calls are there (cases 36 and bash-step-20)" test -s "${M}/${SAMPLE}_chrM_filtered.vcf.gz"
DV_HG=$(awk -F'\t' 'NR == 2 {gsub(/"/, "", $2); print $2}' "${M}/${SAMPLE}_haplogroup.txt" 2>/dev/null)
echo "Haplogroup from the DeepVariant chrM records (case 33): ${DV_HG:-none}"

run_step 12-mito-haplogroup.sh "$SAMPLE"
check_step_exit 12-mito-haplogroup.sh
check "the log names the Mutect2 calls as the input" has '^Input: Mutect2 chrM calls of step 20' "$(cat "$STEP_LOG")"
HG=$(awk -F'\t' 'NR == 2 {gsub(/"/, "", $2); print $2}' "${M}/${SAMPLE}_haplogroup.txt" 2>/dev/null)
check "haplogroup from the Mutect2 calls (${HG:-empty})" test -n "${HG:-}"
check_eq "the same top-level haplogroup as from the DeepVariant records" "${HG:0:1}" "${DV_HG:0:1}"
CHECK="${M}/${SAMPLE}_haplocheck.txt"
head -n 2 "$CHECK" 2>/dev/null
STATUS=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) {h = $i; gsub(/"/, "", h); if (h == "Contamination Status") c = i}; next}
  NR == 2 && c {v = $c; gsub(/"/, "", v); print v}' "$CHECK" 2>/dev/null)
check "haplocheck gives a contamination status (${STATUS:-none})" has '^(YES|NO|ND)$' "${STATUS:-}"
check "the log prints it" has '^Contamination \(haplocheck\): (YES|NO|ND) ' "$(cat "$STEP_LOG")"

run_step generate-report.sh "$SAMPLE"
check_step_exit generate-report.sh
TXT="${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt"
check "the text report has the haplocheck status line" grep -qE '^  Contamination \(haplocheck\): (YES|no|not determined)' "$TXT"

finish
