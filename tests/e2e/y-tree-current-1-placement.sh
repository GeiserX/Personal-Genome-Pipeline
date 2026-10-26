#!/usr/bin/env bash
# Step 37 places the male copy of HG002 on the current YFull YTree beside
# Yleaf's own tree, from the fixture's six chrY stretches:
#   - setup.sh --yleaf-data installs the tree and its two position tables
#     (pinned commits, sha256) in <yleaf data>/ytree;
#   - the step writes <sample>_y_tree.txt naming the tree, and a per-level
#     table whose rows are the path from the root to the terminal branch;
#     the path runs through CT, as it does for every man outside haplogroups
#     A and B (HG002's line is one of them);
#   - the Nextflow run of case nextflow-hardening-2 and the bash step 37 in
#     its GENOME_DIR wrote the same two files.
# Red controls: the same judgement on a copy of the tables with the ancestral
# and derived bases swapped must fail (the walk then leaves CT), and a
# minimum depth no read reaches must give Hg NA.
. "$(dirname "$0")/lib.sh"

YDATA="${GENOME_DIR}/reference/yleaf-${YLEAF_DATA_VERSION}/data"
YT="${YDATA}/ytree"
check "Yleaf's tables and the current tree install" "${REPO}/scripts/setup.sh" --yleaf-data "$GENOME_DIR"
check "the tree folder names YFull YTree ${YTREE_VERSION}" \
  grep -qxF "$(printf 'tree\tYFull YTree %s (CC BY 4.0)' "$YTREE_VERSION")" "${YT}/version.txt"
python3 "${REPO}/bin/y_tree_place.py" positions --tree-dir "$YT" "${CASE_TMP}/positions.txt" > "${CASE_TMP}/positions.log" 2>&1
cat "${CASE_TMP}/positions.log"
check_ge "tree SNP positions the tables place" "$(grep -c . "${CASE_TMP}/positions.txt" 2>/dev/null || true)" 300000

# A male copy of the sample (the BAM hard-linked), as in the Yleaf case.
M="${SAMPLE}yt"
mkdir -p "${GENOME_DIR}/${M}/aligned" "${GENOME_DIR}/${M}/indexcov"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${M}/aligned/${M}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${M}/aligned/${M}_sorted.bam.bai"
printf '#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\n%s\t%s\t-9\t-9\t1\t-9\t1.0\t1.0\n' \
  "$M" "$M" > "${GENOME_DIR}/${M}/indexcov/indexcov-indexcov.ped"

run_step 37-y-haplogroup.sh "$M"
check_step_exit 37-y-haplogroup.sh
OUT="${GENOME_DIR}/${M}/y_haplogroup"
cat "${OUT}/${M}_y_tree.txt" "${OUT}/${M}_y_tree_path.tsv" 2>/dev/null

# col FILE NAME: the value of column NAME in the one data row of FILE
col() { awk -F'\t' -v n="$2" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} NR == 2 {print $c[n]}' "$1" 2>/dev/null; }

# judge PREFIX: 0 when PREFIX_y_tree.txt and PREFIX_y_tree_path.tsv hold a
# placement this sample can have: a terminal branch at level 10 or deeper,
# at least 10 derived SNPs on the path, a path table of exactly Level rows
# that runs through CT. Each failed condition is printed.
judge() {
  local t="${1}_y_tree.txt" p="${1}_y_tree_path.tsv" ok=0 hg lvl pd rows
  hg=$(col "$t" Hg); lvl=$(col "$t" Level); pd=$(col "$t" Path_derived)
  rows=$(awk 'NR > 1' "$p" 2>/dev/null | grep -c . || true)
  [ -n "$hg" ] && [ "$hg" != NA ] || { echo "  judge: no terminal branch (${hg:-none})"; ok=1; }
  [[ "$lvl" =~ ^[0-9]+$ ]] && [ "$lvl" -ge 10 ] || { echo "  judge: level ${lvl:-none} < 10"; ok=1; }
  [[ "$pd" =~ ^[0-9]+$ ]] && [ "$pd" -ge 10 ] || { echo "  judge: ${pd:-no} derived SNPs on the path"; ok=1; }
  [ "$rows" = "$lvl" ] || { echo "  judge: ${rows} path rows for level ${lvl}"; ok=1; }
  awk -F'\t' 'NR > 1 && $2 == "CT" {f = 1} END {exit !f}' "$p" 2>/dev/null || { echo "  judge: the path skips CT"; ok=1; }
  return $ok
}

check "the result names the tree and the tables it used" \
  test "$(col "${OUT}/${M}_y_tree.txt" Tree)|$(col "${OUT}/${M}_y_tree.txt" Positions)" = \
       "YFull YTree ${YTREE_VERSION} (CC BY 4.0)|Yleaf ${YTREE_POSITIONS_VERSION} GRCh38 SNP tables (GPL-3.0)"
check "the placement on the current tree holds (terminal branch, level, path through CT)" judge "${OUT}/${M}"
check "the step prints the current-tree line" \
  has "^Current tree \\(YFull YTree ${YTREE_VERSION//./\\.} \\(CC BY 4\\.0\\)\\): [A-Z]" "$(cat "$STEP_LOG")"
check "and still prints Yleaf's line, with the tree it comes from" \
  has '^Y haplogroup: ' "$(cat "$STEP_LOG")"
check "Yleaf's own tree is named" has "^  \\(Yleaf ${YLEAF_DATA_VERSION//./\\.}, on its own tree, YFull v10\\.01\\)" "$(cat "$STEP_LOG")"
check "the step removed its tree pileup" test ! -e "${OUT}/tree_pileup.txt"

# --- Red controls: the same pileup, judged on swapped tables and at a depth no read reaches ---
mkdir -p "${GENOME_DIR}/${M}/tree-control"
cp "${CASE_TMP}/positions.txt" "${GENOME_DIR}/${M}/tree-control/positions.txt"
sam mpileup -l "${M}/tree-control/positions.txt" -AQ20q1 "${M}/aligned/${M}_sorted.bam" > "${CASE_TMP}/tree_pileup.txt" 2>/dev/null
check_ge "pileup rows at the tree's positions" "$(grep -c . "${CASE_TMP}/tree_pileup.txt" || true)" 1000
python3 "${REPO}/bin/y_tree_place.py" place --tree-dir "$YT" --pileup "${CASE_TMP}/tree_pileup.txt" \
  --sample "$M" --out-prefix "${CASE_TMP}/same"
check "the placement from the case's own pileup is the step's" cmp -s "${CASE_TMP}/same_y_tree.txt" "${OUT}/${M}_y_tree.txt"
SW="${CASE_TMP}/swapped"
mkdir -p "$SW"
cp "${YT}/current_tree.json" "${YT}/version.txt" "$SW/"
for f in positions_ftdna_hg38.txt positions_yfull_hg38.txt; do
  awk -F'\t' -v OFS='\t' '{t = $6; $6 = $7; $7 = t; print}' "${YT}/${f}" > "${SW}/${f}"
done
python3 "${REPO}/bin/y_tree_place.py" place --tree-dir "$SW" --pileup "${CASE_TMP}/tree_pileup.txt" \
  --sample "$M" --out-prefix "${CASE_TMP}/swapped"
not() { ! "$@"; }
check "red control: the judgement fails with ancestral and derived swapped" not judge "${CASE_TMP}/swapped"
python3 "${REPO}/bin/y_tree_place.py" place --tree-dir "$YT" --pileup "${CASE_TMP}/tree_pileup.txt" \
  --sample "$M" --out-prefix "${CASE_TMP}/deep" --min-depth 100000
check_eq "red control: no SNP is called at depth 100000" "$(col "${CASE_TMP}/deep_y_tree.txt" Hg)|$(col "${CASE_TMP}/deep_y_tree.txt" SNPs_called)" "NA|0"
rm -rf "${GENOME_DIR:?}/${M}/tree-control" "${SW:?}"

# --- The Nextflow module wrote what the bash step wrote (case nextflow-hardening-2) ---
NFY="${GENOME_DIR}/nf-hardening/${SAMPLE}/y_haplogroup"
BY="${E2E_WORK}/genome-hardening/${SAMPLE}/y_haplogroup"
for f in "${SAMPLE}_y_tree.txt" "${SAMPLE}_y_tree_path.tsv"; do
  check "Nextflow wrote ${f}" test -s "${NFY}/${f}"
  check "the bash step wrote ${f} on the same BAM" test -s "${BY}/${f}"
  check "and the two are the same" cmp -s "${NFY}/${f}" "${BY}/${f}"
done

finish
