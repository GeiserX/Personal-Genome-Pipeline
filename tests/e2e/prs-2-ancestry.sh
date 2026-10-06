#!/usr/bin/env bash
# Step 26 and step 25's percentiles with an ancestry panel installed (cases 21
# and 03 left HG002's VCF and gVCF on the fixture's slices):
#   - without a panel, step 26 says so in one line and exits 0;
#   - setup.sh --ancestry-panel installs a panel and its GRCh38 site list.
#     On a pull request the panel is the PGS Catalog's small synthetic one
#     (GRCh38_HAPNEST_reference, 268 MB), so the run fits the job; the 1000
#     Genomes panel users install is measured in case prs-3 (monthly and
#     dispatched runs only, see docs/25-prs.md);
#   - step 26 projects HG002 onto it: the ancestry table has a population
#     label, its probabilities and the principal components, and step 25's
#     summary has a percentile for the score, with that population as its group;
#   - both reports show the percentile with its group, and not the raw-score line.
# The score is synthetic: panel SNVs inside the fixture's slices (so the panel's
# own samples carry them too), effect allele = the panel's ALT, no
# strand-ambiguous pair (pgsc_calc drops those).
# PANEL_NAME (from prs-3) picks another panel. Removes what it made.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
D="${GENOME_DIR}/${SAMPLE}"
PANEL_NAME=${PANEL_NAME:-GRCh38_HAPNEST_reference}
PANEL="${GENOME_DIR}/reference/pgsc_calc/${PANEL_NAME}.tar.zst"
SITES="${PANEL%.tar.zst}_GRCh38_sites.tsv"

# --- no panel: one line, exit 0 ---------------------------------------------------
ANCESTRY_PANEL=none run_step 26-ancestry.sh "$SAMPLE"
check_step_exit "26-ancestry.sh (no panel)"
check_eq "step 26 without a panel prints one line" "$(grep -c . "$STEP_LOG")" 1
check "that line says it was skipped and how to install the panel" has '^Step 26 skipped: .*setup\.sh --ancestry-panel' "$(cat "$STEP_LOG")"
check "step 26 without a panel wrote no ancestry table" test ! -e "${D}/ancestry/${SAMPLE}_ancestry.tsv"

# --- the panel ------------------------------------------------------------------
echo "+ ANCESTRY_PANEL_NAME=${PANEL_NAME} scripts/setup.sh --ancestry-panel"
ANCESTRY_PANEL_NAME="$PANEL_NAME" "${REPO}/scripts/setup.sh" --ancestry-panel "$GENOME_DIR" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --ancestry-panel exits 0" "${PIPESTATUS[0]}" 0
check "the panel is installed" test -s "$PANEL"
check_ge "GRCh38 SNVs in the panel's site list" "$(wc -l < "$SITES" 2>/dev/null | tr -d ' ')" 1000
check "the site list is chrN, position, REF, ALT" \
  awk -F'\t' 'NR > 1000 {exit} !($1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[ACGT]$/ && $4 ~ /^[ACGT]$/) {bad = 1} END {exit bad}' "$SITES"

# --- a score of panel SNVs inside the fixture's slices ---------------------------
SCORES="${GENOME_DIR}/prs_scores"
rm -rf "$SCORES" "${D}/prs" "${D}/ancestry"
mkdir -p "$SCORES"
# regions.bed is 0-based; the site list is chr-prefixed, sorted by chrom then position.
awk -F'\t' -v OFS='\t' 'NR == FNR {s[$1] = s[$1] " " $2 + 1 ":" $3; next}
  ($1 in s) {n = split(s[$1], r, " "); for (i = 1; i <= n; i++) {split(r[i], b, ":"); if ($2 >= b[1] && $2 <= b[2]) {print; break}}}' \
  "${FIXTURE_DIR}/regions.bed" "$SITES" \
  | awk -F'\t' '($3 $4) !~ /^(AT|TA|CG|GC)$/' | awk 'NR % 5 == 1' | head -n 300 > "${CASE_TMP}/score_sites.tsv"
N_SCORE=$(wc -l < "${CASE_TMP}/score_sites.tsv" | tr -d ' ')
check_ge "panel SNVs in the fixture's slices used as score sites" "$N_SCORE" 20
for id in $(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${REPO}/assets/pgs_scores.tsv"); do
  { printf '#pgs_id=%s\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n' "$id"
    awk -F'\t' -v OFS='\t' '{sub(/^chr/, "", $1); print $1, $2, $4, $3, ((NR % 7) - 3) / 10}' "${CASE_TMP}/score_sites.tsv"; } \
    | gzip -c > "${SCORES}/${id}.txt.gz"
done

# --- step 26 with the panel (it runs step 25) -------------------------------------
ANCESTRY_PANEL="$PANEL" run_step 26-ancestry.sh "$SAMPLE"
check_step_exit "26-ancestry.sh (with ${PANEL_NAME})"
LOG26=$(cat "$STEP_LOG")
ANC="${D}/ancestry/${SAMPLE}_ancestry.tsv"
kv() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "$ANC" 2>/dev/null; }
POP=$(kv population)
echo "population: ${POP:-none}"; cat "$ANC" 2>/dev/null
check "the ancestry table has a population label" has '^[A-Z]{2,4}$' "${POP:-}"
check "it names the panel" test "$(kv reference_panel)" = "$PANEL_NAME"
check_ge "principal components in the table" "$(awk -F'\t' '$1 ~ /^PC[0-9]+$/ && $2 ~ /^-?[0-9.e+-]+$/' "$ANC" 2>/dev/null | wc -l | tr -d ' ')" 5
check "the probability of the label is in the table" test -n "$(kv "probability_${POP:-none}")"
check "step 26 prints the population" has "Population most similar to ${SAMPLE}: ${POP:-none}" "$LOG26"
SUMMARY="${D}/prs/${SAMPLE}_prs_summary.tsv"
col() { awk -F'\t' -v k="$1" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["PGS_ID"] == "PGS000018" {print $c[k]}' "$SUMMARY" 2>/dev/null; }
PCT=$(col Percentile)
check "step 25's summary has a percentile from 0 to 100 (${PCT:-none})" \
  awk -v p="${PCT:-x}" 'BEGIN {exit !(p ~ /^[0-9.]+$/ && p >= 0 && p <= 100)}'
check_eq "its group is the population of the ancestry table" "$(col Ancestry_Group)" "${POP:-none}"
check_eq "input: genotypes from the gVCF (score and panel sites)" "$(col Input)" gvcf
check "no raw-score line in the step log" lacks 'Raw score only' "$LOG26"

# --- the reports ------------------------------------------------------------------
python3 "${REPO}/bin/render_report.py" --sample "$SAMPLE" --sample-dir "$D" --json "${CASE_TMP}/summary.json" \
  -o "${CASE_TMP}/report.txt" -o "${CASE_TMP}/report.html" >/dev/null
check "text report: the percentile with its group" grep -q "percentile ${PCT:-none} (${POP:-none})" "${CASE_TMP}/report.txt"
check "HTML report: the percentile with its group" grep -q "<td>${PCT:-none} (${POP:-none})</td>" "${CASE_TMP}/report.html"
check "neither report has the raw-score line" bash -c '! grep -q "Raw score only" "$1" "$2"' _ "${CASE_TMP}/report.txt" "${CASE_TMP}/report.html"
printf '#### Ancestry with %s (%s score sites in the fixture slices)\n\npopulation %s, PC1 %s, PRS percentile %s\n\n' \
  "$PANEL_NAME" "$N_SCORE" "${POP:-none}" "$(kv PC1)" "${PCT:-none}" >> "$E2E_NOTES"

rm -rf "$SCORES" "${D:?}/prs" "${D:?}/ancestry" "${PANEL}" "${SITES}"
finish
