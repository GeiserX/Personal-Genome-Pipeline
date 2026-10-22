#!/usr/bin/env bash
# Do a score's chrX rows survive pgsc_calc's ancestry-panel intersection?
# Step 25 on HG002's gVCF (case 21), declared male, with the PGS Catalog's
# synthetic panel (GRCh38_HAPNEST_reference, as case prs-2): every score of
# the list is panel SNVs inside the fixture's slices plus three chrX rows at
# hom-ref positions of the gVCF. With a panel pgsc_calc scores only the
# variants it finds in the panel too, so the case reads pgsc_calc's match log
# for PGS000662 (prostate cancer, which has chrX rows in the real score) and
# checks that the summary's ChrX column says what the log shows: "3 scored"
# when the three chrX rows were matched, "3 left out (not in the ancestry
# panel)" when none were, both counts when some were, and Variants_Matched
# counts only the matched ones. It also counts the panel's own chrX variants and writes both to
# the notes.
#   - the declared sex with the panel: PGS000004 (breast cancer) gets no
#     percentile ("not applicable to the declared sex"), PGS000018 gets one;
#   - control: the same summary read as a no-panel run (prs-table without
#     --panel) says "3 scored" whatever the log says, so the check above can
#     tell the two apart.
# Removes what it made, the panel too.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
D="${GENOME_DIR}/${SAMPLE}"
REF="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta"
PANEL_NAME=GRCh38_HAPNEST_reference
PANEL="${GENOME_DIR}/reference/pgsc_calc/${PANEL_NAME}.tar.zst"
PSITES="${PANEL%.tar.zst}_GRCh38_sites.tsv"

echo "+ ANCESTRY_PANEL_NAME=${PANEL_NAME} scripts/setup.sh --ancestry-panel"
ANCESTRY_PANEL_NAME="$PANEL_NAME" "${REPO}/scripts/setup.sh" --ancestry-panel "$GENOME_DIR" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --ancestry-panel exits 0" "${PIPESTATUS[0]}" 0
check_ge "GRCh38 SNVs in the panel's site list" "$(wc -l < "$PSITES" 2>/dev/null | tr -d ' ')" 1000
# The panel's own chrX variants (its .pvar or .pvar.zst files: chromosome X or 23).
PVARS=$(tar --zstd -tf "$PANEL" 2>/dev/null | grep -E '\.pvar(\.zst)?$' | tr '\n' ' ')
echo "variant files in the panel: ${PVARS:-none}"
PANEL_X=unknown
if [ -n "$PVARS" ]; then
  # shellcheck disable=SC2086  # the member names have no spaces
  PANEL_X=$(tar --zstd -xOf "$PANEL" $PVARS 2>/dev/null | zstd -dcf 2>/dev/null \
    | awk '!/^#/ && ($1 == "X" || $1 == "23" || $1 == "chrX") {n++} END {print n + 0}')
fi
echo "chrX variants in the panel's variant files: ${PANEL_X}"

# --- scores: panel SNVs in the slices, and three chrX hom-ref rows ------------------
awk -F'\t' -v OFS='\t' 'NR == FNR {s[$1] = s[$1] " " $2 + 1 ":" $3; next}
  ($1 in s) {n = split(s[$1], r, " "); for (i = 1; i <= n; i++) {split(r[i], b, ":"); if ($2 >= b[1] && $2 <= b[2]) {print; break}}}' \
  "${FIXTURE_DIR}/regions.bed" "$PSITES" \
  | awk -F'\t' '($3 $4) !~ /^(AT|TA|CG|GC)$/' | awk 'NR % 5 == 1' | head -n 300 > "${CASE_TMP}/score_sites.tsv"
N_SCORE=$(wc -l < "${CASE_TMP}/score_sites.tsv" | tr -d ' ')
check_ge "panel SNVs in the fixture's slices used as score rows" "$N_SCORE" 20
XR=chrX:73700001-74000000
bcf query -r "$XR" -f '%POS\t%INFO/END\t[%GT]\n' "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" > "${CASE_TMP}/x_blocks.tsv"
bcf query -r "$XR" -i 'GT="alt"' -f '%POS\n' "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" > "${CASE_TMP}/x_vars.tsv"
X_SITES="${CASE_TMP}/x_sites.tsv"
python3 - "${CASE_TMP}/x_blocks.tsv" "${CASE_TMP}/x_vars.tsv" "$REF" > "$X_SITES" <<'PY'
import sys
blocks, variants, fasta = sys.argv[1:4]
var = [int(l) for l in open(variants) if l.strip()]
fai = {l.split()[0]: [int(x) for x in l.split()[1:5]] for l in open(fasta + ".fai")}
length, offset, bases, width = fai["chrX"]
fh = open(fasta, "rb")
def base(pos):
    i = pos - 1
    fh.seek(offset + (i // bases) * width + i % bases)
    return fh.read(1).decode().upper()
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
rows = [l.rstrip("\n").split("\t") for l in open(blocks) if l.strip()]
refs = [(int(r[0]), int(r[1])) for r in rows if len(r) == 3 and r[1] not in ("", ".") and r[2] in ("0/0", "0|0")]
out, last = [], 0
for start, end in refs:
    p = (start + end) // 2
    if end - start >= 50 and p - last > 2000 and all(abs(v - p) > 100 for v in var) and base(p) in comp:
        out.append((p, base(p)))
        last = p
for p, b in out[:3]:
    other = next(x for x in "ACGT" if x not in (b, comp[b]))
    print(f"X\t{p}\t{b}\t{other}\t0.1")
PY
check_eq "chrX hom-ref rows" "$(wc -l < "$X_SITES" | tr -d ' ')" 3
SCORES="${GENOME_DIR}/prs_scores"
rm -rf "$SCORES" "${D}/prs" "${D}/ancestry"
mkdir -p "$SCORES"
for id in $(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${REPO}/assets/pgs_scores.tsv"); do
  { printf '#pgs_id=%s\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n' "$id"
    awk -F'\t' -v OFS='\t' '{sub(/^chr/, "", $1); print $1, $2, $4, $3, ((NR % 7) - 3) / 10}' "${CASE_TMP}/score_sites.tsv"
    cat "$X_SITES"; } | gzip -c > "${SCORES}/${id}.txt.gz"
done

# --- step 25 with the panel, male -----------------------------------------------------------
ANCESTRY_PANEL="$PANEL" run_step 25-prs.sh "$SAMPLE" male
check_step_exit "25-prs.sh male, with ${PANEL_NAME}"
SUMMARY="${D}/prs/${SAMPLE}_prs_summary.tsv"
cat "$SUMMARY" 2>/dev/null
col() { awk -F'\t' -v id="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["PGS_ID"] == id {print $c[k]}' "$1" 2>/dev/null; }
MATCH_LOG=$(find "${D}/prs/pgsc_calc/results" -name '*_log.csv.gz' 2>/dev/null | head -n 1)
X_MATCHED=$(python3 - "${MATCH_LOG:-/dev/null}" <<'PY2'
import csv, gzip, sys
try:
    rows = list(csv.DictReader(gzip.open(sys.argv[1], "rt")))
except (OSError, EOFError):
    print("no match log"); sys.exit()
print(sum(1 for r in rows if r.get("accession") == "PGS000662" and r.get("chr_name") in ("X", "23")
          and r.get("match_status") == "matched"))
PY2
)
CHRX=$(col "$SUMMARY" PGS000662 ChrX)
MATCHED=$(col "$SUMMARY" PGS000662 Variants_Matched)
echo "PGS000662: ChrX '${CHRX}', ${MATCHED} of $(col "$SUMMARY" PGS000662 Variants_Total) matched, ${X_MATCHED} chrX rows matched in the log"
check_eq "PGS000662: every row is in the total (${N_SCORE} panel SNVs + 3 chrX)" "$(col "$SUMMARY" PGS000662 Variants_Total)" "$((N_SCORE + 3))"
check "pgsc_calc's match log was read (${MATCH_LOG:-none})" has '^[0-9]+$' "$X_MATCHED"
case "$X_MATCHED" in
  3) check_eq "chrX rows matched through the panel: the summary says scored" "$CHRX" "3 scored" ;;
  0) check_eq "chrX rows not in the panel: the summary says so" "$CHRX" "3 left out (not in the ancestry panel)" ;;
  1|2) check_eq "some chrX rows in the panel: the summary names both counts" "$CHRX" \
         "${X_MATCHED} scored, $((3 - X_MATCHED)) left out (not in the ancestry panel)" ;;
  *) fail "chrX rows matched in the log: ${X_MATCHED}, expected 0 to 3" ;;
esac
# What docs/25-prs.md says of this panel: none of the chrX rows is scored.
check_eq "${PANEL_NAME}: no chrX row of PGS000662 is matched (as docs/25-prs.md says)" "$X_MATCHED" 0
AUTO_MATCHED=$(python3 - "${MATCH_LOG:-/dev/null}" <<'PY3'
import csv, gzip, sys
try:
    rows = list(csv.DictReader(gzip.open(sys.argv[1], "rt")))
except (OSError, EOFError):
    print(-1); sys.exit()
print(sum(1 for r in rows if r.get("accession") == "PGS000662" and r.get("chr_name") not in ("X", "23")
          and r.get("match_status") == "matched"))
PY3
)
check_eq "Variants_Matched is the log's matched rows, chrX included only when matched" "$MATCHED" "$((AUTO_MATCHED + X_MATCHED))"
PCT=$(col "$SUMMARY" PGS000018 Percentile)
check "male with the panel: PGS000018 has a percentile (${PCT:-none})" awk -v p="${PCT:-x}" 'BEGIN {exit !(p ~ /^[0-9.]+$/)}'
check_eq "male with the panel: PGS000004 has no percentile" "$(col "$SUMMARY" PGS000004 Percentile)" NA
check_eq "male with the panel: PGS000004 is not applicable" "$(col "$SUMMARY" PGS000004 Sex_Check)" "not applicable to the declared sex"
# Control: the same pgsc_calc output read as a run without the panel.
docker run --rm --network none --user "$(id -u):$(id -g)" -v "${D}/prs:/p" -v "${REPO}/bin:/b:ro" "$PYTHON_IMAGE" \
  python3 /b/collect_summary.py prs-table --sample "$SAMPLE" --results /p/pgsc_calc/results --sampleset sample \
    --scores /p/pgsc_calc/scores --input-kind gvcf --sex male --out /p/control.tsv > /dev/null 2>&1
check_eq "control: read without the panel the ChrX column says 3 scored" "$(col "${D}/prs/control.tsv" PGS000662 ChrX)" "3 scored"
printf '#### PRS: chrX rows through the ancestry panel (%s)\n\nchrX variants in the panel: %s; PGS000662 chrX rows matched %s of 3; ChrX column "%s"\n\n' \
  "$PANEL_NAME" "${PANEL_X:-unknown}" "${X_MATCHED:-none}" "${CHRX:-none}" >> "$E2E_NOTES"

rm -rf "$SCORES" "${D:?}/prs" "${D:?}/ancestry" "$PANEL" "$PSITES"
finish
