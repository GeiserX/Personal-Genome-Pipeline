#!/usr/bin/env bash
# The multi-SNP rules of step 42 that came with APOL1 and SERPINA1, end to end:
#   - scripts/ci/check-rule-alleles.sh passes the rules as shipped; controls:
#     a copy with SERPINA1 Z (rs28929474) written ref and alt the wrong way
#     round, and APOL1 G1 S342G (rs73885319) written as its minus-strand bases,
#     fails, naming each with its reason;
#   - scripts/42-traits.sh, on a copy of the HG002 VCF and gVCF that case 21
#     wrote, writes one rules line per rule, the APOL1 and SERPINA1 lines among them;
#   - bin/render_report.py on that folder puts both rules in the HTML genotype
#     rules table, each with its trait, gene and sites; control: the same
#     folder with the SERPINA1 line taken out of the rules file renders a table
#     without it, and still with APOL1.
# Writes one sample folder under GENOME_DIR, removed at the end; no other case
# reads it.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
RULES="${REPO}/assets/genotype_rules.tsv"
CHECK="${REPO}/scripts/ci/check-rule-alleles.sh"

# --- the orientation check on the new rows -----------------------------------------
bash "$CHECK" > "${CASE_TMP}/rules.log" 2>&1
check_eq "check-rule-alleles.sh on the shipped rules exits 0" "$?" 0
check "it checked the SERPINA1 Z site" has '^ok: rs28929474 chr14:94378610 C>T$' "$(cat "${CASE_TMP}/rules.log")"
check "it checked the APOL1 G1 site" has '^ok: rs73885319 chr22:36265860 A>G$' "$(cat "${CASE_TMP}/rules.log")"
awk -F'\t' -v OFS='\t' '$1 == "site" && $3 == "rs28929474" { t = $6; $6 = $7; $7 = t }
  $1 == "site" && $3 == "rs73885319" { $6 = "T"; $7 = "C" } { print }' "$RULES" > "${CASE_TMP}/flipped.tsv"
check_eq "control: the flipped copy differs from the rules in two lines" \
  "$(diff "$RULES" "${CASE_TMP}/flipped.tsv" | grep -c '^>' || true)" 2
bash "$CHECK" "${CASE_TMP}/flipped.tsv" > "${CASE_TMP}/flipped.log" 2>&1
FRC=$?
FLOG=$(cat "${CASE_TMP}/flipped.log")
check "control: the flipped copy fails (exit ${FRC})" test "$FRC" -ne 0
check "control: SERPINA1 Z named as swapped" has '^FAIL: rs28929474: ref and alt are swapped' "$FLOG"
check "control: APOL1 G1 named as minus-strand bases" has '^FAIL: rs73885319: the alleles are the minus-strand bases' "$FLOG"

# --- step 42 and the HTML rules table -------------------------------------------------
S="${SAMPLE}R"
rm -rf "${G:?}/${S}"
mkdir -p "${G}/${S}/vcf"
for f in vcf.gz vcf.gz.tbi g.vcf.gz g.vcf.gz.tbi; do
  cp "${G}/${SAMPLE}/vcf/${SAMPLE}.${f}" "${G}/${S}/vcf/${S}.${f}"
done
run_step 42-traits.sh "$S"
check_step_exit 42-traits.sh
RT="${G}/${S}/traits/${S}_genotype_rules.tsv"
check_eq "one rules line per rule" "$(awk 'NR > 1' "$RT" 2>/dev/null | wc -l | tr -d ' ')" \
  "$(awk -F'\t' '$1 == "rule"' "$RULES" | wc -l | tr -d ' ')"
check_eq "the rules file has the APOL1 and SERPINA1 lines" \
  "$(awk -F'\t' '$1 == "APOL1" || $1 == "SERPINA1" { print $1 }' "$RT" 2>/dev/null | sort | tr '\n' ' ')" "APOL1 SERPINA1 "
awk -F'\t' '$1 == "APOL1" || $1 == "SERPINA1"' "$RT" 2>/dev/null | cut -f1-6 | sed 's/^/  /'

render() {  # render OUT.html: the HTML report of the folder
  python3 "${REPO}/bin/render_report.py" --sample "$S" --sample-dir "${G}/${S}" --json "${CASE_TMP}/summary.json" \
    -o "${CASE_TMP}/report.txt" -o "$1" > /dev/null 2>"${CASE_TMP}/render.err"
}
render "${CASE_TMP}/report.html"
HTML=$(cat "${CASE_TMP}/report.html" 2>/dev/null)
check "HTML report: the genotype rules table" has '<h3>Genotype rules \(several sites read together\)</h3>' "$HTML"
check "HTML report: the SERPINA1 row with its trait and sites" \
  has '<td>Alpha-1 antitrypsin Pi type \(M, S, Z\)</td><td>SERPINA1<br><span[^>]*>rs17580, rs28929474</span>' "$HTML"
check "HTML report: the APOL1 row with its trait and sites" \
  has '<td>APOL1 kidney risk variants \(G1 read, G2 not read\)</td><td>APOL1<br><span[^>]*>rs73885319, rs60910145</span>' "$HTML"

# Each row's Genotypes and Result cells are the rules file's: its genotypes, and
# its result when called, else its status (the report adds how it was phased).
cell() {  # cell RULE COLUMN
  awk -F'\t' -v r="$1" -v k="$2" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    $1 == r { print $c[k]; exit }' "$RT" 2>/dev/null
}
for r in SERPINA1 APOL1; do
  geno=$(cell "$r" genotypes | sed 's/,/, /g')
  if [ "$(cell "$r" status)" = called ]; then res=$(cell "$r" result); else res=$(cell "$r" status | tr _ ' '); fi
  row="<td>${r}<br><span style=\"font-size:12px;color:#666\">$(cell "$r" rsids | sed 's/,/, /g')</span></td><td>${geno}</td><td>${res}"
  check "HTML report: the ${r} row's genotypes (${geno}) and result (${res}) are the rules file's" \
    grep -qF -- "$row" "${CASE_TMP}/report.html"
done

awk -F'\t' '$1 != "SERPINA1"' "$RT" > "${CASE_TMP}/rules-no-serpina1.tsv" && cp "${CASE_TMP}/rules-no-serpina1.tsv" "$RT"
render "${CASE_TMP}/report-control.html"
CHTML=$(cat "${CASE_TMP}/report-control.html" 2>/dev/null)
check "control: without its rules line, no SERPINA1 row" lacks '<td>SERPINA1<br>' "$CHTML"
check "control: the APOL1 row is still there" has '<td>APOL1<br>' "$CHTML"

in_genome "$BCFTOOLS_IMAGE" rm -rf "$S"
finish
