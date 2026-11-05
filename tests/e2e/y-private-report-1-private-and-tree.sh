#!/usr/bin/env bash
# Step 37's two additions, on the fixture's chrY stretches and the real
# images and tables:
#   1. the private-SNP list (bin/y_private_snps.py): the chrY calls of the
#      sample's own VCF that no public SNP table names. The list must drop
#      calls the tables do name (the header counts them, and at least one
#      real call is dropped, so the subtraction is doing work on real data),
#      and no listed position may be in either table.
#      Only calls inside the chrY callable mask are listed.
#      Red controls, on a copy of the VCF seeded with two synthetic chrY
#      calls at positions no table names (so they do not depend on what the
#      fixture happens to call): the one in the mask is listed and the one
#      outside it is not; adding the listed one to a copy of the tables
#      removes that row and nothing else; tables with no row stop the tool.
#   2. the placement on the current tree in both reports. The bash report
#      (steps 37 and 24) and the Nextflow report of case
#      nextflow-hardening-2-nf both get a "Y Placement on the Current Public
#      Tree" section whose branch is the one in <id>_y_tree.txt.
#      Red control: the same report code over a copy of the Nextflow results
#      without that file leaves the section missing and the card out, which
#      is what a missing main.nf join or html_report arm would give.
# What this cannot show: whether a private SNP is real. That needs another
# man who carries it, so the list is a starting point, never a result.
. "$(dirname "$0")/lib.sh"
# A real tab for the patterns: has and lacks run grep -E, and GNU grep reads \t
# there as a plain t, so a \t pattern never matches a TSV line on Linux.
T=$'\t'

YDATA="${GENOME_DIR}/reference/yleaf-${YLEAF_DATA_VERSION}/data"
YT="${YDATA}/ytree"
check "Yleaf's tables and the current tree install" "${REPO}/scripts/setup.sh" --yleaf-data "$GENOME_DIR"
check "the tree folder has both position tables and the callable mask" \
  test -s "${YT}/positions_ftdna_hg38.txt" -a -s "${YT}/positions_yfull_hg38.txt" \
  -a -s "${YT}/callable_mask_chrY.fasta.gz"
MASK="${YT}/callable_mask_chrY.fasta.gz"

# A male copy of the sample: the BAM and the VCF hard-linked, a male indexcov
# row (indexcov reads the fixture's slices as female).
M="${SAMPLE}yp"
mkdir -p "${GENOME_DIR}/${M}/aligned" "${GENOME_DIR}/${M}/indexcov" "${GENOME_DIR}/${M}/vcf"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${M}/aligned/${M}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${M}/aligned/${M}_sorted.bam.bai"
ln -f "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${GENOME_DIR}/${M}/vcf/${M}.vcf.gz"
ln -f "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${GENOME_DIR}/${M}/vcf/${M}.vcf.gz.tbi"
printf '#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\n%s\t%s\t-9\t-9\t1\t-9\t1.0\t1.0\n' \
  "$M" "$M" > "${GENOME_DIR}/${M}/indexcov/indexcov-indexcov.ped"

run_step 37-y-haplogroup.sh "$M"
check_step_exit 37-y-haplogroup.sh
OUT="${GENOME_DIR}/${M}/y_haplogroup"
PRIV="${OUT}/${M}_y_private_snps.tsv"
head -4 "$PRIV" 2>/dev/null
check "the step wrote the private-SNP list" test -s "$PRIV"
check "the step line names the file" has "^Private SNPs: ${OUT}/${M}_y_private_snps\\.tsv\$" "$(cat "$STEP_LOG")"
check "its header records the gates and the sample's median chrY depth" \
  has "^# gates${T}min_depth=10 min_gq=20 min_fraction=0\.9 max_depth_factor=2\.0 median_chrY_depth=[0-9.]+ callable_mask=callable_mask_chrY\.fasta\.gz\$" \
  "$(cat "$PRIV")"
check "its header has the columns of a row" has "^Pos${T}Ref${T}Alt${T}Depth${T}Alt_fraction${T}GQ${T}Genotype\$" "$(cat "$PRIV")"

# dropped FIELD: the count of one gate from the header's dropped line
dropped() { sed -n 's/^# dropped\t//p' "$PRIV" | tr ' ' '\n' | sed -n "s/^${1}=//p"; }
check_ge "the subtraction dropped real calls the public tables name" "$(dropped on_public_lists)" 1

# Every listed position must be absent from both tables: the subtraction, on
# the real tables, checked from the file itself rather than from its counter.
awk -F'\t' 'FNR == NR { if ($0 !~ /^(#|Pos\t)/) want[$1] = 1; next }
  $4 in want { print $4; found = 1 } END { exit found }' \
  "$PRIV" "${YT}/positions_ftdna_hg38.txt" "${YT}/positions_yfull_hg38.txt" > "${CASE_TMP}/leaked.txt"
LEAKED=$(grep -c . "${CASE_TMP}/leaked.txt" 2>/dev/null || true)
check_eq "no listed position is on either public table" "${LEAKED:-unknown}" 0
ROWS=$(awk -F'\t' '$0 !~ /^(#|Pos\t)/' "$PRIV" 2>/dev/null | grep -c . || true)
echo "private chrY SNVs listed: ${ROWS}"

# --- Red controls on a seeded copy of the VCF -----------------------------------------------
# Two synthetic haploid chrY calls at positions no table names and the VCF
# does not have: IN where the mask passes the base, OUT where it does not.
# Depth is the sample's own median (at least 10), so the depth gates keep both.
SEED="${CASE_TMP}/seeded.vcf"
python3 - "${GENOME_DIR}/${M}/vcf/${M}.vcf.gz" "$MASK" "$PRIV" "$SEED" \
  "${YT}/positions_ftdna_hg38.txt" "${YT}/positions_yfull_hg38.txt" > "${CASE_TMP}/seed.txt" <<'PY'
import gzip, re, sys
vcf, mask_path, priv, out, *tables = sys.argv[1:]
mask = "".join(l.strip() for l in gzip.open(mask_path, "rt") if not l.startswith(">"))
taken = set()
for t in tables:
    for line in open(t):
        c = line.split("\t")
        if len(c) >= 4 and c[3].isdigit():
            taken.add(int(c[3]))
head, body = [], []
for line in gzip.open(vcf, "rt"):
    (head if line.startswith("#") else body).append(line)
    c = line.split("\t")
    if not line.startswith("#") and c[0] in ("chrY", "Y") and c[1].isdigit():
        taken.add(int(c[1]))
med = re.search(r"median_chrY_depth=([0-9.]+)", open(priv).read())
depth = max(10, int(float(med.group(1)))) if med else 20
# The fixture's TBL1Y stretch (scripts/ci/build-fixture.sh).
free = [p for p in range(6800000, 7000000) if p not in taken]
pin = next(p for p in free if mask[p - 1] == "P")
pout = next(p for p in free if mask[p - 1] != "P")
for p in (pin, pout):
    body.append(f"chrY\t{p}\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP:GQ\t1:0,{depth}:{depth}:60\n")
open(out, "w").writelines(head + body)
print(pin, pout)
PY
read -r P_IN P_OUT < "${CASE_TMP}/seed.txt"
echo "seeded: ${P_IN:-none} (callable), ${P_OUT:-none} (outside the mask)"
python3 "${REPO}/bin/y_private_snps.py" --vcf "$SEED" --tree-dir "$YT" --sample "$M" \
  --out-prefix "${CASE_TMP}/seeded" > "${CASE_TMP}/seeded.log" 2>&1
SEEDED="${CASE_TMP}/seeded_y_private_snps.tsv"
check "the seeded call inside the callable mask is listed" has "^${P_IN:-none}${T}" "$(cat "$SEEDED" 2>/dev/null)"
check "red control: the seeded call outside the callable mask is not listed" \
  lacks "^${P_OUT:-none}${T}" "$(cat "$SEEDED" 2>/dev/null)"
OUTSIDE_BASE=$(dropped outside_callable)
OUTSIDE_SEED=$(sed -n 's/^# dropped\t//p' "$SEEDED" 2>/dev/null | tr ' ' '\n' | sed -n 's/^outside_callable=//p')
check_eq "red control: and it is counted as outside_callable" "${OUTSIDE_SEED:-unknown}" "$(( ${OUTSIDE_BASE:-0} + 1 ))"

CT="${CASE_TMP}/tables"
rm -rf "$CT"; mkdir -p "$CT"
cp "${YT}/positions_yfull_hg38.txt" "${CT}/positions_yfull_hg38.txt"
{ cat "${YT}/positions_ftdna_hg38.txt"; printf 'chrY\tE2E_CONTROL\tR\t%s\tC->T\tC\tT\n' "${P_IN:-0}"; } \
  > "${CT}/positions_ftdna_hg38.txt"
python3 "${REPO}/bin/y_private_snps.py" --vcf "$SEED" --tree-dir "$CT" --callable-mask "$MASK" \
  --sample "$M" --out-prefix "${CASE_TMP}/control" > "${CASE_TMP}/control.log" 2>&1
cat "${CASE_TMP}/control.log"
check "red control: the seeded position added to a table leaves the list" \
  lacks "^${P_IN:-none}${T}" "$(cat "${CASE_TMP}/control_y_private_snps.tsv" 2>/dev/null)"
check "red control: and the rest of the seeded list is unchanged" \
  test "$(awk -F'\t' -v p="${P_IN:-0}" '$0 !~ /^(#|Pos\t)/ && $1 != p' "$SEEDED" | md5sum | cut -d' ' -f1)" = \
       "$(awk -F'\t' '$0 !~ /^(#|Pos\t)/' "${CASE_TMP}/control_y_private_snps.tsv" 2>/dev/null | md5sum | cut -d' ' -f1)"
: > "${CASE_TMP}/empty_ftdna.txt"
mkdir -p "${CASE_TMP}/notables"
cp "${CASE_TMP}/empty_ftdna.txt" "${CASE_TMP}/notables/positions_ftdna_hg38.txt"
cp "${CASE_TMP}/empty_ftdna.txt" "${CASE_TMP}/notables/positions_yfull_hg38.txt"
if python3 "${REPO}/bin/y_private_snps.py" --vcf "${GENOME_DIR}/${M}/vcf/${M}.vcf.gz" \
     --tree-dir "${CASE_TMP}/notables" --sample "$M" --out-prefix "${CASE_TMP}/none" \
     > "${CASE_TMP}/none.log" 2>&1; then
  fail "red control: tables with no position should stop the tool"
else
  check "red control: tables with no position stop the tool" \
    has 'name no position' "$(cat "${CASE_TMP}/none.log")"
fi

# --- The current-tree placement in the bash report --------------------------------------------
TREE_HG=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
  NR == 2 { print $c["Hg"] }' "${OUT}/${M}_y_tree.txt" 2>/dev/null)
check "step 37 wrote the current-tree placement" test -n "${TREE_HG:-}" -a "${TREE_HG:-NA}" != NA
printf 'run\tscript\te2e\n' > "${GENOME_DIR}/${M}/run_manifest.tsv"
run_step 24-html-report.sh "$M"
check_step_exit 24-html-report.sh
HTML=$(cat "${GENOME_DIR}/${M}/${M}_report.html" 2>/dev/null)
check "the HTML report has the current-tree card with the branch" \
  has "Y Placement on the Current Public Tree" "$HTML"
check "and the card names the branch step 37 placed" has "${TREE_HG}" "$HTML"
check "Yleaf's own card is still there, so the two results stand side by side" \
  has 'Y-Chromosome Haplogroup' "$HTML"
# section DIR KEY: the state of a summary.json section
section() { python3 -c "import json, sys; print(json.load(open(sys.argv[1]))['sections'].get(sys.argv[2], {}).get('state', 'absent'))" \
  "$1/summary.json" "$2" 2>/dev/null; }
check_eq "the bash summary.json section y_tree is filled" "$(section "${GENOME_DIR}/${M}" y_tree)" ok
run_step generate-report.sh "$M"
check_step_exit generate-report.sh
TXT=$(cat "${GENOME_DIR}/${M}/${M}_report.txt" 2>/dev/null)
sed -n '/Y Placement on the Current Public Tree/,/^$/p' <<< "$TXT"
check "the text report shows the branch and its level" \
  has "^  Branch: ${TREE_HG} \\(level [0-9]+, " "$TXT"
check "and the path counts" has '^  Path SNPs: [0-9]+ derived, [0-9]+ ancestral, [0-9]+ no call' "$TXT"

# --- The same section in the Nextflow report (case nextflow-hardening-2-nf) -------------------
R="${GENOME_DIR}/nf-hardening/${SAMPLE}"
NF_TREE="${R}/y_haplogroup/${SAMPLE}_y_tree.txt"
if [ -s "$NF_TREE" ] && [ -s "${R}/${SAMPLE}_report.html" ]; then
  NF_HG=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } NR == 2 { print $c["Hg"] }' "$NF_TREE")
  NF_PRIV="${R}/y_haplogroup/${SAMPLE}_y_private_snps.tsv"
  check "the Nextflow run published the private-SNP list" test -s "$NF_PRIV"
  python3 "${REPO}/bin/y_private_snps.py" --vcf "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" --tree-dir "$YT" \
    --sample "$SAMPLE" --out-prefix "${CASE_TMP}/nf-same" > /dev/null 2>&1
  check "and its rows are the tool's own on the same VCF and tables" \
    test "$(awk -F'\t' '$0 !~ /^(#|Pos\t)/' "$NF_PRIV" 2>/dev/null | md5sum | cut -d' ' -f1)" = \
         "$(awk -F'\t' '$0 !~ /^(#|Pos\t)/' "${CASE_TMP}/nf-same_y_private_snps.tsv" 2>/dev/null | md5sum | cut -d' ' -f1)"
  check "and that list has its header, so an empty file is not taken for a match" \
    has "^Pos${T}Ref${T}Alt" "$(cat "$NF_PRIV" 2>/dev/null)"
  check_eq "the Nextflow summary.json section y_tree is filled" "$(section "$R" y_tree)" ok
  check "the Nextflow HTML report has the current-tree card" \
    has "Y Placement on the Current Public Tree" "$(cat "${R}/${SAMPLE}_report.html")"
  check "and it names the branch the Nextflow run placed" \
    has "${NF_HG}" "$(cat "${R}/${SAMPLE}_report.html")"
  # Red control: the report code over a copy of those results without the
  # tree file. It shows the renderer leaves a missing section out; the
  # positive checks above are what a missing main.nf join would turn red.
  CC="${CASE_TMP}/nf-no-tree"
  rm -rf "$CC"; mkdir -p "${CC}/${SAMPLE}"
  cp -R "${R}/y_haplogroup" "${CC}/${SAMPLE}/y_haplogroup"
  rm -f "${CC}/${SAMPLE}/y_haplogroup/${SAMPLE}_y_tree.txt" "${CC}/${SAMPLE}/y_haplogroup/${SAMPLE}_y_tree_path.tsv"
  python3 "${REPO}/bin/render_report.py" --sample "$SAMPLE" --sample-dir "${CC}/${SAMPLE}" \
    --assets "${REPO}/assets" --json "${CC}/${SAMPLE}/summary.json" -o "${CC}/${SAMPLE}/report.html" \
    > "${CASE_TMP}/nf-control.log" 2>&1 || true
  tail -3 "${CASE_TMP}/nf-control.log"
  check_eq "red control: without the tree file the section is missing" "$(section "${CC}/${SAMPLE}" y_tree)" missing
  check "red control: and the card is out of the report" \
    lacks "Y Placement on the Current Public Tree" "$(cat "${CC}/${SAMPLE}/report.html" 2>/dev/null)"
else
  fail "case nextflow-hardening-2-nf left no Nextflow report with a tree file (${NF_TREE})"
fi

finish
