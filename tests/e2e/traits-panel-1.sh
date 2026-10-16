#!/usr/bin/env bash
# Step 42, the traits panel, on the HG002 gVCF case 21 wrote, and its
# Nextflow mirror:
#   - assets/traits.tsv: every reference base is the fixture reference's base
#     (on the contigs it has), every PubMed id resolves and every rsID's GRCh38
#     position and alleles are dbSNP's (tests/test_traits_table.py --fasta
#     --online); control: a copy with one wrong reference base and one PubMed
#     id that does not exist must fail, naming both;
#   - scripts/42-traits.sh: one row per table row plus the APOE line; the
#     table positions inside the fixture's slices are read (called) with the
#     genotype of the GIAB v4.2.1 truth, every position outside them is
#     not_covered (no reads, never taken for the reference);
#   - a planted alternate allele: a copy of the gVCF with CYP2A6*2 (rs1801272)
#     homozygous gives that row the table's hom_alt reading, not HG002's;
#   - both reports (bin/render_report.py on the sample folder) carry the
#     section, the HTML one with the APOE lines folded away;
#   - nextflow run main.nf --tools traits,html_report on the same VCF and
#     gVCF publishes the same two files as the bash step, byte for byte, and
#     its HTML report has the traits card.
# Writes traits/ in the sample folder, which no later case reads; removes the
# rest of what it made.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
TABLE="${REPO}/assets/traits.tsv"
GVCF="${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
check "case 21 wrote a gVCF" vcf_ok "$GVCF"

# --- the table: reference bases, PubMed ids, dbSNP --------------------------------
python3 "${REPO}/tests/test_traits_table.py" --fasta "$REF" --min-fasta-rows 20 --online > "${CASE_TMP}/table.log" 2>&1
TRC=$?
cat "${CASE_TMP}/table.log"
check_eq "tests/test_traits_table.py --fasta --online exits 0" "$TRC" 0
check "the reference bases were checked on the fixture contigs" has 'ref checked against .*: [0-9]+ rows' "$(cat "${CASE_TMP}/table.log")"

# Control: one wrong reference base (rs1801133, chr1, a fixture contig) and one
# PubMed id that does not exist must each fail.
BAD="${CASE_TMP}/traits-planted.tsv"
awk -F'\t' -v OFS='\t' '$3 == "rs1801133" { $5 = "C" } $3 == "rs671" { $13 = "99999999" } { print }' "$TABLE" > "$BAD"
check_eq "control: the planted copy differs from the table in two lines" \
  "$(diff "$TABLE" "$BAD" | grep -c '^>' || true)" 2
python3 "${REPO}/tests/test_traits_table.py" --table "$BAD" --fasta "$REF" --online > "${CASE_TMP}/table-bad.log" 2>&1
BRC=$?
BADLOG=$(cat "${CASE_TMP}/table-bad.log")
check "control: the planted copy fails (exit ${BRC})" test "$BRC" -ne 0
check "control: the wrong reference base is named" has 'rs1801133 chr1:11796321: table ref C, FASTA G' "$BADLOG"
check "control: the PubMed id that does not exist is named" has 'PMID 99999999 does not resolve' "$BADLOG"

# --- the bash step ------------------------------------------------------------------
run_step 42-traits.sh "$SAMPLE"
check_step_exit 42-traits.sh
OUT="${G}/${SAMPLE}/traits/${SAMPLE}_traits.tsv"
ROWS=$(grep -v '^#' "$TABLE" | awk 'NR > 1' | wc -l | tr -d ' ')
check_eq "one output row per table row plus the APOE line" "$(awk 'NR > 1' "$OUT" 2>/dev/null | wc -l | tr -d ' ')" "$((ROWS + 1))"

# Which table positions lie in the fixture's slices (regions.bed is 0-based).
awk -F'\t' 'NR == FNR { s[NR] = $1 "\t" $2 "\t" $3; n = NR; next }
  FNR > 1 && $1 !~ /^#/ && $1 != "chrom" { inside = 0
    for (i = 1; i <= n; i++) { split(s[i], r, "\t"); if (r[1] == $1 && $2 > r[2] && $2 <= r[3]) inside = 1 }
    print $3 "\t" $1 "\t" $2 "\t" $5 "\t" $6 "\t" (inside ? "in" : "out") }' \
  "${FIXTURE_DIR}/regions.bed" <(grep -v '^#' "$TABLE") > "${CASE_TMP}/where.tsv"
IN=$(awk -F'\t' '$6 == "in"' "${CASE_TMP}/where.tsv" | wc -l | tr -d ' ')
check_ge "table positions inside the fixture's slices" "$IN" 2
status_of() { awk -F'\t' -v r="$1" '$1 == r { print $7 }' "$OUT"; }
geno_of()   { awk -F'\t' -v r="$1" '$1 == r { print $6 }' "$OUT"; }
OUTSIDE_NC=0 OUTSIDE=0
while IFS=$'\t' read -r rs c p ref _alt where; do
  if [ "$where" = out ]; then
    OUTSIDE=$((OUTSIDE + 1)); [ "$(status_of "$rs")" != not_covered ] || OUTSIDE_NC=$((OUTSIDE_NC + 1))
    continue
  fi
  # The GIAB v4.2.1 truth: its record at the position, or ref/ref inside the benchmark regions.
  T=$(docker run --rm -v "${FIXTURE_DIR}:/fx:ro" "$BCFTOOLS_IMAGE" bcftools query -r "${c}:${p}" -f '%POS\t%REF\t[%TGT]\n' \
        /fx/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz 2>/dev/null | awk -F'\t' -v p="$p" '$1 == p && length($2) == 1 { print $3; exit }')
  INBED=$(awk -F'\t' -v c="$c" -v p="$p" '$1 == c && p > $2 && p <= $3 { print 1; exit }' \
        "${FIXTURE_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")
  if [ -n "$T" ]; then
    WANT=$(tr '|' '/' <<< "$T" | awk -F/ -v r="$ref" '{ if ($2 == r) print $2 "/" $1; else print $1 "/" $2 }')
  elif [ -n "$INBED" ]; then WANT="${ref}/${ref}"
  else WANT="outside the GIAB benchmark regions"; fi
  check_eq "${rs} (${c}:${p}) is read" "$(status_of "$rs")" called
  check_eq "${rs} genotype is the GIAB truth's" "$(geno_of "$rs")" "$WANT"
done < "${CASE_TMP}/where.tsv"
check_eq "table positions outside the slices that are not_covered" "$OUTSIDE_NC" "$OUTSIDE"
check_ge "control: positions outside the slices exist" "$OUTSIDE" 10
check_eq "the APOE line (its two sites are outside the slices) is not read" "$(status_of rs429358+rs7412)" not_covered
echo "Fixture genotypes at the table positions inside the slices:"
awk -F'\t' '$7 == "called" { print "  " $1 " " $6 " " $8 }' "$OUT"

# --- a planted alternate allele -----------------------------------------------------
P2="${SAMPLE}T"
rm -rf "${G:?}/${P2}"
mkdir -p "${G}/${P2}/vcf"
PL_POS=$(awk -F'\t' '$1 == "rs1801272" { print $3 }' "${CASE_TMP}/where.tsv")
check "rs1801272 lies in a fixture slice" test "$(awk -F'\t' '$1 == "rs1801272" { print $6 }' "${CASE_TMP}/where.tsv")" = in
# The gVCF without the records over the position, plus one hom-alt record there.
python3 - "${G}/${GVCF}" "${G}/${P2}/vcf/planted.vcf" "$PL_POS" "$SAMPLE" <<'PY'
import gzip, sys
src, dst, pos, sample = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
rec = f"chr19\t{pos}\t.\tA\tT,<*>\t40\tPASS\t.\tGT:GQ:DP:AD:VAF:PL\t1/1:40:30:0,30,0:1,0:990,90,0,990,990,990\n"
done = False
with gzip.open(src, "rt") as f, open(dst, "w") as o:
    for line in f:
        if line.startswith("#"):
            o.write(line)
            continue
        c = line.split("\t", 8)
        p = int(c[1])
        end = p + len(c[3]) - 1
        for kv in c[7].split(";"):
            if kv.startswith("END="):
                end = int(kv[4:])
        if c[0] == "chr19" and p <= pos <= end:
            continue
        if not done and c[0] == "chr19" and p > pos:
            o.write(rec)
            done = True
        o.write(line)
    if not done:
        sys.exit("no chr19 record after the planted position")
PY
in_genome "$BCFTOOLS_IMAGE" bash -c "bcftools view -Oz -o '${P2}/vcf/${P2}.g.vcf.gz' '${P2}/vcf/planted.vcf' && bcftools index -t '${P2}/vcf/${P2}.g.vcf.gz'"
check_eq "the planted gVCF has the hom-alt record" \
  "$(bcf query -r "chr19:${PL_POS}" -f '[%GT]\n' "${P2}/vcf/${P2}.g.vcf.gz" 2>/dev/null | head -n 1)" 1/1
run_step 42-traits.sh "$P2"
check_step_exit 42-traits.sh
OUT2="${G}/${P2}/traits/${P2}_traits.tsv"
HOMALT=$(awk -F'\t' '$3 == "rs1801272" { print $10 }' "$TABLE")
R1=$(awk -F'\t' '$1 == "rs1801272" { print $8 }' "$OUT")
R2=$(awk -F'\t' '$1 == "rs1801272" { print $8 }' "$OUT2" 2>/dev/null)
check_eq "planted: rs1801272 genotype" "$(awk -F'\t' '$1 == "rs1801272" { print $6 }' "$OUT2" 2>/dev/null)" T/T
check_eq "planted: rs1801272 reading is the table's hom_alt" "$R2" "$HOMALT"
check "planted: the reading differs from HG002's (${R1})" test "$R1" != "$R2"
check_eq "planted: the other rows read as in HG002" \
  "$(awk -F'\t' '$1 != "rs1801272"' "$OUT2" 2>/dev/null | md5sum)" "$(awk -F'\t' '$1 != "rs1801272"' "$OUT" | md5sum)"
cp "$OUT2" "${CASE_TMP}/planted_traits.tsv" 2>/dev/null
in_genome "$BCFTOOLS_IMAGE" rm -rf "$P2"

# --- both reports ---------------------------------------------------------------------
python3 "${REPO}/bin/render_report.py" --sample "$SAMPLE" --sample-dir "${G}/${SAMPLE}" --json "${CASE_TMP}/summary.json" \
  -o "${CASE_TMP}/report.txt" -o "${CASE_TMP}/report.html" > /dev/null
TXT=$(cat "${CASE_TMP}/report.txt" 2>/dev/null) HTML=$(cat "${CASE_TMP}/report.html" 2>/dev/null)
ONE=$(awk -F'\t' '$7 == "called" { print $1; exit }' "$OUT")
check "text report: the traits section" has '^## Traits Panel \(step 42\)$' "$TXT"
check "text report: a read row (${ONE}) with its genotype" has "${ONE}\\): $(geno_of "$ONE")" "$TXT"
check "HTML report: the traits card" has '<h2>Traits</h2>' "$HTML"
check "HTML report: the row of ${ONE}" has "${ONE}</span></td><td>$(geno_of "$ONE")" "$HTML"
check "HTML report: the APOE lines are folded away" has '<details><summary>APOE and Alzheimer' "$HTML"

# --- the Nextflow mirror ----------------------------------------------------------------
if command -v nextflow >/dev/null; then
  NF_OUT="${G}/nf-traits"
  SHEET="${CASE_TMP}/samplesheet.csv"
  printf 'sample,vcf,vcf_index,gvcf,gvcf_index,sex\n%s,%s,%s,%s,%s,male\n' "$SAMPLE" \
    "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${G}/${GVCF}" "${G}/${GVCF}.tbi" > "$SHEET"
  (
    cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-traits" \
      --input "$SHEET" --reference "$REF" --tools traits,html_report \
      --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB
  ) 2>&1 | tee "$STEP_LOG"
  NF_RC=${PIPESTATUS[0]}
  check_eq "nextflow run --tools traits,html_report exits 0" "$NF_RC" 0
  for f in "${SAMPLE}_traits.tsv" "${SAMPLE}_traits_genotypes.tsv"; do
    check "Nextflow published traits/${f}" test -s "${NF_OUT}/${SAMPLE}/traits/${f}"
    check "traits/${f}: Nextflow and the bash step agree byte for byte" \
      cmp -s "${NF_OUT}/${SAMPLE}/traits/${f}" "${G}/${SAMPLE}/traits/${f}"
  done
  check "control: the comparison sees the planted sample's table as different" \
    test "$(cmp -s "${NF_OUT}/${SAMPLE}/traits/${SAMPLE}_traits.tsv" "${CASE_TMP}/planted_traits.tsv"; echo $?)" -ne 0
  NFHTML=$(cat "${NF_OUT}/${SAMPLE}/${SAMPLE}_report.html" 2>/dev/null)
  check "Nextflow HTML report: the traits card with the row of ${ONE}" has "${ONE}</span></td><td>$(geno_of "$ONE")" "$NFHTML"
else
  fail "nextflow is not on PATH"
fi

finish
