#!/usr/bin/env bash
# The Nextflow side of step 06's cDNA read check: CLINVAR_CDNA_CHECK runs
# bin/clip_check.py on CLINVAR_SCREEN's hits with the sample's BAM and
# --gencode_exons, as scripts/06 does. Case clinvar-clip-check-1 ran first and
# left the copy <sample>cdna, whose reads at the planted hit are spliced alt
# reads, and step 06's table for it. Its planted BAM holds one region only, so
# here it is merged with the real BAM's other reads: the BAM's sex check
# (INDEXCOV) then sees a whole sample. The exon file is built as that case
# builds it. One run, two rows:
#   1. <sample>cdnanf: the VCF with the merged BAM. The planted hit gets the
#      verdict and the read counts step 06 gave it (flagged, 12 of 12 alt
#      reads, 0 of 10 ref reads clipped), and no other hit is flagged.
#   2. <sample>nobam: the VCF alone. The module still writes the table, every
#      hit not checked with the reason "no BAM", instead of passing silently.
# Then a control: the same run with -resume and no --gencode_exons. The planted
# hit is then not checked, so the "flagged" check above can go red.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

C="${SAMPLE}cdna"
BASH_T="${GENOME_DIR}/${C}/clinvar/${C}_clinvar_cdna_check.tsv"
PLANTED_BAM="${GENOME_DIR}/${C}/aligned/${C}_sorted.bam"
HITS_TSV="${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.tsv"
check "case clinvar-clip-check-1 left its planted BAM" test -s "$PLANTED_BAM"
check "case clinvar-clip-check-1 left step 06's table for it" test -s "$BASH_T"
check "case 30 left a hits table" test -s "$HITS_TSV"
CHR=$(planted chrom)
POS=$(planted pos)
P0=$((POS - 1))
FROM=$((P0 - 2000))
TO=$((P0 + 5000))

# The exon file of case clinvar-clip-check-1.
EXONS="${CASE_TMP}/exons.bed"
{
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $1, $2 - 31, $2 + 29}' "$HITS_TSV"
  printf '%s\t%s\t%s\n' "$CHR" "$((P0 + 3000))" "$((P0 + 3150))"
} | sort -k1,1 -k2,2n > "$EXONS"

# The planted region plus the real BAM's reads outside it.
N="${SAMPLE}cdnanf"
D="${GENOME_DIR}/${N}/aligned"
rm -rf "${GENOME_DIR:?}/${N}"
mkdir -p "$D"
printf '%s\t%s\t%s\n' "$CHR" "$FROM" "$TO" > "${GENOME_DIR}/${N}/region.bed"
in_user() { docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$SAMTOOLS_IMAGE" samtools "$@"; }
in_user view -b -L "${N}/region.bed" -o /dev/null -U "${N}/aligned/outside.bam" "${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
in_user merge -f -c -p -o "${N}/aligned/${N}_sorted.bam" "${N}/aligned/outside.bam" "${C}/aligned/${C}_sorted.bam"
in_user index "${N}/aligned/${N}_sorted.bam"
rm -f "${D}/outside.bam"
check_ge "the merged BAM has reads" "$(in_user view -c "${N}/aligned/${N}_sorted.bam" 2>/dev/null || echo 0)" 1000

VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
SHEET="${CASE_TMP}/sheet.csv"
{
  printf 'sample,vcf,vcf_index,bam,bam_index\n'
  printf '%s,%s,%s,%s,%s\n' "$N" "$VCF" "${VCF}.tbi" "${D}/${N}_sorted.bam" "${D}/${N}_sorted.bam.bai"
  printf '%s,%s,%s,,\n' "${SAMPLE}nobam" "$VCF" "${VCF}.tbi"
} > "$SHEET"
NF_OUT="${CASE_TMP}/nf"
echo "+ nextflow run main.nf --tools clinvar --gencode_exons"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "${E2E_WORK}/nf-work-cdna" --input "$SHEET" \
    --reference "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" --tools clinvar \
    --clinvar "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
    --clinvar_index "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
    --gencode_exons "$EXONS" \
    --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB
) > "$STEP_LOG" 2>&1
NF_RC=$?
sed 's/^/    | /' "$STEP_LOG" | grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists'
check_eq "nextflow run exits 0" "$NF_RC" 0

# --- 1. the merged BAM: the verdict step 06 gave -------------------------------------
T1="${NF_OUT}/${N}/clinvar/${N}_clinvar_cdna_check.tsv"
row() { awk -F'\t' -v c="$CHR" -v p="$POS" '$1 == c && $2 == p' "$1" 2>/dev/null; }
NF_ROW=$(row "$T1")
BASH_ROW=$(row "$BASH_T")
echo "Nextflow planted row: ${NF_ROW}"
echo "step 06 planted row:  ${BASH_ROW}"
check_eq "CLINVAR_CDNA_CHECK writes the cdna_check column" \
  "$(head -n 1 "$T1" 2>/dev/null | tr '\t' '\n' | grep -cx cdna_check)" 1
check_eq "one row per hit" "$(awk 'NR > 1' "$T1" 2>/dev/null | wc -l | tr -d ' ')" \
  "$(awk 'NR > 1' "$HITS_TSV" | wc -l | tr -d ' ')"
check_eq "the planted hit is flagged" "$(cut -f9 <<< "$NF_ROW")" flagged
check_eq "the planted hit's counts and verdict match step 06's" "$(cut -f5-10 <<< "$NF_ROW")" "$(cut -f5-10 <<< "$BASH_ROW")"
check_eq "12 of 12 alt reads and 0 of 10 ref reads clipped at the exon edge" "$(cut -f5-8 <<< "$NF_ROW" | tr '\t' ' ')" "12 12 10 0"
check_eq "no other hit is flagged" \
  "$(awk -F'\t' -v c="$CHR" -v p="$POS" 'NR > 1 && $9 == "flagged" && !($1 == c && $2 == p)' "$T1" 2>/dev/null | wc -l | tr -d ' ')" 0
check "the task log names the flagged hit" grep -rqs "FLAGGED ${CHR}:${POS} " "${E2E_WORK}/nf-work-cdna" --include=.command.out

# --- 2. no BAM: every hit not checked, with the reason ---------------------------------
T2="${NF_OUT}/${SAMPLE}nobam/clinvar/${SAMPLE}nobam_clinvar_cdna_check.tsv"
check "CLINVAR_CDNA_CHECK writes a table for the sample without a BAM" test -s "$T2"
check_eq "every hit of the sample without a BAM is not checked" \
  "$(awk -F'\t' 'NR > 1 && $9 != "not_checked"' "$T2" 2>/dev/null | wc -l | tr -d ' ')" 0
check_eq "rows that give no BAM as the reason" \
  "$(awk -F'\t' 'NR > 1 && $10 == "no BAM for this sample in the samplesheet"' "$T2" 2>/dev/null | wc -l | tr -d ' ')" \
  "$(awk 'NR > 1' "$HITS_TSV" | wc -l | tr -d ' ')"

# --- 3. control: the same run without --gencode_exons -------------------------------
# Proves the checks of section 1 can go red: with no exon file the planted hit is
# not checked, so "the planted hit is flagged" would fail. -resume reuses every
# task but CLINVAR_CDNA_CHECK.
NF_OFF="${CASE_TMP}/nf-no-exons"
echo "+ nextflow run main.nf --tools clinvar -resume (no --gencode_exons)"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -resume \
    -work-dir "${E2E_WORK}/nf-work-cdna" --input "$SHEET" \
    --reference "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" --tools clinvar \
    --clinvar "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
    --clinvar_index "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
    --outdir "$NF_OFF" --max_cpus 4 --max_memory 14.GB
) > "${CASE_TMP}/no-exons.log" 2>&1
check_eq "control: nextflow run without --gencode_exons exits 0" "$?" 0
OFF_ROW=$(row "${NF_OFF}/${N}/clinvar/${N}_clinvar_cdna_check.tsv")
echo "control planted row:  ${OFF_ROW}"
check_eq "control: without --gencode_exons the planted hit is not checked, so the flagged check can fail" \
  "$(cut -f9 <<< "$OFF_ROW")" not_checked
check_eq "control: the reason names the missing exon file" "$(cut -f10 <<< "$OFF_ROW")" "no GENCODE exons (--gencode_exons)"

rm -rf "${GENOME_DIR:?}/${N}"
finish
