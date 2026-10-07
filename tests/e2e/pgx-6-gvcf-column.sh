#!/usr/bin/env bash
# A samplesheet row with a BAM, a VCF and the gVCF step 03 wrote beside it
# (the bash leg of the parity check, sample HG002P): the pipeline calls
# nothing again (DEEPVARIANT runs no task), and PharmCAT and PRS read the
# given gVCF, so a PGx or score position where the sample matches the
# reference is a 0/0 call instead of missing.
#
# --sex_check warn: on the fixture's slices indexcov reads HG002 (male) as
# female (see case nextflow-from-fastq-2-nextflow).
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

P=HG002P
G="$GENOME_DIR"
D="${G}/${P}"
for f in "aligned/${P}_sorted.bam" "aligned/${P}_sorted.bam.bai" "vcf/${P}.vcf.gz" "vcf/${P}.vcf.gz.tbi" \
         "vcf/${P}.g.vcf.gz" "vcf/${P}.g.vcf.gz.tbi"; do
  check "case nextflow-from-fastq-1-bash wrote ${f}" nonempty "${P}/${f}"
done
PGS="${E2E_WORK}/parity-pgs"
check "case nextflow-from-fastq-1-bash wrote the score files" test -s "${PGS}/PGS000018.txt.gz"

NF_OUT="${G}/nf-gvcf"
WORK="${E2E_WORK}/nf-work-gvcf"
SHEET="${CASE_TMP}/samplesheet.csv"
rm -rf "$NF_OUT" "$WORK"
printf 'sample,bam,bam_index,vcf,vcf_index,gvcf,gvcf_index,sex\n%s,%s,%s,%s,%s,%s,%s,male\n' "$P" \
  "${D}/aligned/${P}_sorted.bam" "${D}/aligned/${P}_sorted.bam.bai" "${D}/vcf/${P}.vcf.gz" "${D}/vcf/${P}.vcf.gz.tbi" \
  "${D}/vcf/${P}.g.vcf.gz" "${D}/vcf/${P}.g.vcf.gz.tbi" > "$SHEET"
cat "$SHEET"

echo "+ nextflow run main.nf -profile docker (BAM + VCF + gVCF row)"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "$WORK" \
    --input "$SHEET" \
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --sex_check warn \
    --tools pharmcat,prs \
    --pgs_scoring "$PGS" \
    --outdir "$NF_OUT" \
    --max_cpus 4 --max_memory 14.GB
) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
TRACE=$(find "${NF_OUT}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
check "the trace lists the tasks" test -s "${TRACE:-/dev/null}"

# tasks NAME: the tasks of process NAME in the trace, whatever their status
tasks() {
  awk -F'\t' -v p="$1" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    { n = $c["name"]; sub(/ \(.*/, "", n); sub(/.*:/, "", n); if (n == p) k++ }
    END {print k + 0}' "${TRACE:-/dev/null}" 2>/dev/null
}
check_eq "DEEPVARIANT tasks (the row is not called again)" "$(tasks DEEPVARIANT)" 0
check_eq "PHARMCAT_PREPROCESS tasks" "$(tasks PHARMCAT_PREPROCESS)" 1

# What PHARMCAT_PREPROCESS read: its task directory, found by the trace hash.
HASH=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["name"] ~ /PHARMCAT_PREPROCESS/ {print $c["hash"]; exit}' \
  "${TRACE:-/dev/null}" 2>/dev/null)
OUT=$(cat "${WORK}/${HASH:-none}"*/.command.out 2>/dev/null)
echo "PHARMCAT_PREPROCESS (${HASH:-no task}): $(head -n 1 <<< "$OUT")"
check "PharmCAT's input is the gVCF, its reference blocks expanded" has "^Input: ${P}\.g\.vcf\.gz \(reference blocks expanded" "$OUT"
JSON="${NF_OUT}/${P}/pharmcat/${P}.report.json"
check_ge "PHARMCAT genes called" "$( [ -s "$JSON" ] && pharmcat_called "$JSON" | grep -c . || echo 0)" 1
check_eq "PRS input (hom-ref sites from the gVCF)" \
  "$(awk -F'\t' 'NR == 2 {print $NF}' "${NF_OUT}/${P}/prs/${P}_prs_summary.tsv" 2>/dev/null)" gvcf

finish
