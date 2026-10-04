#!/usr/bin/env bash
# The Nextflow pipeline with real containers on the bash leg's VCF and BAM,
# through today's VCF+BAM samplesheet, with the tools that need no large database.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

NF_OUT="${GENOME_DIR}/nf-results"
SHEET="${E2E_WORK}/nf-samplesheet.csv"
G="$GENOME_DIR"
printf 'sample,vcf,vcf_index,bam,bam_index\n%s,%s,%s,%s,%s\n' "$SAMPLE" \
  "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" \
  "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" \
  > "$SHEET"

# Keep going after a failed task, so one run shows every broken module; the
# failed tasks are counted from the trace below instead.
cat > "${CASE_TMP}/e2e.config" <<'NFCONF'
process.errorStrategy = 'ignore'
NFCONF

echo "+ nextflow run main.nf -profile docker"
(
  cd "$E2E_WORK" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -c "${CASE_TMP}/e2e.config" \
    -work-dir "${E2E_WORK}/nf-work" \
    --input "$SHEET" \
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools clinvar,mosdepth,delly,vcfanno,roh,pharmcat,cpic \
    --clinvar "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
    --clinvar_index "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
    --revel "${G}/annotations/revel_grch38.tsv.gz" \
    --revel_index "${G}/annotations/revel_grch38.tsv.gz.tbi" \
    --outdir "$NF_OUT" \
    --max_cpus 4 --max_memory 14.GB
) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${E2E_WORK}/.nextflow.log" "${E2E_WORK}/logs/60-nextflow.nextflow.log" 2>/dev/null
TRACE=$(find "${NF_OUT}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
FAILED=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
                      $c["status"] != "COMPLETED" && $c["status"] != "CACHED" {print $c["name"] " (" $c["status"] ", exit " $c["exit"] ")"}' \
  "${TRACE:-/dev/null}" 2>/dev/null)
echo "tasks that did not complete: ${FAILED:-none}"
check "the trace lists the tasks" test -s "${TRACE:-/dev/null}"
check_eq "tasks that did not complete" "$(grep -c . <<< "$FAILED" || true)" 0
if [ -n "$FAILED" ]; then
  grep -E 'Error executing process|Command exit status|Command error' -A 6 "${E2E_WORK}/.nextflow.log" \
    | grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' | head -80
fi

R="nf-results/${SAMPLE}"
GENE=$(planted gene)
check_ge "CLINVAR_SCREEN reports the planted hit" "$(sample_side_hits "${GENOME_DIR}/${R}/clinvar" "")" 1
check_ge "CLINVAR_SCREEN names its gene (${GENE})" "$(sample_side_hits "${GENOME_DIR}/${R}/clinvar" "$GENE")" 1
check "MOSDEPTH summary has a chr20 row" grep -q '^chr20' "${GENOME_DIR}/${R}/coverage/${SAMPLE}.mosdepth.summary.txt"
check "DELLY writes a readable VCF" vcf_ok "${R}/delly/${SAMPLE}_sv.vcf.gz"
check "DELLY VCF is indexed" nonempty "${R}/delly/${SAMPLE}_sv.vcf.gz.tbi"
check "VCFANNO output is indexed" nonempty "${R}/vep/${SAMPLE}_annotated.vcf.gz.tbi"
check_ge "VCFANNO records carrying REVEL" "$(vcf_count -i 'INFO/REVEL!="."' "${R}/vep/${SAMPLE}_annotated.vcf.gz")" 1
check_ge "ROH output lines" "$(cat "${GENOME_DIR}/${R}"/roh/*_roh.txt 2>/dev/null | grep -vc '^#' || true)" 1
JSON=$(find "${GENOME_DIR}/${R}/pharmcat" -name '*.report.json' 2>/dev/null | awk 'NR == 1')
check_ge "PHARMCAT genes called" "$( [ -n "$JSON" ] && pharmcat_called "$JSON" | grep -c . || echo 0)" 1
check_ge "CPIC_LOOKUP gene rows" \
  "$(awk -F'\t' 'NF >= 3' "${GENOME_DIR}/${R}/cpic/${SAMPLE}_phenotypes.tsv" 2>/dev/null | wc -l | tr -d ' ')" 1

finish
