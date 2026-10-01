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

echo "+ nextflow run main.nf -profile docker"
(
  cd "$E2E_WORK" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "${E2E_WORK}/nf-work" \
    --input "$SHEET" \
    --reference "${G}/reference/Homo_sapiens_assembly38.fasta" \
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
if [ "$NF_RC" -ne 0 ] && [ -f "${E2E_WORK}/.nextflow.log" ]; then
  echo "--- failed tasks (from .nextflow.log) ---"
  grep -E 'Error executing process|Command exit status|Command error' -A 6 "${E2E_WORK}/.nextflow.log" | head -80
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
