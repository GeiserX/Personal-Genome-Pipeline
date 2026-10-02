#!/usr/bin/env bash
# Step 30 on the fixture's VEP output (written as step 13 writes it): the input
# gets compressed and indexed, and the score track lands in INFO.
. "$(dirname "$0")/lib.sh"

IN_N=$(grep -vc '^#' "${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf" || true)
run_step 30-vcfanno.sh "$SAMPLE"
check_step_exit 30-vcfanno.sh

GZ="${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
check "${SAMPLE}_vep.vcf.gz is a readable VCF" vcf_ok "$GZ"
check_eq "${SAMPLE}_vep.vcf.gz keeps every VEP record" "$(vcf_count "$GZ")" "$IN_N"
check "${SAMPLE}_vep.vcf.gz is indexed" nonempty "${GZ}.tbi"
OUTV="${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz"
check "annotated VCF is indexed" nonempty "${OUTV}.tbi"
check_ge "annotated records carrying REVEL" "$(vcf_count -i 'INFO/REVEL!="."' "$OUTV")" 1

finish
