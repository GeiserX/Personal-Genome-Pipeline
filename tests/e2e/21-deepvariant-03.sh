#!/usr/bin/env bash
# Step 03 calls small variants; the VCF names the sample and holds the planted SNV.
. "$(dirname "$0")/lib.sh"

# INTERVALS limits DeepVariant to the fixture's slices (the rest of the
# reference has no reads). It takes effect once step 03 honours INTERVALS the
# way step 03a does; until then DeepVariant walks all 1.8 Gb of reference,
# about 45 minutes on a 4-CPU runner.
INTERVALS=$(awk '{printf "%s%s:%d-%d", (NR > 1 ? " " : ""), $1, $2 + 1, $3}' "${FIXTURE_DIR}/regions.bed")
export INTERVALS
run_step 03-deepvariant.sh "$SAMPLE"
check_step_exit 03-deepvariant.sh

VCF="${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
check "VCF is readable" vcf_ok "$VCF"
check "VCF index exists" nonempty "${VCF}.tbi"
check_ge "PASS records on the chr20 slice" "$(vcf_count -f PASS -r chr20:10000000-10500000 "$VCF")" 200
check_eq "VCF sample column" "$(bcf query -l "$VCF" 2>/dev/null)" "$SAMPLE"
GT=$(bcf query -r "$(planted chrom):$(planted pos)" -f '[%GT]\n' "$VCF" 2>/dev/null | awk 'NR == 1')
check "planted SNV $(planted chrom):$(planted pos) is called non-reference (GT ${GT:-none})" has '1' "${GT:-}"

finish
