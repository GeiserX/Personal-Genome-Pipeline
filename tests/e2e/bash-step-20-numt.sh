#!/usr/bin/env bash
# Step 20 runs NuMTFilterTool after FilterMutectCalls. The fixture covers a few
# slices of each autosome, so the median autosomal depth mosdepth reports is
# 0 and nothing is marked; AUTOSOMAL_COVERAGE set far above the chrM depth
# marks low-depth alleles as possible NuMTs, which shows the filter is live.
. "$(dirname "$0")/lib.sh"

VCF="${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz"
run_step 20-mtoolbox.sh "$SAMPLE"
check_step_exit 20-mtoolbox.sh
check "the log gives the median autosomal depth from mosdepth" has 'Median autosomal coverage: [0-9]+ \(from .*mosdepth\.global\.dist\.txt\)' "$(cat "$STEP_LOG")"
check "the filtered VCF declares possible_numt" has 'ID=possible_numt' "$(bcf view -h "$VCF" 2>/dev/null)"
check_ge "chrM records" "$(vcf_count "$VCF")" 5
BASE=$(vcf_count -i 'FILTER~"possible_numt"' "$VCF")

AUTOSOMAL_COVERAGE=100000 run_step 20-mtoolbox.sh "$SAMPLE"
check_step_exit 20-mtoolbox.sh
HIGH=$(vcf_count -i 'FILTER~"possible_numt"' "$VCF")
echo "possible_numt records: ${BASE} at the measured depth, ${HIGH} at 100000x"
check "more records are marked at 100000x than at the measured depth" test "${HIGH:-0}" -gt "${BASE:-0}"

# Leave the outputs of the measured depth for the cases after this one.
run_step 20-mtoolbox.sh "$SAMPLE"
check_step_exit 20-mtoolbox.sh

finish
