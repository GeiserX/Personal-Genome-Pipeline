#!/usr/bin/env bash
# Step 19 (Delly) with the exclude map setup.sh installs.
. "$(dirname "$0")/lib.sh"

check "Delly's exclude map installs" bash -c '. "$1/scripts/lib/common.sh" && install_data_file delly_exclude' _ "$REPO"
MAP="${GENOME_DIR}/reference/delly_human.hg38.excl.tsv"
check_ge "exclude map lines" "$(grep -c . "$MAP" 2>/dev/null || true)" 1000
# A whole contig is a line with its name alone.
check "the map excludes chr22_KI270879v1_alt" grep -qx 'chr22_KI270879v1_alt' "$MAP"

run_step 19-delly.sh "$SAMPLE"
check_step_exit 19-delly.sh
check "the log names the exclude map" has "Exclude map: ${MAP}" "$(cat "$STEP_LOG")"
VCF="${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz"
check "Delly's VCF is readable" vcf_ok "$VCF"
check "Delly's VCF is indexed" nonempty "${VCF}.tbi"
# The map's chr20 intervals as a BED under GENOME_DIR (bcftools reads it in
# its container): no call may start inside them. Lines that name a whole
# contig have no coordinates and are left out.
EXCL="${SAMPLE}/delly/chr20_exclude.bed"
awk 'BEGIN {OFS = "\t"} $1 == "chr20" && NF >= 3 {print $1, $2, $3}' "$MAP" > "${GENOME_DIR}/${EXCL}"
check_ge "exclude intervals on chr20" "$(grep -c . "${GENOME_DIR}/${EXCL}" || true)" 1
echo "Delly calls: $(vcf_count "$VCF") in all, $(vcf_count -r chr20 "$VCF") on chr20"
check_eq "calls inside the chr20 exclude intervals" "$(vcf_count -R "$EXCL" "$VCF")" 0

finish
