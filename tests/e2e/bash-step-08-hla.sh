#!/usr/bin/env bash
# Step 08 (T1K) on a pinned IPD-IMGT/HLA release, installed the way setup.sh
# installs it (install_data_file in scripts/lib/common.sh). The output names
# the release it was typed against, a typed gene gets an allele, and another
# release builds another index.
. "$(dirname "$0")/lib.sh"

# install RELEASE: hla.dat of that release and the GENCODE gene lines.
install() {
  ( HLA_DB_RELEASE=$1
    # shellcheck source=../../scripts/lib/common.sh
    . "${REPO}/scripts/lib/common.sh"
    install_data_file hla_dat && install_data_file gencode_genes )
}
check "IPD-IMGT/HLA 3.65.0 and the GENCODE genes install" install 3.65.0

run_step 08-hla-typing.sh "$SAMPLE"
check_step_exit 08-hla-typing.sh
REL="${GENOME_DIR}/${SAMPLE}/hla_t1k/database_release.txt"
cat "$REL" 2>/dev/null
check "database_release.txt names IPD-IMGT/HLA 3.65.0" grep -q '^database: IPD-IMGT/HLA Release Version 3\.65\.0$' "$REL"
GT="${GENOME_DIR}/${SAMPLE}/hla_t1k/${SAMPLE}_hla_genotype.tsv"
head -n 12 "$GT" 2>/dev/null
for g in HLA-A HLA-B HLA-C; do
  check_ge "${g} rows with a called allele" "$(awk -F'\t' -v g="$g" '$1 == g && $3 ~ /^HLA-/' "$GT" 2>/dev/null | wc -l | tr -d ' ')" 1
done
check_ge "index directories for release 3.65.0" \
  "$(find "${GENOME_DIR}/t1k_idx" -maxdepth 1 -type d -name 't1k-*_imgt-3.65.0_gencode-*' | wc -l | tr -d ' ')" 1

# Another release: a new index, and the output says which release it used.
check "IPD-IMGT/HLA 3.64.0 installs" install 3.64.0
HLA_DB_RELEASE=3.64.0 run_step 08-hla-typing.sh "$SAMPLE"
check_step_exit 08-hla-typing.sh
check "the log says it built an index for 3.64.0" has 'Building the T1K .* index for IPD-IMGT/HLA 3\.64\.0' "$(cat "$STEP_LOG")"
check_ge "index directories for release 3.64.0" \
  "$(find "${GENOME_DIR}/t1k_idx" -maxdepth 1 -type d -name 't1k-*_imgt-3.64.0_gencode-*' | wc -l | tr -d ' ')" 1
check "database_release.txt names IPD-IMGT/HLA 3.64.0" grep -q '^database: IPD-IMGT/HLA Release Version 3\.64\.0$' "$REL"

finish
