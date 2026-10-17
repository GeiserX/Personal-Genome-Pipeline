#!/usr/bin/env bash
# Step 08 (T1K) on a pinned IPD-IMGT/HLA release, installed the way setup.sh
# installs it (install_data_file in scripts/lib/common.sh). The output names
# the release it was typed against, a typed gene gets an allele, and another
# release builds another index. The default release is read from common.sh,
# never written here, so a pin bump runs this case on the new release; the
# calls of the two releases are printed side by side for the bump's review.
. "$(dirname "$0")/lib.sh"

# install RELEASE: hla.dat of that release and the GENCODE gene lines.
install() {
  ( HLA_DB_RELEASE=$1
    # shellcheck source=../../scripts/lib/common.sh
    . "${REPO}/scripts/lib/common.sh"
    install_data_file hla_dat && install_data_file gencode_genes )
}
# The pin: HLA_DB_RELEASE's default in common.sh.
PIN=$( unset HLA_DB_RELEASE
       # shellcheck source=../../scripts/lib/common.sh
       . "${REPO}/scripts/lib/common.sh"
       printf '%s' "$HLA_DB_RELEASE" )
PIN_RE=${PIN//./\\.}
# The release before the pin, typed second: a new index, and the calls to
# compare with. Change it only when the pin itself becomes 3.65.0.
OTHER=3.65.0
OTHER_RE=${OTHER//./\\.}
check "the pin (${PIN}) is a release and differs from ${OTHER}" \
  test -n "$PIN" -a "$PIN" != "$OTHER"

check "IPD-IMGT/HLA ${PIN} and the GENCODE genes install" install "$PIN"

# Explicit, so a HLA_DB_RELEASE in the environment cannot pick another release.
HLA_DB_RELEASE=$PIN run_step 08-hla-typing.sh "$SAMPLE"
check_step_exit 08-hla-typing.sh
REL="${GENOME_DIR}/${SAMPLE}/hla_t1k/database_release.txt"
cat "$REL" 2>/dev/null
check "database_release.txt names IPD-IMGT/HLA ${PIN}" grep -q "^database: IPD-IMGT/HLA Release Version ${PIN_RE}\$" "$REL"
GT="${GENOME_DIR}/${SAMPLE}/hla_t1k/${SAMPLE}_hla_genotype.tsv"
head -n 12 "$GT" 2>/dev/null
for g in HLA-A HLA-B HLA-C; do
  check_ge "${g} rows with a called allele" "$(awk -F'\t' -v g="$g" '$1 == g && $3 ~ /^HLA-/' "$GT" 2>/dev/null | wc -l | tr -d ' ')" 1
done
check_ge "index directories for release ${PIN}" \
  "$(find "${GENOME_DIR}/t1k_idx" -maxdepth 1 -type d -name "t1k-*_imgt-${PIN}_gencode-*" | wc -l | tr -d ' ')" 1
cp "$GT" "${CASE_TMP}/genotype_${PIN}.tsv" 2>/dev/null

# Another release: a new index, and the output says which release it used.
check "IPD-IMGT/HLA ${OTHER} installs" install "$OTHER"
HLA_DB_RELEASE=$OTHER run_step 08-hla-typing.sh "$SAMPLE"
check_step_exit 08-hla-typing.sh
check "the log says it built an index for ${OTHER}" has "Building the T1K .* index for IPD-IMGT/HLA ${OTHER_RE}" "$(cat "$STEP_LOG")"
check_ge "index directories for release ${OTHER}" \
  "$(find "${GENOME_DIR}/t1k_idx" -maxdepth 1 -type d -name "t1k-*_imgt-${OTHER}_gencode-*" | wc -l | tr -d ' ')" 1
check "database_release.txt names IPD-IMGT/HLA ${OTHER}" grep -q "^database: IPD-IMGT/HLA Release Version ${OTHER_RE}\$" "$REL"

# The calls of the two releases (gene and the two alleles, T1K columns 1, 3
# and 6), for the reviewer of a pin bump: a difference is printed, not
# failed, since a new release may rightly change a call.
calls() { awk -F'\t' '$1 ~ /^HLA-/ {print $1 "\t" $3 "\t" $6}' "$1" 2>/dev/null | sort; }
echo "=== HLA calls: IPD-IMGT/HLA ${OTHER} (left) and ${PIN} (right) ==="
if diff <(calls "$GT") <(calls "${CASE_TMP}/genotype_${PIN}.tsv"); then
  echo "identical calls on the fixture"
fi

finish
