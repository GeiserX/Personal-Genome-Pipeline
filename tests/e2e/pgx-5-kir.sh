#!/usr/bin/env bash
# Step 08 with KIR=true (opt-in): a second T1K pass over the KIR genes
# against IPD-KIR release KIR_DB_RELEASE, which `setup.sh --kir-data`
# installs. The fixture has no slice of the KIR region (chr19 leukocyte
# receptor complex), so the step may find too few reads: then the genotype
# file says so. Either way the release is written beside it.
. "$(dirname "$0")/lib.sh"

echo "+ scripts/setup.sh --kir-data"
"${REPO}/scripts/setup.sh" --kir-data "$GENOME_DIR"
check "setup.sh --kir-data installed IPD-KIR ${KIR_DB_RELEASE}" test -s "${GENOME_DIR}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat"

KIR=true run_step 08-hla-typing.sh "$SAMPLE"
check_step_exit 08-hla-typing.sh
K="${GENOME_DIR}/${SAMPLE}/kir_t1k"
cat "${K}/database_release.txt" "${K}/${SAMPLE}_kir_genotype.tsv" 2>/dev/null | head -n 20
check "database_release.txt names IPD-KIR ${KIR_DB_RELEASE}" \
  grep -q "^database: IPD-KIR Release Version ${KIR_DB_RELEASE//./\\.}$" "${K}/database_release.txt"
check "the KIR genotype file is written" test -s "${K}/${SAMPLE}_kir_genotype.tsv"
check "it holds KIR genotypes or says too few reads" \
  grep -Eq $'^KIR[0-9A-Z]+\t|^# KIR not typed: T1K found too few reads' "${K}/${SAMPLE}_kir_genotype.tsv"
check "the HLA pass still wrote its genotypes" test -s "${GENOME_DIR}/${SAMPLE}/hla_t1k/${SAMPLE}_hla_genotype.tsv"
echo "- KIR on the fixture: $(grep -c $'^KIR' "${K}/${SAMPLE}_kir_genotype.tsv" 2>/dev/null || echo 0) gene rows typed" >> "$E2E_NOTES"

finish
