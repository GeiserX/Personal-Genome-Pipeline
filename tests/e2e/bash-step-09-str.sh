#!/usr/bin/env bash
# Step 09 (ExpansionHunter) and step 09b (Stranger with its GRCh38 catalog).
# Stranger writes its output through a temporary name and annotates every
# ExpansionHunter locus it finds in the catalog.
. "$(dirname "$0")/lib.sh"

run_step 09-expansion-hunter.sh "$SAMPLE" male
check_step_exit 09-expansion-hunter.sh
EH="${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf"
N_EH=$(grep -vc '^#' "${GENOME_DIR}/${EH}" 2>/dev/null || true)
check_ge "ExpansionHunter loci" "$N_EH" 20

run_step 09b-stranger.sh "$SAMPLE"
check_step_exit 09b-stranger.sh
check "the log names the GRCh38 catalog" has 'variant_catalog_grch38' "$(cat "$STEP_LOG")"
OUT="${GENOME_DIR}/${SAMPLE}/expansion_hunter/${SAMPLE}_eh_stranger.vcf"
check "no temporary file is left" test ! -e "${OUT}.tmp"
check "the output declares STR_STATUS" grep -q '^##INFO=<ID=STR_STATUS' "$OUT"
check_eq "records in = records out" "$(grep -vc '^#' "$OUT" 2>/dev/null || true)" "$N_EH"
check_ge "records with the catalog's pathologic threshold" "$(grep -v '^#' "$OUT" 2>/dev/null | grep -c 'STR_PATHOLOGIC_MIN=' || true)" 20

finish
