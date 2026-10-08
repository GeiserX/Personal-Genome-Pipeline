#!/usr/bin/env bash
# Step 10 (TelomereHunter) with the GRCh38 chromosome bands setup.sh installs.
# TelomereHunter runs as the caller (no run_in --root). It runs on a
# copy of the sample under another name (the BAM hard-linked), removed at the
# end, so the cases after it see the sample directory as before.
. "$(dirname "$0")/lib.sh"

check "the GRCh38 bands install" bash -c '. "$1/scripts/lib/common.sh" && install_data_file cytoband' _ "$REPO"
BANDS="${GENOME_DIR}/reference/cytoBand.hg38.txt"
check_ge "band lines on chr1-22, X and Y" "$(grep -cE '^chr([0-9]+|X|Y)[[:space:]]' "$BANDS" 2>/dev/null || true)" 800
check_eq "band lines on other contigs" "$(grep -cvE '^chr([0-9]+|X|Y)[[:space:]]' "$BANDS" 2>/dev/null || true)" 0

T="${SAMPLE}tel"
mkdir -p "${GENOME_DIR}/${T}/aligned"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${T}/aligned/${T}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${T}/aligned/${T}_sorted.bam.bai"
run_step 10-telomere-hunter.sh "$T"
check_step_exit 10-telomere-hunter.sh
check "the log names the GRCh38 bands" has "Chromosome bands: ${BANDS}" "$(cat "$STEP_LOG")"
SUM="${GENOME_DIR}/${T}/telomere/${T}/${T}/${T}_summary.tsv"
head -n 3 "$SUM" 2>/dev/null
check_ge "summary rows" "$(grep -c . "$SUM" 2>/dev/null || true)" 2
check "the summary has a tel_content column" has 'tel_content' "$(head -n 1 "$SUM" 2>/dev/null)"
in_genome "$BCFTOOLS_IMAGE" rm -rf "$T"

finish
