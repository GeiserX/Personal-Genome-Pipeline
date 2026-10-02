#!/usr/bin/env bash
# Step 21 (Cyrius) runs and writes its TSV, whatever the CYP2D6 call is.
#
# Cyrius normalises depth over 3,000 bins on chr1-chr22 before it calls, and
# stops on the first bin whose contig the BAM lacks. The BAM step 02 writes
# has only the fixture's contigs, so this case gives step 21 the fixture's
# HG002_cyrius.bam instead: GIAB's alignment of exactly the regions Cyrius
# reads, under its own sample name.
. "$(dirname "$0")/lib.sh"

CY="${SAMPLE}cyrius"
mkdir -p "${GENOME_DIR}/${CY}/aligned"
cp "${FIXTURE_DIR}/${SAMPLE}_cyrius.bam" "${GENOME_DIR}/${CY}/aligned/${CY}_sorted.bam"
cp "${FIXTURE_DIR}/${SAMPLE}_cyrius.bam.bai" "${GENOME_DIR}/${CY}/aligned/${CY}_sorted.bam.bai"

run_step 21-cyrius.sh "$CY"
check_step_exit 21-cyrius.sh

TSV="${GENOME_DIR}/${CY}/cyrius/${CY}_cyp2d6.tsv"
check_ge "lines in the Cyrius TSV (header and sample)" "$(grep -c . "$TSV" 2>/dev/null || true)" 2
check "the Cyrius TSV has a row for ${CY}" has "^${CY}" "$(awk 'NR > 1' "$TSV" 2>/dev/null)"
[ -s "$TSV" ] && cat "$TSV"

finish
