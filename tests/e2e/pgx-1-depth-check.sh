#!/usr/bin/env bash
# The CYP2D6 depth check (bin/cyp2d6_depth_check.py) before pypgx and Cyrius:
# step 32 (case 38) and step 21 (case 37) wrote it. On the fixture BAM step 02
# aligned to the no-ALT reference the reads at CYP2D6 map uniquely, so step
# 32's check passes and pypgx's CYP2D6 row keeps its call. (The with-ALT side is
# the ALT depth A/B job: scripts/ci/alt-depth-ab.sh runs the same check on its
# with-ALT mapping and requires it flagged.)
. "$(dirname "$0")/lib.sh"

CHECK="${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_cyp2d6_depth_check.tsv"
cat "$CHECK" 2>/dev/null
val() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "$CHECK" 2>/dev/null; }
check_eq "step 32's CYP2D6 depth check on the no-ALT BAM" "$(val status)" ok
check "CYP2D6 has reads (all reads)" awk -v d="$(val gene_depth_all)" 'BEGIN {exit !(d + 0 > 1)}'
check "the flanks have reads (MAPQ >= 1)" awk -v d="$(val flank_depth_mapq1)" 'BEGIN {exit !(d + 0 > 1)}'
CYP2D6=$(awk -F'\t' '$1 == "CYP2D6" {print $2; exit}' "${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv" 2>/dev/null)
check "pypgx's CYP2D6 row keeps its call (${CYP2D6:-none})" lacks '^(Indeterminate|)$' "${CYP2D6:-}"

# Step 21 (case 37) runs on the fixture's HG002_cyrius.bam: GIAB's alignment
# of the regions Cyrius reads plus the two 50 kb flanks the check compares
# CYP2D6 with (since fixture-v6). GIAB aligned it to a reference without ALT
# contigs, so the check passes and step 21 keeps the Cyrius call unmarked.
CY="${SAMPLE}cyrius"
CHECK21="${GENOME_DIR}/${CY}/cyrius/${CY}_cyp2d6_depth_check.tsv"
cat "$CHECK21" 2>/dev/null
val21() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "$CHECK21" 2>/dev/null; }
check_eq "step 21's CYP2D6 depth check on the Cyrius BAM" "$(val21 status)" ok
check "the Cyrius BAM's flanks have reads (MAPQ >= 1)" awk -v d="$(val21 flank_depth_mapq1)" 'BEGIN {exit !(d + 0 > 1)}'
FILTER21=$(awk -F'\t' 'NR > 1 {print $3; exit}' "${GENOME_DIR}/${CY}/cyrius/${CY}_cyp2d6.tsv" 2>/dev/null)
check "step 21 has a Cyrius row (Filter ${FILTER21:-none})" test -n "${FILTER21:-}"
check "step 21 did not mark the Cyrius call (Filter ${FILTER21:-none})" lacks 'CYP2D6_depth_unreliable' "${FILTER21:-}"

finish
