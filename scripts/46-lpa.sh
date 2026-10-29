#!/usr/bin/env bash
# 46-lpa.sh — Lp(a): the LPA KIV-2 repeat count estimated from read depth
# Usage: ./scripts/46-lpa.sh <sample_name>
#
# Opt-in. The number of kringle IV type 2 (KIV-2) repeats in LPA sets the size
# of apo(a); a low total goes with a higher Lp(a) level, a mostly inherited
# cardiovascular risk factor. The step counts the repeats from the depth over
# the KIV-2 exons against the depth over single-copy LPA exons
# (bin/copy_number.py on the LPA_KIV2 row of bin/cn_loci.tsv, laid out by
# bin/lpa_kiv2.py; the same code as the LPA Nextflow processes), with a 95%
# interval. The two common risk variants, rs10455872 and rs3798220, are rows
# of step 42's traits panel; the report's Lp(a) section reads both. A research
# estimate, not a diagnostic test: the actionable test is a blood Lp(a)
# measurement.
#
# Input:  aligned/<sample>_sorted.bam with its .bai (ALIGN_DIR as for the other
#         BAM steps), from an aligner that reports reads with several equally
#         good places (minimap2, BWA-MEM): with no read of MAPQ 0 over KIV-2
#         the count is not assessed.
# Output: lpa/<sample>_lpa_kiv2.tsv          metric<TAB>value: status, the count,
#                                            its 95% interval, the depths, message
#         lpa/<sample>_lpa_kiv2_depth.tsv    samtools bedcov -j over the regions (all MAPQ)
#         lpa/<sample>_lpa_kiv2_depth_q1.tsv the same, reads of MAPQ 1 or more
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image PYTHON_IMAGE SAMTOOLS_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
OUTDIR="${S}/lpa"
TABLE="${OUTDIR}/${SAMPLE}_lpa_kiv2.tsv"
DEPTH="${OUTDIR}/${SAMPLE}_lpa_kiv2_depth.tsv"
DEPTH_Q1="${OUTDIR}/${SAMPLE}_lpa_kiv2_depth_q1.tsv"
BED="${OUTDIR}/regions.bed"

echo "=== Lp(a): LPA KIV-2 repeat count from read depth: ${SAMPLE} ==="
[ -f "$BAM" ] || { echo "ERROR: File not found: ${BAM}" >&2; exit 1; }
[ -f "${BAM}.bai" ] || { echo "ERROR: File not found: ${BAM}.bai" >&2; exit 1; }
echo "BAM: ${BAM}"
mkdir -p "$OUTDIR"
rm -f "$TABLE" "$DEPTH" "$DEPTH_Q1" "$BED"

echo "[1/3] The regions (bin/copy_number.py bed --locus LPA_KIV2)..."
# shellcheck disable=SC2016  # $1 belongs to the inner shell
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g "$PYTHON_IMAGE" \
  sh -c 'python3 /pgp-bin/copy_number.py bed --locus LPA_KIV2 > "$1"' _ "$(cpath "$BED")"

echo "[2/3] Read depth over the KIV-2 and single-copy LPA exons (samtools bedcov -j, all MAPQ and MAPQ >= 1)..."
# bedcov stops on a contig the BAM has no @SQ line for: without chr6 both
# tables stay empty and the count is not assessed. The LPA_DEPTH process runs
# the same commands.
# shellcheck disable=SC2016  # $1 to $5 belong to the inner bash, the program to awk
run_in --cpus 1 --memory 2g "$SAMTOOLS_IMAGE" bash -euo pipefail -c '
    if samtools view -H "$2" | awk -F "\t" "$5"; then
      samtools bedcov -j "$1" "$2" > "$3"
      samtools bedcov -j -Q 1 "$1" "$2" > "$4"
    else : > "$3"; : > "$4"; fi' \
  _ "$(cpath "$BED")" "$(cpath "$BAM")" "$(cpath "$DEPTH")" "$(cpath "$DEPTH_Q1")" \
  '$1 == "@SQ" {for (i = 2; i <= NF; i++) if ($i == "SN:chr6") f = 1} END {exit !f}'
rm -f "$BED"

echo "[3/3] The count and its interval (bin/copy_number.py call --locus LPA_KIV2)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g "$PYTHON_IMAGE" \
  python3 /pgp-bin/copy_number.py call --locus LPA_KIV2 --bedcov "$(cpath "$DEPTH")" \
    --bedcov-mapq1 "$(cpath "$DEPTH_Q1")" --table "$(cpath "$TABLE")"

[ -s "$TABLE" ] || { echo "ERROR: no table written (${TABLE})" >&2; exit 1; }
echo ""
echo "Table: ${TABLE}"
