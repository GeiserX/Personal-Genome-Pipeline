#!/usr/bin/env bash
# 52-vntyper.sh — MUC1 VNTR frameshifts (ADTKD-MUC1) with VNtyper
# Usage: ./scripts/52-vntyper.sh <sample_name>
#
# Opt-in. ADTKD-MUC1, a dominant kidney disease, is caused by frameshifts
# inside the MUC1 coding VNTR, a repeat of 60 bp units, most often one extra C
# (dupC). Short reads cannot be placed inside the repeat, so DeepVariant and
# ExpansionHunter do not see it. VNtyper (VNTYPER_IMAGE, its reference data
# inside) genotypes it with k-mers; bin/vntyper_report.py, the same code as
# the VNTYPER Nextflow process, turns its result into positive, possible,
# negative or not assessed, with the share of known positives VNtyper found
# at the sample's depth over the VNTR. A research call, not a diagnosis.
#
# Input:  aligned/<sample>_sorted.bam with its .bai (ALIGN_DIR as for the other
#         BAM steps), GRCh38 with chr contig names.
# Output: vntyper/<sample>_vntyper.tsv   metric<TAB>value: status, the call,
#                                        the depth, the sensitivity line, message
#         vntyper/run/                   VNtyper's own folder (summary_report.html,
#                                        kestrel/, coverage/, pipeline.log)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image VNTYPER_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
OUTDIR="${S}/vntyper"
RUN="${OUTDIR}/run"
TABLE="${OUTDIR}/${SAMPLE}_vntyper.tsv"
LOG="${OUTDIR}/${SAMPLE}_vntyper.log"
THREADS=${THREADS:-4}

echo "=== MUC1 VNTR frameshifts (VNtyper): ${SAMPLE} ==="
[ -f "$BAM" ] || { echo "ERROR: File not found: ${BAM}" >&2; exit 1; }
[ -f "${BAM}.bai" ] || { echo "ERROR: File not found: ${BAM}.bai" >&2; exit 1; }
echo "BAM: ${BAM}"
mkdir -p "$OUTDIR"
rm -rf "$RUN" "$TABLE" "$LOG"

echo "[1/2] VNtyper on the MUC1 region (fast mode, GRCh38)..."
# The image's entrypoint takes only its own command words, so it is replaced.
# VNtyper reads its references through paths relative to /opt/vntyper, so it
# runs from there. Its exit code goes to the table, not to this script: a run
# that stops says "not assessed" with the log beside it.
rc=0
# shellcheck disable=SC2016  # $1 to $5 belong to the inner bash
run_in --cpus "$THREADS" --memory 8g --entrypoint "" "$VNTYPER_IMAGE" bash -c '
    cd /opt/vntyper && vntyper pipeline --bam "$1" -o "$2" --sample-name "$3" --threads "$4" \
      --fast-mode --reference-assembly hg38 > "$5" 2>&1' \
  _ "$(cpath "$BAM")" "$(cpath "$RUN")" "$SAMPLE" "$THREADS" "$(cpath "$LOG")" || rc=$?
[ "$rc" -eq 0 ] || { echo "WARNING: VNtyper exited ${rc}; last lines of ${LOG}:" >&2; tail -n 20 "$LOG" >&2 || true; }

echo "[2/2] The reading (bin/vntyper_report.py)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g --entrypoint "" "$VNTYPER_IMAGE" \
  python3 /pgp-bin/vntyper_report.py --run-dir "$(cpath "$RUN")" --exit-code "$rc" --log "$(cpath "$LOG")" --out "$(cpath "$TABLE")"

[ -s "$TABLE" ] || { echo "ERROR: no table written (${TABLE})" >&2; exit 1; }
echo ""
awk -F'\t' '$1 == "status" || $1 == "message" || $1 == "sensitivity" {print "  " $2}' "$TABLE"
echo "Table: ${TABLE}"
echo "VNtyper's report: ${RUN}/summary_report.html"
