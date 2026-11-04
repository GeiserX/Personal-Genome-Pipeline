#!/usr/bin/env bash
# 53-locityper.sh — A second opinion at hard polymorphic genes with Locityper
# Usage: ./scripts/53-locityper.sh <sample_name>
#
# Opt-in. HLA class II, KIR, CYP2D6 and FCGR each have one caller here or
# none, so a wrong call there cannot be seen. Locityper (Prodanov et al. 2025,
# Nat Genet, PMID 41107551; MIT licence) genotypes a whole locus at once: it
# picks the two haplotypes of a pangenome panel (here the HPRC v2.1 panel of
# 232 people) that best explain the reads there and their depth. Its call is
# shown beside the primary caller's, never as the call. bin/locityper_report.py,
# the same code as the LOCITYPER Nextflow process, writes one row per locus of
# assets/locityper_loci.bed.
#
# Input:  aligned/<sample>_sorted.bam with its .bai (ALIGN_DIR as for the other
#         BAM steps), aligned to REF_FASTA; Locityper's database
#         (setup.sh --locityper-data: reference/locityper-<LOCITYPER_PANEL_VERSION>/,
#         LOCITYPER_DATA overrides), built from the same REF_FASTA.
# Output: locityper/<sample>_locityper.tsv   one row per locus: the two
#                                            haplotypes, quality, reads, the
#                                            authors' filter
#         locityper/<sample>_locityper.log   Locityper's log
#         locityper/work/                    preproc/ (the sample's read depth,
#                                            insert size and error profile) and
#                                            gt/loci/<locus>/res.json.gz
# LOCITYPER_BG_REGION: the background region preproc measures depth on; by
# default Locityper's own for GRCh38 (chr17:72062001-76562000).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image LOCITYPER_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
DATA=${LOCITYPER_DATA:-${GENOME_DIR}/reference/locityper-${LOCITYPER_PANEL_VERSION}}
case "$DATA" in /*) ;; *) DATA="${GENOME_DIR}/${DATA}" ;; esac
OUTDIR="${S}/locityper"
WORK="${OUTDIR}/work"
TABLE="${OUTDIR}/${SAMPLE}_locityper.tsv"
LOG="${OUTDIR}/${SAMPLE}_locityper.log"
THREADS=${THREADS:-8}

echo "=== Locityper second opinion at hard loci: ${SAMPLE} ==="
for f in "$BAM" "${BAM}.bai" "$REF_FASTA" "${REF_FASTA}.fai"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done
for f in "${DATA}/DONE" "${DATA}/counts.jf" "${DATA}/loci.bed"; do
  [ -s "$f" ] || { echo "ERROR: Locityper's database not found: ${f} (run: ./scripts/setup.sh --locityper-data ${GENOME_DIR})" >&2; exit 1; }
done
# The database's k-mer counts and haplotypes come from one reference: DONE names its .fai.
FAI_SHA=$(_digest sha256 "${REF_FASTA}.fai")
grep -q "fai=${FAI_SHA} " "${DATA}/DONE" || {
  echo "ERROR: Locityper's database ${DATA} was built from another reference than ${REF_FASTA} (run: ./scripts/setup.sh --locityper-data ${GENOME_DIR})" >&2; exit 1; }
echo "BAM: ${BAM}"
echo "Database: ${DATA} ($(grep -vc '^#' "${DATA}/loci.bed") loci)"
mkdir -p "$OUTDIR"
rm -rf "$WORK" "$TABLE" "$LOG"
mkdir -p "$WORK"

echo "[1/2] Locityper preproc and genotype (${THREADS} threads)..."
# Locityper's exit code goes to the table, not to this script: a run that
# stops leaves its loci "not called", with the log beside the table.
rc=0
# shellcheck disable=SC2016  # $1 to $7 belong to the inner bash
run_in --cpus "$THREADS" --memory 16g "$LOCITYPER_IMAGE" bash -c '
    { locityper preproc -a "$1" -r "$2" -j "$3/counts.jf" -o "$4/preproc" -@ "$5" ${6:+-b "$6"} \
        && locityper genotype -a "$1" -d "$3/db" -p "$4/preproc" -o "$4/gt" -@ "$5"; } > "$7" 2>&1' \
  _ "$(cpath "$BAM")" "$REF_FASTA_C" "$(cpath "$DATA")" "$(cpath "$WORK")" "$THREADS" \
  "${LOCITYPER_BG_REGION:-}" "$(cpath "$LOG")" || rc=$?
[ "$rc" -eq 0 ] || { echo "WARNING: Locityper exited ${rc}; last lines of ${LOG}:" >&2; tail -n 20 "$LOG" >&2 || true; }

echo "[2/2] The table (bin/locityper_report.py)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g "$LOCITYPER_IMAGE" \
  python3 /pgp-bin/locityper_report.py --gt-dir "$(cpath "${WORK}/gt")" --loci "$(cpath "${DATA}/loci.bed")" \
    --exit-code "$rc" --out "$(cpath "$TABLE")"

[ -s "$TABLE" ] || { echo "ERROR: no table written (${TABLE})" >&2; exit 1; }
echo ""
awk -F'\t' 'NR > 1 {printf "  %-9s %-7s %s %s\n", $1, $3, ($4 == "." ? "" : $4), ($12 == "." ? "" : "(filter " $12 ")")}' "$TABLE"
echo "Table: ${TABLE}"
