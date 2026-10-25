#!/usr/bin/env bash
# 47-viral.sh — Viral DNA in the reads that do not map to the human genome
# Usage: ./scripts/47-viral.sh <sample_name>
#
# Opt-in. Takes the reads the human alignment left unmapped, and the reads on
# chrEBV (GRCh38's analysis set carries the EBV genome, so EBV reads map there),
# realigns them with minimap2 to the human viral genomes setup.sh --viral-data
# installs (HHV-6A and -6B, the other herpesviruses, HPV, HBV, HCV, HIV-1,
# HTLV-1, polyomaviruses, adenovirus 5, an anellovirus, parvovirus B19 and
# Illumina's phiX control), and reports per genome the reads, the breadth and
# the depth against the human autosomal depth. Inherited chromosomally
# integrated HHV-6 (in about 1 in 100 people) shows as the whole HHV-6 genome at
# about half the human depth. bin/viral_report.py applies the rules; it is the
# same code as the VIRAL_REPORT Nextflow process. Not a test for an infection.
#
# Input:  aligned/<sample>_sorted.bam with its .bai (ALIGN_DIR as for the other
#         BAM steps); reference/viral-<VIRAL_REF_VERSION>/viral_genomes.fa
#         (setup.sh --viral-data; VIRAL_REF overrides, inside GENOME_DIR)
# Output: viral/<sample>_viral.tsv            one row per viral genome, with the call
#         viral/<sample>_viral_reads.bam(.bai) the realigned reads that hit a viral genome
#         viral/<sample>_viral_coverage.tsv   samtools coverage of that BAM
#         viral/<sample>_human_coverage.tsv   samtools coverage of the sample's BAM
# The sample's BAM is read twice end to end: once for its unmapped reads, once
# for the human depth.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image SAMTOOLS_IMAGE MINIMAP2_IMAGE PYTHON_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
VREF=${VIRAL_REF:-${GENOME_DIR}/reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa}
case "$VREF" in /*) ;; *) VREF="${GENOME_DIR}/${VREF}" ;; esac
OUTDIR="${S}/viral"
P="${OUTDIR}/${SAMPLE}"
TABLE="${P}_viral.tsv"

echo "=== Viral DNA in the unmapped reads: ${SAMPLE} ==="
for f in "$BAM" "${BAM}.bai"; do [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }; done
[ -s "$VREF" ] || { echo "ERROR: viral genomes not found: ${VREF} (run: ./scripts/setup.sh --viral-data ${GENOME_DIR})" >&2; exit 1; }
echo "BAM: ${BAM}"
echo "Viral genomes: ${VREF}"
mkdir -p "$OUTDIR"
rm -f "$TABLE" "${P}_viral_reads.bam" "${P}_viral_reads.bam.bai" "${P}_viral_coverage.tsv" "${P}_human_coverage.tsv" \
  "${P}_unmapped.fq" "${P}_viral.sam" "${P}_viral_depth.tsv"
B=$(cpath "$BAM") O=$(cpath "$P") V=$(cpath "$VREF")

# The same commands as the VIRAL_READS process. samtools fastq leaves out
# secondary and supplementary records, so each read is realigned once; a read
# unmapped next to a mapped mate is taken too (an integration junction). The
# reads mapped to chrEBV are taken mapped (-F 4): the unmapped ones are in the
# first set already. Both coverages count duplicates (--ff without DUP), so the
# viral and the human depth are measured alike.
echo "[1/4] Unmapped reads and chrEBV reads (samtools), and the human depth (samtools coverage)..."
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
run_in --cpus "$THREADS" --memory 4g "$SAMTOOLS_IMAGE" bash -euo pipefail -c '
    samtools view -@ "$3" -u -f 4 "$1" | samtools fastq - > "$2_unmapped.fq"
    if samtools view -H "$1" | awk -F "\t" "\$1 == \"@SQ\" && \$2 == \"SN:chrEBV\" {f = 1} END {exit !f}"; then
      samtools view -@ "$3" -u -F 4 "$1" chrEBV | samtools fastq - >> "$2_unmapped.fq"
    fi
    samtools coverage --ff UNMAP,SECONDARY,QCFAIL "$1" > "$2_human_coverage.tsv"' _ "$B" "$O" "$THREADS"
echo "  reads to realign: $(( $(wc -l < "${P}_unmapped.fq") / 4 ))"

echo "[2/4] minimap2 against the viral genomes..."
run_in --cpus "$THREADS" --memory 4g "$MINIMAP2_IMAGE" \
  minimap2 -ax sr --frag=no --secondary=no -t "$THREADS" -o "${O}_viral.sam" "$V" "${O}_unmapped.fq"

echo "[3/4] Sort, coverage and depth over each viral genome (samtools)..."
# shellcheck disable=SC2016
run_in --cpus 2 --memory 4g "$SAMTOOLS_IMAGE" bash -euo pipefail -c '
    samtools sort -o "$1_viral_reads.bam" "$1_viral.sam"
    samtools index "$1_viral_reads.bam"
    samtools coverage --ff UNMAP,SECONDARY,QCFAIL "$1_viral_reads.bam" > "$1_viral_coverage.tsv"
    samtools depth -aa -G UNMAP,SECONDARY,QCFAIL "$1_viral_reads.bam" > "$1_viral_depth.tsv"' _ "$O"

echo "[4/4] Calls (bin/viral_report.py)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
  python3 /pgp-bin/viral_report.py --sample "$SAMPLE" --viral-fasta "$V" \
    --viral-coverage "${O}_viral_coverage.tsv" --viral-depth "${O}_viral_depth.tsv" \
    --human-coverage "${O}_human_coverage.tsv" --out "$(cpath "$TABLE")"
rm -f "${P}_unmapped.fq" "${P}_viral.sam" "${P}_viral_depth.tsv"

[ -s "$TABLE" ] || { echo "ERROR: no table written (${TABLE})" >&2; exit 1; }
echo ""
echo "Table: ${TABLE}"
echo "Viral reads: ${P}_viral_reads.bam"
