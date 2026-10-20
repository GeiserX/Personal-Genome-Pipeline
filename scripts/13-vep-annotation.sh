#!/usr/bin/env bash
# VEP — Ensembl Variant Effect Predictor
# Full functional annotation: consequence, SIFT, PolyPhen, regulatory, etc.
# Requires: VEP cache (~26 GB download, one-time)
# Output: ${SAMPLE}_vep.vcf.gz and its .tbi in $GENOME_DIR/<sample>/vep/
#
# Same annotation as the Nextflow VEP module: --everything with the reference
# FASTA, without which VEP turns HGVS off. A finished run removes the files the
# later steps build from an older _vep output (vcfanno's _annotated.vcf.gz,
# the uncompressed _vep.vcf of earlier versions), so steps 30, 23 and 31
# rebuild from this one instead of reading the old annotation.
#
# PASS only: VEP reads the records with FILTER PASS or '.' (a caller that
# writes no FILTER), from a temporary copy made with bcftools. Steps 23 and
# 31, which read this output directly or through vcfanno (30), keep PASS
# records alone, so the others were annotated for nothing; the Nextflow VEP
# module selects the same records.
#
# Normalised: the PASS copy goes through bcftools norm -m -any -f <reference>
# (multiallelic records split, indels left-aligned) before VEP, because VEP
# matches an allele to gnomAD and ClinVar as those files write it. A deletion
# written right-shifted, or inside a multiallelic record, otherwise finds no
# known variant and gets no frequency, which step 23 would show as
# "frequency unknown" even for a common allele. Records on a contig the
# reference lacks are left out (bcftools norm cannot read their bases), and
# the step says how many contigs that is. The Nextflow VEP module does not
# normalise yet (docs/13-vep-annotation.md).
#
# ClinVar: when clinvar/clinvar_pathogenic_chr.vcf.gz is installed (setup.sh,
# the file step 06 screens against), VEP also annotates it with --custom: its
# CLNSIG, CLNREVSTAT and CLNDN become the CSQ fields ClinVar_CLNSIG,
# ClinVar_CLNREVSTAT and ClinVar_CLNDN. Step 23's ClinVar tier then follows a
# refresh of that file, not the ClinVar of the cache release (CLIN_SIG).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF_DIR=${VCF_DIR:-vcf}
VCF="${GENOME_DIR}/${SAMPLE}/${VCF_DIR}/${SAMPLE}.vcf.gz"
CACHE_DIR="${GENOME_DIR}/vep_cache"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/vep"

echo "=== VEP Annotation: ${SAMPLE} ==="

# The index too: a VCF without one may be half written.
for f in "$VCF" "${VCF}.tbi" "$REF_FASTA" "${REF_FASTA}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# The cache must be the release of VEP_IMAGE (VEP_CACHE_RELEASE in
# versions.env). Another release in the same directory, such as the one CPSR
# uses, does not count.
if [ ! -f "${CACHE_DIR}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
  echo "VEP ${VEP_CACHE_RELEASE} cache not found in ${CACHE_DIR}. Installing (26 GB download)..."
  install_vep_cache "$CACHE_DIR" "$VEP_CACHE_RELEASE"
  echo "Cache installed at ${CACHE_DIR}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/"
fi

OUT="${OUTPUT_DIR}/${SAMPLE}_vep.vcf.gz"
CLINVAR="${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz"
CUSTOM=()
if [ -f "$CLINVAR" ] && [ -f "${CLINVAR}.tbi" ]; then
  CUSTOM=(--custom "file=$(cpath "$CLINVAR"),short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN")
  echo "ClinVar: ${CLINVAR} (as ClinVar_CLNSIG; the file step 06 screens against)"
else
  echo "ClinVar: ${CLINVAR} is not installed; only the cache's CLIN_SIG (run setup.sh for the current file)"
fi
OUT_C="/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"

# Run VEP into a temporary name, so a run that fails or is killed leaves the
# previous result in place. The cache stays writable as before: VEP builds an
# index for a FASTA it finds there without one, and no CI run shows it never
# writes. --compress_output bgzip: the whole-genome output is written
# compressed, ready for tabix.
rm -f "${OUT}.tmp.vcf.gz"
PASS_VCF="${OUTPUT_DIR}/${SAMPLE}.pass.tmp.vcf.gz"
NORM_VCF="${OUTPUT_DIR}/${SAMPLE}.norm.tmp.vcf.gz"
REF_TARGETS="${OUTPUT_DIR}/${SAMPLE}.ref_contigs.tmp.tsv"
trap 'rm -f "$PASS_VCF" "$NORM_VCF" "$REF_TARGETS"' EXIT
run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  bcftools view -f PASS,. -Oz -o "$(cpath "$PASS_VCF")" "/genome/${SAMPLE}/${VCF_DIR}/${SAMPLE}.vcf.gz"
# Split and left-align against the reference (see the header). -T keeps the
# contigs the reference holds: bcftools norm stops at the first record whose
# contig the FASTA lacks. -c w warns on a REF base that differs, as step 06.
awk -v OFS='\t' '{print $1, 1, 2147483647}' "${REF_FASTA}.fai" > "$REF_TARGETS"
OFF_REF=$(run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
    bcftools index -s "/genome/${SAMPLE}/${VCF_DIR}/${SAMPLE}.vcf.gz" \
  | awk -F'\t' 'NR == FNR {ref[$1] = 1; next} $3 > 0 && !($1 in ref) {n++} END {print n + 0}' "${REF_FASTA}.fai" -)
if [ "${OFF_REF:-0}" -gt 0 ]; then
  echo "NOTICE: ${OFF_REF} contig(s) of ${VCF} are not in ${REF_FASTA}; their records are not annotated."
fi
run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  bcftools norm -m -any -c w -f "$REF_FASTA_C" -T "$(cpath "$REF_TARGETS")" -Oz -o "$(cpath "$NORM_VCF")" "$(cpath "$PASS_VCF")"
mv -f "$NORM_VCF" "$PASS_VCF"
# Each VEP fork holds its own copy of the cache index: 2 GB per thread, 8 GB at least.
VEP_MEM=$(( THREADS * 2 > 8 ? THREADS * 2 : 8 ))
run_in \
  --cpus "${THREADS}" --memory "${VEP_MEM}g" \
  -v "${CACHE_DIR}:/opt/vep/.vep" \
  "${VEP_IMAGE}" \
  vep \
    --input_file "$(cpath "$PASS_VCF")" \
    -o "${OUT_C}.tmp.vcf.gz" \
    --vcf \
    --compress_output bgzip \
    --cache \
    --cache_version "${VEP_CACHE_RELEASE}" \
    --dir_cache /opt/vep/.vep \
    --offline \
    --assembly GRCh38 \
    --fasta "${REF_FASTA_C}" \
    --everything \
    --force_overwrite \
    --stats_file "/genome/${SAMPLE}/vep/${SAMPLE}_vep_summary.html" \
    --warning_file "/genome/${SAMPLE}/vep/${SAMPLE}_vep_warnings.txt" \
    --fork "${THREADS}" \
    ${CUSTOM[@]+"${CUSTOM[@]}"}

if ! wrote_vcf "${OUT}.tmp.vcf.gz"; then
  rm -f "${OUT}.tmp.vcf.gz"
  echo "ERROR: VEP exited without a complete output." >&2
  exit 1
fi
rm -f "${OUT}.tbi"
mv -f "${OUT}.tmp.vcf.gz" "$OUT"
run_in "$BCFTOOLS_IMAGE" bcftools index -f -t "$OUT_C"

# Outputs derived from an older VEP run: steps 30, 23 and 31 would read them
# instead of this one.
rm -f "${OUTPUT_DIR}/${SAMPLE}_annotated.vcf.gz" "${OUTPUT_DIR}/${SAMPLE}_annotated.vcf.gz.tbi" \
  "${OUTPUT_DIR}/${SAMPLE}_vep.vcf"

echo "=== VEP complete ==="
echo "Results: ${OUT}"
echo ""
echo "Filter HIGH impact variants:"
echo "  gzip -dc ${OUT} | grep 'HIGH' | head"
