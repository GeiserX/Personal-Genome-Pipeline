#!/usr/bin/env bash
# pgp-lookup.sh — what is my genotype at an rsID, a position, or in a gene?
# Usage: ./scripts/pgp-lookup.sh <sample_name> <rsID|chr:pos|GENE> [more queries...]
#
# Answers from the files the pipeline already wrote, without opening a VCF
# by hand (docs/lookup.md):
#   rs4988235        the genotype at that rsID
#   chr2:135851076   the genotype at that GRCh38 position
#   CFTR             every variant called in the gene
# A site gets its genotype from the gVCF of step 3: a position inside a
# reference block is ref/ref, a position no record covers is a no call. Then
# its depth and quality, and, when the annotated VCF of step 13 (or 30) holds
# a variant there, its gene, consequence, ClinVar significance and gnomAD
# frequency. A gene gets every variant called in it with the same columns,
# and how many of its bases have a call. Both end with what the reports
# (summary.json, from generate-report.sh) say about the gene or rsID: PGx,
# ClinVar, carrier screen, traits.
#
# rsIDs are found in the consumer site table (setup.sh --consumer-sites,
# 614,138 chip rsIDs with dbSNP b157 GRCh38 positions), then in
# assets/traits.tsv and assets/bloodtype_markers.tsv, then in the annotated
# VCF (only the sample's own variants carry an rsID there). Genes are found in
# GENCODE's gene lines, which setup.sh installs.
#
# Nextflow output has the same per-sample folders: run it with
# GENOME_DIR=<outdir>. The rsID table and the gene lines live in the
# reference folder of the data directory; when it is elsewhere, point
# REF_DIR at it.
#
# Read-only: it prints the answer and writes nothing but a temporary slice of
# the VCFs, removed when it ends.
#
# Input:  vcf/<sample>.g.vcf.gz (step 3; the VCF if there is no gVCF),
#         vep/<sample>_annotated.vcf.gz or vep/<sample>_vep.vcf.gz (optional),
#         summary.json (optional)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <rsID|chr:pos|GENE> [more queries...]}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory (the Nextflow --outdir for Nextflow output)}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
shift
if [ $# -eq 0 ]; then
  echo "Usage: $0 <sample_name> <rsID|chr:pos|GENE> [more queries...]" >&2
  exit 2
fi
require_image BCFTOOLS_IMAGE PYTHON_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
REF_DIR=${REF_DIR:-${GENOME_DIR}/reference}
REF_DIR=${REF_DIR%/}

GVCF="${S}/vcf/${SAMPLE}.g.vcf.gz"
if [ ! -s "$GVCF" ]; then
  if [ -s "${S}/vcf/${SAMPLE}.vcf.gz" ]; then
    GVCF="${S}/vcf/${SAMPLE}.vcf.gz"
    echo "NOTE: no gVCF for ${SAMPLE}; reading ${GVCF}. A VCF has no reference blocks, so a site where you match"
    echo "  the reference reads as a no call here. Step 3 writes the gVCF."
  else
    echo "ERROR: no gVCF or VCF for ${SAMPLE} under ${S}/vcf/: run scripts/03-deepvariant.sh ${SAMPLE} first." >&2
    exit 1
  fi
fi
ANN=""
for f in "${S}/vep/${SAMPLE}_annotated.vcf.gz" "${S}/vep/${SAMPLE}_vep.vcf.gz"; do
  if [ -s "$f" ]; then ANN=$f; break; fi
done
SUMMARY="${S}/summary.json"
[ -s "$SUMMARY" ] || SUMMARY=""
SITES="${REF_DIR}/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
[ -s "$SITES" ] || SITES=""
GENES="${REF_DIR}/gencode.v${GENCODE_RELEASE}.basic.genes.gtf"
[ -s "$GENES" ] || GENES=""

echo "=== Genome lookup: ${SAMPLE} ==="
echo "Genotypes:  ${GVCF}"
echo "Annotation: ${ANN:-none (step 13 writes vep/${SAMPLE}_vep.vcf.gz)}"
echo "Reports:    ${SUMMARY:-none (generate-report.sh writes summary.json)}"
echo ""

# The containers see GENOME_DIR at /genome and, when REF_DIR is elsewhere,
# REF_DIR read-only at /pgp-ref.
MOUNTS=(-v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro")
ref_c() {
  case "$1" in
    "${GENOME_DIR}"/*) cpath "$1" ;;
    *) printf '/pgp-ref%s' "${1#"$REF_DIR"}" ;;
  esac
}
case "$REF_DIR/" in
  "${GENOME_DIR}"/*) ;;
  *) [ -d "$REF_DIR" ] && MOUNTS+=(-v "${REF_DIR}:/pgp-ref:ro") ;;
esac

TMP=$(mktemp -d "${S}/.lookup.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# slice VCF OUT REGION: the records of VCF that overlap REGION (reference
# blocks by their END), with the header. Without a .tbi, or when bcftools
# cannot read the region, the whole file is passed on instead: slower, the
# same answer, because the report skips records outside the region.
slice() {
  local vcf=$1 out=$2 region=$3
  if [ -s "${vcf}.tbi" ] || [ -s "${vcf}.csi" ]; then
    if run_in --cpus 1 --memory 1g "$BCFTOOLS_IMAGE" \
         bcftools view -r "$region" -o "$(cpath "$out")" "$(cpath "$vcf")" 2> "${TMP}/bcftools.err"; then
      cpath "$out"
      return 0
    fi
    echo "NOTE: bcftools could not read ${region} from ${vcf} ($(head -n 1 "${TMP}/bcftools.err")); reading the whole file." >&2
  else
    echo "NOTE: ${vcf} has no .tbi index; reading the whole file (bcftools index -t makes this fast)." >&2
  fi
  cpath "$vcf"
}

rc=0
for QUERY in "$@"; do
  RES=()
  [ -z "$SITES" ] || RES+=(--sites "$(ref_c "$SITES")")
  [ -z "$GENES" ] || RES+=(--genes "$(ref_c "$GENES")")
  [ -z "$ANN" ] || RES+=(--annotated "$(cpath "$ANN")")
  qrc=0
  run_in "${MOUNTS[@]}" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
    python3 /pgp-bin/genome_lookup.py resolve --query "$QUERY" \
      --table /pgp-assets/traits.tsv --table /pgp-assets/bloodtype_markers.tsv \
      ${RES[@]+"${RES[@]}"} > "${TMP}/regions.tsv" || qrc=$?
  if [ "$qrc" -ne 0 ]; then
    rc=$qrc
    echo ""
    continue
  fi
  while IFS=$'\t' read -r -u 3 kind chrom start end ref source; do
    region="${chrom}:${start}-${end}"
    G_IN=$(slice "$GVCF" "${TMP}/gvcf.vcf" "$region")
    REP=(--gvcf "$G_IN")
    if [ -n "$ANN" ]; then
      REP+=(--annotated "$(slice "$ANN" "${TMP}/annotated.vcf" "$region")")
    fi
    [ -z "$SUMMARY" ] || REP+=(--summary "$(cpath "$SUMMARY")")
    run_in "${MOUNTS[@]}" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
      python3 /pgp-bin/genome_lookup.py report --query "$QUERY" --kind "$kind" \
        --chrom "$chrom" --start "$start" --end "$end" --ref "$ref" --source "$source" \
        "${REP[@]}" || rc=$?
    echo ""
  done 3< "${TMP}/regions.tsv"
done
exit "$rc"
