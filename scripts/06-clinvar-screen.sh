#!/usr/bin/env bash
# ClinVar Pathogenic Screen — intersect sample VCF with ClinVar pathogenic/LP variants
# Finds known disease-causing variants the person carries
#
# Output: clinvar/${SAMPLE}_clinvar_hits.vcf — the sample's records that match a
# ClinVar Pathogenic/Likely_pathogenic allele, annotated with ClinVar's ID,
# GENEINFO, CLNSIG and CLNREVSTAT. Both reports (steps 24 and generate-report.sh)
# read this file.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF_DIR=${VCF_DIR:-vcf}
VCF="${GENOME_DIR}/${SAMPLE}/${VCF_DIR}/${SAMPLE}.vcf.gz"
REF="$REF_FASTA"
# Source file built by setup.sh; the normalised copy beside it is built from it.
CLINVAR="${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz"
CLINVAR_NORM="${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.norm.vcf.gz"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/clinvar"
HITS="${OUTPUT_DIR}/${SAMPLE}_clinvar_hits.vcf"

echo "=== ClinVar Pathogenic Screen: ${SAMPLE} ==="

for f in "$VCF" "${VCF}.tbi" "$CLINVAR" "${CLINVAR}.tbi" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# bcftools_run [--rw DIR] ARGS...: bcftools in its container.
bcftools_run() {
  local -a opt=()
  if [ "$1" = "--rw" ]; then opt=(--rw "$2"); shift 2; fi
  run_in ${opt[@]+"${opt[@]}"} --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" "$@"
}

# Step 0: the two files must share contig names. A ClinVar file with 1,2,...
# against a chr-prefixed VCF would otherwise give zero hits and look like a clean result.
SAMPLE_CONTIGS=$(bcftools_run bcftools index -s "$(cpath "$VCF")" | cut -f1 | sort -u)
CLINVAR_CONTIGS=$(bcftools_run bcftools index -s "$(cpath "$CLINVAR")" | cut -f1 | sort -u)
SHARED_CONTIGS=$(comm -12 <(printf '%s\n' "$SAMPLE_CONTIGS") <(printf '%s\n' "$CLINVAR_CONTIGS") | grep -c . || true)
if [ "$SHARED_CONTIGS" -eq 0 ]; then
  echo "ERROR: ${VCF} and ${CLINVAR} have no contig name in common." >&2
  echo "  Sample contigs:  $(printf '%s\n' "$SAMPLE_CONTIGS" | head -3 | tr '\n' ' ')..." >&2
  echo "  ClinVar contigs: $(printf '%s\n' "$CLINVAR_CONTIGS" | head -3 | tr '\n' ' ')..." >&2
  echo "  Rebuild the ClinVar file with chr-prefixed names (scripts/setup.sh does this)." >&2
  exit 1
fi

# Step 1: normalise ClinVar once (split multiallelics, left-align), rebuilt when the
# source is newer. isec matches on identical REF/ALT, so both sides must be split
# and left-aligned the same way.
if [ ! -f "$CLINVAR_NORM" ] || [ ! -f "${CLINVAR_NORM}.tbi" ] || [ "$CLINVAR" -nt "$CLINVAR_NORM" ]; then
  echo "Normalising ClinVar into ${CLINVAR_NORM} ..."
  rm -f "${CLINVAR_NORM}.part.vcf.gz" "${CLINVAR_NORM}.part.vcf.gz.tbi"
  # The normalised copy is shared by every sample, so clinvar/ is writable here.
  bcftools_run --rw "$(dirname "$CLINVAR_NORM")" bcftools norm -m -any -c w -f "$(cpath "$REF")" "$(cpath "$CLINVAR")" \
    -Oz -o "$(cpath "${CLINVAR_NORM}.part.vcf.gz")"
  bcftools_run --rw "$(dirname "$CLINVAR_NORM")" bcftools index -t "$(cpath "${CLINVAR_NORM}.part.vcf.gz")"
  mv "${CLINVAR_NORM}.part.vcf.gz" "$CLINVAR_NORM"
  mv "${CLINVAR_NORM}.part.vcf.gz.tbi" "${CLINVAR_NORM}.tbi"
fi

# Step 2: filter the sample. Callers that never write PASS (FILTER '.') would lose
# every record under -f PASS, so fall back to '.,PASS' only when no record is PASS.
HAS_PASS=$(run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  sh -c "bcftools view -H -f PASS '$(cpath "$VCF")' | head -n 1 | wc -l")
if [ "$HAS_PASS" -gt 0 ]; then
  FILTER="PASS"
  echo "Filter mode: PASS records only"
else
  FILTER=".,PASS"
  echo "NOTICE: ${VCF} has no PASS record; using records with FILTER '.' or PASS"
fi

PASS_VCF="${OUTPUT_DIR}/${SAMPLE}_pass.vcf.gz"
run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  sh -c "set -e; bcftools view -f '${FILTER}' -Ou '$(cpath "$VCF")' \
    | bcftools norm -m -any -c w -f '$(cpath "$REF")' -Oz -o '$(cpath "$PASS_VCF")' -
    bcftools index -f -t '$(cpath "$PASS_VCF")'"

PASS_COUNT=$(bcftools_run bcftools index -n "$(cpath "$PASS_VCF")")
if [ "$PASS_COUNT" -eq 0 ]; then
  echo "ERROR: no records left in ${VCF} after the '${FILTER}' filter; nothing to screen." >&2
  exit 1
fi

# Step 3: keep the sample's records whose allele is in ClinVar, then copy ClinVar's
# ID, gene and significance onto them. isec -w1 writes the sample's side, which has
# no GENEINFO/CLNSIG of its own. (annotate -a needs an indexed target, so the
# shared records go to a file first.)
SHARED="${OUTPUT_DIR}/${SAMPLE}_shared.vcf.gz"
run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  sh -c "set -e
    bcftools isec -n=2 -w1 -Oz -o '$(cpath "$SHARED")' '$(cpath "$PASS_VCF")' '$(cpath "$CLINVAR_NORM")'
    bcftools index -f -t '$(cpath "$SHARED")'
    bcftools annotate -a '$(cpath "$CLINVAR_NORM")' --pair-logic exact \
      -c ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT -Ov -o '$(cpath "${HITS}.part")' '$(cpath "$SHARED")'"
mv "${HITS}.part" "$HITS"
rm -f "$SHARED" "${SHARED}.tbi"

HIT_COUNT=$(grep -c -v '^#' "$HITS" || true)
echo "=== ClinVar screen complete ==="
echo "Hits: ${HITS} (sample records matching a ClinVar Pathogenic/Likely_pathogenic allele)"
echo "Count: ${HIT_COUNT} pathogenic hits"
