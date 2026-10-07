#!/usr/bin/env bash
# Imputation Prep — Split VCF by chromosome for Michigan Imputation Server upload
# Creates one PASS (and unfiltered) bgzipped, tabix-indexed VCF per chromosome,
# chr1-22 and chrX, in one container and one pass per chromosome.
# NOTE: MIS requires 20+ samples per job. Single WGS = mainly useful for phasing.
# Input: the variant-only VCF, where a site that matches the reference is
# absent, not a confirmed 0/0. With IMPUTATION_SITES (a CHROM<TAB>POS file of
# the reference panel's sites, or a VCF of them, inside GENOME_DIR) and step
# 03's gVCF, the input is instead every panel site genotyped from the gVCF:
# variant calls and 0/0 calls, with uncovered sites left out. A panel of tens
# of millions of sites needs a few GB of memory here.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
GVCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
IMPUTATION_SITES=${IMPUTATION_SITES:-}
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/imputation"
MIS_DIR="${OUTPUT_DIR}/mis_ready"

echo "=== Imputation Prep: ${SAMPLE} ==="

for f in "$VCF" "${VCF}.tbi"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$MIS_DIR"

INPUT_C="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
if [ -n "$IMPUTATION_SITES" ]; then
  for f in "$IMPUTATION_SITES" "$GVCF" "${GVCF}.tbi"; do
    if [ ! -f "$f" ]; then
      echo "ERROR: IMPUTATION_SITES needs ${f}; step 03 writes the gVCF." >&2
      exit 1
    fi
  done
  SITES_C=$(cpath "$IMPUTATION_SITES") || exit 2
  echo "Input: the sites in ${IMPUTATION_SITES}, genotyped from ${GVCF}"
  # gvcf2vcf expands the reference blocks over the panel sites (one 0/0 record
  # per base, the base from the reference); -T keeps the panel sites,
  # --trim-alt-alleles drops the <*> allele, and a no-call (./.) is dropped.
  # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
  run_in --cpus 2 --memory 8g \
    "${BCFTOOLS_IMAGE}" \
    bash -euo pipefail -c '
      bcftools convert --gvcf2vcf -f "$1" -R "$2" -Ou "$3" \
        | bcftools view -T "$2" --trim-alt-alleles -i "GT!=\"mis\"" -Oz -o "$4.part" --write-index=tbi
      mv -f "$4.part" "$4"
      mv -f "$4.part.tbi" "$4.tbi"' \
    _ "${REF_FASTA_C}" "$SITES_C" "/genome/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" \
      "/genome/${SAMPLE}/imputation/${SAMPLE}.panel_sites.vcf.gz"
  INPUT_C="/genome/${SAMPLE}/imputation/${SAMPLE}.panel_sites.vcf.gz"
elif [ -f "$GVCF" ]; then
  echo "Input: ${VCF} (variant sites only). Set IMPUTATION_SITES to the panel's site list"
  echo "  to add the sites where you match the reference, genotyped from ${GVCF}."
fi

# One container for all chromosomes; each file is written under a .part name
# with its index (--write-index) and renamed when both are complete: the old
# index goes first, then the VCF, then its index, so an index never sits
# beside a VCF it was not built from. A chromosome the VCF's index lists no
# record on gets no file (a file left by an earlier run on another input is
# removed), and the log names it (bcftools would stop on a region it cannot
# place).
CHROMS=()
for i in $(seq 1 22) X; do CHROMS+=("chr${i}"); done
run_in --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bash -euo pipefail -c '
    in=$1 out=$2 sample=$3; shift 3
    present=" $(bcftools index -s "$in" | cut -f1 | tr "\n" " ") "
    for chr in "$@"; do
      case "$present" in
        *" ${chr} "*) ;;
        *)
          echo "No records on ${chr}: no file."
          rm -f "${out}/${sample}_${chr}.vcf.gz" "${out}/${sample}_${chr}.vcf.gz.tbi"
          continue ;;
      esac
      echo "MIS-ready ${chr}..."
      part="${out}/${sample}_${chr}.part.vcf.gz"
      final="${out}/${sample}_${chr}.vcf.gz"
      bcftools view -f PASS,. -r "$chr" -Oz --write-index=tbi -o "$part" "$in"
      rm -f "${final}.tbi"
      mv -f "$part" "$final"
      mv -f "${part}.tbi" "${final}.tbi"
    done' \
  _ "$INPUT_C" "/genome/${SAMPLE}/imputation/mis_ready" "$SAMPLE" "${CHROMS[@]}"
rm -f "${OUTPUT_DIR}/${SAMPLE}.panel_sites.vcf.gz" "${OUTPUT_DIR}/${SAMPLE}.panel_sites.vcf.gz.tbi"

echo "=== Imputation prep complete ==="
echo "MIS-ready VCFs: ${MIS_DIR}/${SAMPLE}_chr{1-22,X}.vcf.gz"
echo ""
echo "Next steps:"
echo "  1. Register at https://imputationserver.sph.umich.edu"
echo "  2. Upload the 23 chromosome VCFs (chr1-22 and chrX)"
echo "  3. Select: refpanel=TOPMed r2, population=eur, build=hg38"
echo "  4. Results expire after 7 days — download immediately"
