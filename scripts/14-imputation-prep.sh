#!/usr/bin/env bash
# Imputation Prep — Split VCF by chromosome for Michigan Imputation Server upload
# Creates one PASS (and unfiltered) bgzipped, tabix-indexed VCF per chromosome,
# chr1-22 and chrX, in one container and one pass per chromosome.
# NOTE: MIS requires 20+ samples per job. Single WGS = mainly useful for phasing.
# NOTE: the input is the variant-only VCF, so a site where the sample matches
# the reference is absent, not a confirmed 0/0, until a gVCF-based input exists.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
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

# One container for all chromosomes; each file is written under a .part name
# with its index (--write-index) and renamed when both are complete. A
# chromosome the VCF's index lists no record on gets no file, and the log
# names it (bcftools would stop on a region it cannot place).
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
        *) echo "No records on ${chr}: no file."; continue ;;
      esac
      echo "MIS-ready ${chr}..."
      part="${out}/${sample}_${chr}.part.vcf.gz"
      bcftools view -f PASS,. -r "$chr" -Oz --write-index=tbi -o "$part" "$in"
      mv -f "${part}.tbi" "${out}/${sample}_${chr}.vcf.gz.tbi"
      mv -f "$part" "${out}/${sample}_${chr}.vcf.gz"
    done' \
  _ "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "/genome/${SAMPLE}/imputation/mis_ready" "$SAMPLE" "${CHROMS[@]}"

echo "=== Imputation prep complete ==="
echo "MIS-ready VCFs: ${MIS_DIR}/${SAMPLE}_chr{1-22,X}.vcf.gz"
echo ""
echo "Next steps:"
echo "  1. Register at https://imputationserver.sph.umich.edu"
echo "  2. Upload the 23 chromosome VCFs (chr1-22 and chrX)"
echo "  3. Select: refpanel=TOPMed r2, population=eur, build=hg38"
echo "  4. Results expire after 7 days — download immediately"
