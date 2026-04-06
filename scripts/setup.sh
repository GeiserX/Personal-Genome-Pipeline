#!/usr/bin/env bash
# setup.sh — One-stop setup: download references, pull Docker images, validate
# Usage: ./scripts/setup.sh <genome_dir>
#
# This script downloads everything needed to run the pipeline:
#   1. GRCh38 reference genome + index (~3.5 GB)
#   2. ClinVar database (~200 MB)
#   3. All Docker images (~10-15 GB)
#   4. Validates the setup
#
# VEP cache (~26 GB) and PCGR ref data (~5 GB) are downloaded separately
# because they are only needed for specific steps and take a long time.
set -euo pipefail

GENOME_DIR=${1:-${GENOME_DIR:-""}}
if [ -z "$GENOME_DIR" ]; then
  echo "Usage: $0 <genome_dir>"
  echo ""
  echo "  <genome_dir>  Where to store reference data and sample outputs."
  echo "                Needs at least 500 GB free space per sample."
  echo ""
  echo "Example:"
  echo "  ./scripts/setup.sh /data/genomics"
  echo "  ./scripts/setup.sh ~/genome_data"
  exit 1
fi

export GENOME_DIR

# Source image versions from the canonical manifest
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

# Download helper with retry and resume
_download() {
  local url="$1" dest="$2" attempts="${3:-3}"
  for i in $(seq 1 "$attempts"); do
    if wget -c -O "$dest" "$url"; then
      return 0
    fi
    echo "  Attempt ${i}/${attempts} failed for $(basename "$dest")..."
    sleep $((i * 5))
  done
  echo "ERROR: Failed to download after ${attempts} attempts: $(basename "$dest")"
  return 1
}

echo "============================================"
echo "  Genomics Pipeline — Setup"
echo "  Data directory: ${GENOME_DIR}"
echo "============================================"
echo ""

# Check Docker
if ! command -v docker &>/dev/null; then
  echo "ERROR: Docker is not installed."
  echo "  Install Docker: https://docs.docker.com/get-docker/"
  exit 1
fi
if ! docker info &>/dev/null; then
  echo "ERROR: Docker daemon is not running."
  echo "  Start Docker Desktop or run: sudo systemctl start docker"
  exit 1
fi
echo "[OK] Docker is running."

# Check disk space
AVAIL_GB=$(df -BG "$GENOME_DIR" 2>/dev/null | awk 'NR==2 {gsub("G",""); print $4}' || echo "0")
if [ -z "$AVAIL_GB" ]; then
  # macOS df format
  AVAIL_KB=$(df -k "$GENOME_DIR" 2>/dev/null | awk 'NR==2 {print $4}' || echo "0")
  AVAIL_GB=$(( AVAIL_KB / 1048576 ))
fi
if [ "$AVAIL_GB" -lt 50 ]; then
  echo "WARNING: Only ${AVAIL_GB} GB free in ${GENOME_DIR}. Need at least 50 GB for references."
fi
echo "[OK] ${AVAIL_GB} GB free in ${GENOME_DIR}."

###############################################################################
# Phase 1: Reference Genome
###############################################################################
echo ""
echo "=== Phase 1: Reference Genome (~3.5 GB) ==="

REFDIR="${GENOME_DIR}/reference"
mkdir -p "$REFDIR"

FASTA="${REFDIR}/Homo_sapiens_assembly38.fasta"
FAI="${REFDIR}/Homo_sapiens_assembly38.fasta.fai"

if [ -f "$FASTA" ] && [ -f "$FAI" ]; then
  echo "[OK] Reference genome already downloaded."
else
  echo "Downloading GRCh38 reference genome..."
  echo "  Source: GATK resource bundle (Google Cloud)"
  echo "  Size: ~3.1 GB (FASTA) + ~2 MB (index)"

  if [ ! -f "$FASTA" ]; then
    _download "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" "$FASTA" || {
      echo "  Try manually: wget -c -O ${FASTA} 'https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta'"
      exit 1
    }
  fi

  if [ ! -f "$FAI" ]; then
    _download "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" "$FAI" || {
      echo "  Generating index with samtools..."
      docker run --rm --user root \
        -v "${GENOME_DIR}:/genome" \
        "$SAMTOOLS_IMAGE" \
        samtools faidx /genome/reference/Homo_sapiens_assembly38.fasta
    }
  fi
  echo "[OK] Reference genome downloaded."
fi

###############################################################################
# Phase 2: ClinVar Database
###############################################################################
echo ""
echo "=== Phase 2: ClinVar Database (~200 MB) ==="

CLINVARDIR="${GENOME_DIR}/clinvar"
mkdir -p "$CLINVARDIR"

CLINVAR="${CLINVARDIR}/clinvar.vcf.gz"
CLINVAR_TBI="${CLINVARDIR}/clinvar.vcf.gz.tbi"

if [ -f "$CLINVAR" ] && [ -f "$CLINVAR_TBI" ]; then
  echo "[OK] ClinVar database already downloaded."
else
  echo "Downloading ClinVar database..."
  if [ ! -f "$CLINVAR" ]; then
    _download "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" "$CLINVAR" || exit 1
  fi
  if [ ! -f "$CLINVAR_TBI" ]; then
    _download "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" "$CLINVAR_TBI" || exit 1
  fi

  # Step A: Chr-rename ClinVar (NCBI uses "1,2,3", pipeline uses "chr1,chr2,chr3")
  CLINVAR_CHR="${CLINVARDIR}/clinvar_chr.vcf.gz"
  if [ ! -f "$CLINVAR_CHR" ]; then
    echo "Creating chr-prefixed ClinVar..."
    docker run --rm --user root \
      -v "${GENOME_DIR}:/genome" \
      "$BCFTOOLS_IMAGE" \
      bash -c 'echo -e "1 chr1\n2 chr2\n3 chr3\n4 chr4\n5 chr5\n6 chr6\n7 chr7\n8 chr8\n9 chr9\n10 chr10\n11 chr11\n12 chr12\n13 chr13\n14 chr14\n15 chr15\n16 chr16\n17 chr17\n18 chr18\n19 chr19\n20 chr20\n21 chr21\n22 chr22\nX chrX\nY chrY\nMT chrM" > /genome/clinvar/chr_rename.txt &&
        bcftools annotate --rename-chrs /genome/clinvar/chr_rename.txt /genome/clinvar/clinvar.vcf.gz -Oz -o /genome/clinvar/clinvar_chr.vcf.gz &&
        bcftools index -t /genome/clinvar/clinvar_chr.vcf.gz'
  fi

  # Step B: Extract pathogenic + likely pathogenic subset from chr-renamed ClinVar
  CLINVAR_PATH="${CLINVARDIR}/clinvar_pathogenic_chr.vcf.gz"
  if [ ! -f "$CLINVAR_PATH" ]; then
    echo "Creating pathogenic/likely pathogenic subset..."
    docker run --rm --user root \
      -v "${GENOME_DIR}:/genome" \
      "$BCFTOOLS_IMAGE" \
      bash -c 'bcftools view -i "CLNSIG~\"Pathogenic\" || CLNSIG~\"Likely_pathogenic\"" /genome/clinvar/clinvar_chr.vcf.gz -Oz \
        -o /genome/clinvar/clinvar_pathogenic_chr.vcf.gz &&
        bcftools index -t /genome/clinvar/clinvar_pathogenic_chr.vcf.gz'
  fi
  echo "[OK] ClinVar database downloaded and indexed."
fi

###############################################################################
# Phase 3: Docker Images
###############################################################################
echo ""
echo "=== Phase 3: Docker Images (~10-15 GB total) ==="

# Image list sourced from versions.env
IMAGES=(
  "$MINIMAP2_IMAGE"
  "$SAMTOOLS_IMAGE"
  "$BCFTOOLS_IMAGE"
  "$DEEPVARIANT_IMAGE"
  "$MANTA_IMAGE"
  "$ANNOTSV_IMAGE"
  "$PHARMCAT_IMAGE"
  "$T1K_IMAGE"
  "$TELOMEREHUNTER_IMAGE"
  "$HAPLOGREP3_IMAGE"
  "$VEP_IMAGE"
  "$PCGR_IMAGE"
  "$DUPHOLD_IMAGE"
  "$GOLEFT_IMAGE"
  "$CNVNATOR_IMAGE"
  "$DELLY_IMAGE"
  "$GATK_IMAGE"
  "$PYTHON_IMAGE"
  "$PLINK2_IMAGE"
  "$FASTP_IMAGE"
  "$MOSDEPTH_IMAGE"
  "$MULTIQC_IMAGE"
  "$EXPANSIONHUNTER_IMAGE"
  "$GRIDSS_IMAGE"
  "$OCTOPUS_IMAGE"
  "$CLAIR3_IMAGE"
  "$SNIFFLES_IMAGE"
)

PULLED=0
SKIPPED=0
FAILED=0

for IMG in "${IMAGES[@]}"; do
  if docker image inspect "$IMG" &>/dev/null; then
    SKIPPED=$((SKIPPED + 1))
  else
    echo "  Pulling: ${IMG}..."
    if docker pull "$IMG" 2>/dev/null; then
      PULLED=$((PULLED + 1))
    else
      echo "  WARNING: Failed to pull ${IMG}. Check the image name/tag."
      FAILED=$((FAILED + 1))
    fi
  fi
done

echo "[OK] Docker images: ${PULLED} pulled, ${SKIPPED} already present, ${FAILED} failed."

if [ "$FAILED" -gt 0 ]; then
  echo ""
  echo "ERROR: ${FAILED} Docker image(s) failed to pull. Fix the issues above and re-run setup."
  exit 1
fi

###############################################################################
# Phase 4: Optional Downloads (instructions only)
###############################################################################
echo ""
echo "=== Phase 4: Optional Downloads (manual) ==="
echo ""
echo "The following are only needed for specific steps and are large downloads:"
echo ""

VEPDIR="${GENOME_DIR}/vep_cache"
if [ -d "$VEPDIR" ] && [ "$(find "$VEPDIR" -name "info.txt" 2>/dev/null | head -1)" ]; then
  echo "[OK] VEP cache already present."
else
  echo "[SKIP] VEP cache (~26 GB) — needed for step 13 (VEP annotation)"
  echo "  Download:"
  echo "    mkdir -p ${VEPDIR}"
  echo "    wget -c -P ${VEPDIR} https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz"
  echo "    tar xzf ${VEPDIR}/homo_sapiens_vep_112_GRCh38.tar.gz -C ${VEPDIR}"
fi
echo ""

PCGRDIR="${GENOME_DIR}/pcgr_data"
if [ -d "$PCGRDIR" ] && [ -d "${PCGRDIR}/20250314/data" ]; then
  echo "[OK] PCGR/CPSR ref data bundle already present."
else
  echo "[SKIP] PCGR/CPSR ref data (~5 GB) — needed for step 17 (cancer predisposition)"
  echo "  Download:"
  echo "    mkdir -p ${PCGRDIR} && cd ${PCGRDIR}"
  echo "    wget -c https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz"
  echo "    tar xzf pcgr_ref_data.20250314.grch38.tgz"
  echo "    mkdir -p 20250314 && mv data/ 20250314/"
fi
echo ""

# VEP 113 cache for CPSR (separate from step 13's release-112 cache)
VEP113="${VEPDIR}/homo_sapiens/113_GRCh38"
if [ -d "$VEP113" ]; then
  echo "[OK] VEP 113 cache (for CPSR) already present."
else
  echo "[SKIP] VEP 113 cache (~26 GB) — needed for step 17 (CPSR uses PCGR 2.2.5 which bundles VEP 113)"
  echo "  This is separate from the release-112 cache used by step 13. Both coexist in vep_cache/."
  echo "  Download:"
  echo "    mkdir -p ${VEPDIR}"
  echo "    wget -c -P ${VEPDIR} https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz"
  echo "    tar xzf ${VEPDIR}/homo_sapiens_vep_113_GRCh38.tar.gz -C ${VEPDIR}"
fi

###############################################################################
# Summary
###############################################################################
echo ""
echo "============================================"
echo "  Setup complete!"
echo ""
echo "  Reference genome:  ${REFDIR}/"
echo "  ClinVar database:  ${CLINVARDIR}/"
echo "  Docker images:     ${PULLED} pulled, ${SKIPPED} cached"
echo "============================================"
echo ""
echo "Next steps:"
echo "  1. Place your FASTQ/BAM/VCF in: ${GENOME_DIR}/<sample_name>/"
echo "  2. Run: export GENOME_DIR=${GENOME_DIR}"
echo "  3. Run: ./scripts/validate-setup.sh <sample_name>"
echo "  4. Run: ./scripts/run-all.sh <sample_name> <male|female>"
echo ""
echo "For optional VEP and CPSR setup, see the instructions above."
