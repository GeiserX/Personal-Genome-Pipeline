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
# VEP cache (~26 GB) and PCGR data (~21 GB) are downloaded separately
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
    wget -c -O "$FASTA" \
      "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" || {
      echo "ERROR: Failed to download reference genome."
      echo "  Try manually: wget -c -O ${FASTA} 'https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta'"
      exit 1
    }
  fi

  if [ ! -f "$FAI" ]; then
    wget -c -O "$FAI" \
      "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" || {
      echo "  Generating index with samtools..."
      docker run --rm --user root \
        -v "${GENOME_DIR}:/genome" \
        staphb/samtools:1.20 \
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
    wget -c -O "$CLINVAR" \
      "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
  fi
  if [ ! -f "$CLINVAR_TBI" ]; then
    wget -c -O "$CLINVAR_TBI" \
      "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
  fi

  # Step A: Chr-rename ClinVar (NCBI uses "1,2,3", pipeline uses "chr1,chr2,chr3")
  CLINVAR_CHR="${CLINVARDIR}/clinvar_chr.vcf.gz"
  if [ ! -f "$CLINVAR_CHR" ]; then
    echo "Creating chr-prefixed ClinVar..."
    docker run --rm --user root \
      -v "${GENOME_DIR}:/genome" \
      staphb/bcftools:1.21 \
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
      staphb/bcftools:1.21 \
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

IMAGES=(
  "staphb/samtools:1.20"
  "staphb/bcftools:1.21"
  "google/deepvariant:1.6.0"
  "quay.io/biocontainers/manta:1.6.0--h9ee0642_2"
  "getwilds/annotsv:latest"
  "pgkb/pharmcat:2.15.5"
  "quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0"
  "weisburd/expansionhunter:latest"
  "lgalarno/telomerehunter:latest"
  "genepi/haplogrep3:latest"
  "ensemblorg/ensembl-vep:release_112.0"
  "sigven/pcgr:1.4.1"
  "brentp/duphold:latest"
  "quay.io/biocontainers/goleft:0.2.4--h9ee0642_1"
  "quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11"
  "quay.io/biocontainers/delly:1.7.3--hd6466ae_0"
  "broadinstitute/gatk:4.6.1.0"
  "python:3.11-slim"
  "pgscatalog/plink2:2.00a5.10"
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
if [ -d "$PCGRDIR" ] && [ -d "${PCGRDIR}/data" ]; then
  echo "[OK] PCGR/CPSR data bundle already present."
else
  echo "[SKIP] PCGR/CPSR data (~21 GB) — needed for step 17 (cancer predisposition)"
  echo "  Download from: https://github.com/sigven/pcgr/releases"
  echo "    mkdir -p ${PCGRDIR}"
  echo "    wget -c -P ${PCGRDIR} https://github.com/sigven/pcgr/releases/download/v1.4.1/pcgr_data_grch38.tar.gz"
  echo "    tar xzf ${PCGRDIR}/pcgr_data_grch38.tar.gz -C ${PCGRDIR}"
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
