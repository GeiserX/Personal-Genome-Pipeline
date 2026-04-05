#!/usr/bin/env bash
# validate-setup.sh — Verify all prerequisites before running the genomics pipeline
# Usage: ./scripts/validate-setup.sh [sample_name]
#
# Checks system requirements, reference data, Docker images, and (optionally)
# sample data readiness. Exits 0 if all critical checks pass, 1 otherwise.
set -euo pipefail

###############################################################################
# Color helpers (gracefully degrade if terminal does not support colors)
###############################################################################
if [ -t 1 ] && command -v tput &>/dev/null && [ "$(tput colors 2>/dev/null || echo 0)" -ge 8 ]; then
  GREEN=$(tput setaf 2)
  YELLOW=$(tput setaf 3)
  RED=$(tput setaf 1)
  BOLD=$(tput bold)
  RESET=$(tput sgr0)
else
  GREEN="" YELLOW="" RED="" BOLD="" RESET=""
fi

pass()  { echo "  ${GREEN}[OK]${RESET}    $1"; }
warn()  { echo "  ${YELLOW}[WARN]${RESET}  $1"; WARNINGS=$((WARNINGS + 1)); }
fail()  { echo "  ${RED}[FAIL]${RESET}  $1"; FAILURES=$((FAILURES + 1)); }
info()  { echo "  ${BOLD}[INFO]${RESET}  $1"; }
header(){ echo ""; echo "${BOLD}=== $1 ===${RESET}"; }

FAILURES=0
WARNINGS=0
MISSING_IMAGES=()

###############################################################################
# 1. System Requirements
###############################################################################
header "System Requirements"

# --- bash version ---
BASH_MAJOR="${BASH_VERSINFO[0]}"
BASH_MINOR="${BASH_VERSINFO[1]}"
if [ "$BASH_MAJOR" -ge 4 ]; then
  pass "bash ${BASH_MAJOR}.${BASH_MINOR} (>= 4.0 required)"
else
  fail "bash ${BASH_MAJOR}.${BASH_MINOR} — version 4.0+ is required. Install a newer bash."
fi

# --- Docker installed ---
if command -v docker &>/dev/null; then
  DOCKER_VERSION=$(docker --version 2>/dev/null | head -1)
  pass "Docker installed: ${DOCKER_VERSION}"
else
  fail "Docker is not installed. Install from https://docs.docker.com/get-docker/"
fi

# --- Docker daemon running ---
if command -v docker &>/dev/null && docker info &>/dev/null 2>&1; then
  pass "Docker daemon is running"

  # --- Docker memory ---
  DOCKER_MEM_BYTES=$(docker info --format '{{.MemTotal}}' 2>/dev/null || echo 0)
  if [ "$DOCKER_MEM_BYTES" -gt 0 ] 2>/dev/null; then
    DOCKER_MEM_GB=$(awk "BEGIN {printf \"%.1f\", ${DOCKER_MEM_BYTES}/1073741824}")
    if awk "BEGIN {exit !(${DOCKER_MEM_GB} >= 16)}" 2>/dev/null; then
      pass "Docker memory: ${DOCKER_MEM_GB} GB (>= 16 GB recommended)"
    elif awk "BEGIN {exit !(${DOCKER_MEM_GB} >= 8)}" 2>/dev/null; then
      warn "Docker memory: ${DOCKER_MEM_GB} GB — 16 GB+ recommended. Some steps may OOM."
    else
      fail "Docker memory: ${DOCKER_MEM_GB} GB — far too low. Increase to at least 16 GB."
    fi
  else
    warn "Could not detect Docker memory allocation"
  fi
else
  if command -v docker &>/dev/null; then
    fail "Docker daemon is not running. Start Docker Desktop or run: sudo systemctl start docker"
  fi
fi

# --- wget or curl ---
if command -v wget &>/dev/null; then
  pass "wget available (used for downloading reference data)"
elif command -v curl &>/dev/null; then
  pass "curl available (can substitute for wget in downloads)"
else
  fail "Neither wget nor curl found. Install one: apt install wget OR brew install wget"
fi

# --- CPU count ---
CPU_COUNT=0
if [ -f /proc/cpuinfo ]; then
  CPU_COUNT=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || echo 0)
elif command -v sysctl &>/dev/null; then
  CPU_COUNT=$(sysctl -n hw.ncpu 2>/dev/null || echo 0)
elif command -v nproc &>/dev/null; then
  CPU_COUNT=$(nproc 2>/dev/null || echo 0)
fi
if [ "$CPU_COUNT" -gt 0 ] 2>/dev/null; then
  if [ "$CPU_COUNT" -ge 16 ]; then
    pass "CPU cores: ${CPU_COUNT} (excellent for parallel steps)"
  elif [ "$CPU_COUNT" -ge 4 ]; then
    pass "CPU cores: ${CPU_COUNT} (adequate; 16+ recommended for faster runs)"
  else
    warn "CPU cores: ${CPU_COUNT} — minimum is 4. Pipeline will be very slow."
  fi
else
  info "Could not detect CPU count"
fi

# --- RAM ---
RAM_KB=0
if [ -f /proc/meminfo ]; then
  RAM_KB=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo 2>/dev/null || echo 0)
elif command -v sysctl &>/dev/null; then
  RAM_BYTES=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
  RAM_KB=$((RAM_BYTES / 1024))
fi
if [ "$RAM_KB" -gt 0 ] 2>/dev/null; then
  RAM_GB=$(awk "BEGIN {printf \"%.0f\", ${RAM_KB}/1048576}")
  if [ "$RAM_GB" -ge 32 ]; then
    pass "System RAM: ${RAM_GB} GB (can run multiple steps in parallel)"
  elif [ "$RAM_GB" -ge 16 ]; then
    pass "System RAM: ${RAM_GB} GB (adequate; 32 GB+ allows more parallelism)"
  else
    warn "System RAM: ${RAM_GB} GB — 16 GB minimum, 32 GB recommended. Reduce --memory flags in scripts."
  fi
else
  info "Could not detect system RAM"
fi

# --- Disk space ---
# Check free space in GENOME_DIR if set, otherwise CWD
CHECK_DIR="${GENOME_DIR:-$(pwd)}"
if command -v df &>/dev/null; then
  # Use 1K blocks for portability (works on Linux and macOS)
  FREE_KB=$(df -Pk "$CHECK_DIR" 2>/dev/null | awk 'NR==2 {print $4}')
  if [ -n "$FREE_KB" ] && [ "$FREE_KB" -gt 0 ] 2>/dev/null; then
    FREE_GB=$(awk "BEGIN {printf \"%.0f\", ${FREE_KB}/1048576}")
    if [ "$FREE_GB" -ge 500 ]; then
      pass "Free disk space: ${FREE_GB} GB in $(df -Pk "$CHECK_DIR" | awk 'NR==2 {print $6}')"
    elif [ "$FREE_GB" -ge 200 ]; then
      warn "Free disk space: ${FREE_GB} GB — 500 GB+ recommended for full pipeline per sample"
    else
      fail "Free disk space: ${FREE_GB} GB — critically low. Need 500 GB+ per sample."
    fi
  fi
fi

# --- Platform note ---
ARCH=$(uname -m 2>/dev/null || echo unknown)
OS=$(uname -s 2>/dev/null || echo unknown)
if [ "$ARCH" = "arm64" ] || [ "$ARCH" = "aarch64" ]; then
  warn "Architecture: ${ARCH} (${OS}) — all pipeline Docker images are amd64. Expect 2-5x slower due to emulation."
else
  info "Architecture: ${ARCH} (${OS})"
fi

###############################################################################
# 2. Environment Variables
###############################################################################
header "Environment Variables"

if [ -n "${GENOME_DIR:-}" ]; then
  pass "GENOME_DIR is set: ${GENOME_DIR}"
  if [ -d "$GENOME_DIR" ]; then
    pass "GENOME_DIR directory exists"
  else
    fail "GENOME_DIR directory does not exist: ${GENOME_DIR}"
    echo "       Create it with: mkdir -p \"${GENOME_DIR}\""
  fi
else
  fail "GENOME_DIR is not set. Export it before running the pipeline:"
  echo "       export GENOME_DIR=/path/to/your/data"
fi

###############################################################################
# 3. Reference Data
###############################################################################
header "Reference Data"

if [ -z "${GENOME_DIR:-}" ]; then
  info "Skipping reference data checks (GENOME_DIR not set)"
else
  # --- GRCh38 FASTA ---
  FASTA="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
  if [ -f "$FASTA" ]; then
    FASTA_SIZE=$(wc -c < "$FASTA" 2>/dev/null || echo 0)
    # 3 GB = 3221225472 bytes
    if [ "$FASTA_SIZE" -gt 3000000000 ] 2>/dev/null; then
      pass "GRCh38 FASTA: $(awk "BEGIN {printf \"%.1f\", ${FASTA_SIZE}/1073741824}") GB"
    else
      fail "GRCh38 FASTA exists but is too small ($(awk "BEGIN {printf \"%.1f\", ${FASTA_SIZE}/1073741824}") GB). Expected > 3 GB. Re-download it."
      echo "       See docs/00-reference-setup.md for download instructions."
    fi
  else
    fail "GRCh38 FASTA not found at: ${FASTA}"
    echo "       Download it:"
    echo "       mkdir -p ${GENOME_DIR}/reference"
    echo "       wget -P ${GENOME_DIR}/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
  fi

  # --- FASTA index (.fai) ---
  FAI="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta.fai"
  if [ -f "$FAI" ]; then
    pass "FASTA index (.fai) present"
  else
    fail "FASTA index not found at: ${FAI}"
    echo "       Download it:"
    echo "       wget -P ${GENOME_DIR}/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
  fi

  # --- ClinVar chr-prefixed VCF ---
  CLINVAR="${GENOME_DIR}/clinvar/clinvar_chr.vcf.gz"
  if [ -f "$CLINVAR" ]; then
    pass "ClinVar VCF (chr-prefixed): present"
  else
    # Check if raw ClinVar exists but chr-prefixed version is missing
    if [ -f "${GENOME_DIR}/clinvar/clinvar.vcf.gz" ]; then
      warn "ClinVar raw VCF found but chr-prefixed version missing."
      echo "       Run: ./scripts/setup.sh ${GENOME_DIR}  OR  see docs/00-reference-setup.md"
    else
      fail "ClinVar VCF not found at: ${CLINVAR}"
      echo "       Download and prepare it — see docs/00-reference-setup.md"
    fi
  fi

  # --- ClinVar pathogenic subset (required by step 6) ---
  CLINVAR_PATH="${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz"
  if [ -f "$CLINVAR_PATH" ]; then
    pass "ClinVar pathogenic subset: present"
  else
    if [ -f "$CLINVAR" ]; then
      warn "ClinVar chr-prefixed found but pathogenic subset missing."
      echo "       Run: ./scripts/setup.sh ${GENOME_DIR}  OR  see docs/00-reference-setup.md"
    fi
  fi

  # --- VEP cache (optional) ---
  VEP_DIR="${GENOME_DIR}/vep_cache/homo_sapiens"
  if [ -d "$VEP_DIR" ]; then
    pass "VEP cache directory: present"
  else
    warn "VEP cache not found at: ${VEP_DIR}"
    echo "       Required for step 13 (VEP annotation). Download ~26 GB:"
    echo "       wget -c https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz"
    echo "       See docs/00-reference-setup.md for full instructions."
  fi

  # --- PCGR data bundle (optional) ---
  PCGR_DIR="${GENOME_DIR}/pcgr_data/20250314/data"
  if [ -d "$PCGR_DIR" ]; then
    pass "PCGR/CPSR ref data bundle: present"
  else
    warn "PCGR ref data bundle not found at: ${PCGR_DIR}"
    echo "       Required for step 17 (CPSR cancer predisposition). Download ~5 GB:"
    echo "       cd ${GENOME_DIR}/pcgr_data"
    echo "       wget -c https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz"
    echo "       tar xzf pcgr_ref_data.20250314.grch38.tgz && mkdir -p 20250314 && mv data/ 20250314/"
    echo "       See docs/17-cpsr.md for full instructions."
  fi

  # --- GATK sequence dictionary (optional) ---
  DICT="${GENOME_DIR}/reference/Homo_sapiens_assembly38.dict"
  if [ -f "$DICT" ]; then
    pass "GATK sequence dictionary (.dict): present"
  else
    warn "GATK sequence dictionary not found at: ${DICT}"
    echo "       Some tools (GATK Mutect2/step 20) require it. Generate with:"
    echo "       docker run --rm -v \"\${GENOME_DIR}:/genome\" broadinstitute/gatk:4.6.1.0 \\"
    echo "         gatk CreateSequenceDictionary -R /genome/reference/Homo_sapiens_assembly38.fasta"
  fi
fi

###############################################################################
# 4. Docker Images
###############################################################################
header "Docker Images"

if ! command -v docker &>/dev/null || ! docker info &>/dev/null 2>&1; then
  info "Skipping Docker image checks (Docker not available)"
else
  # All images from docs/00-reference-setup.md
  IMAGES=(
    "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
    "staphb/samtools:1.20"
    "staphb/bcftools:1.21"
    "google/deepvariant:1.6.0"
    "quay.io/biocontainers/manta:1.6.0--h9ee0642_2"
    "quay.io/biocontainers/delly:1.7.3--hd6466ae_0"
    "quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11"
    "brentp/duphold:latest"
    "getwilds/annotsv:latest"
    "ensemblorg/ensembl-vep:release_112.0"
    "sigven/pcgr:2.2.5"
    "pgkb/pharmcat:3.2.0"
    "weisburd/expansionhunter:latest"
    "lgalarno/telomerehunter:latest"
    "genepi/haplogrep3:latest"
    "quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0"
    "quay.io/biocontainers/goleft:0.2.4--h9ee0642_1"
    "broadinstitute/gatk:4.6.1.0"
    "python:3.11-slim"
    "pgscatalog/plink2:2.00a5.10"
    "quay.io/biocontainers/fastp:1.3.1--h43da1c4_0"
    "quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0"
    "quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"
    "quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5"
    "quay.io/biocontainers/gridss:2.13.2--h96c455f_6"
    "dancooke/octopus:0.7.4"
    "hkubal/clair3:v2.0.0"
    "quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0"
  )

  PULLED=0
  for img in "${IMAGES[@]}"; do
    if docker image inspect "$img" &>/dev/null; then
      PULLED=$((PULLED + 1))
    else
      MISSING_IMAGES+=("$img")
    fi
  done

  if [ ${#MISSING_IMAGES[@]} -eq 0 ]; then
    pass "All ${#IMAGES[@]} Docker images are pulled"
  else
    pass "${PULLED}/${#IMAGES[@]} Docker images already pulled"
    warn "${#MISSING_IMAGES[@]} Docker image(s) not yet pulled:"
    for img in "${MISSING_IMAGES[@]}"; do
      echo "         docker pull ${img}"
    done
    echo ""
    echo "       Pull all missing images at once:"
    echo "         $(printf 'docker pull %s && ' "${MISSING_IMAGES[@]}" | sed 's/ && $//')  "
  fi
fi

###############################################################################
# 5. Sample Data (optional — only if $1 is provided)
###############################################################################
SAMPLE="${1:-}"

if [ -n "$SAMPLE" ]; then
  header "Sample Data: ${SAMPLE}"

  if [ -z "${GENOME_DIR:-}" ]; then
    info "Skipping sample checks (GENOME_DIR not set)"
  elif [ ! -d "${GENOME_DIR}" ]; then
    info "Skipping sample checks (GENOME_DIR does not exist)"
  else
    SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
    HAS_ORA=false
    HAS_FASTQ=false
    HAS_BAM=false
    HAS_VCF=false

    # Check for ORA files
    ORA_COUNT=0
    if [ -d "${SAMPLE_DIR}/fastq" ]; then
      ORA_COUNT=$(find "${SAMPLE_DIR}/fastq" -maxdepth 1 -name "*.ora" 2>/dev/null | wc -l | tr -d ' ')
    fi
    if [ "$ORA_COUNT" -gt 0 ]; then
      HAS_ORA=true
      pass "ORA files found: ${ORA_COUNT} file(s) in ${SAMPLE_DIR}/fastq/"
    fi

    # Check for FASTQ files
    R1="${SAMPLE_DIR}/fastq/${SAMPLE}_R1.fastq.gz"
    R2="${SAMPLE_DIR}/fastq/${SAMPLE}_R2.fastq.gz"
    if [ -f "$R1" ] && [ -f "$R2" ]; then
      HAS_FASTQ=true
      R1_SIZE=$(wc -c < "$R1" 2>/dev/null || echo 0)
      R1_GB=$(awk "BEGIN {printf \"%.1f\", ${R1_SIZE}/1073741824}")
      pass "FASTQ files found: R1 (${R1_GB} GB) + R2"
    elif [ -f "$R1" ]; then
      warn "Only R1 found at ${R1} — missing R2. Pipeline expects paired-end reads."
    else
      # Check for any FASTQ-like files with non-standard names
      FASTQ_COUNT=0
      if [ -d "${SAMPLE_DIR}/fastq" ]; then
        FASTQ_COUNT=$(find "${SAMPLE_DIR}/fastq" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) 2>/dev/null | wc -l | tr -d ' ')
      fi
      if [ "$FASTQ_COUNT" -gt 0 ]; then
        warn "Found ${FASTQ_COUNT} FASTQ file(s) in ${SAMPLE_DIR}/fastq/ but not named ${SAMPLE}_R1.fastq.gz / ${SAMPLE}_R2.fastq.gz"
        echo "       Rename or symlink them:"
        echo "         ln -s your_file_R1.fastq.gz ${SAMPLE_DIR}/fastq/${SAMPLE}_R1.fastq.gz"
        echo "         ln -s your_file_R2.fastq.gz ${SAMPLE_DIR}/fastq/${SAMPLE}_R2.fastq.gz"
      fi
    fi

    # Check for BAM
    BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
    BAI="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam.bai"
    if [ -f "$BAM" ]; then
      HAS_BAM=true
      BAM_SIZE=$(wc -c < "$BAM" 2>/dev/null || echo 0)
      BAM_GB=$(awk "BEGIN {printf \"%.1f\", ${BAM_SIZE}/1073741824}")
      pass "BAM file found: ${BAM_GB} GB"
      if [ -f "$BAI" ]; then
        pass "BAM index (.bai) present"
      else
        warn "BAM index not found. Create it before running BAM-dependent steps:"
        echo "       docker run --rm -v \"\${GENOME_DIR}:/genome\" staphb/samtools:1.20 \\"
        echo "         samtools index /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
      fi
    fi

    # Check for VCF
    VCF="${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz"
    TBI="${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz.tbi"
    if [ -f "$VCF" ]; then
      HAS_VCF=true
      VCF_SIZE=$(wc -c < "$VCF" 2>/dev/null || echo 0)
      VCF_MB=$(awk "BEGIN {printf \"%.0f\", ${VCF_SIZE}/1048576}")
      pass "VCF file found: ${VCF_MB} MB"
      if [ -f "$TBI" ]; then
        pass "VCF index (.tbi) present"
      else
        warn "VCF index not found. Create it before running VCF-dependent steps:"
        echo "       docker run --rm -v \"\${GENOME_DIR}:/genome\" staphb/bcftools:1.21 \\"
        echo "         bcftools index -t /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
      fi
    fi

    # Validate genome build (GRCh38) if BAM or VCF exists
    if $HAS_BAM && command -v docker >/dev/null 2>&1; then
      echo ""
      info "Checking genome build of BAM..."
      BAM_CHR1_LEN=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/samtools:1.20 \
        samtools view -H "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" 2>/dev/null | \
        grep "^@SQ" | grep "SN:chr1" | head -1 | sed 's/.*LN://' | cut -f1 || echo "0")
      if [ -z "$BAM_CHR1_LEN" ] || [ "$BAM_CHR1_LEN" = "0" ]; then
        # Try without chr prefix (hg19 style)
        BAM_CHR1_LEN=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/samtools:1.20 \
          samtools view -H "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" 2>/dev/null | \
          grep "^@SQ" | grep "SN:1[[:space:]]" | head -1 | sed 's/.*LN://' | cut -f1 || echo "0")
        if [ -n "$BAM_CHR1_LEN" ] && [ "$BAM_CHR1_LEN" != "0" ]; then
          fail "BAM uses chromosome names WITHOUT 'chr' prefix (hg19/GRCh37 style)"
          echo "       This pipeline requires GRCh38 (hg38) with 'chr' prefix."
          echo "       Extract FASTQ and re-align: samtools fastq -> step 2 (alignment)"
          echo "       See docs/vendor-guide.md for build conversion instructions."
        fi
      elif [ "$BAM_CHR1_LEN" = "248956422" ]; then
        pass "BAM genome build: GRCh38 (chr1 length = 248,956,422)"
      elif [ "$BAM_CHR1_LEN" = "249250621" ]; then
        fail "BAM genome build: GRCh37/hg19 (chr1 length = 249,250,621)"
        echo "       This pipeline requires GRCh38. Re-align from FASTQ."
      else
        warn "BAM chr1 length (${BAM_CHR1_LEN}) does not match known builds"
        echo "       Expected: 248956422 (GRCh38) or 249250621 (GRCh37)"
      fi
    fi

    if $HAS_VCF && command -v docker >/dev/null 2>&1; then
      VCF_CONTIG=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
        bcftools view -h "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null | \
        grep "^##contig=<ID=chr1," | head -1 || echo "")
      if [ -n "$VCF_CONTIG" ]; then
        VCF_CHR1_LEN=$(echo "$VCF_CONTIG" | sed 's/.*length=//' | tr -d '>' || echo "0")
        if [ "$VCF_CHR1_LEN" = "248956422" ]; then
          pass "VCF genome build: GRCh38"
        elif [ "$VCF_CHR1_LEN" = "249250621" ]; then
          fail "VCF genome build: GRCh37/hg19 — not compatible with this pipeline"
        fi
      else
        # Check for non-chr prefix
        VCF_NO_CHR=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
          bcftools view -h "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null | \
          grep "^##contig=<ID=1," | head -1 || echo "")
        if [ -n "$VCF_NO_CHR" ]; then
          fail "VCF uses chromosome names WITHOUT 'chr' prefix (likely GRCh37)"
          echo "       This pipeline requires GRCh38 with 'chr' prefix."
        fi
      fi
    fi

    # Suggest pipeline entry path
    echo ""
    if $HAS_ORA && ! $HAS_FASTQ && ! $HAS_BAM && ! $HAS_VCF; then
      info "Suggested: Path D (ORA -> FASTQ -> BAM -> VCF)"
      echo "       Start with:  ./scripts/01-ora-to-fastq.sh ${SAMPLE}"
    elif $HAS_FASTQ && ! $HAS_BAM && ! $HAS_VCF; then
      info "Suggested: Path A (FASTQ -> BAM -> VCF)"
      echo "       Start with:  ./scripts/02-alignment.sh ${SAMPLE}"
    elif $HAS_BAM && ! $HAS_VCF; then
      info "Suggested: Path B (BAM -> VCF)"
      echo "       Start with:  ./scripts/03-deepvariant.sh ${SAMPLE}"
    elif $HAS_VCF; then
      info "Suggested: Path C (VCF already available)"
      echo "       Start with:  ./scripts/06-clinvar-screen.sh ${SAMPLE}"
      if $HAS_BAM; then
        echo "       BAM also available — all 27 steps can run."
      else
        echo "       No BAM found — BAM-dependent steps (4, 10, 15, 16, 18, 19, 20) will be skipped."
      fi
    elif ! $HAS_ORA && ! $HAS_FASTQ && ! $HAS_BAM && ! $HAS_VCF; then
      if [ -d "$SAMPLE_DIR" ]; then
        fail "No usable data found for sample '${SAMPLE}' in ${SAMPLE_DIR}/"
        echo "       Expected one of:"
        echo "         ${SAMPLE_DIR}/fastq/${SAMPLE}_R1.fastq.gz + ${SAMPLE}_R2.fastq.gz  (Path A)"
        echo "         ${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam                         (Path B)"
        echo "         ${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz                                 (Path C)"
        echo "         ${SAMPLE_DIR}/fastq/*.ora                                          (Path D)"
      else
        fail "Sample directory does not exist: ${SAMPLE_DIR}/"
        echo "       Create it and place your data files inside:"
        echo "         mkdir -p ${SAMPLE_DIR}/fastq"
        echo "         # Copy your FASTQ/BAM/VCF files into the appropriate subdirectory"
      fi
    fi
  fi
fi

###############################################################################
# Summary
###############################################################################
header "Summary"

echo ""
if [ "$FAILURES" -eq 0 ] && [ "$WARNINGS" -eq 0 ]; then
  echo "  ${GREEN}${BOLD}All checks passed.${RESET} Your setup is ready to run the pipeline."
elif [ "$FAILURES" -eq 0 ]; then
  echo "  ${GREEN}${BOLD}All critical checks passed${RESET} with ${YELLOW}${WARNINGS} warning(s)${RESET}."
  echo "  The pipeline can run, but review the warnings above for best results."
else
  echo "  ${RED}${BOLD}${FAILURES} critical issue(s)${RESET} and ${YELLOW}${WARNINGS} warning(s)${RESET} found."
  echo "  Fix the ${RED}[FAIL]${RESET} items above before running the pipeline."
fi

if [ ${#MISSING_IMAGES[@]} -gt 0 ]; then
  echo ""
  echo "  To pull all missing Docker images (~10-15 GB total):"
  echo "    $(printf 'docker pull %s && ' "${MISSING_IMAGES[@]}" | sed 's/ && $//')"
fi

echo ""
echo "  Documentation:  docs/00-reference-setup.md   (download reference data)"
echo "                  docs/hardware-requirements.md (detailed requirements)"
echo "                  docs/vendor-guide.md          (data format help)"
echo ""

exit "$( [ "$FAILURES" -eq 0 ] && echo 0 || echo 1 )"
