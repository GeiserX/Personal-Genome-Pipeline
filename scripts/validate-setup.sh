#!/usr/bin/env bash
# validate-setup.sh — Verify all prerequisites before running the Personal Genome Pipeline
# Usage: ./scripts/validate-setup.sh [sample_name]
#
# Checks system requirements, reference data, Docker images, and (optionally)
# sample data readiness. Exits 0 if all critical checks pass, 1 otherwise.
set -euo pipefail

# Image versions, data versions and the docker wrapper. They are needed by the
# image list and by the sample checks, which also run when the Docker daemon
# is not up.
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
SAMPLE="${1:-}"
if [ -n "$SAMPLE" ]; then validate_sample "$SAMPLE"; fi

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

# check_bam_reference HEADER: the sample BAM was aligned to REF_FASTA. Its @SQ
# names and lengths must be the .fai's, in the same order (callers look reads
# up by contig index), and it must be coordinate-sorted. Tools give wrong
# results on a BAM from another reference long before any of them fails, so
# this stops the run here. A header with no @SQ line at all is left to
# check_bam_quickcheck, which fails a BAM without one.
check_bam_reference() {
  local header=$1 fai="${REF_FASTA}.fai" diff so bam
  if [ -f "$fai" ] && grep -q '^@SQ' <<< "$header"; then
    diff=$(printf '%s\n' "$header" | awk -F'\t' '
      FNR == 1 { file++ }
      file == 1 && /^@SQ/ {
        sn = ""; ln = ""
        for (i = 2; i <= NF; i++) {
          if ($i ~ /^SN:/) sn = substr($i, 4)
          else if ($i ~ /^LN:/) ln = substr($i, 4)
        }
        bam[++nb] = sn " (" ln " bp)"
      }
      file == 2 && NF >= 2 { ref[++nr] = $1 " (" $2 " bp)" }
      END {
        n = nb > nr ? nb : nr
        for (i = 1; i <= n; i++) {
          if (bam[i] != ref[i]) {
            printf "sequence %d is %s in the BAM and %s in the reference", i,
              (i <= nb ? bam[i] : "missing"), (i <= nr ? ref[i] : "missing")
            exit
          }
        }
      }' - "$fai")
    if [ -z "$diff" ]; then
      pass "BAM header matches the reference: the same $(grep -c . "$fai") sequences in the same order"
    else
      fail "this BAM was aligned to a different reference: realign (docs/realignment.md)"
      echo "       First difference: ${diff}."
      echo "       Reference: ${REF_FASTA}"
    fi
  fi
  so=$(awk -F'\t' '/^@HD/ { for (i = 2; i <= NF; i++) if ($i ~ /^SO:/) print substr($i, 4) }' <<< "$header")
  case "$so" in
    coordinate) pass "BAM is coordinate-sorted (@HD SO:coordinate)" ;;
    "")
      # samtools index refuses an unsorted BAM, so an index made from this
      # BAM shows the order the header does not state. One older than the BAM
      # may belong to a BAM it replaced (htslib warns about that too).
      bam="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
      if [ -f "${bam}.bai" ] && [ ! "${bam}.bai" -ot "$bam" ]; then
        pass "BAM is coordinate-sorted (no @HD SO: tag, but samtools only indexes a sorted BAM and its .bai is not older than it)"
      elif [ -f "${bam}.bai" ]; then
        fail "BAM header does not say how it is sorted (no @HD SO: tag) and its .bai is older than the BAM: index it again (samtools index fails on an unsorted BAM) or realign with step 02"
      else
        fail "BAM header does not say how it is sorted (no @HD SO: tag) and it has no .bai: sort and index it (samtools sort, samtools index) or realign with step 02"
      fi ;;
    *) fail "BAM is sorted by ${so}, not by coordinate: sort it (samtools sort) or realign with step 02" ;;
  esac
}

# check_bam_quickcheck: the sample BAM passes samtools quickcheck (a header
# with at least one sequence, and the end-of-file block). Runs whether or not
# its header could be read.
check_bam_quickcheck() {
  if run_in "${SAMTOOLS_IMAGE}" samtools quickcheck "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" >/dev/null 2>&1; then
    pass "BAM passes samtools quickcheck"
  else
    fail "BAM fails samtools quickcheck (truncated, no sequences in its header, or not a BAM): realign with step 02"
  fi
}

###############################################################################
# 1. System Requirements
###############################################################################
header "System Requirements"

# --- bash version ---
# 4.4: run-all.sh refuses anything older (it expands empty arrays under set -u).
BASH_MAJOR="${BASH_VERSINFO[0]}"
BASH_MINOR="${BASH_VERSINFO[1]}"
if [ $((BASH_MAJOR * 100 + BASH_MINOR)) -ge 404 ]; then
  pass "bash ${BASH_MAJOR}.${BASH_MINOR} (>= 4.4 required)"
else
  fail "bash ${BASH_MAJOR}.${BASH_MINOR} — version 4.4+ is required. Install a newer bash (on macOS: brew install bash)."
fi

# --- Docker installed ---
if command -v "$CONTAINER_ENGINE" &>/dev/null; then
  DOCKER_VERSION=$("$CONTAINER_ENGINE" --version 2>/dev/null | head -1)
  pass "Docker installed: ${DOCKER_VERSION}"
else
  fail "Docker is not installed. Install from https://docs.docker.com/get-docker/"
fi

# --- Docker daemon running ---
if command -v "$CONTAINER_ENGINE" &>/dev/null && "$CONTAINER_ENGINE" info &>/dev/null 2>&1; then
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
  if command -v "$CONTAINER_ENGINE" &>/dev/null; then
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
  warn "Architecture: ${ARCH} (${OS}) — most pipeline Docker images are amd64 only. Expect those to run 2-5x slower under emulation."
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
  FASTA="$REF_FASTA"
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
    echo "       Download it: ./scripts/setup.sh ${GENOME_DIR}"
  fi

  # --- FASTA index (.fai) ---
  FAI="${REF_FASTA}.fai"
  if [ -f "$FAI" ]; then
    pass "FASTA index (.fai) present"
  else
    fail "FASTA index not found at: ${FAI}"
    echo "       Download it: ./scripts/setup.sh ${GENOME_DIR}"
  fi

  # --- No ALT or HLA contigs ---
  # No aligner here runs ALT-aware: a read that matches the primary assembly
  # and an ALT or HLA contig equally gets MAPQ 0, and callers drop it. Depth
  # then thins at CYP2D6, the MHC, KIR and other loci those contigs copy.
  if [ -f "$FAI" ]; then
    ALT_CONTIGS=$(awk -F'\t' '$1 ~ /_alt$/ || $1 ~ /^HLA-/' "$FAI" | wc -l | tr -d ' ')
    if [ "$ALT_CONTIGS" -eq 0 ]; then
      pass "Reference has no ALT or HLA contigs ($(grep -c . "$FAI") sequences)"
    elif [ "${ALLOW_ALT_REFERENCE:-false}" = "true" ]; then
      warn "Reference has ${ALT_CONTIGS} ALT/HLA contigs (allowed by ALLOW_ALT_REFERENCE=true)."
      echo "       Reads that match a primary locus and its ALT copy equally get MAPQ 0, so depth"
      echo "       and calls thin at CYP2D6, the MHC and KIR. See docs/realignment.md."
    else
      fail "Reference has ${ALT_CONTIGS} ALT/HLA contigs (first: $(awk -F'\t' '$1 ~ /_alt$/ || $1 ~ /^HLA-/ {print $1; exit}' "$FAI"))."
      echo "       No step aligns ALT-aware, so reads that match a primary locus and its ALT copy"
      echo "       equally get MAPQ 0 and callers ignore them: depth and calls thin at CYP2D6, the MHC"
      echo "       and KIR. Use the default no-ALT reference (./scripts/setup.sh ${GENOME_DIR}),"
      echo "       or set ALLOW_ALT_REFERENCE=true to keep this one on purpose. See docs/realignment.md."
    fi
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
      fail "ClinVar chr-prefixed found but pathogenic subset missing (step 6 will fail)."
      echo "       Run: ./scripts/setup.sh ${GENOME_DIR}  OR  see docs/00-reference-setup.md"
    else
      fail "ClinVar pathogenic subset not found (step 6 will fail)."
      echo "       Download and prepare it — see docs/00-reference-setup.md"
    fi
  fi

  # --- ClinVar release date (clinvar/RELEASE, written by setup.sh) ---
  # NCBI publishes ClinVar monthly; a copy over 35 days old misses the last release.
  if [ -f "${GENOME_DIR}/clinvar/clinvar.vcf.gz" ]; then
    CLINVAR_RELEASE=$(head -n 1 "${GENOME_DIR}/clinvar/RELEASE" 2>/dev/null || true)
    if [[ "$CLINVAR_RELEASE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
      RELEASE_EPOCH=$(date -u -d "$CLINVAR_RELEASE" +%s 2>/dev/null \
        || date -u -j -f %Y-%m-%d "$CLINVAR_RELEASE" +%s 2>/dev/null || echo "")
      if [ -n "$RELEASE_EPOCH" ]; then
        CLINVAR_AGE=$(( ($(date -u +%s) - RELEASE_EPOCH) / 86400 ))
        if [ "$CLINVAR_AGE" -gt 35 ]; then
          warn "ClinVar release ${CLINVAR_RELEASE} is ${CLINVAR_AGE} days old (over 35)."
          echo "       Replace it with the current release: ./scripts/setup.sh --refresh clinvar ${GENOME_DIR}"
        else
          pass "ClinVar release ${CLINVAR_RELEASE} (${CLINVAR_AGE} days old)"
        fi
      else
        warn "ClinVar release ${CLINVAR_RELEASE}: could not work out its age on this system"
      fi
    else
      warn "ClinVar release date unknown (no ${GENOME_DIR}/clinvar/RELEASE)."
      echo "       Record it with the current release: ./scripts/setup.sh --refresh clinvar ${GENOME_DIR}"
    fi
  fi

  # --- Small pinned data files (setup.sh installs them) ---
  for name in $DATA_FILES; do
    case "$name" in
      delly_exclude) what="Delly exclude map (step 19 runs without it, slower and with more artefacts)" ;;
      cytoband) what="GRCh38 chromosome bands (step 10 falls back to hg19 bands)" ;;
      hla_dat) what="IPD-IMGT/HLA ${HLA_DB_RELEASE} (step 08 is skipped without it)" ;;
      gencode_genes) what="GENCODE ${GENCODE_RELEASE} gene coordinates (step 08 is skipped without them)" ;;
    esac
    if DATA_PATH=$(data_file "$name"); then
      pass "${what%% (*}: present"
    else
      warn "${what} not found at: ${DATA_PATH}"
      echo "       Install it: ./scripts/setup.sh ${GENOME_DIR}"
    fi
  done

  # --- somalier's sites and VerifyBamID2's panel, for step 33 (setup.sh installs them) ---
  QC_MISSING=""
  for f in somalier/sites.hg38.vcf.gz verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.UD \
           verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.mu verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.bed; do
    [ -s "${GENOME_DIR}/reference/${f}" ] || QC_MISSING="${QC_MISSING} reference/${f}"
  done
  if [ -z "$QC_MISSING" ]; then
    pass "somalier sites and VerifyBamID2 panel (step 33): present"
  else
    warn "Step 33 (sample identity and contamination) data not found:${QC_MISSING}"
    echo "       Install it: ./scripts/setup.sh --sample-qc-data ${GENOME_DIR}"
  fi

  # --- VEP cache of the VEP image's release, for step 13 (optional) ---
  VEP_DIR="${GENOME_DIR}/vep_cache/homo_sapiens"
  if [ -f "${VEP_DIR}/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
    pass "VEP ${VEP_CACHE_RELEASE} cache (step 13): present"
  else
    warn "VEP ${VEP_CACHE_RELEASE} cache not found at: ${VEP_DIR}/${VEP_CACHE_RELEASE}_GRCh38"
    echo "       Required for step 13 (VEP annotation), which downloads it (~26 GB) the first time it runs."
    echo "       See docs/00-reference-setup.md for full instructions."
  fi

  # --- VEP cache of the release inside the PCGR image, for step 17/CPSR (optional) ---
  if [ -f "${VEP_DIR}/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
    pass "VEP ${PCGR_VEP_CACHE_RELEASE} cache (step 17 CPSR): present"
  else
    warn "VEP ${PCGR_VEP_CACHE_RELEASE} cache not found at: ${VEP_DIR}/${PCGR_VEP_CACHE_RELEASE}_GRCh38"
    echo "       Required for step 17 (CPSR). ${PCGR_IMAGE} needs VEP ${PCGR_VEP_CACHE_RELEASE}, separate from step 13's VEP ${VEP_CACHE_RELEASE}."
    echo "       curl -fL -C - -O $(vep_cache_url "$PCGR_VEP_CACHE_RELEASE")"
    echo "       See docs/17-cpsr.md for full instructions."
  fi

  # --- PCGR data bundle (optional) ---
  PCGR_DIR="${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}/data"
  if [ -d "$PCGR_DIR" ]; then
    pass "PCGR/CPSR ref data bundle: present"
  else
    warn "PCGR ref data bundle not found at: ${PCGR_DIR}"
    echo "       Required for step 17 (CPSR cancer predisposition). Download ~7 GB:"
    echo "       cd ${GENOME_DIR}/pcgr_data"
    echo "       curl -fL -C - -O https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.${PCGR_DATA_BUNDLE}.grch38.tgz"
    echo "       tar xzf pcgr_ref_data.${PCGR_DATA_BUNDLE}.grch38.tgz && mkdir -p ${PCGR_DATA_BUNDLE} && mv data/ ${PCGR_DATA_BUNDLE}/"
    echo "       See docs/17-cpsr.md for full instructions."
  fi

  # --- AnnotSV annotation data (step 5) ---
  ANNOTSV_DIR="${GENOME_DIR}/annotsv_annotations"
  if [ -d "${ANNOTSV_DIR}/Annotations_Human/Genes/GRCh38" ]; then
    pass "AnnotSV annotation data (step 5): present"
  else
    warn "AnnotSV annotation data not found at: ${ANNOTSV_DIR}/Annotations_Human/Genes/GRCh38"
    echo "       Step 5 (AnnotSV) will be skipped. Download it (~5.3 GB, ~20 GB unpacked) with: ./scripts/setup.sh ${GENOME_DIR}"
  fi

  # --- Annotation databases (optional, for steps 30-31) ---
  ANNOT_DIR="${GENOME_DIR}/annotations"
  ANNOT_COUNT=0
  ANNOT_TOTAL=7
  ANNOT_MISSING=()
  if [ -f "${ANNOT_DIR}/whole_genome_SNVs.tsv.gz" ] && [ -f "${ANNOT_DIR}/whole_genome_SNVs.tsv.gz.tbi" ]; then
    ANNOT_COUNT=$((ANNOT_COUNT + 1))
  else
    ANNOT_MISSING+=("CADD SNVs (whole_genome_SNVs.tsv.gz + .tbi)")
  fi
  if [ -f "${ANNOT_DIR}/gnomad.genomes.r4.0.indel.tsv.gz" ] && [ -f "${ANNOT_DIR}/gnomad.genomes.r4.0.indel.tsv.gz.tbi" ]; then
    ANNOT_COUNT=$((ANNOT_COUNT + 1))
  else
    ANNOT_MISSING+=("CADD indels (gnomad.genomes.r4.0.indel.tsv.gz + .tbi)")
  fi
  # SpliceAI: step 30 accepts the raw or the masked score files
  for kind in snv indel; do
    if { [ -f "${ANNOT_DIR}/spliceai_scores.raw.${kind}.hg38.vcf.gz" ] && [ -f "${ANNOT_DIR}/spliceai_scores.raw.${kind}.hg38.vcf.gz.tbi" ]; } || \
       { [ -f "${ANNOT_DIR}/spliceai_scores.masked.${kind}.hg38.vcf.gz" ] && [ -f "${ANNOT_DIR}/spliceai_scores.masked.${kind}.hg38.vcf.gz.tbi" ]; }; then
      ANNOT_COUNT=$((ANNOT_COUNT + 1))
    else
      ANNOT_MISSING+=("SpliceAI ${kind}s (spliceai_scores.raw.${kind}.hg38.vcf.gz or spliceai_scores.masked.${kind}.hg38.vcf.gz, + .tbi)")
    fi
  done
  if [ -f "${ANNOT_DIR}/revel_grch38.tsv.gz" ] && [ -f "${ANNOT_DIR}/revel_grch38.tsv.gz.tbi" ]; then
    ANNOT_COUNT=$((ANNOT_COUNT + 1))
  else
    ANNOT_MISSING+=("REVEL (revel_grch38.tsv.gz + .tbi)")
  fi
  if [ -f "${ANNOT_DIR}/AlphaMissense_hg38.tsv.gz" ] && [ -f "${ANNOT_DIR}/AlphaMissense_hg38.tsv.gz.tbi" ]; then
    ANNOT_COUNT=$((ANNOT_COUNT + 1))
  else
    ANNOT_MISSING+=("AlphaMissense (AlphaMissense_hg38.tsv.gz + .tbi)")
  fi
  [ -f "${ANNOT_DIR}/gnomad_v4.1_constraint.tsv" ] && ANNOT_COUNT=$((ANNOT_COUNT + 1)) || ANNOT_MISSING+=("gnomAD constraint (gnomad_v4.1_constraint.tsv)")
  if [ "$ANNOT_COUNT" -eq "$ANNOT_TOTAL" ]; then
    pass "Annotation databases: all ${ANNOT_TOTAL} present (CADD SNVs+indels, SpliceAI SNVs+indels, REVEL, AlphaMissense, gnomAD constraint)"
  elif [ "$ANNOT_COUNT" -gt 0 ]; then
    warn "Annotation databases: ${ANNOT_COUNT}/${ANNOT_TOTAL} present. Steps 30-31 will degrade gracefully for missing tracks."
    for missing in "${ANNOT_MISSING[@]}"; do
      echo "         Missing: ${missing}"
    done
    echo "       See docs/00-reference-setup.md for download instructions."
  else
    info "Annotation databases not downloaded (steps 30-31 will be skipped). See docs/00-reference-setup.md"
  fi

  # --- Opt-in steps: Cyrius (21), Parascopy (35), IPD-KIR (08 with KIR=true) ---
  if [ "$(cat "${GENOME_DIR}/tools/cyrius-${CYRIUS_VERSION}/INSTALLED" 2>/dev/null)" = \
       "python=${PYTHON_IMAGE} lock=$(_digest sha256 "${PGP_ROOT}/scripts/cyrius-constraints.txt")" ]; then
    pass "Cyrius ${CYRIUS_VERSION} (opt-in step 21): installed"
  else
    info "Cyrius ${CYRIUS_VERSION} not installed (opt-in step 21; non-commercial licence): ./scripts/setup.sh --cyrius ${GENOME_DIR}"
  fi
  if [ -s "${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}/homology_table/GRCh38.bed.gz" ]; then
    pass "Parascopy ${PARASCOPY_DATA_VERSION} homology table and models (opt-in step 35): present"
  else
    info "Parascopy data not installed (opt-in step 35): ./scripts/setup.sh --parascopy-data ${GENOME_DIR}"
  fi
  if [ -s "${GENOME_DIR}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat" ]; then
    pass "IPD-KIR ${KIR_DB_RELEASE} (step 08 with KIR=true): present"
  else
    info "IPD-KIR ${KIR_DB_RELEASE} not installed (KIR=true in step 08): ./scripts/setup.sh --kir-data ${GENOME_DIR}"
  fi

  # --- GATK sequence dictionary (optional) ---
  DICT="$REF_DICT"
  if [ -f "$DICT" ]; then
    pass "GATK sequence dictionary (.dict): present"
  else
    warn "GATK sequence dictionary not found at: ${DICT}"
    echo "       Some tools (GATK Mutect2/step 20) require it. Generate with:"
    echo "       docker run --rm -v \"${GENOME_DIR}:/genome\" ${GATK_IMAGE} \\"
    echo "         gatk CreateSequenceDictionary -R ${REF_FASTA_C}"
  fi
fi

###############################################################################
# 4. Docker Images
###############################################################################
header "Docker Images"

if ! command -v "$CONTAINER_ENGINE" &>/dev/null || ! "$CONTAINER_ENGINE" info &>/dev/null 2>&1; then
  info "Skipping Docker image checks (Docker not available)"
else
  # Every *_IMAGE line of versions.env not marked `# optional`: the same list
  # setup.sh pulls.
  IMAGES=()
  while IFS= read -r img; do
    IMAGES+=("$img")
  done < <(pipeline_images)

  PULLED=0
  for img in "${IMAGES[@]}"; do
    if "$CONTAINER_ENGINE" image inspect "$img" &>/dev/null; then
      PULLED=$((PULLED + 1))
    else
      MISSING_IMAGES+=("$img")
    fi
  done

  if [ ${#MISSING_IMAGES[@]} -eq 0 ]; then
    pass "All ${#IMAGES[@]} Docker images are pulled"
  else
    pass "${PULLED}/${#IMAGES[@]} Docker images already pulled"
    fail "${#MISSING_IMAGES[@]} Docker image(s) not yet pulled:"
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
        echo "       docker run --rm -v \"${GENOME_DIR}:/genome\" ${SAMTOOLS_IMAGE} \\"
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
        echo "       docker run --rm -v \"${GENOME_DIR}:/genome\" ${BCFTOOLS_IMAGE} \\"
        echo "         bcftools index -t /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
      fi
    fi

    # Validate genome build (GRCh38) if BAM or VCF exists
    if $HAS_BAM && command -v "$CONTAINER_ENGINE" >/dev/null 2>&1; then
      echo ""
      info "Checking genome build of BAM..."
      BAM_CHR1_LEN=$(run_in "${SAMTOOLS_IMAGE}" \
        samtools view -H "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" 2>/dev/null | \
        grep "^@SQ" | grep "SN:chr1" | head -1 | sed 's/.*LN://' | cut -f1 || echo "0")
      if [ -z "$BAM_CHR1_LEN" ] || [ "$BAM_CHR1_LEN" = "0" ]; then
        # Try without chr prefix (hg19 style)
        BAM_CHR1_LEN=$(run_in "${SAMTOOLS_IMAGE}" \
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

      # Read group: GATK steps (20, 03a, 29) reject reads without one, and
      # DeepVariant takes the sample name from it.
      if BAM_HEADER=$(run_in "${SAMTOOLS_IMAGE}" \
          samtools view -H "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" 2>/dev/null); then
        check_bam_reference "$BAM_HEADER"
        if grep -q '^@RG' <<< "$BAM_HEADER"; then
          pass "BAM has a read group (@RG)"
        else
          warn "BAM header has no @RG read group line. GATK steps (20, 03a, 29) will reject its reads."
          echo "       Add one, then replace the BAM and re-index it:"
          echo "       docker run --rm -v \"${GENOME_DIR}:/genome\" ${SAMTOOLS_IMAGE} \\"
          echo "         samtools addreplacerg -r ID:${SAMPLE} -r SM:${SAMPLE} -r PL:ILLUMINA -r LB:${SAMPLE} \\"
          echo "         -o /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.rg.bam /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
          echo "       mv \"${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.rg.bam\" \"${BAM}\""
          echo "       docker run --rm -v \"${GENOME_DIR}:/genome\" ${SAMTOOLS_IMAGE} \\"
          echo "         samtools index /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
        fi
      else
        fail "Could not read the BAM header (samtools view -H failed): the BAM is unreadable or not a BAM"
      fi
      check_bam_quickcheck
    fi

    if $HAS_VCF && command -v "$CONTAINER_ENGINE" >/dev/null 2>&1; then
      VCF_CONTIG=$(run_in "${BCFTOOLS_IMAGE}" \
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
        VCF_NO_CHR=$(run_in "${BCFTOOLS_IMAGE}" \
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
        echo "       BAM also available — all pipeline steps can run."
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
