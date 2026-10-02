#!/usr/bin/env bash
# setup.sh — One-stop setup: download references, pull Docker images, validate
# Usage: ./scripts/setup.sh <genome_dir>
#        ./scripts/setup.sh --pull-only      pull the Docker images and exit
#        ./scripts/setup.sh --refresh clinvar <genome_dir>
#                                            replace ClinVar with NCBI's current
#                                            release and rebuild the files made from it
#
# This script downloads everything needed to run the pipeline:
#   1. GRCh38 reference genome + index (~3.5 GB)
#   2. ClinVar database (~200 MB), with its release date in clinvar/RELEASE
#   3. All Docker images (~10-15 GB)
#   4. The reference's sequence dictionary and small pinned data files: Delly's
#      exclude map, GRCh38 chromosome bands, an IPD-IMGT/HLA release and the
#      GENCODE gene coordinates T1K needs (~350 MB)
#   5. AnnotSV annotation data for step 5 (~5.3 GB download, ~20 GB unpacked)
#
# VEP cache (~26 GB) and PCGR ref data (~5 GB) are downloaded separately
# because they are only needed for specific steps and take a long time.
#
# Every download goes through fetch (scripts/lib/common.sh): it is written to
# <file>.part, checked, and only then renamed, so a file that exists is whole.
#
# The reference can be replaced by another GRCh38 build: set REF_FASTA (where
# it is stored), REF_FASTA_URL and REF_FASTA_MD5 (and REF_FAI_MD5, or leave it
# empty to build the index with samtools).
set -euo pipefail

PULL_ONLY=false
REFRESH=""
case "${1:-}" in
  --pull-only)
    PULL_ONLY=true
    shift ;;
  --refresh)
    REFRESH=${2:-}
    if [ "$REFRESH" != clinvar ]; then
      echo "Usage: $0 --refresh clinvar <genome_dir>" >&2
      echo "  (clinvar is the one database this flag refreshes)" >&2
      exit 1
    fi
    shift 2 ;;
esac

GENOME_DIR=${1:-${GENOME_DIR:-""}}
if [ -z "$GENOME_DIR" ] && ! $PULL_ONLY; then
  echo "Usage: $0 <genome_dir>"
  echo "       $0 --pull-only"
  echo "       $0 --refresh clinvar <genome_dir>"
  echo ""
  echo "  <genome_dir>  Where to store reference data and sample outputs."
  echo "                Needs at least 500 GB free space per sample."
  echo "  --pull-only   Pull the Docker images listed in versions.env and exit."
  echo "  --refresh clinvar"
  echo "                Download NCBI's current ClinVar, rebuild the files made from it"
  echo "                and record its release date in <genome_dir>/clinvar/RELEASE."
  echo ""
  echo "Example:"
  echo "  ./scripts/setup.sh /data/genomics"
  echo "  ./scripts/setup.sh ~/genome_data"
  exit 1
fi

export GENOME_DIR

# Image tags, data versions and the helpers (run_in, fetch, pipeline_images)
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"

REF_FASTA_URL=${REF_FASTA_URL:-https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta}
# md5 of the two files as the bucket reports them (x-goog-hash).
if [ -z "${REF_FASTA_MD5+x}" ]; then REF_FASTA_MD5=7ff134953dcca8c8997453bbb80b6b5e; fi
if [ -z "${REF_FAI_MD5+x}" ]; then REF_FAI_MD5=f76371b113734a56cde236bc0372de0a; fi

# pull_images: pull every image setup pre-pulls (versions.env, minus the
# lines marked `# optional`). Sets PULLED, SKIPPED and FAILED.
pull_images() {
  local img
  PULLED=0
  SKIPPED=0
  FAILED=0
  while IFS= read -r img; do
    if "$CONTAINER_ENGINE" image inspect "$img" &>/dev/null; then
      SKIPPED=$((SKIPPED + 1))
    else
      echo "  Pulling: ${img}..."
      if "$CONTAINER_ENGINE" pull "$img" 2>/dev/null; then
        PULLED=$((PULLED + 1))
      else
        echo "  WARNING: Failed to pull ${img}. Check the image name/tag."
        FAILED=$((FAILED + 1))
      fi
    fi
  done < <(pipeline_images)
  echo "[OK] Docker images: ${PULLED} pulled, ${SKIPPED} already present, ${FAILED} failed."
  if [ "$FAILED" -gt 0 ]; then
    echo ""
    echo "ERROR: ${FAILED} Docker image(s) failed to pull. Fix the issues above and re-run setup."
    return 1
  fi
}

# --- ClinVar ---------------------------------------------------------------------
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVARDIR="${GENOME_DIR:+${GENOME_DIR}/clinvar}"
# The raw download and the two files built from it, each with its index.
CLINVAR_SET="clinvar.vcf.gz clinvar.vcf.gz.tbi clinvar_chr.vcf.gz clinvar_chr.vcf.gz.tbi
  clinvar_pathogenic_chr.vcf.gz clinvar_pathogenic_chr.vcf.gz.tbi"

# clinvar_build_derived RAW_DIR OUT_DIR: from RAW_DIR/clinvar.vcf.gz build, in
# OUT_DIR (both under GENOME_DIR), clinvar_chr.vcf.gz (NCBI's 1, 2 ... MT
# renamed to chr1, chr2 ... chrM) and clinvar_pathogenic_chr.vcf.gz (the
# Pathogenic and Likely_pathogenic records), each with its .tbi.
clinvar_build_derived() {
  local raw out
  raw=$(cpath "$1") && out=$(cpath "$2") || return 2
  awk 'BEGIN { for (i = 1; i <= 22; i++) print i, "chr" i; print "X chrX"; print "Y chrY"; print "MT chrM" }' \
    > "${2}/chr_rename.txt"
  echo "  Creating chr-prefixed ClinVar..."
  run_in --rw "$2" "$BCFTOOLS_IMAGE" \
    bcftools annotate --rename-chrs "${out}/chr_rename.txt" "${raw}/clinvar.vcf.gz" \
      -Oz -o "${out}/clinvar_chr.vcf.gz" || return 1
  run_in --rw "$2" "$BCFTOOLS_IMAGE" bcftools index -f -t "${out}/clinvar_chr.vcf.gz" || return 1
  echo "  Creating pathogenic/likely pathogenic subset..."
  run_in --rw "$2" "$BCFTOOLS_IMAGE" \
    bcftools view -i 'CLNSIG~"Pathogenic" || CLNSIG~"Likely_pathogenic"' "${out}/clinvar_chr.vcf.gz" \
      -Oz -o "${out}/clinvar_pathogenic_chr.vcf.gz" || return 1
  run_in --rw "$2" "$BCFTOOLS_IMAGE" bcftools index -f -t "${out}/clinvar_pathogenic_chr.vcf.gz" || return 1
  rm -f "${2}/chr_rename.txt"
}

# clinvar_install DIR: move every file of CLINVAR_SET that DIR holds into
# clinvar/, then remove DIR. The data files go first and their indexes after.
clinvar_install() {
  local f
  for f in $CLINVAR_SET; do
    case "$f" in *.tbi) continue ;; esac
    if [ -f "${1}/${f}" ]; then mv -f "${1}/${f}" "${CLINVARDIR}/${f}"; fi
  done
  for f in $CLINVAR_SET; do
    case "$f" in *.tbi) ;; *) continue ;; esac
    if [ -f "${1}/${f}" ]; then mv -f "${1}/${f}" "${CLINVARDIR}/${f}"; fi
  done
  rm -rf "$1"
}

# clinvar_record_release: write the ##fileDate of clinvar/clinvar.vcf.gz
# (YYYY-MM-DD) to clinvar/RELEASE, which validate-setup.sh and step 06 print.
clinvar_record_release() {
  local d
  d=$(gzip -cd "${CLINVARDIR}/clinvar.vcf.gz" 2>/dev/null | head -n 100 \
      | awk -F= '/^##fileDate=/ { print $2; exit }') || true
  case "$d" in
    [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]) d="${d:0:4}-${d:4:2}-${d:6:2}" ;;
  esac
  case "$d" in
    [0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9])
      printf '%s\n' "$d" > "${CLINVARDIR}/RELEASE.tmp"
      mv -f "${CLINVARDIR}/RELEASE.tmp" "${CLINVARDIR}/RELEASE"
      echo "[OK] ClinVar release ${d} (recorded in ${CLINVARDIR}/RELEASE)." ;;
    *)
      echo "[WARN] ${CLINVARDIR}/clinvar.vcf.gz has no ##fileDate line; its release date is unknown." ;;
  esac
}

# refresh_clinvar: download NCBI's current ClinVar into clinvar/.refresh/,
# check it against NCBI's md5, build both derived files from it there, and only
# then replace all six files. The normalised copy step 06 keeps is removed, so
# step 06 builds it again from the new release.
refresh_clinvar() {
  local new="${CLINVARDIR}/.refresh"
  # A download cut short last time (*.part) is resumed; anything else is redone.
  mkdir -p "$new"
  find "$new" -mindepth 1 ! -name '*.part' -exec rm -rf {} +
  echo "Downloading the current ClinVar release..."
  fetch "$CLINVAR_URL" "${new}/clinvar.vcf.gz" md5 "${CLINVAR_URL}.md5" || return 1
  fetch "${CLINVAR_URL}.tbi" "${new}/clinvar.vcf.gz.tbi" || return 1
  clinvar_build_derived "$new" "$new" || {
    echo "ERROR: could not build the ClinVar files from the new release; the old ones are kept." >&2
    return 1
  }
  clinvar_install "$new"
  rm -f "${CLINVARDIR}/clinvar_pathogenic_chr.norm.vcf.gz" "${CLINVARDIR}/clinvar_pathogenic_chr.norm.vcf.gz.tbi"
  clinvar_record_release
}

echo "============================================"
echo "  Personal Genome Pipeline — Setup"
echo "  Data directory: ${GENOME_DIR:-(none: --pull-only)}"
echo "============================================"
echo ""

# Check Docker
if ! command -v "$CONTAINER_ENGINE" &>/dev/null; then
  echo "ERROR: Docker is not installed."
  echo "  Install Docker: https://docs.docker.com/get-docker/"
  exit 1
fi
if ! "$CONTAINER_ENGINE" info &>/dev/null; then
  echo "ERROR: Docker daemon is not running."
  echo "  Start Docker Desktop or run: sudo systemctl start docker"
  exit 1
fi
echo "[OK] Docker is running."

if $PULL_ONLY; then
  echo ""
  echo "=== Docker Images (~10-15 GB total) ==="
  pull_images || exit 1
  exit 0
fi

if [ "$REFRESH" = clinvar ]; then
  mkdir -p "$CLINVARDIR"
  echo ""
  echo "=== Refreshing ClinVar ==="
  if [ -f "${CLINVARDIR}/RELEASE" ]; then echo "  Current release: $(cat "${CLINVARDIR}/RELEASE")"; fi
  refresh_clinvar || exit 1
  echo "[OK] ClinVar refreshed: raw file, chr-prefixed file and pathogenic subset replaced."
  exit 0
fi

# Genomes are private data: only the owner may enter the data directory.
mkdir -p "$GENOME_DIR"
chmod 700 "$GENOME_DIR"

# Check disk space (df -Pk works on Linux and macOS; GENOME_DIR may not
# exist yet, so measure its nearest existing parent)
DF_DIR="$GENOME_DIR"
while [ ! -d "$DF_DIR" ]; do
  DF_DIR=$(dirname "$DF_DIR")
done
AVAIL_KB=$(df -Pk "$DF_DIR" | awk 'NR==2 {print $4}')
AVAIL_GB=$(( ${AVAIL_KB:-0} / 1048576 ))
if [ "$AVAIL_GB" -lt 50 ]; then
  echo "WARNING: Only ${AVAIL_GB} GB free in ${GENOME_DIR}. Need at least 50 GB for references."
fi
echo "[OK] ${AVAIL_GB} GB free in ${GENOME_DIR}."

###############################################################################
# Phase 1: Reference Genome
###############################################################################
echo ""
echo "=== Phase 1: Reference Genome (~3.5 GB) ==="

FASTA="$REF_FASTA"
FAI="${REF_FASTA}.fai"
REFDIR=$(dirname "$FASTA")
mkdir -p "$REFDIR"

if [ -f "$FASTA" ] && [ -f "$FAI" ]; then
  echo "[OK] Reference genome already downloaded."
else
  echo "Downloading GRCh38 reference genome..."
  echo "  Source: ${REF_FASTA_URL}"
  echo "  Size: ~3.1 GB (FASTA) + ~2 MB (index)"

  if [ ! -f "$FASTA" ]; then
    # shellcheck disable=SC2046  # the checksum arguments are dropped when REF_FASTA_MD5 is empty
    fetch "$REF_FASTA_URL" "$FASTA" $([ -z "$REF_FASTA_MD5" ] || printf 'md5 %s' "$REF_FASTA_MD5") || {
      echo "  The download is resumed when setup.sh runs again."
      exit 1
    }
  fi

  if [ ! -f "$FAI" ]; then
    if [ -z "$REF_FAI_MD5" ] || ! fetch "${REF_FASTA_URL}.fai" "$FAI" md5 "$REF_FAI_MD5"; then
      echo "  Generating index with samtools..."
      # The index goes next to the FASTA, so its directory is writable here.
      run_in --rw "$REFDIR" "$SAMTOOLS_IMAGE" \
        samtools faidx "$REF_FASTA_C"
    fi
  fi
  echo "[OK] Reference genome downloaded."
fi

###############################################################################
# Phase 2: ClinVar Database
###############################################################################
echo ""
echo "=== Phase 2: ClinVar Database (~200 MB) ==="

mkdir -p "$CLINVARDIR"

CLINVAR="${CLINVARDIR}/clinvar.vcf.gz"
CLINVAR_TBI="${CLINVARDIR}/clinvar.vcf.gz.tbi"

if [ -f "$CLINVAR" ] && [ -f "$CLINVAR_TBI" ]; then
  echo "[OK] ClinVar database already downloaded. To replace it with the current release:"
  echo "  ./scripts/setup.sh --refresh clinvar ${GENOME_DIR}"
else
  echo "Downloading ClinVar database..."
  if [ ! -f "$CLINVAR" ]; then
    # NCBI publishes the md5 next to the file.
    fetch "$CLINVAR_URL" "$CLINVAR" md5 "${CLINVAR_URL}.md5" || exit 1
  fi
  if [ ! -f "$CLINVAR_TBI" ]; then
    fetch "${CLINVAR_URL}.tbi" "$CLINVAR_TBI" || exit 1
  fi
fi

# The two derived files are built whenever one is missing, also on a rerun
# with the raw download already present. They are built in clinvar/.build/
# and moved in when both are complete, so a half-built file is never in place.
if [ ! -f "${CLINVARDIR}/clinvar_chr.vcf.gz.tbi" ] || [ ! -f "${CLINVARDIR}/clinvar_pathogenic_chr.vcf.gz.tbi" ]; then
  rm -rf "${CLINVARDIR}/.build"
  mkdir -p "${CLINVARDIR}/.build"
  clinvar_build_derived "$CLINVARDIR" "${CLINVARDIR}/.build" || exit 1
  clinvar_install "${CLINVARDIR}/.build"
fi
if [ ! -f "${CLINVARDIR}/RELEASE" ]; then clinvar_record_release; fi
echo "[OK] ClinVar database and its chr-prefixed and pathogenic subsets are ready."

###############################################################################
# Phase 3: Docker Images
###############################################################################
echo ""
echo "=== Phase 3: Docker Images (~10-15 GB total) ==="

# The list is every *_IMAGE line of versions.env not marked `# optional`.
pull_images || exit 1

###############################################################################
# Phase 4: Sequence dictionary and small pinned data files
###############################################################################
echo ""
echo "=== Phase 4: Sequence dictionary and pinned data files (~350 MB) ==="

# GATK and Picard (steps 03a, 20, 29, chip-to-vcf) need the reference's .dict.
if [ -f "$REF_DICT" ]; then
  echo "[OK] Sequence dictionary already present."
else
  echo "Creating the sequence dictionary ${REF_DICT}..."
  DICT_PART="${REF_DICT%.dict}.part.dict"
  rm -f "$DICT_PART"
  # The dictionary goes next to the FASTA, so its directory is writable here.
  if run_in --rw "$REFDIR" "$GATK_IMAGE" \
       gatk CreateSequenceDictionary -R "$REF_FASTA_C" -O "$(cpath "$DICT_PART")" \
     && [ -s "$DICT_PART" ]; then
    mv -f "$DICT_PART" "$REF_DICT"
    echo "[OK] Sequence dictionary created."
  else
    rm -f "$DICT_PART"
    echo "[WARN] Could not create ${REF_DICT}. Steps 03a, 20, 29 and chip-to-vcf need it; re-run setup.sh."
  fi
fi

# Each from a fixed commit or release and checked before it is stored
# (scripts/lib/common.sh, install_data_file). A failed download does not stop
# setup: the step that needs the file says so, and setup.sh fetches it on its
# next run.
for name in $DATA_FILES; do
  case "$name" in
    delly_exclude) what="Delly exclude map (step 19)" ;;
    cytoband) what="GRCh38 chromosome bands (step 10)" ;;
    hla_dat) what="IPD-IMGT/HLA ${HLA_DB_RELEASE} (step 08)" ;;
    gencode_genes) what="GENCODE ${GENCODE_RELEASE} gene coordinates (step 08)" ;;
  esac
  if dest=$(data_file "$name"); then
    echo "[OK] ${what} already present."
  elif install_data_file "$name"; then
    echo "[OK] ${what}: ${dest}"
  else
    echo "[WARN] Could not install the ${what}. Re-run setup.sh to try again."
  fi
done

###############################################################################
# Phase 5: AnnotSV annotation data (step 5, a default step)
###############################################################################
# The AnnotSV image holds code only; without this data AnnotSV exits with an
# error. Downloaded under a .part name, extracted into a temporary directory
# and moved into place only when complete.
echo ""
echo "=== Phase 5: AnnotSV annotation data (~5.3 GB download, ~20 GB unpacked) ==="

ANNOTSV_DIR="${GENOME_DIR}/annotsv_annotations"
ANNOTSV_NAME="Annotations_Human_${ANNOTSV_ANNOTATIONS_VERSION}.tar.gz"
ANNOTSV_URL="https://www.lbgi.fr/~geoffroy/Annotations/${ANNOTSV_NAME}"
ANNOTSV_TARBALL="${GENOME_DIR}/${ANNOTSV_NAME}"
if [ -d "${ANNOTSV_DIR}/Annotations_Human/Genes/GRCh38" ]; then
  echo "[OK] AnnotSV annotation data already present."
else
  echo "Downloading AnnotSV ${ANNOTSV_ANNOTATIONS_VERSION} annotation data (~5.3 GB)..."
  echo "  The AnnotSV server is slow (about 0.8 MB/s measured from a GitHub runner), so this can take"
  echo "  1-2 hours. An interrupted download resumes when setup.sh runs again."
  if fetch "$ANNOTSV_URL" "$ANNOTSV_TARBALL"; then
    echo "  Extracting..."
    rm -rf "${ANNOTSV_DIR}.part"
    mkdir -p "${ANNOTSV_DIR}.part"
    if tar -xzf "$ANNOTSV_TARBALL" -C "${ANNOTSV_DIR}.part"; then
      rm -rf "$ANNOTSV_DIR"
      mv "${ANNOTSV_DIR}.part" "$ANNOTSV_DIR"
      rm -f "$ANNOTSV_TARBALL"
      echo "[OK] AnnotSV annotation data: ${ANNOTSV_DIR}/ ($(du -sh "$ANNOTSV_DIR" | cut -f1))"
    else
      rm -rf "${ANNOTSV_DIR}.part"
      echo "[WARN] Could not extract ${ANNOTSV_TARBALL}. Step 5 (AnnotSV) will be skipped until it is in place."
    fi
  else
    echo "[WARN] AnnotSV annotation download failed. Step 5 (AnnotSV) will be skipped until it is in place."
    echo "  Re-run setup.sh, or download it by hand:"
    echo "    curl -fL -C - -o ${ANNOTSV_TARBALL} ${ANNOTSV_URL}"
    echo "    mkdir -p ${ANNOTSV_DIR} && tar -xzf ${ANNOTSV_TARBALL} -C ${ANNOTSV_DIR}"
  fi
fi

###############################################################################
# Phase 6: Optional Downloads (instructions only)
###############################################################################
echo ""
echo "=== Phase 6: Optional Downloads (manual) ==="
echo ""
echo "The following are only needed for specific steps and are large downloads:"
echo ""

VEPDIR="${GENOME_DIR}/vep_cache"
if [ -f "${VEPDIR}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
  echo "[OK] VEP ${VEP_CACHE_RELEASE} cache already present."
else
  echo "[SKIP] VEP ${VEP_CACHE_RELEASE} cache (~26 GB) — needed for step 13 (VEP annotation)"
  echo "  Step 13 downloads, checks and unpacks it the first time it runs:"
  echo "    ./scripts/13-vep-annotation.sh <sample_name>"
fi
echo ""

PCGRDIR="${GENOME_DIR}/pcgr_data"
if [ -d "${PCGRDIR}/${PCGR_DATA_BUNDLE}/data" ]; then
  echo "[OK] PCGR/CPSR ref data bundle already present."
else
  echo "[SKIP] PCGR/CPSR ref data (~5 GB download) — needed for step 17 (cancer predisposition)"
  echo "  Download:"
  echo "    mkdir -p ${PCGRDIR} && cd ${PCGRDIR}"
  echo "    curl -fL -C - -O https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.${PCGR_DATA_BUNDLE}.grch38.tgz"
  echo "    tar xzf pcgr_ref_data.${PCGR_DATA_BUNDLE}.grch38.tgz"
  echo "    mkdir -p ${PCGR_DATA_BUNDLE} && mv data/ ${PCGR_DATA_BUNDLE}/"
fi
echo ""

# CPSR runs the VEP release inside the PCGR image, not the one step 13 uses.
if [ -f "${VEPDIR}/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
  echo "[OK] VEP ${PCGR_VEP_CACHE_RELEASE} cache (for CPSR) already present."
else
  echo "[SKIP] VEP ${PCGR_VEP_CACHE_RELEASE} cache (~26 GB) — needed for step 17 (CPSR runs VEP ${PCGR_VEP_CACHE_RELEASE} inside the PCGR image)"
  echo "  This is separate from the release-${VEP_CACHE_RELEASE} cache used by step 13. Both coexist in vep_cache/."
  echo "  Download:"
  echo "    mkdir -p ${VEPDIR}"
  echo "    curl -fL -C - -o ${VEPDIR}/$(basename "$(vep_cache_url "$PCGR_VEP_CACHE_RELEASE")") $(vep_cache_url "$PCGR_VEP_CACHE_RELEASE")"
  echo "    tar xzf ${VEPDIR}/$(basename "$(vep_cache_url "$PCGR_VEP_CACHE_RELEASE")") -C ${VEPDIR}"
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
