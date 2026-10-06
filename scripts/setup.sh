#!/usr/bin/env bash
# setup.sh — One-stop setup: download references, pull Docker images, validate
# Usage: ./scripts/setup.sh <genome_dir>
#        ./scripts/setup.sh --pull-only      pull the Docker images and exit
#        ./scripts/setup.sh --refresh clinvar <genome_dir>
#                                            replace ClinVar with NCBI's current
#                                            release and rebuild the files made from it
#        ./scripts/setup.sh --cyrius <genome_dir>
#        ./scripts/setup.sh --parascopy-data <genome_dir>
#        ./scripts/setup.sh --kir-data <genome_dir>
#        ./scripts/setup.sh --ancestry-panel <genome_dir>
#                                            what an opt-in step needs, and nothing
#                                            else: Cyrius (step 21), Parascopy's
#                                            homology table and models (step 35),
#                                            IPD-KIR (step 08 with KIR=true),
#                                            pgsc_calc's ancestry panel (step 26,
#                                            and percentiles in step 25)
#
# This script downloads everything needed to run the pipeline:
#   1. GRCh38 reference genome + index: NCBI's GRCh38 no-ALT analysis set
#      (~0.9 GB download, ~3.2 GB unpacked)
#   2. ClinVar database (~200 MB), with its release date in clinvar/RELEASE
#   3. All Docker images (~10-15 GB)
#   4. The reference's sequence dictionary and small pinned data files: Delly's
#      exclude map, GRCh38 chromosome bands, an IPD-IMGT/HLA release, the
#      GENCODE gene coordinates T1K needs (~350 MB), and somalier's sites and
#      VerifyBamID2's marker panel for step 33 (~10 MB)
#   5. AnnotSV annotation data for step 5 (~5.3 GB download, ~20 GB unpacked)
#   6. Step 25: the PGS Catalog scores of assets/pgs_scores.tsv (~400 MB), a
#      checkout of pgsc_calc and the Nextflow plugin it needs, so step 25 and
#      the PRS process run without the network
#
# VEP cache (~26 GB) and PCGR ref data (~7 GB) are downloaded separately
# because they are only needed for specific steps and take a long time.
#
# Every download goes through fetch (scripts/lib/common.sh): it is written to
# <file>.part, checked, and only then renamed, so a file that exists is whole.
#
# The reference can be replaced by another GRCh38 build: set REF_FASTA (where
# it is stored), REF_FASTA_URL (a .gz URL is unpacked) and REF_FASTA_MD5 (an
# md5, the URL of a checksum file that lists the download, or empty for no
# check), and REF_FAI_MD5 the same way for the .fai published beside it (empty
# builds the index with samtools). docs/00-reference-setup.md lists the builds
# that work and docs/realignment.md what a change of reference means for
# existing samples.
set -euo pipefail

PULL_ONLY=false
REFRESH=""
SAMPLE_QC_ONLY=false
OPT_IN=""
case "${1:-}" in
  --cyrius|--parascopy-data|--kir-data|--ancestry-panel)
    OPT_IN=${1#--}
    shift ;;
  --pull-only)
    PULL_ONLY=true
    shift ;;
  --sample-qc-data)
    SAMPLE_QC_ONLY=true
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
  echo "       $0 --sample-qc-data <genome_dir>"
  echo "       $0 --cyrius | --parascopy-data | --kir-data | --ancestry-panel <genome_dir>"
  echo ""
  echo "  <genome_dir>  Where to store reference data and sample outputs."
  echo "                Needs at least 500 GB free space per sample."
  echo "  --pull-only   Pull the Docker images listed in versions.env and exit."
  echo "  --refresh clinvar"
  echo "                Download NCBI's current ClinVar, rebuild the files made from it"
  echo "                and record its release date in <genome_dir>/clinvar/RELEASE."
  echo "  --sample-qc-data"
  echo "                Install only somalier's sites and VerifyBamID2's panel (step 33) and exit."
  echo "  --cyrius      Install Cyrius for the opt-in step 21 (PyPI, hash-locked) and exit."
  echo "                Cyrius is under the PolyForm Strict licence: non-commercial use only."
  echo "  --parascopy-data"
  echo "                Install Parascopy's GRCh38 homology table and models (step 35, ~50 MB) and exit."
  echo "  --kir-data    Install the IPD-KIR ${KIR_DB_RELEASE:-} database (KIR=true in step 08, ~40 MB) and exit."
  echo "  --ancestry-panel"
  echo "                Install pgsc_calc's ancestry reference panel ${PGSC_PANEL:-} (~7 GB download): step 26, and"
  echo "                percentiles instead of raw scores in step 25; then exit."
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

# NCBI's GRCh38 analysis set without ALT contigs: chr1-22, X, Y, M, the
# unplaced and unlocalized scaffolds and chrEBV, 195 sequences. NCBI lists the
# md5 of every file of the directory in md5checksums.txt, which is read at
# download time, and publishes the .fai beside the FASTA.
REF_FASTA_URL=${REF_FASTA_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz}
REF_FAI_URL="${REF_FASTA_URL%.gz}.fai"
if [ -z "${REF_FASTA_MD5+x}" ]; then REF_FASTA_MD5="$(dirname "$REF_FASTA_URL")/md5checksums.txt"; fi
if [ -z "${REF_FAI_MD5+x}" ]; then REF_FAI_MD5="$(dirname "$REF_FASTA_URL")/md5checksums.txt"; fi

# --- Step 33: somalier's sites and VerifyBamID2's marker panel ----------------------
# Pinned files, each checked against the sha256 recorded here. Step 33
# (scripts/33-sample-qc.sh) and validate-setup.sh read them under these names.
SOMALIER_SITES_URL=https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
SOMALIER_SITES_SHA256=d1a853b8bb2e5f1a520bc67c8303be699543d225d785bce51e8335a2420489b5
# The panel of the VerifyBamID release VERIFYBAMID2_IMAGE packages (v2.0.3):
# one line per file, its URL and its sha256.
VB2_PANEL_FILES="https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD 259e320123756bb702542e3b3c4d766b481426d96cd63abd4b9f9fced035f375
https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.mu f7e8b4fad17cc433d887f18bf183fa4b89597d5a375f22133f44f39c5e35ffc8
https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.bed 0025d782137e5906bbc7afc553c49d2c82c5b8f0bfa0649bc415887442f42869"

# install_sample_qc_data: install both under GENOME_DIR/reference unless they
# are there. Returns 1 when a download fails its check.
install_sample_qc_data() {
  local dest="${GENOME_DIR}/reference" url sha rc=0
  if [ -s "${dest}/somalier/sites.hg38.vcf.gz" ]; then
    echo "[OK] somalier sites (step 33) already present."
  elif fetch "$SOMALIER_SITES_URL" "${dest}/somalier/sites.hg38.vcf.gz" sha256 "$SOMALIER_SITES_SHA256"; then
    echo "[OK] somalier sites (step 33): ${dest}/somalier/sites.hg38.vcf.gz"
  else
    echo "[WARN] Could not install somalier's sites. Step 33 needs them; re-run setup.sh."
    rc=1
  fi
  while read -r url sha; do
    if [ -s "${dest}/verifybamid2/$(basename "$url")" ]; then
      continue
    fi
    if ! fetch "$url" "${dest}/verifybamid2/$(basename "$url")" sha256 "$sha"; then
      echo "[WARN] Could not install $(basename "$url") of VerifyBamID2's panel. Step 33 needs it; re-run setup.sh."
      rc=1
    fi
  done <<<"$VB2_PANEL_FILES"
  [ "$rc" -eq 0 ] && echo "[OK] VerifyBamID2 panel (step 33): ${dest}/verifybamid2/"
  return "$rc"
}

if $SAMPLE_QC_ONLY; then
  install_sample_qc_data
  exit $?
fi

# --- Opt-in steps: Cyrius, Parascopy, KIR ------------------------------------------
# None of these is installed by a plain setup.sh run: each belongs to a step a
# default run leaves out, and each is asked for by its own flag.

# install_cyrius: Cyrius (step 21) into GENOME_DIR/tools/cyrius-<version>,
# installed by pip in PYTHON_IMAGE, the image step 21 runs it in, from
# scripts/cyrius-constraints.txt: every package pinned with the sha256 of its
# wheels (--require-hashes), nothing resolved (--no-deps) and nothing built
# (--only-binary). The network is used here, never when step 21 runs. Cyrius
# is under the PolyForm Strict licence 1.0.0 (non-commercial use only).
install_cyrius() {
  local dest="${GENOME_DIR}/tools/cyrius-${CYRIUS_VERSION}" lock="${PGP_ROOT}/scripts/cyrius-constraints.txt" stamp
  stamp="python=${PYTHON_IMAGE} lock=$(_digest sha256 "$lock")"
  if [ "$(cat "${dest}/INSTALLED" 2>/dev/null)" = "$stamp" ]; then
    echo "[OK] Cyrius ${CYRIUS_VERSION} (step 21) already installed: ${dest}"
    return 0
  fi
  echo "Installing Cyrius ${CYRIUS_VERSION} from ${lock} (hash-locked) into ${dest}"
  echo "  Cyrius is under the PolyForm Strict licence 1.0.0: non-commercial use only."
  rm -rf "${dest}.part"
  mkdir -p "${dest}.part"
  # --net: pip downloads the locked wheels from PyPI.
  if ! run_in --net --rw "${GENOME_DIR}/tools" -v "${lock}:/lock.txt:ro" "${PYTHON_IMAGE}" \
      pip install --no-cache-dir --disable-pip-version-check -q --require-hashes --no-deps \
        --only-binary :all: --target "$(cpath "${dest}.part")" -r /lock.txt; then
    rm -rf "${dest}.part"
    echo "[WARN] Could not install Cyrius. Step 21 needs it; run: $0 --cyrius ${GENOME_DIR}"
    return 1
  fi
  echo "$stamp" > "${dest}.part/INSTALLED"
  rm -rf "$dest"
  mv "${dest}.part" "$dest"
  echo "[OK] Cyrius ${CYRIUS_VERSION} (step 21): ${dest}"
}

# install_parascopy_data: Parascopy's precomputed GRCh38 homology table and
# the model parameters of 1000 Genomes populations (step 35), Zenodo record
# PARASCOPY_DATA_RECORD, each archive checked against its md5 there.
PARASCOPY_DATA_RECORD=15019940
PARASCOPY_FILES="GRCh38_v${PARASCOPY_DATA_VERSION}.tar.gz a95bf674f43317d3a4c1b8ddbb140945
models_GRCh38_1KGP_v${PARASCOPY_DATA_VERSION}.tar.gz 244110d8fa883cf11334527ec2383498"
install_parascopy_data() {
  local dest="${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}" name md5
  if [ -s "${dest}/homology_table/GRCh38.bed.gz" ] && [ -s "${dest}/models_GRCh38_1KGP/EUR/SMN1.gz" ]; then
    echo "[OK] Parascopy ${PARASCOPY_DATA_VERSION} homology table and models (step 35) already present."
    return 0
  fi
  rm -rf "${dest}.part"
  mkdir -p "${dest}.part"
  while read -r name md5; do
    if ! fetch "https://zenodo.org/records/${PARASCOPY_DATA_RECORD}/files/${name}" "${dest}.part/${name}" md5 "$md5" \
        || ! tar -xzf "${dest}.part/${name}" -C "${dest}.part"; then
      rm -rf "${dest}.part"
      echo "[WARN] Could not install Parascopy's data. Step 35 needs it; run: $0 --parascopy-data ${GENOME_DIR}"
      return 1
    fi
    rm -f "${dest}.part/${name}"
  done <<<"$PARASCOPY_FILES"
  rm -rf "$dest"
  mv "${dest}.part" "$dest"
  echo "[OK] Parascopy ${PARASCOPY_DATA_VERSION} homology table and models (step 35): ${dest}"
}

# install_kir_data: kir.dat of IPD-KIR release KIR_DB_RELEASE (step 08 with
# KIR=true), from the release branch of the IPD-KIR repository (2.15.0 is
# branch 2150), checked against the md5 that release lists for it.
install_kir_data() {
  local branch=${KIR_DB_RELEASE//./} dest="${GENOME_DIR}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat" url md5
  if [ -s "$dest" ]; then
    echo "[OK] IPD-KIR ${KIR_DB_RELEASE} (step 08, KIR=true) already present."
    return 0
  fi
  url="https://raw.githubusercontent.com/ANHIG/IPDKIR/${branch}"
  # md5checksum.txt lines read "MD5 (kir.dat) = <md5>".
  md5=$(_get "${url}/md5checksum.txt" - 2>/dev/null | awk '$2 == "(kir.dat)" {print $NF}') || md5=""
  if [ -z "$md5" ] || ! fetch "${url}/kir.dat" "$dest" md5 "$md5"; then
    echo "[WARN] Could not install IPD-KIR ${KIR_DB_RELEASE} (is it a release with a branch ${branch}?)."
    return 1
  fi
  if ! grep -q "IPD-KIR Release Version ${KIR_DB_RELEASE}" "$dest"; then
    rm -f "$dest"
    echo "[WARN] ${url}/kir.dat does not say IPD-KIR Release Version ${KIR_DB_RELEASE}; removed it."
    return 1
  fi
  echo "[OK] IPD-KIR ${KIR_DB_RELEASE} (step 08, KIR=true): ${dest}"
}

# --- Step 25 and 26: pgsc_calc, its scores and its ancestry panel ---------------------
PGS_BASE_URL=${PGS_BASE_URL:-https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores}
PGSC_RESOURCES=https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources

# install_pgsc_calc DIR: GitHub's archive of pgsc_calc PGSC_CALC_VERSION,
# checked against PGSC_CALC_SHA256, unpacked into DIR. Step 25 holds the same.
install_pgsc_calc() {
  local dir=$1 tgz="${1}.tar.gz"
  fetch "https://github.com/PGScatalog/pgsc_calc/archive/refs/tags/${PGSC_CALC_VERSION}.tar.gz" "$tgz" \
      sha256 "$PGSC_CALC_SHA256" || return 1
  rm -rf "${dir}.part"
  mkdir -p "${dir}.part"
  tar -xzf "$tgz" -C "${dir}.part" --strip-components 1 || { rm -rf "${dir}.part"; return 1; }
  rm -f "$tgz"
  # An incomplete folder left by an earlier run would receive the new one inside it.
  rm -rf "$dir"
  mv "${dir}.part" "$dir"
}

# install_prs: the GRCh38-harmonised file of every score in
# assets/pgs_scores.tsv (each checked against the md5 the PGS Catalog
# publishes beside it), pgsc_calc PGSC_CALC_VERSION in tools/, and the
# nf-schema plugin its nextflow.config names. With them
# step 25 and the PRS process need no network. A failure is a warning: step
# 25 fetches what is missing the first time it runs.
install_prs() {
  local id url dest rc=0 calc="${GENOME_DIR}/tools/pgsc_calc-${PGSC_CALC_VERSION}"
  while read -r id; do
    dest="${GENOME_DIR}/prs_scores/${id}.txt.gz"
    [ -s "$dest" ] && continue
    url="${PGS_BASE_URL}/${id}/ScoringFiles/Harmonized/${id}_hmPOS_GRCh38.txt.gz"
    if ! fetch "$url" "$dest" md5 "${url}.md5"; then
      echo "[WARN] Could not download ${id} (step 25 tries again when it runs): ${url}"
      rc=1
    fi
  done < <(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${PGP_ROOT}/assets/pgs_scores.tsv")
  [ "$rc" -eq 0 ] && echo "[OK] PGS Catalog scores (step 25): ${GENOME_DIR}/prs_scores/"
  if [ -f "${calc}/main.nf" ]; then
    echo "[OK] pgsc_calc ${PGSC_CALC_VERSION} (step 25) already present."
  elif install_pgsc_calc "$calc"; then
    echo "[OK] pgsc_calc ${PGSC_CALC_VERSION} (step 25): ${calc}"
  else
    echo "[WARN] Could not fetch pgsc_calc ${PGSC_CALC_VERSION}; step 25 tries again when it runs."
  fi
  if ! command -v nextflow >/dev/null 2>&1; then
    echo "[WARN] Nextflow is not installed: run-all.sh and step 25 need it (docs/nextflow.md)."
  elif [ -d "${NXF_HOME:-${HOME}/.nextflow}/plugins/nf-schema-${PGSC_CALC_NF_SCHEMA}" ]; then
    echo "[OK] Nextflow plugin nf-schema ${PGSC_CALC_NF_SCHEMA} (pgsc_calc) already present."
  elif nextflow plugin install "nf-schema@${PGSC_CALC_NF_SCHEMA}" >/dev/null; then
    echo "[OK] Nextflow plugin nf-schema ${PGSC_CALC_NF_SCHEMA} (pgsc_calc)."
  else
    echo "[WARN] Could not install the nf-schema ${PGSC_CALC_NF_SCHEMA} plugin; step 25 tries again when it runs."
  fi
  return 0
}

# install_ancestry_panel: pgsc_calc's reference panel (PGSC_PANEL, the 1000
# Genomes database; ANCESTRY_PANEL_NAME picks another file of the PGS
# Catalog's resources folder, the tests use its small synthetic one), checked
# against the md5 that folder lists, and beside it the panel's GRCh38
# biallelic SNVs (chrN, position, REF, ALT), read from the panel's own .pvar:
# step 25 genotypes them from the gVCF, so the projection counts the sites
# where the sample matches the reference.
install_ancestry_panel() {
  local name=${ANCESTRY_PANEL_NAME:-$PGSC_PANEL} dir="${GENOME_DIR}/reference/pgsc_calc" panel sites
  panel="${dir}/${name}.tar.zst"
  sites="${dir}/${name}_GRCh38_sites.tsv"
  if [ -s "$panel" ] && [ -s "$sites" ]; then
    echo "[OK] ancestry panel ${name} (steps 25 and 26) already present."
    return 0
  fi
  echo "Downloading pgsc_calc's ancestry panel ${name} (the 1000 Genomes panel is about 7 GB; an interrupted download resumes)..."
  if ! fetch "${PGSC_RESOURCES}/${name}.tar.zst" "$panel" md5 "${PGSC_RESOURCES}/md5s.txt"; then
    echo "[WARN] Could not download the ancestry panel; run: $0 --ancestry-panel ${GENOME_DIR}"
    return 1
  fi
  echo "  Listing the panel's GRCh38 SNVs..."
  # shellcheck disable=SC2016  # $1 belongs to the inner sh
  if ! run_in "$PGSC_ZSTD_IMAGE" sh -c 'tar -xOf "$1" --wildcards "GRCh38_*_ALL.pvar.zst" | zstd -dc' _ "$(cpath "$panel")" \
      | awk -F'\t' -v OFS='\t' '
          /^##/ { next }
          /^#/ { for (i = 1; i <= NF; i++) c[$i] = i; next }
          {
            chr = $c["#CHROM"]; sub(/^chr/, "", chr)
            if (chr !~ /^([1-9]|1[0-9]|2[0-2])$/) next
            r = $c["REF"]; a = $c["ALT"]
            if (r ~ /^[ACGT]$/ && a ~ /^[ACGT]$/) print "chr" chr, $c["POS"], r, a
          }' | LC_ALL=C sort -u -k1,1 -k2,2n -k3,3 -k4,4 > "${sites}.part" || [ ! -s "${sites}.part" ]; then
    rm -f "${sites}.part"
    echo "[WARN] Could not read the GRCh38 variants of ${panel}; run: $0 --ancestry-panel ${GENOME_DIR}"
    return 1
  fi
  mv "${sites}.part" "$sites"
  echo "[OK] ancestry panel ${name} (steps 25 and 26): ${panel} ($(wc -l < "$sites" | tr -d ' ') GRCh38 SNVs in ${sites##*/})"
}

case "$OPT_IN" in
  ancestry-panel) install_ancestry_panel; exit $? ;;
  cyrius) install_cyrius; exit $? ;;
  parascopy-data) install_parascopy_data; exit $? ;;
  kir-data) install_kir_data; exit $? ;;
esac

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
# Returns 1 at the first move that fails (set -e does not apply inside a
# function called as `f || ...`), keeping DIR.
clinvar_install() {
  local f
  for f in $CLINVAR_SET; do
    case "$f" in *.tbi) continue ;; esac
    if [ -f "${1}/${f}" ]; then mv -f "${1}/${f}" "${CLINVARDIR}/${f}" || return 1; fi
  done
  for f in $CLINVAR_SET; do
    case "$f" in *.tbi) ;; *) continue ;; esac
    if [ -f "${1}/${f}" ]; then mv -f "${1}/${f}" "${CLINVARDIR}/${f}" || return 1; fi
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
# step 06 builds it again from the new release. RELEASE and that copy go before
# the files move: a refresh that stops between two moves leaves no release date
# on a set that may mix two releases (validate-setup.sh then says the date is
# unknown) and fails, so the next refresh starts again.
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
  rm -f "${CLINVARDIR}/RELEASE" \
    "${CLINVARDIR}/clinvar_pathogenic_chr.norm.vcf.gz" "${CLINVARDIR}/clinvar_pathogenic_chr.norm.vcf.gz.tbi"
  clinvar_install "$new" || {
    echo "ERROR: could not move the new ClinVar files into ${CLINVARDIR}; they may mix two releases. Run setup.sh --refresh clinvar again." >&2
    return 1
  }
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
echo "=== Phase 1: Reference Genome (~0.9 GB download, ~3.2 GB unpacked) ==="

FASTA="$REF_FASTA"
FAI="${REF_FASTA}.fai"
REFDIR=$(dirname "$FASTA")
mkdir -p "$REFDIR"

if [ -f "$FASTA" ] && [ -f "$FAI" ]; then
  echo "[OK] Reference genome already downloaded."
else
  echo "Downloading GRCh38 reference genome..."
  echo "  Source: ${REF_FASTA_URL}"

  if [ ! -f "$FASTA" ]; then
    case "$REF_FASTA_URL" in
      *.gz) DL="${FASTA}.download.gz" ;;
      *) DL="$FASTA" ;;
    esac
    # shellcheck disable=SC2046  # the checksum arguments are dropped when REF_FASTA_MD5 is empty
    fetch "$REF_FASTA_URL" "$DL" $([ -z "$REF_FASTA_MD5" ] || printf 'md5 %s' "$REF_FASTA_MD5") || {
      echo "  The download is resumed when setup.sh runs again."
      exit 1
    }
    if [ "$DL" != "$FASTA" ]; then
      echo "  Unpacking $(basename "$REF_FASTA_URL")..."
      if gzip -dc "$DL" > "${FASTA}.tmp"; then
        mv -f "${FASTA}.tmp" "$FASTA"
        rm -f "$DL"
      else
        rm -f "${FASTA}.tmp" "$DL"
        echo "ERROR: could not unpack $(basename "$REF_FASTA_URL"); removed it, run setup.sh again." >&2
        exit 1
      fi
    fi
  fi

  if [ ! -f "$FAI" ]; then
    if [ -z "$REF_FAI_MD5" ] || ! fetch "$REF_FAI_URL" "$FAI" md5 "$REF_FAI_MD5"; then
      echo "  Generating index with samtools..."
      # The index goes next to the FASTA, so its directory is writable here.
      run_in --rw "$REFDIR" "$SAMTOOLS_IMAGE" \
        samtools faidx "$REF_FASTA_C"
    fi
  fi
  echo "[OK] Reference genome downloaded."
fi
if [ -f "$FAI" ]; then echo "[OK] Reference: ${FASTA} ($(grep -c . "$FAI") sequences)"; fi

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
# A failed download does not stop setup here either.
install_sample_qc_data || true
# Step 25: its scores, pgsc_calc and its plugin (warnings only).
install_prs

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
  echo "[SKIP] PCGR/CPSR ref data (~7 GB download) — needed for step 17 (cancer predisposition)"
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
  echo "[SKIP] VEP ${PCGR_VEP_CACHE_RELEASE} cache (~24 GB) — needed for step 17 (CPSR runs VEP ${PCGR_VEP_CACHE_RELEASE} inside the PCGR image)"
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
