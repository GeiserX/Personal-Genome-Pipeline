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
#        ./scripts/setup.sh --ancestry-fine-panel <genome_dir>
#        ./scripts/setup.sh --yleaf-data <genome_dir>
#        ./scripts/setup.sh --consumer-sites <genome_dir>
#        ./scripts/setup.sh --viral-data <genome_dir>
#                                            what an opt-in step needs, and nothing
#                                            else: Cyrius (step 21), Parascopy's
#                                            homology table and models plus its
#                                            image (step 35), IPD-KIR (step 08
#                                            with KIR=true), pgsc_calc's ancestry
#                                            panel (step 26, and percentiles in
#                                            step 25), the fine-grained ancestry
#                                            panel plus its images (step 26b),
#                                            Yleaf's marker tables plus its image
#                                            (step 37), the site table of the
#                                            raw-data export (step 40), the
#                                            viral genomes of step 47
#        ./scripts/setup.sh --annotsv-data <genome_dir>
#                                            AnnotSV's data for step 5 (the bundle
#                                            when missing, gnomAD-SV's common SVs
#                                            and the MANE Select list) and exit
#        ./scripts/setup.sh --vep-cache <genome_dir>
#                                            the VEP cache of VEP_CACHE_RELEASE
#                                            (VEP, vcfanno, clinical filter and
#                                            slivar need it), before the first run
#
# This script downloads everything needed to run the pipeline:
#   1. GRCh38 reference genome + index: NCBI's GRCh38 no-ALT analysis set
#      (~0.9 GB download, ~3.2 GB unpacked)
#   2. ClinVar database (~200 MB), with its release date in clinvar/RELEASE
#   3. All Docker images (~10-15 GB)
#   4. The reference's sequence dictionary and small pinned data files: Delly's
#      exclude map, GRCh38 chromosome bands, an IPD-IMGT/HLA release, the
#      GENCODE gene coordinates T1K needs (~350 MB) and GENCODE's exons for
#      the cDNA read check of steps 06 and 17, both from one ~80 MB download,
#      and somalier's sites and VerifyBamID2's marker panel for step 33 (~10 MB),
#      and gnomAD v3.1's mitochondrial sites for step 20 (~5 MB)
#   5. AnnotSV annotation data for step 5 (~5.3 GB download, ~20 GB unpacked),
#      plus gnomAD-SV's common SVs (~1.7 GB download, kept as a few MB) and the
#      MANE Select transcript list
#   6. Step 25: the PGS Catalog scores of assets/pgs_scores.tsv (~400 MB), a
#      checkout of pgsc_calc and the Nextflow plugin it needs, so step 25 and
#      the PRS process run without the network
#
# The VEP cache (~26 GB, --vep-cache) and PCGR ref data (~7 GB) are
# downloaded separately because they take a long time.
#
# Every download goes through fetch (scripts/lib/common.sh): it is written to
# <file>.part, checked, and only then renamed, so a file that exists is whole.
# A failed image pull is tried FETCH_TRIES times too, FETCH_WAIT apart, and
# Docker's own message (a rate limit, a missing tag) is printed as it comes.
#
# The reference can be replaced by another GRCh38 build: set REF_FASTA (where
# it is stored), REF_FASTA_URL (a .gz URL is unpacked) and REF_FASTA_MD5 (an
# md5, the URL of a checksum file that lists the download, or empty for no
# check), and REF_FAI_MD5 the same way for the .fai published beside it (empty
# builds the index with samtools). docs/00-reference-setup.md lists the builds
# that work and docs/realignment.md what a change of reference means for
# existing samples.
set -euo pipefail

# The usage text below names data releases, so read them before it prints.
# shellcheck source=../versions.env
. "$(dirname "$0")/../versions.env"

PULL_ONLY=false
REFRESH=""
SAMPLE_QC_ONLY=false
OPT_IN=""
case "${1:-}" in
  --cyrius|--parascopy-data|--kir-data|--ancestry-panel|--ancestry-fine-panel|--yleaf-data|--consumer-sites|--viral-data|--vep-cache|--annotsv-data)
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
  echo "       $0 --cyrius | --parascopy-data | --kir-data | --ancestry-panel | --ancestry-fine-panel | --yleaf-data | --consumer-sites | --viral-data | --vep-cache | --annotsv-data <genome_dir>"
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
  echo "                Install Parascopy's GRCh38 homology table and models (step 35, ~50 MB), pull its image and exit."
  echo "  --kir-data    Install the IPD-KIR ${KIR_DB_RELEASE:-} database (KIR=true in step 08, ~40 MB) and exit."
  echo "  --ancestry-panel"
  echo "                Install pgsc_calc's ancestry reference panel ${PGSC_PANEL:-} (~7 GB download): step 26, and"
  echo "                percentiles instead of raw scores in step 25; then exit."
  echo "  --ancestry-fine-panel"
  echo "                Install the reference of step 26b, gnomAD's HGDP + 1000 Genomes panel (~19 GB is read, a few GB"
  echo "                stored, hours on a home connection), and its images; then exit. Step 26b runs ADMIXTURE,"
  echo "                free for academic use, and RFMix, free for academic research use only: any other use,"
  echo "                a personal one included, needs a licence from Stanford."
  echo "  --yleaf-data  Install Yleaf ${YLEAF_DATA_VERSION:-}'s Y marker tables and tree and the current YFull tree (step 37, ~75 MB), pull its image and exit."
  echo "  --consumer-sites"
  echo "                Build the site table of the raw-data export (step 40): the 23andMe v5 site list joined with"
  echo "                dbSNP ${CONSUMER_DBSNP_BUILD:-}'s GRCh38 VCF (~30 GB download, deleted once the table is built); then exit."
  echo "  --viral-data  Install the human viral genomes and phiX of step 47 (NCBI RefSeq, VIRAL_REF_ACCESSIONS, ~2 MB) and exit."
  echo "  --annotsv-data"
  echo "                Install AnnotSV's data for step 5 and exit: the bundle when it is missing (~5.3 GB), then"
  echo "                gnomAD-SV ${GNOMAD_SV_VERSION:-}'s common SVs (~1.7 GB download) and MANE ${MANE_RELEASE:-}'s Select transcripts."
  echo "  --vep-cache   Install the VEP ${VEP_CACHE_RELEASE:-} cache (step 13 and the Nextflow VEP: ~26 GB download, ~30 GB"
  echo "                unpacked, both on disk during the install) and exit. Without it a run skips VEP, vcfanno,"
  echo "                the clinical filter and slivar."
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

# install_yleaf_data: the yleaf/data folder of Yleaf YLEAF_DATA_VERSION (step
# 37): its GRCh38 marker positions and haplogroup tree. The biocontainer of
# YLEAF_IMAGE installs the code only. GitHub's archive of the tag, checked
# against YLEAF_DATA_SHA256. Then the current Y tree beside it, in data/ytree
# (bin/y_tree_place.py): YFull's YTree YTREE_VERSION and the GRCh38 SNP tables
# of Yleaf YTREE_POSITIONS_VERSION, each at a pinned commit and sha256. A
# failed download of the current tree is a warning, not a failure.
install_yleaf_data() {
  local dest="${GENOME_DIR}/reference/yleaf-${YLEAF_DATA_VERSION}" tgz
  if [ -s "${dest}/data/hg38/new_positions.txt" ] && [ -s "${dest}/data/hg_prediction_tables/tree.json" ]; then
    echo "[OK] Yleaf ${YLEAF_DATA_VERSION} marker tables (step 37) already present."
  else
    tgz="${dest}.tar.gz"
    if ! fetch "https://github.com/genid/Yleaf/archive/refs/tags/${YLEAF_DATA_VERSION}.tar.gz" "$tgz" sha256 "$YLEAF_DATA_SHA256"; then
      echo "[WARN] Could not download Yleaf ${YLEAF_DATA_VERSION}."
      return 1
    fi
    rm -rf "${dest}.part"
    mkdir -p "${dest}.part"
    tar -xzf "$tgz" -C "${dest}.part" --strip-components 2 "Yleaf-${YLEAF_DATA_VERSION}/yleaf/data"
    rm -f "$tgz"
    if [ ! -s "${dest}.part/data/hg38/new_positions.txt" ]; then
      echo "[WARN] The Yleaf ${YLEAF_DATA_VERSION} archive has no yleaf/data/hg38/new_positions.txt."
      rm -rf "${dest}.part"
      return 1
    fi
    rm -rf "$dest"
    mv "${dest}.part" "$dest"
    echo "[OK] Yleaf ${YLEAF_DATA_VERSION} marker tables and tree (step 37): ${dest}/data"
  fi

  # version.txt names every pin (releases, commits, sha256), so a re-pin of
  # any of them, at the same tree release too, installs the tables again.
  local yt="${dest}/data/ytree" label="YFull YTree ${YTREE_VERSION} (CC BY 4.0)" want
  want="$(printf 'tree\t%s\npositions\tYleaf %s GRCh38 SNP tables (GPL-3.0)\npins\t%s %s %s %s %s' \
    "$label" "$YTREE_POSITIONS_VERSION" "$YTREE_COMMIT" "$YTREE_SHA256" \
    "$YTREE_POSITIONS_COMMIT" "$YTREE_FTDNA_SHA256" "$YTREE_YFULL_SHA256")"
  if [ -s "${yt}/version.txt" ] && [ "$(cat "${yt}/version.txt")" = "$want" ]; then
    echo "[OK] Current Y tree (${label}, step 37) already present."
    return 0
  fi
  rm -rf "${yt}.part"
  mkdir -p "${yt}.part"
  if ! fetch "https://raw.githubusercontent.com/YFullTeam/YTree/${YTREE_COMMIT}/current_tree.json" \
         "${yt}.part/current_tree.json" sha256 "$YTREE_SHA256" \
     || ! fetch "https://raw.githubusercontent.com/genid/Yleaf/${YTREE_POSITIONS_COMMIT}/yleaf/data/hg38/ftdna_positions_hg38.txt" \
         "${yt}.part/positions_ftdna_hg38.txt" sha256 "$YTREE_FTDNA_SHA256" \
     || ! fetch "https://raw.githubusercontent.com/genid/Yleaf/${YTREE_POSITIONS_COMMIT}/yleaf/data/hg38/new_positions.txt" \
         "${yt}.part/positions_yfull_hg38.txt" sha256 "$YTREE_YFULL_SHA256"; then
    # Yleaf's own data is in place, so step 37 still runs: it says the
    # current tree is not installed and writes Yleaf's result only. An older
    # tree from an earlier pin goes too, so no step places on tables that
    # versions.env no longer names.
    echo "[WARN] Could not download the current Y tree (${label}); step 37 will place the sample on Yleaf's tree only. Run this again to retry."
    rm -rf "${yt}.part" "$yt"
    return 0
  fi
  printf '%s\n' "$want" > "${yt}.part/version.txt"
  rm -rf "$yt"
  mv "${yt}.part" "$yt"
  echo "[OK] Current Y tree (${label}) and Yleaf ${YTREE_POSITIONS_VERSION}'s SNP positions (step 37): ${yt}"
}

# install_viral_data: the viral genomes of step 47, VIRAL_REF_ACCESSIONS in
# that order as one FASTA from NCBI's E-utilities, checked against
# VIRAL_REF_SHA256. minimap2 indexes a file this small in a second at each run.
install_viral_data() {
  local dest="${GENOME_DIR}/reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa"
  if [ -s "$dest" ] && [ "$(_digest sha256 "$dest")" = "$VIRAL_REF_SHA256" ]; then
    echo "[OK] Viral genomes ${VIRAL_REF_VERSION} (step 47) already present."
    return 0
  fi
  rm -f "$dest"
  if ! fetch "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=${VIRAL_REF_ACCESSIONS}" \
      "$dest" sha256 "$VIRAL_REF_SHA256"; then
    echo "[WARN] Could not install the viral genomes of step 47 (NCBI E-utilities)."
    return 1
  fi
  echo "[OK] Viral genomes ${VIRAL_REF_VERSION} (step 47, $(grep -c '^>' "$dest") sequences): ${dest}"
}

# install_consumer_sites: the site table of the raw-data export (step 40). The
# 23andMe v5 site list (GRCh37 positions and reference bases) of
# CONSUMER_SITES_COMMIT, checked against CONSUMER_SITES_SHA256, is joined by
# rsID with dbSNP CONSUMER_DBSNP_BUILD's GRCh38 VCF (about 30 GB, checked
# against CONSUMER_DBSNP_MD5) in bin/export_consumer.py build-sites, which
# takes the GRCh38 position from dbSNP instead of a liftover and lists every
# site it drops with the reason. grep first keeps the lines of the list's
# rsIDs, so Python reads about 600,000 lines instead of more than a billion.
# Each rsID is matched with the tabs around it, which is the whole ID column;
# never with -w, whose word check made GNU grep read 30 MB in over 20 minutes
# with 620,000 patterns.
# The dbSNP file is deleted once the table is built. The table built is then
# compared with CONSUMER_SITES_TABLE_SHA256: a difference is a warning.
install_consumer_sites() {
  local dir="${GENOME_DIR}/reference/consumer" tpl db out dropped
  out="${dir}/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
  dropped="${dir}/23andme_v5_dropped_dbsnp${CONSUMER_DBSNP_BUILD}.tsv"
  if [ -s "$out" ]; then
    consumer_sites_pin_check "$out" "already present"
    return 0
  fi
  tpl="${dir}/23andme_v5_hg19_ref.txt.gz"
  if ! fetch "https://raw.githubusercontent.com/arrogantrobot/23andme2vcf/${CONSUMER_SITES_COMMIT}/23andme_v5_hg19_ref.txt.gz" \
       "$tpl" sha256 "$CONSUMER_SITES_SHA256"; then
    echo "[WARN] Could not download the 23andMe v5 site list."
    return 1
  fi
  db="${dir}/dbsnp_${CONSUMER_DBSNP_BUILD}_GRCh38.vcf.gz"
  # A run after a failed build reuses the download, once its md5 matches the
  # pin again: a file under that name that is not the pinned one is removed
  # and downloaded anew, so the table is never built from another file.
  if [ -s "$db" ]; then
    echo "Checking the md5 of the dbSNP ${CONSUMER_DBSNP_BUILD} file already downloaded: ${db}"
    if [ "$(_digest md5 "$db")" != "$CONSUMER_DBSNP_MD5" ]; then
      echo "[WARN] ${db} is not the pinned dbSNP ${CONSUMER_DBSNP_BUILD} (md5 differs): removing it and downloading it again."
      rm -f "$db"
    fi
  fi
  if [ -s "$db" ]; then
    echo "dbSNP ${CONSUMER_DBSNP_BUILD} for GRCh38 already downloaded and checked: ${db}"
  else
    echo "Downloading dbSNP ${CONSUMER_DBSNP_BUILD} for GRCh38 (about 30 GB; deleted once the table is built)..."
    if ! fetch "https://ftp.ncbi.nih.gov/snp/archive/${CONSUMER_DBSNP_BUILD}/VCF/GCF_000001405.40.gz" "$db" md5 "$CONSUMER_DBSNP_MD5"; then
      echo "[WARN] Could not download dbSNP ${CONSUMER_DBSNP_BUILD}. Run this command again: the download resumes."
      return 1
    fi
  fi
  echo "Building the site table (reads the whole dbSNP file once)..."
  # The table is written to reference/consumer/, outside any sample directory.
  # shellcheck disable=SC2016  # $1 to $5 expand in the container's shell
  if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 2 --memory 4g "$PYTHON_IMAGE" bash -euo pipefail -c '
      gzip -dc "$1" | cut -f3 | grep "^rs" | sed "s/.*/\t&\t/" > "$3.ids"
      gzip -dc "$2" | LC_ALL=C grep -Ff "$3.ids" \
        | python3 /pgp-bin/export_consumer.py build-sites --template "$1" --dbsnp - --out "$3" --dropped "$4" --source "$5"
      rm -f "$3.ids"' _ "$(cpath "$tpl")" "$(cpath "$db")" "$(cpath "$out")" "$(cpath "$dropped")" \
      "23andme2vcf ${CONSUMER_SITES_COMMIT} 23andme_v5_hg19_ref.txt.gz; dbSNP ${CONSUMER_DBSNP_BUILD} GCF_000001405.40"; then
    rm -f "$out" "${out}.ids" "${out}.tmp" "${dropped}.tmp"
    echo "[WARN] Could not build the site table of step 40; the dbSNP download is kept: run this command again."
    return 1
  fi
  rm -f "$db"
  consumer_sites_pin_check "$out" "built (dropped sites and their reasons: ${dropped})"
}

# consumer_sites_pin_check TABLE WHAT: compare the site table with
# CONSUMER_SITES_TABLE_SHA256 and say so. A difference is a warning: both
# inputs are checked against their own pins, so it means build-sites changed
# since the pin was taken, and step 40 can still use the table.
consumer_sites_pin_check() {
  local got
  got=$(consumer_sites_digest "$1")
  if [ "$got" = "$CONSUMER_SITES_TABLE_SHA256" ]; then
    echo "[OK] Site table of step 40 $2: $1, the same as the pinned build"
  else
    echo "[WARN] Site table of step 40 $2: $1 differs from the pinned build (content sha256 ${got}, pinned ${CONSUMER_SITES_TABLE_SHA256})."
    echo "       It was built by another version of bin/export_consumer.py build-sites. Step 40 uses it; to rebuild it, remove it and run this command again."
  fi
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
# against the md5 that folder lists, and beside it the panel's common GRCh38
# biallelic SNVs (chrN, position, REF, ALT), from the panel's own genotypes:
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
  [ -s "$panel" ] || echo "Downloading pgsc_calc's ancestry panel ${name} (the 1000 Genomes panel is about 7 GB; an interrupted download resumes)..."
  # fetch stores a download only once its md5 matches, so a panel that is there is whole.
  if [ ! -s "$panel" ] && ! fetch "${PGSC_RESOURCES}/${name}.tar.zst" "$panel" md5 "${PGSC_RESOURCES}/md5s.txt"; then
    echo "[WARN] Could not download the ancestry panel; run: $0 --ancestry-panel ${GENOME_DIR}"
    return 1
  fi
  # Only the SNVs pgsc_calc can use for the projection: its FILTER_VARIANTS
  # keeps biallelic ACGT SNVs on the autosomes with a panel frequency of 5%
  # or more (its maf_ref). The whole GRCh38 table holds about 62 million
  # SNVs; genotyping all of them from a gVCF would take step 25 hours.
  echo "  Listing the panel's common GRCh38 SNVs (plink2 on the panel's genotypes)..."
  local tmp="${dir}/.${name}.sites.tmp" prefix
  rm -rf "$tmp"
  mkdir -p "$tmp"
  # shellcheck disable=SC2016  # $1 and $2 belong to the inner sh
  if ! run_in --rw "$dir" "$PGSC_ZSTD_IMAGE" sh -c 'cd "$2" && tar -xf "$1" --wildcards "GRCh38_*_ALL.pgen" "GRCh38_*_ALL.psam" "GRCh38_*_ALL.pvar.zst"' \
      _ "$(cpath "$panel")" "$(cpath "$tmp")"; then
    rm -rf "$tmp"
    echo "[WARN] Could not unpack the GRCh38 genotypes of ${panel}; run: $0 --ancestry-panel ${GENOME_DIR}"
    return 1
  fi
  prefix=$(find "$tmp" -maxdepth 1 -name 'GRCh38_*_ALL.pgen' | head -1)
  prefix=${prefix%.pgen}
  # One chromosome at a time: plink2 then loads a 22nd of the panel's variant
  # table (the whole table, 61.6 million variants, did not fit a 16 GB runner).
  # --allow-extra-chr: the panel names contigs plink2 does not know, as
  # pgsc_calc's own plink2 calls on it allow.
  # plink2 still indexes the whole table: 4 GB of workspace ran out of memory
  # even for chromosome 1. 8 GB, the most pgsc_calc's own call on the panel
  # gets here, or three quarters of a smaller machine's RAM.
  local c ok=true mem_mb
  if [ -r /proc/meminfo ]; then
    mem_mb=$(awk '/^MemTotal:/ {print int($2 * 3 / 4 / 1024)}' /proc/meminfo)
  else
    mem_mb=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 17179869184) * 3 / 4 / 1048576 ))
  fi
  [ "$mem_mb" -le 8000 ] || mem_mb=8000
  [ -n "$prefix" ] || ok=false
  for c in $(seq 1 22); do
    $ok || break
    # A panel without a common SNV on a chromosome ("No variants remaining") is not an error.
    if ! run_in --rw "$dir" "$PLINK2_IMAGE" plink2 --pfile "$(cpath "$prefix")" vzs --chr "$c" \
          --allow-extra-chr --snps-only just-acgt --max-alleles 2 --maf 0.05 \
          --make-just-pvar --threads "$THREADS" --memory "$mem_mb" --out "$(cpath "${tmp}/common_${c}")" > "${tmp}/plink2.log" 2>&1 \
       && ! grep -q 'No variants remaining' "${tmp}/plink2.log"; then
      head -n 40 "${tmp}/plink2.log"
      ok=false
    fi
  done
  if ! $ok || ! awk -F'\t' -v OFS='\t' '
          /^##/ { next }
          /^#/ { for (i = 1; i <= NF; i++) c[$i] = i; next }
          {
            chr = $c["#CHROM"]; sub(/^chr/, "", chr)
            r = $c["REF"]; a = $c["ALT"]
            if (r ~ /^[ACGT]$/ && a ~ /^[ACGT]$/) print "chr" chr, $c["POS"], r, a
          }' "${tmp}"/common_*.pvar | LC_ALL=C sort -u -k1,1 -k2,2n -k3,3 -k4,4 > "${sites}.part" || [ ! -s "${sites}.part" ]; then
    rm -rf "$tmp" "${sites}.part"
    echo "[WARN] Could not list the common GRCh38 SNVs of ${panel}; run: $0 --ancestry-panel ${GENOME_DIR}"
    return 1
  fi
  rm -rf "$tmp"
  mv "${sites}.part" "$sites"
  echo "[OK] ancestry panel ${name} (steps 25 and 26): ${panel} ($(wc -l < "$sites" | tr -d ' ') common GRCh38 SNVs in ${sites##*/})"
}

# --- Step 26b: the fine-grained ancestry panel (opt-in) -----------------------------
# gnomAD's HGDP + 1000 Genomes callset, phased with SHAPEIT5 (Koenig et al.
# 2024): 4,091 people of 80 populations, open access. setup.sh streams the
# common biallelic SNVs of chr1 to chr22 from gnomAD's public bucket (about
# 19 GB is read; only the kept SNVs are stored) and keeps one per 2 kb.
# One of the 22 files, named in full, so scripts/ci/check-links.sh can request
# it; the others differ only in the chromosome.
ANCESTRY_FINE_CHR22_URL=https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.bcf
ANCESTRY_FINE_URL=${ANCESTRY_FINE_CHR22_URL%/*}
ANCESTRY_FINE_SUFFIX=${ANCESTRY_FINE_CHR22_URL#*hgdp1kgp_chr22}
ANCESTRY_FINE_META_URL=https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv
# Beagle's GRCh38 PLINK genetic map (HapMap, lifted from GRCh37), for phasing and painting.
# Beagle's server is its only host and was unreachable for hours on 2026-10-15;
# the Internet Archive's copy is byte for byte the same file (same sha256), so
# the install falls back to it, checked against the same pin.
ANCESTRY_FINE_MAP_URL=https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
ANCESTRY_FINE_MAP_MIRROR_URL=https://web.archive.org/web/20260306143720id_/https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip

# ancestry_fine_pins: one line per chromosome file of the panel: chrN, the
# crc32c Google reports for it and its size, from a HEAD request each. Their
# sha256 is ANCESTRY_FINE_PANEL_SHA256, so a file the bucket replaces changes
# it. bcftools reads the panel a region at a time, so the bytes cannot be
# checked as they stream; the bucket's own checksum of each file can.
ancestry_fine_pins() {
  local c
  command -v curl >/dev/null 2>&1 || { echo "ERROR: --ancestry-fine-panel needs curl to check the panel's files" >&2; return 1; }
  for c in $(seq 1 22); do
    curl -fsI --connect-timeout 20 --max-time 60 --retry 3 --retry-delay 5 \
        "${ANCESTRY_FINE_URL}/hgdp1kgp_chr${c}${ANCESTRY_FINE_SUFFIX}" | tr -d '\r' \
      | awk -v c="chr${c}" 'tolower($1) == "x-goog-hash:" && $2 ~ /^crc32c=/ {h = substr($2, 8)}
                            tolower($1) == "content-length:" {n = $2}
                            END {if (h == "" || n == "") exit 1; print c "\t" h "\t" n}' \
      || { echo "ERROR: no crc32c or size for the panel's chr${c} file (${ANCESTRY_FINE_URL})" >&2; return 1; }
  done
}

# install_ancestry_fine_panel: the reference of step 26b, in
# reference/ancestry_fine/ANCESTRY_FINE_PANEL (ANCESTRY_FINE_DIR: another
# folder inside GENOME_DIR). ADMIXTURE is free for academic use; RFMix for
# academic research use only, and any other use, a personal one included,
# needs a licence from Stanford. Four parts, each redone only when its inputs change:
#   samples.tsv  every person of the panel with population, continental and
#                regional group, reference and painting flags (panel-map of
#                bin/ancestry_fine_report.py); ANCESTRY_FINE_SAMPLE_MAP gives
#                a table of your own instead;
#   panel/       chrN.vcf.gz: their phased genotypes at the biallelic SNVs with a
#                minor allele frequency of 5% or more, one per 2 kb (bcftools
#                +prune); sites.tsv (chrN, position, REF, ALT) and positions.tsv;
#   maps/        chrN.map (Beagle) and chrN.rfmix.map (RFMix's column order);
#   model/       ref.bed of the reference individuals on the SNVs left after LD
#                pruning (r2 above 0.1 within 500 kb removed); for each level
#                (continent, region) the group allele frequencies of a supervised
#                ADMIXTURE run, <level>.P, with their column order in levels.tsv;
#                and paint_map.tsv, the reference individuals RFMix paints with.
# ANCESTRY_FINE_REGIONS (a BED file) keeps only those regions, and
# ANCESTRY_FINE_EXCLUDE (sample ids) leaves people out of the panel: the e2e
# case uses both, for the fixture's slices and the two people of its
# synthetic mix. ANCESTRY_FINE_THIN (default 2kb) is the window that keeps
# one SNV (the e2e case keeps more of its few megabases with 200bp), and
# ANCESTRY_FINE_PAINT_PER_GROUP (default 50, at least 10) caps the painting
# references of each continental group.
install_ancestry_fine_panel() {
  local dir=${ANCESTRY_FINE_DIR:-${GENOME_DIR}/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}}
  local regions=${ANCESTRY_FINE_REGIONS:-} exclude=${ANCESTRY_FINE_EXCLUDE:-} smap=${ANCESTRY_FINE_SAMPLE_MAP:-}
  local thin=${ANCESTRY_FINE_THIN:-2kb}
  local v c e level col k order rstamp="" pstamp mstamp tmp chroms=""
  cpath "$dir" >/dev/null || return 1
  echo "Installing the fine-grained ancestry panel ${ANCESTRY_FINE_PANEL} (step 26b) in ${dir}"
  echo "  Step 26b runs ADMIXTURE, free for academic use, and RFMix, free for academic research use only:"
  echo "  any other use of RFMix, a personal one included, needs a licence from Stanford."
  for v in BCFTOOLS_HTTPS_IMAGE ADMIXTURE_IMAGE BEAGLE_IMAGE RFMIX_IMAGE; do ensure_image "$v" 26b || return 1; done
  mkdir -p "${dir}/panel" "${dir}/maps" "${dir}/model"
  echo "  Checking the panel's files against ANCESTRY_FINE_PANEL_SHA256..."
  if ! ancestry_fine_pins > "${dir}/pins.tsv"; then
    echo "[WARN] Could not read the panel's files at ${ANCESTRY_FINE_URL}; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
    return 1
  fi
  if [ "$(_digest sha256 "${dir}/pins.tsv")" != "$ANCESTRY_FINE_PANEL_SHA256" ]; then
    echo "[WARN] The panel's files at ${ANCESTRY_FINE_URL} are not the ones ANCESTRY_FINE_PANEL_SHA256 pins (see ${dir}/pins.tsv)."
    echo "       The bucket replaced a file: read what changed before moving the pin in versions.env."
    return 1
  fi
  # The metadata and the map, unless they are there with their pinned sha256;
  # the map's mirror (the fourth field) is tried when its own host fails.
  for v in "${ANCESTRY_FINE_META_URL}|${dir}/meta.tsv|${ANCESTRY_FINE_META_SHA256}|" \
           "${ANCESTRY_FINE_MAP_URL}|${dir}/maps/plink.GRCh38.map.zip|${ANCESTRY_FINE_MAP_SHA256}|${ANCESTRY_FINE_MAP_MIRROR_URL}"; do
    IFS='|' read -r c e k m <<< "$v"
    [ -s "$e" ] && [ "$(_digest sha256 "$e")" = "$k" ] && continue
    fetch "$c" "$e" sha256 "$k" && continue
    if [ -n "$m" ]; then
      rm -f "${e}.part"
      echo "  Trying the mirror of ${c##*/}: ${m}"
      fetch "$m" "$e" sha256 "$k" && continue
    fi
    echo "[WARN] Could not download ${c}; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
    return 1
  done
  if [ -n "$regions" ]; then
    [ -s "$regions" ] || { echo "[WARN] ANCESTRY_FINE_REGIONS=${regions} is not a file"; return 1; }
    cp "$regions" "${dir}/regions.bed"
    rstamp=$(_digest sha256 "${dir}/regions.bed")
    chroms=$(awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {sub(/^chr/, "", $1); print $1}' "${dir}/regions.bed" | sort -un | sed 's/^/chr/' | tr '\n' ' ')
  else
    rm -f "${dir}/regions.bed"
    for c in $(seq 1 22); do chroms+="chr${c} "; done
  fi

  # samples.tsv: from the metadata, or the table given.
  # --net: the sample list is the header of a panel file in gnomAD's bucket.
  # shellcheck disable=SC2016  # $1 belongs to the inner sh
  if ! run_in --net --rw "$dir" "$BCFTOOLS_HTTPS_IMAGE" sh -c 'bcftools query -l "$1" 2>/dev/null' \
        _ "$ANCESTRY_FINE_CHR22_URL" > "${dir}/bcf_samples.txt.part" \
     || [ ! -s "${dir}/bcf_samples.txt.part" ]; then
    rm -f "${dir}/bcf_samples.txt.part"
    echo "[WARN] Could not read the panel's sample list; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
    return 1
  fi
  mv "${dir}/bcf_samples.txt.part" "${dir}/bcf_samples.txt"
  if [ -n "$smap" ]; then
    head -n 1 "$smap" | grep -q $'^sample\tpopulation\tcontinent\tregion\treference\tpaint$' \
      || { echo "[WARN] ANCESTRY_FINE_SAMPLE_MAP=${smap} does not start with the header of samples.tsv"; return 1; }
    cp "$smap" "${dir}/samples.tsv.part"
    if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
          python3 /pgp-bin/ancestry_fine_report.py check-map --map "$(cpath "${dir}/samples.tsv.part")"; then
      rm -f "${dir}/samples.tsv.part"
      echo "[WARN] ANCESTRY_FINE_SAMPLE_MAP=${smap} has a group too small to model"
      return 1
    fi
  elif ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
        python3 /pgp-bin/ancestry_fine_report.py panel-map --meta "$(cpath "${dir}/meta.tsv")" \
          --bcf-samples "$(cpath "${dir}/bcf_samples.txt")" --exclude "$exclude" \
          --paint-per-group "${ANCESTRY_FINE_PAINT_PER_GROUP:-50}" --out "$(cpath "${dir}/samples.tsv.part")"; then
    rm -f "${dir}/samples.tsv.part"
    echo "[WARN] Could not write the panel's sample table"
    return 1
  fi
  mv "${dir}/samples.tsv.part" "${dir}/samples.tsv"
  awk -F'\t' 'NR > 1 {print $1}' "${dir}/samples.tsv" > "${dir}/panel/keep.txt"

  # panel/: redone when the pins, the regions or the people kept change; a
  # chromosome done by an interrupted run with the same inputs is kept.
  pstamp="pins=${ANCESTRY_FINE_PANEL_SHA256} regions=${rstamp:-all} thin=${thin} keep=$(_digest sha256 "${dir}/panel/keep.txt") bcftools=${BCFTOOLS_HTTPS_IMAGE}"
  if [ "$(cat "${dir}/panel/DONE" 2>/dev/null)" != "$pstamp" ]; then
    rm -f "${dir}/panel/DONE" "${dir}/model/DONE"
    for c in $chroms; do
      [ ! -s "${dir}/panel/${c}.vcf.gz.tbi" ] || [ "$(cat "${dir}/panel/${c}.stamp" 2>/dev/null)" != "$pstamp" ] || continue
      echo "  ${c}: reading the panel's common SNVs..."
      # --net: the panel is read from gnomAD's bucket over https. htslib
      # stores the remote index in the working directory, so the container
      # works in panel/ and deletes it after. bcftools warns about the files'
      # PP header line on every read; those lines go.
      # shellcheck disable=SC2016  # $1 to $5 belong to the inner bash
      if ! run_in --net --rw "$dir" "$BCFTOOLS_HTTPS_IMAGE" bash -euo pipefail -c '
            cd "$(dirname "$3")"; trap "rm -f ./*.bcf.csi" EXIT
            r=(); [ -z "$4" ] || r=(-R "$4" --regions-overlap pos)
            { bcftools view "${r[@]}" --force-samples -S "$2" -m2 -M2 -v snps -Ou "$1" \
                | bcftools view -q 0.05:minor -Ou \
                | bcftools annotate -x ID,INFO,^FORMAT/GT -Ou \
                | bcftools +prune -n 1 -w "$5" -N maxAF -Oz -o "$3.part"; } 2> >(grep -v bcf_hdr_check_sanity >&2)
            bcftools index -t "$3.part"
            mv "$3.part.tbi" "$3.tbi"; mv "$3.part" "$3"' \
            _ "${ANCESTRY_FINE_URL}/hgdp1kgp_${c}${ANCESTRY_FINE_SUFFIX}" \
            "$(cpath "${dir}/panel/keep.txt")" "$(cpath "${dir}/panel/${c}.vcf.gz")" \
            "$([ -z "$rstamp" ] || cpath "${dir}/regions.bed")" "$thin"; then
        echo "[WARN] Could not read ${c} of the panel; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
        return 1
      fi
      echo "$pstamp" > "${dir}/panel/${c}.stamp"
    done
    # sites.tsv: chrN, position, REF and ALT of every SNV; positions.tsv: the
    # first two columns, the form bcftools -R and -T read.
    # shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
    if ! run_in --rw "$dir" "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
          'cd "$1"; for c in $2; do bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" "$c.vcf.gz"; done' \
          _ "$(cpath "${dir}/panel")" "$chroms" > "${dir}/panel/sites.tsv.part" || [ ! -s "${dir}/panel/sites.tsv.part" ]; then
      echo "[WARN] Could not list the panel's sites"
      return 1
    fi
    cut -f1,2 "${dir}/panel/sites.tsv.part" > "${dir}/panel/positions.tsv"
    mv "${dir}/panel/sites.tsv.part" "${dir}/panel/sites.tsv"
    echo "$pstamp" > "${dir}/panel/DONE"
  fi
  echo "  panel: $(wc -l < "${dir}/panel/sites.tsv" | tr -d ' ') SNVs of $(wc -l < "${dir}/panel/keep.txt" | tr -d ' ') people"

  # maps/: Beagle's file of each chromosome, and the same in RFMix's column
  # order (chromosome, position, cM).
  for c in $chroms; do
    [ -s "${dir}/maps/${c}.rfmix.map" ] && continue
    # shellcheck disable=SC2016  # no shell variables: the python reads sys.argv
    if ! run_in --rw "$dir" "$PYTHON_IMAGE" python3 -c '
import sys, zipfile
z, c, out = sys.argv[1:]
with zipfile.ZipFile(z) as f:
    text = f.read(f"chr_in_chrom_field/plink.chr{c}.GRCh38.map").decode()
rows = [line.split() for line in text.splitlines() if line.strip()]
open(f"{out}/{c}.map", "w").write(text)
open(f"{out}/{c}.rfmix.map", "w").write("".join(f"{r[0]}\t{r[3]}\t{r[2]}\n" for r in rows))
' "$(cpath "${dir}/maps/plink.GRCh38.map.zip")" "$c" "$(cpath "${dir}/maps")"; then
      echo "[WARN] Could not unpack the genetic map of ${c}"
      return 1
    fi
  done

  # model/: redone when samples.tsv or the panel changes.
  mstamp="panel=$(_digest sha256 "${dir}/panel/sites.tsv") samples=$(_digest sha256 "${dir}/samples.tsv") admixture=${ADMIXTURE_IMAGE}"
  if [ "$(cat "${dir}/model/DONE" 2>/dev/null)" != "$mstamp" ]; then
    tmp="${dir}/model/tmp"
    rm -rf "$tmp" "${dir}/model/DONE"
    mkdir -p "$tmp"
    awk -F'\t' 'NR > 1 && $5 == 1 {print $1}' "${dir}/samples.tsv" > "${tmp}/ref.txt"
    echo "  model: LD pruning on $(wc -l < "${tmp}/ref.txt" | tr -d ' ') reference individuals, then ADMIXTURE's group frequencies..."
    # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
    if ! run_in --rw "$dir" "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
          for c in $4; do bcftools view -S "$2" -Ou "$1/$c.vcf.gz" | bcftools +prune -m 0.1 -w 500kb -Ob -o "$3/$c.bcf"; done
          for c in $4; do echo "$3/$c.bcf"; done > "$3/bcfs.txt"
          bcftools concat -f "$3/bcfs.txt" -Oz -o "$3/adm.vcf.gz"' \
          _ "$(cpath "${dir}/panel")" "$(cpath "${tmp}/ref.txt")" "$(cpath "$tmp")" "$chroms" \
       || ! run_in --rw "$dir" "$PLINK2_IMAGE" plink2 --vcf "$(cpath "${tmp}/adm.vcf.gz")" --make-bed \
              --set-all-var-ids '@:#' --threads "$THREADS" --out "$(cpath "${tmp}/ref")" > "${tmp}/plink2.log" 2>&1; then
      tail -n 20 "${tmp}/plink2.log" 2>/dev/null || true
      echo "[WARN] Could not make the reference's ADMIXTURE input"
      return 1
    fi
    : > "${tmp}/levels.tsv"
    for v in continent:3 region:4; do
      level=${v%:*} col=${v#*:}
      awk -F'\t' -v c="$col" 'NR == FNR {if (FNR > 1) g[$1] = $c; next} {print g[$2]}' "${dir}/samples.tsv" "${tmp}/ref.fam" > "${tmp}/ref_${level}.pop"
      k=$(sort -u "${tmp}/ref_${level}.pop" | grep -cv '^-\{0,1\}$' || true)
      for e in bed bim fam; do ln -sf "ref.${e}" "${tmp}/ref_${level}.${e}"; done
      # shellcheck disable=SC2016  # $1 to $4 belong to the inner sh
      if ! run_in --rw "$dir" "$ADMIXTURE_IMAGE" sh -c 'cd "$1" && admixture --supervised -j"$4" "ref_$2.bed" "$3" > "admixture_$2.log" 2>&1' \
            _ "$(cpath "$tmp")" "$level" "$k" "$THREADS" \
         || ! order=$(run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" python3 /pgp-bin/ancestry_fine_report.py order \
              --pop "$(cpath "${tmp}/ref_${level}.pop")" --q "$(cpath "${tmp}/ref_${level}.${k}.Q")"); then
        tail -n 20 "${tmp}/admixture_${level}.log" 2>/dev/null || true
        echo "[WARN] ADMIXTURE could not estimate the frequencies of the ${level} groups"
        return 1
      fi
      mv "${tmp}/ref_${level}.${k}.P" "${tmp}/${level}.P"
      printf '%s\t%s\t%s\n' "$level" "$k" "$order" >> "${tmp}/levels.tsv"
    done
    awk -F'\t' -v OFS='\t' 'NR > 1 && $6 == 1 {print $1, $3}' "${dir}/samples.tsv" > "${tmp}/paint_map.tsv"
    for e in ref.bed ref.bim ref.fam continent.P region.P levels.tsv paint_map.tsv; do mv -f "${tmp}/${e}" "${dir}/model/${e}"; done
    rm -rf "$tmp"
    echo "$mstamp" > "${dir}/model/DONE"
  fi
  echo "[OK] fine-grained ancestry panel ${ANCESTRY_FINE_PANEL} (step 26b): ${dir} ($(wc -l < "${dir}/model/ref.bim" | tr -d ' ') SNVs for ADMIXTURE, $(wc -l < "${dir}/model/paint_map.tsv" | tr -d ' ') painting references)"
}

# pull_image IMAGE: docker pull, tried FETCH_TRIES times (default 3),
# FETCH_WAIT seconds apart (default 5, then 10, ...), the knobs fetch uses.
# Docker's own message goes to the terminal: a rate limit (toomanyrequests)
# or a registry error passes with a retry, a wrong tag (manifest unknown)
# does not.
pull_image() {
  local img=$1 tries=${FETCH_TRIES:-3} wait=${FETCH_WAIT:-} i
  case "$tries" in ''|*[!0-9]*) tries=0 ;; *) tries=$((10#$tries)) ;; esac
  [ "$tries" -gt 0 ] || { echo "ERROR: FETCH_TRIES must be a whole number above 0, got '${FETCH_TRIES}'" >&2; return 1; }
  case "$wait" in *[!0-9]*) echo "ERROR: FETCH_WAIT must be a whole number of seconds, got '${wait}'" >&2; return 1 ;; esac
  for ((i = 1; i <= tries; i++)); do
    "$CONTAINER_ENGINE" pull "$img" && return 0
    echo "  Pull attempt ${i}/${tries} failed: ${img} (Docker's message is above)" >&2
    [ "$i" -lt "$tries" ] && sleep "${wait:-$((i * 5))}"
  done
  return 1
}

# keep_image IMAGE: tag a digest-pinned image (tag_pinned_image in
# scripts/lib/common.sh), so a dangling-image cleanup does not delete it.
# When the tag fails, what Docker lists decides the warning, as in
# validate-setup.sh (image_is_dangling): on the classic image store an image
# with no tag is dangling and the cleanup deletes it; on the containerd image
# store, or with another local tag, it is not; or Docker did not say.
keep_image() {
  local rc=0 msg
  tag_pinned_image "$1" && return 0
  image_is_dangling "$1" || rc=$?
  case $rc in
    0) msg="Docker lists it as dangling, so a dangling-image cleanup would delete it" ;;
    1) msg="it has no pgp-pinned tag, but Docker does not list it as dangling, so a dangling-image cleanup keeps it" ;;
    *) msg="Docker did not list its dangling images, so a dangling-image cleanup may delete it" ;;
  esac
  echo "  WARNING: Could not tag ${1} as $(pinned_tag "$1"); ${msg}." >&2
  return 1
}

# --- Step 5: AnnotSV's annotation data, gnomAD-SV and MANE ----------------------------
# The AnnotSV image holds code only; without its data AnnotSV exits with an
# error. Two files are added to the bundle, each pinned in versions.env:
#   - gnomAD-SV GNOMAD_SV_VERSION's sites (GNOMAD_SV_SHA256): its common
#     deletions, duplications, insertions and inversions, written where AnnotSV
#     reads a gnomAD-SV v4 source of its benign set (ACMG 2F/4O, B_*_AFmax);
#   - MANE MANE_RELEASE's summary (MANE_SUMMARY_SHA256): the MANE Select RefSeq
#     transcripts, the list step 05 passes as -txFile, so the Tx, Frameshift,
#     Location and exon columns follow MANE Select.
ANNOTSV_DIR="${GENOME_DIR:+${GENOME_DIR}/annotsv_annotations}"
ANNOTSV_NAME="Annotations_Human_${ANNOTSV_ANNOTATIONS_VERSION}.tar.gz"
ANNOTSV_URL="https://www.lbgi.fr/~geoffroy/Annotations/${ANNOTSV_NAME}"
ANNOTSV_TARBALL="${GENOME_DIR:+${GENOME_DIR}/${ANNOTSV_NAME}}"
GNOMAD_SV_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/${GNOMAD_SV_VERSION}/genome_sv/gnomad.v${GNOMAD_SV_VERSION}.sv.sites.vcf.gz"
MANE_SUMMARY_URL="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANE_RELEASE}/MANE.GRCh38.v${MANE_RELEASE}.summary.txt.gz"

# install_annotsv_bundle: AnnotSV's annotation bundle, downloaded under a .part
# name, extracted into a temporary directory and moved into place only when
# complete.
install_annotsv_bundle() {
  if [ -d "${ANNOTSV_DIR}/Annotations_Human/Genes/GRCh38" ]; then
    echo "[OK] AnnotSV annotation data already present."
    return 0
  fi
  echo "Downloading AnnotSV ${ANNOTSV_ANNOTATIONS_VERSION} annotation data (~5.3 GB)..."
  echo "  The AnnotSV server is slow (about 0.8 MB/s measured from a GitHub runner), so this can take"
  echo "  1-2 hours. An interrupted download resumes when setup.sh runs again."
  if ! fetch "$ANNOTSV_URL" "$ANNOTSV_TARBALL"; then
    echo "[WARN] AnnotSV annotation download failed. Step 5 (AnnotSV) will be skipped until it is in place."
    echo "  Re-run setup.sh, or download it by hand:"
    echo "    curl -fL -C - -o ${ANNOTSV_TARBALL} ${ANNOTSV_URL}"
    echo "    mkdir -p ${ANNOTSV_DIR} && tar -xzf ${ANNOTSV_TARBALL} -C ${ANNOTSV_DIR}"
    return 1
  fi
  echo "  Extracting..."
  rm -rf "${ANNOTSV_DIR}.part"
  mkdir -p "${ANNOTSV_DIR}.part"
  if ! tar -xzf "$ANNOTSV_TARBALL" -C "${ANNOTSV_DIR}.part"; then
    rm -rf "${ANNOTSV_DIR}.part"
    echo "[WARN] Could not extract ${ANNOTSV_TARBALL}. Step 5 (AnnotSV) will be skipped until it is in place."
    return 1
  fi
  if [ ! -d "${ANNOTSV_DIR}.part/Annotations_Human/Genes/GRCh38" ]; then
    rm -rf "${ANNOTSV_DIR}.part" "$ANNOTSV_TARBALL"
    echo "[WARN] ${ANNOTSV_TARBALL} holds no Annotations_Human/Genes/GRCh38; removed it. Step 5 (AnnotSV) will be skipped until it is in place."
    return 1
  fi
  rm -rf "$ANNOTSV_DIR"
  mv "${ANNOTSV_DIR}.part" "$ANNOTSV_DIR"
  rm -f "$ANNOTSV_TARBALL"
  echo "[OK] AnnotSV annotation data: ${ANNOTSV_DIR}/ ($(du -sh "$ANNOTSV_DIR" | cut -f1))"
}

# install_annotsv_gnomad_sv: gnomAD-SV's sites, reduced to the records AnnotSV
# keeps for its benign set (FILTER PASS; DEL, DUP, INS or INV; AN >= 1000,
# AF >= 0.001 and N_HOMALT >= 5; AnnotSV then applies its -benignAF, 1%, when
# it annotates). AnnotSV 3.5 reads gnomAD-SV v4 from files named
# gnomad.v4.0.sv.chr*.vcf.gz in BenignSV/GRCh38, and from the INFO keys END,
# SVTYPE, controls_and_biobanks_AN, _AF and _N_HOMALT in that order. The v4.1
# sites file has no controls-and-biobanks subset, so its AN, AF and N_HOMALT of
# all genomes are written under those keys, in a file of that name. One AnnotSV
# run (annotsv_merge_benign) merges the file into its benign_*_SV_GRCh38.sorted.bed
# and deletes it; only then does gnomad-sv.installed record what was added.
install_annotsv_gnomad_sv() {
  local benign="${ANNOTSV_DIR}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38"
  local stamp="${benign}/gnomad-sv.installed"
  local want="gnomAD-SV ${GNOMAD_SV_VERSION} sha256 ${GNOMAD_SV_SHA256}"
  local src="${ANNOTSV_DIR}/gnomad.v${GNOMAD_SV_VERSION}.sv.sites.vcf.gz"
  local out="${benign}/gnomad.v4.0.sv.chrAll.from-v${GNOMAD_SV_VERSION}-sites.vcf.gz"
  local n
  if [ "$(cat "$stamp" 2>/dev/null)" = "$want" ]; then
    echo "[OK] gnomAD-SV ${GNOMAD_SV_VERSION} frequencies for AnnotSV already present."
    return 0
  fi
  if [ -s "$stamp" ]; then
    echo "  ${stamp} records '$(cat "$stamp")'; adding ${want}."
    echo "  The rows of the earlier release stay in AnnotSV's benign files; reinstall the bundle to drop them."
  fi
  mkdir -p "$benign"
  # fetch keeps a download only once its checksum matches, so a file there is whole.
  if [ ! -s "$src" ]; then
    echo "Downloading gnomAD-SV ${GNOMAD_SV_VERSION} sites (~1.7 GB)..."
    if ! fetch "$GNOMAD_SV_URL" "$src" sha256 "$GNOMAD_SV_SHA256"; then
      echo "[WARN] Could not download gnomAD-SV ${GNOMAD_SV_VERSION}. Step 5 stops until it is in place; run setup.sh again."
      return 1
    fi
  fi
  echo "  Keeping the records AnnotSV uses for its benign set..."
  if ! gzip -dc "$src" | awk -F'\t' -v OFS='\t' -v ver="$GNOMAD_SV_VERSION" '
      BEGIN {
        print "##fileformat=VCFv4.2"
        print "##source=gnomAD-SV v" ver " sites, all genomes: AN, AF and N_HOMALT under the controls_and_biobanks_* keys AnnotSV 3.5 reads"
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
      }
      /^#/ || $7 != "PASS" { next }
      {
        end = ""; t = ""; an = ""; af = ""; hom = ""
        n = split($8, kv, ";")
        for (i = 1; i <= n; i++) {
          eq = index(kv[i], "=")
          if (!eq) continue
          k = substr(kv[i], 1, eq - 1); v = substr(kv[i], eq + 1)
          if (k == "END") end = v
          else if (k == "SVTYPE") t = v
          else if (k == "AN") an = v
          else if (k == "AF") af = v
          else if (k == "N_HOMALT") hom = v
        }
        if (t != "DEL" && t != "DUP" && t != "INS" && t != "INV") next
        if (end == "" || an == "" || af == "" || hom == "") next
        if (an + 0 < 1000 || af + 0 < 0.001 || hom + 0 < 5) next
        print $1, $2, $3, "N", "<" t ">", ".", "PASS", "END=" end ";SVTYPE=" t ";controls_and_biobanks_AN=" an ";controls_and_biobanks_AF=" af ";controls_and_biobanks_N_HOMALT=" hom
      }' | gzip -c > "${out}.part"; then
    rm -f "${out}.part"
    echo "[WARN] Could not read ${src}; removed nothing. Run setup.sh again."
    return 1
  fi
  n=$(gzip -dc "${out}.part" | grep -vc '^#' || true)
  if [ "$n" -eq 0 ]; then
    rm -f "${out}.part" "$src"
    echo "[WARN] ${src} gave no record AnnotSV can use; removed it. Run setup.sh again."
    return 1
  fi
  mv "${out}.part" "$out"
  if ! annotsv_merge_benign "$out"; then
    rm -f "$out"
    echo "[WARN] AnnotSV did not merge gnomAD-SV ${GNOMAD_SV_VERSION} into its benign set; kept ${src}. Run setup.sh again."
    return 1
  fi
  printf '%s\n' "$want" > "$stamp"
  rm -f "$src"
  echo "[OK] gnomAD-SV ${GNOMAD_SV_VERSION}: ${n} common SVs merged into AnnotSV's benign set (${benign})"
}

# annotsv_merge_benign ADAPTER: AnnotSV merges a gnomad.v4.0.sv.chr*.vcf.gz it
# finds in BenignSV/GRCh38 into its shared benign_*_SV_GRCh38.sorted.bed at
# the start of a run, then deletes it. One AnnotSV run here, as the caller, on
# a one-deletion VCF does that merge at setup: sample runs then find nothing to
# merge, so parallel Nextflow tasks never rewrite those files at once, and
# tasks running as root never take them over. Succeeds when AnnotSV exits 0,
# ADAPTER is gone and the benign loss file exists.
annotsv_merge_benign() {
  local adapter=$1 probe="${ANNOTSV_DIR}/.setup-merge" rc=0
  local loss="${ANNOTSV_DIR}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/benign_Loss_SV_GRCh38.sorted.bed"
  rm -rf "$probe"
  mkdir -p "${probe}/out"
  printf '%s\n' '##fileformat=VCFv4.2' \
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">' \
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">' \
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">' \
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' > "${probe}/probe.vcf"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tprobe\n' >> "${probe}/probe.vcf"
  printf 'chr20\t1000001\tsetup_probe\tN\t<DEL>\t100\tPASS\tEND=1001000;SVTYPE=DEL;SVLEN=-1000\tGT\t0/1\n' >> "${probe}/probe.vcf"
  echo "  AnnotSV merges the records into its benign set (one run, a few minutes)..."
  # The annotations directory is writable: AnnotSV rewrites its benign files there.
  run_in --rw "$ANNOTSV_DIR" --cpus 2 --memory 8g "$ANNOTSV_IMAGE" \
    AnnotSV -SVinputFile "$(cpath "${probe}/probe.vcf")" \
      -outputFile "$(cpath "${probe}/out/probe.tsv")" -outputDir "$(cpath "${probe}/out")" \
      -genomeBuild GRCh38 -annotationsDir "$(cpath "$ANNOTSV_DIR")" > "${probe}.log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ] || [ -e "$adapter" ] || [ ! -s "$loss" ]; then
    echo "[WARN] AnnotSV's merge run: exit ${rc}, $([ -e "$adapter" ] && echo "${adapter##*/} still there" || echo "${adapter##*/} consumed"). Its log: ${probe}.log"
    tail -n 20 "${probe}.log" >&2
    rm -rf "$probe"
    return 1
  fi
  rm -rf "$probe" "${probe}.log"
}

# install_annotsv_mane: the MANE Select RefSeq transcripts of MANE MANE_RELEASE,
# one per line, read by step 05 and the Nextflow ANNOTSV as -txFile.
# mane_select_refseq.release records the release and checksum.
install_annotsv_mane() {
  local out="${ANNOTSV_DIR}/mane_select_refseq.txt"
  local stamp="${ANNOTSV_DIR}/mane_select_refseq.release"
  local want="MANE ${MANE_RELEASE} sha256 ${MANE_SUMMARY_SHA256}"
  local src="${ANNOTSV_DIR}/MANE.GRCh38.v${MANE_RELEASE}.summary.txt.gz"
  local n
  if [ -s "$out" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$want" ]; then
    echo "[OK] MANE ${MANE_RELEASE} Select transcripts for AnnotSV already present."
    return 0
  fi
  if [ ! -s "$src" ] && ! fetch "$MANE_SUMMARY_URL" "$src" sha256 "$MANE_SUMMARY_SHA256"; then
    echo "[WARN] Could not download the MANE ${MANE_RELEASE} summary. Step 5 stops until it is in place; run setup.sh again."
    return 1
  fi
  # Columns by name: RefSeq_nuc holds NM_ ids with a version, MANE_status
  # "MANE Select" or "MANE Plus Clinical".
  if ! gzip -dc "$src" | awk -F'\t' '
      NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
      $col["MANE_status"] == "MANE Select" && $col["RefSeq_nuc"] ~ /^NM_/ { print $col["RefSeq_nuc"] }' > "${out}.part"; then
    rm -f "${out}.part"
    echo "[WARN] Could not read ${src}. Run setup.sh again."
    return 1
  fi
  n=$(wc -l < "${out}.part" | tr -d ' ')
  if [ "$n" -eq 0 ]; then
    rm -f "${out}.part" "$src"
    echo "[WARN] ${src} lists no MANE Select RefSeq transcript; removed it. Run setup.sh again."
    return 1
  fi
  mv "${out}.part" "$out"
  printf '%s\n' "$want" > "$stamp"
  rm -f "$src"
  echo "[OK] MANE ${MANE_RELEASE}: ${n} MANE Select transcripts for AnnotSV (${out})"
}

# install_annotsv_data: the bundle, then the two files step 05 needs with it.
install_annotsv_data() {
  install_annotsv_bundle || return 1
  local rc=0
  install_annotsv_gnomad_sv || rc=1
  install_annotsv_mane || rc=1
  return "$rc"
}

# ensure_image VAR STEP: the image of an opt-in step (`# optional` in
# versions.env, so the main setup leaves it out), pulled with its data, so
# the run does not pull it in the middle.
ensure_image() {
  local img=${!1}
  if "$CONTAINER_ENGINE" image inspect "$img" &>/dev/null; then
    echo "[OK] ${img} (step $2) already pulled."
  elif pull_image "$img"; then
    echo "[OK] ${img} (step $2) pulled."
  else
    echo "[WARN] Could not pull ${img}. Step $2 needs it; run this command again."
    return 1
  fi
  keep_image "$img"
}

# install_vep_cache_opt_in: the VEP cache step 13 and the Nextflow VEP read,
# installed before the first run; without it run-all.sh skips VEP, vcfanno,
# the clinical filter and slivar. The tarball (~26 GB) and the unpacked tree
# (~30 GB) are both on disk until the install ends.
install_vep_cache_opt_in() {
  local dir="${GENOME_DIR}/vep_cache"
  if [ -f "${dir}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
    echo "[OK] VEP ${VEP_CACHE_RELEASE} cache already present: ${dir}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
    return 0
  fi
  echo "Installing the VEP ${VEP_CACHE_RELEASE} cache: ~26 GB download, ~30 GB unpacked, both on disk until it ends."
  if ! install_vep_cache "$dir" "$VEP_CACHE_RELEASE"; then
    echo "[WARN] Could not install the VEP ${VEP_CACHE_RELEASE} cache. Run this command again: the download resumes."
    return 1
  fi
  echo "[OK] VEP ${VEP_CACHE_RELEASE} cache: ${dir}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
}

case "$OPT_IN" in
  ancestry-panel) install_ancestry_panel; exit $? ;;
  ancestry-fine-panel) install_ancestry_fine_panel; exit $? ;;
  cyrius) install_cyrius; exit $? ;;
  parascopy-data) install_parascopy_data && ensure_image PARASCOPY_IMAGE 35; exit $? ;;
  kir-data) install_kir_data; exit $? ;;
  yleaf-data) install_yleaf_data && ensure_image YLEAF_IMAGE 37; exit $? ;;
  consumer-sites) install_consumer_sites; exit $? ;;
  viral-data) install_viral_data; exit $? ;;
  vep-cache) install_vep_cache_opt_in; exit $? ;;
  annotsv-data) install_annotsv_data; exit $? ;;
esac

# pull_images: pull every image setup pre-pulls (versions.env, minus the
# lines marked `# optional`), and tag the digest-pinned ones, those already
# present included. Sets PULLED, SKIPPED, FAILED and TAGGED (pinned tags added).
pull_images() {
  local img
  PULLED=0
  SKIPPED=0
  FAILED=0
  TAGGED=0
  while IFS= read -r img; do
    if "$CONTAINER_ENGINE" image inspect "$img" &>/dev/null; then
      SKIPPED=$((SKIPPED + 1))
    else
      echo "  Pulling: ${img}..."
      if pull_image "$img"; then
        PULLED=$((PULLED + 1))
      else
        echo "  WARNING: Failed to pull ${img}; Docker's message is above."
        FAILED=$((FAILED + 1))
        continue
      fi
    fi
    pinned_tag_present "$img" && continue
    if keep_image "$img"; then TAGGED=$((TAGGED + 1)); else FAILED=$((FAILED + 1)); fi
  done < <(pipeline_images)
  echo "[OK] Docker images: ${PULLED} pulled, ${SKIPPED} already present, ${FAILED} failed, ${TAGGED} pgp-pinned tag(s) added."
  if [ "$FAILED" -gt 0 ]; then
    echo ""
    echo "ERROR: ${FAILED} Docker image(s) failed to pull or tag. Fix the issues above and re-run setup."
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
    gencode_exons) what="GENCODE ${GENCODE_RELEASE} exons (steps 06 and 17)" ;;
  esac
  if dest=$(data_file "$name"); then
    echo "[OK] ${what} already present."
  elif install_data_file "$name"; then
    echo "[OK] ${what}: ${dest}"
  else
    echo "[WARN] Could not install the ${what}. Re-run setup.sh to try again."
  fi
done
# gnomAD v3.1's chrM sites (step 20 annotates the calls with them), checked
# against GNOMAD_MITO_SHA256. A failed download does not stop setup either.
GNOMAD_MITO="${GENOME_DIR}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
if [ -s "$GNOMAD_MITO" ] && [ "$(_digest sha256 "$GNOMAD_MITO")" = "$GNOMAD_MITO_SHA256" ]; then
  echo "[OK] gnomAD v3.1 mitochondrial sites (step 20) already present (sha256 matches GNOMAD_MITO_SHA256)."
else
  if [ -s "$GNOMAD_MITO" ]; then
    echo "[WARN] ${GNOMAD_MITO} is not the pinned file (sha256 differs from GNOMAD_MITO_SHA256): downloading it again."
  fi
  if fetch "$GNOMAD_MITO_URL" "$GNOMAD_MITO" sha256 "$GNOMAD_MITO_SHA256"; then
    echo "[OK] gnomAD v3.1 mitochondrial sites (step 20): ${GNOMAD_MITO}"
  else
    echo "[WARN] Could not install gnomAD's mitochondrial sites. Step 20 runs without its gnomAD columns; re-run setup.sh to try again."
  fi
fi
# A failed download does not stop setup here either.
install_sample_qc_data || true
# Step 25: its scores, pgsc_calc and its plugin (warnings only).
install_prs

###############################################################################
# Phase 5: AnnotSV annotation data (step 5, a default step)
###############################################################################
# The bundle, gnomAD-SV's common SVs and the MANE Select transcripts
# (install_annotsv_data above). A failure is a warning: the rest of setup goes on.
echo ""
echo "=== Phase 5: AnnotSV annotation data (~5.3 GB + ~1.7 GB download, ~20 GB unpacked) ==="

install_annotsv_data || true

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
  echo "[SKIP] VEP ${VEP_CACHE_RELEASE} cache (~26 GB) — without it a run skips VEP, vcfanno, the clinical filter and slivar"
  echo "  Install it before the first run (downloads, checks and unpacks it; ~30 GB unpacked):"
  echo "    ./scripts/setup.sh --vep-cache ${GENOME_DIR}"
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
