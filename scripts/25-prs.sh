#!/usr/bin/env bash
# 25-prs.sh — Polygenic scores with pgsc_calc, the PGS Catalog's calculator
# Usage: ./scripts/25-prs.sh <sample_name>
#
# Scores the sample with every score in assets/pgs_scores.tsv. pgsc_calc
# (PGSC_CALC_VERSION in versions.env, a Nextflow pipeline of the PGS Catalog
# team) matches each score's variants to your genotypes and sums them. With
# the ancestry reference panel installed (scripts/setup.sh --ancestry-panel)
# it also places you among the panel's samples and reports where your score
# falls among the group whose genetic ancestry is most similar to yours: a
# percentile. Without the panel the summary says "raw score only": a raw sum
# cannot be compared with anyone. A score is not a diagnosis either way.
#
# Requires: the VCF of step 3, Java 17+ and Nextflow (as for run-all.sh). With
# step 3's gVCF next to the VCF, the score positions (and the panel's) are
# genotyped from the gVCF, so a site where you match the reference is a real
# 0/0 instead of a missing site. Without it, only your variant sites are
# scored and the summary's Input column says vcf.
#
# Env: ANCESTRY_PANEL  the panel (default reference/pgsc_calc/${PGSC_PANEL}.tar.zst;
#                      "none" scores without it even when it is installed)
#      PGSC_MAX_MEMORY memory pgsc_calc may use, e.g. 12.GB (default: 3/4 of the RAM)
#      PGSC_CALC_DIR   pgsc_calc's code (default tools/pgsc_calc-${PGSC_CALC_VERSION},
#                      which setup.sh or this step installs)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image PYTHON_IMAGE BCFTOOLS_IMAGE PLINK2_IMAGE PGSC_UTILS_IMAGE PGSC_FRAPOSA_IMAGE \
  PGSC_PYYAML_IMAGE PGSC_ZSTD_IMAGE PGSC_REPORT_IMAGE

VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
GVCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
OUTDIR="${GENOME_DIR}/${SAMPLE}/prs"
WORK="${OUTDIR}/pgsc_calc"
SCORING_DIR="${GENOME_DIR}/prs_scores"
SCORE_LIST="${PGP_ROOT}/assets/pgs_scores.tsv"
PANEL=${ANCESTRY_PANEL:-${GENOME_DIR}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst}
# setup.sh --ancestry-panel writes the panel's GRCh38 SNVs beside it.
PANEL_SITES="${PANEL%.tar.zst}_GRCh38_sites.tsv"
PGSC_CALC_DIR=${PGSC_CALC_DIR:-${GENOME_DIR}/tools/pgsc_calc-${PGSC_CALC_VERSION}}
# pgsc_calc names a sampleset with letters and digits only; the sample keeps its name inside the VCF.
SAMPLESET=sample
mkdir -p "$OUTDIR" "$SCORING_DIR"

if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF not found: ${VCF}"
  echo "  Run step 3 (DeepVariant) first."
  exit 1
fi
for t in nextflow java; do
  if ! command -v "$t" >/dev/null 2>&1; then
    echo "ERROR: ${t} is not on PATH. Step 25 runs pgsc_calc, a Nextflow pipeline: install Java 17+ and Nextflow ${NEXTFLOW_VERSION} (docs/nextflow.md)." >&2
    exit 1
  fi
done

USE_PANEL=false
if [ "${ANCESTRY_PANEL:-}" != none ] && [ -f "$PANEL" ]; then
  if [ ! -s "$PANEL_SITES" ]; then
    echo "ERROR: the ancestry panel ${PANEL} has no site list beside it (${PANEL_SITES})." >&2
    echo "  Run scripts/setup.sh --ancestry-panel ${GENOME_DIR} again, or ANCESTRY_PANEL=none for raw scores." >&2
    exit 1
  fi
  USE_PANEL=true
fi

echo "============================================"
echo "  Step 25: Polygenic scores"
echo "  Tool: pgsc_calc ${PGSC_CALC_VERSION}"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${VCF}"
if $USE_PANEL; then
  echo "  Ancestry panel: ${PANEL}"
else
  echo "  Ancestry panel: none (raw scores only; scripts/setup.sh --ancestry-panel installs it)"
fi
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# --- 1. Scoring files -------------------------------------------------------------
# Only GRCh38-harmonised files are used. The author-reported files are often
# GRCh37 or rsID-only, and scoring them against a GRCh38 VCF gives a number
# that looks valid but is not, so there is no fallback to them.
PGS_BASE_URL=${PGS_BASE_URL:-https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores}
mapfile -t PGS_IDS < <(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "$SCORE_LIST")
if [ "${#PGS_IDS[@]}" -eq 0 ]; then
  echo "ERROR: no PGS id in ${SCORE_LIST}" >&2
  exit 1
fi
# Prints the genome build recorded in a scoring file's #HmPOS_build header (empty if none).
hm_build() {
  { gzip -cd "$1" 2>/dev/null || true; } | awk -F= '/^#HmPOS_build=/ {b=$2} !/^#/ {exit} END {print b}'
}
echo "[1/5] PGS Catalog scoring files (${#PGS_IDS[@]} scores, ${SCORE_LIST#"${PGP_ROOT}/"})..."
for PGS_ID in "${PGS_IDS[@]}"; do
  SCORE_FILE="${SCORING_DIR}/${PGS_ID}.txt.gz"
  # A cached file from an older version may be the author-reported build.
  if [ -f "$SCORE_FILE" ] && [ "$(hm_build "$SCORE_FILE")" != "GRCh38" ]; then
    echo "  Cached ${PGS_ID} is not a GRCh38-harmonised file; downloading it again."
    rm -f "$SCORE_FILE"
  fi
  if [ -f "$SCORE_FILE" ]; then
    continue
  fi
  echo "  Downloading ${PGS_ID}..."
  URL="${PGS_BASE_URL}/${PGS_ID}/ScoringFiles/Harmonized/${PGS_ID}_hmPOS_GRCh38.txt.gz"
  # The PGS Catalog publishes an md5 next to every scoring file.
  if ! fetch "$URL" "$SCORE_FILE" md5 "${URL}.md5"; then
    echo "ERROR: Could not download the GRCh38-harmonised scoring file for ${PGS_ID}:" >&2
    echo "  ${URL}" >&2
    exit 1
  fi
done

# --- 2. pgsc_calc's scoring files --------------------------------------------------
# bin/collect_summary.py prs-format (the Nextflow PRS_PREPARE process runs the
# same) writes each file as a custom GRCh38 file: chr_name and chr_position
# from the harmonised hm_chr and hm_pos, labelled with the catalog's trait,
# so pgsc_calc never needs the network or a liftover. It refuses a file that
# is not harmonised to GRCh38 or not additive.
echo ""
echo "[2/5] Writing the scores as pgsc_calc reads them..."
rm -rf "$WORK"
mkdir -p "${WORK}/scores"
IDS=$(IFS=,; echo "${PGS_IDS[*]}")
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${SCORE_LIST}:/pgs_scores.tsv:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/collect_summary.py prs-format \
    --scores "$(cpath "$SCORING_DIR")" --ids "$IDS" --labels /pgs_scores.tsv \
    --out "$(cpath "${WORK}/scores")" --alleles "$(cpath "${WORK}/score_alleles.tsv")"

# --- 3. The genotypes pgsc_calc scores ---------------------------------------------
# Both inputs end as target.vcf.gz: the score positions (and with the panel the
# panel's SNVs) only. Both lists are autosomal (prs-format drops a score's
# chrX, chrY and MT rows; the panel list is chromosomes 1 to 22), which keeps
# chrX out (plink2 refuses it without the sample's sex) and gives pgsc_calc a
# small file to convert.
echo ""
TARGET="${WORK}/target.vcf.gz"
# ALLELES: every candidate ALT of each position (score effect and other
# alleles, the panel's ALT); SITES: each position once, for bcftools -R/-T.
ALLELES="${WORK}/alleles.tsv"
SITES="${WORK}/sites.tsv"
if $USE_PANEL; then
  { cat "${WORK}/score_alleles.tsv"; cut -f1,2,4 "$PANEL_SITES"; } | LC_ALL=C sort -u -k1,1 -k2,2n -k3,3 > "$ALLELES"
else
  cp "${WORK}/score_alleles.tsv" "$ALLELES"
fi
cut -f1,2 "$ALLELES" | uniq > "$SITES"
if [ -f "$GVCF" ] && [ -f "${GVCF}.tbi" ]; then
  INPUT_KIND=gvcf
  echo "[3/5] Genotyping the score positions$($USE_PANEL && echo " and the panel's") from the gVCF (${GVCF})..."
  # gvcf2vcf expands each reference block overlapping a position into one 0/0
  # record per base, with the base from the reference; -T keeps the positions,
  # --trim-alt-alleles drops the <*> allele, and a no-call (./.) is dropped,
  # so a position without coverage stays missing. A 0/0 record then has ALT
  # '.', which no allele can match: the awk sets ALT to the position's first
  # candidate allele that is not the reference, so a score allele other than
  # the reference matches it with a dosage of 0 and counts as matched, and
  # the panel's variant is found in the sample.
  # shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
  run_in --cpus 2 --memory 8g \
    "${BCFTOOLS_IMAGE}" \
    bash -euo pipefail -c '
      bcftools convert --gvcf2vcf -f "$1" -R "$2" -Ou "$3" \
        | bcftools view -T "$2" --trim-alt-alleles -i "GT!=\"mis\"" -Ov' \
    _ "${REF_FASTA_C}" "$(cpath "$SITES")" "/genome/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" \
  | awk -F'\t' -v OFS='\t' -v alleles="$ALLELES" '
      BEGIN { while ((getline l < alleles) > 0) { split(l, f, "\t"); k = f[1] ":" f[2]; ea[k] = ea[k] " " f[3] } }
      /^#/ { print; next }
      $5 == "." {
        n = split(ea[$1 ":" $2], c, " ")
        for (i = 1; i <= n; i++) if (c[i] != $4) { $5 = c[i]; break }
      }
      { print }' | gzip -c > "${TARGET}.tmp"
else
  INPUT_KIND=vcf
  echo "[3/5] No gVCF beside the VCF: scoring the variant sites only (sites where you match the reference are missing)..."
  # shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
  run_in --cpus 1 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -euo pipefail -c 'bcftools view -T "$1" -Ov "$2"' \
    _ "$(cpath "$SITES")" "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | gzip -c > "${TARGET}.tmp"
fi
mv -f "${TARGET}.tmp" "$TARGET"
rm -f "$ALLELES" "$SITES"
N_TARGET=$({ gzip -dc "$TARGET" | grep -vc '^#'; } || true)
echo "  ${N_TARGET} genotyped positions for pgsc_calc"

# --- 4. pgsc_calc ------------------------------------------------------------------
echo ""
echo "[4/5] pgsc_calc ${PGSC_CALC_VERSION}..."
if [ ! -f "${PGSC_CALC_DIR}/main.nf" ]; then
  # GitHub's archive of the tag, checked against PGSC_CALC_SHA256 (setup.sh holds the same).
  echo "  Installing pgsc_calc ${PGSC_CALC_VERSION} into ${PGSC_CALC_DIR} (setup.sh does this once)..."
  fetch "https://github.com/PGScatalog/pgsc_calc/archive/refs/tags/${PGSC_CALC_VERSION}.tar.gz" "${PGSC_CALC_DIR}.tar.gz" \
    sha256 "$PGSC_CALC_SHA256"
  rm -rf "${PGSC_CALC_DIR}.part"
  mkdir -p "${PGSC_CALC_DIR}.part"
  tar -xzf "${PGSC_CALC_DIR}.tar.gz" -C "${PGSC_CALC_DIR}.part" --strip-components 1
  rm -f "${PGSC_CALC_DIR}.tar.gz"
  # An incomplete folder left by an earlier run would receive the new one inside it.
  rm -rf "$PGSC_CALC_DIR"
  mv "${PGSC_CALC_DIR}.part" "$PGSC_CALC_DIR"
fi
NXF_HOME=${NXF_HOME:-${HOME}/.nextflow}
if [ ! -d "${NXF_HOME}/plugins/nf-schema-${PGSC_CALC_NF_SCHEMA}" ]; then
  echo "  Installing the nf-schema ${PGSC_CALC_NF_SCHEMA} plugin pgsc_calc needs..."
  nextflow plugin install "nf-schema@${PGSC_CALC_NF_SCHEMA}"
fi
# Nextflow refuses a task that asks for more CPUs than the machine has.
PGSC_CPUS=$THREADS
NCPU=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo "$THREADS")
[ "$NCPU" -ge "$PGSC_CPUS" ] 2>/dev/null || PGSC_CPUS=$NCPU
if [ -z "${PGSC_MAX_MEMORY:-}" ]; then
  if [ -r /proc/meminfo ]; then
    kb=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
  else
    kb=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 17179869184) / 1024 ))
  fi
  PGSC_MAX_MEMORY="$(( kb * 3 / 4 / 1048576 )).GB"
fi
# docker_ref IMAGE: IMAGE with docker.io/ in front when it names no registry:
# pgsc_calc sets docker.registry to quay.io, which would be put there instead.
docker_ref() {
  case "${1%%/*}" in
    *.*|*:*|localhost) printf '%s' "$1" ;;
    *) printf 'docker.io/%s' "$1" ;;
  esac
}
# The images of versions.env for pgsc_calc's process labels (the same table
# as scripts/ci/gen-containers-config.sh's NATIVE row), and no network for its
# containers: the scores, the genotypes and the panel are all local.
{
  echo 'process {'
  for pair in pgscatalog_utils=PGSC_UTILS_IMAGE plink2=PLINK2_IMAGE zstd=PGSC_ZSTD_IMAGE \
              report=PGSC_REPORT_IMAGE pyyaml=PGSC_PYYAML_IMAGE fraposa=PGSC_FRAPOSA_IMAGE; do
    v=${pair#*=}
    printf "    withLabel: '%s' { ext.docker = '%s'; ext.docker_version = '' }\n" "${pair%%=*}" "$(docker_ref "${!v}")"
  done
  echo '}'
  # shellcheck disable=SC2016  # $(id -u) is expanded by the task's shell
  echo 'docker.runOptions = '"'"'-u $(id -u):$(id -g) --network none'"'"
} > "${WORK}/images.config"
printf 'sampleset,path_prefix,chrom,format\n%s,%s,,vcf\n' "$SAMPLESET" "${WORK}/target" > "${WORK}/samplesheet.csv"
PANEL_ARGS=()
if $USE_PANEL; then PANEL_ARGS=(--run_ancestry "$PANEL"); fi
LOG="${WORK}/pgsc_calc.log"
rc=0
ZERO=()
if [ "$N_TARGET" -eq 0 ]; then
  # Not one score position is in the input: there is nothing for pgsc_calc to match.
  echo "  None of the score positions is in this input; no score." | tee "$LOG"
  ZERO=(--zero-matches)
else
  # A fresh work folder each run: pgsc_calc keeps converted genotypes with
  # storeDir, and an old run's would be scored instead of this VCF. The
  # nf-core institutional configs are not fetched (NXF_OFFLINE).
  ( cd "$WORK" && NXF_OFFLINE=true nextflow -log "${WORK}/nextflow.log" run "${PGSC_CALC_DIR}/main.nf" \
      -profile docker -c "${WORK}/images.config" -work-dir "${WORK}/work" -ansi-log false \
      --input "${WORK}/samplesheet.csv" --target_build GRCh38 \
      --scorefile "${WORK}/scores/*.txt.gz" \
      ${PANEL_ARGS[@]+"${PANEL_ARGS[@]}"} \
      --outdir "${WORK}/results" \
      --max_cpus "$PGSC_CPUS" --max_memory "$PGSC_MAX_MEMORY" ) > "$LOG" 2>&1 || rc=$?
  cat "$LOG"
fi
if [ "$rc" -ne 0 ]; then
  # Every score under pgsc_calc's minimum overlap: no sum, the match rates
  # from its log. Not one score variant in the genotypes: every score unmatched.
  if grep -q 'All scores fail to meet match threshold' "$LOG" "${WORK}/nextflow.log" 2>/dev/null; then
    echo "  Every score matched under pgsc_calc's minimum overlap of its variants in this input; no score."
    ZERO=(--below-threshold "$(cpath "$LOG")")
  elif grep -qE 'ZeroMatchesError|No match candidates found for any scoring files' "$LOG" "${WORK}/nextflow.log" 2>/dev/null; then
    echo "  None of the score variants is in this input; no score."
    ZERO=(--zero-matches)
  else
    echo "ERROR: pgsc_calc failed (exit ${rc}); see ${LOG} and ${WORK}/nextflow.log" >&2
    exit 1
  fi
fi

# --- 5. Summary ---------------------------------------------------------------------
echo ""
echo "[5/5] Summary..."
RESULTS_FILE="${OUTDIR}/${SAMPLE}_prs_summary.tsv"
ANCESTRY_DIR="${GENOME_DIR}/${SAMPLE}/ancestry"
ANCESTRY_FILE="${ANCESTRY_DIR}/${SAMPLE}_ancestry.tsv"
rm -f "$RESULTS_FILE"
ANC_ARGS=()
if $USE_PANEL; then
  mkdir -p "$ANCESTRY_DIR"
  rm -f "$ANCESTRY_FILE"
  ANC_ARGS=(--panel "$(basename "${PANEL%.tar.zst}")" --ancestry-out "$(cpath "$ANCESTRY_FILE")")
fi
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/collect_summary.py prs-table \
    --sample "$SAMPLE" --results "$(cpath "${WORK}/results")" --sampleset "$SAMPLESET" \
    --scores "$(cpath "${WORK}/scores")" --input-kind "$INPUT_KIND" \
    ${ZERO[@]+"${ZERO[@]}"} ${ANC_ARGS[@]+"${ANC_ARGS[@]}"} --out "$(cpath "$RESULTS_FILE")"
# Kept: pgsc_calc's results (its HTML report, the match log); dropped: its work
# folder and its run reports, whose names carry the time of the run.
rm -rf "${WORK}/work" "${WORK}/.nextflow" "${WORK}/results/pipeline_info"

if awk -F'\t' 'NR > 1 && $6 + 0 < 50 {found = 1} END {exit !found}' "$RESULTS_FILE"; then
  echo "WARNING: under half of a score's variants were genotyped; that sum is not comparable to published distributions."
fi
if [ "$INPUT_KIND" = vcf ]; then
  echo "NOTE: hom-ref sites are absent from this VCF, so the score is biased; not comparable to published distributions"
fi

echo ""
echo "============================================"
echo "  Polygenic scores complete: ${SAMPLE}"
echo ""
echo "  Summary: ${RESULTS_FILE}"
column -t -s $'\t' "$RESULTS_FILE" 2>/dev/null || cat "$RESULTS_FILE"
if [ -s "$ANCESTRY_FILE" ]; then
  echo "  Ancestry (step 26's table): ${ANCESTRY_FILE}"
fi
echo "  pgsc_calc's report: ${WORK}/results/${SAMPLESET}/score/report.html"
echo "============================================"
echo ""
echo "IMPORTANT: a polygenic score is NOT diagnostic. It estimates relative genetic"
echo "predisposition; lifestyle, environment and other genes are not captured."
if $USE_PANEL; then
  echo "The percentile compares your score with the reference samples whose genetic"
  echo "ancestry is most similar to yours. See docs/25-prs.md."
else
  echo "Raw score only: without the ancestry reference panel there is no percentile,"
  echo "and a raw sum cannot be compared with anyone (scripts/setup.sh --ancestry-panel)."
  echo "See docs/25-prs.md."
fi
