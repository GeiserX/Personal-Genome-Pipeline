#!/usr/bin/env bash
# What the 1000 Genomes ancestry panel costs on a GitHub-hosted runner, for
# docs/25-prs.md: the download, the disk at its fullest and the memory at its
# highest while pgsc_calc projects a genome-wide target onto the panel and
# scores it. Runs on the monthly and dispatched E2E runs, and with
# PGSC_MEASURE=1; a pull request run skips it (about 7 GB more to download).
#
#   1. setup.sh --ancestry-panel installs pgsc_1000G_v1 and its site list;
#   2. pgsc_calc (the release step 25 installed) runs with --run_ancestry on
#      the PGS Catalog's synthetic genome-wide test target (HAPNEST, 600
#      samples) and PGS000018 (1.7 million variants), so the panel side of
#      the work (intersection, QC, PCA, projection) is the size a real genome
#      gives, which the fixture's slices are not;
#   3. a sampler records used memory and free disk every 5 s; the trace adds
#      each task's peak RSS.
# Leaves the panel installed for case prs-4, which projects HG002 onto it
# and then removes it.
. "$(dirname "$0")/lib.sh"

if [ "${GITHUB_EVENT_NAME:-}" = pull_request ] && [ "${PGSC_MEASURE:-}" != 1 ]; then
  echo "pull request run: the 1000 Genomes panel measurement runs on the monthly and dispatched runs only."
  finish
fi
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
. "${REPO}/versions.env"
CALC="${GENOME_DIR}/tools/pgsc_calc-${PGSC_CALC_VERSION}"
check "step 25 (case prs-1) installed pgsc_calc" test -f "${CALC}/main.nf"
RES=https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources
PANEL="${GENOME_DIR}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst"
MEASURE="${E2E_WORK}/logs/pgsc-measure.tsv"
mkdir -p "$(dirname "$MEASURE")"

# --- sampler: used memory (MB) and free disk of the work area (KB) -----------------
SAMPLES="${CASE_TMP}/samples.tsv"
: > "$SAMPLES"
( trap 'exit 0' TERM
  while :; do
    printf '%s\t%s\t%s\n' "$(date +%s)" "$(free -m | awk '/^Mem:/ {print $3}')" "$(df -k --output=avail "$E2E_WORK" | tail -1)" >> "$SAMPLES"
    sleep 5
  done ) &
SAMPLER=$!
trap 'kill "$SAMPLER" 2>/dev/null || true' EXIT
sleep 6
BASE_MEM=$(awk -F'\t' 'NR == 1 {print $2}' "$SAMPLES")
BASE_DISK=$(awk -F'\t' 'NR == 1 {print $3}' "$SAMPLES")

# --- 1. the panel -------------------------------------------------------------------
t0=$(date +%s)
"${REPO}/scripts/setup.sh" --ancestry-panel "$GENOME_DIR" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --ancestry-panel exits 0" "${PIPESTATUS[0]}" 0
T_SETUP=$(( $(date +%s) - t0 ))
PANEL_BYTES=$(stat -c %s "$PANEL" 2>/dev/null || echo 0)
N_SITES=$(wc -l < "${PANEL%.tar.zst}_GRCh38_sites.tsv" 2>/dev/null | tr -d ' ')
check_ge "panel size in bytes" "$PANEL_BYTES" 1000000000
check_ge "GRCh38 SNVs in the site list" "${N_SITES:-0}" 100000

# --- 2. a genome-wide target projected and scored ---------------------------------
G="${CASE_TMP}/genomewide"
rm -rf "$G"
mkdir -p "$G"
curl -fsSL -o "${G}/target.tar.zst" "${RES}/GRCh38_HAPNEST_target.tar.zst"
tar -C "$G" --zstd -xf "${G}/target.tar.zst" && rm -f "${G}/target.tar.zst"
PREFIX=$(find "$G" -name '*.pgen' | head -1); PREFIX=${PREFIX%.pgen}
check "the HAPNEST target is unpacked" test -s "${PREFIX}.pgen"
curl -fsSL -o "${G}/PGS000018_hmPOS_GRCh38.txt.gz" \
  https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000018/ScoringFiles/Harmonized/PGS000018_hmPOS_GRCh38.txt.gz
printf 'sampleset,path_prefix,chrom,format\ntest,%s,,pfile\n' "$PREFIX" > "${G}/samplesheet.csv"
t0=$(date +%s)
( cd "$G" && nextflow -log "${G}/nextflow.log" run "${CALC}/main.nf" -profile docker -ansi-log false -work-dir "${G}/work" \
    --input "${G}/samplesheet.csv" --target_build GRCh38 --scorefile "${G}/PGS000018_hmPOS_GRCh38.txt.gz" \
    --run_ancestry "$PANEL" --min_overlap 0.5 --outdir "${G}/results" \
    --max_cpus "$(nproc)" --max_memory "$(( $(free -g | awk '/^Mem:/ {print $2}') - 1 )).GB" ) 2>&1 | tee "$STEP_LOG"
RC=${PIPESTATUS[0]}
T_RUN=$(( $(date +%s) - t0 ))
check_eq "pgsc_calc with the 1000 Genomes panel on the genome-wide target exits 0" "$RC" 0
check "it wrote ancestry-adjusted scores" test -s "${G}/results/test/score/test_pgs.txt.gz"
sleep 6
kill "$SAMPLER" 2>/dev/null || true
WORK_BYTES=$(du -sb "${G}/work" 2>/dev/null | cut -f1)
TRACE=$(find "${G}/results/pipeline_info" -name 'execution_trace_*.txt' 2>/dev/null | head -1)
# The five tasks with the highest peak RSS, as Nextflow's trace gives it (e.g. "3.2 GB").
TOP=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
  { v = $c["peak_rss"]; n = v + 0; if (v ~ /GB/) n *= 1024; else if (v ~ /KB/) n /= 1024; else if (v !~ /MB/) n /= 1048576
    printf "%.0f\t%s\n", n, $c["name"] }' "${TRACE:-/dev/null}" 2>/dev/null | sort -rn | head -5)
PEAK_TASK_MB=$(head -1 <<< "$TOP" | cut -f1)
PEAK_MEM=$(awk -F'\t' 'NR == 1 || $2 > m {m = $2} END {print m}' "$SAMPLES")
MIN_DISK=$(awk -F'\t' 'NR == 1 || $3 < m {m = $3} END {print m}' "$SAMPLES")
DISK_GB=$(awk -v b="$BASE_DISK" -v m="$MIN_DISK" 'BEGIN {printf "%.1f", (b - m) / 1048576}')
MEM_GB=$(awk -v b="$BASE_MEM" -v p="$PEAK_MEM" 'BEGIN {printf "%.1f", (p - b) / 1024}')
{
  printf 'item\tvalue\n'
  printf 'panel\t%s\n' "$PGSC_PANEL"
  printf 'panel_download_bytes\t%s\n' "$PANEL_BYTES"
  printf 'panel_grch38_snvs\t%s\n' "$N_SITES"
  printf 'setup_seconds\t%s\n' "$T_SETUP"
  printf 'genomewide_run_seconds\t%s\n' "$T_RUN"
  printf 'disk_peak_gb_above_start\t%s\n' "$DISK_GB"
  printf 'pgsc_calc_work_bytes_at_end\t%s\n' "${WORK_BYTES:-0}"
  printf 'memory_peak_gb_above_start\t%s\n' "$MEM_GB"
  printf 'task_peak_rss_mb\t%s\n' "${PEAK_TASK_MB:-0}"
  printf 'runner_memory_gb\t%s\n' "$(free -g | awk '/^Mem:/ {print $2}')"
  printf 'runner_cpus\t%s\n' "$(nproc)"
} > "$MEASURE"
cat "$MEASURE"; echo "tasks with the highest peak RSS (MB, name):"; echo "$TOP"
{
  printf '#### pgsc_calc with the 1000 Genomes panel (%s), genome-wide target\n\n```\n' "$PGSC_PANEL"
  cat "$MEASURE"; echo "top tasks by peak RSS (MB):"; echo "$TOP"
  printf '```\n\n'
} >> "$E2E_NOTES"
check_ge "a task's peak RSS was recorded (MB)" "${PEAK_TASK_MB:-0}" 1
rm -rf "$G"
finish
