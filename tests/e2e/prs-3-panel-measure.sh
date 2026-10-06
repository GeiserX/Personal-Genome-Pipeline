#!/usr/bin/env bash
# What the 1000 Genomes ancestry panel costs on a GitHub-hosted runner, for
# docs/25-prs.md: the download, the disk at its fullest and the memory at its
# highest while pgsc_calc projects a genome-wide target onto the panel and
# scores it. Runs on a dispatched E2E run (or with PGSC_MEASURE=1) only: the
# panel is 7.4 GB, and the monthly and pull request runs have no room for it
# in their time budget.
#
#   --download   (started in the background by case prs-1, so the download
#                runs while other cases do) fetches the panel with 8 byte
#                ranges at once: EBI serves one connection at about 1.2 MB/s
#                to a runner, 100 minutes for the whole file. The parts are
#                joined and checked against the md5 the PGS Catalog lists.
#   (no option)  waits for that download; setup.sh --ancestry-panel then
#                writes the panel's site list beside it; pgsc_calc (the
#                release step 25 installed) runs with --run_ancestry on the
#                PGS Catalog's synthetic genome-wide test target (HAPNEST, 600
#                samples) and PGS000018 (1.7 million variants), so the panel
#                side of the work (intersection, QC, PCA, projection) is the
#                size a real genome gives, which the fixture's slices are not.
#                A sampler records used memory and free disk every 5 s; the
#                trace adds each task's peak RSS. pgsc_calc gets 35 minutes:
#                a run that does not finish is reported as such, with what it
#                used until then.
# Case prs-2 projects HG002 onto the small synthetic panel on every run.
# Removes the panel and everything else it made.
. "$(dirname "$0")/lib.sh"

measure_on() { [ "${GITHUB_EVENT_NAME:-}" = workflow_dispatch ] || [ "${PGSC_MEASURE:-}" = 1 ]; }
RES=https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources
DIR="${GENOME_DIR}/reference/pgsc_calc"
PANEL="${DIR}/${PGSC_PANEL}.tar.zst"
STATUS="${E2E_WORK}/panel-download.status"
MEASURE="${E2E_WORK}/logs/pgsc-measure.tsv"

# --- background download ----------------------------------------------------------
if [ "${1:-}" = --download ]; then
  mkdir -p "$DIR"
  t0=$(date +%s)
  url="${RES}/${PGSC_PANEL}.tar.zst"
  size=$(curl -fsSIL "$url" | tr -d '\r' | awk 'tolower($1) == "content-length:" {n = $2} END {print n}')
  want=$(curl -fsSL "${RES}/md5s.txt" | awk -v f="${PGSC_PANEL}.tar.zst" '$2 == f {print $1}')
  if [ -z "$size" ] || [ -z "$want" ]; then echo "fail no-size-or-md5" > "$STATUS"; exit 1; fi
  n=8 chunk=$(( (size + 7) / 8 ))
  pids=()
  for i in $(seq 0 $((n - 1))); do
    a=$((i * chunk)); b=$(( (i + 1) * chunk - 1 )); [ "$b" -lt "$size" ] || b=$((size - 1))
    curl -fsSL --retry 5 --retry-delay 10 -r "${a}-${b}" -o "${PANEL}.part${i}" "$url" & pids+=($!)
  done
  ok=true
  for p in "${pids[@]}"; do wait "$p" || ok=false; done
  if $ok; then
    for i in $(seq 0 $((n - 1))); do cat "${PANEL}.part${i}"; done > "${PANEL}.part"
    got=$(md5sum "${PANEL}.part" | cut -d' ' -f1)
    if [ "$got" = "$want" ]; then mv "${PANEL}.part" "$PANEL"; else ok=false; fi
  fi
  rm -f "${PANEL}".part*
  if $ok; then
    echo "ok $(( $(date +%s) - t0 )) ${size}" > "$STATUS"
  else
    echo "fail download-or-md5 $(( $(date +%s) - t0 ))" > "$STATUS"
  fi
  exit 0
fi

if ! measure_on; then
  echo "not a dispatched run: the 1000 Genomes panel measurement runs on a dispatched E2E run only (or PGSC_MEASURE=1)."
  finish
fi
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
CALC="${GENOME_DIR}/tools/pgsc_calc-${PGSC_CALC_VERSION}"
check "step 25 (case prs-1) installed pgsc_calc" test -f "${CALC}/main.nf"
mkdir -p "$(dirname "$MEASURE")"

# --- the download prs-1 started -------------------------------------------------------
deadline=$(( $(date +%s) + 30 * 60 ))
until [ -s "$STATUS" ] || [ "$(date +%s)" -ge "$deadline" ]; do sleep 20; done
echo "panel download: $(cat "$STATUS" 2>/dev/null || echo 'not finished after 30 minutes')"
read -r D_STATE D_SECONDS D_BYTES < "$STATUS" 2>/dev/null || true
check_eq "the panel download (8 ranges, md5 checked)" "${D_STATE:-none}" ok
if [ "${D_STATE:-}" != ok ]; then
  printf 'item\tvalue\npanel\t%s\ndownload\t%s\n' "$PGSC_PANEL" "$(cat "$STATUS" 2>/dev/null || echo unfinished)" > "$MEASURE"
  finish
fi

# --- sampler: used memory (MB) and free disk of the work area (KB) -----------------
SAMPLES="${CASE_TMP}/samples.tsv"
: > "$SAMPLES"
( trap 'exit 0' TERM
  while :; do
    printf '%s\t%s\t%s\n' "$(date +%s)" "$(free -m | awk '/^Mem:/ {print $3}')" "$(df -k --output=avail "$E2E_WORK" | tail -1)" >> "$SAMPLES"
    sleep 5
  done ) &
SAMPLER=$!
trap 'kill "$SAMPLER" 2>/dev/null || true; rm -rf "${CASE_TMP}/genomewide" "$DIR"' EXIT
sleep 6
BASE_MEM=$(awk -F'\t' 'NR == 1 {print $2}' "$SAMPLES")
BASE_DISK=$(awk -F'\t' 'NR == 1 {print $3}' "$SAMPLES")

# --- the site list ----------------------------------------------------------------------
t0=$(date +%s)
"${REPO}/scripts/setup.sh" --ancestry-panel "$GENOME_DIR" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --ancestry-panel exits 0" "${PIPESTATUS[0]}" 0
T_SITES=$(( $(date +%s) - t0 ))
N_SITES=$(wc -l < "${PANEL%.tar.zst}_GRCh38_sites.tsv" 2>/dev/null | tr -d ' ')
check_ge "GRCh38 SNVs in the site list" "${N_SITES:-0}" 100000

# --- a genome-wide target projected and scored ---------------------------------
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
( cd "$G" && timeout -k 60 2100 nextflow -log "${G}/nextflow.log" run "${CALC}/main.nf" -profile docker -ansi-log false -work-dir "${G}/work" \
    --input "${G}/samplesheet.csv" --target_build GRCh38 --scorefile "${G}/PGS000018_hmPOS_GRCh38.txt.gz" \
    --run_ancestry "$PANEL" --min_overlap 0.5 --outdir "${G}/results" \
    --max_cpus "$(nproc)" --max_memory "$(( $(free -g | awk '/^Mem:/ {print $2}') - 1 )).GB" ) > "${CASE_TMP}/genomewide.log" 2>&1
RC=$?
tail -n 40 "${CASE_TMP}/genomewide.log"
T_RUN=$(( $(date +%s) - t0 ))
case "$RC" in
  0) RESULT=finished ;;
  124|137) RESULT="stopped after 35 minutes" ;;
  *) RESULT="failed (exit ${RC})" ;;
esac
check_eq "pgsc_calc with the 1000 Genomes panel on the genome-wide target" "$RESULT" finished
sleep 6
kill "$SAMPLER" 2>/dev/null || true
WORK_BYTES=$(du -sb "${G}/work" 2>/dev/null | cut -f1)
TRACE=$(find "${G}/results/pipeline_info" -name 'execution_trace_*.txt' 2>/dev/null | head -1)
# The tasks with the highest peak RSS, as Nextflow's trace gives it (e.g. "3.2 GB").
TOP=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
  { v = $c["peak_rss"]; n = v + 0; if (v ~ /GB/) n *= 1024; else if (v ~ /KB/) n /= 1024; else if (v !~ /MB/) n /= 1048576
    printf "%.0f\t%s\t%s\n", n, $c["name"], $c["realtime"] }' "${TRACE:-/dev/null}" 2>/dev/null | sort -rn | head -6)
PEAK_TASK_MB=$(head -1 <<< "$TOP" | cut -f1)
PEAK_MEM=$(awk -F'\t' 'NR == 1 || $2 > m {m = $2} END {print m}' "$SAMPLES")
MIN_DISK=$(awk -F'\t' 'NR == 1 || $3 < m {m = $3} END {print m}' "$SAMPLES")
DISK_GB=$(awk -v b="$BASE_DISK" -v m="$MIN_DISK" 'BEGIN {printf "%.1f", (b - m) / 1048576}')
MEM_GB=$(awk -v b="$BASE_MEM" -v p="$PEAK_MEM" 'BEGIN {printf "%.1f", (p - b) / 1024}')
{
  printf 'item\tvalue\n'
  printf 'panel\t%s\n' "$PGSC_PANEL"
  printf 'panel_download_bytes\t%s\n' "$D_BYTES"
  printf 'panel_download_seconds_8_ranges\t%s\n' "$D_SECONDS"
  printf 'panel_grch38_snvs\t%s\n' "$N_SITES"
  printf 'site_list_seconds\t%s\n' "$T_SITES"
  printf 'genomewide_run\t%s\n' "$RESULT"
  printf 'genomewide_run_seconds\t%s\n' "$T_RUN"
  printf 'disk_peak_gb_beyond_the_panel\t%s\n' "$DISK_GB"
  printf 'pgsc_calc_work_bytes_at_end\t%s\n' "${WORK_BYTES:-0}"
  printf 'memory_peak_gb_above_start\t%s\n' "$MEM_GB"
  printf 'task_peak_rss_mb\t%s\n' "${PEAK_TASK_MB:-0}"
  printf 'runner_memory_gb\t%s\n' "$(free -g | awk '/^Mem:/ {print $2}')"
  printf 'runner_cpus\t%s\n' "$(nproc)"
} > "$MEASURE"
cat "$MEASURE"; echo "tasks with the highest peak RSS (MB, name, time):"; echo "$TOP"
{
  printf '#### pgsc_calc with the 1000 Genomes panel (%s), genome-wide target\n\n```\n' "$PGSC_PANEL"
  cat "$MEASURE"; echo "top tasks by peak RSS (MB, name, time):"; echo "$TOP"
  printf '```\n\n'
} >> "$E2E_NOTES"
check_ge "a task's peak RSS was recorded (MB)" "${PEAK_TASK_MB:-0}" 1
finish
