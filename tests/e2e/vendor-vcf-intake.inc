# shellcheck shell=bash
# vendor-vcf-intake.inc — helpers for the vendor-vcf-intake-*.sh cases, sourced
# after lib.sh. Not a case itself: the runner only runs *.sh files.

G="$GENOME_DIR"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
CHR_VCF="${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"   # case 21's DeepVariant call on the slices
INTAKE="${G}/intake"                              # inputs these cases build
mkdir -p "$INTAKE"

# Every Nextflow run in these cases writes a trace with start and end times,
# so a case can tell which task waited for which.
cat > "${CASE_TMP}/intake.config" <<'NFCONF'
trace.fields = 'task_id,name,status,exit,submit,start,complete'
NFCONF

# nf_run NAME SHEET TOOLS [ARGS...]: one Nextflow run, outputs in
# GENOME_DIR/nf-NAME, console log in NF_LOG, exit code in NF_RC, trace in NF_TRACE.
nf_run() {
  local name=$1 sheet=$2 tools=$3
  shift 3
  NF_OUT="${G}/nf-${name}"
  NF_LOG="${CASE_TMP}/${name}.log"
  mkdir -p "${CASE_TMP}/${name}"
  echo "+ nextflow run main.nf (${name}): --tools ${tools} $*"
  (
    cd "${CASE_TMP}/${name}" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -c "${CASE_TMP}/intake.config" \
      -work-dir "${E2E_WORK}/nf-work-intake" \
      --input "$sheet" --reference "$REF" --tools "$tools" \
      --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB "$@"
  ) > "$NF_LOG" 2>&1
  NF_RC=$?
  cp "${CASE_TMP}/${name}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.${name}.nextflow.log" 2>/dev/null
  NF_TRACE=$(find "${NF_OUT}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
  sed 's/^/    | /' "$NF_LOG" | grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists'
  echo "+ exit ${NF_RC}"
}

# sheet FILE VCF: a VCF-only samplesheet for SAMPLE.
sheet() {
  printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$SAMPLE" "$2" "$2.tbi" > "$1"
}

# trace_col PROCESS COLUMN: COLUMN of the trace row whose name contains
# PROCESS (first match), empty when there is none.
trace_col() {
  awk -F'\t' -v p="$1" -v k="$2" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    index($c["name"], p) {print $c[k]; exit}' "${NF_TRACE:-/dev/null}" 2>/dev/null
}

# trace_names: the process name of every task in the trace, one per line.
trace_names() {
  awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} {print $c["name"]}' \
    "${NF_TRACE:-/dev/null}" 2>/dev/null
}

# json_calls FILE: PharmCAT's calls in a match.json or phenotype.json, without
# the run's own metadata (time, version, file names), as sorted canonical JSON.
# Two runs that called the same thing print the same text.
json_calls() {
  python3 - "$1" <<'PY'
import json, re, sys
try:
    data = json.load(open(sys.argv[1]))
except Exception as e:
    print(f"unreadable: {e}")
    sys.exit(0)
SKIP = re.compile(r"(?i)(metadata|timestamp|date|time|file|path|version|title)")
def strip(x):
    if isinstance(x, dict):
        return {k: strip(v) for k, v in x.items() if not SKIP.search(k)}
    if isinstance(x, list):
        return [strip(v) for v in x]
    return x
print(json.dumps(strip(data), sort_keys=True))
PY
}

# html_stat FILE HEADING LABEL: the value of the stat LABEL in the card titled HEADING.
html_stat() {
  awk -v h="<h2>$2</h2>" -v l="<span class=\"label\">$3</span>" '
    index($0, h) {card = 1; next}
    card && /<h2>/ {exit}
    card && index($0, l) {
      line = $0
      if (line !~ /class="value"/) { getline nxt; line = line nxt }
      sub(/.*class="value">/, "", line); gsub(/<[^>]*>/, "", line); sub(/[ \t]+$/, "", line)
      print line; exit
    }' "$1" 2>/dev/null
}

# bcftools_here: bcftools from the pinned image on the current directory, as
# docs/vcf-first.md defines it, so the printed and documented commands run as written.
bcftools_fn='bcftools() { docker run --rm -i -u "$(id -u):$(id -g)" -v "$PWD:$PWD" -w "$PWD" "$BCFTOOLS_IMAGE" bcftools "$@"; }'
