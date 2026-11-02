#!/usr/bin/env bash
# The live site (bin/pgp_live.py) on a real Nextflow run: the pinned Nextflow
# loads the nf-weblog plugin from the config run-all.sh writes when
# PGP_LIVE_URL is set (the same text, the pin read from run-all.sh), posts each
# task to the site's ingest port, and the site stores the messages and shows
# step 11 (ROH) done in state.json. A control reads state.json before the run:
# step 11 is 'not planned' then, so the state check can fail.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
PLUGIN=$(sed -n 's/^NF_WEBLOG_PLUGIN=\([^ ]*\)$/\1/p' "${REPO}/scripts/run-all.sh")
check "run-all.sh pins nf-weblog with a version" has '^nf-weblog@[0-9]+\.[0-9]+\.[0-9]+$' "$PLUGIN"

G="$GENOME_DIR"
ROOT="${CASE_TMP}/root" STATE="${CASE_TMP}/state"
mkdir -p "${ROOT}/${SAMPLE}/nextflow"   # a run folder the site lists (it has nextflow/)
port() { python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1", 0)); print(s.getsockname()[1])'; }
WEB=$(port) ING=$(port)
python3 "${REPO}/bin/pgp_live.py" --root "$ROOT" --state-dir "$STATE" --port "$WEB" --ingest-port "$ING" --quiet \
  2> "${CASE_TMP}/site.err" &
SITE=$!
trap 'kill "$SITE" 2>/dev/null' EXIT
for _ in $(seq 50); do curl -s -o /dev/null "http://127.0.0.1:${WEB}/" && break; sleep 0.2; done
state11() { curl -s "http://127.0.0.1:${WEB}/runs/${SAMPLE}/state.json" \
  | python3 -c 'import json, sys; print(next(s["state"] for s in json.load(sys.stdin)["steps"] if s["step"] == "11"))' 2>/dev/null \
  || echo "no state.json"; }
check_eq "control: before the run, step 11 is not planned" "$(state11)" "not planned"

SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$SAMPLE" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
  "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" > "$SHEET"
printf "plugins {\n    id '%s'\n}\nweblog {\n    enabled = true\n    url     = '%s'\n}\n" "$PLUGIN" \
  "http://127.0.0.1:${ING}/ingest/${SAMPLE}" > "${CASE_TMP}/live.config"
echo "+ nextflow run main.nf --tools roh -c live.config"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -c "${CASE_TMP}/live.config" \
    -work-dir "${CASE_TMP}/work" --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools roh --outdir "${CASE_TMP}/out" --max_cpus 2 --max_memory 6.GB
) 2>&1 | tee "$STEP_LOG"
check_eq "nextflow run exits 0" "${PIPESTATUS[0]}" 0
check "the run loaded the plugin" grep -q "nf-weblog" "${CASE_TMP}/.nextflow.log"

EV="${STATE}/runs/${SAMPLE}/events.jsonl"
count() { python3 -c 'import json, sys
n = 0
for l in open(sys.argv[1]):
    e = json.loads(l); t = e.get("trace") or {}
    n += e.get("event") == sys.argv[2] and (sys.argv[3] == "" or (t.get("process") or "").endswith(sys.argv[3]))
print(n)' "$EV" "$@" 2>/dev/null || echo 0; }
check_ge "the site stored the run's 'started' message" "$(count started '')" 1
check_ge "the site stored ROH's process_completed" "$(count process_completed ROH)" 1
check_ge "the site stored the run's 'completed' message" "$(count completed '')" 1
check "no stored message holds the samplesheet or the VCF path" lacks "${SHEET}|${SAMPLE}\\.vcf\\.gz" "$(cat "$EV" 2>/dev/null)"
check_eq "state.json shows step 11 done" "$(state11)" "done"
finish
