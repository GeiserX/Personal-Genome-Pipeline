#!/usr/bin/env bash
# generate-report.sh three times on one sample: what changed since the
# previous run (bin/run_delta.py) and the history it keeps.
#   1. the first run writes "no previous run" and keeps its summary.json and
#      run_manifest.tsv as history/<UTC time>/;
#   2. a second run with the same outputs writes an empty delta;
#   3. with one ClinVar hit reclassified and one image digest changed in the
#      newest history copy, the delta names exactly those two changes.
# The report renderer and run_delta.py run with the host's python3 on the host
# side of the call's mounts, so both are the real code.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
S="${GENOME_DIR}/sample1"
mkdir -p "${S}/clinvar"
REV='criteria_provided,_multiple_submitters,_no_conflicts'
printf '%s\n' '##fileformat=VCFv4.2' \
  "$(printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1')" \
  "$(printf 'chr1\t1000\t.\tA\tG\t.\tPASS\tCLNSIG=Pathogenic;CLNREVSTAT=%s;GENEINFO=GENEA:1\tGT\t0/1' "$REV")" \
  "$(printf 'chr2\t2000\t.\tC\tT\t.\tPASS\tCLNSIG=Likely_pathogenic;CLNREVSTAT=%s;GENEINFO=GENEB:2\tGT\t0/1' "$REV")" \
  > "${S}/clinvar/sample1_clinvar_hits.vcf"

use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
case "${*:2}" in
  "python3 /pgp-bin/"*)
    # The real script, on the host paths of every /... argument.
    shift
    mapped=()
    for a in "$@"; do
      if [[ "$a" == /* ]] && h=$(host_path "$a"); then mapped+=("$h"); else mapped+=("$a"); fi
    done
    exec "${mapped[@]}" ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

DELTA="${S}/sample1_delta"
history_dirs() { find "${S}/history" -mindepth 1 -maxdepth 1 -type d | LC_ALL=C sort; }
changes() { python3 -c 'import json, sys; print(json.dumps(json.load(open(sys.argv[1]))["changes"]))' "${DELTA}.json"; }

# 1. no previous run
run_expect 0 first "${SCRIPTS}/generate-report.sh" sample1
for ext in txt html json; do [ -s "${DELTA}.${ext}" ] || fail "the first run wrote no ${DELTA}.${ext}"; done
grep -q 'No previous run of this sample on record' "${DELTA}.txt" || fail "the first delta does not say there is no previous run: $(cat "${DELTA}.txt")"
[ "$(history_dirs | wc -l)" -eq 1 ] || fail "after one run history/ holds $(history_dirs | wc -l) folders, expected 1"
first=$(history_dirs | tail -n 1)
[ -s "${first}/summary.json" ] && [ -s "${first}/run_manifest.tsv" ] \
  || fail "history/$(basename "$first") lacks summary.json or run_manifest.tsv: $(ls "$first")"
output_has first 'Delta saved to: .*sample1_delta.txt \(this run kept as history/'

# 2. the same outputs again: no change
run_expect 0 second "${SCRIPTS}/generate-report.sh" sample1
[ "$(changes)" = "[]" ] || fail "two runs on the same outputs list changes: $(changes)"
grep -q '^No change: the same inputs and the same results' "${DELTA}.txt" || fail "the second delta does not say no change: $(cat "${DELTA}.txt")"
[ "$(history_dirs | wc -l)" -eq 2 ] || fail "after two runs history/ holds $(history_dirs | wc -l) folders, expected 2"

# 3. plant a reclassification and a new digest in the previous run
prev=$(history_dirs | tail -n 1)
python3 - "${prev}/summary.json" <<'PY'
import json, sys
p = sys.argv[1]
s = json.load(open(p))
hit = next(h for h in s["sections"]["clinvar"]["data"]["hits"] if h["pos"] == 2000)
hit["significance"] = "Uncertain significance"
json.dump(s, open(p, "w"))
PY
PY_REF=$(awk -F'\t' '$1 == "version" && $2 == "PYTHON_IMAGE" {print $3}' "${prev}/run_manifest.tsv")
awk -F'\t' -v OFS='\t' -v ref="$PY_REF" '$1 == "image" && $2 == ref {$3 = "python@sha256:0000"} {print}' \
  "${prev}/run_manifest.tsv" > "${CASE_WORK}/m.tsv" && mv "${CASE_WORK}/m.tsv" "${prev}/run_manifest.tsv"
grep -q 'python@sha256:0000' "${prev}/run_manifest.tsv" || fail "could not plant the digest of ${PY_REF}"
run_expect 0 third "${SCRIPTS}/generate-report.sh" sample1
got=$(python3 -c 'import json, sys
for c in json.load(open(sys.argv[1]))["changes"]:
    print(c["area"], c["kind"], c["item"], sep="|")' "${DELTA}.json" | LC_ALL=C sort)
want=$(printf '%s\n' 'ClinVar|reclassified|chr2:2000 C>T (GENEB)' 'inputs|input|PYTHON_IMAGE' | LC_ALL=C sort)
[ "$got" = "$want" ] || fail "the delta should name exactly the two planted changes; got: ${got//$'\n'/ ; }"
grep -q 'reclassified: chr2:2000 C>T (GENEB): Uncertain significance -> Likely pathogenic' "${DELTA}.txt" \
  || fail "the text delta does not show the reclassification: $(cat "${DELTA}.txt")"
grep -q 'python@sha256:0000' "${DELTA}.html" || fail "the HTML delta does not show the old digest"
[ "$(history_dirs | wc -l)" -eq 3 ] || fail "after three runs history/ holds $(history_dirs | wc -l) folders, expected 3"
echo "The delta says no previous run, then no change, then exactly the planted reclassification and digest."
