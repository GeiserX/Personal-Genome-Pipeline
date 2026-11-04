#!/usr/bin/env bash
# bin/owngenome_recheck.py --rerun replays the last run-all.sh command exactly:
#   1. run-all.sh writes 'meta command' (the switches it reads that are set,
#      then its arguments, shell-quoted) and 'meta switches' into
#      logs/run_status.tsv;
#   2. --rerun, from a shell that sets another TOOLS and THREADS, starts
#      nextflow with the same arguments, the same launch directory and the
#      same NXF_VER as the run did, and records the rerun in recheck/state.json
#      and one more row in recheck/log.tsv (the log only grows);
#   3. control: with one argument dropped from the 'meta command' row the
#      replayed nextflow line differs, so check 2 can fail;
#   4. --rerun refuses (exit 3) while a run is going on (no 'meta finished_utc')
#      and starts no nextflow.
# The upstream answers come from an offline file: nothing here asks the network.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
S="${G}/sample1"
ST="${S}/logs/run_status.tsv"
seed_reference "$G"
seed_clinvar "$G"
use_output_hook   # the reports at the end read the files their containers write
seed_sample "$G" sample1
# A ClinVar file with its release date, as step 06 reads it.
printf '##fileformat=VCFv4.1\n##fileDate=20260901\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
  | gzip -c > "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz"
ANSWERS="${REPO_ROOT}/tests/fixtures/recheck/offline-nothing-newer.json"
RECHECK=(python3 "${REPO_ROOT}/bin/owngenome_recheck.py" --sample sample1 --offline "$ANSWERS")
nf_lines() { grep '^nextflow :: .*/main\.nf ' "$FAKE_DOCKER_LOG" || true; }
meta() { awk -F'\t' -v k="$1" '$1 == "meta" && $2 == k {print $3}' "$ST"; }

# 1. The run, with switches and arguments that need quoting
run_expect 0 run env TOOLS=pharmcat,cpic THREADS=3 INTERVALS="chr20 chr22" "${SCRIPTS}/run-all.sh" sample1 male --sex_check warn
first=$(nf_lines | tail -1)
[ -n "$first" ] || fail "run-all.sh started no nextflow"
cmd=$(meta command)
[ -n "$cmd" ] || fail "no 'meta command' row: $(cat "$ST")"
for want in "TOOLS=pharmcat\\,cpic" "THREADS=3" "INTERVALS=chr20\\ chr22" "SKIP_VALIDATION=true" "/run-all.sh sample1 male --sex_check warn"; do
  [[ "$cmd" == *"$want"* ]] || fail "the command row lacks '${want}': ${cmd}"
done
[[ " $(meta switches) " == *" TOOLS "* && " $(meta switches) " == *" NXF_VER "* ]] || fail "'meta switches' does not name the switches: $(meta switches)"
[ -s "${S}/run_manifest.tsv" ] || fail "run-all.sh wrote no run manifest"

# 2. The replay, from a shell with other switches set
run_expect 0 rerun env TOOLS=vep THREADS=9 "${RECHECK[@]}" --rerun
second=$(nf_lines | tail -1)
[ "$(nf_lines | wc -l)" -eq 2 ] || fail "want two nextflow runs, got: $(nf_lines)"
[ "$second" = "$first" ] || fail "the rerun did not replay the run: was ${first}, now ${second}"
output_has rerun '^Rerunning: '
python3 - "${S}/recheck/state.json" <<'PY' || fail "recheck/state.json has no rerun record"
import json, sys
s = json.load(open(sys.argv[1]))
r = s["rerun"]
assert r["rc"] == 0 and r["command"] and r["finished_utc"] and "predicted_vs_redone" in r, r
assert all(i["class"] == "no-newer-release" for i in s["items"]), s["items"]
PY
rows1=$(wc -l < "${S}/recheck/log.tsv")
for e in rerun-start rerun-end; do grep -q $'\t'"$e"$'\t' "${S}/recheck/log.tsv" || fail "log.tsv lacks the ${e} row"; done

# 3. Control: a command row with one argument dropped replays something else
awk -F'\t' -v OFS='\t' '$1 == "meta" && $2 == "command" {sub(/ --sex_check warn$/, "", $3)} {print}' "$ST" > "${ST}.new" && mv "${ST}.new" "$ST"
[[ "$(meta command)" != *"--sex_check"* ]] || fail "the control did not drop the argument"
run_expect 0 dropped "${RECHECK[@]}" --rerun
third=$(nf_lines | tail -1)
[ "$third" != "$first" ] || fail "control: a command row without --sex_check replayed the same line, so check 2 cannot fail"
[ "$(wc -l < "${S}/recheck/log.tsv")" -gt "$rows1" ] || fail "log.tsv did not grow"

# 4. A run going on: no rerun
grep -v $'^meta\tfinished_utc\t' "$ST" > "${ST}.new" && mv "${ST}.new" "$ST"
n=$(nf_lines | wc -l)
run_expect 3 running "${RECHECK[@]}" --rerun
output_has running 'is going on'
[ "$(nf_lines | wc -l)" -eq "$n" ] || fail "--rerun started nextflow while a run was going on"
echo "OK: --rerun replays the recorded command and refuses during a run"
