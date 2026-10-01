#!/usr/bin/env bash
# Scratch helper (deleted before merge): the input-check controls on
# origin/main, to show each new check would have caught something there.
set -uo pipefail

D=${1:?}; MAIN=${2:?}
cd "$MAIN" || exit 1
test -e "$D/ref.dict" || samtools dict "$D/ref.fa" -o "$D/ref.dict"
common=(-profile docker -c "$D/ci.config" --reference "$D/ref.fa")

run() {
  local label=$1; shift
  if nextflow run main.nf "$@" > ctl.log 2>&1; then
    echo "main, ${label}: run SUCCEEDED"
  else
    echo "main, ${label}: run failed with:"
    grep -m3 -E 'ERROR|Error|error' ctl.log | sed 's/^/    /'
  fi
}

run "--tools roh,clinvar_screen" "${common[@]}" --input "$D/samplesheet.csv" --outdir m1 --tools roh,clinvar_screen
run "duplicate sample id" "${common[@]}" --input "$D/samplesheet_dup.csv" --outdir m2 --tools roh
run "--tools cnvpytor without resources" "${common[@]}" --input "$D/samplesheet.csv" --outdir m3 --tools cnvpytor
run "--tools manta,survivor_merge" "${common[@]}" --input "$D/samplesheet.csv" --outdir m4 --tools manta,survivor_merge
echo "    consensus records: $(bcftools view -H m4/HG002/sv_merged/HG002_sv_consensus.vcf.gz 2>/dev/null | wc -l)"
run "VCF with FILTER '.' everywhere" "${common[@]}" --input "$D/samplesheet_nofilter.csv" --outdir m5 --tools roh
run "expansion_hunter, samplesheet with sex=male" "${common[@]}" --input "$D/samplesheet.csv" --outdir m6 --tools expansion_hunter --expansion_catalog "$D/eh_catalog_chr20.json"
echo "    ExpansionHunter command on main:"
grep -h -A7 '^ExpansionHunter' work/*/*/.command.sh | sed 's/^/    /' | head -8
