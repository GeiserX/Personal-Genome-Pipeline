#!/usr/bin/env bash
# Scratch helper (deleted before merge): the same real run against origin/main,
# which is expected to fail in VCFANNO and DELLY.
set -euo pipefail

D=${1:?}; MAIN=${2:?}
cd "$MAIN"
# main requires the .dict on every run
samtools dict "$D/ref.fa" -o "$D/ref.dict"
echo "process.errorStrategy = 'finish'" > finish.config

if nextflow run main.nf -profile docker -c "$D/ci.config" -c finish.config \
    --input "$D/samplesheet.csv" --reference "$D/ref.fa" --outdir results_main \
    --tools vcfanno,delly \
    --cadd_snv "$D/cadd_tiny.tsv.gz" --cadd_snv_index "$D/cadd_tiny.tsv.gz.tbi" > main.log 2>&1; then
  echo "UNEXPECTED: main's real run succeeded"
  tail -30 main.log
  exit 1
fi
echo "main's real run failed, as expected. Failing tasks:"
grep -E "Process .* terminated|Error executing process" main.log | sort -u
for f in work/*/*/.command.err; do
  if grep -q 'command not found' "$f"; then
    echo "--- $(head -c 400 "$(dirname "$f")/.command.sh" | grep -m1 -oE 'vcfanno|delly call|bcftools annotate' || true) :: $f"
    grep 'command not found' "$f"
  fi
done
n=$(cat work/*/*/.command.err | grep -c 'command not found' || true)
echo "command-not-found lines: ${n}"
test "$n" -gt 0
