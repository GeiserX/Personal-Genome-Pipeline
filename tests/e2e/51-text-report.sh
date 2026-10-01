#!/usr/bin/env bash
# generate-report.sh writes the text summary; the ClinVar section names the gene.
. "$(dirname "$0")/lib.sh"

run_step generate-report.sh "$SAMPLE"
check_step_exit generate-report.sh

check "no shell arithmetic error" lacks 'integer expression expected' "$(cat "$STEP_LOG")"
TXT="${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt"
check "text report exists" test -s "$TXT"
check "the ClinVar section names the planted gene ($(planted gene))" grep -q "$(planted gene)" "$TXT"

finish
