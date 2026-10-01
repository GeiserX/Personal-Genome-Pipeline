#!/usr/bin/env bash
# Step 24's HTML report shows the planted ClinVar hit with its gene and no empty cells.
. "$(dirname "$0")/lib.sh"

run_step 24-html-report.sh "$SAMPLE"
check_step_exit 24-html-report.sh

HTML="${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html"
check "report exists" test -s "$HTML"
REPORT=$(cat "$HTML" 2>/dev/null)
check_eq "cells that read '.|.' (no gene, no significance)" "$(grep -c '<td>\.|\.</td>' <<< "$REPORT" || true)" 0
# The planted hit's table row (found by its position) has a cell that starts
# with the gene, whether the gene has its own column or shares one.
ROW=$(grep -o "<tr>.*<td>$(planted pos)</td>.*</tr>" <<< "$REPORT" | awk 'NR == 1')
echo "planted row: ${ROW:-none}"
check "the ClinVar row of the planted hit names its gene ($(planted gene))" \
  has "<td>$(planted gene)(</td>|[:|,])" "${ROW:-}"

finish
