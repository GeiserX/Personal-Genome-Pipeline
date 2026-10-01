#!/usr/bin/env bash
# Step 24's HTML report shows the planted ClinVar hit with its gene, its
# significance and no empty cell.
. "$(dirname "$0")/lib.sh"

run_step 24-html-report.sh "$SAMPLE"
check_step_exit 24-html-report.sh

HTML="${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html"
check "report exists" test -s "$HTML"
REPORT=$(cat "$HTML" 2>/dev/null)
# The planted hit's table row, found by its position. The report writes '.'
# for a missing gene, significance or review status. Without a row the checks
# below read '<td>.</td>' and fail.
ROW=$(grep -o "<tr>.*<td>$(planted pos)</td>.*</tr>" <<< "$REPORT" | awk 'NR == 1')
# SCRATCH, reverted in the next commit: blank the gene cell so the checks below must fail.
ROW=$(sed "s#<td>$(planted gene)</td>#<td>.</td>#" <<< "$ROW")
echo "planted row: ${ROW:-none}"
check "the planted row has no empty '.' cell" lacks '<td>\.</td>' "${ROW:-<td>.</td>}"
check "the planted row names its gene ($(planted gene))" \
  has "<td>$(planted gene)(</td>|[:|,])" "${ROW:-}"
check "the planted row shows its significance (Pathogenic)" has '<td>Pathogenic</td>' "${ROW:-}"

finish
