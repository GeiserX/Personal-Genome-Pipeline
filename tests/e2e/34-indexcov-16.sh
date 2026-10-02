#!/usr/bin/env bash
# Step 16 (goleft indexcov) writes a .ped row for the sample. The slices are too
# small for a trusted sex call, so this only checks the row is there and prints
# it (column 5 is the inferred sex).
. "$(dirname "$0")/lib.sh"

SEX_CHECK=warn run_step 16-indexcov.sh "$SAMPLE" male
check_step_exit 16-indexcov.sh

PED="${GENOME_DIR}/${SAMPLE}/indexcov/indexcov-indexcov.ped"
ROW=$(awk '!/^#/ && NF >= 6' "$PED" 2>/dev/null | awk 'NR == 1')
echo "ped header: $(awk 'NR == 1' "$PED" 2>/dev/null)"
echo "ped row:    ${ROW:-none}"
check "indexcov .ped has a sample row with at least 6 columns" test -n "${ROW:-}"

finish
