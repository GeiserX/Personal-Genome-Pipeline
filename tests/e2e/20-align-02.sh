#!/usr/bin/env bash
# Step 02 aligns the fixture reads; GATK and DeepVariant need a read group.
. "$(dirname "$0")/lib.sh"

run_step 02-alignment.sh "$SAMPLE"
check_step_exit 02-alignment.sh

BAM="${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
check "BAM passes samtools quickcheck" sam quickcheck -v "$BAM"
check "BAM index exists" nonempty "${BAM}.bai"
HDR=$(sam view -H "$BAM")
check "header has an @RG line" has '^@RG' "$HDR"
TAB=$'	'
check "@RG SM equals the sample name (${SAMPLE})" has "^@RG.*${TAB}SM:${SAMPLE}(${TAB}|\$)" "$HDR"
IDX=$(sam idxstats "$BAM")
for c in chr1 chr2 chr4 chr5 chr6 chr10 chr12 chr16 chr19 chr20 chr22 chrX chrY chrM; do
  check_ge "mapped reads on ${c}" "$(awk -v c="$c" '$1 == c {print $3}' <<< "$IDX")" 100
done

finish
