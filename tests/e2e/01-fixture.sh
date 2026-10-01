#!/usr/bin/env bash
# The downloaded fixture is whole: the BAM reads, every .gz file decompresses,
# every sliced contig has reads, and the VEP subset has CSQ annotations.
. "$(dirname "$0")/lib.sh"

fx() { docker run --rm -v "${FIXTURE_DIR}:/f:ro" -w /f "$SAMTOOLS_IMAGE" samtools "$@"; }

check "HG002_slice.bam passes samtools quickcheck" fx quickcheck -v HG002_slice.bam
for f in "${FIXTURE_DIR}"/*.gz; do
  check "gzip -t $(basename "$f")" gzip -t "$f"
done

IDX=$(fx idxstats HG002_slice.bam)
for c in chr1 chr2 chr4 chr5 chr6 chr10 chr12 chr16 chr19 chr20 chr22 chrX chrY chrM; do
  check_ge "reads on ${c} in HG002_slice.bam" "$(awk -v c="$c" '$1 == c {print $3}' <<< "$IDX")" 1
done

R1=$(( $(gzip -dc "${FIXTURE_DIR}/${SAMPLE}_R1.fastq.gz" | wc -l) / 4 ))
R2=$(( $(gzip -dc "${FIXTURE_DIR}/${SAMPLE}_R2.fastq.gz" | wc -l) / 4 ))
check_eq "R2 holds as many reads as R1" "$R2" "$R1"
check_ge "read pairs in the FASTQ" "$R1" 100000

VEP_N=$(grep -vc '^#' "${FIXTURE_DIR}/${SAMPLE}_vep.vcf" || true)
check_ge "records in ${SAMPLE}_vep.vcf" "$VEP_N" 20
check "at most 200 records in ${SAMPLE}_vep.vcf" test "$VEP_N" -le 200
check "${SAMPLE}_vep.vcf declares CSQ" grep -q '^##INFO=<ID=CSQ' "${FIXTURE_DIR}/${SAMPLE}_vep.vcf"

TOTAL=$(find "$FIXTURE_DIR" -maxdepth 1 -type f -printf '%s\n' | awk '{s += $1} END {print s}')
check "fixture is under 1.5 GB (${TOTAL} bytes)" test "$TOTAL" -lt $((1500 * 1024 * 1024))

finish
