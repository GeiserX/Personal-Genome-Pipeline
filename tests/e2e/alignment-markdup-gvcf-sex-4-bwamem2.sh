#!/usr/bin/env bash
# Step 02a (BWA-MEM2) on a mini reference: the fixture's regions cut out of the
# fixture reference, about 7 Mb, small enough for bwa-mem2 index on a runner
# (GRCh38 needs about 90 GB of RAM). Checks that the index is built under a
# temporary prefix and renamed whole, and that the BAM is sorted, has its
# duplicates marked and is renamed into place with nothing left behind.
. "$(dirname "$0")/lib.sh"

MINI="${GENOME_DIR}/bwamini"
MREF="${MINI}/mini.fa"
BS="${SAMPLE}B"
rm -rf "$MINI" "${GENOME_DIR:?}/${BS}"
mkdir -p "$MINI" "${GENOME_DIR}/${BS}/fastq"

# One contig per fixture region, named chr_start_end (1-based, inclusive).
mapfile -t REGS < <(awk '{printf "%s:%d-%d\n", $1, $2 + 1, $3}' "${FIXTURE_DIR}/regions.bed")
sam faidx reference/Homo_sapiens_assembly38.fasta "${REGS[@]}" | sed -E '/^>/ s/[:-]/_/g' > "$MREF"
check_eq "mini reference contigs (one per fixture region)" "$(grep -c '^>' "$MREF" || true)" "${#REGS[@]}"
check_ge "mini reference bases" "$(grep -v '^>' "$MREF" | tr -d '\n' | wc -c | tr -d ' ')" 5000000

for r in R1 R2; do
  ln -f "${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_${r}.fastq.gz" "${GENOME_DIR}/${BS}/fastq/${BS}_${r}.fastq.gz"
done

REF_FASTA="$MREF" THREADS=2 run_step 02a-alignment-bwamem2.sh "$BS"
check_step_exit 02a-alignment-bwamem2.sh
check "step 02a built the index" grep -q '^BWA-MEM2 index built\.' "$STEP_LOG"

for e in 0123 amb ann pac bwt.2bit.64; do
  check "index file mini.fa.${e} is next to the FASTA" nonempty "bwamini/mini.fa.${e}"
done
check_eq "temporary index files left next to the FASTA" \
  "$(find "$MINI" -maxdepth 1 -name 'mini.fa.tmp.*' | wc -l | tr -d ' ')" 0

D="${BS}/aligned_bwamem2"
BAM="${D}/${BS}_sorted.bam"
check "BAM passes samtools quickcheck" sam quickcheck -v "$BAM"
check "BAM index exists" nonempty "${BAM}.bai"
FS=$(sam flagstat "$BAM" 2>/dev/null)
echo "$FS"
TOTAL=$(awk '/ in total / {print $1; exit}' <<< "$FS")
MAPPED=$(awk '/ mapped \(/ {print $1; exit}' <<< "$FS")
DUPS=$(awk '/ duplicates$/ {print $1; exit}' <<< "$FS")
check_ge "reads in the BAM" "${TOTAL:-0}" 100000
check_ge "mapped reads, percent of all" "$(( ${MAPPED:-0} * 100 / ${TOTAL:-1} ))" 50
check_ge "reads flagged as duplicates (samtools flagstat)" "${DUPS:-0}" 1
HDR=$(sam view -H "$BAM" 2>/dev/null)
check "the BAM header records samtools markdup" has '^@PG.*ID:samtools.*markdup' "$HDR"
check "the BAM header records bwa-mem2" has '^@PG.*ID:bwa-mem2' "$HDR"
TAB=$'	'
check "@RG SM equals the sample name (${BS})" has "^@RG.*${TAB}SM:${BS}(${TAB}|\$)" "$HDR"
check "the BAM is coordinate-sorted" has '^@HD.*SO:coordinate' "$HDR"
check_eq "files in aligned_bwamem2/ other than the BAM and its index" \
  "$(find "${GENOME_DIR}/${D}" -mindepth 1 ! -name "${BS}_sorted.bam" ! -name "${BS}_sorted.bam.bai" | wc -l | tr -d ' ')" 0

rm -rf "$MINI" "${GENOME_DIR:?}/${BS}"
finish
