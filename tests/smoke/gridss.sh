# shellcheck shell=sh
# GRIDSS_IMAGE row: GRIDSS indexes the reference beside the FASTA (step 04b
# mounts the reference directory writable for this), so it gets its own copy.
set -e
mkdir -p ref
cp /in/mini.fa /in/mini.fa.fai /in/mini.dict ref/
gridss -r ref/mini.fa -o gridss.vcf.gz -a assembly.bam -t 4 --jvmheap 6g /in/mini.bam
