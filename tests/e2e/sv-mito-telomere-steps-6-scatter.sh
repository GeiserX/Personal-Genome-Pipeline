#!/usr/bin/env bash
# Steps 03a (GATK HaplotypeCaller) and 03b (FreeBayes) on chr20 and chr22 of
# the fixture, scattered (one container per chromosome, two at a time) and in
# one process (SCATTER=false): both runs must give the same records. A copy of
# the sample (the BAM hard-linked) keeps case 22's vcf_gatk/ as it is.
. "$(dirname "$0")/lib.sh"

T="${SAMPLE}sc"
mkdir -p "${GENOME_DIR}/${T}/aligned"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${T}/aligned/${T}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${T}/aligned/${T}_sorted.bam.bai"
records() { bcf query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT]\n' "$1" 2>/dev/null; }

for step in 03a-gatk-haplotypecaller.sh:vcf_gatk 03b-freebayes.sh:vcf_freebayes; do
  s=${step%%:*} dir=${step#*:}
  VCF="${T}/${dir}/${T}.vcf.gz"
  INTERVALS="chr20 chr22" SCATTER_JOBS=2 run_step "$s" "$T"
  check_step_exit "$s (scattered)"
  check "${s}: two units, two at a time" has '2 unit\(s\), 2 at a time' "$(cat "$STEP_LOG")"
  records "$VCF" > "${CASE_TMP}/${dir}-scattered.tsv"
  INTERVALS="chr20 chr22" SCATTER=false run_step "$s" "$T"
  check_step_exit "$s (one process)"
  check "${s}: one unit" has '1 unit\(s\), 1 at a time' "$(cat "$STEP_LOG")"
  records "$VCF" > "${CASE_TMP}/${dir}-single.tsv"
  N=$(grep -c . "${CASE_TMP}/${dir}-single.tsv" || true)
  check_ge "${s}: records on chr20 and chr22" "$N" 200
  check_ge "${s}: of them on chr22" "$(grep -c '^chr22' "${CASE_TMP}/${dir}-single.tsv" || true)" 50
  check "${s}: the scattered run has the same records" cmp -s "${CASE_TMP}/${dir}-scattered.tsv" "${CASE_TMP}/${dir}-single.tsv"
  diff "${CASE_TMP}/${dir}-scattered.tsv" "${CASE_TMP}/${dir}-single.tsv" | head -n 10
done

finish
