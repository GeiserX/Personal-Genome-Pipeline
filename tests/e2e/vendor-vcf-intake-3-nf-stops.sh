#!/usr/bin/env bash
# VCF_PRECHECK stops a vendor VCF the steps would misread, before any analysis,
# naming the sample and the fix:
#   - Ensembl contig names (1, MT): without the stop, a run without clinvar
#     ended with exit 0, an empty haplogroup file and a wrong ROH summary. The
#     printed rename command, run as printed, gives the chr slice's ROH and
#     haplogroup outputs (case 2) exactly;
#   - a gVCF (reference blocks, ALT <*> or <NON_REF> with INFO/END) with
#     pharmcat selected: it used to die inside PHARMCAT_PREPROCESS;
#   - a plain VCF named *.g.vcf.gz with pharmcat selected: PharmCAT refuses the
#     name, so the message says to rename the file.
. "$(dirname "$0")/lib.sh"
. "$(dirname "$0")/vendor-vcf-intake.inc"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

# only_precheck: true when VCF_PRECHECK is the only task the run started.
only_precheck() {
  local others
  others=$(trace_names | grep -v 'VCF_PRECHECK' | paste -sd' ' -)
  [ -z "$others" ] || { echo "    tasks besides VCF_PRECHECK: ${others}"; return 1; }
  trace_names | grep -q VCF_PRECHECK
}

# --- Ensembl contig names --------------------------------------------------------
mkdir -p "${INTAKE}/ensembl"
{ for c in $(seq 1 22) X Y; do echo "chr$c $c"; done; echo "chrM MT"; } > "${INTAKE}/to_ensembl.txt"
bcf annotate --no-version --rename-chrs intake/to_ensembl.txt -Oz -o "intake/ensembl/${SAMPLE}.vcf.gz" "${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
bcf index -f -t "intake/ensembl/${SAMPLE}.vcf.gz"
check "the Ensembl copy has no chr-named contig with records" \
  lacks '^chr' "$(bcf index -s "intake/ensembl/${SAMPLE}.vcf.gz" | cut -f1)"

sheet "${CASE_TMP}/ensembl.csv" "${INTAKE}/ensembl/${SAMPLE}.vcf.gz"
nf_run ensembl "${CASE_TMP}/ensembl.csv" roh,mito_haplogroup
LOG=$(cat "$NF_LOG")
check "Ensembl names: the run fails" test "$NF_RC" -ne 0
check "Ensembl names: the message names the sample and the file" \
  has "Sample '${SAMPLE}': no contig in ${SAMPLE}.vcf.gz is named the chr way" "$LOG"
check "Ensembl names: the message prints the rename command" has 'bcftools annotate --rename-chrs chr_map.txt' "$LOG"
check "Ensembl names: no analysis task started" only_precheck
check "Ensembl names: no haplogroup file published" test ! -e "${NF_OUT}/${SAMPLE}/mito"
check "a failed run prints the failure message" has 'Pipeline failed' "$LOG"
check "no onComplete handler error on a failed run" lacks 'Failed to invoke .workflow.onComplete. event handler' "$LOG"

# The printed commands, run as printed in the file's folder.
# (awk: Nextflow may print the error twice; the commands count once)
grep -E '^    (for c in |echo "MT chrM"|bcftools )' "$NF_LOG" | sed 's/^    //' | awk '!seen[$0]++' > "${CASE_TMP}/printed.sh"
echo "printed commands:"; sed 's/^/    /' "${CASE_TMP}/printed.sh"
check_eq "four command lines printed" "$(grep -c . "${CASE_TMP}/printed.sh" || true)" 4
(
  cd "${INTAKE}/ensembl" || exit 1
  eval "$bcftools_fn"
  set -e
  # shellcheck source=/dev/null
  . "${CASE_TMP}/printed.sh"
) > "${CASE_TMP}/printed.log" 2>&1
check_eq "the printed commands exit 0" "$?" 0
cat "${CASE_TMP}/printed.log"
RENAMED="${INTAKE}/ensembl/${SAMPLE}.chr.vcf.gz"
check "the printed commands wrote ${SAMPLE}.chr.vcf.gz and its index" test -s "$RENAMED" -a -s "${RENAMED}.tbi"

sheet "${CASE_TMP}/renamed.csv" "$RENAMED"
nf_run renamed "${CASE_TMP}/renamed.csv" roh,mito_haplogroup
check_eq "renamed: nextflow run exits 0" "$NF_RC" 0
R="${NF_OUT}/${SAMPLE}"
check_eq "renamed: ROH segments match the chr slice's" \
  "$(grep -v '^#' "${R}/roh/${SAMPLE}_roh.txt" 2>/dev/null | md5sum)" \
  "$(grep -v '^#' "${INTAKE}/chr-run/${SAMPLE}_roh.txt" 2>/dev/null | md5sum)"
check_ge "renamed: ROH file has segments" "$(grep -c '^RG' "${R}/roh/${SAMPLE}_roh.txt" 2>/dev/null || true)" 1
check_eq "renamed: haplogroup file matches the chr slice's" \
  "$(md5sum < "${R}/mito/${SAMPLE}_haplogroup.txt" 2>/dev/null)" \
  "$(md5sum < "${INTAKE}/chr-run/${SAMPLE}_haplogroup.txt" 2>/dev/null)"
check_ge "renamed: haplogroup file has a call" "$(awk 'NR > 1' "${R}/mito/${SAMPLE}_haplogroup.txt" 2>/dev/null | grep -c . || true)" 1

# --- gVCF with pharmcat --------------------------------------------------------------
# The chr slice plus two reference blocks, one per common style.
mkdir -p "${INTAKE}/gvcf"
{
  bcf view -h "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | grep '^##' | grep -v '^##INFO=<ID=END,'
  echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the reference block">'
  bcf view -h "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | grep '^#CHROM'
  bcf view -H "${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
  printf 'chr20\t10000001\t.\tN\t<*>\t0\t.\tEND=10000050\tGT\t0/0\n'
  printf 'chr10\t94700001\t.\tN\t<NON_REF>\t.\t.\tEND=94700100\tGT\t0/0\n'
} > "${INTAKE}/gvcf/unsorted.vcf"
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools sort -Oz -o intake/gvcf/${SAMPLE}.g.vcf.gz intake/gvcf/unsorted.vcf
  bcftools index -f -t intake/gvcf/${SAMPLE}.g.vcf.gz"
check_eq "the gVCF has two reference blocks" \
  "$(vcf_count -i 'INFO/END!="."' "intake/gvcf/${SAMPLE}.g.vcf.gz")" 2

sheet "${CASE_TMP}/gvcf.csv" "${INTAKE}/gvcf/${SAMPLE}.g.vcf.gz"
nf_run gvcf "${CASE_TMP}/gvcf.csv" pharmcat,roh
LOG=$(cat "$NF_LOG")
check "gVCF: the run fails" test "$NF_RC" -ne 0
check "gVCF: the message names the sample and says PharmCAT refuses a gVCF" \
  has "Sample '${SAMPLE}': ${SAMPLE}.g.vcf.gz is a gVCF .*PharmCAT refuses a gVCF" "$LOG"
check "gVCF: the message says what to do" has 'Remove the reference blocks' "$LOG"
check "gVCF: stopped in the precheck, before PHARMCAT_PREPROCESS" only_precheck

# --- A plain VCF named like a gVCF ---------------------------------------------------
mkdir -p "${INTAKE}/named"
cp "$CHR_VCF" "${INTAKE}/named/${SAMPLE}.g.vcf.gz"
cp "${CHR_VCF}.tbi" "${INTAKE}/named/${SAMPLE}.g.vcf.gz.tbi"
sheet "${CASE_TMP}/named.csv" "${INTAKE}/named/${SAMPLE}.g.vcf.gz"
nf_run named "${CASE_TMP}/named.csv" pharmcat
LOG=$(cat "$NF_LOG")
check "gVCF name only: the run fails" test "$NF_RC" -ne 0
check "gVCF name only: the message says to rename the file" \
  has "${SAMPLE}.g.vcf.gz has no reference blocks, but its name contains .g.vcf.*Rename the file" "$LOG"
check "gVCF name only: stopped in the precheck" only_precheck

finish
