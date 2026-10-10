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
#     name, so the message says to rename the file;
#   - two samples in one VCF, merged with real bcftools: the steps would mix
#     their genotypes;
#   - a ##contig chr1 length that is GRCh37's: every position would be read on
#     the wrong build.
# And validate-setup.sh, on a VCF whose header has no ##contig lines: its REF
# spot-check runs real bcftools norm against the reference, passes the
# fixture VCF and fails a copy whose positions are shifted off its REF bases.
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

# --- Two samples in one VCF ---------------------------------------------------------
# The fixture call plus a copy of it under another sample name, merged.
mkdir -p "${INTAKE}/joint"
FIRST=$(bcf query -l "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | head -n 1)
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  echo '${SAMPLE}_B' > intake/joint/rename.txt
  bcftools reheader -s intake/joint/rename.txt -o intake/joint/b.vcf.gz ${SAMPLE}/vcf/${SAMPLE}.vcf.gz
  bcftools index -f -t intake/joint/b.vcf.gz
  bcftools merge -Oz -o intake/joint/${SAMPLE}.vcf.gz ${SAMPLE}/vcf/${SAMPLE}.vcf.gz intake/joint/b.vcf.gz
  bcftools index -f -t intake/joint/${SAMPLE}.vcf.gz"
check_eq "the merged VCF holds two samples" \
  "$(bcf query -l "intake/joint/${SAMPLE}.vcf.gz" | paste -sd, -)" "${FIRST},${SAMPLE}_B"

sheet "${CASE_TMP}/joint.csv" "${INTAKE}/joint/${SAMPLE}.vcf.gz"
nf_run joint "${CASE_TMP}/joint.csv" roh,mito_haplogroup
LOG=$(cat "$NF_LOG")
check "two samples: the run fails" test "$NF_RC" -ne 0
check "two samples: the message names the sample, the count and both names" \
  has "Sample '${SAMPLE}': ${SAMPLE}\.vcf\.gz holds 2 samples \(${FIRST}, ${SAMPLE}_B\)" "$LOG"
check "two samples: the message prints the bcftools command that keeps one" \
  has "bcftools view -s ${FIRST} -a -c 1 -Oz -o ${FIRST}\.vcf\.gz ${SAMPLE}\.vcf\.gz" "$LOG"
check "two samples: no analysis task started" only_precheck

# --- A GRCh37 chr1 length in the header ------------------------------------------------
# The fixture call with only its chr1 ##contig line changed to GRCh37's length.
mkdir -p "${INTAKE}/grch37"
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools view -h ${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    | sed 's/^##contig=<ID=chr1,length=248956422/##contig=<ID=chr1,length=249250621/' > intake/grch37/header.txt
  bcftools reheader -h intake/grch37/header.txt -o intake/grch37/${SAMPLE}.vcf.gz ${SAMPLE}/vcf/${SAMPLE}.vcf.gz
  bcftools index -f -t intake/grch37/${SAMPLE}.vcf.gz"
check "the copy's header gives GRCh37's chr1 length" \
  has '^##contig=<ID=chr1,length=249250621' "$(bcf view -h "intake/grch37/${SAMPLE}.vcf.gz")"

sheet "${CASE_TMP}/grch37.csv" "${INTAKE}/grch37/${SAMPLE}.vcf.gz"
nf_run grch37 "${CASE_TMP}/grch37.csv" roh,mito_haplogroup
LOG=$(cat "$NF_LOG")
check "GRCh37 length: the run fails" test "$NF_RC" -ne 0
check "GRCh37 length: the message names the sample and the build" \
  has "Sample '${SAMPLE}': ${SAMPLE}\.vcf\.gz is not on GRCh38" "$LOG"
check "GRCh37 length: the message gives both chr1 lengths" has 'chr1 length 249250621 \(GRCh38: 248956422\)' "$LOG"
check "GRCh37 length: no analysis task started" only_precheck
check "GRCh37 length: no ROH file published" test ! -e "${NF_OUT}/${SAMPLE}/roh"

# --- validate-setup.sh on VCFs without ##contig lines ----------------------------------
# nocontig38: the fixture call with its ##contig lines removed (same records).
# nocontig_shift: the same, every position moved 7 bases, so most REF bases
# no longer match the reference, as on another build.
for s in nocontig38 nocontig_shift; do mkdir -p "${G}/${s}/vcf"; done
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools view -h ${SAMPLE}/vcf/${SAMPLE}.vcf.gz | grep -v '^##contig=' > intake/nocontig_header.txt
  bcftools reheader -h intake/nocontig_header.txt -o nocontig38/vcf/nocontig38.vcf.gz ${SAMPLE}/vcf/${SAMPLE}.vcf.gz
  bcftools index -f -t nocontig38/vcf/nocontig38.vcf.gz
  bcftools view ${SAMPLE}/vcf/${SAMPLE}.vcf.gz | awk -F'\t' -v OFS='\t' '/^#/ { print; next } { \$2 = \$2 + 7; print }' \
    | bcftools view -Oz -o intake/shifted.vcf.gz -
  bcftools reheader -h intake/nocontig_header.txt -o nocontig_shift/vcf/nocontig_shift.vcf.gz intake/shifted.vcf.gz
  bcftools index -f -t nocontig_shift/vcf/nocontig_shift.vcf.gz"
check "nocontig38 has no ##contig line" lacks '^##contig' "$(bcf view -h nocontig38/vcf/nocontig38.vcf.gz 2>/dev/null)"
check_ge "nocontig38 has records" "$(vcf_count nocontig38/vcf/nocontig38.vcf.gz)" 100

for s in nocontig38 nocontig_shift; do
  echo "+ scripts/validate-setup.sh ${s}"
  "${REPO}/scripts/validate-setup.sh" "$s" > "${CASE_TMP}/validate-${s}.log" 2>&1
  grep -E 'VCF|##contig' "${CASE_TMP}/validate-${s}.log" | sed 's/^/    | /'
done
V38=$(cat "${CASE_TMP}/validate-nocontig38.log")
VSH=$(cat "${CASE_TMP}/validate-nocontig_shift.log")
check "no ##contig lines: validate-setup says the header does not show the build" \
  has '\[WARN\].*VCF header has no ##contig lines' "$V38"
check "no ##contig lines, GRCh38 records: the REF spot-check passes (real bcftools norm)" \
  has '\[OK\].*VCF REF bases match the reference in all of the first [0-9]+ records' "$V38"
check "no ##contig lines, shifted records: the REF spot-check fails" \
  has '\[FAIL\].*VCF REF bases differ from the reference in [0-9]+ of the first [0-9]+ records: the VCF is not on GRCh38' "$VSH"
check "the spot-check ran on both (never a 'could not compare')" \
  lacks 'Could not compare the VCF' "${V38}${VSH}"

finish
