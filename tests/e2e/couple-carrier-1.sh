#!/usr/bin/env bash
# Step 39 and the Nextflow COUPLE_CARRIER process on a couple: HG002 (the
# fixture's sample, his DeepVariant calls from the base cases) and a synthetic
# female partner, PARTNER, in a GENOME_DIR of their own. ClinVar is the
# fixture's file plus four planted pathogenic records (2 stars): two at real
# heterozygous SNVs of HG002, in NAGA and in VDR, and two at the partner's
# SNVs, in NAGA (another position) and in UGT1A1. The partner's Delly calls
# hold one heterozygous deletion over three VDR exons (GENCODE 50).
#   - step 06 finds the planted hits of both; step 39 reports NAGA (SNV and
#     SNV) and VDR (SNV and deletion) as shared, 1 in 4 per pregnancy, and
#     UGT1A1 as the partner's only; HG002's planted SNAP25 hit (homozygous) is
#     his only;
#   - red controls: with the exon list gone the deletion does not count and
#     VDR is the partner's only; with the gene list emptied the step stops
#     and writes no table;
#   - both partners' reports (step 24) show the couple card;
#   - the pipeline with --pedigree and both as VCF rows (clinvar,html_report)
#     writes the same couple table and report as step 39 given the same
#     inputs (no Delly calls: the pipeline calls Delly from a BAM), and the
#     partner's HTML report has the card. somalier gets 20 sites, so the pair
#     is "not checked" rather than measured.
#   - a stub run with partners and without clinvar in --tools schedules no
#     COUPLE_CARRIER (red control: the same run with the gate taken out).
# Cost: about 5 minutes, most of it the Nextflow runs.
. "$(dirname "$0")/lib.sh"

G0="$GENOME_DIR"
C="${E2E_WORK}/genome-couple"
REF_NAME=GRCh38_no_alt_analysis_set
GENES="${REPO}/assets/carrier_genes.tsv"
in_e2e() { docker run --rm -v "${E2E_WORK}:/w" "$BCFTOOLS_IMAGE" "$@"; }
in_e2e rm -rf /w/genome-couple
mkdir -p "${C}/reference" "${C}/clinvar" "${C}/HG002/vcf" "${C}/PARTNER/vcf" "${C}/PARTNER/delly"
for f in "${G0}/reference/${REF_NAME}".fasta "${G0}/reference/${REF_NAME}".fasta.fai "${G0}/reference/${REF_NAME}".dict \
         "${G0}/HG002/vcf/HG002.vcf.gz" "${G0}/HG002/vcf/HG002.vcf.gz.tbi"; do
  d="${C}/reference"; case "$f" in */HG002/*) d="${C}/HG002/vcf" ;; esac
  ln -f "$f" "${d}/" 2>/dev/null || cp -f "$f" "${d}/"
done
export GENOME_DIR="$C"

span() { awk -F'\t' -v g="$1" '$1 == g {print $2 ":" $3 "-" $4}' "$GENES"; }
# het_snv GENE: the first PASS heterozygous SNV of HG002 in GENE (GQ 20 or more): POS REF ALT
het_snv() {
  bcf query -i 'TYPE="snp" && FILTER="PASS" && GT="het" && FMT/GQ>=20' -r "$(span "$1")" \
    -f '%POS\t%REF\t%ALT\n' HG002/vcf/HG002.vcf.gz 2>/dev/null | awk -F'\t' '$3 !~ /,/ {print; exit}'
}
refbase() { sam faidx "reference/${REF_NAME}.fasta" "${1}:${2}-${2}" 2>/dev/null | sed 1d | tr -d '\n' | tr '[:lower:]' '[:upper:]'; }
other() { case "$1" in A) echo C ;; *) echo A ;; esac; }

read -r N_POS N_REF N_ALT <<< "$(het_snv NAGA)"
read -r V_POS V_REF V_ALT <<< "$(het_snv VDR)"
check "HG002 has a PASS het SNV in NAGA (chr22:${N_POS:-none})" test -n "${N_POS:-}"
check "HG002 has a PASS het SNV in VDR (chr12:${V_POS:-none})" test -n "${V_POS:-}"
NAGA_END=$(span NAGA | cut -d- -f2)
P_POS=$(( ${N_POS:-0} + 101 )); [ "$P_POS" -le "${NAGA_END:-0}" ] || P_POS=$(( ${N_POS:-0} - 101 ))
P_REF=$(refbase chr22 "$P_POS"); P_ALT=$(other "$P_REF")
U_RANGE=$(span UGT1A1)
U_POS=$(( ($(cut -d: -f2 <<< "$U_RANGE" | cut -d- -f1) + $(cut -d- -f2 <<< "$U_RANGE")) / 2 ))
U_REF=$(refbase chr2 "$U_POS"); U_ALT=$(other "$U_REF")
check "the partner's NAGA and UGT1A1 positions have a reference base (${P_REF:-?}, ${U_REF:-?})" \
  sh -c "echo '${P_REF}${U_REF}' | grep -Eq '^[ACGT]{2}$'"

# --- ClinVar: the fixture's file plus the four planted records -------------------
REV="criteria_provided,_multiple_submitters,_no_conflicts"
plant() {  # CHROM POS ID REF ALT GENEINFO
  printf '%s\t%s\t%s\t%s\t%s\t.\t.\tALLELEID=%s;CLNDN=Synthetic_e2e_couple_record;CLNREVSTAT=%s;CLNSIG=Pathogenic;CLNVC=single_nucleotide_variant;GENEINFO=%s;ORIGIN=1\n' \
    "$1" "$2" "$3" "$4" "$5" "$3" "$REV" "$6"
}
{
  plant chr22 "$N_POS" 900000101 "$N_REF" "$N_ALT" NAGA:4668
  plant chr12 "$V_POS" 900000102 "$V_REF" "$V_ALT" VDR:7421
  plant chr22 "$P_POS" 900000103 "$P_REF" "$P_ALT" NAGA:4668
  plant chr2 "$U_POS" 900000104 "$U_REF" "$U_ALT" UGT1A1:54658
} > "${C}/clinvar/planted.txt"
cp "${G0}/clinvar/clinvar_pathogenic_chr.vcf.gz" "${C}/clinvar/fixture.vcf.gz"
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  { bcftools view -h clinvar/fixture.vcf.gz
    { bcftools view -H clinvar/fixture.vcf.gz; cat clinvar/planted.txt; } | sort -t\"\$(printf '\\t')\" -k1,1V -k2,2n
  } | bcftools view -Oz -o clinvar/clinvar_pathogenic_chr.vcf.gz -
  bcftools index -f -t clinvar/clinvar_pathogenic_chr.vcf.gz
  rm -f clinvar/fixture.vcf.gz"
check_eq "planted records in the couple's ClinVar file" \
  "$(bcf view -H clinvar/clinvar_pathogenic_chr.vcf.gz 2>/dev/null | cut -f3 | grep -cE '^90000010[1-4]$')" 4

# --- PARTNER: VCF and gVCF with HG002's header, and Delly's deletion -------------------
# header FILE: the header of one of HG002's files, with PARTNER as its sample
header() { gzip -dc "$1" | awk '/^#/ {print; next} {exit}' | sed '$ s/\tHG002$/\tPARTNER/'; }
{
  header "${G0}/HG002/vcf/HG002.g.vcf.gz"
  printf 'chr2\t%s\t.\t%s\t%s,<*>\t40\tPASS\t.\tGT:GQ:DP:AD:VAF:PL\t0/1:40:30:15,15,0:0.5,0:40,0,40,990,990,990\n' "$U_POS" "$U_REF" "$U_ALT"
  printf 'chr22\t%s\t.\t%s\t%s,<*>\t40\tPASS\t.\tGT:GQ:DP:AD:VAF:PL\t0/1:40:30:15,15,0:0.5,0:40,0,40,990,990,990\n' "$P_POS" "$P_REF" "$P_ALT"
} > "${C}/PARTNER/vcf/PARTNER.g.vcf"
{
  header "${G0}/HG002/vcf/HG002.vcf.gz"
  printf 'chr2\t%s\t.\t%s\t%s\t40\tPASS\t.\tGT:GQ:DP:AD:VAF:PL\t0/1:40:30:15,15:0.5:40,0,40\n' "$U_POS" "$U_REF" "$U_ALT"
  printf 'chr22\t%s\t.\t%s\t%s\t40\tPASS\t.\tGT:GQ:DP:AD:VAF:PL\t0/1:40:30:15,15:0.5:40,0,40\n' "$P_POS" "$P_REF" "$P_ALT"
} > "${C}/PARTNER/vcf/PARTNER.vcf"
check_eq "the partner's files name her as their sample" \
  "$(grep -h '^#CHROM' "${C}/PARTNER/vcf/PARTNER.vcf" "${C}/PARTNER/vcf/PARTNER.g.vcf" | awk -F'\t' '{print $NF}' | sort -u)" PARTNER
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  for f in PARTNER/vcf/PARTNER.vcf PARTNER/vcf/PARTNER.g.vcf; do
    bcftools sort -Oz -o \"\${f}.gz\" \"\$f\"; bcftools index -f -t \"\${f}.gz\"; rm -f \"\$f\"
  done"
check_eq "the partner's VCF holds her two SNVs" "$(vcf_count PARTNER/vcf/PARTNER.vcf.gz)" 2
# One heterozygous deletion over three VDR exons: chr12:47855001-47858000
printf '##fileformat=VCFv4.2\n##source=Synthetic_e2e_couple\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPARTNER\nchr12\t47855000\tDEL00000001\tN\t<DEL>\t60\tPASS\tPRECISE;SVTYPE=DEL;END=47858000\tGT:FT\t0/1:PASS\n' \
  | gzip -n > "${C}/PARTNER/delly/PARTNER_sv.vcf.gz"
# The merged GENCODE 50 exons of VDR
printf 'chr12\t%s\t%s\n' 47841536 47845005 47846334 47846451 47846656 47846808 47855629 47855801 47857128 47857249 \
  47857503 47857688 47865046 47865177 47878967 47879115 47882693 47882774 47899836 47899958 47904563 47905003 \
  47942862 47943048 > "${C}/reference/vdr_exons.bed"

printf 'sample\tsex\tfather\tmother\tpartner\nHG002\tmale\t.\t.\tPARTNER\nPARTNER\tfemale\t.\t.\tHG002\n' > "${C}/pedigree.tsv"

# --- step 06 for both, then step 39 -------------------------------------------------
for s in HG002 PARTNER; do
  run_step 06-clinvar-screen.sh "$s"
  check_step_exit 06-clinvar-screen.sh
done
check "HG002's planted NAGA and VDR hits" \
  sh -c "grep -q 900000101 '${C}/HG002/clinvar/HG002_clinvar_hits.vcf' && grep -q 900000102 '${C}/HG002/clinvar/HG002_clinvar_hits.vcf'"
check "the partner's planted NAGA and UGT1A1 hits" \
  sh -c "grep -q 900000103 '${C}/PARTNER/clinvar/PARTNER_clinvar_hits.vcf' && grep -q 900000104 '${C}/PARTNER/clinvar/PARTNER_clinvar_hits.vcf'"

# col TABLE GENE COLUMN: COLUMN of GENE's row in a couple table
col() { awk -F'\t' -v g="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $1 == g {print $c[k]; exit}' "$1" 2>/dev/null; }
genes_of() { awk -F'\t' -v r="$2" 'NR > 1 && $5 == r {print $1}' "$1" 2>/dev/null | sort | paste -sd, -; }

export CARRIER_EXONS="${C}/reference/vdr_exons.bed"
run_step 39-couple.sh couple "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
T="${C}/couple/couple_HG002_PARTNER.tsv"
cat "$T" 2>/dev/null | cut -f1,5,6,7,8,10
check_eq "shared genes: NAGA (SNV and SNV) and VDR (SNV and deletion)" "$(genes_of "$T" shared)" NAGA,VDR
check_eq "NAGA: 1 in 4 per pregnancy" "$(col "$T" NAGA offspring_risk)" "1 in 4 per pregnancy"
check "VDR: the partner's side is the deletion over three exons" \
  has '^deletion chr12:47855001-47858000 \(3.0 kb\) het, 3 exons of the gene \(Delly\)$' "$(col "$T" VDR findings_b)"
check "UGT1A1 is the partner's only" test "$(col "$T" UGT1A1 result)" = one_partner -a "$(col "$T" UGT1A1 findings_a)" = .
check "SNAP25 (the fixture's planted hit, homozygous) is HG002's only" \
  test "$(col "$T" SNAP25 result)" = one_partner -a "$(col "$T" SNAP25 findings_b)" = .
check "each partner's copy is the couple table" \
  sh -c "cmp -s '$T' '${C}/HG002/couple/HG002_couple_carrier.tsv' && cmp -s '$T' '${C}/PARTNER/couple/PARTNER_couple_carrier.tsv'"
MD="${C}/couple/couple_HG002_PARTNER.md"
check "the report names what was not assessed and points to a genetic counsellor" \
  sh -c "grep -q 'Deletions of HG002: no Delly calls' '$MD' && grep -q 'SMN1 copy number of PARTNER: no Parascopy' '$MD' && grep -q 'genetic counsellor' '$MD'"

# Red controls: no exon list, then an emptied gene list
cp "$T" "${CASE_TMP}/with_exons.tsv"
CARRIER_EXONS="${C}/reference/none.bed" run_step 39-couple.sh couple "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
check_eq "red: without the exon list VDR is not shared" "$(genes_of "$T" shared)" NAGA
head -n 15 "$GENES" > "${C}/reference/empty_genes.tsv"
CARRIER_GENES="${C}/reference/empty_genes.tsv" run_step 39-couple.sh couple "${C}/pedigree.tsv"
check "red: with the gene list emptied step 39 stops (exit ${STEP_RC})" test "$STEP_RC" -ne 0
check "and says why, and writes no couple table" \
  sh -c "grep -q 'holds no gene' '$STEP_LOG' && test ! -e '$T'"
check "and leaves no stale copy in either partner's folder" \
  sh -c "test ! -e '${C}/HG002/couple/HG002_couple_carrier.tsv' && test ! -e '${C}/PARTNER/couple/PARTNER_couple_carrier.tsv'"
run_step 39-couple.sh couple "${C}/pedigree.tsv"
check "the step back on the full inputs writes the first table again" cmp -s "$T" "${CASE_TMP}/with_exons.tsv"

# Both partners' reports
for s in HG002 PARTNER; do
  run_step 24-html-report.sh "$s"
  check_step_exit 24-html-report.sh
  check "${s}'s report has the couple card with NAGA shared" \
    has 'Shared: NAGA \(Schindler disease\): 1 in 4 per pregnancy' "$(cat "${C}/${s}/${s}_report.html" 2>/dev/null)"
done

# --- the pipeline -----------------------------------------------------------------
# The same inputs as the pipeline gets: no Delly calls (it calls Delly from a BAM).
mv "${C}/PARTNER/delly" "${CASE_TMP}/partner_delly"
unset CARRIER_EXONS
run_step 39-couple.sh couple_nf "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
mv "${CASE_TMP}/partner_delly" "${C}/PARTNER/delly"
B="${C}/couple_nf"
SITES="${CASE_TMP}/sites.vcf"
{
  printf '##fileformat=VCFv4.2\n'
  awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${C}/reference/${REF_NAME}.fasta.fai"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  bcf query -i 'TYPE="snp" && GT="het"' -r chr20 -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\tPASS\tAF=0.5\n' HG002/vcf/HG002.vcf.gz 2>/dev/null \
    | awk -F'\t' '$5 !~ /,/' | head -n 20
} > "$SITES"
SHEET="${CASE_TMP}/samplesheet.csv"
{
  echo 'sample,vcf,vcf_index,gvcf,gvcf_index,sex'
  echo "HG002,${G0}/HG002/vcf/HG002.vcf.gz,${G0}/HG002/vcf/HG002.vcf.gz.tbi,${G0}/HG002/vcf/HG002.g.vcf.gz,${G0}/HG002/vcf/HG002.g.vcf.gz.tbi,male"
  echo "PARTNER,${C}/PARTNER/vcf/PARTNER.vcf.gz,${C}/PARTNER/vcf/PARTNER.vcf.gz.tbi,${C}/PARTNER/vcf/PARTNER.g.vcf.gz,${C}/PARTNER/vcf/PARTNER.g.vcf.gz.tbi,female"
} > "$SHEET"
OUT="${CASE_TMP}/nf-out"
rm -rf "$OUT"
echo "+ nextflow run main.nf --pedigree (HG002 and PARTNER) --tools clinvar,html_report"
( cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-couple" \
    --input "$SHEET" --reference "${C}/reference/${REF_NAME}.fasta" --tools clinvar,html_report \
    --clinvar "${C}/clinvar/clinvar_pathogenic_chr.vcf.gz" --clinvar_index "${C}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
    --pedigree "${C}/pedigree.tsv" --somalier_sites "$SITES" --outdir "$OUT" --max_cpus 4 --max_memory 14.GB ) \
  > "${CASE_TMP}/nf.log" 2>&1
NF_RC=$?
grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "${CASE_TMP}/nf.log" | tail -n 40
cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
check_eq "pipeline with --pedigree: exit" "$NF_RC" 0
check "COUPLE_CARRIER ran" has 'COUPLE_CARRIER' "$(cat "${CASE_TMP}/nf.log")"
NT="${OUT}/family/couple_HG002_PARTNER.tsv"
check_eq "the pipeline's shared gene: NAGA" "$(genes_of "$NT" shared)" NAGA
for f in couple_HG002_PARTNER.tsv couple_HG002_PARTNER.md; do
  check "parity: ${f} of step 39 and of the pipeline are the same bytes" cmp -s "${B}/${f}" "${OUT}/family/${f}"
done
check "parity check sees a difference: a planted one is caught" \
  sh -c "sed '2s/\$/x/' '${B}/couple_HG002_PARTNER.tsv' > '${CASE_TMP}/planted.tsv'; ! cmp -s '${CASE_TMP}/planted.tsv' '$NT'"
check "the pipeline published each partner's copy" \
  sh -c "cmp -s '$NT' '${OUT}/HG002/couple/HG002_couple_carrier.tsv' && cmp -s '$NT' '${OUT}/PARTNER/couple/PARTNER_couple_carrier.tsv'"
check "the pipeline's report of PARTNER has the couple card" \
  has 'Shared: NAGA \(Schindler disease\): 1 in 4 per pregnancy' "$(cat "${OUT}/PARTNER/PARTNER_report.html" 2>/dev/null)"
echo "couple: shared $(genes_of "$T" shared); pipeline shared $(genes_of "$NT" shared)" >> "$E2E_NOTES"

# --- without clinvar in --tools -------------------------------------------------
# A stub run (no containers) of the stub samplesheet, whose pedigree marks two
# samples as partners, with Delly, ExpansionHunter and Parascopy but not
# clinvar: COUPLE_CARRIER must not be scheduled, because without ClinVar hits
# the screen stops with exit 1 and would fail the whole run. Red control: the
# same run on a copy of the repository with the clinvar gate taken out of
# main.nf schedules it.
tasks() { find "$1" -name .command.run -exec sed -n "s/^### name: '\(.*\)'\$/\1/p" {} + \
  | grep -c "$2" || true; }
stub_run() {  # stub_run DIR: copy the tracked files to DIR and run the stub there
  rm -rf "$1"; mkdir -p "$1"
  git -C "$REPO" ls-files -z | (cd "$REPO" && tar --null -T - -cf -) | tar -xf - -C "$1"
  [ "${2:-}" = ungated ] && sed -i "s/couple_ids && tools_list.contains('clinvar') ?/couple_ids ?/" "${1}/main.nf"
  ( cd "$1" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${1}/out" -work-dir "${1}/work" \
      --pedigree assets/stub/pedigree.tsv --tools delly,expansion_hunter,parascopy,html_report ) > "${1}.log" 2>&1
}
stub_run "${CASE_TMP}/stub-ungated" ungated
check "red control: the copy's main.nf has the gate taken out" \
  sh -c "grep -q 'def ch_couple_outputs = couple_ids ?' '${CASE_TMP}/stub-ungated/main.nf'"
check_ge "red: with the gate taken out, the stub run without clinvar schedules COUPLE_CARRIER" \
  "$(tasks "${CASE_TMP}/stub-ungated/work" COUPLE_CARRIER)" 1
stub_run "${CASE_TMP}/stub-gated"
STUB_RC=$?
tail -n 15 "${CASE_TMP}/stub-gated.log"
check_eq "stub run with partners and without clinvar: exit" "$STUB_RC" 0
check "it warns that the couple screen does not run" grep -q 'clinvar is not in --tools' "${CASE_TMP}/stub-gated.log"
check_ge "it ran the partners' Delly tasks" "$(tasks "${CASE_TMP}/stub-gated/work" DELLY)" 2
check_eq "and schedules no COUPLE_CARRIER task" "$(tasks "${CASE_TMP}/stub-gated/work" COUPLE_CARRIER)" 0
rm -rf "${CASE_TMP}/stub-ungated" "${CASE_TMP}/stub-gated"
in_e2e rm -rf /w/genome-couple /w/nf-work-couple
finish
