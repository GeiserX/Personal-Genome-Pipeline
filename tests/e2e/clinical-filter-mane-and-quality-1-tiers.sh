#!/usr/bin/env bash
# Steps 23 and 31 on a synthetic VEP sample (MANEQ) whose records test the
# transcript the impact is read from, the genotype check and the frequency tier:
#   1000 HIGH on a non-MANE transcript, intronic on the MANE Select one, rare,
#        ClinVar pathogenic          -> not in step 23's HIGH tier; its summary
#                                       row says MODIFIER on MANE, HIGH as worst
#   2000 HIGH, 0/1 at GQ 12, AD 3,2  -> HIGH tier, CALL_CHECK weak (GQ12)
#   2500 HIGH, 1/1 at GQ 11, AD 12,5 -> HIGH tier, weak (GQ11, AB0.29)
#   2600 HIGH, 1/1 at GQ 50, AD 0,30 -> HIGH tier, ok
#   3000 HIGH 3-base deletion, no MAX_AF or gnomAD AF
#                                    -> freq_unknown tier, not HIGH, in both steps
#   4000 HIGH, MAX_AF 0.3            -> in no tier
#   5000 no MANE transcript: MODERATE on the canonical one, HIGH on another
#                                    -> step 23's MODERATE tier, TRANSCRIPT CANONICAL
#   6000 HIGH and 6100 MODERATE in one gene, both het and rare
#                                    -> step 31's compound-het TSV, with the
#                                       gnomAD AF and the MANE consequence
# Then step 23 again with CLINICAL_MIN_GQ=0: 2000 is no longer weak, 2500 still
# is (its alt reads), so the floor is the one the variable sets.
. "$(dirname "$0")/lib.sh"

S=MANEQ
FIELDS='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|CANONICAL|MANE_SELECT|Existing_variation|MAX_AF|gnomADe_AF|gnomADg_AF|CLIN_SIG'
mkdir -p "${GENOME_DIR}/${S}/vep"
{
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n'
  printf '##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: %s">\n' "$FIELDS"
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
  printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
  printf '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$S"
  # pos ref alt sample; then one line per CSQ entry:
  # consequence impact symbol canonical mane MAX_AF gnomADe gnomADg clin_sig ('-' = empty)
  while read -r pos ref alt sample n; do
    csq=""
    for _ in $(seq "$n"); do
      read -r cons impact sym can mane max e g clin
      for v in can mane max e g clin; do [ "${!v}" = - ] && printf -v "$v" '%s' ''; done
      csq+="${csq:+,}${alt}|${cons}|${impact}|${sym}|ENSG${sym}|Transcript|ENST${pos}${#csq}|protein_coding|${can}|${mane}|rs${pos}|${max}|${e}|${g}|${clin}"
    done
    printf 'chr20\t%s\t.\t%s\t%s\t50\tPASS\tCSQ=%s\tGT:GQ:DP:AD\t%s\n' "$pos" "$ref" "$alt" "$csq" "$sample"
  done <<'RECORDS'
1000 A G 0/1:50:30:15,15 2
stop_gained HIGH GENEA - - 0.0001 0.0001 0.0001 pathogenic
intron_variant MODIFIER GENEA YES NM_0001.1 0.0001 0.0001 0.0001 pathogenic
2000 C T 0/1:12:.:3,2 1
frameshift_variant HIGH GENEB YES NM_0002.1 0.0001 0.0001 0.0001 -
2500 G A 1/1:11:17:12,5 1
stop_gained HIGH GENEC YES NM_0003.1 0.0001 0.0001 0.0001 -
2600 G A 1/1:50:30:0,30 1
stop_gained HIGH GENED YES NM_0004.1 0.0001 0.0001 0.0001 -
3000 CTTT C 0/1:50:30:15,15 1
frameshift_variant HIGH GENEE YES NM_0005.1 - - - -
4000 T C 0/1:50:30:15,15 1
stop_gained HIGH GENEF YES NM_0006.1 0.3 0.3 0.3 -
5000 A C 0/1:50:30:15,15 2
missense_variant MODERATE GENEG YES - 0.0001 0.0001 0.0001 -
stop_gained HIGH GENEG - - 0.0001 0.0001 0.0001 -
6000 C T 0/1:50:30:15,15 1
stop_gained HIGH GENEH YES NM_0008.1 0.0002 0.0001 0.0002 -
6100 G A 0/1:50:30:15,15 1
missense_variant MODERATE GENEH YES NM_0008.1 0.0002 0.0001 0.0002 -
RECORDS
} > "${GENOME_DIR}/${S}/vep/${S}_vep.vcf"
grep -v '^##' "${GENOME_DIR}/${S}/vep/${S}_vep.vcf"

positions() { bcf query -f '%POS\n' "$1" 2>/dev/null | paste -sd, -; }
# cell FILE POS COLUMN: the value of COLUMN in the row of POS.
cell() {
  awk -F'\t' -v p="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    $c["POS"] == p && (k in c) {print $c[k]; exit}' "$1" 2>/dev/null
}

# --- step 23 -------------------------------------------------------------------
run_step 23-clinical-filter.sh "$S"
check_step_exit 23-clinical-filter.sh
LOG23=$(cat "$STEP_LOG")
C="${S}/clinical/${S}"
SUM="${GENOME_DIR}/${C}_clinical_summary.tsv"
cat "$SUM" 2>/dev/null
check_eq "step 23: rare HIGH tier on the MANE Select transcript" "$(positions "${C}_high_impact.vcf.gz")" 2000,2500,2600,6000
check_eq "step 23: rare MODERATE tier (5000 on its canonical transcript)" "$(positions "${C}_rare_moderate.vcf.gz")" 5000,6100
check_eq "step 23: the deletion with no frequency is in its own tier" "$(positions "${C}_freq_unknown.vcf.gz")" 3000
check "step 23: the log counts it apart" has 'Frequency unknown \(not rare\): 1$' "$LOG23"
check_eq "step 23: 1000 (HIGH off MANE) keeps its MANE consequence in the summary" \
  "$(cell "$SUM" 1000 IMPACT) $(cell "$SUM" 1000 Consequence) $(cell "$SUM" 1000 TRANSCRIPT)" "MODIFIER intron_variant MANE"
check_eq "step 23: and shows the worst one beside it" \
  "$(cell "$SUM" 1000 WORST_IMPACT) $(cell "$SUM" 1000 WORST_Consequence)" "HIGH stop_gained"
check_eq "step 23: 5000 is ranked on its canonical transcript" "$(cell "$SUM" 5000 TRANSCRIPT) $(cell "$SUM" 5000 IMPACT)" "CANONICAL MODERATE"
check_eq "step 23: 3000's frequency tier" "$(cell "$SUM" 3000 FREQ_TIER)" unknown
check_eq "step 23: 2000 carries its GQ and AD" "$(cell "$SUM" 2000 GQ) $(cell "$SUM" 2000 AD)" "12 3,2"
check_eq "step 23: the GQ 12 heterozygote is a weak call" "$(cell "$SUM" 2000 CALL_CHECK) $(cell "$SUM" 2000 CALL_WHY)" "weak GQ12"
check_eq "step 23: the GQ 11 1/1 on AD 12,5 is weak by GQ and alt reads" \
  "$(cell "$SUM" 2500 CALL_CHECK) $(cell "$SUM" 2500 CALL_WHY)" "weak GQ11,AB0.29"
check_eq "step 23: the GQ 50 1/1 on AD 0,30 is not" "$(cell "$SUM" 2600 CALL_CHECK)" ok
check "step 23: the log counts the weak calls" has 'Weak genotype calls \(CALL_CHECK=weak, check the reads\): 2$' "$LOG23"
check "step 23: the clinical VCF names the floor it used" \
  has '^##PGP_CALL_FLOOR=GQ=30,DP=10,HET_AB=0.2-0.8,HOM_AB=0.8$' "$(bcf view -h "${C}_clinical.vcf.gz" 2>/dev/null)"

CLINICAL_MIN_GQ=0 run_step 23-clinical-filter.sh "$S"
check_step_exit 23-clinical-filter.sh
check_eq "step 23 with CLINICAL_MIN_GQ=0: 2000 is no longer weak" "$(cell "$SUM" 2000 CALL_CHECK)" ok
check_eq "step 23 with CLINICAL_MIN_GQ=0: 2500 is still weak on its alt reads" "$(cell "$SUM" 2500 CALL_WHY)" AB0.29

# --- step 31 -------------------------------------------------------------------
run_step 31-slivar.sh "$S"
check_step_exit 31-slivar.sh
LOG31=$(cat "$STEP_LOG")
L="${S}/slivar/${S}"
check_eq "step 31: the deletion with no frequency is not rare_high" \
  "$(positions "${L}_freq_unknown.vcf.gz")" 3000
check "step 31: the log counts it apart" has 'freq_unknown \(HIGH/MODERATE, not rare\): 1$' "$LOG31"
check_eq "step 31: prioritized keeps it for the compound-het search" \
  "$(positions "${L}_prioritized.vcf.gz")" 1000,2000,2500,2600,3000,5000,6000,6100
SUM31="${GENOME_DIR}/${L}_slivar_summary.tsv"
CH="${GENOME_DIR}/${L}_compound_hets.tsv"
cat "$SUM31" "$CH" 2>/dev/null
check_eq "step 31 summary: 3000's frequency tier" "$(cell "$SUM31" 3000 FREQ_TIER)" unknown
check_eq "step 31 summary: 1000's MANE consequence beside its worst one" \
  "$(cell "$SUM31" 1000 Consequence) $(cell "$SUM31" 1000 MANE_Consequence)" "stop_gained intron_variant"
check_eq "step 31 summary: no MANE consequence for 5000 (no MANE transcript)" \
  "$(cell "$SUM31" 5000 TRANSCRIPT) $(cell "$SUM31" 5000 MANE_Consequence)" "CANONICAL ."
check "step 31: the compound-het TSV has the AF and MANE columns" \
  has $'\tMAX_AF\tgnomADe_AF\tgnomADg_AF\tFREQ_TIER\tTRANSCRIPT\tMANE_Consequence$' "$(head -1 "$CH" 2>/dev/null)"
check_eq "step 31: the GENEH pair's first row carries its gnomAD AF and MANE consequence" \
  "$(cell "$CH" 6000 gnomADe_AF) $(cell "$CH" 6000 gnomADg_AF) $(cell "$CH" 6000 FREQ_TIER) $(cell "$CH" 6000 MANE_Consequence)" \
  "0.0001 0.0002 rare stop_gained"
check_eq "step 31: and its partner" "$(cell "$CH" 6100 GENE) $(cell "$CH" 6100 MANE_Consequence)" "GENEH missense_variant"

finish
