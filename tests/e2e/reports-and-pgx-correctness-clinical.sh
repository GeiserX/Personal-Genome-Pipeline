#!/usr/bin/env bash
# Steps 23 and 31: real gene symbols, one rarity rule on every tier, both
# SpliceAI genes, and one gnomAD constraint loader.
#   1. HG002 (case 41's output): no HIGH or MODERATE row of the clinical summary
#      has '.' as its gene (the old step read a SYMBOL= tag nothing writes).
#   2. With a gnomAD-style constraint table for the fixture's genes (a
#      non-canonical decoy row before each canonical one), steps 23 and 31 print
#      the same LOEUF for a gene and a numeric mis_z.
#   3. Two synthetic CSQ samples, one with MAX_AF and one with only gnomADe_AF
#      and gnomADg_AF: a missense common in genomes but absent from exomes, and
#      a common stop-gain, are not in the rare tiers of either step; a rare
#      missense and a stop-gain with no frequency are; a SpliceAI value whose
#      second gene scores 0.5 reaches the splice tier.
. "$(dirname "$0")/lib.sh"

SUM23="${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv"
head -5 "$SUM23" 2>/dev/null
check "the clinical summary has a GENE column" has '(^|	)GENE(	|$)' "$(head -1 "$SUM23" 2>/dev/null)"
NOGENE=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
  ($c["IMPACT"] == "HIGH" || $c["IMPACT"] == "MODERATE") && $c["GENE"] == "." {n++} END {print n + 0}' "$SUM23" 2>/dev/null)
check_eq "HIGH/MODERATE rows whose gene is '.'" "${NOGENE:-missing}" 0
check_ge "rows with a gene symbol" "$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["GENE"] != "." {n++} END {print n + 0}' "$SUM23" 2>/dev/null)" 1

# --- 2. constraint join in steps 23 and 31 -------------------------------------
CONSTRAINT="${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv"
trap 'rm -f "$CONSTRAINT"' EXIT
SUM31="${GENOME_DIR}/${SAMPLE}/slivar/${SAMPLE}_slivar_summary.tsv"
{ awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["GENE"] != "." {print $c["GENE"]}' "$SUM23"
  awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["SYMBOL"] != "." && $c["SYMBOL"] != "" {print $c["SYMBOL"]}' "$SUM31"
} 2>/dev/null | LC_ALL=C sort -u > "${CASE_TMP}/genes.txt"
echo "genes for the constraint table: $(paste -sd' ' "${CASE_TMP}/genes.txt")"
awk 'BEGIN {OFS = "\t"; print "gene", "gene_id", "transcript", "canonical", "mane_select", "lof.oe_ci.upper", "lof.pLI", "mis.z_score"}
     {n++; print $1, "ENSGX" n, "ENST9" n, "false", "false", "0.99", "0.00", "0.10"
            print $1, "ENSGX" n, "ENST1" n, "true", "true", sprintf("0.%02d", n), "0.95", sprintf("2.%02d", n)}' \
  "${CASE_TMP}/genes.txt" > "$CONSTRAINT"
run_step 23-clinical-filter.sh "$SAMPLE"
check_step_exit 23-clinical-filter.sh
run_step 31-slivar.sh "$SAMPLE"
check_step_exit 31-slivar.sh
col() {  # col FILE GENE_COLUMN VALUE_COLUMN: "gene<TAB>value" per row
  awk -F'\t' -v g="$2" -v v="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} (g in c) && (v in c) {print $c[g] "\t" $c[v]}' "$1" 2>/dev/null | LC_ALL=C sort -u
}
col "$SUM23" GENE LOEUF > "${CASE_TMP}/loeuf23.tsv"
col "$SUM31" SYMBOL LOEUF > "${CASE_TMP}/loeuf31.tsv"
SHARED=$(LC_ALL=C join -t $'\t' "${CASE_TMP}/loeuf23.tsv" "${CASE_TMP}/loeuf31.tsv")
echo "gene, LOEUF in step 23, LOEUF in step 31:"; printf '%s\n' "${SHARED:-none}"
check_ge "genes in both summaries" "$(grep -c . <<< "$SHARED" || true)" 1
check_eq "genes whose LOEUF differs between steps 23 and 31" "$(awk -F'\t' '$2 != $3' <<< "$SHARED" | grep -c . || true)" 0
check_eq "shared genes given the decoy (non-canonical) LOEUF 0.99" "$(awk -F'\t' '$2 == "0.99"' <<< "$SHARED" | grep -c . || true)" 0
check_ge "step 23 rows with a numeric mis_z" \
  "$(col "$SUM23" GENE mis_z | awk -F'\t' '$2 ~ /^[0-9.]+$/' | grep -c . || true)" 1
check_ge "step 31 rows with a numeric mis_z" \
  "$(col "$SUM31" SYMBOL mis_z | awk -F'\t' '$2 ~ /^[0-9.]+$/' | grep -c . || true)" 1
rm -f "$CONSTRAINT"

# --- 3. synthetic frequencies --------------------------------------------------
# write_sample NAME WITH_MAX_AF: a VEP-style VCF at NAME/vep/NAME_vep.vcf. With
# WITH_MAX_AF=no the CSQ has gnomADe_AF and gnomADg_AF only.
#   1000 missense, absent from exomes, 20% in genomes      -> not rare
#   2000 missense, rare everywhere                         -> rare MODERATE
#   3000 stop-gain, 30% everywhere                         -> not rare
#   4000 stop-gain, no frequency at all                    -> rare HIGH
#   5000 intron, SpliceAI 0.50 for its second gene, rare   -> splice tier
write_sample() {
  local name=$1 with_max=$2 fields='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|CANONICAL|gnomADe_AF|gnomADg_AF'
  [ "$with_max" = yes ] && fields+='|MAX_AF'
  fields+='|CLIN_SIG'
  mkdir -p "${GENOME_DIR}/${name}/vep"
  {
    printf '##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n'
    printf '##FILTER=<ID=PASS,Description="All filters passed">\n'
    printf '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: %s">\n' "$fields"
    printf '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI scores">\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$name"
    # pos ref alt consequence impact symbol gnomADe gnomADg MAX_AF extra-INFO
    while read -r pos ref alt cons impact sym e g max extra; do
      [ "$e" = - ] && e=""; [ "$g" = - ] && g=""; [ "$max" = - ] && max=""; [ "$extra" = - ] && extra=""
      csq="${alt}|${cons}|${impact}|${sym}|ENSG${pos}|Transcript|ENST${pos}|protein_coding|YES|${e}|${g}"
      [ "$with_max" = yes ] && csq+="|${max}"
      csq+="|"
      printf 'chr20\t%s\t.\t%s\t%s\t50\tPASS\tCSQ=%s%s\tGT\t0/1\n' "$pos" "$ref" "$alt" "$csq" "$extra"
    done <<'RECORDS'
1000 A G missense_variant MODERATE GENEX - 0.2 0.2 -
2000 C T missense_variant MODERATE GENEY 0.0001 0.0002 0.0002 -
3000 G A stop_gained HIGH GENEZ 0.3 0.3 0.3 -
4000 T C stop_gained HIGH GENEW - - - -
5000 A C intron_variant MODIFIER GENEA - - - ;SpliceAI=C|GENEA|0.00|0.00|0.01|0.00|1|2|3|4,C|GENEB|0.50|0.00|0.00|0.00|1|2|3|4
RECORDS
  } > "${GENOME_DIR}/${name}/vep/${name}_vep.vcf"
}
write_sample FREQMAX yes
write_sample FREQEG no
grep -v '^##' "${GENOME_DIR}/FREQEG/vep/FREQEG_vep.vcf"

positions() { bcf query -f '%POS\n' "$1" 2>/dev/null | paste -sd, -; }
for S in FREQMAX FREQEG; do
  run_step 23-clinical-filter.sh "$S"
  check_step_exit 23-clinical-filter.sh
  check "${S} step 23: says which frequency it filtered on" has 'Population frequency: (MAX_AF|gnomADe_AF and gnomADg_AF)' "$(cat "$STEP_LOG")"
  MOD=$(positions "${S}/clinical/${S}_rare_moderate.vcf.gz")
  HIGH=$(positions "${S}/clinical/${S}_high_impact.vcf.gz")
  SPLICE=$(positions "${S}/clinical/${S}_spliceai_high.vcf.gz")
  echo "${S} step 23: MODERATE ${MOD:-none}; HIGH ${HIGH:-none}; SpliceAI ${SPLICE:-none}"
  check_eq "${S} step 23: rare MODERATE is the rare missense only" "${MOD:-none}" 2000
  check_eq "${S} step 23: rare HIGH is the stop-gain with no frequency only" "${HIGH:-none}" 4000
  check_eq "${S} step 23: the SpliceAI tier tests the value's second gene" "${SPLICE:-none}" 5000
  check_eq "${S} step 23: summary rows with gene '.'" \
    "$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["GENE"] == "." {n++} END {print n + 0}' \
        "${GENOME_DIR}/${S}/clinical/${S}_clinical_summary.tsv" 2>/dev/null)" 0

  run_step 31-slivar.sh "$S"
  check_step_exit 31-slivar.sh
  PRI=$(positions "${S}/slivar/${S}_prioritized.vcf.gz")
  echo "${S} step 31: prioritized ${PRI:-none}"
  check_eq "${S} step 31: prioritized are the rare missense and the rare stop-gain" "${PRI:-none}" 2000,4000
done

finish
