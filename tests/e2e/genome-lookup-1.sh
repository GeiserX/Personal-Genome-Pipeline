#!/usr/bin/env bash
# scripts/pgp-lookup.sh on the HG002 gVCF case 21 wrote and the annotated VCF
# case 40 wrote, with the real bcftools and python images:
#   - an rsID the annotated VCF carries resolves to its position, and the
#     answer has the genotype bcftools reads from the gVCF there, the gene
#     and the rsID from VEP's CSQ;
#   - a position inside a reference block of the gVCF reads as homozygous
#     reference, naming the block;
#   - a position no gVCF record covers reads as a no call;
#   - a gene (a one-line GENCODE file in a REF_DIR outside GENOME_DIR, so the
#     /pgp-ref mount is used) lists the annotated variant inside it and the
#     share of its bases with a call;
#   - an rsID found nowhere exits 3 with the reason.
# Writes nothing in the sample folder (the slices are removed on exit).
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
GVCF="${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
ANN="${G}/${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz"
check "case 21 wrote a gVCF" vcf_ok "$GVCF"
check "case 40 wrote the annotated VCF" nonempty "${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz.tbi"

lookup() {
  echo "+ scripts/pgp-lookup.sh ${SAMPLE} $*"
  "${REPO}/scripts/pgp-lookup.sh" "$SAMPLE" "$@" > "$STEP_LOG" 2>&1
  STEP_RC=$?
  cat "$STEP_LOG"
  echo "+ exit ${STEP_RC}"
  OUT=$(cat "$STEP_LOG")
}

# --- an rsID of the annotated VCF ----------------------------------------------------
# The first record whose CSQ has a gene symbol and an rsID in Existing_variation.
PICK=$(gzip -dc "$ANN" | python3 -c '
import re, sys
fields = []
for line in sys.stdin:
    if line.startswith("##INFO=<ID=CSQ,"):
        fields = re.search(r"Format: ([^\"]+)\"", line).group(1).split("|")
        continue
    if line.startswith("#"):
        continue
    c = line.rstrip("\n").split("\t")
    m = re.search(r"(?:^|;)CSQ=([^;]+)", c[7])
    if not m or not fields:
        continue
    for e in m.group(1).split(","):
        d = dict(zip(fields, e.split("|")))
        rs = [x for x in d.get("Existing_variation", "").split("&") if re.fullmatch(r"rs\d+", x)]
        if d.get("SYMBOL") and rs and len(c[3]) == 1 and len(c[4]) == 1:
            print(c[0], c[1], rs[0], d["SYMBOL"], sep="\t")
            sys.exit(0)
')
IFS=$'\t' read -r V_CHR V_POS V_RS V_SYM <<< "$PICK"
check "the annotated VCF has an SNV with a gene and an rsID (${PICK:-none})" test -n "${V_RS:-}"

lookup "$V_RS"
check_eq "rsID lookup exits 0" "$STEP_RC" 0
check "the rsID resolves to ${V_CHR}:${V_POS}" has "^== ${V_RS}: ${V_CHR}:${V_POS} " "$OUT"
check "the answer names the rsID among the variant's ids" has "known as: .*${V_RS}" "$OUT"
check "the answer has a gene line" has "^    gene: " "$OUT"
# What bcftools reads from the gVCF there: a variant record at the position,
# else the block covering it, else nothing.
GT_AT=$(bcf query -r "${V_CHR}:${V_POS}" -f '%POS\t%ALT\t[%TGT]\n' "$GVCF" 2>/dev/null \
  | awk -F'\t' -v p="$V_POS" '$1 == p && $2 != "<*>" { print $3; exit }')
if [ -n "$GT_AT" ]; then
  check "the genotype is bcftools' ${GT_AT}" has "^  genotype: ${GT_AT//|/\\|} \(" "$OUT"
else
  check "no variant record in the gVCF at ${V_POS}: a block or a no call" has '^  genotype: (ref/ref|[ACGT]/[ACGT] \(homozygous reference|no call)' "$OUT"
fi

# --- a position inside a reference block -------------------------------------------
BLOCK=$(bcf query -i 'ALT="<*>" && GT="0/0"' -f '%CHROM\t%POS\t%INFO/END\n' "$GVCF" 2>/dev/null \
  | awk -F'\t' '$3 != "." && $3 - $2 >= 10 { print; exit }')
IFS=$'\t' read -r B_CHR B_POS B_END <<< "$BLOCK"
check "the gVCF has a reference block of 10 bases or more (${BLOCK:-none})" test -n "${B_END:-}"
lookup "${B_CHR}:$((B_POS + 5))"
check_eq "position lookup exits 0" "$STEP_RC" 0
check "inside the block ${B_CHR}:${B_POS}-${B_END}: homozygous reference" \
  has "^  genotype: ref/ref \(homozygous reference, inside the reference block ${B_CHR}:${B_POS}-${B_END}\)" "$OUT"

# --- a position no record covers ------------------------------------------------------
# chr1:1 is in the telomere N run, outside every fixture slice.
check "control: chr1:1 lies outside the fixture slices" test -z "$(awk -F'\t' '$1 == "chr1" && $2 < 1' "${FIXTURE_DIR}/regions.bed")"
check_eq "control: bcftools finds no gVCF record over chr1:1" "$(vcf_count -r chr1:1 "$GVCF")" 0
lookup chr1:1
check "chr1:1 is a no call, not the reference" has '^  genotype: no call \(no gVCF record covers this position' "$OUT"

# --- a gene, from a REF_DIR outside GENOME_DIR -----------------------------------------
REL=$(sed -n 's/^GENCODE_RELEASE=//p' "${REPO}/scripts/lib/common.sh")
RD="${CASE_TMP}/ref"
mkdir -p "$RD"
printf '%s\tHAVANA\tgene\t%s\t%s\t.\t+\t.\tgene_id "ENSG_E2E"; gene_type "protein_coding"; gene_name "%s"; level 2;\n' \
  "$V_CHR" "$((V_POS - 500))" "$((V_POS + 500))" "$V_SYM" > "${RD}/gencode.v${REL}.basic.genes.gtf"
REF_DIR="$RD" lookup "$V_SYM"
check_eq "gene lookup exits 0" "$STEP_RC" 0
check "the gene span comes from the REF_DIR file" has "^== ${V_SYM}: ${V_CHR}:$((V_POS - 500))-$((V_POS + 500)) " "$OUT"
check "the variant ${V_CHR}:${V_POS} is listed" has "^  - ${V_CHR}:${V_POS} " "$OUT"
check "the share of bases with a call" has 'of 1,001 bases have a genotype call' "$OUT"

# --- an rsID found nowhere --------------------------------------------------------------
lookup rs999999999999
check_eq "an unknown rsID exits 3" "$STEP_RC" 3
check "and says where it looked" has 'ERROR: rs999999999999 was not found in' "$OUT"

check "no slice is left in the sample folder" test -z "$(find "${G}/${SAMPLE}" -maxdepth 1 -name '.lookup.*')"

finish
