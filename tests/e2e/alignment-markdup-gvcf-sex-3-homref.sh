#!/usr/bin/env bash
# Steps 25, 07 and 14 read hom-ref genotypes from the gVCF step 03 writes
# (case 21 called it; case 31 ran PharmCAT on it):
#   PRS: a synthetic five-site score whose effect alleles are the reference
#     bases at sites where HG002 matches the reference (GIAB v4.2.1 truth, no
#     variant within 200 bp) sums to 2 x (1+2+4+8+16) = 62 from the gVCF, and
#     to less from the variant-only VCF;
#   PharmCAT: PGx positions inside the fixture's CYP2C19/CYP2C9 slice are
#     reference calls, not missing;
#   imputation prep: panel sites become 0/0 calls.
. "$(dirname "$0")/lib.sh"

D="${GENOME_DIR}/${SAMPLE}"
GVCF="${D}/vcf/${SAMPLE}.g.vcf.gz"
check "case 21 wrote a gVCF" vcf_ok "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
check "gVCF index exists" nonempty "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz.tbi"

# --- five hom-ref sites on the chr20 slice, from the truth set -------------------
SITES="${CASE_TMP}/homref_sites.tsv"
python3 - "${FIXTURE_DIR}/HG002_truth_chr20.bed" "${FIXTURE_DIR}/HG002_truth_chr20.vcf.gz" \
  "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" > "$SITES" <<'PY'
import gzip, sys
bed, truth, fasta = sys.argv[1:4]
iv = [tuple(int(x) for x in l.split()[1:3]) for l in open(bed) if l.startswith("chr20\t")]
var = [int(l.split("\t")[1]) for l in gzip.open(truth, "rt") if l.startswith("chr20\t")]
fai = {l.split()[0]: [int(x) for x in l.split()[1:5]] for l in open(fasta + ".fai")}
length, offset, bases, width = fai["chr20"]
def base(pos):
    i = pos - 1
    with open(fasta, "rb") as f:
        f.seek(offset + (i // bases) * width + i % bases)
        return f.read(1).decode().upper()
out, p = [], 10_050_000
while len(out) < 5 and p < 10_450_000:
    inside = any(s + 200 < p <= e - 200 for s, e in iv)
    clear = all(abs(v - p) > 200 for v in var)
    b = base(p)
    if inside and clear and b in "ACGT":
        out.append((p, b))
    p += 7_919
for p, b in out:
    print(f"chr20\t{p}\t{b}")
PY
cat "$SITES"
check_eq "hom-ref sites chosen on the chr20 slice" "$(wc -l < "$SITES" | tr -d ' ')" 5

# --- PRS ---------------------------------------------------------------------
PRS_IDS="PGS000018 PGS000014 PGS000004 PGS000662 PGS000016 PGS000334 PGS000027 PGS000017 PGS000055"
mkdir -p "${GENOME_DIR}/prs_scores"
for id in $PRS_IDS; do
  { printf '#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\n'
    awk -F'\t' -v OFS='\t' '{print "20", $2, $3, 2 ^ (NR - 1)}' "$SITES"; } | gzip -c > "${GENOME_DIR}/prs_scores/${id}.txt.gz"
done
SUMMARY="${D}/prs/${SAMPLE}_prs_summary.tsv"
row() { awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i} $2 == "PGS000018" {print $c["Score_SUM"], $c["Variants_Matched"], $c["Input"]}' "$SUMMARY" 2>/dev/null; }

run_step 25-prs.sh "$SAMPLE"
check_step_exit 25-prs.sh
read -r G_SUM G_MATCHED G_INPUT <<< "$(row)"
G_LOG=$(cat "$STEP_LOG")
check_eq "PRS input with the gVCF present" "${G_INPUT:-}" gvcf
check_eq "PRS sum from the gVCF (2 x 31)" "${G_SUM:-}" 62
check_eq "PRS variants matched from the gVCF" "${G_MATCHED:-}" 5
check "no 'hom-ref sites are absent' banner on the gVCF path" lacks 'hom-ref sites are absent' "$G_LOG"

mkdir -p "${CASE_TMP}/hidden"
mv "$GVCF" "${GVCF}.tbi" "${CASE_TMP}/hidden/"
run_step 25-prs.sh "$SAMPLE"
check_step_exit "25-prs.sh (variant-only VCF)"
read -r V_SUM V_MATCHED V_INPUT <<< "$(row)"
V_LOG=$(cat "$STEP_LOG")
mv "${CASE_TMP}/hidden/${SAMPLE}.g.vcf.gz" "${CASE_TMP}/hidden/${SAMPLE}.g.vcf.gz.tbi" "${D}/vcf/"
check_eq "PRS input without a gVCF" "${V_INPUT:-}" vcf
check "the variant-only sum (${V_SUM:-none}) is lower than the gVCF sum (${G_SUM:-none})" \
  awk -v v="${V_SUM:-NA}" -v g="${G_SUM:-NA}" 'BEGIN { if (v == "NA") v = 0; exit !(g != "NA" && v + 0 < g + 0) }'
check "the 'hom-ref sites are absent' banner on the variant-only path" has 'hom-ref sites are absent' "$V_LOG"
printf '#### PRS synthetic control (five sites where HG002 is hom-ref, effect allele = reference)\n\nfrom the gVCF: sum %s, %s of 5 matched; from the variant-only VCF: sum %s, %s of 5 matched\n\n' \
  "${G_SUM:-none}" "${G_MATCHED:-none}" "${V_SUM:-none}" "${V_MATCHED:-none}" >> "$E2E_NOTES"
rm -rf "${GENOME_DIR:?}/prs_scores" "${D:?}/prs"

# --- PharmCAT (case 31) ---------------------------------------------------------
P31_LOG="${E2E_WORK}/logs/31-pharmcat-07.log"
check "case 31 gave PharmCAT the gVCF" grep -q "^Input: .*${SAMPLE}\.g\.vcf\.gz" "$P31_LOG"
# chr10:94700000-95000000 is the fixture's CYP2C19 + CYP2C9 slice.
POSITIONS=$(docker run --rm "$PHARMCAT_IMAGE" sh -c 'gzip -dc /pharmcat/pharmcat_positions.vcf.bgz' \
  | awk -F'\t' '$1 == "chr10" && $2 >= 94700000 && $2 <= 95000000 {print $2}' | sort -u | wc -l | tr -d ' ')
MISSING_VCF="${D}/vcf/${SAMPLE}.missing_pgx_var.vcf"
check "PharmCAT's missing-positions report exists" test -f "$MISSING_VCF"
MISSING=$(awk -F'\t' '$1 == "chr10" && $2 >= 94700000 && $2 <= 95000000 {print $2}' "$MISSING_VCF" 2>/dev/null | sort -u | wc -l | tr -d ' ')
echo "PharmCAT positions in the CYP2C slice: ${POSITIONS}; reported missing: ${MISSING}"
check_ge "PharmCAT positions in the CYP2C slice" "$POSITIONS" 20
check "under 10% of them reported missing (${MISSING} of ${POSITIONS})" test "$((MISSING * 10))" -lt "$POSITIONS"
printf '#### PharmCAT on the CYP2C19/CYP2C9 slice\n\n%s of %s PGx positions reported missing\n\n' "$MISSING" "$POSITIONS" >> "$E2E_NOTES"

# --- imputation prep with panel sites --------------------------------------------
cut -f1,2 "$SITES" > "${D}/imputation_sites.tsv"
IMPUTATION_SITES="${D}/imputation_sites.tsv" run_step 14-imputation-prep.sh "$SAMPLE"
check_step_exit 14-imputation-prep.sh
CHR20="${SAMPLE}/imputation/mis_ready/${SAMPLE}_chr20.vcf.gz"
check_eq "0/0 calls at the panel sites on chr20" "$(bcf view -H -g hom -i 'GT="ref"' "$CHR20" 2>/dev/null | wc -l | tr -d ' ')" 5
check_eq "records on chr20 (the panel sites only)" "$(vcf_count "$CHR20")" 5
rm -rf "${D:?}/imputation" "${D}/imputation_sites.tsv"

finish
