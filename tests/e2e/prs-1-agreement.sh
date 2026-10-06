#!/usr/bin/env bash
# Step 25 runs pgsc_calc; this case ties its sum to a known answer and to the
# plink2 scoring step 25 did before it (cases 21 and 03 left a VCF and a gVCF
# of HG002):
#   - a synthetic score on the chr20 slice: five sites where HG002 matches the
#     reference (GIAB v4.2.1 truth, no variant within 200 bp) with the
#     reference base as effect allele (dosage 2), one such site whose effect
#     allele is another base (dosage 0), and three truth SNVs with their ALT
#     as effect allele, each with its own weight;
#   - pgsc_calc's sum is compared with plink2 --score (cols=+scoresums,
#     no-mean-imputation, the old path) on the genotypes step 25 handed
#     pgsc_calc, and the hom-ref part alone must be 2 x (1+2+4+8+16) = 62;
#   - control: plink2 on the same score with one weight changed must differ
#     from pgsc_calc's sum, so the comparison can fail.
# Writes nothing other cases read; removes what it made.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
D="${GENOME_DIR}/${SAMPLE}"
check "case 21 wrote a gVCF" vcf_ok "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"

# --- sites -----------------------------------------------------------------------
SITES="${CASE_TMP}/sites.tsv"
python3 - "${FIXTURE_DIR}/HG002_truth_chr20.bed" "${FIXTURE_DIR}/HG002_truth_chr20.vcf.gz" \
  "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" > "$SITES" <<'PY'
import gzip, sys
bed, truth, fasta = sys.argv[1:4]
iv = [tuple(int(x) for x in l.split()[1:3]) for l in open(bed) if l.startswith("chr20\t")]
var, snv = [], []
for l in gzip.open(truth, "rt"):
    if not l.startswith("chr20\t"):
        continue
    f = l.rstrip("\n").split("\t")
    var.append(int(f[1]))
    gt = f[9].split(":")[0].replace("|", "/")
    # pgsc_calc drops strand-ambiguous pairs (A/T, C/G), plink2 would not.
    amb = {f[3], f[4]} in ({"A", "T"}, {"C", "G"})
    if len(f[3]) == 1 and len(f[4]) == 1 and not amb and 10_050_000 < int(f[1]) < 10_450_000 and gt in ("0/1", "1/1"):
        snv.append((int(f[1]), f[4], f[3]))
fai = {l.split()[0]: [int(x) for x in l.split()[1:5]] for l in open(fasta + ".fai")}
length, offset, bases, width = fai["chr20"]
def base(pos):
    i = pos - 1
    with open(fasta, "rb") as fh:
        fh.seek(offset + (i // bases) * width + i % bases)
        return fh.read(1).decode().upper()
homref, p = [], 10_050_000
while len(homref) < 6 and p < 10_450_000:
    inside = any(s + 200 < p <= e - 200 for s, e in iv)
    clear = all(abs(v - p) > 200 for v in var)
    b = base(p)
    if inside and clear and b in "ACGT":
        homref.append((p, b))
    p += 7_919
# chrom, pos, effect allele, other allele, weight, kind
for i, (p, b) in enumerate(homref[:5]):
    print(f"20\t{p}\t{b}\t\t{2 ** i}\thomref")
p, b = homref[5]
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
other = next(x for x in "ACGT" if x not in (b, comp[b]))
print(f"20\t{p}\t{other}\t{b}\t3.5\thomref_effect_not_ref")
for (p, alt, ref), w in zip(snv[:: max(1, len(snv) // 3)][:3], (0.25, -1.5, 0.125)):
    print(f"20\t{p}\t{alt}\t{ref}\t{w}\tsnv")
PY
cat "$SITES"
check_eq "score sites (5 hom-ref, 1 hom-ref with another effect allele, 3 truth SNVs)" "$(wc -l < "$SITES" | tr -d ' ')" 9

# The nine ids of assets/pgs_scores.tsv, each the synthetic score (step 25 scores the whole list).
SCORES="${GENOME_DIR}/prs_scores"
rm -rf "$SCORES" "${D}/prs"
mkdir -p "$SCORES"
for id in $(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${REPO}/assets/pgs_scores.tsv"); do
  { printf '#pgs_id=%s\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n' "$id"
    cut -f1-5 "$SITES"; } | gzip -c > "${SCORES}/${id}.txt.gz"
done

# --- pgsc_calc (step 25) ---------------------------------------------------------
ANCESTRY_PANEL=none run_step 25-prs.sh "$SAMPLE"
check_step_exit 25-prs.sh
SUMMARY="${D}/prs/${SAMPLE}_prs_summary.tsv"
col() { awk -F'\t' -v k="$1" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["PGS_ID"] == "PGS000018" {print $c[k]}' "$SUMMARY" 2>/dev/null; }
PGSC_SUM=$(col Score_SUM)
check_eq "pgsc_calc matched every score site" "$(col Variants_Matched)" 9
check_eq "input: genotypes from the gVCF" "$(col Input)" gvcf
check_eq "no percentile without a panel" "$(col Percentile)" NA
check "the raw-score line in the step log" has 'Raw score only' "$(cat "$STEP_LOG")"
check "the label comes from assets/pgs_scores.tsv" test "$(col Condition)" = "Coronary artery disease"
TARGET="${D}/prs/pgsc_calc/target.vcf.gz"
check "step 25 kept the genotypes it gave pgsc_calc" test -s "$TARGET"

# --- the old path: plink2 --score on the same genotypes ---------------------------
plink_sum() {  # plink_sum SCORE_TSV: SCORE1_SUM of chr:pos / allele / weight
  docker run --rm --network none --user "$(id -u):$(id -g)" -v "${CASE_TMP}:/w" -v "$(dirname "$TARGET"):/t:ro" -w /w "$PLINK2_IMAGE" \
    sh -c 'plink2 --vcf /t/target.vcf.gz --make-pgen --out g --set-all-var-ids "@:#" --new-id-max-allele-len 100 \
             --chr 1-22 --allow-extra-chr --output-chr chrM >/dev/null &&
           plink2 --pfile g --score "$1" 1 2 3 ignore-dup-ids no-mean-imputation cols=+scoresums --out s --allow-extra-chr >/dev/null' \
    _ "/w/$1" >/dev/null 2>&1 || { echo "plink2 failed"; return; }
  awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} NR == 2 {print $c["SCORE1_SUM"]}' "${CASE_TMP}/s.sscore"
}
awk -F'\t' -v OFS='\t' '{print "chr" $1 ":" $2, $3, $5}' "$SITES" > "${CASE_TMP}/score.tsv"
PLINK_SUM=$(plink_sum score.tsv)
awk -F'\t' -v OFS='\t' '$6 == "homref" {print "chr" $1 ":" $2, $3, $5}' "$SITES" > "${CASE_TMP}/homref.tsv"
HOMREF_SUM=$(plink_sum homref.tsv)
echo "pgsc_calc sum ${PGSC_SUM:-none}; plink2 sum ${PLINK_SUM:-none}; plink2 hom-ref part ${HOMREF_SUM:-none}"
close() { awk -v a="$1" -v b="$2" 'BEGIN { if (a == "" || b == "" || a !~ /^-?[0-9.e+-]+$/ || b !~ /^-?[0-9.e+-]+$/) exit 1; d = a - b; exit !(d < 1e-6 && d > -1e-6) }'; }
check "the hom-ref sites alone sum to 2 x (1+2+4+8+16) = 62 (plink2: ${HOMREF_SUM:-none})" close "${HOMREF_SUM:-}" 62
check "pgsc_calc's sum (${PGSC_SUM:-none}) equals plink2's (${PLINK_SUM:-none})" close "${PGSC_SUM:-}" "${PLINK_SUM:-}"
# Control: one weight changed, so the comparison above is shown able to fail.
awk -F'\t' -v OFS='\t' 'NR == 7 {$3 = $3 + 10} {print}' "${CASE_TMP}/score.tsv" > "${CASE_TMP}/score_changed.tsv"
CHANGED_SUM=$(plink_sum score_changed.tsv)
check "control: a changed weight gives a different sum (${CHANGED_SUM:-none} vs ${PGSC_SUM:-none})" \
  bash -c '! awk -v a="$1" -v b="$2" "BEGIN { d = a - b; exit !(d < 1e-6 && d > -1e-6) }"' _ "${CHANGED_SUM:-0}" "${PGSC_SUM:-0}"
printf '#### PRS: pgsc_calc against the old plink2 path (synthetic nine-site score on chr20)\n\npgsc_calc %s, plink2 %s, hom-ref part %s (want 62), changed-weight control %s\n\n' \
  "${PGSC_SUM:-none}" "${PLINK_SUM:-none}" "${HOMREF_SUM:-none}" "${CHANGED_SUM:-none}" >> "$E2E_NOTES"

rm -rf "$SCORES" "${D:?}/prs"
finish
