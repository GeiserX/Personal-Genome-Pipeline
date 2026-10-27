#!/usr/bin/env bash
# Step 41's platelet (HPA) and neutrophil (HNA) antigens, through the real
# step and its pinned images. The fixture's slices hold none of these genes,
# so the case writes a synthetic VCF with planted genotypes at the sites of
# assets/bloodtype_markers.tsv (no BAM: every site with no record is read as
# the reference allele, confidence medium) and checks the rows:
#   - HPA-1 1/1 is HPA-1(a-b+), and the note says the person lacks HPA-1a,
#     with FNAIT and anti-HPA-1a;
#   - HPA-6 0/1, HNA-3 0/1 and HNA-4 1/1 flip those systems;
#   - four FCGR3B markers 0/1 (with c.194 0/1) give HNA-1a/1b at confidence low;
#   - control: the same step on a VCF with none of these records gives the
#     reference phenotypes, so the checks above see the planted genotypes.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
MARKERS="${REPO}/assets/bloodtype_markers.tsv"

# vcf OUT NAME:GT...: a one-sample VCF with those genotypes at the table's sites.
vcf() {
  local out=$1; shift
  { printf '##fileformat=VCFv4.2\n'
    for c in chr1 chr5 chr6 chr7 chr9 chr16 chr17 chr18 chr19; do printf '##contig=<ID=%s>\n' "$c"; done
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
    awk -F'\t' -v want="$*" 'BEGIN {n = split(want, w, " "); for (i = 1; i <= n; i++) {split(w[i], p, ":"); gt[p[1]] = p[2]}}
      !/^#/ && $4 == "site" && ($6 in gt) {print $1 "\t" $3 "\t.\t" $7 "\t" $8 "\t50\tPASS\t.\tGT\t" gt[$6]}' "$MARKERS" |
      sort -k1,1V -k2,2n
  } > "${T}/in.vcf"
  docker run --rm -u "$(id -u):$(id -g)" -v "${E2E_WORK}:${E2E_WORK}" "$BCFTOOLS_IMAGE" \
    bcftools view -Oz -o "$out" "${T}/in.vcf"
}
# sample NAME VCF: a sample folder with that VCF and no BAM.
sample() { rm -rf "${G:?}/$1"; mkdir -p "${G}/$1/vcf"; cp "$2" "${G}/$1/vcf/$1.vcf.gz"; }
cell() {
  awk -F'\t' -v s="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $1 == s {print $c[k]}' \
    "${G}/$1/bloodtype/$1_bloodtype.tsv" 2>/dev/null
}

# --- planted -----------------------------------------------------------------------
vcf "${T}/planted.vcf.gz" HPA1_176:1/1 HPA6_1544:0/1 HNA3_461:0/1 HNA4_230:1/1 \
  FCGR3B_108:0/1 FCGR3B_114:0/1 FCGR3B_194:0/1 FCGR3B_244:0/1 FCGR3B_316:0/1
check_eq "planted records in the VCF" \
  "$(docker run --rm -v "${E2E_WORK}:${E2E_WORK}" "$BCFTOOLS_IMAGE" bcftools view -H "${T}/planted.vcf.gz" | wc -l | tr -d ' ')" 9
P="${SAMPLE}plt"
sample "$P" "${T}/planted.vcf.gz"
run_step 41-bloodtype.sh "$P"; check_step_exit 41-bloodtype.sh
grep -E '^(HPA|HNA)' "${G}/${P}/bloodtype/${P}_bloodtype.tsv" 2>/dev/null | cut -f1-4
check_eq "a row per system, the HPA and HNA ones included" \
  "$(awk 'END {print NR - 1}' "${G}/${P}/bloodtype/${P}_bloodtype.tsv" 2>/dev/null)" 22
check_eq "HPA-1 1/1: HPA-1(a-b+), 1b/1b" "$(cell "$P" HPA-1 phenotype) $(cell "$P" HPA-1 genotype)" "HPA-1(a-b+) 1b/1b"
check "and the note says the person lacks HPA-1a, with FNAIT and anti-HPA-1a" \
  has '^lacks HPA-1a: .*FNAIT.*anti-HPA-1a' "$(cell "$P" HPA-1 note)"
check_eq "HPA-6 0/1: HPA-6(a+b+)" "$(cell "$P" HPA-6 phenotype)" "HPA-6(a+b+)"
check_eq "HNA-3 0/1: HNA-3(a+b+)" "$(cell "$P" HNA-3 phenotype)" "HNA-3(a+b+)"
check_eq "HNA-4 1/1: HNA-4(a-b+)" "$(cell "$P" HNA-4 phenotype)" "HNA-4(a-b+)"
check "and the note gives the neutrophil meaning" has '^lacks HNA-4a: .*neonatal alloimmune neutropenia' \
  "$(cell "$P" HNA-4 note)"
check_eq "HNA-1 from FCGR3B: a+b+, FCGR3B*01/FCGR3B*02, low" \
  "$(cell "$P" HNA-1 phenotype) $(cell "$P" HNA-1 genotype) $(cell "$P" HNA-1 confidence)" \
  "HNA-1(a+b+c-) FCGR3B*01/FCGR3B*02 low"
check "the evidence file has the HNA-1 sites as called" \
  test "$(awk -F'\t' '$2 ~ /^FCGR3B_/ && $9 == "called"' "${G}/${P}/bloodtype/${P}_bloodtype_evidence.tsv" | wc -l)" -eq 5

# --- control: none of these records ------------------------------------------------
vcf "${T}/none.vcf.gz" KEL_578:0/1
R="${SAMPLE}plr"
sample "$R" "${T}/none.vcf.gz"
run_step 41-bloodtype.sh "$R"; check_step_exit 41-bloodtype.sh
check_eq "control: HPA-1 is HPA-1(a+b-) without the planted record" "$(cell "$R" HPA-1 phenotype)" "HPA-1(a+b-)"
check_eq "control: HNA-4 is HNA-4(a+b-)" "$(cell "$R" HNA-4 phenotype)" "HNA-4(a+b-)"
check_eq "control: HNA-3 is GRCh38's HNA-3(a-b+)" "$(cell "$R" HNA-3 phenotype)" "HNA-3(a-b+)"

finish
