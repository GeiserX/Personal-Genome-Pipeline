#!/usr/bin/env bash
# APOL1 G2 (rs71785313) beside step 42's APOL1 rule (bin/genotype_rules_indels.py),
# with the real bcftools images:
#   - dbSNP orientation, online: the record pinned in the module (DBSNP) is what
#     dbSNP's GRCh38 VCF of the build versions.env pins holds for rs71785313 (same
#     position and REF, the module's ALT among dbSNP's); control: the same compare on
#     a copy with REF and ALT swapped fails;
#   - the report's reader through a real .tbi: a bgzipped, indexed gVCF (bcftools
#     view -Oz, bcftools index -t) with a G2 deletion, read by collect_summary's
#     genotype_rules, gives G1/G2 beside a G0/G1 rule line; control: the same
#     folder with a reference block over G2 gives G0/G1.
# Writes one folder under GENOME_DIR, removed at the end; no other case reads it.
. "$(dirname "$0")/lib.sh"

URL="https://ftp.ncbi.nih.gov/snp/archive/${CONSUMER_DBSNP_BUILD:?}/VCF/GCF_000001405.40.gz"
PINNED=$(python3 -c "import sys; sys.path.insert(0, '${REPO}/bin'); import genotype_rules_indels as m; \
print(m.DBSNP[0], m.DBSNP[2], m.DBSNP[3], m.DBSNP[4], sep='\t')")
echo "pinned: ${PINNED}"

# same_record PINNED LIVE: the pinned rsid, POS and REF equal dbSNP's, and the pinned ALT is one of dbSNP's.
same_record() {
  awk -F'\t' 'NR == 1 { id = $1; pos = $2; ref = $3; alt = $4; next }
              $1 == id && $2 == pos && $3 == ref && index("," $4 ",", "," alt ",") { ok = 1 }
              END { exit ok ? 0 : 1 }' <(echo "$1") <(echo "$2")
}

LIVE=$(docker run --rm "${BCFTOOLS_HTTPS_IMAGE:?}" bcftools query -r NC_000022.11:36265990-36266010 \
  -i 'ID="rs71785313"' -f '%ID\t%POS\t%REF\t%ALT\n' "$URL" 2>"${CASE_TMP}/dbsnp.err")
echo "dbSNP ${CONSUMER_DBSNP_BUILD}: ${LIVE:-<nothing>}"
check_eq "dbSNP's VCF has one rs71785313 record near chr22:36265995" "$(grep -c '^rs71785313' <<< "$LIVE")" 1
check "the module's pinned G2 record is dbSNP's" same_record "$PINNED" "$LIVE"
SWAPPED=$(awk -F'\t' -v OFS='\t' '{ t = $3; $3 = $4; $4 = t; print }' <<< "$PINNED")
check "control: a pinned record with REF and ALT swapped does not match dbSNP" \
  test "$(same_record "$SWAPPED" "$LIVE"; echo $?)" -ne 0

# --- the reader through a real .tbi ---------------------------------------------------
D="${GENOME_DIR}/genotype-rules-g2"
in_genome "$BCFTOOLS_IMAGE" rm -rf genotype-rules-g2
mkdir -p "$D"
for s in G2HET G2REF; do
  mkdir -p "${D}/${s}/vcf" "${D}/${s}/traits"
  python3 - "$REPO" "${D}/${s}" "$s" <<'PY'
import os, sys
repo, d, s = sys.argv[1:]
sys.path.insert(0, os.path.join(repo, "bin"))
import traits_report as tr
rules = tr.read_rules(os.path.join(repo, "assets", "genotype_rules.tsv"))
g = os.path.join(d, "g.tsv")
with open(g, "w") as f:   # G0/G1 at the two G1 sites
    f.write("chr22\t36265860\tA\tG\tA/G\t40\nchr22\t36265988\tT\tG\tT/G\t40\n")
tr.write_tsv(os.path.join(d, "traits", f"{s}_genotype_rules.tsv"), tr.rule_lines(rules, tr.read_genotypes(g)),
             tr.RULE_OUT_COLUMNS)
os.remove(g)
body = ["chr22\t36200000\t.\tA\t<*>\t0\t.\tEND=36265994\tGT:GQ\t0/0:50"]
if s == "G2HET":
    body.append("chr22\t36265995\trs71785313\tAATAATT\tA,<*>\t30\tPASS\t.\tGT:GQ\t0/1:30")
    body.append("chr22\t36266002\t.\tT\t<*>\t0\t.\tEND=36300000\tGT:GQ\t0/0:50")
else:
    body.append("chr22\t36265995\t.\tA\t<*>\t0\t.\tEND=36300000\tGT:GQ\t0/0:50")
with open(os.path.join(d, "vcf", f"{s}.g.vcf"), "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=chr22,length=50818468>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n"
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{s}\n" + "\n".join(body) + "\n")
PY
  in_genome "$BCFTOOLS_IMAGE" bash -c "bcftools view -Oz -o 'genotype-rules-g2/${s}/vcf/${s}.g.vcf.gz' \
    'genotype-rules-g2/${s}/vcf/${s}.g.vcf' && bcftools index -t 'genotype-rules-g2/${s}/vcf/${s}.g.vcf.gz' \
    && rm 'genotype-rules-g2/${s}/vcf/${s}.g.vcf'"
  check "${s}: bcftools wrote the gVCF and its .tbi" test -s "${D}/${s}/vcf/${s}.g.vcf.gz.tbi"
done

read_result() {
  python3 - "$REPO" "${D}/$1" "$1" <<'PY'
import os, sys
repo, d, s = sys.argv[1:]
sys.path.insert(0, os.path.join(repo, "bin"))
import collect_summary as cs
r = next(x for x in cs.genotype_rules(d, s) if x["rule"] == "APOL1")
print(r["result"], r["g2"]["result"], r["g2"]["source"], r["g2"]["status"], sep="\t")
PY
}
HET=$(read_result G2HET)
REF=$(read_result G2REF)
echo "G2HET: ${HET}"
echo "G2REF: ${REF}"
check_eq "G2HET: the rule's G0/G1 with one G2 copy reads G1/G2 from the indexed gVCF" "$HET" \
  "$(printf 'G0/G1\tG1/G2\tvcf/G2HET.g.vcf.gz\tcalled')"
check_eq "control: G2REF, a reference block over G2, reads G0/G1" "$REF" \
  "$(printf 'G0/G1\tG0/G1\tvcf/G2REF.g.vcf.gz\tcalled')"

in_genome "$BCFTOOLS_IMAGE" rm -rf genotype-rules-g2
finish
