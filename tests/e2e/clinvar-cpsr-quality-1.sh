#!/usr/bin/env bash
# Step 06's call floor, the gnomAD frequency the reports put beside each hit,
# and the cDNA read check on CPSR's P/LP rows, with the real tools.
#   1. The sample as case 30 left it: the hits file carries the default floor
#      line, every hit has its GQ and DP, and bin/clinvar_hits.awk and the
#      summary mark the same hits low confidence.
#   2. Step 06 again on a copy of the sample with CLINVAR_MIN_GQ=99 and
#      CLINVAR_MIN_DP=999: the floor line says so, every hit is marked, the
#      step counts them, and the summary's headline count leaves them out.
#   3. A VEP file written for that copy, bgzipped and indexed by bcftools
#      (htslib), with one record per hit and filler records so the file spans
#      many BGZF blocks: the summary reads each hit's gnomAD AF through the
#      .tbi, and the 3% hit is marked common in both reports.
#   4. bin/clip_check.py --cpsr in PYPGX_IMAGE on a CPSR table holding the
#      planted hit as a P/LP row (CPSR's contig name, no "chr") and a benign
#      row: on the copy with spliced alt reads that clinvar-clip-check-1 left
#      the row is flagged, on the sample's own reads it is clear, and the
#      benign row is not checked.
. "$(dirname "$0")/lib.sh"

CHR=$(planted chrom)
POS=$(planted pos)
REF=$(planted ref)
ALT=$(planted alt)
in_user() {
  local image=$1; shift
  docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$image" "$@"
}
summary_py() {  # summary_py <sample> <python expression over d, the summary>
  python3 - "$1" "$2" "${REPO}" "${GENOME_DIR}" <<'PY'
import os, sys
sys.path.insert(0, os.path.join(sys.argv[3], "bin"))
import collect_summary
d = collect_summary.collect(sys.argv[1], os.path.join(sys.argv[4], sys.argv[1]))
print(eval(sys.argv[2]))
PY
}

# --- 1. the sample as case 30 left it ------------------------------------------
H="${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf"
check "case 30 left a hits file" test -s "$H"
check_eq "the hits file carries the default call floor" "$(grep -c '^##PGP_CALL_FLOOR=GQ=30,DP=10$' "$H")" 1
N=$(grep -vc '^#' "$H" || true)
AWK_ROWS=$(awk -f "${REPO}/bin/clinvar_hits.awk" "$H")
check_eq "bin/clinvar_hits.awk prints 13 columns on every hit" \
  "$(awk -F'\t' 'NF == 13' <<< "$AWK_ROWS" | wc -l | tr -d ' ')" "$N"
check_eq "every hit has a GQ and a DP from DeepVariant" \
  "$(awk -F'\t' '$11 ~ /^[0-9]+$/ && $12 ~ /^[0-9]+$/' <<< "$AWK_ROWS" | wc -l | tr -d ' ')" "$N"
check_eq "the summary marks the same hits low confidence as the awk" \
  "$(summary_py "$SAMPLE" "' '.join(f\"{h['chrom']}:{h['pos']}\" for h in sorted(d['sections']['clinvar']['data']['hits'], key=lambda h: (h['chrom'], h['pos'])) if h['low_confidence'])")" \
  "$(awk -F'\t' '$13 == "yes" {print $2 ":" $3}' <<< "$AWK_ROWS" | sort -t: -k1,1 -k2,2n | tr '\n' ' ' | sed 's/ $//')"

# --- 2. step 06 on a copy with a floor no call reaches -------------------------
T="${SAMPLE}floor"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/vcf"
ln "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${GENOME_DIR}/${T}/vcf/${T}.vcf.gz"
ln "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${GENOME_DIR}/${T}/vcf/${T}.vcf.gz.tbi"
CLINVAR_MIN_GQ=99 CLINVAR_MIN_DP=999 run_step 06-clinvar-screen.sh "$T"
check_step_exit 06-clinvar-screen.sh
TH="${GENOME_DIR}/${T}/clinvar/${T}_clinvar_hits.vcf"
TN=$(grep -vc '^#' "$TH" || true)
check_ge "the copy has the sample's hits" "$TN" 1
check_eq "the floor line says GQ 99 and DP 999" "$(grep -c '^##PGP_CALL_FLOOR=GQ=99,DP=999$' "$TH")" 1
check "the step counts every hit low confidence" \
  has "Low-confidence calls among them \(GQ below 99 or DP below 999\): ${TN}$" "$(cat "$STEP_LOG")"
check "the step marks the planted hit" has "${CHR}:${POS} ${REF}>${ALT} .*LOW-CONFIDENCE CALL" "$(cat "$STEP_LOG")"
check_eq "the summary leaves every hit out of the headline count and counts them apart" \
  "$(summary_py "$T" "(d['sections']['clinvar']['data']['count'], d['sections']['clinvar']['data']['low_confidence'])")" \
  "(0, ${TN})"
CLINVAR_MIN_GQ=x run_step 06-clinvar-screen.sh "$T"
check "a floor that is not a whole number stops the step" test "$STEP_RC" -ne 0
check "and says which" has "CLINVAR_MIN_GQ must be a whole number" "$(cat "$STEP_LOG")"

# --- 3. the gnomAD AF through an htslib index ----------------------------------
mkdir -p "${GENOME_DIR}/${T}/vep"
python3 - "${GENOME_DIR}/${T}/clinvar/${T}_clinvar_hits.tsv" "${GENOME_DIR}/${T}/vep/plain.vcf" \
  "${CASE_TMP}/want.tsv" "$CHR" "$POS" <<'PY'
import sys
hits, out, want, pchrom, ppos = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
rows = [l.rstrip("\n").split("\t") for l in open(hits)][1:]
recs = {}
expect = []
for i, r in enumerate(rows):
    chrom, pos, ref, alt = r[0], int(r[1]), r[2], r[3]
    af = 0.03 if (chrom, pos) == (pchrom, ppos) else round(0.0001 * (i + 1), 6)
    snv = len(ref) == 1 and len(alt) == 1
    allele = alt if snv or alt[:1] != ref[:1] else (alt[1:] or "-")
    recs[(chrom, pos, ref, alt)] = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tCSQ={allele}|x|LOW|G|{af}|{af / 2}\tGT\t0/1"
    expect.append(f"{chrom}:{pos}\t{af}")
# Filler records around the planted hit (2 Mb each side) and the others
# (50 kb), so the file spans many BGZF blocks and many 16 kb windows of the
# linear index.
at = {(c, p) for c, p, _, _ in recs}
for chrom, pos in sorted(at):
    span = 2_000_000 if (chrom, pos) == (pchrom, ppos) else 50_000
    for p in range(max(1, pos - span), pos + span, 97):
        if (chrom, p) not in at:
            recs.setdefault((chrom, p, "A", "C"), f"{chrom}\t{p}\t.\tA\tC\t50\tPASS\tCSQ=C|x|LOW|G|0.4|0.4\tGT\t0/1")
order = {c: i for i, c in enumerate(sorted({k[0] for k in recs}, key=lambda c: (len(c), c)))}
with open(out, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for c in order:
        f.write(f"##contig=<ID={c}>\n")
    f.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
            'Format: Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF|gnomADg_AF">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    for k in sorted(recs, key=lambda k: (order[k[0]], k[1])):
        f.write(recs[k] + "\n")
open(want, "w").write("\n".join(sorted(expect)) + "\n")
print(f"{len(recs)} VEP records, {len(expect)} at hits")
PY
in_user "$BCFTOOLS_IMAGE" bcftools view --no-version -Oz -o "${T}/vep/${T}_vep.vcf.gz" "${T}/vep/plain.vcf"
in_user "$BCFTOOLS_IMAGE" bcftools index -t "${T}/vep/${T}_vep.vcf.gz"
VEPGZ="${GENOME_DIR}/${T}/vep/${T}_vep.vcf.gz"
check "bcftools wrote the VEP file and its .tbi" test -s "${VEPGZ}.tbi"
# A BGZF block holds at most 64 KiB of text, so 2 MB of text is 30 blocks or more.
check_ge "the VEP text spans 30 BGZF blocks or more (bytes)" \
  "$(wc -c < "${GENOME_DIR}/${T}/vep/plain.vcf" | tr -d ' ')" 2000000
rm -f "${GENOME_DIR}/${T}/vep/plain.vcf"
GOT=$(summary_py "$T" "'\n'.join(sorted(f\"{h['chrom']}:{h['pos']}\t{h['af']}\" for h in d['sections']['clinvar']['data']['hits']))")
check_eq "every hit's gnomAD AF is read through the .tbi" "$GOT" "$(cat "${CASE_TMP}/want.tsv")"
check_eq "the planted 3% hit is the one marked common" \
  "$(summary_py "$T" "[f\"{h['chrom']}:{h['pos']}\" for h in d['sections']['clinvar']['data']['hits'] if h['common']]")" \
  "['${CHR}:${POS}']"
python3 "${REPO}/bin/render_report.py" --sample "$T" --sample-dir "${GENOME_DIR}/${T}" --json "${CASE_TMP}/t.json" \
  -o "${CASE_TMP}/t.txt" -o "${CASE_TMP}/t.html"
check "the text report prints the AF and the common mark" \
  has "frequency: gnomAD AF 0.03; common: unlikely to be a high-penetrance finding" "$(cat "${CASE_TMP}/t.txt" 2>/dev/null)"
check "the HTML hits table has the frequency column" has '<th>gnomAD frequency</th>' "$(cat "${CASE_TMP}/t.html" 2>/dev/null)"
check "the summary validates against the schema" \
  python3 "${REPO}/tests/schema/validate.py" "${REPO}/tests/schema/summary.schema.json" "${CASE_TMP}/t.json"

# --- 4. the cDNA read check on CPSR rows ----------------------------------------
C="${SAMPLE}cdna"
SPLICED="${GENOME_DIR}/${C}/aligned/${C}_sorted.bam"
check "clinvar-clip-check-1 left the copy with spliced alt reads" test -s "${SPLICED}.bai"
EXONS="${CASE_TMP}/exons.bed"
P0=$((POS - 1))
B0=$((P0 + 3000))
{
  printf '%s\t%s\t%s\n' "$CHR" "$((P0 - 30))" "$((P0 + 30))"
  printf '%s\t%s\t%s\n' "$CHR" "$B0" "$((B0 + 150))"
} > "$EXONS"
mkdir -p "${GENOME_DIR}/${C}/cpsr"
CLS="${GENOME_DIR}/${C}/cpsr/planted.classification.tsv.gz"
{
  printf 'VAR_ID\tGENOMIC_CHANGE\tSYMBOL\tCLINVAR_CLASSIFICATION\tCPSR_CLASSIFICATION\tCLASSIFICATION\n'
  printf '%s_%s_%s_%s\t%s:g.%s%s>%s\t%s\tPathogenic\tPathogenic\tPathogenic\n' "${CHR#chr}" "$POS" "$REF" "$ALT" \
    "${CHR#chr}" "$POS" "$REF" "$ALT" "$(planted gene)"
  printf '%s_%s_A_G\t%s:g.%sA>G\tOTHER\tBenign\tBenign\tBenign\n' "${CHR#chr}" "$((POS + 500))" "${CHR#chr}" "$((POS + 500))"
} | gzip -c > "$CLS"
cp "$EXONS" "${GENOME_DIR}/${C}/cpsr/exons.bed"
clip() {  # clip <bam relative to GENOME_DIR> <out relative to GENOME_DIR>
  docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w /genome \
    "$PYPGX_IMAGE" python3 /pgp-bin/clip_check.py --cpsr "${C}/cpsr/planted.classification.tsv.gz" \
    --bam "$1" --exons "${C}/cpsr/exons.bed" --out "$2"
}
clip "${C}/aligned/${C}_sorted.bam" "${C}/cpsr/spliced.tsv"
clip "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${C}/cpsr/genomic.tsv"
ROW_S=$(awk -F'\t' 'NR > 1' "${GENOME_DIR}/${C}/cpsr/spliced.tsv" 2>/dev/null)
ROW_G=$(awk -F'\t' 'NR > 1' "${GENOME_DIR}/${C}/cpsr/genomic.tsv" 2>/dev/null)
echo "spliced: ${ROW_S}"
echo "genomic: ${ROW_G}"
check_eq "only the P/LP row is checked, under CPSR's contig name" "$(cut -f1-4 <<< "$ROW_S" | tr '\t' ' ')" \
  "${CHR#chr} ${POS} ${REF} ${ALT}"
check_eq "the planted P/LP row on the spliced reads is flagged" "$(cut -f9 <<< "$ROW_S")" flagged
check_eq "the same row on the sample's own reads is clear" "$(cut -f9 <<< "$ROW_G")" clear
rm -rf "${GENOME_DIR:?}/${C}/cpsr"

finish
