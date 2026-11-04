#!/usr/bin/env bash
# Step 06's genome-wide cDNA scan (bin/cdna_scan.py) and its always-report
# list of known deep intronic sites, in bash and in Nextflow, on copies of the
# fixture sample (its VCF, its gVCF from case 21 and its BAM).
#   1. A copy whose BAM also holds spliced reads across three made-up exons
#      of one stretch of chr20 (exon 1 to 2, 2 to 3 and 1 to 3, as a DNA
#      aligner writes cDNA: clipped at the exon end, the rest an SA alignment
#      at the next exon's start). The scan lists that one locus with 3
#      junctions. The default list (assets/deep_intronic_sites.tsv) gets a row
#      per site, each with a call and the reason when there is none, and both
#      reports print every row.
#   2. The same copy with a list of its own (CLINVAR_ALWAYS_REPORT): the
#      planted ClinVar SNV gives its genotype, a site inside a gVCF reference
#      block gives "ref" with the block's GQ ("ref (low confidence)" below
#      the GQ floor), a site no record covers "no call".
#   3. Red control: the real BAM alone with the same exon file lists nothing,
#      and without the gVCF the block site is "no call" with that reason.
#   4. Nextflow (--tools clinvar,html_report, the row's gvcf column, the same
#      exons): CLINVAR_SCREEN writes the same always-report table as step 06,
#      CLINVAR_CDNA_CHECK the same scan table, and HTML_REPORT shows both.
. "$(dirname "$0")/lib.sh"

GENCODE_RELEASE=$(sed -n 's/^GENCODE_RELEASE=//p' "${REPO}/scripts/lib/common.sh")
EXONS="${GENOME_DIR}/reference/gencode.v${GENCODE_RELEASE}.basic.exons.bed"
trap 'rm -f "$EXONS"' EXIT
CHR=$(planted chrom)
POS=$(planted pos)
GVCF="${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
check "case 21 left the sample's gVCF" nonempty "$GVCF"
in_user() { docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$SAMTOOLS_IMAGE" samtools "$@"; }

# Three exons 10 kb past the planted hit, away from its reads.
X=$((POS + 10000))
E1S=$X; E1E=$((X + 150)); E2S=$((X + 2000)); E2E=$((X + 2150)); E3S=$((X + 4000)); E3E=$((X + 4200))
printf '%s\t%s\t%s\n' "$CHR" "$E1S" "$E1E" "$CHR" "$E2S" "$E2E" "$CHR" "$E3S" "$E3E" > "$EXONS"

# copy NAME [nogvcf]: a sample folder of hard links to the fixture sample's files.
copy() {
  rm -rf "${GENOME_DIR:?}/$1"
  mkdir -p "${GENOME_DIR}/$1/vcf" "${GENOME_DIR}/$1/aligned"
  ln "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${GENOME_DIR}/$1/vcf/$1.vcf.gz"
  ln "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${GENOME_DIR}/$1/vcf/$1.vcf.gz.tbi"
  if [ "${2:-}" != nogvcf ]; then
    ln "${GENOME_DIR}/${GVCF}" "${GENOME_DIR}/$1/vcf/$1.g.vcf.gz"
    ln "${GENOME_DIR}/${GVCF}.tbi" "${GENOME_DIR}/$1/vcf/$1.g.vcf.gz.tbi"
  fi
}

# --- the spliced copy: the real BAM plus 18 spliced reads ----------------------------
T="${SAMPLE}cs"
copy "$T"
in_user faidx reference/GRCh38_no_alt_analysis_set.fasta "${CHR}:$((X + 1))-$((E3E + 200))" > "${CASE_TMP}/region.fa"
in_user view -H "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" > "${CASE_TMP}/spliced.sam"
python3 - "${CASE_TMP}/region.fa" "$CHR" "$X" "$E1E" "$E2S" "$E2E" "$E3S" >> "${CASE_TMP}/spliced.sam" <<'PY'
import sys
fa, chrom = sys.argv[1], sys.argv[2]
x, e1e, e2s, e2e, e3s = (int(v) for v in sys.argv[3:8])
g = "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()
base = lambda a, b: g[a - x:b - x]
for name, end, start in (("j12", e1e, e2s), ("j23", e2e, e3s), ("j13", e1e, e3s)):
    for i in range(6):
        a = 60 + 4 * i
        c = 100 - a
        print(f"{name}_{i}\t0\t{chrom}\t{end - a + 1}\t60\t{a}M{c}S\t*\t0\t0\t{base(end - a, end) + base(start, start + c)}"
              f"\t{'I' * 100}\tSA:Z:{chrom},{start + 1},+,{a}S{c}M,60,0;")
PY
cp "${CASE_TMP}/spliced.sam" "${GENOME_DIR}/${T}/aligned/spliced.sam"
in_user sort -o "${T}/aligned/spliced.bam" "${T}/aligned/spliced.sam"
in_user merge -f -c -p -o "${T}/aligned/${T}_sorted.bam" "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${T}/aligned/spliced.bam"
in_user index "${T}/aligned/${T}_sorted.bam"
rm -f "${GENOME_DIR}/${T}/aligned/spliced.sam" "${GENOME_DIR}/${T}/aligned/spliced.bam"
check_eq "the spliced copy has the 18 planted reads" \
  "$(in_user view -c "${T}/aligned/${T}_sorted.bam" "${CHR}:$((E1E - 10))-$((E1E))" 2>/dev/null | tr -d ' ')" \
  "$(( $(in_user view -c "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${CHR}:$((E1E - 10))-$((E1E))" 2>/dev/null | tr -d ' ') + 12 ))"

# A site inside a gVCF reference block of the chr20 slice, and the block's GQ.
BLOCK=$(bcf query -r "${CHR}:$((POS - 50000))-$((POS - 1000))" -f '%POS\t%END\t[%GT]\t[%GQ]\t%ALT\n' "$GVCF" 2>/dev/null \
  | awk -F'\t' '$3 == "0/0" && $2 - $1 >= 20 && $5 ~ /^<(\*|NON_REF)>$/ {print $1 + 10 "\t" $4; exit}')
BPOS=$(cut -f1 <<< "$BLOCK")
BGQ=$(cut -f2 <<< "$BLOCK")
check "a reference block of 20 bases or more on ${CHR} before the planted hit (${BPOS:-none}, GQ ${BGQ:-none})" test -n "$BPOS"
BREF=$(sam faidx reference/GRCh38_no_alt_analysis_set.fasta "${CHR}:${BPOS}-${BPOS}" 2>/dev/null | sed 1d | tr -d '\n' | tr '[:lower:]' '[:upper:]')
BALT=$( [ "$BREF" = A ] && echo C || echo A )
SITES="${CASE_TMP}/sites.tsv"
{
  printf '# the case list\nchrom\tpos\tref\talt\tgene\thgvs\tclinvar_id\n'
  printf '%s\t%s\t%s\t%s\tPLANTED\tc.1\t1\n' "$CHR" "$POS" "$(planted ref)" "$(planted alt)"
  printf '%s\t%s\t%s\t%s\tBLOCK\tc.2\t2\n' "$CHR" "$BPOS" "$BREF" "$BALT"
  printf 'chr2\t1000\tA\tG\tNOWHERE\tc.3\t3\n'
} > "$SITES"

# --- 2. the case's own list ----------------------------------------------------------
CLINVAR_ALWAYS_REPORT="$SITES" run_step 06-clinvar-screen.sh "$T"
check_step_exit 06-clinvar-screen.sh
DEEP="${GENOME_DIR}/${T}/clinvar/${T}_clinvar_deep_intronic.tsv"
call() { awk -F'\t' -v g="$1" '$5 == g {print $8 "|" $9}' "$DEEP" 2>/dev/null; }
check "the planted SNV gives its genotype with GQ and DP" has '^(het|hom)\|GQ [0-9]+, DP [0-9]+$' "$(call PLANTED)"
# Below step 6's call floor (CLINVAR_MIN_GQ, default 30) the block is a low-confidence ref.
if [ "${BGQ:-0}" -ge 30 ] 2>/dev/null; then WANT="ref|GQ ${BGQ} (gVCF block)"
else WANT="ref (low confidence)|GQ ${BGQ} (gVCF block), below the call floor of GQ 30"; fi
check_eq "the site inside a gVCF reference block gives ref with the block's GQ" "$(call BLOCK)" "$WANT"
check_eq "a site no record covers gives no call, and why" "$(call NOWHERE)" "no call|no record in the VCF or the gVCF covers the site"

# --- 1. the default list and the scan -------------------------------------------------
run_step 06-clinvar-screen.sh "$T"
check_step_exit 06-clinvar-screen.sh
SCAN="${GENOME_DIR}/${T}/clinvar/${T}_clinvar_cdna_scan.tsv"
echo "scan table:"; sed 's/^/    | /' "$SCAN" 2>/dev/null
check_eq "the scan lists one locus" "$(awk -F'\t' '!/^#/ && $1 != "chrom"' "$SCAN" 2>/dev/null | wc -l | tr -d ' ')" 1
check_eq "the locus: exon 1's end to exon 3's start, 3 junctions, 18 reads" \
  "$(awk -F'\t' -v OFS=' ' '!/^#/ && $1 != "chrom" {print $1, $2, $3, $4, $5}' "$SCAN" 2>/dev/null)" \
  "${CHR} ${E1E} $((E3S + 1)) 3 18"
check "the step's log names the locus" has "SPLICED READS ${CHR}:${E1E}-$((E3S + 1)): 3 exon junctions, 18 reads" "$(cat "$STEP_LOG")"
LIST="${REPO}/assets/deep_intronic_sites.tsv"
N_SITES=$(awk -F'\t' '!/^#/ && $1 != "chrom"' "$LIST" | wc -l | tr -d ' ')
check_eq "one always-report row per site of the default list" "$(awk 'NR > 1' "$DEEP" 2>/dev/null | wc -l | tr -d ' ')" "$N_SITES"
check_eq "the rows are the list's sites, in its order" \
  "$(awk -F'\t' 'NR > 1 {print $5 ":" $6}' "$DEEP" 2>/dev/null | paste -sd, -)" \
  "$(awk -F'\t' '!/^#/ && $1 != "chrom" {print $5 ":" $6}' "$LIST" | paste -sd, -)"
check_eq "every row has a call and a detail" "$(awk -F'\t' 'NR > 1 && ($8 == "" || $9 == "")' "$DEEP" 2>/dev/null | wc -l | tr -d ' ')" 0
python3 "${REPO}/bin/render_report.py" --sample "$T" --sample-dir "${GENOME_DIR}/${T}" --json "${CASE_TMP}/t.json" \
  -o "${CASE_TMP}/t.txt" -o "${CASE_TMP}/t.html"
check "the text report has the always-report list" has "Known deep intronic sites \(always reported\)" "$(cat "${CASE_TMP}/t.txt" 2>/dev/null)"
check_eq "the text report prints every site" \
  "$(grep -cE '^    (USH2A|CFTR|ATM|CEP290|GLA) +NM_' "${CASE_TMP}/t.txt" 2>/dev/null)" "$N_SITES"
check "the text report prints the scan" has "cDNA scan of every exon: 1 loci" "$(cat "${CASE_TMP}/t.txt" 2>/dev/null)"
check_eq "the HTML report has a row per site" "$(grep -c '<tr><td>[A-Z0-9]*</td><td>NM_' "${CASE_TMP}/t.html" 2>/dev/null)" "$N_SITES"
check "the HTML report has the scan's locus" has "<td>${CHR}:${E1E}-$((E3S + 1))</td><td>3</td><td>18</td>" "$(cat "${CASE_TMP}/t.html" 2>/dev/null)"
cp "$DEEP" "${CASE_TMP}/bash_deep.tsv"
cp "$SCAN" "${CASE_TMP}/bash_scan.tsv"

# --- 3. red control: the real BAM alone, and no gVCF -----------------------------------
R="${SAMPLE}csr"
copy "$R" nogvcf
ln "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${R}/aligned/${R}_sorted.bam"
ln "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${R}/aligned/${R}_sorted.bam.bai"
CLINVAR_ALWAYS_REPORT="$SITES" run_step 06-clinvar-screen.sh "$R"
check_step_exit 06-clinvar-screen.sh
check_eq "control: the genomic reads alone list no locus" \
  "$(awk -F'\t' '!/^#/ && $1 != "chrom"' "${GENOME_DIR}/${R}/clinvar/${R}_clinvar_cdna_scan.tsv" 2>/dev/null | wc -l | tr -d ' ')" 0
check "control: the scan did read the exons" has '^# cdna_scan: read 3 exon regions on 1 contigs; 0 loci listed$' \
  "$(head -n 1 "${GENOME_DIR}/${R}/clinvar/${R}_clinvar_cdna_scan.tsv" 2>/dev/null)"
DEEP="${GENOME_DIR}/${R}/clinvar/${R}_clinvar_deep_intronic.tsv"
check_eq "control: without the gVCF the block site is no call, and why" "$(call BLOCK)" \
  "no call|no record in the VCF, and no gVCF to tell a reference call from no call"

# --- 4. Nextflow ----------------------------------------------------------------------
if command -v nextflow >/dev/null; then
  VCF="${GENOME_DIR}/${T}/vcf/${T}.vcf.gz"
  G="${GENOME_DIR}/${T}/vcf/${T}.g.vcf.gz"
  B="${GENOME_DIR}/${T}/aligned/${T}_sorted.bam"
  SHEET="${CASE_TMP}/sheet.csv"
  printf 'sample,vcf,vcf_index,bam,bam_index,gvcf,gvcf_index\n%s,%s,%s,%s,%s,%s,%s\n' \
    "$T" "$VCF" "${VCF}.tbi" "$B" "${B}.bai" "$G" "${G}.tbi" > "$SHEET"
  cp "$EXONS" "${CASE_TMP}/exons.bed"
  NF_OUT="${CASE_TMP}/nf"
  echo "+ nextflow run main.nf --tools clinvar,html_report --gencode_exons"
  (
    cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-construct" --input "$SHEET" \
      --reference "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" --tools clinvar,html_report \
      --clinvar "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
      --clinvar_index "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
      --gencode_exons "${CASE_TMP}/exons.bed" \
      --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB
  ) > "${CASE_TMP}/nf.log" 2>&1
  NF_RC=$?
  sed 's/^/    | /' "${CASE_TMP}/nf.log" | grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists'
  check_eq "nextflow run exits 0" "$NF_RC" 0
  check "CLINVAR_SCREEN writes the same always-report table as step 06" \
    cmp -s "${NF_OUT}/${T}/clinvar/${T}_clinvar_deep_intronic.tsv" "${CASE_TMP}/bash_deep.tsv"
  check "CLINVAR_CDNA_CHECK writes the same scan table as step 06" \
    cmp -s "${NF_OUT}/${T}/clinvar/${T}_clinvar_cdna_scan.tsv" "${CASE_TMP}/bash_scan.tsv"
  HTML="${NF_OUT}/${T}/${T}_report.html"
  check_eq "the Nextflow report has a row per site" "$(grep -c '<tr><td>[A-Z0-9]*</td><td>NM_' "$HTML" 2>/dev/null)" "$N_SITES"
  check "the Nextflow report has the scan's locus" has "<td>${CHR}:${E1E}-$((E3S + 1))</td><td>3</td><td>18</td>" \
    "$(cat "$HTML" 2>/dev/null)"
else
  fail "nextflow is not on PATH"
fi

rm -rf "${GENOME_DIR:?}/${T}" "${GENOME_DIR:?}/${R}"
finish
