#!/usr/bin/env bash
# Step 26b (fine-grained ancestry) against a reduced panel, read from gnomAD's
# bucket by setup.sh --ancestry-fine-panel for the fixture's slices plus
# chr22:20-35 Mb (the slices are too short to paint: RFMix needs 5 cM):
#   - without the panel, step 26b says so in one line and exits 0;
#   - the panel installs: both ADMIXTURE levels, painting references, and
#     neither person of the synthetic mix in it;
#   - HG002 from its gVCF on the slices (cases 21 and 03): the continental
#     fractions sum to 1, West Eurasia comes first with 80% or more, and with
#     no chromosome long enough the step says it painted nothing;
#   - HG002T, HG002's GIAB v4.2.1 genotypes at the panel's SNVs inside GIAB's
#     confident regions (0/0 where GIAB lists no variant): West Eurasia 85% or
#     more, and 80% or more of the painted chr22 painted West Eurasia;
#   - MIX, a planted person: chr22 below 27.5 Mb and the first half of each
#     slice from a Yoruba (NA18486), the rest from a Han Chinese (NA18525),
#     both left out of the panel: Africa and East Asia each 35% to 65%, and
#     the painting switches at the junction (Africa left of 25.5 Mb, East Asia
#     right of 29.5 Mb, 80% or more of each side);
#   - red control: the same panel with every reference person's labels moved
#     to the next group (Africa's people called West Eurasia, and so on) and
#     the ADMIXTURE frequencies fitted again: the HG002T checks must fail on
#     it (sample HG002R, the same genotypes);
#   - the Nextflow module on HG002T (with html_report) writes the fractions of
#     the bash step, and the HTML report has its card.
# The time and disk of the panel install go to the job summary.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
G="$GENOME_DIR"
PANEL="${G}/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}"
A=NA18486 B=NA18525
LONG_CHR=chr22 LONG_START=20000000 LONG_END=35000000 JUNCTION=27500000

frac() {  # frac SAMPLE LEVEL CODE: the fraction column of that row
  awk -F'\t' -v l="$2" -v c="$3" '$1 == l && $2 == c {print $4}' "${G}/$1/ancestry_fine/$1_ancestry_fine.tsv" 2>/dev/null
}
ge() { awk -v x="${1:-}" -v y="$2" 'BEGIN {exit !(x ~ /^[0-9.]+$/ && x + 0 >= y)}'; }
between() { awk -v x="${1:-}" -v lo="$2" -v hi="$3" 'BEGIN {exit !(x ~ /^[0-9.]+$/ && x + 0 >= lo && x + 0 <= hi)}'; }
# side SAMPLE CHROM FROM TO CODE: share of the painted bp of CHROM between FROM and TO painted CODE
side() {
  awk -F'\t' -v c="$2" -v a="$3" -v b="$4" -v g="$5" '
    $1 == c { s = ($2 > a ? $2 : a); e = ($3 < b ? $3 : b); if (e > s) { t += e - s; if ($5 == g) h += e - s } }
    END { if (t > 0) printf "%.3f\n", h / t }' "${G}/$1/ancestry_fine/$1_painting.bed" 2>/dev/null
}

# --- no panel: one line, exit 0 ---------------------------------------------------
ANCESTRY_FINE_DIR=none run_step 26b-ancestry-fine.sh "$SAMPLE"
check_step_exit "26b-ancestry-fine.sh (no panel)"
check_eq "step 26b without a panel prints one line" "$(grep -c . "$STEP_LOG")" 1
check "that line says it was skipped and how to install the panel" has '^Step 26b skipped: .*setup\.sh --ancestry-fine-panel' "$(cat "$STEP_LOG")"

# --- the reduced panel ----------------------------------------------------------------
REGIONS="${CASE_TMP}/regions.bed"
{ cat "${FIXTURE_DIR}/regions.bed"; printf '%s\t%s\t%s\n' "$LONG_CHR" "$LONG_START" "$LONG_END"; } > "$REGIONS"
echo "+ ANCESTRY_FINE_REGIONS=... ANCESTRY_FINE_EXCLUDE=${A},${B} scripts/setup.sh --ancestry-fine-panel"
T0=$(date +%s)
ANCESTRY_FINE_REGIONS="$REGIONS" ANCESTRY_FINE_EXCLUDE="${A},${B}" \
  "${REPO}/scripts/setup.sh" --ancestry-fine-panel "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --ancestry-fine-panel exits 0" "${PIPESTATUS[0]}" 0
T_SETUP=$(( $(date +%s) - T0 ))
check "the panel is finished (model/DONE)" test -s "${PANEL}/model/DONE"
check_eq "continental level: 6 groups" "$(awk -F'\t' '$1 == "continent" {print $2}' "${PANEL}/model/levels.tsv" 2>/dev/null)" 6
check_eq "regional level: 18 groups" "$(awk -F'\t' '$1 == "region" {print $2}' "${PANEL}/model/levels.tsv" 2>/dev/null)" 18
check_ge "SNVs for ADMIXTURE" "$(wc -l < "${PANEL}/model/ref.bim" 2>/dev/null | tr -d ' ')" 500
check_ge "painting references" "$(wc -l < "${PANEL}/model/paint_map.tsv" 2>/dev/null | tr -d ' ')" 150
check "the two people of the mix are not in the panel" \
  awk -F'\t' -v a="$A" -v b="$B" '$1 == a || $1 == b {bad = 1} END {exit bad}' "${PANEL}/samples.tsv"
printf '#### Step 26b reduced panel (fixture slices plus %s:%s-%s)\n\n%s s to install, %s on disk, %s SNVs, %s for ADMIXTURE\n\n' \
  "$LONG_CHR" "$LONG_START" "$LONG_END" "$T_SETUP" "$(du -sh "$PANEL" | cut -f1)" \
  "$(wc -l < "${PANEL}/panel/sites.tsv" | tr -d ' ')" "$(wc -l < "${PANEL}/model/ref.bim" | tr -d ' ')" >> "$E2E_NOTES"

# --- HG002 from its gVCF on the fixture's slices ---------------------------------------
T0=$(date +%s)
run_step 26b-ancestry-fine.sh "$SAMPLE"
check_step_exit "26b-ancestry-fine.sh ${SAMPLE}"
echo "step 26b on ${SAMPLE}: $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
cat "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.tsv" 2>/dev/null
check "it read the gVCF" has 'Input: the gvcf' "$(cat "$STEP_LOG")"
SUM=$(awk -F'\t' '$1 == "continent" {s += $4} END {printf "%.3f", s}' "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.tsv" 2>/dev/null)
check "the continental fractions sum to 1 (${SUM})" between "$SUM" 0.99 1.01
check_eq "West Eurasia comes first" "$(awk -F'\t' '$1 == "continent" {print $2; exit}' "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.tsv" 2>/dev/null)" WEA
check "with 80% or more ($(frac "$SAMPLE" continent WEA))" ge "$(frac "$SAMPLE" continent WEA)" 0.8
check "each continental fraction has a bootstrap range" \
  awk -F'\t' '$1 == "continent" && ($6 == "." || $7 == "." || $6 > $4 + 1e-9 || $7 < $4 - 1e-9) {bad = 1} END {exit bad}' \
  "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.tsv"
check "no slice spans 5 cM, so nothing was painted, and the summary says so" \
  has '^Chromosome painting: none' "$(cat "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.txt" 2>/dev/null)"
check "the PNG is a PNG" has '^89504e47' "$(head -c 4 "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_painting.png" 2>/dev/null | od -An -tx1 | tr -d ' \n')"

# --- HG002T: GIAB's genotypes at the panel's SNVs ------------------------------------
# Inside GIAB's confident regions a site GIAB lists no variant at is 0/0.
mkdir -p "${G}/HG002T/vcf"
docker run --rm -u "$(id -u):$(id -g)" -v "${G}:/genome:ro" -v "${FIXTURE_DIR}:/fx:ro" "$BCFTOOLS_IMAGE" \
  bcftools view -R "/genome/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/panel/positions.tsv" -Ov \
  /fx/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > "${CASE_TMP}/truth_at_sites.vcf"
python3 -I - "${PANEL}/panel/sites.tsv" "${CASE_TMP}/truth_at_sites.vcf" \
    "${FIXTURE_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" "${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai" \
    > "${G}/HG002T/vcf/HG002T.vcf" <<'PY'
import bisect, sys
sites, truth, bed, fai = sys.argv[1:]
iv = {}
for line in open(bed):
    c = line.split()
    if len(c) >= 3:
        iv.setdefault(c[0], []).append((int(c[1]), int(c[2])))
starts = {k: [s for s, _ in sorted(v)] for k, v in iv.items()}
ends = {k: [e for _, e in sorted(v)] for k, v in iv.items()}
def confident(c, p):
    i = bisect.bisect_right(starts.get(c, []), p - 1) - 1
    return i >= 0 and ends[c][i] >= p
gt = {}
for line in open(truth):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    g = f[9].split(":")[f[8].split(":").index("GT")].replace("|", "/")
    gt[(f[0], int(f[1]))] = (f[3], f[4], {"1/0": "0/1"}.get(g, g))
chroms = sorted({l.split("\t")[0] for l in open(sites)})
print("##fileformat=VCFv4.2")
print('##FILTER=<ID=PASS,Description="All filters passed">')
for l in open(fai):
    c = l.split("\t")
    if c[0] in chroms:
        print(f"##contig=<ID={c[0]},length={c[1]}>")
print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002T")
for l in open(sites):
    c, p, r, a = l.rstrip("\n").split("\t")
    p = int(p)
    if not confident(c, p):
        continue
    t = gt.get((c, p))
    if t is None:
        g = "0/0"
    elif t[0] == r and t[1] == a and t[2] in ("0/1", "1/1", "0/0"):
        g = t[2]
    else:
        continue
    print(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{g}")
PY
in_genome "$BCFTOOLS_IMAGE" bash -c 'bcftools view -Oz -o HG002T/vcf/HG002T.vcf.gz HG002T/vcf/HG002T.vcf && bcftools index -t HG002T/vcf/HG002T.vcf.gz'
check_ge "HG002T: GIAB genotypes at panel SNVs" "$(vcf_count HG002T/vcf/HG002T.vcf.gz)" 3000
T0=$(date +%s)
run_step 26b-ancestry-fine.sh HG002T
check_step_exit "26b-ancestry-fine.sh HG002T"
echo "step 26b on HG002T (painting ${LONG_CHR}): $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
cat "${G}/HG002T/ancestry_fine/HG002T_ancestry_fine.txt" 2>/dev/null
WEA_T=$(frac HG002T continent WEA) PAINT_T=$(side HG002T "$LONG_CHR" "$LONG_START" "$LONG_END" WEA)
check "HG002T: West Eurasia 85% or more (${WEA_T})" ge "$WEA_T" 0.85
check "HG002T: 80% or more of ${LONG_CHR} painted West Eurasia (${PAINT_T:-none})" ge "$PAINT_T" 0.8
check_ge "HG002T: painted Mb of ${LONG_CHR} on two haplotypes" \
  "$(awk -F'\t' -v c="$LONG_CHR" '$1 == c {t += $3 - $2} END {print int(t / 1e6)}' "${G}/HG002T/ancestry_fine/HG002T_painting.bed" 2>/dev/null)" 20

# --- MIX: half a Yoruba, half a Han Chinese ---------------------------------------------
# Their genotypes at the panel's SNVs, from the same bucket (the panel holds
# neither), unphased: the step phases them itself.
mkdir -p "${G}/MIX/vcf" "${G}/MIX/tmp"
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
docker run --rm -u "$(id -u):$(id -g)" -v "${G}:/genome" -w /genome/MIX/tmp "$BCFTOOLS_HTTPS_IMAGE" bash -c '
  for c in $(cut -f1 "$1" | uniq); do
    r=$(awk -v c="$c" '"'"'$1 == c {printf "%s%s:%d-%d", s, $1, $2 + 1, $3; s = ","}'"'"' "$3")
    bcftools view -s "$2" -T "$1" -r "$r" -m2 -M2 -v snps \
      "https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_${c}.filtered.SNV_INDEL.phased.shapeit5.bcf" 2>/dev/null \
      | bcftools query -H -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n"
  done; rm -f ./*.csi' _ "/genome/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/panel/positions.tsv" "${A},${B}" \
  "/genome/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/regions.bed" > "${CASE_TMP}/ab.tsv"
check_ge "genotypes of ${A} and ${B} at the panel's SNVs" "$(grep -vc '^#' "${CASE_TMP}/ab.tsv" || true)" 3000
# The first half of each region (of chr22:20-35 Mb, below 27.5 Mb) from A,
# the rest from B; the columns are found by name in bcftools' -H header.
{ printf '##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n'
  cut -f1 "${PANEL}/panel/sites.tsv" | uniq | awk -F'\t' 'NR == FNR {u[$1] = 1; next} ($1 in u) {printf "##contig=<ID=%s,length=%s>\n", $1, $2}' - "${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai"
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMIX\n'
  awk -F'\t' -v OFS='\t' -v regions="$REGIONS" -v a="$A" -v b="$B" '
    BEGIN { while ((getline l < regions) > 0) { split(l, f, "\t"); n++; rc[n] = f[1]; rs[n] = f[2]; re[n] = f[3] } }
    /^#/ { for (i = 5; i <= NF; i++) { if (index($i, "]" a ":")) ca = i; if (index($i, "]" b ":")) cb = i }; next }
    {
      mid = 0
      for (i = 1; i <= n; i++) if (rc[i] == $1 && $2 > rs[i] && $2 <= re[i]) mid = (rs[i] + re[i]) / 2
      g = ($2 < mid) ? $ca : $cb; sub(/\|/, "/", g); if (g == "1/0") g = "0/1"
      print $1, $2, ".", $3, $4, ".", "PASS", ".", "GT", g
    }' "${CASE_TMP}/ab.tsv"; } > "${G}/MIX/vcf/MIX.vcf"
in_genome "$BCFTOOLS_IMAGE" bash -c 'bcftools sort -Oz -o MIX/vcf/MIX.vcf.gz MIX/vcf/MIX.vcf && bcftools index -t MIX/vcf/MIX.vcf.gz'
T0=$(date +%s)
run_step 26b-ancestry-fine.sh MIX
check_step_exit "26b-ancestry-fine.sh MIX"
echo "step 26b on MIX: $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
cat "${G}/MIX/ancestry_fine/MIX_ancestry_fine.txt" 2>/dev/null
check "MIX: Africa 35% to 65% ($(frac MIX continent AFR))" between "$(frac MIX continent AFR)" 0.35 0.65
check "MIX: East Asia 35% to 65% ($(frac MIX continent EAS))" between "$(frac MIX continent EAS)" 0.35 0.65
L=$(side MIX "$LONG_CHR" "$LONG_START" $((JUNCTION - 2000000)) AFR)
R=$(side MIX "$LONG_CHR" $((JUNCTION + 2000000)) "$LONG_END" EAS)
check "MIX: Africa on 80% or more of ${LONG_CHR} left of the junction (${L:-none})" ge "$L" 0.8
check "MIX: East Asia on 80% or more of ${LONG_CHR} right of the junction (${R:-none})" ge "$R" 0.8

# --- red control: the reference labels moved to the next group ------------------------
RED="${G}/reference/ancestry_fine/red-control"
rm -rf "$RED"
cp -a "$PANEL" "$RED"
rm -rf "${RED}/model"
awk -F'\t' -v OFS='\t' '
  BEGIN { split("AFR WEA CSA EAS AMR OCE", c, " "); for (i = 1; i <= 6; i++) nc[c[i]] = c[i % 6 + 1]
          split("WAF EAF CAF SAF NWE NEE IBE ITA SAR CAU LEV BDN NAF CSA EAS SIB AMR OCE", r, " "); for (i = 1; i <= 18; i++) nr[r[i]] = r[i % 18 + 1] }
  NR == 1 { print; next }
  { if ($3 in nc) $3 = nc[$3]; if ($4 in nr) $4 = nr[$4]; print }' "${PANEL}/samples.tsv" > "${CASE_TMP}/rotated.tsv"
rm -f "${RED}/samples.tsv"
ANCESTRY_FINE_DIR="$RED" ANCESTRY_FINE_REGIONS="$REGIONS" ANCESTRY_FINE_SAMPLE_MAP="${CASE_TMP}/rotated.tsv" \
  "${REPO}/scripts/setup.sh" --ancestry-fine-panel "$G" > "${CASE_TMP}/red_setup.log" 2>&1
RED_RC=$?
check_eq "setup.sh refits the model on the rotated labels" "$RED_RC" 0
check "and reuses the panel's genotypes (no chromosome read again)" lacks 'reading the panel' "$(cat "${CASE_TMP}/red_setup.log")"
mkdir -p "${G}/HG002R/vcf"
cp "${G}/HG002T/vcf/HG002T.vcf.gz" "${G}/HG002R/vcf/HG002R.vcf.gz"
cp "${G}/HG002T/vcf/HG002T.vcf.gz.tbi" "${G}/HG002R/vcf/HG002R.vcf.gz.tbi"
ANCESTRY_FINE_DIR="$RED" run_step 26b-ancestry-fine.sh HG002R
check_step_exit "26b-ancestry-fine.sh HG002R (HG002T against the rotated labels)"
WEA_R=$(frac HG002R continent WEA) PAINT_R=$(side HG002R "$LONG_CHR" "$LONG_START" "$LONG_END" WEA)
check "red: with rotated labels the HG002T check fails (West Eurasia ${WEA_R:-none}, below 0.85)" test "$(ge "$WEA_R" 0.85 && echo pass)" != pass
check "red: and so does the painting check (${PAINT_R:-none}, below 0.8)" test "$(ge "$PAINT_R" 0.8 && echo pass)" != pass
rm -rf "$RED"

# --- the Nextflow module ------------------------------------------------------------------
NF_OUT="${G}/nf-ancestry-fine"
SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,vcf,vcf_index\nHG002T,%s,%s\n' "${G}/HG002T/vcf/HG002T.vcf.gz" "${G}/HG002T/vcf/HG002T.vcf.gz.tbi" > "$SHEET"
T0=$(date +%s)
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "${E2E_WORK}/nf-work-ancestry-fine" \
    --input "$SHEET" \
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools ancestry_fine,html_report \
    --ancestry_fine_panel "$PANEL" \
    --outdir "$NF_OUT" \
    --max_cpus 4 --max_memory 14.GB
) 2>&1 | tee "$STEP_LOG"
check_eq "nextflow run with --tools ancestry_fine,html_report exits 0" "${PIPESTATUS[0]}" 0
echo "Nextflow ancestry_fine on HG002T: $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
NF_TSV="${NF_OUT}/HG002T/ancestry_fine/HG002T_ancestry_fine.tsv"
check "the module published the fractions table" test -s "$NF_TSV"
check "its fractions are the bash step's (within 0.01 each)" python3 -I - "$NF_TSV" "${G}/HG002T/ancestry_fine/HG002T_ancestry_fine.tsv" <<'PY'
import sys
rows = [{(c[0], c[1]): float(c[3]) for c in (l.rstrip("\n").split("\t") for l in open(p)) if c[0] != "level"} for p in sys.argv[1:]]
diff = {k: (rows[0].get(k), rows[1].get(k)) for k in set(rows[0]) | set(rows[1])
        if rows[0].get(k) is None or rows[1].get(k) is None or abs(rows[0][k] - rows[1][k]) > 0.01}
print("differences:", diff or "none")
sys.exit(1 if diff else 0)
PY
check "the HTML report has the fine-grained ancestry card with its painting" \
  python3 -I -c 'import sys; t = open(sys.argv[1]).read(); sys.exit(0 if "Fine-Grained Ancestry" in t and "data:image/png;base64," in t else 1)' \
  "${NF_OUT}/HG002T/HG002T_report.html"

finish
