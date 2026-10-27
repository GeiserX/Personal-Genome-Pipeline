#!/usr/bin/env bash
# Step 51 (mosaic chromosomal alterations, MoChA) with the real images, and its
# Nextflow mirror:
#   - the pins: the MoChA image of versions.env prints MOCHA_VERSION, and
#     mocha.c at MOCHA_COMMIT defines that version (control: mocha.c at an
#     older commit defines another one, so the check can fail);
#   - without the panel the step stops and names setup.sh --mocha-data;
#   - setup.sh --mocha-data with MOCHA_REGIONS cut to the fixture's chr20 slice:
#     the 1000 Genomes chr20 file is checked against the release's md5 and cut,
#     the genetic map of chr20 unpacked, the folder under 200 MB;
#   - HG002 (case 21's VCF) phased by Beagle on that slice: phased hets, MoChA
#     makes no call and none is reported;
#   - MOSAIC, a synthetic phased VCF on chr16, chr19, chr20 and chr22 at 30x
#     (REF bases from the fixture reference) with a 20 Mb CN-LOH on 20p from
#     the telomere in 15% of the cells (the haplotype balance 0.575 instead of
#     0.5): exactly one call, CN-LOH, on chr20 over that stretch, with a cell
#     fraction of 0.10 to 0.20;
#   - MOSAIC0, the same sites, genotypes and depths with no shift: no call, so
#     the one-call check fails on it (red control);
#   - nextflow run main.nf --tools mocha,html_report on HG002 writes the
#     bash step's call table byte for byte and the same counts, and the HTML
#     report has the card. Control: MOSAIC's table differs.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
MD="${G}/reference/mocha-${MOCHA_PANEL_VERSION}"
summ() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "${G}/$1/mosaic/$1_mocha_summary.tsv" 2>/dev/null; }
reported() { awk -F'\t' 'NR > 1 && $13 == "yes"' "${G}/$1/mosaic/$1_mocha.tsv" 2>/dev/null; }
nrep() { reported "$1" | grep -c . || true; }
between() { awk -v x="${1:-}" -v lo="$2" -v hi="$3" 'BEGIN {exit !(x ~ /^[0-9.]+$/ && x + 0 >= lo && x + 0 <= hi)}'; }
differs() { [ -n "$1" ] && [ "$1" != "$2" ]; }
not_same() { ! cmp -s "$1" "$2"; }
status_ok() { [ "$1" = assessed ] || { [ "$1" = "not assessed" ] && [ -n "$2" ]; }; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && sleep $((i * 15)); done; return 1; }
show() { echo "--- $1"; cat "${G}/$1/mosaic/$1_mocha_summary.tsv" "${G}/$1/mosaic/$1_mocha.tsv" \
  "${G}/$1/mosaic/$1_mocha_calls_raw.tsv" "${G}/$1/mosaic/$1_mocha_stats.tsv" "${G}/$1/mosaic/$1_mocha_sites.tsv" 2>/dev/null; }

# --- the pins ---------------------------------------------------------------------------
check "MOCHA_IMAGE is pinned by digest" has '@sha256:[0-9a-f]{64}$' "$MOCHA_IMAGE"
USAGE=$(docker run --rm "$MOCHA_IMAGE" bcftools +mocha 2>&1 || true)
check "the image's plugin prints MOCHA_VERSION (${MOCHA_VERSION})" has "\(version ${MOCHA_VERSION}" "$USAGE"
define_at() { curl -fsSL "https://raw.githubusercontent.com/freeseek/mocha/$1/mocha.c" | sed -n 's/^#define MOCHA_VERSION "\(.*\)"$/\1/p'; }
retry define_at "$MOCHA_COMMIT" > "${T}/define.txt"
check_eq "mocha.c at MOCHA_COMMIT defines MOCHA_VERSION" "$(cat "${T}/define.txt")" "$MOCHA_VERSION"
retry define_at 3eacae195e6228877ab354b3bef6c6a7ad360700 > "${T}/define_old.txt"
check "control: mocha.c of 2024-09-27 defines another version ($(cat "${T}/define_old.txt"))" \
  differs "$(cat "${T}/define_old.txt")" "$MOCHA_VERSION"

# --- no panel ------------------------------------------------------------------------
MOCHA_DATA=reference/no-such-panel run_step 51-mocha.sh "$SAMPLE"
check_eq "without the panel the step exits 1" "$STEP_RC" 1
check "and names setup.sh --mocha-data" has 'setup\.sh --mocha-data' "$(cat "$STEP_LOG")"

# --- the panel, cut to the fixture's chr20 slice -------------------------------------
awk '$1 == "chr20"' "${FIXTURE_DIR}/regions.bed" > "${T}/regions.bed"
check_eq "the fixture has one chr20 slice" "$(grep -c . "${T}/regions.bed")" 1
echo "+ MOCHA_REGIONS=${T}/regions.bed scripts/setup.sh --mocha-data"
T0=$(date +%s)
MOCHA_REGIONS="${T}/regions.bed" "${REPO}/scripts/setup.sh" --mocha-data "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --mocha-data exits 0" "${PIPESTATUS[0]}" 0
T_SETUP=$(( $(date +%s) - T0 ))
check "the panel is finished (DONE)" test -s "${MD}/DONE"
check_ge "panel SNVs on the chr20 slice" "$(vcf_count "reference/mocha-${MOCHA_PANEL_VERSION}/panel/chr20.vcf.gz")" 1000
check_eq "people in the panel" "$(bcf query -l "reference/mocha-${MOCHA_PANEL_VERSION}/panel/chr20.vcf.gz" 2>/dev/null | grep -c .)" 3202
check_ge "lines of the chr20 genetic map" "$(wc -l < "${MD}/maps/chr20.map" 2>/dev/null | tr -d ' ')" 1000
check "the downloaded chromosome file is deleted" test ! -e "${MD}/raw"
MB=$(du -sm "$MD" | cut -f1)
check "the panel folder is under 200 MB (${MB} MB)" test "$MB" -lt 200
printf '#### Step 51 panel (chr20 slice)\n\n%s s to install, %s MB on disk, %s SNVs\n\n' "$T_SETUP" "$MB" \
  "$(vcf_count "reference/mocha-${MOCHA_PANEL_VERSION}/panel/chr20.vcf.gz")" >> "$E2E_NOTES"

# --- HG002, phased by Beagle ----------------------------------------------------------
check "case 21's VCF is there" test -s "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi"
run_step 51-mocha.sh "$SAMPLE"
check_step_exit 51-mocha.sh
show "$SAMPLE"
check_eq "HG002: phased by Beagle" "$(summ "$SAMPLE" phasing)" beagle
check_eq "HG002: MoChA read chr20, the panel's one chromosome" "$(summ "$SAMPLE" chromosomes)" chr20
check_ge "HG002: phased heterozygous sites" "$(summ "$SAMPLE" phased_het_sites)" 100
# 500 kb hold about 500 hets, too few for the genome-wide BAF autocorrelation
# MoChA's sample filter reads, so the status may be either state (with a
# reason); MoChA's own call count is what must be zero.
check "HG002: a status with its reason ($(summ "$SAMPLE" status), baf_auto $(summ "$SAMPLE" baf_auto))" \
  status_ok "$(summ "$SAMPLE" status)" "$(summ "$SAMPLE" message)"
check_eq "HG002 chr20: MoChA made no call" "$(summ "$SAMPLE" calls_mocha)" 0
check_eq "HG002 chr20: no mosaic call reported" "$(nrep "$SAMPLE")" 0
echo "mocha: HG002 chr20 slice: $(summ "$SAMPLE" het_sites) hets ($(summ "$SAMPLE" phased_het_sites) phased)," \
  "depth $(summ "$SAMPLE" depth)x, floor $(summ "$SAMPLE" floor), $(summ "$SAMPLE" calls_mocha) MoChA calls" >> "$E2E_NOTES"

# --- MOSAIC and MOSAIC0: a planted CN-LOH, and none ---------------------------------
for s in MOSAIC MOSAIC0; do rm -rf "${G:?}/${s}"; mkdir -p "${G}/${s}/vcf"; done
python3 -I - "${G}/reference/GRCh38_no_alt_analysis_set.fasta" "${G}/MOSAIC/vcf/MOSAIC.vcf" "${G}/MOSAIC0/vcf/MOSAIC0.vcf" <<'PY'
import math, random, sys
fa, out_shift, out_flat = sys.argv[1:]
CHROMS = ("chr16", "chr19", "chr20", "chr22")
EVENT = ("chr20", 1, 20_000_000)          # 20p from the telomere
BDEV = 0.075                              # CN-LOH in 15% of the cells: haplotype balance 0.575
fai = {}
for line in open(fa + ".fai"):
    c = line.split("\t")
    fai[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
ALT = {"A": "G", "G": "A", "C": "T", "T": "C"}
def poisson(lam):
    l, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= random.random()
        if p <= l:
            return k
        k += 1
def binom(n, f):
    return sum(random.random() < f for _ in range(n))
outs = {}
for path, shift in ((out_shift, True), (out_flat, False)):
    random.seed(51)
    lines = ["##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="All filters passed">']
    lines += [f"##contig=<ID={c},length={fai[c][0]}>" for c in CHROMS]
    lines += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
              '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + ("MOSAIC" if shift else "MOSAIC0")]
    with open(fa, "rb") as f:
        for c in CHROMS:
            length, off, lb, lw = fai[c]
            pos = 100_000
            while pos < length - 100_000:
                pos += random.randint(800, 1200)
                f.seek(off + (pos - 1) // lb * lw + (pos - 1) % lb)
                ref = f.read(1).decode().upper()
                p, r, depth = random.uniform(0.05, 0.5), random.random(), poisson(30)
                o = random.random() < 0.5
                bal = 0.5 + BDEV if shift and c == EVENT[0] and EVENT[1] <= pos <= EVENT[2] else 0.5
                k_het, k_err = binom(depth, bal), binom(depth, 0.005)
                if ref not in ALT:
                    continue
                if r < 2 * p * (1 - p):
                    gt = "0|1" if o else "1|0"
                    alt = depth - k_het if o else k_het
                elif r < 2 * p * (1 - p) + p * p:
                    gt, alt = "1|1", depth - k_err
                else:
                    continue
                lines.append(f"{c}\t{pos}\t.\t{ref}\t{ALT[ref]}\t50\tPASS\t.\tGT:GQ:DP:AD\t{gt}:50:{depth}:{depth - alt},{alt}")
    open(path, "w").write("\n".join(lines) + "\n")
    print(path, len(lines), "lines")
PY
for s in MOSAIC MOSAIC0; do
  in_genome "$BCFTOOLS_IMAGE" bash -c "bcftools view -Oz -o ${s}/vcf/${s}.vcf.gz ${s}/vcf/${s}.vcf && bcftools index -t ${s}/vcf/${s}.vcf.gz"
  check_ge "${s}: synthetic records" "$(vcf_count "${s}/vcf/${s}.vcf.gz")" 80000
  MOCHA_PHASED_VCF="${s}/vcf/${s}.vcf.gz" run_step 51-mocha.sh "$s"
  check_step_exit "51-mocha.sh ${s}"
  show "$s"
  check_eq "${s}: read as already phased" "$(summ "$s" phasing)" input
  check_eq "${s}: assessed" "$(summ "$s" status)" assessed
done
ROW=$(reported MOSAIC)
check_eq "MOSAIC: exactly one call reported" "$(nrep MOSAIC)" 1
check_eq "MOSAIC: on chr20" "$(cut -f1 <<< "$ROW")" chr20
check_eq "MOSAIC: typed CN-LOH" "$(cut -f6 <<< "$ROW")" CN-LOH
check "MOSAIC: it starts in the first Mb ($(cut -f2 <<< "$ROW"))" between "$(cut -f2 <<< "$ROW")" 0 1000000
check "MOSAIC: and ends at 19 to 21 Mb ($(cut -f3 <<< "$ROW"))" between "$(cut -f3 <<< "$ROW")" 19000000 21000000
check "MOSAIC: cell fraction 0.10 to 0.20 ($(cut -f7 <<< "$ROW"))" between "$(cut -f7 <<< "$ROW")" 0.10 0.20
check_eq "MOSAIC0 (no shift): no call reported" "$(nrep MOSAIC0)" 0
check "red control: the one-call check fails on MOSAIC0" test "$(nrep MOSAIC0)" != 1
check "red control: the no-MoChA-call check on HG002 fails on MOSAIC ($(summ MOSAIC calls_mocha))" differs "$(summ MOSAIC calls_mocha)" 0
echo "mocha: MOSAIC planted CN-LOH chr20:1-20000000 at 0.15: $(cut -f1-7 <<< "$ROW" | tr '\t' ' ');" \
  "floor $(summ MOSAIC floor), depth $(summ MOSAIC depth)x; MOSAIC0 $(summ MOSAIC0 calls_mocha) MoChA calls" >> "$E2E_NOTES"

# --- Nextflow against the bash step ------------------------------------------------------
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
SHEET="${T}/samplesheet.csv"
printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$SAMPLE" \
  "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" > "$SHEET"
NF_OUT="${G}/nf-mocha"
echo "+ nextflow run main.nf -profile docker --tools mocha,html_report"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-mocha" \
    --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools mocha,html_report \
    --mocha_data "$MD" --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
check_eq "nextflow run exits 0" "${PIPESTATUS[0]}" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null || true
# A failed task's own error, for diagnosis (report-only tasks fail without stopping the run).
while read -r x; do
  [ "$(cat "$x")" = 0 ] || { echo "--- failed task ${x%/.exitcode}"; tail -n 30 "${x%/.exitcode}/.command.err"; }
done < <(find "${E2E_WORK}/nf-work-mocha" -name .exitcode 2>/dev/null)
A="${G}/${SAMPLE}/mosaic/${SAMPLE}_mocha.tsv" B="${NF_OUT}/${SAMPLE}/mosaic/${SAMPLE}_mocha.tsv"
check "Nextflow wrote ${SAMPLE}_mocha.tsv" test -s "$B"
check "Nextflow's ${SAMPLE}_mocha.tsv is the bash step's" cmp -s "$A" "$B"
cmp -s "$A" "$B" || diff "$A" "$B" | head -20
for k in status chromosomes het_sites calls_mocha calls_reported floor_min phasing; do
  check_eq "Nextflow's summary: ${k} as the bash step's" \
    "$(awk -F'\t' -v k="$k" '$1 == k {print $2}' "${NF_OUT}/${SAMPLE}/mosaic/${SAMPLE}_mocha_summary.tsv" 2>/dev/null)" "$(summ "$SAMPLE" "$k")"
done
check "control: the comparison sees a difference" not_same "${G}/MOSAIC/mosaic/MOSAIC_mocha.tsv" "$B"
HTML=$(find "$NF_OUT" -name "${SAMPLE}_report.html" | head -n 1)
check "the HTML report has the MoChA card" has "<h2>Mosaic Chromosomal Alterations \(MoChA\)</h2>" "$(cat "$HTML" 2>/dev/null)"

rm -rf "${G:?}/MOSAIC0" "${G:?}/MOSAIC/vcf/MOSAIC.vcf" "${E2E_WORK:?}/nf-work-mocha"
finish
