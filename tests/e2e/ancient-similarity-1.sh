#!/usr/bin/env bash
# Step 44 (ancient populations) against a synthetic AADR-like panel, and the
# real AADR v66.p1 annotation:
#   - without the panel, step 44 says so in one line and exits 0;
#   - the AADR's .anno, .ind and .snp files of AADR_VERSION, downloaded from
#     Harvard Dataverse, have the md5s versions.env pins, and select chooses
#     800 or more present-day people in 35 or more groups and 10,000 or more
#     ancient people from them (the 7.1 GB .geno is not downloaded here);
#   - a synthetic panel in the AADR's own formats (the transposed packed
#     .geno, TGENO, that only EIGENSOFT 9 reads), on made-up GRCh37
#     coordinates 1,000 bp below the real GRCh38 ones with a chain that moves
#     them back: 3 present-day groups (Sardinian drawn from HG002's GIAB
#     genotypes, French and Finnish from random frequencies) and 3 ancient
#     groups, one of them pseudo-haploid copies of HG002's GIAB genotypes
#     (Planted_HG002), with missing sites. setup.sh --ancient-panel builds it,
#     every site lifts back to its GRCh38 position, and the install's own
#     axis check passes;
#   - step 44 on HG002's gVCF (case 21): most sites are used, Planted_HG002 is
#     the nearest group with its date and place, the plot is an SVG with the
#     sample marked, and both reports carry the section;
#   - red control: the same panel with the two alleles of every site swapped
#     (each genotype read as its opposite) must not put Planted_HG002 first;
#   - nextflow run --tools ancient,html_report on the same gVCF ranks the
#     groups as the bash step does, and its HTML report has the card.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
GVCF="${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
PANEL="${G}/reference/ancient/aadr_${AADR_VERSION}"
OUT="${G}/${SAMPLE}/ancient"
check "case 21 wrote a gVCF" vcf_ok "$GVCF"
rank_of() { awk -F'\t' -v g="$2" '$2 == g {print $1}' "$1" 2>/dev/null; }

# --- no panel: one line, exit 0 -----------------------------------------------------
run_step 44-ancient.sh "$SAMPLE"
check_step_exit "44-ancient.sh (no panel)"
check_eq "step 44 without a panel prints one line" "$(grep -c . "$STEP_LOG")" 1
check "that line says it was skipped and how to install the panel" has '^Step 44 skipped: .*setup\.sh --ancient-panel' "$(cat "$STEP_LOG")"

# --- the real AADR annotation -----------------------------------------------------------
REAL="${CASE_TMP}/aadr"
mkdir -p "$REAL"
for e in ANNO:anno IND:ind SNP:snp; do
  id_var="AADR_${e%%:*}_ID" md5_var="AADR_${e%%:*}_MD5"
  curl -fsSL --retry 3 -o "${REAL}/aadr.${e#*:}" "https://dataverse.harvard.edu/api/access/datafile/${!id_var}"
  check_eq "AADR ${AADR_VERSION} .${e#*:}: the md5 versions.env pins" "$(md5sum "${REAL}/aadr.${e#*:}" | cut -d' ' -f1)" "${!md5_var}"
done
python3 "${REPO}/bin/ancient_report.py" select --anno "${REAL}/aadr.anno" --ind "${REAL}/aadr.ind" \
  --out-ind "${REAL}/panel.ind" --poplist "${REAL}/poplist.txt" --ancients "${REAL}/ancients.tsv" | tee "${CASE_TMP}/select.log"
check_ge "real AADR: present-day reference people" "$(awk '$3 != "Ignore" && $3 != "Ancient"' "${REAL}/panel.ind" | wc -l | tr -d ' ')" 800
check_ge "real AADR: present-day groups" "$(wc -l < "${REAL}/poplist.txt" | tr -d ' ')" 35
check_ge "real AADR: ancient people projected" "$(awk '$3 == "Ancient"' "${REAL}/panel.ind" | wc -l | tr -d ' ')" 10000
check_eq "real AADR: the .ind and the .snp have the counts the 1240K .geno's header holds (23089 people, 1233013 sites)" \
  "$(wc -l < "${REAL}/aadr.ind" | tr -d ' ') $(wc -l < "${REAL}/aadr.snp" | tr -d ' ')" "23089 1233013"
printf '#### Step 44, real AADR %s selection\n\n%s\n\n' "$AADR_VERSION" "$(cat "${CASE_TMP}/select.log")" >> "$E2E_NOTES"

# --- the synthetic panel ------------------------------------------------------------------
SRC="${CASE_TMP}/src"
CHAIN="${G}/reference/liftover-e2e/synthetic_hg19ToHg38.chain"
mkdir -p "$SRC" "$(dirname "$CHAIN")"
python3 -I - "${FIXTURE_DIR}/regions.bed" "${FIXTURE_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "${FIXTURE_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" "$REF" "$SRC" "$CHAIN" <<'PY'
import bisect, gzip, math, random, struct, sys
regions, truth, bed, ref, out, chain = sys.argv[1:]
OFF = 1000
rnd = random.Random(44)
fai = {l.split("\t")[0]: [int(x) for x in l.split("\t")[1:5]] for l in open(ref + ".fai")}
fa = open(ref, "rb")
def base(c, p):
    n, off, lb, lw = fai[c]
    fa.seek(off + ((p - 1) // lb) * lw + (p - 1) % lb)
    return fa.read(1).decode().upper()
slices = [(c, int(s), int(e)) for c, s, e, *_ in (l.split("\t") for l in open(regions))
          if c in fai and c[3:].isdigit() and int(c[3:]) <= 22 and int(s) > OFF + 10]
conf = {}
for l in open(bed):
    c, s, e = l.split()[:3]
    conf.setdefault(c, []).append((int(s), int(e)))
for v in conf.values():
    v.sort()
def confident(c, p):
    v = conf.get(c, [])
    i = bisect.bisect_right([s for s, _ in v], p - 1) - 1
    return i >= 0 and v[i][1] >= p
truth_gt = {}
for l in gzip.open(truth, "rt"):
    if l.startswith("#"):
        continue
    f = l.split("\t")
    if len(f[3]) == 1 and len(f[4]) == 1:
        g = f[9].split(":")[0].replace("|", "/")
        truth_gt[(f[0], int(f[1]))] = (f[3], f[4], {"0/1": 1, "1/0": 1, "1/1": 2}.get(g))
sites = []  # (chrom, pos38, ref, alt, HG002 alt count)
for c, s, e in slices:
    for p in range(s + 1, e + 1):
        t = truth_gt.get((c, p))
        if t and t[2] is not None and confident(c, p) and base(c, p) == t[0]:
            sites.append((c, p, t[0], t[1], t[2]))
        elif p % 150 == 0 and (c, p) not in truth_gt and confident(c, p) and base(c, p) in "ACGT":
            r = base(c, p)
            sites.append((c, p, r, {"A": "G", "C": "T", "G": "A", "T": "C"}[r], 0))
sites.sort(key=lambda x: (int(x[0][3:]), x[1]))
assert len(sites) >= 1000, len(sites)
n = len(sites)
# People: (id, sex, group, kind, date_bp, country, lat, lon, dosages of the .snp's first allele (REF))
people = []
freq = {"French": [rnd.uniform(0.05, 0.95) for _ in sites], "Finnish": [rnd.uniform(0.05, 0.95) for _ in sites],
        "Sardinian": [{0: 0.1, 1: 0.5, 2: 0.9}[x[4]] for x in sites]}
for g, f in freq.items():
    for i in range(30):
        people.append((f"{g}{i}.DG", "M", g, "modern", 0, "Europe", 45, 5,
                       [2 - sum(rnd.random() < q for _ in range(2)) for q in f]))
def pseudo(alts):
    return [9 if rnd.random() < 0.3 else (0 if rnd.random() < a / 2 else 2) for a in alts]
for i in range(6):
    people.append((f"PL{i}.AG", "F", "Planted_HG002", "ancient", 5000, "Testland", 48, 10, pseudo([x[4] for x in sites])))
    people.append((f"FR{i}.AG", "M", "Anc_French_like", "ancient", 7000, "Franceland", 46, 2,
                   pseudo([2 * q for q in freq["French"]])))
    people.append((f"FI{i}.AG", "M", "Anc_Finnish_like", "ancient", 3500, "Finland", 62, 25,
                   pseudo([2 * q for q in freq["Finnish"]])))
people.append(("OUT0.AG", "M", "Far_Away", "ancient", 5000, "China", 40, 116, pseudo([1] * n)))  # outside West Eurasia
cols = ["Genetic ID", "Individual ID", "Date mean in BP in years before 1950 CE", "Date standard deviation in BP",
        "Full Date", "Group ID", "Locality", "Political Entity", "Latitude", "Longitude", "Data type",
        "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)", "ASSESSMENT", "ASSESSMENT WARNINGS"]
with open(f"{out}/aadr.anno", "w") as o:
    o.write("\t".join(cols) + "\n")
    for pid, sex, g, kind, bp, country, lat, lon, _ in people:
        o.write("\t".join([pid, pid.split(".")[0], str(bp), "50", "synthetic", g, "Somewhere", country, str(lat), str(lon),
                           "Shotgun.diploid" if kind == "modern" else "1240k", "1100000" if kind == "modern" else "500000",
                           "PASS", ""]) + "\n")
with open(f"{out}/aadr.ind", "w") as o:
    for pid, sex, g, *_ in people:
        o.write(f"{pid:>20} {sex} {g}\n")
with open(f"{out}/aadr.snp", "w") as o:
    for i, (c, p, r, a, _) in enumerate(sites):
        o.write(f"{'syn' + str(i):>20} {int(c[3:]):>2} 0.000000 {p - OFF:>12} {r} {a}\n")
# TGENO: a 48-byte header, then one record per person of ceil(n / 4) bytes,
# two bits per site (first site in the high bits), 3 for missing.
trlen = math.ceil(n * 2 / 8)
with open(f"{out}/aadr.geno", "wb") as o:
    o.write(f"TGENO {len(people):7d} {n:7d} 0 0".encode().ljust(48, b"\0"))
    for *_, d in people:
        rec = bytearray(b"\xff" * trlen)
        for k, x in enumerate(d):
            v = 3 if x == 9 else x
            sh = (3 - k % 4) * 2
            rec[k // 4] = (rec[k // 4] & ~(3 << sh) & 0xFF) | (v << sh)
        o.write(rec)
# The chain: each chromosome's slices, GRCh37 = GRCh38 - OFF.
with open(chain, "w") as o:
    k = 0
    for c in sorted({x[0] for x in sites}, key=lambda c: int(c[3:])):
        sz = fai[c][0]
        for _, s, e in [x for x in slices if x[0] == c]:
            k += 1
            o.write(f"chain 1000 {c} {sz} + {s - OFF} {e - OFF} {c} {sz} + {s} {e} {k}\n{e - s}\n\n")
with open(f"{out}/expected_sites.tsv", "w") as o:
    for i, (c, p, r, a, _) in enumerate(sites):
        o.write(f"syn{i}\t{c}\t{p}\t{r}\t{a}\n")
print(f"synthetic panel: {n} sites, {len(people)} people")
PY
check_ge "synthetic panel: sites in the fixture's slices" "$(wc -l < "${SRC}/aadr.snp" | tr -d ' ')" 1000
echo "+ ANCIENT_SOURCE_DIR=... ANCIENT_CHAIN=... scripts/setup.sh --ancient-panel"
T0=$(date +%s)
ANCIENT_SOURCE_DIR="$SRC" ANCIENT_CHAIN="$CHAIN" "${REPO}/scripts/setup.sh" --ancient-panel "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --ancient-panel exits 0" "${PIPESTATUS[0]}" 0
echo "setup.sh --ancient-panel (synthetic): $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
check "the panel is finished (DONE)" test -s "${PANEL}/DONE"
check "its two smartpca runs gave the same axes" has 'Axes agree with the panel' "$(cat "$STEP_LOG")"
check "every site lifted back to its GRCh38 position, alleles in the .snp's order" \
  cmp -s <(tail -n +2 "${PANEL}/sites_grch38.tsv") "${SRC}/expected_sites.tsv"
check_eq "present-day people (3 groups of 30)" "$(wc -l < "${PANEL}/moderns.ind" 2>/dev/null | tr -d ' ')" 90
check_eq "ancient people projected (Far_Away is outside West Eurasia)" "$(awk 'NR > 1' "${PANEL}/ancients.tsv" 2>/dev/null | wc -l | tr -d ' ')" 18
check_eq "panel.evec: the ancient people projected" "$(awk '$NF == "Ancient"' "${PANEL}/panel.evec" 2>/dev/null | wc -l | tr -d ' ')" 18
check "the AADR files were deleted once the panel was built" test ! -e "${G}/reference/ancient/raw_${AADR_VERSION}"

# --- step 44 on HG002 ------------------------------------------------------------------------
run_step 44-ancient.sh "$SAMPLE"
check_step_exit "44-ancient.sh ${SAMPLE}"
cat "${OUT}/${SAMPLE}_ancient.txt" 2>/dev/null
TSV="${OUT}/${SAMPLE}_ancient.tsv" INFO="${OUT}/${SAMPLE}_ancient_info.tsv"
info() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "$INFO" 2>/dev/null; }
USED=$(info sites_used) ALL=$(info panel_sites)
check_ge "sites used (of ${ALL:-?})" "$USED" "$(( ${ALL:-0} * 8 / 10 ))"
check_eq "the nearest ancient group is Planted_HG002" "$(rank_of "$TSV" Planted_HG002)" 1
check_eq "its row: 6 people, 3050 BCE, Testland" \
  "$(awk -F'\t' '$2 == "Planted_HG002" {print $3, $8, $9}' "$TSV" 2>/dev/null)" "6 3050 BCE Testland"
check_eq "three ancient groups ranked" "$(awk 'NR > 1' "$TSV" 2>/dev/null | wc -l | tr -d ' ')" 3
SVG="${OUT}/${SAMPLE}_ancient_pca.svg"
check "the plot is an SVG with the sample's label" has "^<svg .*>${SAMPLE}</text>" "$(tr -d '\n' < "$SVG" 2>/dev/null)"
check "the text carries the caveat" has 'not descent' "$(cat "${OUT}/${SAMPLE}_ancient.txt" 2>/dev/null)"
python3 "${REPO}/bin/render_report.py" --sample "$SAMPLE" --sample-dir "${G}/${SAMPLE}" --json "${CASE_TMP}/summary.json" \
  -o "${CASE_TMP}/report.txt" -o "${CASE_TMP}/report.html" > /dev/null
check "text report: the section with Planted_HG002 first" has '^ +1\. Planted_HG002 \(3050 BCE, Testland; 6 people\)' "$(cat "${CASE_TMP}/report.txt")"
check "HTML report: the card and the plot" has '<h2>Ancient Populations</h2>.*data:image/svg\+xml;base64,' "$(tr -d '\n' < "${CASE_TMP}/report.html")"

# --- red control: every site's alleles swapped ------------------------------------------------
RED="${G}/reference/ancient/red_swapped"
rm -rf "$RED"
cp -r "$PANEL" "$RED"
awk -F'\t' -v OFS='\t' 'NR > 1 {t = $4; $4 = $5; $5 = t} {print}' "${PANEL}/sites_grch38.tsv" > "${RED}/sites_grch38.tsv"
check "red: the copy differs in every site" test "$(diff "${PANEL}/sites_grch38.tsv" "${RED}/sites_grch38.tsv" | grep -c '^>')" -eq "$(( $(wc -l < "${PANEL}/sites_grch38.tsv") - 1 ))"
cp -r "$OUT" "${CASE_TMP}/bash_ancient"
ANCIENT_PANEL_DIR="$RED" run_step 44-ancient.sh "$SAMPLE"
check_step_exit "44-ancient.sh ${SAMPLE} on the swapped panel"
RED_RANK=$(rank_of "$TSV" Planted_HG002)
check "red: with every genotype read as its opposite, Planted_HG002 is not first (rank ${RED_RANK:-none})" test "${RED_RANK:-0}" -ne 1
rm -rf "$OUT" "$RED"
cp -r "${CASE_TMP}/bash_ancient" "$OUT"

# --- the Nextflow mirror --------------------------------------------------------------------
if command -v nextflow >/dev/null; then
  NF_OUT="${G}/nf-ancient"
  SHEET="${CASE_TMP}/samplesheet.csv"
  printf 'sample,vcf,vcf_index,gvcf,gvcf_index,sex\n%s,%s,%s,%s,%s,male\n' "$SAMPLE" \
    "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${G}/${GVCF}" "${G}/${GVCF}.tbi" > "$SHEET"
  (
    cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-ancient" \
      --input "$SHEET" --reference "$REF" --tools ancient,html_report --ancient_panel "$PANEL" \
      --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB
  ) 2>&1 | tee "$STEP_LOG"
  check_eq "nextflow run --tools ancient,html_report exits 0" "${PIPESTATUS[0]}" 0
  NF_TSV="${NF_OUT}/${SAMPLE}/ancient/${SAMPLE}_ancient.tsv"
  for f in "${SAMPLE}_ancient.tsv" "${SAMPLE}_ancient_info.tsv" "${SAMPLE}_ancient_pca.svg" "${SAMPLE}_ancient.txt"; do
    check "Nextflow published ancient/${f}" test -s "${NF_OUT}/${SAMPLE}/ancient/${f}"
  done
  check_eq "Nextflow ranks the groups as the bash step does" "$(cut -f1-3 "$NF_TSV" 2>/dev/null | tr '\n' ' ')" "$(cut -f1-3 "$TSV" | tr '\n' ' ')"
  check_eq "and uses the same sites" "$(awk -F'\t' '$1 == "sites_used" {print $2}' "${NF_OUT}/${SAMPLE}/ancient/${SAMPLE}_ancient_info.tsv" 2>/dev/null)" "$USED"
  check "Nextflow HTML report: the card" has '<h2>Ancient Populations</h2>' "$(cat "${NF_OUT}/${SAMPLE}/${SAMPLE}_report.html" 2>/dev/null)"
else
  fail "nextflow is not on PATH"
fi

finish
