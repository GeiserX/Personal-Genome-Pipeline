#!/usr/bin/env bash
# Step 38's shared segments (Beagle against step 26b's panel, then hap-ibd) on
# the whole of chr20 of the GIAB Ashkenazi trio: HG002 and his parents HG003
# and HG004. The fixture's slices are far shorter than hap-ibd's 2 cM minimum,
# so the genotypes are GIAB's v4.2.1 benchmark calls of each person on chr20,
# read from GIAB's bucket, at the panel's SNVs inside every member's confident
# regions (0/0 where GIAB lists no variant there). They go in as the family's
# joint call, and step 38 runs with FAMILY_IBD_ONLY=true, the mode that reads
# an earlier joint call.
#   - setup.sh --ancestry-fine-panel installs the panel for chr20 only
#     (ANCESTRY_FINE_REGIONS) in a folder of its own (ANCESTRY_FINE_DIR);
#   - the output schema: the exact headers of ibd_pairs.tsv and
#     ibd_segments.tsv, one pairs row per pair, and the segments rows adding
#     up to each pair's count, total and longest stretch;
#   - each parent-child pair shares a stretch over 95% or more of chr20's
#     analysed genetic map, and is estimated first degree (parent-child);
#   - the two parents share no stretch longer than 5 cM;
#   - red control: HG005, an unrelated GIAB sample (the Chinese trio's son),
#     planted as the mother "HG004": the pair HG002 and HG004 shares under
#     5% of chr20, and the 95% check fails on it;
#   - FAMILY_IBD_ONLY=true without a joint call, or without the panel, stops
#     with the reason.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
PANEL="${G}/reference/ibd_chr20_panel"
REF_FAI="${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai"
GIAB=https://giab.s3.amazonaws.com/release
declare -A TRUTH=(
  [HG002]="${GIAB}/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark"
  [HG003]="${GIAB}/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark"
  [HG004]="${GIAB}/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark"
  [HG005]="${GIAB}/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark"
)
PAIRS_HEAD=$'sample_a\tsample_b\tdeclared\texpected_relatedness\tsegments\tshared_cm\tibd2_cm\tlongest_cm\tanalysed_cm\tshared_fraction\trelatedness\tdegree'
SEGS_HEAD=$'sample_a\tsample_b\tchrom\tstart\tend\tcm\tibd2_cm'

# pval DIR A B COLUMN: COLUMN of the ibd_pairs.tsv row of the pair A, B
pval() {
  awk -F'\t' -v a="$2" -v b="$3" -v k="$4" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    ($1 == a && $2 == b) || ($1 == b && $2 == a) {print $c[k]; exit}' "${1}/ibd_pairs.tsv" 2>/dev/null
}
ge() { awk -v x="${1:-}" -v y="$2" 'BEGIN {exit !(x ~ /^[0-9.]+$/ && x + 0 >= y)}'; }
lt() { awk -v x="${1:-}" -v y="$2" 'BEGIN {exit !(x ~ /^[0-9.]+$/ && x + 0 < y)}'; }
# parent_child DIR A B: the acceptance rule, 95% or more of the analysed map shared
parent_child() { ge "$(pval "$1" "$2" "$3" shared_fraction)" 0.95; }
not() { ! "$@"; }

# --- 1. The panel, chr20 only ------------------------------------------------------------
printf 'chr20\t0\t%s\n' "$(awk '$1 == "chr20" {print $2}' "$REF_FAI")" > "${CASE_TMP}/chr20.bed"
echo "+ ANCESTRY_FINE_DIR=${PANEL#"${G}"/} ANCESTRY_FINE_REGIONS=chr20 scripts/setup.sh --ancestry-fine-panel"
T0=$(date +%s)
ANCESTRY_FINE_DIR="$PANEL" ANCESTRY_FINE_REGIONS="${CASE_TMP}/chr20.bed" \
  "${REPO}/scripts/setup.sh" --ancestry-fine-panel "$G" > "${CASE_TMP}/setup.log" 2>&1
check_eq "setup.sh --ancestry-fine-panel (chr20) exits 0" "$?" 0
tail -n 5 "${CASE_TMP}/setup.log"
echo "ibd-segments: chr20 panel install $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
check "the panel's chr20 SNVs and map are there" test -s "${PANEL}/panel/chr20.vcf.gz" -a -s "${PANEL}/maps/chr20.map" -a -s "${PANEL}/panel/DONE"
check_ge "panel SNVs on chr20" "$(wc -l < "${PANEL}/panel/sites.tsv" | tr -d ' ')" 10000

# --- 2. GIAB's chr20 genotypes of the four, from GIAB's bucket -----------------------------
for s in HG002 HG003 HG004 HG005; do
  bed=${TRUTH[$s]}.bed
  [ "$s" = HG005 ] || bed=${TRUTH[$s]}_noinconsistent.bed
  curl -fsSL --retry 5 "$bed" | awk '$1 == "chr20"' > "${CASE_TMP}/${s}.bed"
  check_ge "${s}: confident regions on chr20" "$(wc -l < "${CASE_TMP}/${s}.bed" | tr -d ' ')" 1000
  docker run --rm -u "$(id -u):$(id -g)" -v "${CASE_TMP}:/w" -w /w "$BCFTOOLS_HTTPS_IMAGE" \
    bcftools view -H -r chr20 -v snps -m2 -M2 "${TRUTH[$s]}.vcf.gz" > "${CASE_TMP}/${s}.vcf" 2> "${CASE_TMP}/${s}.log"
  check_eq "${s}: chr20 SNVs read from GIAB" "$?" 0
  check_ge "${s}: chr20 SNVs" "$(wc -l < "${CASE_TMP}/${s}.vcf" | tr -d ' ')" 30000
done

# joint DIR NAME=SOURCE...: the members' genotypes at the panel's SNVs inside
# every member's confident regions, as DIR/family.vcf.gz
joint() {
  local dir=$1; shift
  mkdir -p "${G}/${dir}"
  python3 -I - "${PANEL}/panel/sites.tsv" "$REF_FAI" "$CASE_TMP" "$@" > "${G}/${dir}/family.vcf" <<'PY'
import bisect, sys
sites, fai, tmp, *members = sys.argv[1:]
names = [m.split("=")[0] for m in members]
srcs = [m.split("=")[1] for m in members]
def bed(path):
    iv = sorted((int(c[1]), int(c[2])) for c in (l.split() for l in open(path)) if len(c) >= 3)
    return [s for s, _ in iv], [e for _, e in iv]
def inside(b, p):
    i = bisect.bisect_right(b[0], p - 1) - 1
    return i >= 0 and b[1][i] >= p
beds = [bed(f"{tmp}/{s}.bed") for s in srcs]
gts = []
for s in srcs:
    g = {}
    for l in open(f"{tmp}/{s}.vcf"):
        f = l.rstrip("\n").split("\t")
        gt = f[9].split(":")[f[8].split(":").index("GT")].replace("|", "/")
        g[int(f[1])] = (f[3], f[4], {"1/0": "0/1"}.get(gt, gt))
    gts.append(g)
print("##fileformat=VCFv4.2")
for l in open(fai):
    c = l.split("\t")
    if c[0] == "chr20":
        print(f"##contig=<ID=chr20,length={c[1]}>")
print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names))
n = 0
for l in open(sites):
    c, p, r, a = l.rstrip("\n").split("\t")
    p = int(p)
    if c != "chr20" or not all(inside(b, p) for b in beds):
        continue
    row = []
    for g in gts:
        t = g.get(p)
        if t is None:
            row.append("0/0")
        elif t[0] == r and t[1] == a and t[2] in ("0/0", "0/1", "1/1"):
            row.append(t[2])
        else:
            break
    if len(row) == len(gts):
        print(f"chr20\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t" + "\t".join(row))
        n += 1
print(n, "sites", file=sys.stderr)
PY
  in_genome "$BCFTOOLS_IMAGE" bash -c "bcftools view -Oz -o ${dir}/family.vcf.gz ${dir}/family.vcf && bcftools index -t ${dir}/family.vcf.gz" \
    && rm -f "${G}/${dir}/family.vcf"
  printf 'sample\tsex\tfather\tmother\tpartner\nHG002\tmale\tHG003\tHG004\t.\nHG003\tmale\t.\t.\tHG004\nHG004\tfemale\t.\t.\tHG003\n' \
    > "${G}/${dir}/pedigree.tsv"
}
rm -rf "${G}/ibd_trio" "${G}/ibd_swap"
joint ibd_trio HG002=HG002 HG003=HG003 HG004=HG004
check_ge "the trio's joint call: panel SNVs on chr20" "$(vcf_count ibd_trio/family.vcf.gz)" 5000
# The red control: HG005's genotypes under the mother's name
joint ibd_swap HG002=HG002 HG003=HG003 HG004=HG005
check_ge "the planted trio's joint call: panel SNVs on chr20" "$(vcf_count ibd_swap/family.vcf.gz)" 5000

# --- 3. Step 38, shared segments only -------------------------------------------------------
export FAMILY_IBD_ONLY=true ANCESTRY_FINE_DIR="$PANEL" THREADS=4
T0=$(date +%s)
run_step 38-family.sh ibd_trio "${G}/ibd_trio/pedigree.tsv"
check_step_exit 38-family.sh
echo "ibd-segments: step 38 shared segments on chr20 of the trio $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
D="${G}/ibd_trio"
cat "${D}/ibd_pairs.tsv" "${D}/ibd_segments.tsv" 2>/dev/null
check_eq "ibd_pairs.tsv header" "$(head -n 1 "${D}/ibd_pairs.tsv" 2>/dev/null)" "$PAIRS_HEAD"
check_eq "ibd_segments.tsv header" "$(head -n 1 "${D}/ibd_segments.tsv" 2>/dev/null)" "$SEGS_HEAD"
check_eq "one pairs row per pair of the three" "$(awk 'NR > 1' "${D}/ibd_pairs.tsv" 2>/dev/null | cut -f1,2 | tr '\t\n' '- ')" "HG002-HG003 HG002-HG004 HG003-HG004 "
check_eq "phased and analysed: chr20" "$(cat "${D}/ibd/chroms.txt" 2>/dev/null)" chr20
check "the analysed map is most of chr20 (90 cM or more: $(pval "$D" HG002 HG003 analysed_cm))" ge "$(pval "$D" HG002 HG003 analysed_cm)" 90
# Each pair's segments rows add up to its pairs row
check "segments rows add up to each pair's count, total and longest" python3 -I - "${D}/ibd_pairs.tsv" "${D}/ibd_segments.tsv" <<'PY'
import sys
with open(sys.argv[1]) as f:
    h = f.readline().rstrip("\n").split("\t")
    P = [dict(zip(h, l.rstrip("\n").split("\t"))) for l in f]
with open(sys.argv[2]) as f:
    h = f.readline().rstrip("\n").split("\t")
    S = [dict(zip(h, l.rstrip("\n").split("\t"))) for l in f]
for p in P:
    s = [float(x["cm"]) for x in S if (x["sample_a"], x["sample_b"]) == (p["sample_a"], p["sample_b"])]
    assert len(s) == int(p["segments"]), (p, len(s))
    assert abs(sum(s) - float(p["shared_cm"])) <= 0.01 * len(s) + 0.01, (p, sum(s))
    assert abs(max(s or [0.0]) - float(p["longest_cm"])) <= 0.01, (p, max(s or [0.0]))
PY
for p in "HG002 HG003" "HG002 HG004"; do
  # shellcheck disable=SC2086  # the pair splits on purpose
  set -- $p
  check "$1 and $2: declared parent-child" test "$(pval "$D" "$1" "$2" declared)" = parent-child
  check "$1 and $2: share 95% or more of chr20's map ($(pval "$D" "$1" "$2" shared_fraction))" parent_child "$D" "$1" "$2"
  check_eq "$1 and $2: estimate" "$(pval "$D" "$1" "$2" degree)" "first degree (parent-child)"
done
check_eq "HG003 and HG004: declared partners" "$(pval "$D" HG003 HG004 declared)" partners
check "HG003 and HG004: no shared stretch longer than 5 cM (longest $(pval "$D" HG003 HG004 longest_cm))" \
  lt "$(pval "$D" HG003 HG004 longest_cm)" 5
check "the report has the section and the parent-child line" \
  has '\| HG002 and HG003 \| parent-child \| [0-9.]+ cM \| (9[5-9]|100)\.[0-9]% \|' "$(cat "${D}/ibd_report.md" 2>/dev/null)"
echo "ibd-segments: HG002-HG003 $(pval "$D" HG002 HG003 shared_fraction), HG002-HG004 $(pval "$D" HG002 HG004 shared_fraction), HG003-HG004 longest $(pval "$D" HG003 HG004 longest_cm) cM" >> "$E2E_NOTES"

# --- 4. Red control: HG005 planted as the mother ----------------------------------------
run_step 38-family.sh ibd_swap "${G}/ibd_swap/pedigree.tsv"
check_step_exit "38-family.sh (HG005 as the mother)"
X="${G}/ibd_swap"
cat "${X}/ibd_pairs.tsv" "${X}/ibd_segments.tsv" 2>/dev/null
check "HG002 and the planted HG004 share under 5% of chr20 ($(pval "$X" HG002 HG004 shared_fraction))" \
  lt "$(pval "$X" HG002 HG004 shared_fraction)" 0.05
check "the 95% parent-child check fails on that pair" not parent_child "$X" HG002 HG004
check "HG002 and HG003 still share 95% or more ($(pval "$X" HG002 HG003 shared_fraction))" parent_child "$X" HG002 HG003

# --- 5. What FAMILY_IBD_ONLY needs -------------------------------------------------------
mkdir -p "${G}/ibd_none"
cp "${G}/ibd_trio/pedigree.tsv" "${G}/ibd_none/"
run_step 38-family.sh ibd_none "${G}/ibd_none/pedigree.tsv"
check "no joint call: it stops (exit ${STEP_RC}) and names the file" \
  sh -c '[ "$1" -eq 1 ] && grep -q "FAMILY_IBD_ONLY=true reads .*family.vcf.gz" "$2"' _ "$STEP_RC" "$STEP_LOG"
ANCESTRY_FINE_DIR=none run_step 38-family.sh ibd_trio "${G}/ibd_trio/pedigree.tsv"
check "no panel: it stops (exit ${STEP_RC}) and says how to install it" \
  sh -c '[ "$1" -eq 1 ] && grep -q "setup.sh --ancestry-fine-panel" "$2"' _ "$STEP_RC" "$STEP_LOG"

rm -rf "${G}/ibd_trio" "${G}/ibd_swap" "${G}/ibd_none" "$PANEL"
finish
