#!/usr/bin/env bash
# Step 33's mosaic loss of chrY or chrX (bin/mosaic_sex_loss.py) on the fixture:
#   - HG002 (male), declared male, with the mosdepth windows of case 35: the
#     estimated loss of Y is near zero, inside its 95% range, not detected;
#   - the same BAM with 30% of its chrY reads removed (samtools view -s on
#     chrY only), step 16b run on it: the loss reads what was removed, detected;
#   - a synthetic woman, through the tool in PYTHON_IMAGE as the step runs
#     it: HG002's windows with chrX at autosomal depth and no chrY, and a
#     made-up VCF of heterozygous SNVs at 0.5 on chr20 and chrX, reads no
#     loss; with chrX depth at 0.85 and the chrX sites at 0.5 +- 0.088 (one X
#     lost in 30% of the cells) it reads about 0.30 from depth and from
#     allele fractions, detected.
# The fixture holds chrY:2.70-3.00 Mb, of which 2.78-3.00 Mb lies in the
# single-copy part the estimate reads (about 21 blocks of 10 kb), and
# chrX:73.7-74.0 Mb. Step 33 needs somalier's sites and VerifyBamID2's panel
# cut to the fixture reference's contigs, as case qc-1 does; this case runs
# before qc-1, so it installs its own copy in a folder of its own and leaves
# reference/somalier and reference/verifybamid2 to qc-1.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
WIN="${G}/${SAMPLE}/mosdepth/${SAMPLE}.regions.bed.gz"
tval() { awk -F'\t' -v k="$2" '$1 == k { print $2; exit }' "$1" 2>/dev/null; }
# within V LO HI: LO <= V <= HI, failing on an empty value.
within() { [ -n "$1" ] && awk -v v="$1" -v a="$2" -v b="$3" 'BEGIN { exit !(v + 0 >= a + 0 && v + 0 <= b + 0) }'; }

check "case 35 left HG002's mosdepth windows" test -s "$WIN"

# --- 1. somalier's sites and VerifyBamID2's panel, cut to the fixture's contigs ---------
SETUP="${G}/mosaic_setup"
"${REPO}/scripts/setup.sh" --sample-qc-data "$SETUP" > "${CASE_TMP}/setup.log" 2>&1
check_eq "setup.sh --sample-qc-data exits 0" "$?" 0
SITES="${G}/reference/mosaic_fixture/sites_fixture.vcf"
PANEL="${G}/reference/mosaic_fixture/1000g.phase3.100k.b38.vcf.gz.dat"
mkdir -p "$(dirname "$SITES")"
python3 - "${REF}.fai" "${SETUP}/reference/somalier/sites.hg38.vcf.gz" "$SITES" \
  "${SETUP}/reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat" "$PANEL" <<'PY'
import gzip, sys
fai, sites, fix_sites, panel, fix_panel = sys.argv[1:]
contigs = {l.split("\t")[0] for l in open(fai)}
with gzip.open(sites, "rt") as f, open(fix_sites, "w") as out:
    for l in f:
        if l.startswith("#") or l.split("\t", 1)[0] in contigs:
            out.write(l)
bed = open(panel + ".bed").read().splitlines()
keep = [i for i, l in enumerate(bed) if l.split("\t", 1)[0] in contigs]
for ext in ("bed", "mu", "UD"):
    lines = open(panel + "." + ext).read().splitlines()
    assert len(lines) == len(bed), ext
    with open(fix_panel + "." + ext, "w") as out:
        out.writelines(lines[i] + "\n" for i in keep)
PY
check_ge "somalier sites on the fixture's contigs" "$(grep -vc '^#' "$SITES" 2>/dev/null || true)" 9000
rm -rf "$SETUP"
export SOMALIER_SITES="$SITES" VERIFYBAMID2_PANEL="$PANEL"

# --- 2. HG002 as it is ----------------------------------------------------------------------
run_step 33-sample-qc.sh "$SAMPLE" male
check_step_exit 33-sample-qc.sh
T="${G}/${SAMPLE}/qc/${SAMPLE}_sample_qc.tsv"
grep '^mosaic_' "$T"
F0=$(tval "$T" mosaic_loss_depth) LO0=$(tval "$T" mosaic_loss_depth_low) HI0=$(tval "$T" mosaic_loss_depth_high)
check_eq "the estimate assumes the declared male (somalier cannot tell on the slices)" \
  "$(tval "$T" mosaic_sex)/$(tval "$T" mosaic_sex_from)/$(tval "$T" mosaic_chrom)" male/declared/chrY
check_ge "10 kb blocks on chrY's single-copy part" "$(tval "$T" mosaic_depth_blocks)" 15
check "HG002: the normalised chrY ratio is near 0.5 ($(tval "$T" mosaic_depth_ratio))" \
  within "$(tval "$T" mosaic_depth_ratio)" 0.45 0.55
check "HG002: the loss of Y is near zero (${F0})" within "$F0" -0.10 0.10
check "HG002: zero lies inside its 95% range (${LO0} to ${HI0})" within 0 "${LO0:-1}" "${HI0:--1}"
check_eq "HG002: not detected" "$(tval "$T" mosaic_call)" not_detected
check "the step prints the estimate" has 'Mosaic loss of chrY: not_detected' "$(cat "$STEP_LOG")"

# For the record, not a check: the raw ratio, mosdepth's whole-chrY mean over
# the autosomal mean of its summary, gaps and unsequenced slices included.
RAW=$(awk -F'\t' '$1 ~ /^chr[0-9]+$/ { b += $3; l += $2 } $1 == "chrY" { y = $4 }
  END { if (l && b) printf "%.3f", y / (b / l) }' "${G}/${SAMPLE}/mosdepth/${SAMPLE}.mosdepth.summary.txt")
echo "raw chrY ratio from mosdepth's summary: ${RAW:-?}; normalised: $(tval "$T" mosaic_depth_ratio)"

# --- 3. HG002 with 30% of its chrY reads removed -------------------------------------------
L="${SAMPLE}loy"
mkdir -p "${G}/${L}/aligned"
in_genome "$SAMTOOLS_IMAGE" sh -c "set -e
  samtools view -b -e 'rname != \"chrY\"' -o ${L}/aligned/rest.bam ${SAMPLE}/aligned/${SAMPLE}_sorted.bam
  samtools view -b -s 7.7 -o ${L}/aligned/y.bam ${SAMPLE}/aligned/${SAMPLE}_sorted.bam chrY
  samtools merge -f -o ${L}/aligned/${L}_sorted.bam ${L}/aligned/rest.bam ${L}/aligned/y.bam
  samtools index ${L}/aligned/${L}_sorted.bam
  rm ${L}/aligned/rest.bam ${L}/aligned/y.bam"
check "the BAM without 30% of chrY passes quickcheck" sam quickcheck "${L}/aligned/${L}_sorted.bam"
Y_ALL=$(sam view -c -F 0x904 "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" chrY 2>/dev/null)
Y_KEPT=$(sam view -c -F 0x904 "${L}/aligned/${L}_sorted.bam" chrY 2>/dev/null)
A_ALL=$(sam view -c -F 0x904 "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" chr20 2>/dev/null)
A_KEPT=$(sam view -c -F 0x904 "${L}/aligned/${L}_sorted.bam" chr20 2>/dev/null)
REMOVED=$(awk -v a="${Y_ALL:-0}" -v k="${Y_KEPT:-0}" 'BEGIN { if (a) printf "%.3f", 1 - k / a }')
echo "chrY reads: ${Y_ALL:-?} before, ${Y_KEPT:-?} after (${REMOVED} removed); chr20: ${A_ALL:-?} and ${A_KEPT:-?}"
check "about 30% of the chrY reads are gone (${REMOVED})" within "$REMOVED" 0.25 0.35
check_eq "no chr20 read is gone" "${A_KEPT:-x}" "${A_ALL:-y}"

run_step 16b-mosdepth.sh "$L"
check_step_exit 16b-mosdepth.sh
run_step 33-sample-qc.sh "$L" male
check_step_exit 33-sample-qc.sh
TL="${G}/${L}/qc/${L}_sample_qc.tsv"
grep '^mosaic_' "$TL"
F30=$(tval "$TL" mosaic_loss_depth)
WANT=$(awk -v r="$REMOVED" -v f="${F0:-0}" 'BEGIN { printf "%.3f", r + (1 - r) * f }')
check "30% of chrY gone: the loss reads ${WANT}, what was removed on top of HG002's own (${F30})" \
  within "$F30" "$(awk -v w="$WANT" 'BEGIN { print w - 0.05 }')" "$(awk -v w="$WANT" 'BEGIN { print w + 0.05 }')"
check "30% of chrY gone: the low end of its range is above 0.10 ($(tval "$TL" mosaic_loss_depth_low))" \
  within "$(tval "$TL" mosaic_loss_depth_low)" 0.10 1
check_eq "30% of chrY gone: detected" "$(tval "$TL" mosaic_call)" detected
check "the step says it is not a diagnosis" has 'A population-level marker, not a diagnosis' "$(cat "$STEP_LOG")"
echo "mosaic loss of Y: HG002 ${F0} (95% range ${LO0} to ${HI0}); with ${REMOVED} of chrY's reads removed ${F30}" >> "$E2E_NOTES"
in_genome "$BCFTOOLS_IMAGE" rm -rf "$L"

# --- 4. a synthetic woman --------------------------------------------------------------------
W="${SAMPLE}fem"
mkdir -p "${G}/${W}/mosdepth" "${G}/${W}/vcf" "${G}/${W}/qc"
# woman LOSS: windows and VCF of a woman who lost one X in LOSS of her cells.
woman() {
  python3 - "$WIN" "${G}/${W}" "$W" "$1" <<'PY'
import gzip, random, sys
win, out, name, loss = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
d = loss / (2 * (2 - loss))
rnd = random.Random(11)
with gzip.open(win, "rt") as f, gzip.open(f"{out}/mosdepth/{name}.regions.bed.gz", "wt") as o:
    for line in f:
        c, s, e, dp = line.rstrip("\n").split("\t")
        if c == "chrY":
            continue
        if c == "chrX":   # HG002's one X, as two, of which a share of the cells lost one
            dp = f"{float(dp) * 2 * (1 - loss / 2):.2f}"
        o.write(f"{c}\t{s}\t{e}\t{dp}\n")
lines = ["##fileformat=VCFv4.2\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n"]
def site(chrom, pos, p):
    n = rnd.randint(20, 40)
    a = sum(rnd.random() < p for _ in range(n))
    lines.append(f"{chrom}\t{pos}\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:{n - a},{a}\n")
for i in range(2000):
    site("chr20", 10_000_000 + 250 * i, 0.5)
for i in range(300):
    site("chrX", 73_700_000 + 1000 * i, 0.5 + (d if rnd.random() < 0.5 else -d))
with gzip.open(f"{out}/vcf/{name}.vcf.gz", "wt") as o:
    o.writelines(lines)
with open(f"{out}/qc/{name}_sample_qc.tsv", "w") as o:
    o.write(f"key\tvalue\nsample\t{name}\ninferred_sex\tunknown\ndeclared_sex\tfemale\n")
PY
  docker run --rm -u "$(id -u):$(id -g)" -v "${REPO}/bin:/pgp-bin:ro" -v "${G}:/genome" "$PYTHON_IMAGE" \
    python3 /pgp-bin/mosaic_sex_loss.py --table "/genome/${W}/qc/${W}_sample_qc.tsv" \
      --regions "/genome/${W}/mosdepth/${W}.regions.bed.gz" --vcf "/genome/${W}/vcf/${W}.vcf.gz"
}
TW="${G}/${W}/qc/${W}_sample_qc.tsv"
woman 0
FW0=$(tval "$TW" mosaic_loss_depth)
check_eq "woman, no loss: the estimate reads chrX for a declared female" \
  "$(tval "$TW" mosaic_sex)/$(tval "$TW" mosaic_chrom)" female/chrX
check "woman, no loss: the depth loss is near zero ($(tval "$TW" mosaic_loss_depth))" \
  within "$(tval "$TW" mosaic_loss_depth)" -0.10 0.10
check "woman, no loss: the allele-fraction loss is below 0.10 ($(tval "$TW" mosaic_loss_af))" \
  within "$(tval "$TW" mosaic_loss_af)" 0 0.0999
check_eq "woman, no loss: 300 chrX sites read" "$(tval "$TW" mosaic_x_het_sites)" 300
check_eq "woman, no loss: not detected" "$(tval "$TW" mosaic_call)" not_detected
woman 0.3
# The chrX ratio of the no-loss woman times 0.85: what HG002's own chrX reads, minus 15%.
WANTW=$(awk -v f="${FW0:-0}" 'BEGIN { printf "%.3f", 2 - 2 * (1 - f / 2) * 0.85 }')
check "woman, one X lost in 30%: the depth loss reads ${WANTW} ($(tval "$TW" mosaic_loss_depth))" \
  within "$(tval "$TW" mosaic_loss_depth)" "$(awk -v w="$WANTW" 'BEGIN { print w - 0.05 }')" \
  "$(awk -v w="$WANTW" 'BEGIN { print w + 0.05 }')"
check "woman, one X lost in 30%: the allele-fraction loss reads about 0.30 ($(tval "$TW" mosaic_loss_af))" \
  within "$(tval "$TW" mosaic_loss_af)" 0.20 0.40
check_eq "woman, one X lost in 30%: detected" "$(tval "$TW" mosaic_call)" detected
rm -rf "${G:?}/${W}" "${G}/reference/mosaic_fixture"

finish
