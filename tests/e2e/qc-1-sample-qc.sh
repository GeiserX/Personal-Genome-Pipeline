#!/usr/bin/env bash
# Step 33 (somalier, VerifyBamID2) on the fixture, with three controls:
#   - setup.sh --sample-qc-data installs somalier's sites and VerifyBamID2's
#     panel, each checked against its pinned sha256;
#   - the clean HG002 BAM of case 20: FREEMIX stays below 0.03. Only 2 of
#     somalier's chrX sites fall inside the slices, too few for a sex call, so
#     the step says it could not check the sex and goes on, also when male is
#     declared: somalier's table then carries male as the pedigree's sex,
#     which the verdict does not take for a call;
#   - sex: with sites at the slice's own chrX calls (case 21; HG002 is male,
#     so somalier sees them homozygous), declared female stops the step,
#     somalier's own table shows the pedigree's female next to the male it
#     set from the reads, declared male passes (and MultiQC, step 28, shows
#     male in its somalier Sex column), and SEX_CHECK=warn goes on;
#   - contamination: reads of HG001, an unrelated GIAB sample streamed from
#     GIAB's GRCh38 BAM, added to HG002's until they are about 10% of the mix:
#     FREEMIX rises above 0.03, which warns and never stops the step.
# The fixture's reference holds 14 of GRCh38's contigs, and somalier stops on
# a site whose contig the FASTA lacks ("sequence chr11 not found in fasta"),
# so every run here reads the installed sites and panel cut to those contigs
# (SOMALIER_SITES, VERIFYBAMID2_PANEL). They and the sites with the slice's
# chrX calls stay under reference/ for case qc-2.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
QC="${G}/${SAMPLE}/qc"
# tval FILE KEY: the value of KEY in a step 33 table.
tval() { awk -F'\t' -v k="$2" '$1 == k { print $2; exit }' "$1" 2>/dev/null; }
# scol COLUMN: COLUMN of somalier's samples.tsv for the sample.
scol() {
  awk -F'\t' -v k="$1" 'NR == 1 {sub(/^#/, ""); for (i = 1; i <= NF; i++) c[$i] = i; next} {print $c[k]; exit}' \
    "${QC}/somalier/${SAMPLE}.samples.tsv" 2>/dev/null
}
# below A B / above A B: numeric comparisons that fail on an empty value.
below() { [ -n "$1" ] && awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 < b + 0) }'; }
above() { [ -n "$1" ] && awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 > b + 0) }'; }

# --- 1. setup.sh installs the data -----------------------------------------------------
"${REPO}/scripts/setup.sh" --sample-qc-data "$G" > "${CASE_TMP}/setup.log" 2>&1
check_eq "setup.sh --sample-qc-data exits 0" "$?" 0
cat "${CASE_TMP}/setup.log"
SITES="${G}/reference/somalier/sites.hg38.vcf.gz"
PANEL="${G}/reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat"
check_ge "somalier sites" "$(gzip -dc "$SITES" 2>/dev/null | grep -vc '^#' || true)" 17000
for e in UD mu bed; do
  check_eq "VerifyBamID2 panel .${e} lines" "$(wc -l < "${PANEL}.${e}" 2>/dev/null | tr -d ' ')" 100000
done

# The installed files, cut to the fixture reference's contigs.
FIX_SITES="${G}/reference/somalier/sites_fixture.vcf"
FIX_PANEL="${G}/reference/verifybamid2_fixture/1000g.phase3.100k.b38.vcf.gz.dat"
mkdir -p "$(dirname "$FIX_PANEL")"
python3 - "${REF}.fai" "$SITES" "$FIX_SITES" "$PANEL" "$FIX_PANEL" <<'PY'
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
print(len(keep), "panel markers on the fixture's contigs")
PY
check_ge "somalier sites on the fixture's contigs" "$(grep -vc '^#' "$FIX_SITES" 2>/dev/null || true)" 9000
check_ge "panel markers on the fixture's contigs" "$(wc -l < "${FIX_PANEL}.bed" 2>/dev/null | tr -d ' ')" 40000
export SOMALIER_SITES="$FIX_SITES" VERIFYBAMID2_PANEL="$FIX_PANEL"

# --- 2. the clean sample -----------------------------------------------------------------
run_step 33-sample-qc.sh "$SAMPLE"
check_step_exit 33-sample-qc.sh
T="${QC}/${SAMPLE}_sample_qc.tsv"
cat "$T" 2>/dev/null
check_eq "somalier reads the sample under its @RG SM" "$(tval "$T" somalier_id)" "$SAMPLE"
check_ge "somalier sites genotyped on the slices" "$(tval "$T" sites_genotyped)" 60
check_eq "no sex call from 2 chrX sites" "$(tval "$T" inferred_sex)" unknown
check_eq "so the sex is not checked" "$(tval "$T" sex_check)" not_checked
check "FREEMIX of the clean sample is below 0.03 ($(tval "$T" freemix))" below "$(tval "$T" freemix)" 0.03
check_eq "contamination verdict" "$(tval "$T" contamination)" ok
check_eq "VerifyBamID2's marker check was skipped (the slices hold fewer than 1,000 markers)" \
  "$(tval "$T" verifybamid2_marker_check)" skipped
check "the log says why it ran again" has 'Fewer than 1,000 panel markers have reads' "$(cat "$STEP_LOG")"
check_eq "no declared sex: somalier's pedigree sex is unknown" "$(scol original_pedigree_sex)" unknown

# Declared male on the same 2 chrX sites: somalier starts its sex column from
# the pedigree's male and keeps it, as the reads cannot tell.
run_step 33-sample-qc.sh "$SAMPLE" male
check_step_exit 33-sample-qc.sh
check_eq "declared male: somalier's table has it as the pedigree's sex (MultiQC's Sex column)" \
  "$(scol original_pedigree_sex)" male
check_eq "declared male on 2 chrX sites: no sex call from the reads" "$(tval "$T" inferred_sex)" unknown
check_eq "so the sex is still not checked" "$(tval "$T" sex_check)" not_checked

# --- 3. sex, from sites at the slice's own chrX calls ----------------------------------------
# HG002 is male: case 21 calls chrX outside the pseudoautosomal regions
# haploid, and somalier finds those sites homozygous. Its rule: male when
# heterozygous / homozygous-ALT chrX sites is below 0.05 over more than 10.
XS="${G}/reference/somalier/sites_slice_chrX.vcf"
# Only calls with no reference read in DeepVariant's AD: a few chrX sites
# whose reads are mixed (paralogous mapping) read as heterozygous to
# somalier, and 2 such sites in 39 already break its male rule (below 0.05).
bcf view -H -f PASS -v snps -i 'FMT/AD[0:0]==0 && QUAL>=30' -r chrX:2781480-155701382 "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null \
  | awk -F'\t' -v OFS='\t' '{ split($5, a, ","); print $1, $2, ".", $4, a[1], ".", "PASS", "AF=0.5" }' \
  > "${CASE_TMP}/x_sites.tsv"
check_ge "chrX SNVs of case 21 outside the PARs, no reference read" "$(wc -l < "${CASE_TMP}/x_sites.tsv" | tr -d ' ')" 15
python3 - "$FIX_SITES" "${CASE_TMP}/x_sites.tsv" "${REF}.fai" "$XS" <<'PY'
import sys
sites, extra, fai, out = sys.argv[1:]
order = {l.split("\t")[0]: i for i, l in enumerate(open(fai))}
head, recs = [], []
for l in open(sites):
    if l.startswith("##"):
        head.append(l)
    elif l.startswith("#"):
        cols = l
    else:
        recs.append(l.rstrip("\n").split("\t")[:8])
seen = {(r[0], r[1]) for r in recs}
recs += [r.rstrip("\n").split("\t") for r in open(extra) if tuple(r.split("\t")[:2]) not in seen]
recs.sort(key=lambda r: (order.get(r[0], 1 << 30), int(r[1])))
with open(out, "w") as f:
    f.writelines(head)
    f.write(cols)
    f.writelines("\t".join(r) + "\n" for r in recs)
print(len(recs), "sites")
PY

SOMALIER_SITES="$XS" run_step 33-sample-qc.sh "$SAMPLE" female
check "declared female: the step stops (exit ${STEP_RC})" test "$STEP_RC" -ne 0
check "the message names both sexes" has 'SEX CHECK MISMATCH: declared female, somalier infers male' "$(cat "$STEP_LOG")"
check "and how to go on" has 'Set SEX_CHECK=warn' "$(cat "$STEP_LOG")"
cat "$T" 2>/dev/null
check_ge "chrX sites somalier genotyped" "$(tval "$T" x_sites)" 11
check_ge "homozygous ALT among them" "$(tval "$T" x_hom_alt)" 10
check_eq "the table records the mismatch" "$(tval "$T" sex_check)" mismatch
check_eq "somalier's table: the declared female as the pedigree's sex" "$(scol original_pedigree_sex)" female
check_eq "somalier's own check: it set the sex to male (1) from the reads" "$(scol sex)" 1
check "somalier's log says it changed the pedigree's sex" has "setting sex to male for ${SAMPLE}" "$(cat "$STEP_LOG")"

SOMALIER_SITES="$XS" run_step 33-sample-qc.sh "$SAMPLE" male
check_step_exit 33-sample-qc.sh
check_eq "declared male: the check passes" "$(tval "$T" sex_check)" ok
check_eq "somalier's table: the declared male as the pedigree's sex" "$(scol original_pedigree_sex)" male
# MultiQC reads that column as its somalier "Sex"; -9 before somalier had a pedigree.
FLAGSTAT="${G}/${SAMPLE}/aligned/${SAMPLE}_flagstat.txt"
HAD_FLAGSTAT=false; [ -f "$FLAGSTAT" ] && HAD_FLAGSTAT=true
run_step 28-multiqc.sh "$SAMPLE"
check_step_exit 28-multiqc.sh
MQC="${G}/${SAMPLE}/multiqc/multiqc_data/multiqc_somalier.txt"
check_eq "MultiQC's somalier table: Sex male" \
  "$(awk -F'\t' -v s="$SAMPLE" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $1 == s {print $c["original_pedigree_sex"]}' "$MQC" 2>/dev/null)" male
in_genome "$BCFTOOLS_IMAGE" rm -rf "${SAMPLE}/multiqc"
$HAD_FLAGSTAT || rm -f "$FLAGSTAT"
check "the log says so" has 'Sex check: OK' "$(cat "$STEP_LOG")"

SEX_CHECK=warn SOMALIER_SITES="$XS" run_step 33-sample-qc.sh "$SAMPLE" female
check_step_exit 33-sample-qc.sh
check "SEX_CHECK=warn: the mismatch is printed and the step goes on" has 'SEX_CHECK=warn: continuing' "$(cat "$STEP_LOG")"

# Leave the table of the default run for the report cases that follow.
run_step 33-sample-qc.sh "$SAMPLE"
check_step_exit 33-sample-qc.sh

# --- 4. contamination: HG001 reads mixed into HG002 ------------------------------------------
# HG001 (NA12878) is not related to HG002. Its GIAB GRCh38 BAM is streamed
# over the slices that hold most panel markers (HLA, chr20, chr2, chr12, chr1,
# chr19) and sampled to about a ninth of HG002's depth, so about 10% of the
# mixed reads are HG001's. Its read groups are dropped, so the mix reads as
# the one sample of HG002's header.
HG001_BAM=https://giab.s3.amazonaws.com/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam
REGIONS="chr6:29900000-33100000 chr20:10000000-10500000 chr2:233600000-233800000 chr12:47800000-47950000 chr1:109600000-109800000 chr19:40800000-41050000"
M="${SAMPLE}mix"
mkdir -p "${G}/${M}/aligned" "${G}/contam"
HG002_DP=$(sam coverage -r chr20:10000000-10500000 "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" 2>/dev/null | awk 'NR == 2 { print $7 }')
echo "HG002 depth on the chr20 slice: ${HG002_DP:-?}x"
# GIAB's HG001 BAM is about 300x; a ninth of HG002's depth from it:
FRACTION=$(awk -v d="${HG002_DP:-30}" 'BEGIN { printf "%.4f", d / 9 / 300 }')
echo "streaming HG001 at a fraction of ${FRACTION}"
curl -fsSL --retry 5 -o "${G}/contam/hg001.bai" "${HG001_BAM}.bai"
# Three attempts: one network hiccup on the NCBI stream must not fail the job.
for attempt in 1 2 3; do
  # shellcheck disable=SC2086  # the regions split on purpose
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
    -v /etc/ssl/certs:/etc/ssl/certs:ro -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt \
    -v "${G}:/genome" -w /genome "$SAMTOOLS_IMAGE" \
    samtools view -M -s "7${FRACTION#0}" -x RG -F 0x900 -X "$HG001_BAM" /genome/contam/hg001.bai $REGIONS \
    > "${G}/contam/hg001.sam" 2> "${CASE_TMP}/stream.log"
  STREAM_RC=$?
  [ "$STREAM_RC" -eq 0 ] && break
  echo "HG001 stream attempt ${attempt} failed (exit ${STREAM_RC}):"
  tail -n 3 "${CASE_TMP}/stream.log"
  [ "$attempt" -lt 3 ] && sleep 30
done
check_eq "HG001 reads streamed" "$STREAM_RC" 0
tail -n 3 "${CASE_TMP}/stream.log"
N_HG001=$(wc -l < "${G}/contam/hg001.sam" | tr -d ' ')
# shellcheck disable=SC2086
N_HG002=$(sam view -c -F 0x900 "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" $REGIONS 2>/dev/null)
echo "reads in those regions: HG002 ${N_HG002:-?}, HG001 ${N_HG001}"
SHARE=$(awk -v a="$N_HG001" -v b="${N_HG002:-0}" 'BEGIN { printf "%.3f", (a + b) ? a / (a + b) : 0 }')
check "HG001's share of the mixed reads is about 10% (${SHARE})" \
  awk -v s="$SHARE" 'BEGIN { exit !(s >= 0.07 && s <= 0.14) }'
in_genome "$SAMTOOLS_IMAGE" sh -c "set -e
  { samtools view -H ${SAMPLE}/aligned/${SAMPLE}_sorted.bam; samtools view ${SAMPLE}/aligned/${SAMPLE}_sorted.bam; cat contam/hg001.sam; } \
    | samtools sort -@ 2 -m 1G -o ${M}/aligned/${M}_sorted.bam -
  samtools index ${M}/aligned/${M}_sorted.bam"
rm -f "${G}/contam/hg001.sam"
check "the mixed BAM passes quickcheck" sam quickcheck "${M}/aligned/${M}_sorted.bam"

run_step 33-sample-qc.sh "$M"
check_step_exit 33-sample-qc.sh
TM="${G}/${M}/qc/${M}_sample_qc.tsv"
cat "$TM" 2>/dev/null
check "FREEMIX of the mix is above 0.03 ($(tval "$TM" freemix))" above "$(tval "$TM" freemix)" 0.03
check_eq "contamination verdict" "$(tval "$TM" contamination)" warn
check "the step warns, and still exits 0" has 'is above 0\.03: about [0-9]+% of the reads may come' "$(cat "$STEP_LOG")"
echo "FREEMIX: clean $(tval "$T" freemix), with ${SHARE} HG001 reads $(tval "$TM" freemix)" >> "$E2E_NOTES"
in_genome "$BCFTOOLS_IMAGE" rm -rf "$M" contam

finish
