#!/usr/bin/env bash
# Step 38 and the Nextflow family workflow (--pedigree) on the GIAB Ashkenazi
# trio: HG002, the fixture's sample, and his parents HG003 and HG004, whose
# reads are streamed from GIAB's GRCh38 BAMs (public, consented samples) over
# ten of the fixture's slices, 2.45 Mb in all, sampled to about 30x, then
# aligned (step 02) and called with their sex (step 03). It all happens in a
# GENOME_DIR of its own, so no later case sees the parents.
#   - the joint call holds the three; somalier measures both parent-child
#     pairs as first degree and the parents as unrelated, and the check
#     passes. The slices hold too few of somalier's own sites, so the sites
#     here are the SNVs of HG001, an unrelated GIAB sample, in those slices;
#   - red control: HG002 and HG003 swapped in the pedigree. Step 38 and the
#     pipeline both stop at the relationship check and name the pair, and
#     nothing is phased;
#   - WhatsHap phases HG002's heterozygous sites by inheritance, and its phase
#     agrees with the parent of origin the genotypes give; on chrX outside the
#     PARs HG002's alleles come from his mother only;
#   - a planted Mendelian violation (one of HG004's heterozygous calls that
#     HG002 inherited, turned into a reference block in a copy of her gVCF)
#     is a de novo candidate of the pipeline's run, and not of the clean run;
#   - parity: step 38 and the pipeline write the same relationships,
#     inheritance table, de novo list and joint and phased records;
#   - a single-sample run without --pedigree publishes the same files as main
#     does, after dropping the lines that hold a date or a command line; a
#     planted difference shows that the comparison sees one.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G0="$GENOME_DIR"
F="${E2E_WORK}/genome-family"
REF_NAME=GRCh38_no_alt_analysis_set
REGIONS=(chr1:109600001-109800000 chr2:233600001-233800000 chr4:68500001-68700000 chr10:94700001-95000000
         chr12:47800001-47950000 chr16:28580001-28630000 chr19:40800001-41050000 chr20:10000001-10500000
         chr22:42000001-42300000 chrX:73700001-74000000)
GIAB=https://giab.s3.amazonaws.com/data/AshkenazimTrio
declare -A BAM_URL=(
  [HG003]="${GIAB}/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.60x.1.bam"
  [HG004]="${GIAB}/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.60x.1.bam"
)
declare -A SEX=([HG002]=male [HG003]=male [HG004]=female)

# --- 1. A GENOME_DIR for the trio ---------------------------------------------------------
in_e2e() { docker run --rm -v "${E2E_WORK}:/w" "$BCFTOOLS_IMAGE" "$@"; }
in_e2e rm -rf /w/genome-family
mkdir -p "${F}/reference" "${F}/HG002/aligned"
for f in "${G0}/reference/${REF_NAME}".*; do ln -f "$f" "${F}/reference/" 2>/dev/null || cp -f "$f" "${F}/reference/"; done
for f in "${G0}/HG002/aligned/HG002_sorted.bam" "${G0}/HG002/aligned/HG002_sorted.bam.bai"; do
  ln -f "$f" "${F}/HG002/aligned/" 2>/dev/null || cp -f "$f" "${F}/HG002/aligned/"
done
check "the reference and its minimap2 index are in the trio's GENOME_DIR" \
  test -s "${F}/reference/${REF_NAME}.fasta" -a -s "${F}/reference/${REF_NAME}.sr.mmi"
export GENOME_DIR="$F"

# --- 2. HG003 and HG004: stream, sample to about 30x, align, call ----------------------------------
for s in HG003 HG004; do
  mkdir -p "${F}/${s}/fastq" "${F}/stream"
  curl -fsSL --retry 5 -o "${F}/stream/${s}.bai" "${BAM_URL[$s]}.bai"
  for attempt in 1 2 3; do
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
      -v /etc/ssl/certs:/etc/ssl/certs:ro -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt \
      -v "${F}:/genome" -w /genome "$SAMTOOLS_IMAGE" sh -c "set -e
        samtools view -M -u -s 7.5 -F 0x900 -X '${BAM_URL[$s]}' stream/${s}.bai ${REGIONS[*]} \
          | samtools collate -u -O - stream/${s}.collate \
          | samtools fastq -n -1 ${s}/fastq/${s}_R1.fastq.gz -2 ${s}/fastq/${s}_R2.fastq.gz -0 /dev/null -s /dev/null -" \
      > "${CASE_TMP}/stream.log" 2>&1
    rc=$?
    [ "$rc" -eq 0 ] && break
    echo "${s} stream attempt ${attempt} failed (exit ${rc}):"
    tail -n 3 "${CASE_TMP}/stream.log"
    [ "$attempt" -lt 3 ] && sleep 30
  done
  check_eq "${s} reads streamed from GIAB" "$rc" 0
  check_ge "${s} read pairs" "$(( $(gzip -dc "${F}/${s}/fastq/${s}_R1.fastq.gz" 2>/dev/null | wc -l) / 4 ))" 100000
  run_step 02-alignment.sh "$s"
  check_step_exit 02-alignment.sh
done
in_e2e rm -rf /w/genome-family/stream
INTERVALS="${REGIONS[*]}"
export INTERVALS
for s in HG002 HG003 HG004; do
  run_step 03-deepvariant.sh "$s" "${SEX[$s]}"
  check_step_exit 03-deepvariant.sh
  check "${s}: gVCF and index" test -s "${F}/${s}/vcf/${s}.g.vcf.gz" -a -s "${F}/${s}/vcf/${s}.g.vcf.gz.tbi"
  echo "${s} depth on the chr20 slice: $(sam coverage -r chr20:10000001-10500000 "${s}/aligned/${s}_sorted.bam" 2>/dev/null | awk 'NR == 2 {print $7}')x"
done
unset INTERVALS

# somalier's sites: GIAB's HG001 benchmark SNVs inside the autosomal slices,
# streamed from GIAB. HG001 is not related to the trio. Sites taken from the
# trio's own calls bias somalier's relatedness down (0.353 and 0.385 for the
# two parent-child pairs, -0.147 for the parents, in the first run of this
# case), because a site is kept only when one of the three carries it.
SITES="${F}/reference/somalier/sites_hg001.vcf"
HG001_VCF=https://giab.s3.amazonaws.com/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
AUTO=$(printf '%s\n' "${REGIONS[@]}" | grep -v '^chrX' | paste -sd, -)
mkdir -p "$(dirname "$SITES")" "${F}/hg001"
# The whole file (126 MB) and its index, then a local query of the slices
curl -fsSL --retry 5 -o "${F}/hg001/hg001.vcf.gz" "$HG001_VCF"
curl -fsSL --retry 5 -o "${F}/hg001/hg001.vcf.gz.tbi" "${HG001_VCF}.tbi"
bcf view -H -v snps -m2 -M2 -r "$AUTO" hg001/hg001.vcf.gz > "${CASE_TMP}/hg001.vcf" 2> "${CASE_TMP}/hg001.log"
check_eq "HG001's benchmark SNVs read from GIAB's file" "$?" 0
cut -f1,2,4,5 "${CASE_TMP}/hg001.vcf" > "${CASE_TMP}/snvs.tsv"
python3 - "${F}/reference/${REF_NAME}.fasta.fai" "${CASE_TMP}/snvs.tsv" "$SITES" <<'PY'
import sys
fai, snvs, out = sys.argv[1:]
order, lens = {}, {}
for i, l in enumerate(open(fai)):
    c = l.split("\t")
    order[c[0]], lens[c[0]] = i, c[1]
keep = {tuple(l.rstrip("\n").split("\t")) for l in open(snvs) if l.strip()}
keep = sorted((k for k in keep if k[0] in order), key=lambda k: (order[k[0]], int(k[1])))
with open(out, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for c in sorted(order, key=order.get):
        f.write("##contig=<ID=%s,length=%s>\n" % (c, lens[c]))
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for c, p, r, a in keep:
        f.write("\t".join([c, p, ".", r, a, ".", "PASS", "AF=0.5"]) + "\n")
print(len(keep), "sites")
PY
check_ge "somalier sites: HG001's SNVs in the autosomal slices" "$(grep -vc '^#' "$SITES")" 1500

printf 'sample\tsex\tfather\tmother\tpartner\nHG002\tmale\tHG003\tHG004\t.\nHG003\tmale\t.\t.\tHG004\nHG004\tfemale\t.\t.\tHG003\n' \
  > "${F}/pedigree.tsv"
# HG002 and HG003 swapped: "HG003" is now the child of "HG002" and HG004
printf 'sample\tsex\tfather\tmother\tpartner\nHG003\tmale\tHG002\tHG004\t.\nHG002\tmale\t.\t.\tHG004\nHG004\tfemale\t.\t.\tHG002\n' \
  > "${F}/pedigree_swapped.tsv"

# tval FILE A B COLUMN: COLUMN of the relationships.tsv row of the pair A, B
tval() {
  awk -F'\t' -v a="$2" -v b="$3" -v k="$4" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    ($1 == a && $2 == b) || ($1 == b && $2 == a) {print $c[k]; exit}' "$1" 2>/dev/null
}
between() { [ -n "$1" ] && awk -v x="$1" -v lo="$2" -v hi="$3" 'BEGIN { exit !(x + 0 >= lo && x + 0 <= hi) }'; }

# --- 3. Step 38 ------------------------------------------------------------------------------
export SOMALIER_SITES="$SITES" THREADS=4 GLNEXUS_MEM=8
FAM="${F}/trio"
run_step 38-family.sh trio "${F}/pedigree.tsv"
check_step_exit 38-family.sh
check_eq "the joint VCF holds the trio" "$(bcf query -l trio/family.vcf.gz 2>/dev/null | tr '\n' ' ')" "HG002 HG003 HG004 "
check_ge "joint VCF records" "$(vcf_count trio/family.vcf.gz)" 1500
R="${FAM}/relationships.tsv"
cat "$R" 2>/dev/null
for p in "HG002 HG003" "HG002 HG004"; do
  # shellcheck disable=SC2086  # the pair splits on purpose
  check "${p/ / and }: declared parent-child, relatedness $(tval "$R" $p relatedness) is first degree" \
    between "$(tval "$R" $p relatedness)" 0.354 0.707
  # shellcheck disable=SC2086
  check_eq "${p/ / and }: verdict" "$(tval "$R" $p verdict)" ok
  # shellcheck disable=SC2086
  check_ge "${p/ / and }: sites compared" "$(tval "$R" $p sites)" 1000
done
check "HG003 and HG004 (partners): relatedness $(tval "$R" HG003 HG004 relatedness) is unrelated" \
  between "$(tval "$R" HG003 HG004 relatedness)" -1 0.0884
check_eq "HG003 and HG004: verdict" "$(tval "$R" HG003 HG004 verdict)" ok
MD="${FAM}/family_report.md"
PH=$(sed -nE 's/.*WhatsHap phased ([0-9]+) heterozygous sites of HG002 by inheritance\. At ([0-9]+) of them the genotypes also give the origin, and the phase agrees with it at ([0-9]+).*/\1 \2 \3/p' "$MD" 2>/dev/null)
read -r PH_N PH_INF PH_AGREE <<< "${PH:-0 0 0}"
check "WhatsHap wrote the phased VCF" test -s "${FAM}/family.phased.vcf.gz"
check_ge "HG002 heterozygous sites WhatsHap phased by inheritance" "${PH_N:-0}" 300
check "the phase agrees with the genotype origin at 98% or more (${PH_AGREE:-?} of ${PH_INF:-?})" \
  awk -v a="${PH_AGREE:-0}" -v n="${PH_INF:-0}" 'BEGIN { exit !(n >= 300 && a >= 0.98 * n) }'
INH="${FAM}/HG002_inheritance.tsv.gz"
icount() { gzip -dc "$INH" 2>/dev/null | awk -F'\t' "NR > 1 && ($1)" | wc -l | tr -d ' '; }
check_ge "HG002 sites with the father's allele known" "$(icount '$8 != "."')" 500
check_ge "HG002 sites with the mother's allele known" "$(icount '$9 != "."')" 500
check_ge "chrX outside the PARs: alleles from the mother" "$(icount '$1 == "chrX" && $10 == "genotype" && $9 != "."')" 20
check_eq "chrX outside the PARs: none from the father (HG002 is a son)" "$(icount '$1 == "chrX" && $8 != "."')" 0
check "the report states the de novo candidates" has 'De novo candidates \(an allele neither parent has\): [0-9]+' "$(cat "$MD" 2>/dev/null)"
grep -A3 '^De novo' "$MD" 2>/dev/null | head -n 3
DN_CLEAN=$(awk 'END {print NR - 1}' "${FAM}/de_novo.tsv" 2>/dev/null)
echo "family: de novo candidates of the clean trio ${DN_CLEAN}, HG002 phased hets ${PH_N}, agreement ${PH_AGREE}/${PH_INF}" >> "$E2E_NOTES"

run_step 38-family.sh trio_swapped "${F}/pedigree_swapped.tsv"
check "step 38, HG002 and HG003 swapped: it stops (exit ${STEP_RC})" test "$STEP_RC" -eq 1
check "and names a pair that disagrees" has 'HG003 and HG004: declared parent-child \(first degree\), but somalier measures (unrelated|third degree)' "$(cat "$STEP_LOG")"
check "and says what it may be" has 'A sample swap' "$(cat "$STEP_LOG")"
check_eq "its relationships.tsv marks the mismatches" "$(grep -c $'\tmismatch\t' "${F}/trio_swapped/relationships.tsv" 2>/dev/null)" 2
check "nothing was phased or reported after the stop" \
  test ! -e "${F}/trio_swapped/family.phased.vcf.gz" -a ! -e "${F}/trio_swapped/family_report.md"

# --- 4. The pipeline ----------------------------------------------------------------------------
# nf_family NAME PEDIGREE [HG004_GVCF]: the trio as VCF rows with their gVCFs.
nf_family() {
  local name=$1 ped=$2 hg004=${3:-${F}/HG004/vcf/HG004.g.vcf.gz} dir="${CASE_TMP}/$1" s g
  rm -rf "$dir"
  mkdir -p "$dir"
  echo 'sample,vcf,vcf_index,gvcf,gvcf_index,sex' > "${dir}/samplesheet.csv"
  for s in HG002 HG003 HG004; do
    g="${F}/${s}/vcf/${s}.g.vcf.gz"
    [ "$s" = HG004 ] && g=$hg004
    printf '%s,%s,%s,%s,%s,%s\n' "$s" "${F}/${s}/vcf/${s}.vcf.gz" "${F}/${s}/vcf/${s}.vcf.gz.tbi" "$g" "${g}.tbi" "${SEX[$s]}" \
      >> "${dir}/samplesheet.csv"
  done
  OUT="${dir}/out"
  echo "+ nextflow run main.nf (${name}) --pedigree ${ped##*/}"
  ( cd "$dir" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${F}/nfwork" --input "${dir}/samplesheet.csv" \
      --reference "${F}/reference/${REF_NAME}.fasta" --tools roh \
      --pedigree "$ped" --somalier_sites "$SITES" \
      --outdir "$OUT" --max_cpus 4 --max_memory 14.GB ) > "${dir}/run.log" 2>&1
  RC=$?
  LOG="${dir}/run.log"
  grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "$LOG"
  cp "${dir}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.${name}.nextflow.log" 2>/dev/null
}

nf_family clean "${F}/pedigree.tsv"
check_eq "pipeline, the trio: exit" "$RC" 0
NF="${OUT}/family"
cmp_items=(relationships.tsv HG002_inheritance.tsv.gz de_novo.tsv)
for f in "${cmp_items[@]}"; do
  check "parity: ${f} of step 38 and of the pipeline are the same bytes" cmp -s "${FAM}/${f}" "${NF}/${f}"
done
records() { gzip -dc "$1" 2>/dev/null | grep -v '^#' | md5sum | cut -d' ' -f1; }
check_ge "the pipeline's joint VCF has records" "$(gzip -dc "${NF}/family.vcf.gz" 2>/dev/null | grep -vc '^#')" 1500
check_ge "the pipeline's phased VCF has records" "$(gzip -dc "${NF}/family.phased.vcf.gz" 2>/dev/null | grep -vc '^#')" 1500
check_eq "parity: the joint VCF records" "$(records "${NF}/family.vcf.gz")" "$(records "${FAM}/family.vcf.gz")"
check_eq "parity: the phased VCF records" "$(records "${NF}/family.phased.vcf.gz")" "$(records "${FAM}/family.phased.vcf.gz")"
check "parity check sees a difference: a planted one is caught" \
  sh -c "gzip -dc '${FAM}/HG002_inheritance.tsv.gz' | sed '2s/\$/x/' | gzip -n > '${CASE_TMP}/planted.tsv.gz'; ! cmp -s '${CASE_TMP}/planted.tsv.gz' '${NF}/HG002_inheritance.tsv.gz'"
check "the pipeline published somalier's tables under family/somalier" test -s "${NF}/somalier/family.pairs.tsv"

nf_family swapped "${F}/pedigree_swapped.tsv"
check "pipeline, HG002 and HG003 swapped: the run stops (exit ${RC})" test "$RC" -ne 0
check "the message names a pair that disagrees" has 'HG003 and HG004: declared parent-child \(first degree\), but somalier measures (unrelated|third degree)' "$(cat "$LOG")"
check "FAMILY_CHECK is the task that failed" has "FAMILY:FAMILY_CHECK" "$(cat "$LOG")"
check "no report and no phasing after the stop" test ! -e "${OUT}/family/family_report.md" -a ! -e "${OUT}/family/family.phased.vcf.gz"

# A planted de novo: a site HG002 inherited from HG004 (HG002 0/1, HG003 0/0,
# HG004 0/1, all with GQ 30 or more), made a reference block in a copy of
# HG004's gVCF.
SITE=$(bcf query -s HG002,HG003,HG004 -i 'TYPE="snp"' -f '%CHROM\t%POS[\t%GT\t%GQ]\n' -r chr20 trio/family.vcf.gz 2>/dev/null \
  | awk -F'\t' '$3 == "0/1" && $5 == "0/0" && $7 == "0/1" && $4 >= 30 && $6 >= 30 && $8 >= 30 {print $1 "\t" $2; exit}')
read -r P_CHR P_POS <<< "${SITE:-none 0}"
check "a site HG002 inherited from HG004 to plant at (${P_CHR}:${P_POS})" test "$P_POS" -gt 0
check_eq "it is not a de novo candidate of the clean run" \
  "$(awk -F'\t' -v c="$P_CHR" -v p="$P_POS" '$2 == c && $3 == p' "${NF}/de_novo.tsv" 2>/dev/null | wc -l | tr -d ' ')" 0
mkdir -p "${F}/planted"
gzip -dc "${F}/HG004/vcf/HG004.g.vcf.gz" | awk -F'\t' -v OFS='\t' -v c="$P_CHR" -v p="$P_POS" '
  !/^#/ && $1 == c && $2 == p { print $1, $2, ".", $4, "<*>", "0", ".", "END=" p, "GT:GQ:MIN_DP:PL", "0/0:50:30:0,90,990"; n++; next }
  { print }
  END { if (n != 1) exit 1 }' > "${F}/planted/HG004.g.vcf"
check "HG004's record at the site became a reference block" test "$?" -eq 0
bcf view -Oz -o planted/HG004.g.vcf.gz planted/HG004.g.vcf
bcf index -t planted/HG004.g.vcf.gz
nf_family planted "${F}/pedigree.tsv" "${F}/planted/HG004.g.vcf.gz"
check_eq "pipeline, the planted trio: exit" "$RC" 0
DN=$(awk -F'\t' -v c="$P_CHR" -v p="$P_POS" '$1 == "HG002" && $2 == c && $3 == p' "${OUT}/family/de_novo.tsv" 2>/dev/null)
echo "${DN:-no row}"
check "the planted site is a de novo candidate of HG002" test -n "$DN"
check_eq "with both parents homozygous reference" "$(cut -f7,8 <<< "$DN" | cut -d: -f1 | tr '\t' ' ')" "0/0 0/0"

# --- 5. A single-sample run is what main publishes ---------------------------------------------
# norm DIR: one line per published file (outside pipeline_info), its path and
# the md5 of its text without the lines that hold a date or a command line.
norm() {
  ( cd "$1" && find . -type f ! -path './pipeline_info/*' | LC_ALL=C sort | while read -r f; do
      printf '%s %s\n' "$f" "$(gzip -cdf "$f" | grep -avE 'Date=|^##[A-Za-z_]*[Cc]ommand(line)?=|^# The command line' | md5sum | cut -d' ' -f1)"
    done )
}
MAIN="${CASE_TMP}/main"
rm -rf "$MAIN"
mkdir -p "$MAIN"
if git -C "$REPO" fetch -q --depth 1 https://github.com/GeiserX/Personal-Genome-Pipeline.git main 2>"${CASE_TMP}/fetch.log"; then
  git -C "$REPO" archive FETCH_HEAD | tar -x -C "$MAIN"
fi
check "main's pipeline is checked out" test -s "${MAIN}/main.nf"
single() {
  local tree=$1 dir="${CASE_TMP}/single-$2"
  rm -rf "$dir"
  mkdir -p "$dir"
  printf 'sample,vcf,vcf_index,sex\nHG002,%s,%s,male\n' "${F}/HG002/vcf/HG002.vcf.gz" "${F}/HG002/vcf/HG002.vcf.gz.tbi" > "${dir}/samplesheet.csv"
  ( cd "$dir" && nextflow run "${tree}/main.nf" -profile docker -ansi-log false -work-dir "${dir}/work" \
      --input "${dir}/samplesheet.csv" --reference "${F}/reference/${REF_NAME}.fasta" --tools roh,clinvar \
      --clinvar "${G0}/clinvar/clinvar_pathogenic_chr.vcf.gz" --clinvar_index "${G0}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
      --outdir "${dir}/out" --max_cpus 4 --max_memory 14.GB ) > "${dir}/run.log" 2>&1
  SRC=$?
  grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "${dir}/run.log" | tail -n 5
  norm "${dir}/out" > "${dir}/files.txt"
}
single "$MAIN" main
check_eq "main, one sample, no --pedigree: exit" "$SRC" 0
single "$REPO" head
check_eq "this branch, one sample, no --pedigree: exit" "$SRC" 0
check_ge "files the single-sample run publishes" "$(wc -l < "${CASE_TMP}/single-head/files.txt" | tr -d ' ')" 3
check "a single-sample run publishes the same files, with the same content, as main" \
  diff "${CASE_TMP}/single-main/files.txt" "${CASE_TMP}/single-head/files.txt"
check "no family/ folder without --pedigree" test ! -e "${CASE_TMP}/single-head/out/family"
F1=$(awk 'NR == 1 {print $1}' "${CASE_TMP}/single-head/files.txt")
cp -r "${CASE_TMP}/single-head/out" "${CASE_TMP}/single-planted"
echo "planted" >> "${CASE_TMP}/single-planted/${F1#./}"
norm "${CASE_TMP}/single-planted" > "${CASE_TMP}/single-planted.txt"
check "the comparison sees a planted difference (a line added to ${F1})" \
  sh -c "! diff -q '${CASE_TMP}/single-head/files.txt' '${CASE_TMP}/single-planted.txt' > /dev/null"
in_e2e rm -rf /w/genome-family
finish
