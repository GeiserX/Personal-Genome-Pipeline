#!/usr/bin/env bash
# build-fixture.sh — build the small HG002 data set the e2e job runs every tool on.
#
# Usage: scripts/ci/build-fixture.sh <out_dir>
#
# Streams a few regions of the public GIAB HG002 60x GRCh38 BAM over HTTPS
# (samtools reads only the byte ranges it needs, through the .bai), samples them
# down to about 30x and writes into <out_dir>:
#
#   HG002_R1.fastq.gz, HG002_R2.fastq.gz   name-sorted read pairs, for step 02
#   HG002_slice.bam (+.bai)                the downsampled GIAB alignments
#   fixture_ref.fa.gz (+.fai .gzi .dict)   whole chr5 chr6 chr10 chr12 chr20
#                                          chr22 chrX chrY chrM of the NCBI
#                                          GRCh38 no-alt analysis set
#   clinvar.vcf.gz, clinvar_chr.vcf.gz,    ClinVar records inside the regions,
#   clinvar_pathogenic_chr.vcf.gz (+.tbi)  plus one planted record (planted.tsv)
#   HG002_vep.vcf                          up to 200 GIAB truth variants annotated
#                                          by VEP --database --everything
#   revel_synthetic.tsv.gz (+.tbi)         a tiny score file in REVEL's layout;
#                                          the values are made up, not REVEL's
#   HG002_sv_manta_style.vcf.gz (+.tbi)    ten Manta-style SV records
#   HG002_truth_chr20.vcf.gz (+.tbi), HG002_truth_chr20.bed
#                                          GIAB v4.2.1 truth for the chr20 slice
#   regions.bed, planted.tsv, MANIFEST.txt, SHA256SUMS
#
# HG002 is a public, consented Genome in a Bottle sample, so nothing here is
# personal data. Needs Docker, curl, bgzip and tabix (htslib), about 10 GB of
# free disk and a network connection (VEP queries Ensembl's public database).
# The release tag the e2e job downloads is tests/fixtures/VERSION; see
# docs/testing.md for how to rebuild and publish a new version.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
# shellcheck source=../../versions.env
. "${REPO}/versions.env"

OUT_ARG=${1:?Usage: $0 <out_dir>}
mkdir -p "$OUT_ARG"
OUT="$(cd "$OUT_ARG" && pwd)"
WORK="${OUT}/.work"
mkdir -p "$WORK"

VERSION=$(tr -d '[:space:]' < "${REPO}/tests/fixtures/VERSION")
THREADS=${THREADS:-4}
TARGET_DEPTH=${TARGET_DEPTH:-30}
SEED=42
SAMPLE=HG002
MAX_TOTAL_BYTES=$((1500 * 1024 * 1024))

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
BAM_URL="${GIAB}/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
TRUTH_BASE="${GIAB}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark"
REF_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
REF_NAME=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
CLINVAR_URL=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

# Slices, chosen so later packages (paralogs, ploidy, sample QC, Y haplogroup)
# need no rebuild. Coordinates are GRCh38, 1-based, inclusive.
REGIONS=(
  "chr20:10000000-10500000"   # small variants; the planted ClinVar record (SNAP25)
  "chr22:42000000-42300000"   # CYP2D6 and CYP2D7 with flanks
  "chr12:47800000-47950000"   # VDR, pypgx's control gene
  "chr10:94700000-95000000"   # CYP2C19 and CYP2C9
  "chr6:29900000-33100000"    # HLA
  "chr5:69900000-71100000"    # SMN1 and SMN2
  "chrX:73700000-74000000"    # non-PAR chrX (XIST)
  "chrX:1000000-1200000"      # PAR1
  "chrY:2700000-3000000"      # non-PAR chrY (SRY, RPS4Y1, ZFY)
  "chrM"                      # whole mitochondrial genome
)
CONTIGS=(chr5 chr6 chr10 chr12 chr20 chr22 chrX chrY chrM)

# Truth variants sent to VEP: gene windows and how many records to take from
# each (200 at most in total). The HLA genes give missense variants, so steps
# 23 and 31 have MODERATE-impact records and compound-het candidates.
VEP_WINDOWS=(
  "chr6:29941260-29945884 60"    # HLA-A
  "chr6:31353872-31357188 60"    # HLA-B
  "chr6:31268749-31272130 30"    # HLA-C
  "chr10:94762681-94855547 20"   # CYP2C19
  "chr10:94938658-94990091 20"   # CYP2C9
  "chr20:10218830-10307418 9"    # SNAP25
)
PLANT_WINDOW="chr20:10230000-10300000"   # inside SNAP25
PLANT_GENE="SNAP25:6616"
PLANT_ID=900000001

for tool in docker curl bgzip tabix md5sum sha256sum awk sort; do
  command -v "$tool" >/dev/null || { echo "ERROR: ${tool} is required" >&2; exit 1; }
done

# Run a tool from a pinned image, as the calling user, with OUT mounted at /w.
in_image() {
  local image=$1; shift
  docker run --rm -i -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${OUT}:/w" -w /w "$image" "$@"
}
sam() { in_image "$SAMTOOLS_IMAGE" samtools "$@"; }
bcf() { in_image "$BCFTOOLS_IMAGE" bcftools "$@"; }
fetch() { curl -fsSL --retry 5 --retry-delay 10 -o "$2" "$1"; }

echo "=== Fixture ${VERSION}: ${OUT} ==="
sam --version | head -3

# --- Reference: whole contigs from the no-alt analysis set ------------------
echo "[1/8] Reference contigs: ${CONTIGS[*]}"
fetch "${REF_BASE}/${REF_NAME}" "${WORK}/${REF_NAME}"
REF_MD5=$(curl -fsSL "${REF_BASE}/md5checksums.txt" | awk -v f="./${REF_NAME}" '$2 == f {print $1}')
GOT_MD5=$(md5sum "${WORK}/${REF_NAME}" | awk '{print $1}')
if [ -z "$REF_MD5" ] || [ "$REF_MD5" != "$GOT_MD5" ]; then
  echo "ERROR: reference md5 ${GOT_MD5} does not match NCBI's '${REF_MD5}'" >&2
  exit 1
fi
gzip -dc "${WORK}/${REF_NAME}" > "${WORK}/full.fa"
rm -f "${WORK}/${REF_NAME}"
sam faidx /w/.work/full.fa
sam faidx -o /w/.work/mini.fa /w/.work/full.fa "${CONTIGS[@]}"
rm -f "${WORK}/full.fa" "${WORK}/full.fa.fai"
bgzip -@ "$THREADS" -c "${WORK}/mini.fa" > "${OUT}/fixture_ref.fa.gz"
rm -f "${WORK}/mini.fa"
sam faidx /w/fixture_ref.fa.gz
sam dict -o /w/fixture_ref.dict /w/fixture_ref.fa.gz

# --- Reads: stream the slices, sample down to TARGET_DEPTH -------------------
echo "[2/8] Streaming ${#REGIONS[@]} regions from the GIAB HG002 60x BAM"
sam view -@ "$THREADS" -b -o /w/.work/slice_full.bam "$BAM_URL" "${REGIONS[@]}"
sam index /w/.work/slice_full.bam
FULL_DEPTH=$(sam coverage -r chr20:10000000-10500000 /w/.work/slice_full.bam | awk 'NR == 2 {print $7}')
# samtools -s takes SEED.FRACTION; mates share a read name, so pairs stay whole.
FRACTION=$(awk -v d="$FULL_DEPTH" -v t="$TARGET_DEPTH" 'BEGIN {f = t / d; if (f > 0.9999) f = 0.9999; printf "%.4f", f}')
echo "  chr20 slice depth ${FULL_DEPTH}x; keeping a fraction of ${FRACTION} (seed ${SEED})"
sam view -@ "$THREADS" -b -s "${SEED}${FRACTION#0}" \
  -o "/w/${SAMPLE}_slice.bam" /w/.work/slice_full.bam
sam index "/w/${SAMPLE}_slice.bam"
rm -f "${WORK}/slice_full.bam" "${WORK}/slice_full.bam.bai" "${WORK}"/*.bai
SLICE_DEPTH=$(sam coverage -r chr20:10000000-10500000 "/w/${SAMPLE}_slice.bam" | awk 'NR == 2 {print $7}')

echo "[3/8] Paired FASTQ (name-collated, secondary and supplementary records dropped)"
sam collate -u -O "/w/${SAMPLE}_slice.bam" /w/.work/collate \
  | sam fastq -n -c 6 -1 "/w/${SAMPLE}_R1.fastq.gz" -2 "/w/${SAMPLE}_R2.fastq.gz" \
      -0 /dev/null -s /dev/null -
R1_READS=$(( $(gzip -dc "${OUT}/${SAMPLE}_R1.fastq.gz" | wc -l) / 4 ))
R2_READS=$(( $(gzip -dc "${OUT}/${SAMPLE}_R2.fastq.gz" | wc -l) / 4 ))
if [ "$R1_READS" -ne "$R2_READS" ] || [ "$R1_READS" -eq 0 ]; then
  echo "ERROR: R1 has ${R1_READS} reads and R2 ${R2_READS}" >&2
  exit 1
fi

# --- GIAB truth: chr20 slice, the planted variant, the VEP input -------------
echo "[4/8] GIAB v4.2.1 truth"
fetch "${TRUTH_BASE}.vcf.gz" "${WORK}/truth.vcf.gz"
fetch "${TRUTH_BASE}.vcf.gz.tbi" "${WORK}/truth.vcf.gz.tbi"
fetch "${TRUTH_BASE}_noinconsistent.bed" "${WORK}/truth.bed"
bcf view -r chr20:10000000-10500000 -Oz -o "/w/${SAMPLE}_truth_chr20.vcf.gz" /w/.work/truth.vcf.gz
bcf index -t "/w/${SAMPLE}_truth_chr20.vcf.gz"
awk 'BEGIN {OFS = "\t"} $1 == "chr20" && $3 > 9999999 && $2 < 10500000 {
       if ($2 < 9999999) $2 = 9999999; if ($3 > 10500000) $3 = 10500000; print }' \
  "${WORK}/truth.bed" > "${OUT}/${SAMPLE}_truth_chr20.bed"

# The planted ClinVar record: a homozygous-alt truth SNV inside SNAP25, so the
# sample carries it at any depth and step 06 has a hit with a known gene.
read -r PLANT_CHROM PLANT_POS PLANT_REF PLANT_ALT < <(
  bcf view -H -v snps -i 'GT="AA"' -r "$PLANT_WINDOW" /w/.work/truth.vcf.gz \
    | awk 'NR == 1 {print $1, $2, $4, $5}')
if [ -z "${PLANT_POS:-}" ]; then
  echo "ERROR: no homozygous-alt truth SNV in ${PLANT_WINDOW}" >&2
  exit 1
fi
printf 'chrom\tpos\tref\talt\tgene\tclinvar_id\n%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$PLANT_CHROM" "$PLANT_POS" "$PLANT_REF" "$PLANT_ALT" "${PLANT_GENE%%:*}" "$PLANT_ID" \
  > "${OUT}/planted.tsv"
echo "  planted: ${PLANT_CHROM}:${PLANT_POS} ${PLANT_REF}>${PLANT_ALT} (${PLANT_GENE%%:*})"

# --- ClinVar subset, built the way scripts/setup.sh builds the full files ----
echo "[5/8] ClinVar subset"
fetch "$CLINVAR_URL" "${WORK}/clinvar_full.vcf.gz"
fetch "${CLINVAR_URL}.tbi" "${WORK}/clinvar_full.vcf.gz.tbi"
CLINVAR_DATE=$(bcf view -h /w/.work/clinvar_full.vcf.gz | awk -F= '/^##fileDate=/ {print $2; exit}')
NOCHR_REGIONS=()
for r in "${REGIONS[@]}"; do
  r="${r#chr}"
  [ "$r" = "M" ] && r="MT"
  NOCHR_REGIONS+=("$r")
done
NOCHR_LIST=$(IFS=,; echo "${NOCHR_REGIONS[*]}")
PLANT_LINE=$(printf '%s\t%s\t%s\t%s\t%s\t.\t.\t%s' "${PLANT_CHROM#chr}" "$PLANT_POS" "$PLANT_ID" \
  "$PLANT_REF" "$PLANT_ALT" \
  "ALLELEID=${PLANT_ID};CLNDN=Synthetic_e2e_fixture_record;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic;CLNVC=single_nucleotide_variant;GENEINFO=${PLANT_GENE};ORIGIN=1")
{
  bcf view -h /w/.work/clinvar_full.vcf.gz
  { bcf view -H -r "$NOCHR_LIST" /w/.work/clinvar_full.vcf.gz; printf '%s\n' "$PLANT_LINE"; } \
    | sort -t$'\t' -k1,1V -k2,2n
} | bcf view -Oz -o /w/clinvar.vcf.gz -
bcf index -t /w/clinvar.vcf.gz
rm -f "${WORK}/clinvar_full.vcf.gz" "${WORK}/clinvar_full.vcf.gz.tbi"
{
  for c in $(seq 1 22) X Y; do echo "${c} chr${c}"; done
  echo "MT chrM"
} > "${WORK}/chr_rename.txt"
bcf annotate --rename-chrs /w/.work/chr_rename.txt /w/clinvar.vcf.gz -Oz -o /w/clinvar_chr.vcf.gz
bcf index -t /w/clinvar_chr.vcf.gz
bcf view -i 'CLNSIG~"Pathogenic" || CLNSIG~"Likely_pathogenic"' /w/clinvar_chr.vcf.gz \
  -Oz -o /w/clinvar_pathogenic_chr.vcf.gz
bcf index -t /w/clinvar_pathogenic_chr.vcf.gz

# --- VEP-annotated subset (the offline cache does not fit a runner) ----------
echo "[6/8] VEP --database on at most 200 truth variants"
{
  bcf view -h /w/.work/truth.vcf.gz | awk -v s="$SAMPLE" 'BEGIN {OFS = "\t"} /^#CHROM/ {$10 = s} {print}'
  for w in "${VEP_WINDOWS[@]}"; do
    read -r region cap <<< "$w"
    bcf view -H -r "$region" /w/.work/truth.vcf.gz | awk -v n="$cap" 'NR <= n'
  done | sort -t$'\t' -k1,1V -k2,2n -u | awk 'BEGIN {OFS = "\t"} {$7 = "PASS"; print}'
} > "${WORK}/vep_input.vcf"
VEP_INPUT_RECORDS=$(grep -vc '^#' "${WORK}/vep_input.vcf" || true)
in_image "$VEP_IMAGE" vep \
  --input_file /w/.work/vep_input.vcf \
  --output_file "/w/${SAMPLE}_vep.vcf" \
  --vcf --database --assembly GRCh38 --everything \
  --fasta /w/fixture_ref.fa.gz \
  --force_overwrite --no_stats
VEP_RECORDS=$(grep -vc '^#' "${OUT}/${SAMPLE}_vep.vcf" || true)

# A score file in REVEL's layout for every SNV of the VEP subset. The values
# are synthetic (0.010 to 0.990, from the position) and mean nothing; they only
# let steps 30 and the VCFANNO module show that a score track is applied.
{
  printf '#chr\tgrch38_pos\tref\talt\tREVEL\n'
  bcf query -i 'TYPE="snp"' -f '%CHROM\t%POS\t%REF\t%ALT{0}\n' "/w/${SAMPLE}_vep.vcf" \
    | awk 'BEGIN {OFS = "\t"} {printf "%s\t%s\t%s\t%s\t%.3f\n", $1, $2, $3, $4, ($2 % 99 + 1) / 100}' \
    | sort -t$'\t' -k1,1V -k2,2n -u
} | bgzip -c > "${OUT}/revel_synthetic.tsv.gz"
tabix -f -s 1 -b 2 -e 2 "${OUT}/revel_synthetic.tsv.gz"

# --- Ten Manta-style SV records (input for AnnotSV and SV readers) ----------
echo "[7/8] Manta-style SV VCF"
SVS=(
  "chr20 10100000 DEL 2000"
  "chr20 10300000 DUP 5000"
  "chr22 42200000 DEL 3000"
  "chr10 94800000 DEL 1500"
  "chr6 30500000 INV 10000"
  "chr6 32000000 DEL 800"
  "chr5 70100000 DUP 20000"
  "chr12 47850000 INS 300"
  "chrX 73900000 DEL 4000"
  "chrY 2850000 DEL 2500"
)
SV_REF_REGIONS=()
for sv in "${SVS[@]}"; do
  read -r c p _ _ <<< "$sv"
  SV_REF_REGIONS+=("${c}:${p}-${p}")
done
mapfile -t SV_BASES < <(sam faidx /w/fixture_ref.fa.gz "${SV_REF_REGIONS[@]}" | awk '!/^>/ {print toupper($0)}')
{
  printf '##fileformat=VCFv4.1\n##source=GenerateSVCandidates 1.6.0\n'
  awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${OUT}/fixture_ref.fa.gz.fai"
  cat <<'HEADER'
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INS,Description="Insertion">
HEADER
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$SAMPLE"
  i=0
  for sv in "${SVS[@]}"; do
    read -r c p t l <<< "$sv"
    base="${SV_BASES[$i]}"
    i=$((i + 1))
    case "$t" in
      DEL) end=$((p + l)); svlen="-${l}" ;;
      INS) end=$p; svlen="$l" ;;
      *)   end=$((p + l)); svlen="$l" ;;
    esac
    printf '%s\t%s\tManta%s:fixture:%s\t%s\t<%s>\t500\tPASS\tEND=%s;SVTYPE=%s;SVLEN=%s;IMPRECISE;CIPOS=-50,50;CIEND=-50,50\tGT:GQ:PR:SR\t0/1:99:20,12:18,9\n' \
      "$c" "$p" "$t" "$i" "$base" "$t" "$end" "$t" "$svlen"
  done | sort -t$'\t' -k1,1V -k2,2n
} | bcf view -Oz -o "/w/${SAMPLE}_sv_manta_style.vcf.gz" -
bcf index -t "/w/${SAMPLE}_sv_manta_style.vcf.gz"

printf '%s\n' "${REGIONS[@]}" | awk -F'[:-]' 'BEGIN {OFS = "\t"} NF == 3 {print $1, $2 - 1, $3} NF == 1 {print $1, 0, 16569}' \
  > "${OUT}/regions.bed"

# --- Self-checks: the same ones the e2e job repeats after download -----------
echo "[8/8] Checks"
rm -rf "$WORK"
sam quickcheck -v "/w/${SAMPLE}_slice.bam"
for f in "${OUT}"/*.gz; do gzip -t "$f"; done
IDXSTATS=$(sam idxstats "/w/${SAMPLE}_slice.bam")
for c in "${CONTIGS[@]}"; do
  n=$(awk -v c="$c" '$1 == c {print $3}' <<< "$IDXSTATS")
  if [ "${n:-0}" -le 0 ]; then
    echo "ERROR: no reads on ${c} in ${SAMPLE}_slice.bam" >&2
    exit 1
  fi
done
if [ "$VEP_RECORDS" -lt 20 ] || [ "$VEP_RECORDS" -gt 200 ] || ! grep -q '^##INFO=<ID=CSQ' "${OUT}/${SAMPLE}_vep.vcf"; then
  echo "ERROR: ${SAMPLE}_vep.vcf has ${VEP_RECORDS} records (want 20-200) or no CSQ header" >&2
  exit 1
fi

{
  echo "fixture_version: ${VERSION}"
  echo "build_script_sha256: $(sha256sum "${REPO}/scripts/ci/build-fixture.sh" | awk '{print $1}')"
  echo "built_utc: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  echo "sample: ${SAMPLE} (GIAB, public)"
  echo "bam_source: ${BAM_URL}"
  echo "reference_source: ${REF_BASE}/${REF_NAME} (md5 ${REF_MD5})"
  echo "reference_contigs: ${CONTIGS[*]}"
  echo "regions: ${REGIONS[*]}"
  echo "chr20_slice_depth_before: ${FULL_DEPTH}"
  echo "subsample_fraction: ${FRACTION} (seed ${SEED})"
  echo "chr20_slice_depth_after: ${SLICE_DEPTH}"
  echo "read_pairs: ${R1_READS}"
  echo "truth_source: ${TRUTH_BASE}.vcf.gz"
  echo "clinvar_source: ${CLINVAR_URL} (fileDate ${CLINVAR_DATE})"
  echo "planted_clinvar_record: ${PLANT_CHROM}:${PLANT_POS} ${PLANT_REF}>${PLANT_ALT} ID ${PLANT_ID} GENEINFO=${PLANT_GENE} CLNSIG=Pathogenic (synthetic)"
  echo "vep: ${VEP_IMAGE} --database --everything; ${VEP_INPUT_RECORDS} records in, ${VEP_RECORDS} out"
  echo "revel_synthetic: made-up scores, not REVEL"
  echo "images: ${SAMTOOLS_IMAGE} ${BCFTOOLS_IMAGE} ${VEP_IMAGE}"
  echo "idxstats:"
  sed 's/^/  /' <<< "$IDXSTATS"
} > "${OUT}/MANIFEST.txt"

(cd "$OUT" && find . -maxdepth 1 -type f ! -name 'SHA256SUMS*' -printf '%f\n' | LC_ALL=C sort | xargs sha256sum) > "${OUT}/SHA256SUMS.tmp"
mv "${OUT}/SHA256SUMS.tmp" "${OUT}/SHA256SUMS"
TOTAL=$(find "$OUT" -maxdepth 1 -type f -printf '%s\n' | awk '{s += $1} END {print s}')
if [ "$TOTAL" -gt "$MAX_TOTAL_BYTES" ]; then
  echo "ERROR: fixture is ${TOTAL} bytes, over the 1.5 GB budget" >&2
  exit 1
fi

echo "=== Fixture ${VERSION} built: $((TOTAL / 1024 / 1024)) MB ==="
cat "${OUT}/MANIFEST.txt"
