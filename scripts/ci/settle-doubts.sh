#!/usr/bin/env bash
# settle-doubts.sh — run the tools to answer questions a code review could not
# answer by reading. Prints one table row per question: the command and what
# was observed. Not a gate: a question whose command fails is answered by that
# failure, and the script exits 0 unless it cannot set itself up.
#
# Usage: scripts/ci/settle-doubts.sh <work_dir>
# Needs Docker, the gh CLI (GH_TOKEN) and about 30 GB of disk. Runs in the E2E
# workflow (job settle-doubts, on dispatch or on a push that changes this file).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
# shellcheck source=../../versions.env
. "${REPO}/versions.env"

W_ARG=${1:?Usage: $0 <work_dir>}
mkdir -p "$W_ARG"
W="$(cd "$W_ARG" && pwd)"
SUMMARY=${GITHUB_STEP_SUMMARY:-/dev/null}
GH_REPO=${GITHUB_REPOSITORY:-GeiserX/Personal-Genome-Pipeline}
TAG=$(tr -d '[:space:]' < "${REPO}/tests/fixtures/VERSION")
SLICE=chr20:10000000-10500000
HAPPY_IMAGE="jmcdani20/hap.py:v0.3.12"            # the image benchmark-variants.sh uses
BWAMEM2_IMAGE="quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5"   # as in 02a-alignment-bwamem2.sh
CLAIR3_IMAGE="hkubal/clair3:v2.0.2"                 # as in 03e-clair3.sh

export PATH="${REPO}/tests/e2e/bin:${PATH}"   # clamps --cpus to this machine
export THREADS=4
export SAMPLE=HG002
FX="${W}/fixture"
G="${W}/genome"            # chr20 + chrM reference, for the alignment questions
export GENOME_DIR="$G"
LOGS="${W}/logs"
mkdir -p "$FX" "$LOGS"

ROWS=()
# clean <text>: one table cell (no newlines, no pipes).
clean() { tr '\n|' ' /' <<< "$1" | sed 's/  */ /g; s/ $//'; }
# row <n> <question> <command> <answer>
row() {
  ROWS+=("| $1 | $(clean "$2") | \`$(clean "$3")\` | $(clean "$4") |")
  echo "== row $1: $4"
}
# last_lines <file> <n>: the tail of a log, for an answer cell.
last_lines() { tail -n "$2" "$1" 2>/dev/null | tr -s '\n' ' '; }
in_g() { local image=$1; shift; docker run --rm -i -v "${G}:/genome" -w /genome "$image" "$@"; }

# --- Setup: fixture, chr20+chrM reference, reads ------------------------------
gh release download "$TAG" -R "$GH_REPO" -D "$FX" --clobber
(cd "$FX" && sha256sum -c --quiet SHA256SUMS)
mkdir -p "${G}/reference" "${G}/${SAMPLE}/fastq"
in_fx() { docker run --rm -i -v "${FX}:/f" -v "${G}:/genome" -w /f "$SAMTOOLS_IMAGE" "$@"; }
in_fx samtools faidx -o /genome/reference/Homo_sapiens_assembly38.fasta /f/fixture_ref.fa.gz chr20 chrM
in_g "$SAMTOOLS_IMAGE" samtools faidx reference/Homo_sapiens_assembly38.fasta
in_g "$SAMTOOLS_IMAGE" samtools dict -o reference/Homo_sapiens_assembly38.dict reference/Homo_sapiens_assembly38.fasta
# Reads the GIAB alignment placed on chr20 or chrM, as pairs.
in_fx samtools view -u /f/HG002_slice.bam chr20 chrM \
  | in_fx samtools collate -u -O - /tmp/c \
  | in_fx samtools fastq -n -1 "/genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz" \
      -2 "/genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz" -0 /dev/null -s /dev/null -
cp "${FX}/HG002_truth_chr20.vcf.gz" "${FX}/HG002_truth_chr20.vcf.gz.tbi" "${FX}/HG002_truth_chr20.bed" "${G}/reference/"

# From here on a failing command is an answer, not a reason to stop.
set +e

# Step 02 as shipped (it builds the minimap2 index on first use).
"${REPO}/scripts/02-alignment.sh" "$SAMPLE" > "${LOGS}/02.log" 2>&1
BAM="${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
HAS_RG=$(in_g "$SAMTOOLS_IMAGE" samtools view -H "$BAM" | grep -c '^@RG' || true)

# --- 1. GATK and DeepVariant on a step-02 BAM without @RG ----------------------
# If step 02 already writes @RG, strip it so the question is still the one asked.
NORG=HG002norg
mkdir -p "${G}/${NORG}/aligned"
in_g "$SAMTOOLS_IMAGE" bash -c "samtools view -h ${BAM} | grep -v '^@RG' \
  | sed 's/\tRG:Z:[^\t]*//' | samtools view -b -o ${NORG}/aligned/${NORG}_sorted.bam - \
  && samtools index ${NORG}/aligned/${NORG}_sorted.bam"
"${REPO}/scripts/20-mtoolbox.sh" "$NORG" > "${LOGS}/q1_step20.log" 2>&1; RC20=$?
INTERVALS="$SLICE" "${REPO}/scripts/03a-gatk-haplotypecaller.sh" "$NORG" > "${LOGS}/q1_step03a.log" 2>&1; RC03A=$?
ERR20=$(grep -m1 -iE 'read group|USER ERROR|sample list|samples cannot|Exception' "${LOGS}/q1_step20.log")
ERR03A=$(grep -m1 -iE 'read group|USER ERROR|sample list|samples cannot|Exception' "${LOGS}/q1_step03a.log")

run_dv() {   # run_dv <sample> <out name>: step 03's command on the chr20 slice only
  docker run --rm --cpus 4 --memory 14g -v "${G}:/genome" "$DEEPVARIANT_IMAGE" \
    /opt/deepvariant/bin/run_deepvariant --model_type=WGS \
      --ref=/genome/reference/Homo_sapiens_assembly38.fasta \
      --reads="/genome/$1/aligned/$1_sorted.bam" \
      --regions="$SLICE" \
      --output_vcf="/genome/$1/$2.vcf.gz" --num_shards=4
}
run_dv "$NORG" dv_norg > "${LOGS}/q1_dv.log" 2>&1; RCDV=$?
DV_NAME=$(in_g "$BCFTOOLS_IMAGE" bcftools query -l "${NORG}/dv_norg.vcf.gz" 2>/dev/null || echo "no VCF (exit ${RCDV})")
row 1 "Do step 20 (Mutect2) and step 03a (HaplotypeCaller, ${SLICE}) fail on a step-02 BAM with no @RG, and what sample name does DeepVariant write?" \
  "scripts/20-mtoolbox.sh; INTERVALS=${SLICE} scripts/03a-gatk-haplotypecaller.sh; run_deepvariant without --sample_name" \
  "step 02 as shipped writes @RG: $([ "$HAS_RG" -gt 0 ] && echo yes || echo no). On the no-@RG BAM: step 20 exit ${RC20} (${ERR20:-no read-group error printed}); step 03a exit ${RC03A} (${ERR03A:-no read-group error printed}); DeepVariant sample name: ${DV_NAME}"

# --- 2. TelomereHunter default banding -----------------------------------------
docker run --rm "$TELOMEREHUNTER_IMAGE" telomerehunter --help > "${LOGS}/q2.log" 2>&1; RC=$?
BAND=$(grep -iE -A3 'band' "${LOGS}/q2.log" | tr -s ' \n' ' ' | cut -c1-600)
row 2 "telomerehunter --help in the pinned image: is the default banding hg19?" \
  "docker run TELOMEREHUNTER_IMAGE telomerehunter --help" \
  "exit ${RC}; banding lines: ${BAND:-none mention banding}"

# --- 3. Clair3 model directories -----------------------------------------------
docker pull -q "$CLAIR3_IMAGE" >/dev/null 2>&1
docker run --rm "$CLAIR3_IMAGE" bash -c 'for d in /opt/models/r1041_e82_400bps_sup_v500 /opt/models/hifi_revio; do
  [ -d "$d" ] && echo "$d: present" || echo "$d: MISSING"; done; echo "models in /opt/models: $(ls /opt/models 2>&1 | tr "\n" " ")"' \
  > "${LOGS}/q3.log" 2>&1
row 3 "Do the two model directories hardcoded in scripts/03e-clair3.sh exist in the pinned image?" \
  "docker run ${CLAIR3_IMAGE} ls -d /opt/models/r1041_e82_400bps_sup_v500 /opt/models/hifi_revio" \
  "$(grep -E 'present|MISSING|models in' "${LOGS}/q3.log" | cut -c1-700)"

# --- 4. TIDDIT with only a BWA-MEM2 index ---------------------------------------
in_g "$BWAMEM2_IMAGE" bwa-mem2 index reference/Homo_sapiens_assembly38.fasta > "${LOGS}/q4_index.log" 2>&1
IDX_FILES=$(cd "${G}/reference" && ls Homo_sapiens_assembly38.fasta.* | tr '\n' ' ')
"${REPO}/scripts/04a-tiddit.sh" "$SAMPLE" > "${LOGS}/q4.log" 2>&1; RC=$?
ASM=$(grep -m1 -E 'BWA index detected|No BWA index' "${LOGS}/q4.log")
ERR=$(grep -E '^[A-Za-z]+Error|Exception' "${LOGS}/q4.log" | tail -n 1)
# Control: the same run with the BWA-MEM2 index moved away (--skip_assembly).
mkdir -p "${G}/bwamem2-aside"
mv "${G}/reference/Homo_sapiens_assembly38.fasta.bwt.2bit.64" "${G}/bwamem2-aside/"
rm -rf "${G}/${SAMPLE}/sv_tiddit"
"${REPO}/scripts/04a-tiddit.sh" "$SAMPLE" > "${LOGS}/q4_control.log" 2>&1; RCC=$?
SVC=$(in_g "$BCFTOOLS_IMAGE" bcftools view -H "${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz" 2>/dev/null | wc -l | tr -d ' ')
row 4 "What does TIDDIT do when only the BWA-MEM2 index exists: crash, or run without assembly?" \
  "bwa-mem2 index ref.fasta; scripts/04a-tiddit.sh (chr20+chrM reference); then again without the .bwt.2bit.64 file" \
  "index files: ${IDX_FILES}; script says: ${ASM:-nothing}; exit ${RC}: ${ERR:-no Python error}. Control without the BWA-MEM2 index (--skip_assembly): exit ${RCC}, ${SVC:-0} SV records"

# --- 5. bcftools convert --tsv2vcf on a five-column AncestryDNA-style file -------
mkdir -p "${G}/chip"
in_g "$BCFTOOLS_IMAGE" bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT{0}\t[%GT]\n' -i 'TYPE="snp"' \
  reference/HG002_truth_chr20.vcf.gz 2>/dev/null | awk 'NR <= 20' \
  | awk 'BEGIN {OFS = "\t"; print "#AncestryDNA raw data download (synthetic, from the GIAB truth)"; print "rsid\tchromosome\tposition\tallele1\tallele2"}
         {a1 = ($6 ~ /^0/) ? $4 : $5; a2 = ($6 ~ /1$/) ? $5 : $4; print "rs" NR, $2, $3, a1, a2}' \
  > "${G}/chip/ancestry_5col.txt"
in_g "$BCFTOOLS_IMAGE" bcftools convert --tsv2vcf chip/ancestry_5col.txt -f reference/Homo_sapiens_assembly38.fasta \
  -s "$SAMPLE" -c ID,CHROM,POS,AA -Oz -o chip/out.vcf.gz > "${LOGS}/q5.log" 2>&1; RC=$?
GTS=$(in_g "$BCFTOOLS_IMAGE" bcftools query -f '[%GT] ' chip/out.vcf.gz 2>/dev/null | cut -c1-120)
NREC=$(in_g "$BCFTOOLS_IMAGE" bcftools view -H chip/out.vcf.gz 2>/dev/null | wc -l | tr -d ' ')
row 5 "What does bcftools convert --tsv2vcf -c ID,CHROM,POS,AA do with a five-column AncestryDNA-style file?" \
  "bcftools convert --tsv2vcf ancestry_5col.txt -f ref -s ${SAMPLE} -c ID,CHROM,POS,AA" \
  "exit ${RC}; ${NREC:-0} of 20 rows became records; genotypes: ${GTS:-none}; messages: $(last_lines "${LOGS}/q5.log" 3 | cut -c1-300)"

# --- 6. AnnotSV without its annotation data ----------------------------------------
mkdir -p "${G}/${SAMPLE}/sv"
cp "${FX}/HG002_sv_manta_style.vcf.gz" "${FX}/HG002_sv_manta_style.vcf.gz.tbi" "${G}/${SAMPLE}/sv/"
docker run --rm --user root -v "${G}:/genome" "$ANNOTSV_IMAGE" AnnotSV \
  -SVinputFile "/genome/${SAMPLE}/sv/HG002_sv_manta_style.vcf.gz" \
  -outputFile "/genome/${SAMPLE}/sv/annotsv_direct.tsv" -genomeBuild GRCh38 -annotationMode both \
  > "${LOGS}/q6.log" 2>&1; RC=$?
SV_VCF="${G}/${SAMPLE}/sv/HG002_sv_manta_style.vcf.gz" "${REPO}/scripts/05-annotsv.sh" "$SAMPLE" > "${LOGS}/q6_script.log" 2>&1; RCS=$?
row 6 "AnnotSV exit code without its annotation data" \
  "AnnotSV -SVinputFile <ten Manta-style SVs> -genomeBuild GRCh38 -annotationMode both; SV_VCF=... scripts/05-annotsv.sh" \
  "AnnotSV exit ${RC}: $(grep -m2 -iE 'error|not found|annotation' "${LOGS}/q6.log" | tr '\n' ' ' | cut -c1-300); scripts/05-annotsv.sh exit ${RCS}"

# --- 7. VEP --everything --offline without --fasta -----------------------------------
row 7 "First lines VEP prints for --everything --offline without --fasta" \
  "vep --offline --everything (no --fasta)" \
  "skipped: needs the 26 GB offline cache, which does not fit a GitHub runner"

# --- 8. minimap2 killed mid-run --------------------------------------------------------
KILL=HG002kill
mkdir -p "${G}/${KILL}/fastq"
# minimap2 writes nothing until its first batch (500 Mbp) is mapped, so the
# input is the fixture reads four times over, and the kill comes right after
# minimap2 logs its first mapped batch: samtools sort has data by then.
for r in R1 R2; do
  for _ in 1 2 3 4; do cat "${FX}/HG002_${r}.fastq.gz"; done > "${G}/${KILL}/fastq/${KILL}_${r}.fastq.gz"
done
INPUT_READS=$(( $(gzip -dc "${G}/${KILL}/fastq/${KILL}_R1.fastq.gz" | wc -l) / 2 ))
"${REPO}/scripts/02-alignment.sh" "$KILL" > "${LOGS}/q8.log" 2>&1 &
PID=$!
CID=""
for _ in $(seq 1 120); do
  CID=$(docker ps -q --filter "ancestor=${MINIMAP2_IMAGE}" | awk 'NR == 1')
  [ -n "$CID" ] && break
  sleep 1
done
KILLED=no
WHEN="never started"
if [ -n "$CID" ]; then
  for t in $(seq 1 600); do
    if docker logs "$CID" 2>&1 | grep -q 'mapped [0-9]* sequences'; then WHEN="after the first mapped batch (${t} s)"; break; fi
    docker ps -q --no-trunc | grep -q "^${CID}" || { WHEN="it had already exited"; break; }
    sleep 1
  done
  docker kill "$CID" >/dev/null 2>&1 && KILLED=yes
fi
wait "$PID"; RC=$?
KBAM="${KILL}/aligned/${KILL}_sorted.bam"
if [ -f "${G}/${KBAM}" ]; then
  QC=$(in_g "$SAMTOOLS_IMAGE" samtools quickcheck -v "$KBAM" 2>&1 && echo "passes quickcheck" || echo "fails quickcheck")
  NREADS=$(in_g "$SAMTOOLS_IMAGE" samtools view -c "$KBAM" 2>/dev/null || echo "unreadable")
  LEFT="a BAM is left ($(stat -c %s "${G}/${KBAM}") bytes, ${QC}, ${NREADS} records; the input has ${INPUT_READS} reads)"
else
  LEFT="no BAM is left"
fi
row 8 "Kill the minimap2 container mid-run: is the BAM left behind valid or partial?" \
  "scripts/02-alignment.sh & ; docker kill <minimap2 container> once it logs a mapped batch" \
  "killed: ${KILLED}, ${WHEN}; script exit ${RC}; ${LEFT}; index left: $([ -f "${G}/${KBAM}.bai" ] && echo yes || echo no); $(grep -m1 -iE 'samtools sort|truncated|EOF' "${LOGS}/q8.log")"

# --- 9. Images that run as non-root, and whether they can write to /genome ------------
mkdir -p "${W}/mount-test"
chmod 755 "${W}/mount-test"
NONROOT=()
CHECKED=0
while IFS='=' read -r var val; do
  img=$(sed -E 's/^"([^"]*)".*/\1/' <<< "$val")
  CHECKED=$((CHECKED + 1))
  user=$(docker buildx imagetools inspect "$img" --format '{{json .Image}}' 2>/dev/null | python3 -c '
import json, sys
d = json.load(sys.stdin)
if "config" not in d:
    d = d.get("linux/amd64") or next(iter(d.values()))
print((d.get("config") or {}).get("User") or "")' 2>/dev/null || echo "?")
  case "$user" in
    ""|root|0|0:0) ;;
    *)
      w=$(docker run --rm -v "${W}/mount-test:/genome" --entrypoint sh "$img" -c \
            'touch /genome/.probe 2>/dev/null && echo writable || echo "not writable"' 2>/dev/null || echo "could not run sh")
      NONROOT+=("${var} (user ${user}: /genome ${w})")
      ;;
  esac
done < <(grep -E '^[A-Z0-9_]+_IMAGE=' "${REPO}/versions.env")
row 9 "Which images run as a non-root user, and can they write to a /genome mount owned by the host user?" \
  "docker buildx imagetools inspect <image> (config User); docker run -v dir:/genome <image> touch /genome/.probe" \
  "${CHECKED} images in versions.env checked; non-root: ${NONROOT[*]:-none}. Every other image runs as root by default."

# --- 10. T1K coordinate file built from the reference ------------------------------------
T1K=${W}/t1k
mkdir -p "${T1K}/ref"
docker run --rm -v "${FX}:/f" -v "${T1K}:/t" "$SAMTOOLS_IMAGE" samtools faidx -o /t/ref/chr6.fasta /f/fixture_ref.fa.gz chr6
docker run --rm -v "${T1K}:/t" -w /t "$T1K_IMAGE" t1k-build.pl -o /t/hlaidx --download IPD-IMGT/HLA > "${LOGS}/q10_download.log" 2>&1
docker run --rm -v "${T1K}:/t" -w /t "$T1K_IMAGE" t1k-build.pl -d /t/hlaidx/hla.dat -g /t/ref/chr6.fasta -o /t/hlaidx_grch38 \
  > "${LOGS}/q10_build.log" 2>&1; RC=$?
COORD="${T1K}/hlaidx_grch38/_dna_coord.fa"
if [ -f "$COORD" ]; then
  TOTAL=$(grep -c '^>' "$COORD" || true)
  NEG=$(grep -c ' -1 -1 ' "$COORD" || true)
  ANS="exit ${RC}; ${NEG} of ${TOTAL} entries have ' -1 -1 ' coordinates; e.g. $(grep -m2 ' -1 -1 ' "$COORD" | tr '\n' ' ' | cut -c1-200). Step 08 passes the FASTA to -g; T1K's README passes a GENCODE GTF there"
else
  ANS="exit ${RC}; no coordinate file; $(last_lines "${LOGS}/q10_build.log" 3 | cut -c1-300)"
fi
row 10 "Does T1K's coordinate file built from the reference contain ' -1 -1 ' rows?" \
  "t1k-build.pl -d hla.dat -g chr6.fasta (as step 08 does, with chr6 only: every HLA gene is on chr6)" "$ANS"

# --- 11. fastp trimming A/B against the GIAB truth ----------------------------------------
TRIM=HG002trim
mkdir -p "${G}/${TRIM}/fastq"
cp "${G}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz" "${G}/${TRIM}/fastq/${TRIM}_R1.fastq.gz"
cp "${G}/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz" "${G}/${TRIM}/fastq/${TRIM}_R2.fastq.gz"
"${REPO}/scripts/01b-fastp-qc.sh" "$TRIM" > "${LOGS}/q11_fastp.log" 2>&1
"${REPO}/scripts/02-alignment.sh" "$TRIM" > "${LOGS}/q11_align.log" 2>&1
run_dv "$SAMPLE" dv_raw > "${LOGS}/q11_dv_raw.log" 2>&1
run_dv "$TRIM" dv_trim > "${LOGS}/q11_dv_trim.log" 2>&1
happy() {   # happy <sample> <vcf name>: SNP and INDEL recall/precision/F1 on the slice
  docker run --rm --user root -v "${G}:/genome" "$HAPPY_IMAGE" /opt/hap.py/bin/hap.py \
    /genome/reference/HG002_truth_chr20.vcf.gz "/genome/$1/$2.vcf.gz" \
    -r /genome/reference/Homo_sapiens_assembly38.fasta -f /genome/reference/HG002_truth_chr20.bed \
    -o "/genome/$1/happy_$2" --engine=vcfeval > "${LOGS}/q11_happy_$2.log" 2>&1 || true
  awk -F',' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
             ($1 == "SNP" || $1 == "INDEL") && $2 == "PASS" {
               printf "%s recall %.4f precision %.4f F1 %.4f (TP %s FN %s FP %s); ", $1,
                 $c["METRIC.Recall"], $c["METRIC.Precision"], $c["METRIC.F1_Score"], $c["TRUTH.TP"], $c["TRUTH.FN"], $c["QUERY.FP"]}' \
    "${G}/$1/happy_$2.summary.csv" 2>/dev/null || true
}
RAW=$(happy "$SAMPLE" dv_raw)
TRM=$(happy "$TRIM" dv_trim)
HAPPY_COLS=$(head -n 1 "${G}/${SAMPLE}/happy_dv_raw.summary.csv" 2>/dev/null | cut -d, -f1-14)
row 11 "fastp A/B: DeepVariant on ${SLICE} with and without step 01b, against GIAB v4.2.1 (hap.py)" \
  "01b-fastp-qc.sh; 02-alignment.sh; run_deepvariant --regions ${SLICE}; hap.py --engine=vcfeval -f truth.bed" \
  "without trimming: ${RAW:-hap.py produced no summary}; with trimming: ${TRM:-hap.py produced no summary}. summary.csv columns 1-14: ${HAPPY_COLS:-none}"

# --- Table -------------------------------------------------------------------------------
{
  echo "### Settle doubts (fixture ${TAG}, $(git -C "$REPO" rev-parse --short HEAD))"
  echo
  echo "| # | Question | Command | Observed |"
  echo "|---|---|---|---|"
  printf '%s\n' "${ROWS[@]}"
  echo
  echo "Full logs: the settle-doubts step output."
} | tee -a "$SUMMARY"

for f in "${LOGS}"/*.log; do
  echo "::group::$(basename "$f")"
  tail -n 60 "$f"
  echo "::endgroup::"
done
