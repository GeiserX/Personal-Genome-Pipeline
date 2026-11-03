#!/usr/bin/env bash
# Mitochondrial analysis — Heteroplasmy detection + variant calling with GATK Mutect2
# Input: Sorted BAM
# Output: Mitochondrial VCF with heteroplasmy fractions
# Runtime: ~15-30 minutes
# Note: Originally planned for MToolBox, but no working Docker image exists.
#       GATK Mutect2 in mitochondrial mode is the standard alternative; the file
#       name is historical.
#
# NuMTs (nuclear copies of mitochondrial DNA) put reads on chrM that come from
# the autosomes. GATK's NuMTFilterTool marks an allele as possible_numt when
# its depth is what such copies could explain at the sample's median autosomal
# coverage. That median is read from step 16b's mosdepth output, or taken from
# AUTOSOMAL_COVERAGE; without either the filter marks nothing.
#
# Last, each call is annotated with its ClinVar record (clinvar/clinvar_chr.vcf.gz,
# which setup.sh installs) and with gnomAD v3.1's chrM sites (GNOMAD_MITO_URL in
# versions.env, which setup.sh installs under gnomad_mito/), into
# mito/<sample>_chrM_annotated.tsv (bin/collect_summary.py mito-annotate). A
# missing file leaves its columns '.' and the table says why.
#
# Then the mtDNA copy number: chrM reads per base over autosomal reads per base
# in the BAM index (samtools idxstats), times 2, into mito/<sample>_mtdna_copies.tsv
# (not assessed under 20 reads per kb of autosome, as on a sliced BAM). And large deletions: eKLIPse
# (EKLIPSE_IMAGE) on the chrM reads, into mito/<sample>_chrM_deletions.tsv. The
# report names the variants and deletions of assets/mito_confirmed.tsv among
# the calls.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/mito"

echo "=== Mitochondrial Analysis (GATK Mutect2): ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"
# Tables of an earlier run go first, so a stage that fails now leaves none behind.
rm -f "${OUTPUT_DIR}/${SAMPLE}_mtdna_copies.tsv" "${OUTPUT_DIR}/${SAMPLE}_chrM_deletions.tsv"


echo "[1/8] Extracting chrM reads..."
run_in --cpus 2 --memory 4g \
  "$SAMTOOLS_IMAGE" \
  bash -c "
    samtools view -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam chrM \
      > /genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam && \
    samtools index /genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam
  "

echo "[2/8] Checking sequence dictionary..."
if [ ! -f "$REF_DICT" ]; then
  echo "  Creating sequence dictionary..."
  # The dictionary goes next to the FASTA, so its directory is writable here.
  run_in --rw "$(dirname "$REF_FASTA")" --cpus 2 --memory 4g \
    "$GATK_IMAGE" \
    gatk CreateSequenceDictionary \
      -R "${REF_FASTA_C}" \
      -O "$(cpath "$REF_DICT")"
else
  echo "  Sequence dictionary already exists, skipping."
fi

echo "[3/8] Running Mutect2 in mitochondrial mode..."
run_in --cpus 4 --memory 8g \
  "$GATK_IMAGE" \
  gatk Mutect2 \
    -R "${REF_FASTA_C}" \
    -I "/genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam" \
    -L chrM \
    --mitochondria-mode \
    --max-mnp-distance 0 \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2.vcf.gz"

echo "[4/8] Filtering variants..."
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  gatk FilterMutectCalls \
    -R "${REF_FASTA_C}" \
    -V "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2.vcf.gz" \
    --mitochondria-mode \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2_filtered.vcf.gz"

echo "[5/8] Marking possible NuMTs..."
# Median autosomal depth: the depth at which half of the chr1-22 bases are
# covered at least that deep, from mosdepth's per-chromosome distribution.
MOSDEPTH_PREFIX="${SAMPLE_DIR}/mosdepth/${SAMPLE}.mosdepth"
if [ -n "${AUTOSOMAL_COVERAGE:-}" ]; then
  echo "  Median autosomal coverage: ${AUTOSOMAL_COVERAGE} (AUTOSOMAL_COVERAGE)"
elif [ -s "${MOSDEPTH_PREFIX}.summary.txt" ] && [ -s "${MOSDEPTH_PREFIX}.global.dist.txt" ]; then
  AUTOSOMAL_COVERAGE=$(awk -F'\t' '
    FNR == NR { if ($1 ~ /^chr[0-9]+$/) len[$1] = $2; next }
    ($1 in len) { at_least[$2] += len[$1] * $3 }
    END {
      for (c in len) total += len[c]
      best = 0
      if (total > 0) for (k in at_least) if (at_least[k] / total >= 0.5 && k + 0 > best) best = k + 0
      print best
    }' "${MOSDEPTH_PREFIX}.summary.txt" "${MOSDEPTH_PREFIX}.global.dist.txt")
  echo "  Median autosomal coverage: ${AUTOSOMAL_COVERAGE} (from ${MOSDEPTH_PREFIX}.global.dist.txt)"
else
  AUTOSOMAL_COVERAGE=0
  echo "  WARNING: no mosdepth output (run scripts/16b-mosdepth.sh first) and no AUTOSOMAL_COVERAGE:"
  echo "  the NuMT filter runs with coverage 0 and marks nothing."
fi
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  gatk NuMTFilterTool \
    -R "${REF_FASTA_C}" \
    -V "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2_filtered.vcf.gz" \
    --autosomal-coverage "$AUTOSOMAL_COVERAGE" \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz"

echo "[6/8] Annotating the calls with ClinVar and gnomAD..."
CLINVAR_VCF="${GENOME_DIR}/clinvar/clinvar_chr.vcf.gz"
GNOMAD_MITO="${GENOME_DIR}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
[ -s "${CLINVAR_VCF}.tbi" ] || echo "  WARNING: no ${CLINVAR_VCF}.tbi (setup.sh installs ClinVar): no ClinVar columns."
[ -s "$GNOMAD_MITO" ] || echo "  WARNING: no ${GNOMAD_MITO} (setup.sh installs it): no gnomAD columns."
if [ ! -s "${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz" ]; then
  echo "ERROR: NuMTFilterTool wrote no ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz" >&2
  exit 1
fi
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/collect_summary.py mito-annotate \
    --vcf "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz" \
    --clinvar "$(cpath "$CLINVAR_VCF")" \
    --gnomad "$(cpath "$GNOMAD_MITO")" \
    --out "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_annotated.tsv"

echo "[7/8] mtDNA copy number (chrM against autosomal reads, from the BAM index)..."
run_in --cpus 1 --memory 1g "$SAMTOOLS_IMAGE" \
  bash -c "samtools idxstats /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam > /genome/${SAMPLE}/mito/${SAMPLE}_idxstats.txt"
run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/collect_summary.py mito-copies \
    --idxstats "/genome/${SAMPLE}/mito/${SAMPLE}_idxstats.txt" \
    --out "/genome/${SAMPLE}/mito/${SAMPLE}_mtdna_copies.tsv"

echo "[8/8] Large deletions (eKLIPse, soft-clipped chrM reads)..."
# eKLIPse writes into a new eKLIPse_<id> folder under -out and -tmp on each run.
EK_DIR="${OUTPUT_DIR}/eklipse"
rm -rf "$EK_DIR"
mkdir -p "${EK_DIR}/out" "${EK_DIR}/tmp"
printf '%s\t%s\n' "/genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam" "$SAMPLE" > "${EK_DIR}/input.tsv"
# GRCh38's chrM is the rCRS: the reference is the GenBank file the image ships beside eKLIPse.py.
# shellcheck disable=SC2016  # expanded inside the container
run_in --cpus 2 --memory 4g "$EKLIPSE_IMAGE" bash -c '
  EK=$(command -v eKLIPse.py)
  python3 "$EK" -in "$1/input.tsv" -ref "$(dirname "$EK")/data/NC_012920.1.gb" -out "$1/out" -tmp "$1/tmp" -thread 2
' _ "/genome/${SAMPLE}/mito/eklipse"
# eKLIPse stops with exit 1 on an error (set -e stops here then). Should it
# still leave no table, there is none (the old one went at the start) and the
# report says the deletions were not looked for.
CSV=$(find "${EK_DIR}/out" -name eKLIPse_deletions.csv | head -n 1)
if [ -z "$CSV" ]; then
  echo "  WARNING: eKLIPse exited 0 but wrote no eKLIPse_deletions.csv under ${EK_DIR}/out: no deletion table."
else
  run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
    python3 /pgp-bin/collect_summary.py mito-deletions \
      --eklipse "/genome/${SAMPLE}/mito/eklipse/out/${CSV#"${EK_DIR}/out/"}" \
      --out "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_deletions.tsv"
fi

echo "=== Mitochondrial analysis complete ==="
echo "Raw calls: ${OUTPUT_DIR}/${SAMPLE}_chrM_mutect2.vcf.gz"
echo "Filtered (FilterMutectCalls, then possible_numt): ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz"
echo "ClinVar and gnomAD per allele: ${OUTPUT_DIR}/${SAMPLE}_chrM_annotated.tsv"
echo "mtDNA copy number: ${OUTPUT_DIR}/${SAMPLE}_mtdna_copies.tsv"
echo "Large deletions (eKLIPse): ${OUTPUT_DIR}/${SAMPLE}_chrM_deletions.tsv (eKLIPse's own files and plot: ${EK_DIR}/out)"
echo ""
echo "PASS variants below AF 0.95 (possible heteroplasmy; calls near the ends of the"
echo "control region, chrM:1-500 and 16000-16569, are less reliable, and length"
echo "variants in the D-loop poly-C tracts m.303-315 and m.16184-16193 are common and"
echo "usually technical, so the reports count them apart):"
echo "  bcftools query -f '%POS\t%REF\t%ALT\t[%AF]\n' -i 'FILTER=\"PASS\"' ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz | awk '\$4 < 0.95'"
