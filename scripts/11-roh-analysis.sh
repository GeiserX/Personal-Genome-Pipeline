#!/usr/bin/env bash
# ROH Analysis — Runs of homozygosity (consanguinity screening)
# Detects long stretches of homozygous DNA that indicate shared ancestry.
# The summary leaves out chrX/chrY and the part of each segment in a centromere
# or heterochromatin band (MASKED_BANDS below), where few markers map.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTPUT="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh.txt"

echo "=== ROH Analysis: ${SAMPLE} ==="

for f in "$VCF" "${VCF}.tbi"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# Auto-detect chip data: if FORMAT/PL is absent, use -G30 (genotype-only mode)
HAS_PL=$(run_in \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -h "/data/${SAMPLE}.vcf.gz" | grep -c '##FORMAT=<ID=PL' || true)

ROH_FLAGS=(--AF-dflt 0.4)
if [ "${HAS_PL}" -eq 0 ]; then
  echo "  No FORMAT/PL tag detected (chip data) — using -G30 genotype-only mode"
  ROH_FLAGS+=(-G30)
fi

# Only PASS calls (and records with no filter, as chip VCFs have): DeepVariant's
# RefCall and other filtered records are not genotypes to count on. The same
# input as the Nextflow ROH module. pipefail: a failed view fails the step.
# (The flags and the validated sample name are plain words, safe inline.)
run_in \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  "${BCFTOOLS_IMAGE}" \
  bash -euo pipefail -c "bcftools view -f PASS,. -Ou /data/${SAMPLE}.vcf.gz | bcftools roh ${ROH_FLAGS[*]} -o /data/${SAMPLE}_roh.txt -"

# GRCh38's centromere and heterochromatin bands: the acen, gvar and stalk
# rows of UCSC's cytoBand.hg38.txt, neighbouring rows merged (0-based start,
# exclusive end). bin/collect_summary.py (MASKED_BANDS, which says how to
# regenerate it) and modules/local/roh/main.nf hold the same line.
MASKED_BANDS="chr1:121700000-143200000 chr2:91800000-96000000 chr3:87800000-98600000 chr4:48200000-51800000 chr5:46100000-51400000 chr6:58500000-62600000 chr7:58100000-62100000 chr8:43200000-47200000 chr9:42200000-61500000 chr10:38000000-41600000 chr11:51000000-55800000 chr12:33200000-37800000 chr13:0-18900000 chr14:0-18200000 chr15:0-20500000 chr16:35300000-47000000 chr17:22700000-27400000 chr18:15400000-21500000 chr19:19900000-31900000 chr20:25700000-30400000 chr21:0-13000000 chr22:0-17400000 chrX:58100000-63800000 chrY:10300000-10600000 chrY:26600000-57227415"

# Summary: autosomal segments with 5 Mb or more outside those bands, the
# threshold of docs/11-roh-analysis.md and of the reports; the Nextflow module
# writes the same file.
SUMMARY="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt"
{
  echo "# ROH Summary for ${SAMPLE}"
  echo "# Segments with >=5MB on autosomes (excludes chrX/chrY) outside centromere and heterochromatin bands; length_bp is that part, masked_bp the part in the bands"
  printf 'chrom\tstart\tend\tlength_bp\tlength_mb\tmasked_bp\n'
  awk -v bands="$MASKED_BANDS" '
    BEGIN { n = split(bands, b, " "); for (i = 1; i <= n; i++) { split(b[i], x, /[:-]/); c[i] = x[1]; s[i] = x[2]; e[i] = x[3] } }
    $1 == "RG" && $3 !~ /chrX|chrY/ {
      m = 0
      for (i = 1; i <= n; i++) if (c[i] == $3) { o = ($5 < e[i] ? $5 : e[i]) - ($4 - 1 > s[i] ? $4 - 1 : s[i]); if (o > 0) m += o }
      if (m > $6) m = $6
      l = $6 - m
      if (l >= 5000000) printf "%s\t%s\t%s\t%d\t%.1f\t%d\n", $3, $4, $5, l, l / 1e6, m
    }' "$OUTPUT"
} > "${SUMMARY}.tmp"
mv -f "${SUMMARY}.tmp" "$SUMMARY"

echo "=== ROH complete ==="
echo "Results: ${OUTPUT}"
echo "Summary: ${SUMMARY}"
echo ""
echo "Autosomal ROH of 5 Mb or more outside centromere and heterochromatin bands (potential consanguinity signal):"
awk -F'\t' 'NR > 3 {printf "%s:%s-%s  %sMB\n", $1, $2, $3, $5}' "$SUMMARY"
echo ""
echo "NOTE: no population allele frequencies were used (--AF-dflt 0.4), so read the total and the short"
echo "      segments as rough; the long segments are the signal."
