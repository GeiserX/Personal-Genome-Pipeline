#!/usr/bin/env bash
# Scratch helper (deleted before merge): real AnnotSV run with the 3.5
# human annotations (GRCh38 part only, to fit the runner disk).
set -euo pipefail

D=${1:?}; A=${2:?}
mkdir -p "$A"
echo "== AnnotSV human annotations 3.5 (5.3 GB download, GRCh37 files skipped)"
curl -fsSL https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz \
  | tar -xz -C "$A" --exclude='*GRCh37*'
du -sh "$A"
find "$A" -maxdepth 3 -type d | sort | head -40
df -h /mnt /

cd "${GITHUB_WORKSPACE:?}"
nextflow run main.nf -profile docker -c "$D/ci.config" \
  --input "$D/samplesheet.csv" --reference "$D/ref.fa" --outdir results_annotsv \
  --tools manta,duphold,annotsv \
  --annotsv_annotations "$A"

T=results_annotsv/HG002/annotsv/HG002_sv_annotated.tsv
test -s "$T"
echo "== AnnotSV output: $(wc -l < "$T") lines"
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "ACMG_class") c = i; if (!c) { print "no ACMG_class column"; exit 1 } next } { n[$c]++ } END { for (k in n) print "ACMG_class=" k ": " n[k] }' "$T"
filled=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "ACMG_class") c = i; next } $c != "" && $c != "NA" && $c != "." { n++ } END { print n + 0 }' "$T")
echo "rows with ACMG_class filled: ${filled}"
test "$filled" -gt 0
