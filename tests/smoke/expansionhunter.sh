# shellcheck shell=sh
# EXPANSIONHUNTER_IMAGE row. Step 09 reads the catalog the image ships; it is
# copied out so the row can check it. The fixture's slices hold none of its
# disease loci, so the genotype comes from a one-locus catalog for a CA repeat
# found in the mini reference (image-smoke.sh writes it).
set -e
CATALOG=/usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json
cp "$CATALOG" bundled_catalog.json
cat /in/eh_catalog.json
ExpansionHunter --reads /in/mini.bam --reference /in/mini.fa --variant-catalog /in/eh_catalog.json \
  --output-prefix eh --threads 4 --sex male
