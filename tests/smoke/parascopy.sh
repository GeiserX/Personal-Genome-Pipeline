# shellcheck shell=sh
# PARASCOPY_IMAGE row: step 35 on the fixture. Parascopy's GRCh38 homology
# table and 1000 Genomes models are downloaded from Zenodo as
# `setup.sh --parascopy-data` does (PARASCOPY_DATA_VERSION, md5 checked). The
# slice BAM holds reads only in the fixture's regions, so the background depth
# comes from 100 bp windows over chr20:10.05-10.45 Mb (PARASCOPY_DEPTH_BED in
# step 35) instead of Parascopy's genome-wide set; then the SMN1/SMN2 locus
# (chr5 slice) gets its copy number from the EUR model.
set -e
: "${PARASCOPY_DATA_VERSION:?versions.env sets it}"
Z=https://zenodo.org/records/15019940/files
for f in "GRCh38_v${PARASCOPY_DATA_VERSION}.tar.gz a95bf674f43317d3a4c1b8ddbb140945" \
         "models_GRCh38_1KGP_v${PARASCOPY_DATA_VERSION}.tar.gz 244110d8fa883cf11334527ec2383498"; do
  set -- $f
  python3 -c 'import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], sys.argv[2])' "${Z}/$1" "$1"
  echo "$2  $1" | md5sum -c -
  tar -xzf "$1"
done
awk 'BEGIN { for (s = 10050000; s < 10450000; s += 100) printf "chr20\t%d\t%d\n", s, s + 100 }' > windows.bed
parascopy depth -i /in/HG002_slice.bam::HG002 -f /in/ref.fa -b windows.bed -o depth -@ 4
parascopy cn-using models_GRCh38_1KGP/EUR/SMN1.gz -i /in/HG002_slice.bam::HG002 -f /in/ref.fa \
  -t homology_table/GRCh38.bed.gz -d depth -o cn -@ 4
gzip -dc cn/res.samples.bed.gz | grep -v '^##' | cut -f1-11 | tee smn.tsv
