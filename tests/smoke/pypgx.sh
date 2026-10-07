# shellcheck shell=sh
# PYPGX_IMAGE row: the CYP2D6 path of step 32 (depth of coverage, VDR control
# statistics, run-ngs-pipeline) on the GIAB slice, with the pypgx-bundle tag
# PYPGX_BUNDLE_VERSION names. pypgx makes the input VCF itself.
set -e
export PYPGX_BUNDLE=/in/pypgx-bundle
ln -sfn /in/pypgx-bundle "${HOME}/pypgx-bundle"
BAM=/in/HG002_slice.bam
pypgx create-input-vcf input.vcf.gz /in/ref.fa "$BAM" --assembly GRCh38 --genes CYP2D6
pypgx prepare-depth-of-coverage depth.zip "$BAM" --assembly GRCh38 --genes CYP2D6
pypgx compute-control-statistics VDR control.zip "$BAM" --assembly GRCh38
pypgx run-ngs-pipeline CYP2D6 CYP2D6 --variants input.vcf.gz --depth-of-coverage depth.zip \
  --control-statistics control.zip --assembly GRCh38
python3 - <<'PY' > CYP2D6.genotype.txt
import csv, io, zipfile
z = zipfile.ZipFile("CYP2D6/results.zip")
name = next(n for n in z.namelist() if n.endswith("data.tsv"))
rows = list(csv.DictReader(io.TextIOWrapper(z.open(name)), delimiter="\t"))
print(rows[0]["Genotype"])
PY
cat CYP2D6.genotype.txt
