# shellcheck shell=sh
# HAPLOGREP3_IMAGE row. haplogrep3 reads haplogrep3.yaml and its trees from
# the working directory; step 12 runs it in the image's own working
# directory, so the row moves to the directory of the binary first.
set -e
BIN=$(readlink -f "$(command -v haplogrep3)")
cd "$(dirname "$BIN")"
echo "running in $(pwd)"
haplogrep3 classify --tree phylotree-fu-rcrs@1.2 --input /in/mito.vcf.gz \
  --output /out/haplogroup.txt --extend-report
cat /out/haplogroup.txt
