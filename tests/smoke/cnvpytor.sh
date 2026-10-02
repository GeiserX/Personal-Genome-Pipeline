# shellcheck shell=sh
# CNVPYTOR_IMAGE row (runs as root). The image ships without the GC and mask
# resource files, and cnvpytor stops every command (exit 0, no calls) unless
# all seven exist. Step 18 mounts the pinned files over the package's data
# directory; the mini reference matches no genome they describe, so empty
# files stand in for them here.
set -e
D=$(python -c 'import cnvpytor, os; print(os.path.dirname(cnvpytor.__file__) + "/data")' | tail -n 1)
echo "resource directory: ${D}"
for f in gc_hg19 mask_hg19 gc_hg38 mask_hg38 gc_chm13v2.0 gc_chm13v1.1 gc_kn99; do
  [ -e "${D}/${f}.pytor" ] || touch "${D}/${f}.pytor"
done
cnvpytor -root mini.pytor -rd /in/mini.bam -chrom chr20
cnvpytor -root mini.pytor -his 1000
cnvpytor -root mini.pytor -partition 1000
cnvpytor -root mini.pytor -call 1000 > calls.tsv
cat calls.tsv
