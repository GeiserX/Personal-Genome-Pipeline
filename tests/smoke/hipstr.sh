# shellcheck shell=sh
# HIPSTR_IMAGE row: what step 45 runs, on the fixture's GIAB slice BAM. Eleven
# Y-STR regions of HipSTR's annotated marker set (the source of
# assets/ystr_panel.tsv) lie in the fixture's chrY slices; three of them
# (DYS570, DYS388, DYS442) are FamilyTreeDNA markers. HipSTR runs with the
# step's flags, and the image's Python, which bin/ystr_report.py runs on,
# reads the bgzipped VCF.
set -e
printf 'chrY\t%s\t%s\t%s\t%s\t%s\n' \
  2954700 2954743 4 11 DYS453   2974166 2974201 4 9 DYS579    6993190 6993257 4 17 DYS570 \
  12635604 12635639 3 12 DYS388 12649170 12649229 4 15 DYS442 14611544 14611591 4 12 DYS558 \
  14635532 14635551 5 4 DYS571.1 14635563 14635602 4 10 DYS571.2 19647924 19647959 3 12 DYS489 \
  20607709 20607741 3 11 DYS470.1 20607773 20607784 3 4 DYS470.2 > regions.bed
HipSTR --bams /in/HG002_slice.bam --fasta /in/ref.fa --regions regions.bed --str-vcf calls.vcf.gz \
  --bam-samps HG002 --bam-libs HG002 --haploid-chrs chrY --def-stutter-model --min-reads 5 \
  --max-str-len 150 --output-filters --log hipstr.log
python3 --version > python_version.txt 2>&1
python3 - <<'PY'
import gzip
called = []
for line in gzip.open("calls.vcf.gz", "rt"):
    if line.startswith("#"):
        continue
    c = line.rstrip("\n").split("\t")
    f = dict(zip(c[8].split(":"), c[9].split(":")))
    if f.get("GT", ".") != "." and "GB" in f and "DP" in f and "Q" in f and "DSTUTTER" in f:
        called.append(f"{c[2]}\t{f['GT']}\t{f['GB']}\t{f['DP']}\t{f['Q']}")
open("called.tsv", "w").write("".join(x + "\n" for x in called))
print("\n".join(called))
PY
