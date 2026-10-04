# shellcheck shell=sh
# TELOMEREHUNTER_IMAGE row. The fixture BAM has no telomeric read, so its
# tel_content is 0 whatever the tool does. This plants 1,100 unmapped
# telomeric reads of 150 bp (600 TTAGGG, 400 CCCTAA, 100 TTAGGG with one
# TCAGGG in the middle) with the image's own samtools, then runs TelomereHunter
# as step 10 does (--plotNone). PLANTED=0 plants none, the control that shows
# the tel_content check fails without them.
set -e
n=${PLANTED:-1}
awk -v n="$n" 'BEGIN {
  q = ""; for (i = 0; i < 150; i++) q = q "I"
  if (!n) exit
  for (i = 1; i <= 1100; i++) {
    if (i <= 600) { s = ""; for (k = 0; k < 25; k++) s = s "TTAGGG" }
    else if (i <= 1000) { s = ""; for (k = 0; k < 25; k++) s = s "CCCTAA" }
    else { s = ""; for (k = 0; k < 12; k++) s = s "TTAGGG"; s = s "TCAGGG"; for (k = 0; k < 12; k++) s = s "TTAGGG" }
    printf "telo%d\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tRG:Z:HG002\n", i, s, q
  }
}' > planted.sam
samtools view -h /in/HG002_slice.bam | cat - planted.sam | samtools view -b -o telo.bam -
samtools index telo.bam
echo "telo.bam: $(samtools view -c telo.bam) reads, $(samtools view -c -f 4 telo.bam) unmapped"
telomerehunter -ibt telo.bam -o th -p HG002 --plotNone
cat th/HG002/HG002_summary.tsv
