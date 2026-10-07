# shellcheck shell=sh
# T1K_IMAGE row: build the HLA index from IPD-IMGT/HLA as step 08 does, then
# genotype the reads of the fixture's HLA slice (chr6:29.9-33.1 Mb).
set -e
t1k-build.pl -o hlaidx --download IPD-IMGT/HLA
SEQ=$(ls hlaidx/*_dna_seq.fa 2>/dev/null | head -n 1)
if [ -z "$SEQ" ]; then
  t1k-build.pl -d hlaidx/hla.dat -o hlaidx_seq
  SEQ=$(ls hlaidx_seq/*_dna_seq.fa | head -n 1)
fi
echo "allele sequences: ${SEQ}"
run-t1k -1 /in/hla_R1.fq.gz -2 /in/hla_R2.fq.gz -f "$SEQ" --preset hla-wgs -t 4 --od t1k -o HG002_hla
cat t1k/HG002_hla_genotype.tsv
