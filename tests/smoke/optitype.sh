# shellcheck shell=sh
# OPTITYPE_IMAGE row: type the read pairs of the fixture's HLA slice
# (chr6:29.9-33.1 Mb) with OptiType, as step 08 does with HLA_QC=true, with
# the solver the step names. The reads are unpacked first: step 08 gives
# OptiType plain FASTQ.
set -e
optitype check-deps
gzip -dc /in/hla_R1.fq.gz > hla_R1.fq
gzip -dc /in/hla_R2.fq.gz > hla_R2.fq
optitype run -i hla_R1.fq -i hla_R2.fq --dna --solver glpk --threads 4 -o ot -p HG002
cp "$(find ot -name '*_result.tsv' | head -n 1)" HG002_result.tsv
cat HG002_result.tsv
rm -f hla_R1.fq hla_R2.fq
