# shellcheck shell=sh
# SLIVAR_IMAGE row: compound-hets exactly as step 31 calls it, then once more
# on a VCF whose heterozygous calls slivar expr tagged, so pairs must come out
# (the fixture's HLA genes hold many heterozygous variants).
set -e
printf 'HG002\tHG002\t0\t0\t0\t-9\n' > sample.ped
slivar compound-hets --allow-non-trios --vcf /in/sample.vcf.gz --ped sample.ped > step31.vcf
slivar expr --vcf /in/sample.vcf.gz --ped sample.ped --pass-only \
  --sample-expr 'comphet_side:sample.het' -o tagged.vcf
slivar compound-hets --allow-non-trios --vcf tagged.vcf --ped sample.ped > comphet.vcf
