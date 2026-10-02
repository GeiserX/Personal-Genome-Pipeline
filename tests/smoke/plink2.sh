# shellcheck shell=sh
# PLINK2_IMAGE row: the two plink2 calls of step 25 (make-pgen, then score) with
# a score file of five SNVs of the sample, effect allele ALT, weight 1.
set -e
plink2 --vcf /in/sample.vcf.gz --make-pgen --out p --threads 4 --memory 2000 \
  --set-all-var-ids '@:#' --new-id-max-allele-len 100 --chr 1-22 --allow-extra-chr --output-chr chrM
awk '!/^#/ && length($4) == 1 && length($5) == 1 && n < 5 {print $3 "\t" $5 "\t1"; n++}' p.pvar > score.tsv
cat score.tsv
plink2 --pfile p --score score.tsv 1 2 3 ignore-dup-ids no-mean-imputation cols=+scoresums \
  --out score --threads 4 --memory 2000 --allow-extra-chr
cat score.sscore
