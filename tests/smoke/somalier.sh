# shellcheck shell=sh
# SOMALIER_IMAGE row: step 33's two somalier calls on the fixture slice, with
# the sites of somalier's GRCh38 file that fall inside the fixture's regions
# (somalier_sites.vcf.in, 92 of its 17,766 sites). extract reads the sample
# name from the @RG SM tag; relate writes the per-sample table step 33 reads
# (with --sites, which its hom-ref and hom-alt counts need).
set -e
somalier extract -d ex --sites /smoke/somalier_sites.vcf.in -f /in/ref.fa /in/HG002_slice.bam
ls ex
somalier relate --sites /smoke/somalier_sites.vcf.in -o rel ex/*.somalier
cat rel.samples.tsv
