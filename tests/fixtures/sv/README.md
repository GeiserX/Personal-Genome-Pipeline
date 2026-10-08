# Synthetic SV calls for step 22

Three caller VCFs on chr20 with a known answer, read by
`tests/e2e/sv-mito-telomere-steps-1-sv-merge.sh` and by the `SURVIVOR_IMAGE`
row of `tests/smoke/commands.tsv`. The records are made up; they describe no
person.

| Event | manta.vcf.in | delly.vcf.in | cnvpytor.vcf.in | Consensus |
|---|---|---|---|---|
| one 2 kb deletion, called 2 bp apart across a 1 kb boundary | 10,100,999-10,102,999 | 10,101,001-10,103,001 | | one record, SUPP=2 |
| two deletions starting in the same 1 kb window, 600 bp and 5 kb long | 10,200,100-10,200,700 | 10,200,300-10,205,300 | | none: the ends are 4.6 kb apart |
| a deletion one caller saw | | | 10,300,000-10,310,000 | none |

The position binning step 22 used before (chromosome, `int(POS/1000)` and
SVTYPE) gets the first two wrong: it splits the first event over two windows
and counts the second pair as agreement.
