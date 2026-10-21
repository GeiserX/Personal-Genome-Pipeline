# Made-up reads for the compound-het phase (step 31)

`make_fixture.py OUT_DIR [--swap]` writes the fixture: read pairs of 100 bp
from 500 bp fragments on a made-up 80 kb contig (`chrP`, seeded random
sequence), a VEP VCF of six heterozygous variants for step 31, and the three
pairs as slivar names them. The reads and variants are made up; they describe
no person. The same arguments give the same files, so nothing binary is kept
in git: the tests build the BAM from `reads.sam` (pysam in
`tests/test_comphet_phase.py`, samtools in
`tests/e2e/comphet-phasing-1-reads.sh`).

| Gene | Variants | Copies | Expected phase |
|---|---|---|---|
| GENEC | SNV at 1000, 2-base deletion at 1003 | both on copy 1 | cis, from 22 single reads |
| GENET | SNVs at 5000 and 5400 | one on each copy | trans, from 20 read pairs (no single read covers both) |
| GENEF | SNVs at 20000 and 70000 | one on each copy | unknown: 50 kb apart, no read pair covers both |

`--swap` moves GENEC's deletion to copy 2: the same positions then read
trans. That is the control that the cis verdict comes from the alleles in the
reads, not from the distance.
