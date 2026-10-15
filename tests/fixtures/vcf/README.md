# Synthetic VCFs

Invented records for `tests/test_vcf_precheck_samples.sh` and the
`validate-multisample-vcf` fake-docker case. No record comes from a real
genome: positions, alleles and the sample names `SAMPLE_A` and `SAMPLE_B`
are made up (except the REF bases of `grch37_no_contigs.vcf`, below, which
are the reference genome's).

- `one_sample.vcf`: one sample column; every intake check accepts it.
- `two_samples.vcf`: the same kind of records with two sample columns, the
  shape of a joint-called file; `VCF_PRECHECK` and `validate-setup.sh`
  refuse it.
- `grch37_no_contigs.vcf`, for the e2e case
  `tests/e2e/vcf-precheck-ref-spotcheck-1.sh`: 60 SNVs on chr1 and chr20
  between 10.0 and 10.1 Mb, one every 3,301 bases, with no `##contig` lines.
  The REF bases are GRCh37's (UCSC hg19) at those positions; the ALT alleles
  and the genotypes are made up. At the same positions GRCh38 has another
  base in 45 of the 60, so `VCF_PRECHECK` stops the file by its REF bases.

The tests compress them into a temporary directory; nothing here is read
by a real run.
