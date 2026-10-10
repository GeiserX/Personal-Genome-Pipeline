# Synthetic VCFs

Invented records for `tests/test_vcf_precheck_samples.sh` and the
`validate-multisample-vcf` fake-docker case. No record comes from a real
genome: positions, alleles and the sample names `SAMPLE_A` and `SAMPLE_B`
are made up.

- `one_sample.vcf`: one sample column; every intake check accepts it.
- `two_samples.vcf`: the same kind of records with two sample columns, the
  shape of a joint-called file; `VCF_PRECHECK` and `validate-setup.sh`
  refuse it.

The tests compress them into a temporary directory; nothing here is read
by a real run.
