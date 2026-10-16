# Synthetic VCFs

Invented records for `tests/test_vcf_precheck_samples.sh`, the
`validate-multisample-vcf` fake-docker case and the combined-allele tests.
No record comes from a real genome: positions, alleles and the sample names
`SAMPLE_A` and `SAMPLE_B` are made up. The exceptions are the REF bases of
`grch37_no_contigs.vcf` and the positions and REF bases of the two
combined-allele files, all described below, which are the reference genome's.

- `one_sample.vcf`: one sample column; every intake check accepts it.
- `two_samples.vcf`: the same kind of records with two sample columns, the
  shape of a joint-called file; `VCF_PRECHECK` and `validate-setup.sh`
  refuse it.
- `combined_allele_sample.vcf` and `combined_allele_clinvar.vcf`: a delins
  written as two records 3 bp apart (a 2 bp deletion and a 3 bp insertion)
  plus a single SNV, and a ClinVar-style file with a 2-star record for the
  combined allele, 1-star records for its two parts and a record for the
  SNV. Positions are on GRCh38 chr22 and the REF bases are the reference's
  there, so step 06 runs on them with the real fixture reference; the genes
  (`SYNTH1`, `SYNTH2`) and ClinVar IDs are made up.
  `tests/test_clinvar_join_adjacent.py` and the e2e case
  `clinvar-combined-allele-1` read them.
- `grch37_no_contigs.vcf`, for the e2e case
  `tests/e2e/vcf-precheck-ref-spotcheck-1.sh`: 60 SNVs on chr1 and chr20
  between 10.0 and 10.1 Mb, one every 3,301 bases, with no `##contig` lines.
  The REF bases are GRCh37's (UCSC hg19) at those positions; the ALT alleles
  and the genotypes are made up. At the same positions GRCh38 has another
  base in 45 of the 60, so `VCF_PRECHECK` stops the file by its REF bases.

The tests compress them into a temporary directory; nothing here is read
by a real run.
