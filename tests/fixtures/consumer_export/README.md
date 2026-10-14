# Site table slice for CI

`sites.tsv.gz` is the site table of step 40 (`bin/export_consumer.py
build-sites`, the format `setup.sh --consumer-sites` writes) cut to the e2e
fixture's regions. The e2e case `consumer-export-1` and the all-tools stub
profile (`conf/test_all.config`) read it.

How it was made:

1. dbSNP b157's GRCh38 VCF (`GCF_000001405.40.gz`, the build
   `CONSUMER_DBSNP_BUILD` names) was read over HTTPS with `tabix` for each
   region of the fixture's `regions.bed` (release `fixture-v6`), keeping the
   records whose rsID is in the 23andMe v5 site list of
   `CONSUMER_SITES_COMMIT`: 7,667 rsIDs.
2. The list's rows of those rsIDs and those dbSNP records went through
   `export_consumer.py build-sites`: 7,623 sites kept; 38 dropped as not an
   SNV in dbSNP, and 6 pseudoautosomal chrX sites whose chrX record lies
   outside the fixture's chrX slices.

Check of the GRCh37 side: the kept sites were looked up in dbSNP b157's
GRCh37 VCF (`GCF_000001405.25.gz`) at their GRCh37 position. Off chrM,
7,337 of 7,338 sites are there with the same reference base (the other,
on chr5, has no GRCh37 record at that position). The chrM sites carry
dbSNP's rCRS position (see docs/40-export-consumer.md).

The slice has rsIDs, positions and reference bases only: no genotype of
anyone.
