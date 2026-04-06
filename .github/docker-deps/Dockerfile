# This file exists solely for Dependabot to track Docker image versions.
# It is NOT used to build any container.
#
# When Dependabot opens a PR updating an image tag below, also update
# versions.env (the canonical source) and re-run ./scripts/setup.sh.
#
# Images pinned to :latest have no versioned tags from their publisher.
# Dependabot cannot suggest updates for these.

# Core tools
FROM quay.io/biocontainers/minimap2:2.28--he4a0461_0 AS minimap2
FROM staphb/samtools:1.20 AS samtools
FROM staphb/bcftools:1.21 AS bcftools
FROM google/deepvariant:1.6.0 AS deepvariant
FROM quay.io/biocontainers/fastp:1.3.1--h43da1c4_0 AS fastp

# Structural variant callers
FROM quay.io/biocontainers/manta:1.6.0--h9ee0642_2 AS manta
FROM quay.io/biocontainers/delly:1.7.3--hd6466ae_0 AS delly
FROM quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11 AS cnvnator
FROM quay.io/biocontainers/gridss:2.13.2--h96c455f_6 AS gridss
FROM brentp/duphold:v0.2.3 AS duphold
FROM getwilds/annotsv:3.4.4 AS annotsv

# Variant annotation
FROM ensemblorg/ensembl-vep:release_112.0 AS vep
FROM sigven/pcgr:2.2.5 AS pcgr

# Pharmacogenomics
FROM pgkb/pharmcat:3.2.0 AS pharmcat

# Specialized analysis
# telomerehunter: only :latest available from publisher
FROM lgalarno/telomerehunter:latest AS telomerehunter
# haplogrep3: only :latest available from publisher
FROM genepi/haplogrep3:latest AS haplogrep3
FROM quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 AS t1k
FROM quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5 AS expansionhunter
FROM quay.io/biocontainers/goleft:0.2.4--h9ee0642_1 AS goleft
FROM quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0 AS mosdepth

# General purpose
FROM broadinstitute/gatk:4.6.1.0 AS gatk
FROM broadinstitute/picard:3.4.0 AS picard
FROM python:3.11-slim AS python
FROM pgscatalog/plink2:2.00a5.10 AS plink2
FROM quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0 AS multiqc

# Alternative variant callers
FROM dancooke/octopus:0.7.4 AS octopus
FROM hkubal/clair3:v2.0.0 AS clair3
FROM quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0 AS sniffles
