# Step 0: Reference Data Setup

One-time downloads required before running the pipeline.

## GRCh38 Reference Genome

```bash
GENOMA_DIR=/path/to/genome/data
mkdir -p ${GENOMA_DIR}/reference
cd ${GENOMA_DIR}/reference

# Download GRCh38 reference (3.1GB)
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai

# Create minimap2 index (7GB, ~30 min)
docker run --rm -v ${GENOMA_DIR}/reference:/ref staphb/samtools:1.20 \
  minimap2 -d /ref/GRCh38.mmi /ref/Homo_sapiens_assembly38.fasta
```

## ClinVar Database

```bash
cd ${GENOMA_DIR}/reference

# Download latest ClinVar (updated monthly)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# ClinVar uses "1,2,3..." chromosomes. Our BAMs use "chr1,chr2,chr3..."
# Create chr-prefixed version:
docker run --rm -v ${GENOMA_DIR}/reference:/ref staphb/bcftools:1.21 bash -c '
  echo -e "1 chr1\n2 chr2\n3 chr3\n4 chr4\n5 chr5\n6 chr6\n7 chr7\n8 chr8\n9 chr9\n10 chr10\n11 chr11\n12 chr12\n13 chr13\n14 chr14\n15 chr15\n16 chr16\n17 chr17\n18 chr18\n19 chr19\n20 chr20\n21 chr21\n22 chr22\nX chrX\nY chrY\nMT chrM" > /ref/chr_rename.txt
  bcftools annotate --rename-chrs /ref/chr_rename.txt /ref/clinvar.vcf.gz -Oz -o /ref/clinvar_chr.vcf.gz
  bcftools index -t /ref/clinvar_chr.vcf.gz
'

# Extract pathogenic/likely pathogenic only (faster for screening):
docker run --rm -v ${GENOMA_DIR}/reference:/ref staphb/bcftools:1.21 bash -c '
  bcftools view -i "CLNSIG~\"Pathogenic\" || CLNSIG~\"Likely_pathogenic\"" /ref/clinvar_chr.vcf.gz -Oz -o /ref/clinvar_pathogenic_chr.vcf.gz
  bcftools index -t /ref/clinvar_pathogenic_chr.vcf.gz
'
```

## VEP Cache (~26GB)

```bash
mkdir -p ${GENOMA_DIR}/vep_cache/tmp
cd ${GENOMA_DIR}/vep_cache/tmp

# Download (manual wget is more reliable than VEP INSTALL.pl)
wget -c https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz

# Extract to parent directory
cd ${GENOMA_DIR}/vep_cache
tar xzf tmp/homo_sapiens_vep_112_GRCh38.tar.gz
# Creates: ${GENOMA_DIR}/vep_cache/homo_sapiens/112_GRCh38/
```

## T1K HLA Reference (~30 min)

```bash
mkdir -p ${GENOMA_DIR}/t1k_idx

# Step 1: Download IPD-IMGT/HLA database (~2 min)
docker run --rm -v ${GENOMA_DIR}/t1k_idx:/idx \
  quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
  t1k-build.pl -o /idx/hlaidx --download IPD-IMGT/HLA

# Step 2: Build coordinate file from genome (~30 min, reads entire 3.1GB FASTA)
# CRITICAL: Use the actual FASTA, NOT the .fai index!
docker run --rm --cpus 4 --memory 8g \
  -v ${GENOMA_DIR}:/genoma \
  quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
  t1k-build.pl \
    -d /genoma/t1k_idx/hlaidx/hla.dat \
    -g /genoma/reference/Homo_sapiens_assembly38.fasta \
    -o /genoma/t1k_idx/hlaidx_grch38
```

## Docker Images to Pre-Pull

```bash
# Core tools
docker pull staphb/bcftools:1.21
docker pull staphb/samtools:1.20
docker pull google/deepvariant:1.6.0

# Analysis tools
docker pull getwilds/annotsv:latest
docker pull weisburd/expansionhunter:latest
docker pull lgalarno/telomerehunter:latest
docker pull quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0
docker pull ensemblorg/ensembl-vep:release_112.0
docker pull pgkb/pharmcat:2.15.5

# Optional (HLA-LA with pre-built graph, 4.5GB)
docker pull jiachenzdocker/hla-la:latest
```

## Disk Space Summary

| Resource | Size | Notes |
|---|---|---|
| GRCh38 FASTA + FAI | 3.2 GB | One-time |
| minimap2 index | 7 GB | One-time |
| ClinVar DB | 200 MB | Update monthly |
| VEP cache | 26 GB | One-time per release |
| T1K HLA index | 450 MB | One-time |
| Per sample (BAM) | 30-40 GB | Keep |
| Per sample (VCF) | 100 MB | Keep |
| Per sample (intermediates) | ~5 GB | Can clean up |
| **Total per sample** | **~35-45 GB** | |
