# Step 0: Reference Data Setup

One-time downloads required before running the pipeline. Core downloads: ~70-75 GB. Including all optional resources (somatic, T1K, GIAB): ~78 GB download, ~96 GB after extraction.

> **Estimated time:** 1-3 hours depending on internet speed. VEP caches (26 GB each for steps 13 and 17) are the largest downloads.

## GRCh38 Reference Genome

The foundation for everything. All tools need this.

```bash
export GENOME_DIR=/path/to/your/data
mkdir -p ${GENOME_DIR}/reference
cd ${GENOME_DIR}/reference

# Download GRCh38 reference (3.1 GB)
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai

# Verify the download
md5sum Homo_sapiens_assembly38.fasta
# Expected: 64b32de2fc934679c16e83a2bc072064
```

### Why GRCh38?

This pipeline uses **GRCh38** (also called hg38) exclusively. It's the current standard genome build with:
- Corrected mitochondrial sequence (rCRS)
- Better representation of centromeres and telomeres
- ALT contigs for highly polymorphic regions (MHC, KIR)
- `chr` prefix naming (chr1, chr2, ..., chrX, chrY, chrM)

If your data is on **GRCh37/hg19**, extract FASTQ from BAM and re-align. See [vendor-guide.md](vendor-guide.md#genome-build-grch37-hg19-vs-grch38-hg38).

## ClinVar Database

Updated monthly by NCBI. Contains known pathogenic/benign variant classifications.

```bash
mkdir -p ${GENOME_DIR}/clinvar
cd ${GENOME_DIR}/clinvar

# Download latest ClinVar (~200 MB)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# ClinVar uses "1,2,3..." chromosomes. Our BAMs use "chr1,chr2,chr3..."
# Create chr-prefixed version:
docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 bash -c '
  echo -e "1 chr1\n2 chr2\n3 chr3\n4 chr4\n5 chr5\n6 chr6\n7 chr7\n8 chr8\n9 chr9\n10 chr10\n11 chr11\n12 chr12\n13 chr13\n14 chr14\n15 chr15\n16 chr16\n17 chr17\n18 chr18\n19 chr19\n20 chr20\n21 chr21\n22 chr22\nX chrX\nY chrY\nMT chrM" > /genome/clinvar/chr_rename.txt
  bcftools annotate --rename-chrs /genome/clinvar/chr_rename.txt /genome/clinvar/clinvar.vcf.gz -Oz -o /genome/clinvar/clinvar_chr.vcf.gz
  bcftools index -t /genome/clinvar/clinvar_chr.vcf.gz
'

# Extract pathogenic/likely pathogenic only (faster for screening):
docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 bash -c '
  bcftools view -i "CLNSIG~\"Pathogenic\" || CLNSIG~\"Likely_pathogenic\"" /genome/clinvar/clinvar_chr.vcf.gz -Oz -o /genome/clinvar/clinvar_pathogenic_chr.vcf.gz
  bcftools index -t /genome/clinvar/clinvar_pathogenic_chr.vcf.gz
'
```

> **Tip:** Re-download ClinVar monthly for the latest classifications. ClinVar adds ~1000 new pathogenic variants per month.

## VEP Cache (~26 GB)

Ensembl Variant Effect Predictor annotation database. Required for step 13.

```bash
mkdir -p ${GENOME_DIR}/vep_cache/tmp
cd ${GENOME_DIR}/vep_cache/tmp

# Download (manual wget is more reliable than VEP INSTALL.pl)
# The -c flag enables resume if the download is interrupted
wget -c https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz

# Extract to parent directory (~30 GB extracted)
cd ${GENOME_DIR}/vep_cache
tar xzf tmp/homo_sapiens_vep_112_GRCh38.tar.gz
# Creates: ${GENOME_DIR}/vep_cache/homo_sapiens/112_GRCh38/

# Optional: delete the tarball to save 26 GB
# rm tmp/homo_sapiens_vep_112_GRCh38.tar.gz
```

> **Warning:** The VEP `INSTALL.pl` script downloads to a temporary directory that may lack write permissions inside Docker. Always download manually with `wget -c`. See [lessons-learned.md](lessons-learned.md) for details.

## PCGR/CPSR Ref Data Bundle (~5 GB)

Required for step 17 (CPSR cancer predisposition screening). Includes ClinVar, gnomAD, CancerMine, and other databases. PCGR 2.x uses a separate, smaller ref data bundle — VEP cache is mounted independently.

> **Important:** PCGR 2.2.5 bundles VEP 113, which requires a **release-113** cache — different from the release-112 cache used by step 13 above. See the next section for the VEP 113 download.

```bash
mkdir -p ${GENOME_DIR}/pcgr_data
cd ${GENOME_DIR}/pcgr_data

# Download (~5 GB)
wget -c https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz

# Extract and organize into version-stamped directory
tar xzf pcgr_ref_data.20250314.grch38.tgz
mkdir -p 20250314 && mv data/ 20250314/
# Creates: ${GENOME_DIR}/pcgr_data/20250314/data/

# Optional: delete the tarball to save 5 GB
# rm pcgr_ref_data.20250314.grch38.tgz
```

## VEP 113 Cache for CPSR (~26 GB)

PCGR 2.2.5 (step 17) bundles VEP 113 internally, which needs the **release-113** cache. This is separate from the release-112 cache used by step 13. Both coexist in the same `vep_cache/` directory under different subdirectories (`112_GRCh38/` and `113_GRCh38/`).

```bash
mkdir -p ${GENOME_DIR}/vep_cache/tmp
cd ${GENOME_DIR}/vep_cache/tmp

# Download VEP 113 cache (~26 GB)
wget -c https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz

# Extract alongside the existing release-112 cache
cd ${GENOME_DIR}/vep_cache
tar xzf tmp/homo_sapiens_vep_113_GRCh38.tar.gz
# Creates: ${GENOME_DIR}/vep_cache/homo_sapiens/113_GRCh38/
```

> If you only run step 13 (VEP annotation) and skip step 17 (CPSR), you only need the release-112 cache. If you only run step 17, you only need release-113.

## T1K HLA Reference (Optional)

Only needed for step 8 (HLA typing). ~30 minutes to build.

```bash
mkdir -p ${GENOME_DIR}/t1k_idx

# Step 1: Download IPD-IMGT/HLA database (~2 min)
docker run --rm -v ${GENOME_DIR}/t1k_idx:/idx \
  quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
  t1k-build.pl -o /idx/hlaidx --download IPD-IMGT/HLA

# Step 2: Build coordinate file from genome (~30 min, reads entire 3.1GB FASTA)
# CRITICAL: Use the actual FASTA, NOT the .fai index!
docker run --rm --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
  t1k-build.pl \
    -d /genome/t1k_idx/hlaidx/hla.dat \
    -g /genome/reference/Homo_sapiens_assembly38.fasta \
    -o /genome/t1k_idx/hlaidx_grch38
```

> **Note:** HLA typing from WGS is challenging. T1K coordinates may have ~50% unmapped alleles. For clinical HLA typing, dedicated lab assays are more reliable. See [lessons-learned.md](lessons-learned.md#t1k-coordinate-file-with-wrong-values).

## Somatic Calling Resources (Optional, for Step 29)

Only needed if you plan to run step 29 (somatic variant calling with Mutect2 tumor-only mode). These resources significantly reduce false positives. See [29-mutect2-somatic.md](29-mutect2-somatic.md) for details.

```bash
mkdir -p ${GENOME_DIR}/somatic

# gnomAD AF-only VCF (~6.5 GB) — germline allele frequencies for filtering
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
  -O ${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi \
  -O ${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz.tbi

# Panel of Normals (~1 GB) — recurrent technical artifacts from 1000 Genomes
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz \
  -O ${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi \
  -O ${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz.tbi
```

## Annotation Databases (Optional, for Steps 30-31)

Deep pathogenicity scoring and variant prioritization. These databases power vcfanno (step 30) and slivar (step 31). All are optional — step 30 gracefully skips any missing track.

> **Total download:** ~104 GB (CADD is 83 GB alone). If disk space is tight, start with REVEL + AlphaMissense (~1.2 GB combined) — they provide the highest value per byte for missense variant interpretation.

### CADD v1.7 Pre-scored (SNVs + Indels) — ~83 GB

The most widely used deleteriousness metric. Scores all possible SNVs genome-wide plus gnomAD indels.

> **License:** Free for non-commercial/academic use. Commercial use requires a license from the University of Washington. The pipeline does not redistribute CADD data — users download directly from the source.

```bash
mkdir -p ${GENOME_DIR}/annotations
cd ${GENOME_DIR}/annotations

# Whole-genome SNV scores (~81.5 GB, pre-indexed)
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi

# gnomAD v4.0 indel scores (~1.2 GB, pre-indexed)
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi
```

> **Note:** CADD TSVs use chromosome names without `chr` prefix (1, 2, 3...). vcfanno handles this via column mapping in the TOML config — no manual renaming needed.

### SpliceAI Pre-scored — ~20 GB

Deep learning splice-site variant prediction. Catches pathogenic intronic variants that VEP's rule-based splice prediction misses.

```bash
cd ${GENOME_DIR}/annotations

# SNV splice scores (~15.5 GB)
wget -c https://download.molgeniscloud.org/downloads/vip/resources/GRCh38/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz
wget -c https://download.molgeniscloud.org/downloads/vip/resources/GRCh38/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz.tbi

# Indel splice scores (~4.5 GB)
wget -c https://download.molgeniscloud.org/downloads/vip/resources/GRCh38/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz
wget -c https://download.molgeniscloud.org/downloads/vip/resources/GRCh38/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz.tbi
```

> **Alternative source:** Illumina BaseSpace at `https://basespace.illumina.com/s/otSPW8hnhaZR` (requires free account).

### REVEL v1.3 — ~526 MB

Ensemble pathogenicity scoring for missense variants, combining 13 individual tools. Recommended by ClinGen for missense variant classification.

```bash
cd ${GENOME_DIR}/annotations

# Download and prepare for vcfanno
wget -c https://zenodo.org/record/7072866/files/revel-v1.3_all_chromosomes.zip
unzip revel-v1.3_all_chromosomes.zip

# Convert to tabix-indexed TSV for GRCh38
# Extract GRCh38 columns, add chr prefix, sort, bgzip, index
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c '
    cd /genome/annotations
    cat revel_with_transcript_ids | tr "," "\t" > revel_tabbed.tsv
    # Header line
    head -1 revel_tabbed.tsv | sed "s/^/#/" > revel_grch38.tsv
    # Extract GRCh38 rows (col3=grch38_pos), skip header, add chr prefix, sort
    tail -n+2 revel_tabbed.tsv | \
      awk -F"\t" "\$3 != \".\" {print \"chr\"\$1\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$8}" | \
      sort -k1,1V -k2,2n >> revel_grch38.tsv
    bgzip revel_grch38.tsv
    tabix -s 1 -b 2 -e 2 revel_grch38.tsv.gz
    # Clean up intermediates
    rm -f revel_tabbed.tsv revel_with_transcript_ids revel-v1.3_all_chromosomes.zip
  '
```

> **Thresholds:** ClinGen recommends REVEL >= 0.644 as supporting evidence of pathogenicity (PP3_Moderate) and >= 0.932 for strong evidence (PP3_Strong) for missense variants.

### AlphaMissense — ~613 MB

DeepMind's protein-structure-informed missense classifier. Complements REVEL with structural context. Covers all ~71 million possible human missense variants.

> **License:** CC BY-NC-SA 4.0 (non-commercial, share-alike). The pipeline does not redistribute these scores.

```bash
cd ${GENOME_DIR}/annotations

# Download pre-scored GRCh38 predictions
wget -c https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz

# Index for vcfanno (skip header lines starting with #)
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  tabix -s 1 -b 2 -e 2 -S 1 /genome/annotations/AlphaMissense_hg38.tsv.gz
```

> **Thresholds:** am_pathogenicity < 0.34 = likely benign, > 0.564 = likely pathogenic.

### gnomAD v4.1 Constraint Metrics — ~91 MB

Per-gene loss-of-function intolerance scores. Essential for interpreting novel variants in constrained genes.

```bash
cd ${GENOME_DIR}/annotations

# Gene-level constraint table
wget -c -O gnomad_v4.1_constraint.tsv \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
```

> **Key metrics:** LOEUF (oe_lof_upper) < 0.35 = constrained for LoF; pLI >= 0.9 = intolerant to LoF; missense Z > 3.09 = constrained for missense.

## Docker Images — Pre-Pull All

Pull all images in advance to avoid download delays during analysis:

```bash
# Core pipeline
docker pull quay.io/biocontainers/minimap2:2.28--he4a0461_0
docker pull staphb/samtools:1.20
docker pull staphb/bcftools:1.21
docker pull google/deepvariant:1.6.0

# SV callers
docker pull quay.io/biocontainers/manta:1.6.0--h9ee0642_2
docker pull quay.io/biocontainers/delly:1.7.3--hd6466ae_0
docker pull quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11
docker pull brentp/duphold:v0.2.3

# Annotation
docker pull getwilds/annotsv:3.4.4
docker pull ensemblorg/ensembl-vep:release_112.0
docker pull sigven/pcgr:2.2.5

# Pharmacogenomics
docker pull pgkb/pharmcat:3.2.0

# Specialized
docker pull quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5
docker pull lgalarno/telomerehunter:latest
docker pull genepi/haplogrep3:latest
docker pull quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0
docker pull quay.io/biocontainers/goleft:0.2.4--h9ee0642_1
docker pull broadinstitute/gatk:4.6.1.0

# QC & reporting
docker pull quay.io/biocontainers/fastp:1.3.1--h43da1c4_0
docker pull quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0
docker pull quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0

# Annotation enrichment (v0.4.0)
docker pull quay.io/biocontainers/vcfanno:0.3.7--he881be0_0
docker pull quay.io/biocontainers/slivar:0.3.3--h5f107b1_0
docker pull quay.io/biocontainers/pypgx:0.26.0--pyh7e72e81_0

# Alternative / new callers
docker pull dancooke/octopus:0.7.4
docker pull quay.io/biocontainers/gridss:2.13.2--h96c455f_6
docker pull hkubal/clair3:v2.0.0
docker pull quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0
```

## Alternative Callers (Optional, for Benchmarking)

Only needed if you plan to run alternative variant callers. See [benchmarking.md](benchmarking.md).

```bash
# Alternative aligners
docker pull quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5

# Alternative variant callers (GATK image already pulled above)
docker pull quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2

# Alternative SV caller
docker pull quay.io/biocontainers/tiddit:3.9.5--py312h6e8b409_0

# Alternative small variant caller (SNVs + indels, complements Manta)
docker pull quay.io/biocontainers/strelka:2.9.10--h9ee0642_1

# Benchmarking (truth set evaluation)
docker pull jmcdani20/hap.py:v0.3.12
```

### GATK Sequence Dictionary

GATK HaplotypeCaller requires a `.dict` file alongside the reference FASTA. If you don't have one:

```bash
docker run --rm --user root \
  -v ${GENOME_DIR}:/genome \
  broadinstitute/gatk:4.6.1.0 \
  gatk CreateSequenceDictionary \
    -R /genome/reference/Homo_sapiens_assembly38.fasta
```

### BWA-MEM2 Index

BWA-MEM2 requires its own index files (different from minimap2's `.mmi`). Build once (~1 hour, ~6 GB):

```bash
docker run --rm --user root \
  --cpus 8 --memory 24g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5 \
  bwa-mem2 index /genome/reference/Homo_sapiens_assembly38.fasta
# Creates: .0123, .amb, .ann, .bwt.2bit.64, .pac alongside the FASTA
```

### GIAB Truth Set (for hap.py Benchmarking)

Download a GIAB truth set for benchmarking variant callers. **HG002** (Ashkenazi Jewish male) is preferred because its truth set covers more difficult genomic regions. HG001 (NA12878) is an alternative used in the [quick test](quick-test.md).

**Important:** Truth set benchmarking only works when the query VCF comes from the **same biological sample** as the truth set. You must sequence HG002 (or HG001) DNA, not your own sample. See [benchmarking.md](benchmarking.md) for details.

```bash
mkdir -p ${GENOME_DIR}/giab
cd ${GENOME_DIR}/giab

# HG002 truth set (recommended, GRCh38 v4.2.1)
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# Alternative: HG001/NA12878 (used in quick-test.md)
# wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
# wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
# wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

**Total Docker image size:** ~10-15 GB (compressed, after layer deduplication). Alternative tools add ~3-5 GB.

## Disk Space Summary

| Resource | Download | Extracted | Notes |
|---|---|---|---|
| GRCh38 FASTA + FAI | 3.1 GB | 3.1 GB | Same size (not compressed) |
| ClinVar DB (all versions) | 200 MB | 400 MB | Including chr-prefixed and pathogenic-only |
| VEP cache | 26 GB | 30 GB | Largest single database |
| PCGR/CPSR data bundle | 5 GB | 5 GB | Smaller in PCGR 2.x |
| VEP 113 cache (CPSR) | 26 GB | 30 GB | Separate from step 13's VEP 112 cache |
| T1K HLA index | 50 MB | 450 MB | Optional |
| Somatic resources (gnomAD + PoN) | 7.5 GB | 7.5 GB | Optional (step 29) |
| Annotation databases (CADD+SpliceAI+REVEL+AM+gnomAD) | ~104 GB | ~104 GB | Optional (step 30-31) |
| Docker images | 10-15 GB | 10-15 GB | Cached by Docker engine |
| **Total (core)** | **~78 GB** | **~96 GB** | Without annotation databases |
| **Total (with annotation enrichment)** | **~182 GB** | **~200 GB** | With all v0.4.0 databases |

> **Tip:** If disk space is tight, you can skip the VEP cache (step 13) and PCGR bundle (step 17) initially. The core pipeline (steps 2-3-6-7) only needs the reference FASTA and ClinVar (~3.5 GB total).

## Verifying Your Setup

After all downloads, verify everything is in place:

```bash
echo "Checking reference setup..."
[ -f "${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta" ] && echo "  GRCh38 FASTA: OK" || echo "  GRCh38 FASTA: MISSING"
[ -f "${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta.fai" ] && echo "  FASTA index: OK" || echo "  FASTA index: MISSING"
[ -f "${GENOME_DIR}/clinvar/clinvar_chr.vcf.gz" ] && echo "  ClinVar (chr): OK" || echo "  ClinVar: MISSING"
[ -d "${GENOME_DIR}/vep_cache/homo_sapiens/112_GRCh38" ] && echo "  VEP cache: OK" || echo "  VEP cache: MISSING"
[ -d "${GENOME_DIR}/pcgr_data/20250314/data" ] && echo "  PCGR data: OK" || echo "  PCGR data: MISSING"
[ -d "${GENOME_DIR}/vep_cache/homo_sapiens/113_GRCh38" ] && echo "  VEP 113 cache (CPSR): OK" || echo "  VEP 113 cache (CPSR): MISSING"

echo "Annotation databases (optional, for steps 30-31):"
for DB_PAIR in \
  "whole_genome_SNVs.tsv.gz:CADD SNVs" \
  "gnomad.genomes.r4.0.indel.tsv.gz:CADD indels" \
  "spliceai_scores.raw.snv.hg38.vcf.gz:SpliceAI SNVs" \
  "spliceai_scores.raw.indel.hg38.vcf.gz:SpliceAI indels" \
  "revel_grch38.tsv.gz:REVEL" \
  "AlphaMissense_hg38.tsv.gz:AlphaMissense" \
  "gnomad_v4.1_constraint.tsv:gnomAD constraint"; do
  DB_FILE="${DB_PAIR%%:*}"
  DB_NAME="${DB_PAIR#*:}"
  if [ -f "${GENOME_DIR}/annotations/${DB_FILE}" ]; then
    # Check for .tbi index (not needed for plain TSV files)
    if [[ "$DB_FILE" == *.vcf.gz ]] || [[ "$DB_FILE" == *.tsv.gz && "$DB_FILE" != *constraint* ]]; then
      [ -f "${GENOME_DIR}/annotations/${DB_FILE}.tbi" ] && echo "  ${DB_NAME}: OK (+index)" || echo "  ${DB_NAME}: OK (WARNING: .tbi index missing)"
    else
      echo "  ${DB_NAME}: OK"
    fi
  else
    echo "  ${DB_NAME}: not downloaded"
  fi
done
echo "Done."
```
