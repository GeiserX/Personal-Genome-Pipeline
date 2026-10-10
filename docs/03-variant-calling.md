# Step 3: Variant Calling (BAM to VCF and gVCF)

## What This Does
Identifies all positions where the sample's DNA differs from the reference genome: SNPs (single nucleotide changes) and small indels (insertions/deletions <50bp).

It writes two files:

- `vcf/${SAMPLE}.vcf.gz`: the variant sites only.
- `vcf/${SAMPLE}.g.vcf.gz` (a gVCF): the same calls plus "reference blocks", stretches where the sample matches the reference with a genotype quality. A position inside a block is a confirmed 0/0. A position outside any covered block was not covered. The variant-only VCF cannot tell those two apart. PharmCAT (step 7), the PRS (step 25) and imputation prep with panel sites (step 14) read hom-ref genotypes from the gVCF.

## Why
The VCF file is the foundation for ALL downstream analyses: ClinVar screening, pharmacogenomics, PRS, ROH, etc.

## Tool
- **DeepVariant** v1.10.0 — Google's deep learning variant caller (state-of-the-art accuracy)

## Docker Image
- `DEEPVARIANT_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- Sorted, indexed BAM file
- GRCh38 reference genome (FASTA + FAI)

## Command
```bash
./scripts/03-deepvariant.sh your_sample male    # or female; see "Sex chromosomes" below
```

What the script runs:

```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
THREADS=8

docker run --rm \
  --cpus ${THREADS} --memory 32g \
  -v ${GENOME_DIR}:/genome \
  -v "$PWD/assets/par_grch38.bed:/pgp/par_grch38.bed:ro" \
  "${DEEPVARIANT_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="/genome/${REF_FASTA}" \
    --reads=/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    --output_vcf=/genome/${SAMPLE}/vcf/${SAMPLE}.part.vcf.gz \
    --output_gvcf=/genome/${SAMPLE}/vcf/${SAMPLE}.part.g.vcf.gz \
    --intermediate_results_dir=/genome/${SAMPLE}/vcf/deepvariant_tmp \
    --sample_name="${SAMPLE}" \
    --num_shards=${THREADS} \
    --haploid_contigs=chrX,chrY \
    --par_regions_bed=/pgp/par_grch38.bed   # these two lines for a male sample only

# The script then renames the .part files (and their .tbi indexes) to
# ${SAMPLE}.vcf.gz and ${SAMPLE}.g.vcf.gz and removes deepvariant_tmp/.

# For WES data, use MODEL_TYPE=WES:
# MODEL_TYPE=WES ./scripts/03-deepvariant.sh your_sample male

# To call only some regions, set INTERVALS (space-separated, passed to --regions):
# INTERVALS="chr20:10000001-10500000" ./scripts/03-deepvariant.sh your_sample male

# A second BAM (BWA-MEM2, long reads) into its own directory, so vcf/ is kept:
# ALIGN_DIR=aligned_bwamem2 VCF_OUT_DIR=vcf_bwamem2 ./scripts/03-deepvariant.sh your_sample male

# Output: ~93MB VCF with ~5.5M total variants (~4.6M PASS), and the gVCF
```

Settings, all optional:

| Variable | Default | What it does |
|---|---|---|
| `THREADS` | 8 | CPUs of the container and DeepVariant's `--num_shards` |
| `DV_MEM` | `32g` | container memory |
| `ALIGN_DIR` | `aligned` | where the BAM is, inside the sample directory |
| `VCF_OUT_DIR` | `vcf` | where the VCF and gVCF go, inside the sample directory |
| `MODEL_TYPE` | `WGS` | `WGS`, `WES`, `PACBIO` or `ONT_R104` |
| `INTERVALS` | whole genome | regions to call |

DeepVariant writes its outputs under `.part` names, and the script renames them only when the VCF, the gVCF and both indexes are complete, so a killed run leaves no VCF that looks finished. Its intermediate files go to `deepvariant_tmp/` in the output directory, not into the container's own disk, and are removed at the end of a run.

## Sex chromosomes

A male sample has one X and one Y. Outside the pseudoautosomal regions (PARs, the ends of X and Y that pair with each other), a heterozygous call on chrX or chrY cannot be real. With `male`, DeepVariant calls chrX and chrY haploid outside the PARs listed in `assets/par_grch38.bed` (the GRCh38 PAR1 and PAR2 of chrX and chrY), so those calls are homozygous. `female` and no sex call every contig diploid; no sex prints a note saying so. `run-all.sh` passes the sex it was given.

The alternative callers (03a GATK, 03b FreeBayes, 03d Octopus) call every contig diploid, chrX and chrY of a male sample too. Clair3 (03e, long reads) takes the same `male` or `female` argument (`--gender`).

## Alternative callers: scatter

GATK HaplotypeCaller (03a) multithreads only its PairHMM and FreeBayes (03b) runs on one thread, so each would leave most of the CPUs idle for a whole genome. Both split the work into units: chr1-22, chrX, chrY and chrM one each and every other contig together, or each region of `INTERVALS` (space-separated, for example `INTERVALS="chr20 chr22"`). `SCATTER_JOBS` units (default THREADS/2) run at once, GATK with 2 CPUs and 8 GB each, FreeBayes with 1 CPU and 8 GB, and the parts are joined in reference order. `SCATTER=false` runs one process over everything (32 GB). On chr20 and chr22 of the e2e fixture the scattered and the single run of each caller give the same records (`tests/e2e/sv-mito-telomere-steps-6-scatter.sh`).

## Resource Requirements
- CPU: `THREADS` (default 8) sets `--cpus` and `--num_shards`. `make_examples` runs one process per shard; `call_variants` on CPU scales sub-linearly ([upstream](https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-details.md#call_variants)). `run-all.sh` runs DeepVariant with 8 shards on any machine; [CPU requirements](hardware-requirements.md#cpu-requirements) shows how to give it more
- RAM: 32GB recommended (`DV_MEM`)
- Disk: the intermediate files in `deepvariant_tmp/` need free space in the sample directory while the step runs
- GPU: neither entry point has a GPU option; only `call_variants` could use one (see [troubleshooting](troubleshooting.md#step-3-deepvariant-and-a-gpu))
- Time: see [Hardware and storage requirements](hardware-requirements.md#runtime-per-step)

## Output Interpretation
- **PASS** variants: high-confidence calls (~4.6M per 30X WGS sample)
- **RefCall**: site looks like reference (not a variant)
- **LowQual**: low confidence
- GQ (Genotype Quality): higher = more confident
- DP (Read Depth): typical 25-35x for a 30X WGS sample

## Notes
- DeepVariant is the gold standard for SNPs/indels but does NOT detect:
  - Structural variants >50bp (use Manta, step 4)
  - Repeat expansions (use ExpansionHunter, step 9)
  - Copy number variants (use CNVpytor or GATK gCNV)
- bcftools can also call variants but is significantly less accurate than DeepVariant
