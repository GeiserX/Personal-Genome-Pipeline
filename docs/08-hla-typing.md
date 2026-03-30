# Step 8: HLA Typing (T1K)

## What This Does
Determines your HLA genotype (Human Leukocyte Antigen) from WGS data — the immune system genes that control tissue compatibility and drug hypersensitivity reactions.

## Why
HLA alleles determine transplant compatibility, predisposition to autoimmune diseases, and severe adverse drug reactions (e.g., HLA-B*57:01 and abacavir, HLA-B*58:01 and allopurinol).

## Tool
- **T1K** v1.0.9 — efficient HLA genotyping from sequencing reads

## Docker Image
```
quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0
```

## Prerequisites
- Aligned BAM from step 2
- Pre-built HLA reference index from step 00 (`t1k_idx/hlaidx_grch38`)

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

mkdir -p ${GENOME_DIR}/${SAMPLE}/hla

docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
  run-t1k \
    -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    -f /genome/t1k_idx/hlaidx_grch38_rna_seq.fa \
    --preset hla \
    -o /genome/${SAMPLE}/hla/${SAMPLE}_t1k \
    -t 4
```

## Output
- `${SAMPLE}_t1k_genotype.tsv` — HLA allele calls per locus (A, B, C, DRB1, DQB1, etc.)
- Two alleles per locus (one per chromosome)

## Alternative: HLA-LA
For a second opinion or when T1K results are ambiguous:

```bash
docker run --rm \
  --cpus 8 --memory 16g \
  -v ${GENOME_DIR}:/genome \
  jiachenzdocker/hla-la:latest \
  HLA-LA.pl \
    --BAM /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    --graph PRG_MHC_GRCh38_withIMGT \
    --sampleID ${SAMPLE} \
    --maxThreads 8 \
    --workingDir /genome/${SAMPLE}/hla_la
```

- Docker image is 4.5GB (includes pre-built graph)
- Slower but uses a different algorithm — useful for validation

## Key HLA Alleles for Drug Safety
| Allele | Drug | Risk |
|---|---|---|
| HLA-B*57:01 | Abacavir (HIV) | Severe hypersensitivity reaction |
| HLA-B*58:01 | Allopurinol (gout) | Stevens-Johnson syndrome / TEN |
| HLA-B*15:02 | Carbamazepine | Stevens-Johnson syndrome (SE Asian) |
| HLA-A*31:01 | Carbamazepine | Drug reaction with eosinophilia |
| HLA-B*57:01 | Flucloxacillin | Drug-induced liver injury |

## Important Notes
- HLA typing from WGS is **approximate** — clinical HLA typing for transplant or critical drug decisions uses dedicated high-resolution panels (sequence-based typing)
- WGS-based HLA is sufficient for pharmacogenomic screening (presence/absence of risk alleles)
- T1K requires the pre-built index from step 00 — do not skip the `t1k-build.pl` step
- Running both T1K and HLA-LA and comparing results increases confidence in the calls
- HLA region is the most polymorphic in the human genome — ambiguous calls are expected for rare alleles
