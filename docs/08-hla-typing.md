# Step 8: HLA Typing (T1K)

## What This Does
Determines your HLA genotype (Human Leukocyte Antigen) from WGS data — the immune system genes that control tissue compatibility and drug hypersensitivity reactions.

## Why
HLA alleles determine transplant compatibility, predisposition to autoimmune diseases, and severe adverse drug reactions (e.g., HLA-B*57:01 and abacavir, HLA-B*58:01 and allopurinol).

## Tool
- **T1K** v1.0.9 — efficient HLA genotyping from sequencing reads

## Docker Image
- `T1K_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- Aligned BAM from step 2 (or set `ALIGN_DIR`, for example `ALIGN_DIR=aligned_bwamem2`)
- The IPD-IMGT/HLA release `HLA_DB_RELEASE` (3.65.0) and the GENCODE gene lines, both installed by `setup.sh` (see [reference setup](00-reference-setup.md#small-pinned-data-files)). Without them the step prints `SKIPPED` and the command that installs them.
- The step builds its T1K index on the first run, into `t1k_idx/t1k-<T1K version>_imgt-<release>_gencode-<release>/`. A new T1K image or another `HLA_DB_RELEASE` builds a new index.

## Command
```bash
export GENOME_DIR=/path/to/your/data
./scripts/08-hla-typing.sh your_sample
```

`THREADS` (default 4) sets the container's CPUs and T1K's `-t`.

The T1K call the script makes. It uses the DNA index with its coordinate file and the `hla-wgs` preset, which are the right inputs for whole-genome DNA reads. `IDX` is the index directory the script built:

```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
IDX=$(ls -d ${GENOME_DIR}/t1k_idx/t1k-*_imgt-3.65.0_gencode-* | head -n 1)
IDX=/genome/${IDX#${GENOME_DIR}/}

mkdir -p ${GENOME_DIR}/${SAMPLE}/hla_t1k

docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  "${T1K_IMAGE}" \
  run-t1k \
    -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    -f ${IDX}/hla_dna_seq.fa \
    -c ${IDX}/hla_dna_coord.fa \
    --preset hla-wgs \
    -t 4 \
    --od /genome/${SAMPLE}/hla_t1k/ \
    -o ${SAMPLE}_hla
```

## Output
- `${SAMPLE}/hla_t1k/${SAMPLE}_hla_genotype.tsv` — HLA allele calls per locus (A, B, C, DRB1, DQB1, etc.)
- Two alleles per locus (one per chromosome)
- `${SAMPLE}/hla_t1k/database_release.txt` — the IPD-IMGT/HLA release the calls come from (read from `hla.dat`), the T1K version and the GENCODE release of the gene positions

HLA-A and HLA-B go on to PharmCAT: [step 36](36-pgx-consensus.md) cuts them to two fields (`*57:01`) and writes them as PharmCAT's outside calls, so step 7 gives the drug guidance for them (abacavir, allopurinol, carbamazepine, oxcarbazepine, phenytoin and others) and step 27 lists it. When either allele of HLA-A or HLA-B has a T1K quality of 0 or below, step 36 passes neither allele of that gene: the whole gene is withheld from PharmCAT. The report shows each allele's quality, marks such a gene low confidence and says it was not passed to PharmCAT.

## KIR Genes (opt-in)

`KIR=true` adds a second T1K pass over the KIR genes (killer-cell immunoglobulin-like receptors, chr19), with T1K's `kir-wgs` preset, against the IPD-KIR release `KIR_DB_RELEASE` names in `versions.env` (2.15.0). Install the database once:

```bash
./scripts/setup.sh --kir-data "$GENOME_DIR"     # kir/IPD-KIR_2.15.0/kir.dat, md5-checked
KIR=true ./scripts/08-hla-typing.sh your_sample
```

The index is built once per T1K version, IPD-KIR release and GENCODE release (`t1k_idx/t1k-<version>_kir-<release>_gencode-<release>`). Several KIR genes (KIR2DL2, KIR2DL5, KIR2DS1-3, KIR2DS5, KIR3DS1) are not on the GRCh38 primary assembly and get no coordinates; T1K types them from the reads of the genes that are. Output, in `${SAMPLE}/kir_t1k/`:

- `${SAMPLE}_kir_genotype.tsv` — the KIR genotypes in T1K's layout; when the BAM has too few reads at the KIR genes to type any of them, one line that says so instead;
- `database_release.txt` — the IPD-KIR release (read from `kir.dat`), the T1K version and the GENCODE release.

With Nextflow: `--kir true --kir_dat ${GENOME_DIR}/kir/IPD-KIR_2.15.0/kir.dat` with `hla_typing` in `--tools` (processes `KIR_BUILD` and `KIR_TYPING`); with `run-all.sh`: `KIR=true`. No step interprets the KIR genotypes yet: they matter for transplant matching and some immune conditions, not for drug guidance.

## Alternative: HLA-LA

> **Known issue:** in this pipeline's tests the image below crashes during graph alignment, and HLA-LA is still unsolved (see [Troubleshooting](troubleshooting.md#hla-typing-step-8-known-difficulties)). The command is kept for reference; use T1K for results.

For a second opinion or when T1K results are ambiguous:

```bash
docker run --rm \
  --cpus 8 --memory 16g \
  -v ${GENOME_DIR}:/genome \
  jiachenzdocker/hla-la@sha256:ecca23de6635aa85e60b4ee39dd4e15341b5febb514e5478f2b2a086f05a447c \
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
- `scripts/08-hla-typing.sh` runs `t1k-build.pl` itself when the index for its T1K version and database release is missing (a few minutes). The gene positions come from GENCODE's annotation, not from the FASTA: built from the FASTA or its `.fai`, every gene gets `-1 -1` coordinates and T1K extracts no reads. The step stops when one of the six typed genes has no coordinates
- Running both T1K and HLA-LA and comparing results increases confidence in the calls
- HLA region is the most polymorphic in the human genome — ambiguous calls are expected for rare alleles
