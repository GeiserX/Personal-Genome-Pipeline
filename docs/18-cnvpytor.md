# Step 18: Depth-Based CNV Calling with CNVpytor

## What This Does
Detects copy number variants (CNVs) using read-depth analysis — complementary to Manta's paired-end/split-read approach. Especially effective for large CNVs (>1 kb) that Manta may miss.

CNVpytor is the maintained Python successor to CNVnator (same lab — Abyzov — and the same read-depth method), replacing the heavy CERN ROOT dependency with a portable HDF5 `.pytor` store.

## Why
Manta (step 4) detects SVs from discordant read pairs and split reads, which works well for balanced SVs (inversions, translocations) and smaller deletions/duplications. CNVpytor uses read-depth signal only, making it better for:
- Large duplications/deletions (>10 kb)
- Tandem duplications where breakpoints are ambiguous
- Validating Manta/Delly calls with orthogonal evidence
- Copy-number-neutral loss of heterozygosity (LOH), which paired-end callers cannot see

## Tool
- **CNVpytor** (Abyzov lab) — Python reimplementation of CNVnator (Abyzov et al., Genome Research 2011)

## Docker Image
```
quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0
```

> **Reference resources required.** The biocontainer ships **without** the GC/mask resource files and its built-in `-download` is broken in 1.3.2. Pinned resource files must be present at `${GENOME_DIR}/reference/cnvpytor/` and are bind-mounted into the container. See **[00-reference-setup.md](00-reference-setup.md)** for the one-time download.

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
IMG=quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0
# The container has no GC/mask data; mount the pinned resource dir onto its data path.
DATA=/usr/local/lib/python3.12/site-packages/cnvpytor/data
MOUNTS="-v ${GENOME_DIR}:/genome -v ${GENOME_DIR}/reference/cnvpytor:${DATA}"
PYTOR=/genome/${SAMPLE}/cnvpytor/${SAMPLE}.pytor
mkdir -p ${GENOME_DIR}/${SAMPLE}/cnvpytor

# 1. Import read depth from the BAM (no reference FASTA needed)
docker run --rm --user root ${MOUNTS} ${IMG} \
  cnvpytor -root ${PYTOR} -rd /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam -j 4

# 2. Read-depth histogram + automatic GC correction
docker run --rm --user root ${MOUNTS} ${IMG} cnvpytor -root ${PYTOR} -his 1000

# 3. Partition (mean-shift segmentation)
docker run --rm --user root ${MOUNTS} ${IMG} cnvpytor -root ${PYTOR} -partition 1000

# 4. Call CNVs (tab-separated)
docker run --rm --user root ${MOUNTS} ${IMG} cnvpytor -root ${PYTOR} -call 1000 \
  > ${GENOME_DIR}/${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.txt

# 5. Export a VCF (deletions / duplications / LOH). `-view` reads commands from stdin, so pass -i.
docker run --rm --user root -i ${MOUNTS} ${IMG} cnvpytor -root ${PYTOR} -view 1000 <<'VIEW'
set print_filename /genome/${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.raw.vcf
print calls
VIEW
```

`scripts/18-cnvpytor.sh` runs all of the above and then normalizes the VCF (reheaders the full contig set from the reference `.fai`, sorts, bgzips, and indexes) so it merges cleanly with the other SV callers in step 22.

## Bin Size
The `1000` parameter is the bin size in base pairs (must be divisible by 100). Use:
- `1000` for 30X WGS (recommended)
- `10000`/`100000` for coarser, faster genome-wide scans
- CNVpytor can compute several bin sizes at once: `-his 1000 10000 100000`

## Output
- `${SAMPLE}.pytor` — HDF5 read-depth store (intermediate; a few GB, can be deleted)
- `${SAMPLE}_cnvs.txt` — tab-separated CNV calls, 11 columns:
  1. type (`deletion` / `duplication`)
  2. region (`chr:start-end`)
  3. size (bp)
  4. read-depth level (normalized; 1.0 ≈ diploid)
  5–8. e-values (t-test / Gaussian tail, full and middle region)
  9. q0 (fraction of reads with mapping quality 0)
  10. pN (fraction of reference N bases)
  11. dG (distance to nearest gap >100 bp)
- `${SAMPLE}_cnvs.vcf.gz` (+ `.tbi`) — normalized VCF with `SVTYPE=DEL/DUP/LOH`, `END`, `SVLEN`, `GT`, and `CN`

## Filtering
```bash
# Keep confident CNVs (e-value in col 5 < 0.05, size in col 3 > 1 kb)
awk -F'\t' '$5 < 0.05 && $3 > 1000' ${SAMPLE}_cnvs.txt
```
The first call is often a low-`level` artifact spanning the chromosome-start N gap (`pN` near 1.0); the `pN`/`q0` columns let you drop such regions.

## Runtime
~1-3 hours per 30X WGS genome (CNVpytor is multi-threaded via `-j`).

## Notes
- Run AFTER alignment (step 2). Independent of Manta/Delly — can run in parallel.
- CNVpytor and Manta overlap on the majority of true large CNVs; variants called by both have lower false-positive rates than single-caller calls.
- No CERN ROOT dependency — the `.pytor` store is standard HDF5, so the container is small and reproducible.
- Reproducibility: the GC/mask resources are pinned to the CNVpytor v1.3.2 tag and mounted offline; no network access is needed at run time (see [00-reference-setup.md](00-reference-setup.md)).
- For maximum sensitivity, intersect CNVpytor calls with Manta and Delly (step 19) for a consensus call set (step 22).
- CNVpytor can also model B-allele frequency (`-snp`/`-baf`) for allele-specific CNV/LOH; this pipeline uses the read-depth path.
- Read-depth import is restricted to the canonical chromosomes (`chr1`–`chr22`, `chrX`, `chrY`). A full-reference GRCh38 BAM also carries hundreds of ALT/HLA/decoy contigs that the GC-correction data does not cover; without this restriction `-rd` chokes on them and produces no calls.
