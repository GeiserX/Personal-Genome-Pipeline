# Realigning After a Reference Change

The pipeline aligns to NCBI's **GRCh38 no-ALT analysis set** (`reference/GRCh38_no_alt_analysis_set.fasta`, 195 sequences). Versions before it used the Broad `hg38` FASTA (`reference/Homo_sapiens_assembly38.fasta`, 3,366 sequences with ALT, HLA and decoy contigs). A BAM aligned to the old file, or to any other GRCh38 file such as a vendor's, does not match the new one, and every result made from that BAM has to be made again.

This page says why, what to redo, what to keep, how to do it one sample at a time on your own machine, and how to check that a BAM is on the new reference. `validate-setup.sh` refuses a BAM from another reference, so a run cannot mix the two by accident.

## Why the reference changed

ALT contigs are second copies of regions that vary a lot between people: the MHC (the HLA genes), KIR, the CYP2D6 locus and others. An aligner that is run ALT-aware treats the two copies as one place. None of the aligners here is (minimap2 in step 02, BWA-MEM2 in step 02a). On a reference with ALT contigs, a read that matches the primary copy and the ALT copy equally well is placed on one of them at random with mapping quality (MAPQ) 0, and the callers skip MAPQ 0 reads. Depth then thins at exactly the loci where copy number and HLA type matter most, and a depth-based caller (pypgx, Cyrius, CNVpytor) can report a deletion that is not there.

### How much depth ALT contigs cost

The `ALT depth A/B` workflow ([`scripts/ci/alt-depth-ab.sh`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/ci/alt-depth-ab.sh)) measures it on the test fixture, a 30x slice of the public GIAB HG002 genome. It maps the fixture's CYP2D6-region and MHC-region read pairs to both whole references with step 02's aligner and preset, and reads the mean depth with mosdepth, for all reads and for MAPQ >= 1, the reads a caller uses:

| Region | Reference | Depth, all reads | Depth, MAPQ >= 1 | MAPQ >= 1 / all |
|---|---|---|---|---|
| CYP2D6 | no-ALT | 18.5 | 16.8 | 0.91 |
| CYP2D6 | with-ALT | 4.1 | 0.0 | 0.00 |
| CYP2D7 | no-ALT | 23.8 | 22.1 | 0.93 |
| CYP2D7 | with-ALT | 5.9 | 0.1 | 0.03 |
| CYP2D flanks | no-ALT | 29.7 | 29.7 | 1.00 |
| CYP2D flanks | with-ALT | 15.0 | 8.1 | 0.54 |
| HLA-A | no-ALT | 26.9 | 26.9 | 1.00 |
| HLA-A | with-ALT | 0.7 | 0.0 | 0.00 |
| HLA-B | no-ALT | 26.9 | 26.9 | 1.00 |
| HLA-B | with-ALT | 1.0 | 0.2 | 0.18 |

| Reference | Read set | Primary alignments | On ALT or HLA contigs |
|---|---|---|---|
| no-ALT | CYP2D slice | 60,823 | 0 |
| no-ALT | MHC slice | 635,211 | 0 |
| with-ALT | CYP2D slice | 60,820 | 24,386 (40%) |
| with-ALT | MHC slice | 635,805 | 542,204 (85%) |

On the old reference a caller saw almost nothing at CYP2D6, CYP2D7, HLA-A and HLA-B (0.0 to 0.2x at MAPQ >= 1), where the new one gives 17 to 27x. The reads were still aligned, but to the ALT and HLA contigs: 85% of the MHC slice's alignments went there, because the old reference has several more copies of the MHC as ALT contigs plus the HLA allele contigs, and a read that fits all of them equally is placed on one at random. Even the flanks lost half their depth, so the old reference has a second copy of that stretch of chr22 too. T1K (step 8) takes a BAM's HLA reads from the positions in its coordinate file, which are on chr6, so on the old reference it saw a small fraction of them. On the new reference CYP2D6 and CYP2D7 keep about 90% of their reads at MAPQ >= 1; the rest match CYP2D6, CYP2D7 and the CYP2D8 pseudogene equally, and no aligner can tell those apart.

CYP2D6 and CYP2D7 are the regions Cyrius (step 21) reads; the flanks are two 50 kb stretches of the same slice outside the genes, chr22:42.05-42.10 Mb and 42.20-42.25 Mb; HLA-A and HLA-B are the gene bodies. The workflow runs from the Actions tab (`ALT depth A/B`, "Run workflow") and writes the same tables to its job summary.

## What has to be redone

Everything made from the BAM: step 02 (alignment) and every step after it. In practice that is the whole sample directory except the FASTQ: the BAM, the VCF and gVCF of step 3, the SV and CNV calls, HLA types, the pharmacogenomic calls (steps 7, 21, 27 and 32), annotations, PRS, ancestry and the reports.

Two kinds of sample need little or nothing:

- **A sample that starts from a vendor VCF** (no reads). Its positions on chr1-22, X, Y and M are the same on both references, so the analysis steps run on it as before.
- **Chip data** only needs `scripts/chip-to-vcf.sh` once more (seconds), so the VCF header lists the new reference's contigs.

## What stays

None of these depend on ALT contigs, so they stay as they are:

- the FASTQ files;
- ClinVar, the VEP caches, the PCGR/CPSR bundle and the AnnotSV data;
- the annotation score files (CADD, SpliceAI, REVEL, AlphaMissense, gnomAD constraint) and the PGS scoring files;
- the pypgx bundle, the ancestry panel and the GIAB truth sets (GIAB defines v4.2.1 on this same reference);
- Delly's exclude map, the chromosome bands, CNVpytor's GC and mask files (chr1-22, X and Y only);
- the IPD-IMGT/HLA release, the GENCODE gene file and the T1K index built from them. T1K's coordinate file comes from those two files only (`t1k-build.pl -d hla.dat -g genes.gtf`, step 8); it never reads the reference FASTA, so it is the same on both references and is not rebuilt.

## Indexes

| Index | On the new reference | Notes |
|---|---|---|
| `.fai` | `setup.sh` downloads NCBI's and checks its md5 | 195 lines |
| `.dict` | `setup.sh` builds it | GATK, Picard and `chip-to-vcf.sh` read it |
| minimap2 `.sr.mmi` | step 02 builds it the first time it runs | about 30 minutes, with a 32 GB container cap; named after the reference, `GRCh38_no_alt_analysis_set.sr.mmi`, so the old index is never used |
| BWA-MEM2 (step 02a only) | by hand | about 90 GB of RAM for the build: [BWA-MEM2 index](00-reference-setup.md#bwa-mem2-index) |
| classic BWA (GRIDSS, TIDDIT's assembly) | by hand | `bwa index`, about an hour: [GRIDSS](04b-gridss.md) |

The old reference's files (`Homo_sapiens_assembly38.fasta`, `.fai`, `.dict`, `.sr.mmi` and any BWA index beside it) are no longer read. Delete them once the realigned samples check out; together they take about 12 GB.

## Step by step, one sample

Run these on the machine that holds your data, in the repository folder. The realignment itself is a normal step 02 run: 1-2 hours for a 30x genome, plus the one-time index build.

**1. Install the new reference** next to the old one:

```bash
export GENOME_DIR=/path/to/your/data
./scripts/setup.sh ${GENOME_DIR}
./scripts/validate-setup.sh
```

`setup.sh` prints `[OK] Reference: .../GRCh38_no_alt_analysis_set.fasta (195 sequences)`, and `validate-setup.sh` prints `[OK] Reference has no ALT or HLA contigs (195 sequences)`.

**2. Move the old results aside.** Nothing is deleted, so you can compare old and new results later:

```bash
export SAMPLE=your_name
OLD="${SAMPLE}.old-reference"
mv "${GENOME_DIR}/${SAMPLE}" "${GENOME_DIR}/${OLD}"
mkdir -p "${GENOME_DIR}/${SAMPLE}/fastq"
```

**3. Get the reads back as FASTQ.** If you still have the FASTQ files, move them in:

```bash
mv "${GENOME_DIR}/${OLD}/fastq/${SAMPLE}_R1.fastq.gz" "${GENOME_DIR}/${OLD}/fastq/${SAMPLE}_R2.fastq.gz" \
  "${GENOME_DIR}/${SAMPLE}/fastq/"
```

### From an existing BAM

If the BAM is all you have (a vendor BAM, or the FASTQ was deleted), turn it back into paired FASTQ. `samtools collate` groups the two reads of each pair, and `samtools fastq` writes them, skipping secondary and supplementary records, so every read appears once. It needs free disk of about the BAM's size for its temporary files, and takes about an hour for a 30x BAM:

```bash
source versions.env   # from the repository root
docker run --rm --user "$(id -u):$(id -g)" --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "${SAMTOOLS_IMAGE}" \
  bash -c "set -euo pipefail
    samtools collate -u -O -@ 4 /genome/${OLD}/aligned/${SAMPLE}_sorted.bam /genome/${SAMPLE}/fastq/collate \
    | samtools fastq -n -@ 4 \
        -1 /genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz \
        -2 /genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz \
        -0 /dev/null -s /dev/null -"
```

Check that no pair was lost: the two numbers must be equal. The first counts the read pairs in the BAM (first reads of primary records), the second the reads in R1 (several minutes for a 30x file):

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" "${SAMTOOLS_IMAGE}" \
  samtools view -c -f 0x40 -F 0x900 /genome/${OLD}/aligned/${SAMPLE}_sorted.bam
echo $(( $(gzip -dc "${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz" | wc -l) / 4 ))
```

A read whose mate is not in the BAM at all is left out (`-s /dev/null`), so a vendor BAM that was filtered can give a slightly smaller second number. The reads keep any trimming done before the old alignment; `SKIP_TRIM=true` in the next step skips fastp's second pass over them.

**4. Run the pipeline** as for a new sample. Step 02 builds the new minimap2 index on its first run:

```bash
./scripts/run-all.sh ${SAMPLE} <male|female>
```

**5. Check that the BAM is on the new reference** (next section). When it is and the results look right, delete `${GENOME_DIR}/${OLD}` and the old reference files.

## Checks that prove a BAM is on the new reference

Each one reads the BAM step 02 wrote, `${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam`:

```bash
source versions.env   # from the repository root
BAM=/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
sam() { docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" "${SAMTOOLS_IMAGE}" samtools "$@"; }

./scripts/validate-setup.sh ${SAMPLE}                     # [OK] BAM header matches the reference: the same 195 sequences in the same order
sam view -H "$BAM" | grep -c '^@SQ'                       # 195
sam view -H "$BAM" | grep -c -E 'SN:(chr[^[:space:]]*_alt|HLA-)'   # 0
sam view -H "$BAM" | grep '^@PG' | grep -o 'GRCh38_no_alt_analysis_set[^ ]*'   # the minimap2 line names GRCh38_no_alt_analysis_set.sr.mmi
```

`validate-setup.sh` is the one that matters: it compares the BAM's sequence names and lengths with the reference's `.fai`, in order, and on a BAM from another reference fails with `this BAM was aligned to a different reference: realign (docs/realignment.md)` and the first sequence that differs. It also checks that the BAM is coordinate-sorted and passes `samtools quickcheck`. `run-all.sh` runs it before any step.

The VCF made from the new BAM names no ALT contig either:

```bash
source versions.env   # from the repository root
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
  bcftools view -h /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz | grep -c '^##contig=<ID=[^,]*_alt,'   # 0
```

## Keeping the old reference on purpose

To keep analysing with the Broad file, for example to finish a comparison, point `REF_FASTA` at it and allow its ALT contigs:

```bash
export REF_FASTA=reference/Homo_sapiens_assembly38.fasta
export ALLOW_ALT_REFERENCE=true
./scripts/validate-setup.sh ${SAMPLE}   # [WARN] Reference has ... ALT/HLA contigs (allowed by ALLOW_ALT_REFERENCE=true)
```

The depth loss above then applies to everything that sample produces. Unset both variables to go back to the default.

## References considered and not chosen

The decoy variant of the same NCBI set and GIAB's GRCh38 file that also masks false duplications were considered. [Reference setup](00-reference-setup.md#why-the-no-alt-analysis-set) says why neither is the default and how to use the decoy variant through `REF_FASTA`; it has no ALT or HLA contigs, so `validate-setup.sh` accepts it as it is.
