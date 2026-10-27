# Step 50: Mobile Element Insertions (Alu, L1, SVA)

## What This Does

Finds the copies of Alu, L1 (LINE-1) and SVA elements that the genome carries and the reference does not, and genotypes each one as one or two copies, with MEGAnE ([Kojima et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37169872/), Nat Genet). Then it reads a short catalog of insertions known to cause disease as seen, not seen, or not callable. Opt-in.

Mobile elements are pieces of DNA that copy themselves into new places. In the 1000 Genomes Project each genome carried a median of about 1,200 insertions that the reference lacks, most of them Alu ([Sudmant et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26432246/)), nearly all common and harmless. A few land in a gene and break it:

| Gene | Insertion | Condition | Inheritance | ClinVar | Source |
|---|---|---|---|---|---|
| FKTN | SVA in the 3' UTR | Fukuyama congenital muscular dystrophy | autosomal recessive | [3199](https://www.ncbi.nlm.nih.gov/clinvar/variation/3199/) | [Kobayashi et al. 1998](https://pubmed.ncbi.nlm.nih.gov/9690476/) |
| TAF1 | SVA in intron 32 | X-linked dystonia-parkinsonism | X-linked | [9862](https://www.ncbi.nlm.nih.gov/clinvar/variation/9862/) | [Aneichyk et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29474918/) |
| BBS1 | SVA-F in exon 13 (c.1214_1215insSVA) | Bardet-Biedl syndrome 1 | autosomal recessive | [981241](https://www.ncbi.nlm.nih.gov/clinvar/variation/981241/) | [Delvallée et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33169370/) |

The FKTN insertion is a founder variant: Kobayashi et al. found it on every patient chromosome with the founder haplotype, 87% of those they studied. Delvallée et al. found the BBS1 insertion in eight families.

## Why

Manta, Delly and the other SV callers see an insertion as a break in the alignment, but they do not say what was inserted, and an SVA, about 2 kb of repeat, is much longer than the reads. A caller built for mobile elements uses what such an insertion leaves in short reads: pairs with one end in the element and the other in unique sequence nearby, and reads that cross the element's edge. MEGAnE does both, and genotypes each insertion from the depth of that evidence.

## Tool

- MEGAnE v1.0.1.beta (MIT licence), `call_genotype_38 -only_ins`: insertions only, not the reference copies a sample lacks
- samtools to list the chromosomes the BAM has and to measure the depth over the catalog windows
- `bin/mei_report.py` (standard library Python) for the calls table and the catalog

## Docker Image

- `MEGANE_IMAGE`: MEGAnE's only image (amd64 only, about 6.4 GB to download, 19.5 GB unpacked), its last release, pinned to its digest. It carries its own GRCh38 RepeatMasker annotation and BLAST database.
- `SAMTOOLS_IMAGE` and `PYTHON_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. `MEGANE_IMAGE` is an `# optional` line: `setup.sh --megane-kmers` pulls it.

## Reference data

`./scripts/setup.sh --megane-kmers <genome_dir>` builds MEGAnE's k-mer set of the reference: every 32-mer seen more than once in it, which MEGAnE uses to drop reads that map to many places. It is written to `reference/megane/<reference name>.mk` with its count file `.mi` beside it (`MEGANE_KMERS` overrides the path, inside `GENOME_DIR`). The builder holds 16 bytes per reference base in memory: about 50 GB for GRCh38, which MEGAnE's authors give as 48 GB and about 10 minutes. `setup.sh` checks the memory Docker has and stops with a message when it is short; then build it on a larger machine and copy `reference/megane/` over. See [Reference data setup](00-reference-setup.md#megane-k-mer-set-optional-for-step-50).

The catalog is `assets/mei_sites.tsv`: one row per insertion, with its GRCh38 site (the base after which the element is inserted, from the ClinVar record), a window of the site -/+ 500 bp, the element, the inheritance, the condition, the ClinVar id and the PubMed id. Checked on 2026-10-27 against ClinVar, PubMed and Ensembl.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps), paired-end short reads. A VCF alone has no reads, so a sample without a BAM is left out.
- the k-mer set above
- the sex (`male` or `female`), optional: MEGAnE uses it to genotype chrX and chrY

## Command

```bash
./scripts/50-mei.sh your_name male
```

With `run-all.sh`: `TOOLS=...,mei`; a run lists step 50 only when `TOOLS` names it. In the Nextflow pipeline: `mei` in `--tools` with `--megane_kmers reference/megane/GRCh38_no_alt_analysis_set.mk` (processes `MEI_DEPTH`, `MEI` and `MEI_REPORT`).

## How it works

1. samtools lists the chromosomes MEGAnE reads (chr1 to chr22, chrX and chrY, as the BAM names them) and measures the depth over the catalog windows on the contigs the BAM has (`samtools bedcov -c`).
2. MEGAnE runs on the BAM with that chromosome list. MEGAnE's own list holds every GRCh38 contig, ALT contigs included, and it stops on a BAM that lacks one; a BAM aligned to the no-ALT analysis set lacks all of them. On a BAM without chrX the step adds `-no_sex_chr`, or MEGAnE stops. MEGAnE checks that the `.fai` is newer than the FASTA, so it reads the reference through a link with a fresh copy of the `.fai`.
3. MEGAnE writes `MEI_final_percentile_genotyped.vcf`, the file its log recommends for one sample. When it finds too few insertions to fit its filter thresholds (a BAM of one chromosome), it writes `MEI_final_failed_genotyped.vcf` with fixed thresholds instead; the step then uses that file and prints a warning.
4. `bin/mei_report.py` keeps the Alu, L1 and SVA insertions that MEGAnE identified (`MEPRED=PASS`) with one or two copies. Other repeat classes in MEGAnE's VCF (DNA transposons, old LTRs) are left out: none is active in human genomes, so a call of one is read as an artefact.
5. Each catalog row is **seen** when a kept insertion of the same element lies in its window, **not callable** when the BAM has no reads of that contig or the mean depth over the window is under 10, and **not seen** otherwise. The 10x is this pipeline's choice, not a published value.

## Output

| File | Contents |
|---|---|
| `mei/${SAMPLE}_mei.tsv` | one row per insertion: `chrom`, `pos`, `end` (MEGAnE's right breakpoint), `element` (Alu, L1, SVA), `subfamily` (for example AluYa5, L1HS, SVA_F), `copies` (1 or 2), `filter` (MEGAnE's FILTER: `PASS`, or flags such as `LC`, low confidence), `confidence` (`H` or `L`, how sure the copy number is) |
| `mei/${SAMPLE}_mei_sites.tsv` | the catalog with `status` (seen, not seen, not callable), `mean_depth` over the window and `call` (the insertion seen: position, subfamily, copies, filter) |
| `mei/${SAMPLE}_megane.vcf` | MEGAnE's VCF as it wrote it |
| `mei/${SAMPLE}_mei_depth.tsv` | `samtools bedcov -c` over the catalog windows |
| `mei/${SAMPLE}_megane.log` | MEGAnE's log |

Both reports show the count per element and the catalog in their "Mobile element insertions" section. Without step 50, the report lists mobile element insertions under what the run did not assess.

## Runtime

On HG002's chr20 at 30x (13 million reads) on a 4-CPU GitHub runner, MEGAnE took 23 minutes with the unmapped reads in the BAM, 21 of them in its search of the unmapped reads against the whole genome. MEGAnE's authors give about an hour for a 30x genome on 4 threads. The step asks for 16 GB of memory for MEGAnE. Not measured here on a whole genome.

## Validation

The e2e case `tests/e2e/mei-1.sh` runs the step on HG002's chr20: the reads of chr20 in GIAB's 60x novoalign BAM, half the pairs (30x), realigned with minimap2 to chr20 alone, and the k-mer set of that chromosome. Its Alu calls are scored against the GIAB v5.0q insertions of HG002 on chr20 (RepeatMasker SINE/Alu, interspersed, 250 to 400 bp, inside the benchmark regions): recall must be at least 0.70 and precision at least 0.80, within 100 bp. With every read pair near a GIAB Alu insertion removed, the recall check fails. `tests/test_mei_report.py` plants an SVA in the FKTN window of a synthetic MEGAnE VCF and checks it is seen, and not seen without it.

## Notes

- Not seen is not proof of absence. An SVA is about 2 kb of repeat, and short reads find its edges, not its body; a low-complexity or repeat-rich window can hide one. A finding is confirmed in a clinical laboratory, usually by long-range PCR across the site.
- Carrying one copy of the FKTN or BBS1 insertion means carrier status for a recessive condition, not the disease. The TAF1 insertion is on chrX: one copy in a male means the X he has carries it.
- The e2e case scores Alu calls only: chr20 has 15 GIAB Alu insertions that pass its filter, and too few L1 or SVA insertions to count. The L1 and SVA calls are listed but not measured here.
- MEGAnE's last release is a 2022 beta. Its image is amd64 only, so on an Apple Silicon Mac Docker runs it under emulation (not tested here).
- The step needs paired-end short reads. It has not been tested on long reads.
