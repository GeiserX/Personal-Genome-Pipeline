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
- AnnotSV, when its data is there, for the genes each insertion falls in, where in them, and the benign set's frequency of an insertion at that place

## Docker Image

- `MEGANE_IMAGE`: MEGAnE's only image (amd64 only, about 6.4 GB to download, 19.5 GB unpacked), its last release, pinned to its digest. It carries its own GRCh38 RepeatMasker annotation and BLAST database.
- `SAMTOOLS_IMAGE` and `PYTHON_IMAGE`
- `ANNOTSV_IMAGE`, the image of [step 05](05-annotsv.md), when AnnotSV's data is there

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. `MEGANE_IMAGE` is an `# optional` line: `setup.sh --megane-kmers` pulls it.

## Reference data

`./scripts/setup.sh --megane-kmers <genome_dir>` builds MEGAnE's k-mer set of the reference: every 32-mer seen more than once in it, which MEGAnE uses to drop reads that map to many places. It is written to `reference/megane/<reference name>.mk` with its count file `.mi` beside it (`MEGANE_KMERS` overrides the path, inside `GENOME_DIR`). The builder holds 16 bytes per reference base in memory: about 50 GB for GRCh38, which MEGAnE's authors give as 48 GB and about 10 minutes. `setup.sh` checks the memory Docker has and stops with a message when it is short; then build it on a larger machine and copy `reference/megane/` over. See [Reference data setup](00-reference-setup.md#megane-k-mer-set-optional-for-step-50).

The catalog is `assets/mei_sites.tsv`: one row per insertion, with its GRCh38 site (the base after which the element is inserted, from the ClinVar record), a window of the site -/+ 500 bp, the element, the inheritance, the condition, the ClinVar id and the PubMed id. Checked on 2026-10-27 against ClinVar, PubMed and Ensembl.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps), paired-end short reads. A VCF alone has no reads, so a sample without a BAM is left out.
- the k-mer set above
- the sex (`male` or `female`), optional: MEGAnE uses it to genotype chrX and chrY
- AnnotSV's data, optional: `annotsv_annotations/` in `GENOME_DIR`, as `setup.sh` installs it for [step 05](05-annotsv.md), with the MANE Select list and gnomAD-SV's common SVs (`setup.sh --annotsv-data`). The step stops when the folder is there without those two files. Without the folder, the calls get no gene or frequency.

## Command

```bash
./scripts/50-mei.sh your_name male
```

With `run-all.sh`: `TOOLS=...,mei`; a run lists step 50 only when `TOOLS` names it. In the Nextflow pipeline: `mei` in `--tools` with `--megane_kmers ${GENOME_DIR}/reference/megane/GRCh38_no_alt_analysis_set.mk` (Nextflow reads a relative path from the launch directory; `run-all.sh` passes it for you) (processes `MEI_DEPTH`, `MEI`, `MEI_ANNOTSV` and `MEI_REPORT`). `MEI_ANNOTSV` runs when `--annotsv_annotations` is set, which `run-all.sh` does when AnnotSV's data is there and `annotsv` is in `TOOLS`.

## How it works

1. samtools lists the chromosomes MEGAnE reads (chr1 to chr22, chrX and chrY, as the BAM names them) and measures the depth over the catalog windows on the contigs the BAM has (`samtools bedcov -c`).
2. MEGAnE runs on the BAM with that chromosome list. MEGAnE's own list holds every GRCh38 contig, ALT contigs included, and it stops on a BAM that lacks one; a BAM aligned to the no-ALT analysis set lacks all of them. On a BAM without chrX the step adds `-no_sex_chr`, or MEGAnE stops. MEGAnE checks that the `.fai` is newer than the FASTA, so it reads the reference through a link with a fresh copy of the `.fai`.
3. MEGAnE writes `MEI_final_percentile_genotyped.vcf`, the file its log recommends for one sample. When it finds too few insertions to fit its filter thresholds (a BAM of one chromosome), it writes `MEI_final_failed_genotyped.vcf` with fixed thresholds instead; the step then uses that file and prints a warning.
4. `bin/mei_report.py` keeps the Alu, L1 and SVA insertions that MEGAnE identified (`MEPRED=PASS`) with one or two copies. Other repeat classes in MEGAnE's VCF (DNA transposons, old LTRs) are left out: none is active in human genomes, so a call of one is read as an artefact.
5. With AnnotSV's data, `bin/mei_report.py --annotsv-input` writes the kept insertions as the VCF AnnotSV reads, one `<INS:ME>` record each with `END` at its position and an ID the table is joined on, and AnnotSV annotates them as [step 05](05-annotsv.md) does the other SVs (`-annotationMode both`, the MANE Select transcripts). MEGAnE gives no inserted length, so AnnotSV takes each insertion as a point and skips its minimum-size check.
6. Each catalog row is **seen** when a kept insertion of the same element lies in its window, **not callable** when the BAM has no reads of that contig or the mean depth over the window is under 10, and **not seen** otherwise. The 10x is this pipeline's choice, not a published value.

## Output

| File | Contents |
|---|---|
| `mei/${SAMPLE}_mei.tsv` | one row per insertion: `chrom`, `pos`, `end` (MEGAnE's right breakpoint), `element` (Alu, L1, SVA), `subfamily` (for example AluYa5, L1HS, SVA_F), `copies` (1 or 2), `filter` (MEGAnE's FILTER: `PASS`, or flags such as `LC`, low confidence), `confidence` (`H` or `L`, how sure the copy number is), and from AnnotSV: `genes` (the genes the insertion falls in), `location` (each gene with where in it: `exon3`, `intron2-intron2`, ...), `benign_af` (`B_ins_AFmax`, the highest allele frequency of a benign insertion at that place in AnnotSV's benign set, gnomAD-SV among them) and `acmg_class` (AnnotSV's `ACMG_class`). All four are `.` without AnnotSV's data, and for a call AnnotSV did not annotate |
| `mei/${SAMPLE}_mei_annotsv.tsv` | AnnotSV's own table of those insertions, with AnnotSV's data only |
| `mei/${SAMPLE}_mei_sites.tsv` | the catalog with `status` (seen, not seen, not callable), `mean_depth` over the window and `call` (the insertion seen: position, subfamily, copies, filter) |
| `mei/${SAMPLE}_megane.vcf` | MEGAnE's VCF as it wrote it |
| `mei/${SAMPLE}_mei_depth.tsv` | `samtools bedcov -c` over the catalog windows |
| `mei/${SAMPLE}_megane.log` | MEGAnE's log |

Both reports show the count per element and the catalog in their "Mobile element insertions" section; the AnnotSV columns are in the calls table only.

AnnotSV classifies deletions and duplications only. In AnnotSV 3.5 the ranking of an insertion is empty (`SVrankingINS` in its `AnnotSV-ranking.tcl` scores nothing), so `acmg_class` is `NA` for every insertion. Read an insertion by its gene and location instead: one inside an exon can break the gene; most in introns or between genes do nothing, though some do (the TAF1 SVA above is intronic). A `benign_af` above a few percent means insertions at that place are common. Without step 50, the report lists mobile element insertions under what the run did not assess.

## Runtime

On HG002's chr20 at 30x (13 million mapped reads) on a 4-CPU GitHub runner, MEGAnE took 23 minutes with the unmapped reads in the BAM, 21 of them in its search of the unmapped reads against the whole genome, and 2.5 to 3 minutes with them dropped (AMD EPYC 7763 and Intel Xeon Platinum 8573C runners). MEGAnE's authors give about an hour for a 30x genome on 4 threads. MEGAnE runs with 16 GB of memory in both entry points: `scripts/50-mei.sh` passes `--memory 16g`, and `conf/base.config` gives the `MEI` process 16 GB (doubled on the one retry after an out-of-memory exit) instead of its `process_high` label's 32 GB. Not measured here on a whole genome. AnnotSV on the insertions gets 4 CPUs and 8 GB, as in step 05.

## Validation

The e2e case `tests/e2e/mei-1.sh` runs the step on HG002's chr20: the reads of chr20 in GIAB's 60x novoalign BAM, half the pairs (30x), realigned with minimap2 to chr20 alone, the unmapped reads dropped (so the case fits a runner's time), and the k-mer set of that chromosome. Its Alu calls are scored against the GIAB v5.0q insertions of HG002 on chr20 (RepeatMasker SINE/Alu, interspersed, 250 to 400 bp, inside the benchmark regions): recall must be at least 0.70 and precision at least 0.80, within 100 bp. On 2026-10-27 it found 12 of the 15 insertions (0.80), and 13 of its 15 Alu calls inside the benchmark regions matched a GIAB insertion (0.87); with the unmapped reads kept, the same 12 of 15 and 13 of 15. With every read pair within 1 kb of a GIAB Alu insertion removed, it found 0 of 15 and the recall check failed. The same case then runs AnnotSV on the calls with made-up AnnotSV data (the bundle of `tests/fixtures/annotsv/class45.sh` plus one made-up gene over all of chr20) and checks that every call comes back inside that gene, with its location there and class `NA`; the same VCF read without AnnotSV's table fails that check. It also runs the Nextflow pipeline with `--tools mei` and `--annotsv_annotations` on the same BAM and checks that its calls table is the bash step's, byte for byte. `tests/test_mei_report.py` plants an SVA in the FKTN window of a synthetic MEGAnE VCF and checks it is seen, and not seen without it, and joins a table in AnnotSV's columns back to the calls.

## Notes

- Not seen is not proof of absence. An SVA is about 2 kb of repeat, and short reads find its edges, not its body; a low-complexity or repeat-rich window can hide one. A finding is confirmed in a clinical laboratory, usually by long-range PCR across the site.
- Carrying one copy of the FKTN or BBS1 insertion means carrier status for a recessive condition, not the disease. The TAF1 insertion is on chrX: one copy in a male means the X he has carries it.
- The e2e case scores Alu calls only: chr20 has 15 GIAB Alu insertions that pass its filter, and too few L1 or SVA insertions to count. The L1 and SVA calls are listed but not measured here.
- MEGAnE's last release is a 2022 beta. Its image is amd64 only, so on an Apple Silicon Mac Docker runs it under emulation (not tested here).
- The step needs paired-end short reads. It has not been tested on long reads.
