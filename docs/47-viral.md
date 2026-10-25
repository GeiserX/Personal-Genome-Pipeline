# Step 47: Viral DNA in the Unmapped Reads

## What This Does

Takes the reads that did not map to the human genome, realigns them to 24 viral genomes, and reports per virus how many reads it got, how much of its genome they cover and how deep, against the human depth of the same sample. Opt-in.

Its main finding is **inherited chromosomally integrated HHV-6** (iciHHV-6). In under 1% of people in the USA and the UK, and about 1% in other studies, a whole HHV-6A or HHV-6B genome sits in a telomere of one chromosome and is inherited like any other part of the genome ([Pellett et al. 2012](https://doi.org/10.1002/rmv.715); [Liu et al. 2020](https://doi.org/10.1371/journal.pgen.1008915)). Every cell then carries one viral copy, so the reads cover the whole viral genome at about half the depth of a human chromosome, which has two copies. Knowing it matters because a blood PCR test for HHV-6 then reads high, which can be taken for an active infection ([Pellett et al. 2012](https://doi.org/10.1002/rmv.715)).

The other viruses are reported as found or not found. **This is not a test for an infection.** Most adults carry herpesviruses that stay in their cells for life, and laboratory reagents carry viral DNA of their own: in 8,240 blood genomes, [Moustafa et al. 2017](https://doi.org/10.1371/journal.ppat.1006292) found sequences of 19 human viruses in 42% of the people, and judged most of the other 75 viruses they saw to be contamination of reagents and the environment.

## Why

The reads a human alignment leaves unmapped are already in the BAM and nothing else in the pipeline reads them. Moustafa et al. 2017 and Liu et al. 2020 both took exactly these reads from whole-genome sequencing of blood to find viral DNA and integrated HHV-6.

## Tool

- samtools to take the reads out of the BAM and to measure the depths
- minimap2 (`-ax sr`) to realign them to the viral genomes
- `bin/viral_report.py` (standard library Python) for the calls

## Docker Image

- `SAMTOOLS_IMAGE`, `MINIMAP2_IMAGE` and `PYTHON_IMAGE` in the bash step; `MINIMAP2_SAMTOOLS_IMAGE` (the same minimap2 and samtools releases in one image) and `PYTHON_IMAGE` in the Nextflow processes

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. The step adds no image.

## Reference data

`./scripts/setup.sh --viral-data <genome_dir>` downloads the 24 NCBI RefSeq records of `VIRAL_REF_ACCESSIONS` (versions.env) as one FASTA from NCBI's E-utilities and checks it against `VIRAL_REF_SHA256`; it is about 1.6 MB, stored as `reference/viral-<VIRAL_REF_VERSION>/viral_genomes.fa`. An accession with its version number names one fixed sequence. The download was checked on 2026-10-25: two fetches gave the same file, with 24 sequences.

| Group | Genomes (RefSeq) |
|---|---|
| HHV-6 | HHV-6A (NC_001664.4), HHV-6B (NC_000898.1) |
| Other herpesviruses | EBV (NC_007605.1), CMV (NC_006273.2), HSV-1 (NC_001806.2), HSV-2 (NC_001798.2), VZV (NC_001348.1), HHV-7 (NC_001716.2), HHV-8 (NC_009333.1) |
| Papillomaviruses | HPV-16 (NC_001526.4), HPV-18 (NC_001357.1), HPV-6b (NC_001355.1), HPV-11 (NC_001525.1) |
| Hepatitis and retroviruses | HBV (NC_003977.2), HCV (NC_004102.1), HIV-1 (NC_001802.1), HTLV-1 (NC_001436.1) |
| Polyomaviruses | BK (NC_001538.1), JC (NC_001699.1), Merkel cell (NC_010277.2) |
| Others | adenovirus 5 (AC_000008.1), Torque teno virus 1 (NC_002076.2), parvovirus B19 (NC_000883.2) |
| Lab control | phiX174 (NC_001422.1), which Illumina adds to sequencing runs |

HCV is an RNA virus with no DNA stage; its row stays because Moustafa et al. 2017 saw HCV DNA once, which they read as an integration.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps). A VCF alone has no reads, so a sample without a BAM is left out.
- the viral FASTA above (`VIRAL_REF` overrides its path, inside `GENOME_DIR`)

## Command

```bash
./scripts/47-viral.sh your_name
```

With `run-all.sh`: `TOOLS=...,viral`; a run lists step 47 only when `TOOLS` names it. In the Nextflow pipeline: `viral` in `--tools` with `--viral_ref reference/viral-1/viral_genomes.fa` (processes `VIRAL_READS` and `VIRAL_REPORT`).

## How it works

1. The reads with the unmapped flag (`samtools view -f 4`), which includes a read unmapped beside a mapped mate, and the reads mapped to chrEBV, the EBV genome the GRCh38 analysis set carries. `samtools fastq` leaves out secondary and supplementary records, so each read is realigned once.
2. minimap2 aligns them to the viral genomes as single reads, with no secondary alignments.
3. `samtools coverage` gives the reads per viral genome, and `samtools depth -aa` the depth at every base, zero included.
4. `samtools coverage` over the sample's BAM gives the human depth: chr1 to chr22, over the bases with at least one read. Both coverages count duplicates and leave out unmapped, secondary and QC-failed reads, so the two depths are measured the same way.

Per genome, `bin/viral_report.py` writes:

| Column | Meaning |
|---|---|
| `reads` | reads aligned to the genome |
| `breadth` | the fraction of the genome with at least one read |
| `mean_depth`, `mode_depth` | the mean depth, and the most common depth over the genome's bases |
| `relative_depth` | `mode_depth` / human depth |
| `copies_per_cell` | 2 x `mean_depth` / human depth (a cell holds two copies of each autosome) |
| `call` | see below |

The calls:

- **inherited integrated HHV-6 likely**: HHV-6A or HHV-6B with `relative_depth` of 0.30 or more. This is the rule of [Liu et al. 2020](https://doi.org/10.1371/journal.pgen.1008915): the mode of the per-base depth over the viral genome against the genome's sequencing depth, with integration at 0.30 or more. One integrated copy gives about 0.5; viraemia gives 10 to 1000 times less, by the studies they cite. The mode is used, not the mean, so a pile of reads on a short repeat (HHV-6's ends carry the same TTAGGG repeat as human telomeres) cannot lift a genome that has no reads elsewhere. When both species pass, the deeper one is called and the other row is read as reads the two genomes share.
- **lab control (phiX)**: any phiX read. It was added in the lab and says nothing about the sample.
- **viral DNA found**: at least 10 reads, over at least 10% of the genome.
- **a few reads**: at least 10 reads, over less than 10% of the genome.
- **not found**: fewer than 10 reads.

The 10 reads and 10% are this pipeline's choices, not published values.

## Output

| File | Contents |
|---|---|
| `viral/${SAMPLE}_viral.tsv` | one row per viral genome: `virus`, `accession`, `length`, `reads`, `breadth`, `mean_depth`, `mode_depth`, `human_depth`, `relative_depth`, `copies_per_cell`, `call`, `note` |
| `viral/${SAMPLE}_viral_reads.bam` (`.bai`) | the realigned reads, on the viral genomes; open it in IGV with the viral FASTA |
| `viral/${SAMPLE}_viral_coverage.tsv` | `samtools coverage` of that BAM |
| `viral/${SAMPLE}_human_coverage.tsv` | `samtools coverage` of the sample's BAM |

Both reports show the HHV-6 call and every genome with reads in their "Viral DNA in the unmapped reads" section.

## Runtime

Not measured on a whole genome. samtools reads the BAM twice from end to end, once for the unmapped reads and once for the human depth; minimap2 then aligns only the unmapped reads, to 24 small genomes.

## Notes

- Not a diagnosis. A virus found means its DNA was in the sequenced sample, which can come from latent virus in blood cells, a past or current infection, or a reagent.
- What integrated HHV-6 means for health is not settled. [Pellett et al. 2012](https://doi.org/10.1002/rmv.715) list what is known: it can make blood tests read as an active infection, transplant recipients with it may be at higher risk, and its role in disease needs more research. Each child of a carrier has about a one-in-two chance of inheriting it, as for any variant on one chromosome.
- A finding is confirmed by a test in a clinical laboratory; Liu et al. 2020 used droplet digital PCR.
- [Liu et al. 2020](https://doi.org/10.1371/journal.pgen.1008915) also describe people who carry only one direct repeat (DR) of HHV-6B, a "solo-DR": reads on the genome's ends only. Such a pattern gives a low mode depth here and is not called integrated; look at the BAM when HHV-6 has reads but no call.
- The human depth is measured where reads are, so a BAM cut to a few regions still gives the depth of those regions.
- On the e2e fixture the case plants synthetic reads drawn from the HHV-6B genome at half the human depth as unmapped reads, and the call is integration; at a tenth of the human depth over the same genome it is not. It plants CMV at a low level and phiX too. Nothing in the case is a finding about the fixture's own sample.
