# Step 52: MUC1 VNTR Frameshifts (VNtyper)

## What This Does

Looks for the frameshifts inside the MUC1 coding VNTR that cause autosomal dominant tubulointerstitial kidney disease (ADTKD-MUC1), with VNtyper. The step says **positive**, **possible**, **negative** or **not assessed**, and gives the share of known positives VNtyper found at this sample's depth over the repeat. Opt-in.

It is a research call, not a diagnosis. VNtyper's authors ask for any call to be confirmed by another method: SNaPshot for dupC, long-read sequencing for other variants.

## Why

ADTKD-MUC1 is a dominant kidney disease that leads to kidney failure in adulthood. It is caused by frameshifts inside a large coding VNTR in exon 2 of MUC1: a repeat of 60 bp units, 20 to over 125 of them per allele. The most common change is one extra C in a run of seven C (dupC). The causal variants were found by other means because massively parallel sequencing missed them (Kirby et al., Nat Genet 2013, PMID 23396133).

Short reads from inside the repeat fit every unit equally well and map with mapping quality 0, so DeepVariant and GATK do not call these frameshifts. ExpansionHunter counts repeat units and cannot see a one-base change inside one of them. VNtyper (Saei et al., iScience 2023, PMID 37456840) reads the region's k-mers with Kestrel instead of alignments and calls the frameshift directly. This step owns MUC1: a MUC1 value from the repeat-expansion steps is not a call.

## Tool

- VNtyper, the release `VNTYPER_IMAGE` pins (2.0.47 today; BSD-3-Clause, https://github.com/hassansaei/VNtyper), `vntyper pipeline` on the BAM with `--fast-mode` and `--reference-assembly hg38`
- `bin/vntyper_report.py` (standard library Python, run with the image's Python) for the reading

## Docker Image

- `VNTYPER_IMAGE`: the publisher's image, `ghcr.io/hassansaei/vntyper`, at the release `versions.env` pins. It carries VNtyper's own reference data (its MUC1 motifs and the indexes it realigns against), so the step installs nothing. The image is about 3 GB compressed. Docker Hub's `saei/vntyper` stops at 2.0.3, and Biocontainers has no build.

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. The line is marked optional, so `setup.sh` does not pull it: the first run of the step does.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps), aligned to GRCh38 with `chr` contig names. VNtyper checks the length of chr1 in the header against the assembly it is told and stops when they differ.

## Command

```bash
./scripts/52-vntyper.sh your_name
```

With `run-all.sh`: `TOOLS=...,vntyper`. It is opt-in, so a run lists step 52 only when `TOOLS` names it, and it is skipped for a sample with no BAM. In the Nextflow pipeline: `vntyper` in `--tools` (process `VNTYPER`, for each sample with a BAM). VNtyper reads its references through paths relative to `/opt/vntyper`, so both run it from there, and both replace the image's entrypoint, which accepts only its own command words.

## How the reading is made

VNtyper slices the reads of chr1:155,184,000-155,194,000, genotypes the VNTR with Kestrel and grades the call. Its documentation describes the grades. `bin/vntyper_report.py` reads its `kestrel/kestrel_result.tsv` and `coverage/coverage_summary.tsv`:

| Result | When |
|---|---|
| `positive` | Kestrel calls a frameshift with confidence `High_Precision` or `High_Precision*`, VNtyper did not flag it, and at most half of the VNTR has no read |
| `possible` | a call with confidence `Low_Precision`, a call VNtyper flagged, a call over a VNTR more than half without reads, or a call where VNtyper reported no share of the VNTR without reads |
| `negative` | no call, and VNtyper's coverage check is `PASS` or `REDUCED` (on GRCh38, a mean of 100x or more, or 50x to 100x) |
| `not assessed` | no call and a coverage check of `FAIL` or `NOT_EVALUATED` (under 50x, or no depth measured), VNtyper stopped with an error, or it wrote no result (a missing file, or a Kestrel file with no result row) |

### Sensitivity at the observed depth

The table's `sensitivity` row gives the mean depth over the VNTR and what VNtyper's authors measured at that depth. They downsampled confirmed positive GRCh38 exomes and re-ran them: 28 exomes at nine fractions with VNtyper 2.0.3 and 14 at nine fractions with 2.0.42 (VNtyper 2.0.47 documentation, "Coverage thresholds", read 2026-10-27):

| Mean depth over the VNTR | Runs with the variant found |
|---|---|
| 100x or more | 160 of 161 (99%) |
| 50x to 100x | 95 of 108 (88%) |
| under 50x | 48 of 109 (44%) |

These are repeated runs of 42 samples, not independent tests, and they are exomes, not genomes. The authors saw no false positive in 1,200 simulated negative runs. VNtyper's FAQ says a whole genome at 30x typically gives enough depth over the VNTR.

## Output

| File | Contents |
|---|---|
| `vntyper/${SAMPLE}_vntyper.tsv` | `metric`, `value`: `status`, `confidence`, `variant`, `flag`, `alt_depth`, `region_depth`, `depth_score`, `vntr_mean_depth`, `vntr_uncovered_pct`, `coverage_qc`, `sensitivity`, `message`, `advice` |
| `vntyper/${SAMPLE}_vntyper.log` | VNtyper's console output |
| `vntyper/run/` | VNtyper's own folder: `summary_report.html` (with an IGV view of the reads), `kestrel/`, `coverage/`, `pipeline.log`, `pipeline_summary.json` |

Both reports have a MUC1 VNTR section with the result, the variant and its confidence, VNtyper's coverage check, the sensitivity line and the advice. Without a positive, possible or negative reading, the reports list MUC1 VNTR frameshifts under what the pipeline did not assess.

## Runtime

A few minutes: VNtyper reads only the MUC1 region through the BAM index.

## Notes

- Research call, not a diagnosis. A negative lowers the chance of the frameshifts VNtyper looks for and does not rule ADTKD-MUC1 out.
- `--fast-mode` skips VNtyper's search of unmapped and partly mapped reads, as VNtyper's own examples do.
- VNtyper's optional adVNTR module (`--extra-modules advntr`) is not run.
- The e2e fixture has no reads over MUC1, so on it the step says not assessed. The container test and the e2e case `vntyper-1` run it on two BAMs from VNtyper's own public test data (Zenodo record 19181821): a sample with a frameshift, which must read positive, and one without, which must read negative.
