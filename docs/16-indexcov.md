# Step 16: Coverage QC and Sex Chromosome Verification

## What This Does
Ultra-fast whole-genome coverage profiling directly from the BAM index file. Infers sex chromosome copy number (CNchrX, CNchrY) and detects sex chromosome aneuploidies (XXY, XYY, X0). Produces per-chromosome depth uniformity plots.

## Why
Coverage QC catches alignment problems, sample swaps, and sequencing artifacts early — before spending hours on variant calling. Comparing the sex inferred from X/Y coverage with the sex you declare is a cheap sample-identity check: a swapped sample usually shows up as a mismatch. It can also reveal sex-chromosome aneuploidies like Klinefelter syndrome (XXY). It cannot tell two samples of the same sex apart; [step 33](33-sample-qc.md) checks the sex again from the reads, compares the samples of a run with each other, and estimates contamination.

## Tool
- **goleft indexcov** (Brent Pedersen)

## Docker Image
- `GOLEFT_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
./scripts/16-indexcov.sh your_name male     # or female; run-all.sh passes the sex you give it
./scripts/16-indexcov.sh your_name          # no declared sex: report only, no check
```

The script runs:
```bash
source versions.env   # from the repository root
docker run --rm \
  --cpus 1 --memory 1g \
  -v ${GENOME_DIR}:/genome \
  "${GOLEFT_IMAGE}" \
  goleft indexcov \
  --directory /genome/${SAMPLE}/indexcov \
  /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
```

### Sex check

goleft writes the inferred sex to `indexcov-indexcov.ped`, whose columns are `#family_id sample_id paternal_id maternal_id sex phenotype CNchrX CNchrY ...`. The `sex` column uses PED coding: `1` male, `2` female, anything else unknown. The `phenotype` column is always `-9`. The script finds the columns by header name, prints the raw row, CNchrX, CNchrY and the inferred sex, and then compares it with the sex you declared:

- **Match:** prints `Sex check: OK`.
- **Mismatch:** prints both values and exits non-zero. A mismatch means a sample swap, a wrong declared sex, or a sex-chromosome aneuploidy. Steps that take the declared sex (ExpansionHunter, step 9) would otherwise use a wrong value without warning.
- `SEX_CHECK=warn ./scripts/16-indexcov.sh your_name female` prints the mismatch and exits 0, for when you know why they differ (for example 47,XXY).

## Output Files
| File | Description |
|---|---|
| `indexcov-indexcov.ped` | PED file with CN values for chrX, chrY, and autosomes |
| `indexcov-indexcov.roc` | ROC-like data for each chromosome |
| `indexcov-indexcov.bed.gz` | Per-16KB-bin normalized depth across all chromosomes |
| `index.html` | Interactive HTML report with all plots |

## Interpretation
### Sex Chromosome Copy Number
| Karyotype | CNchrX | CNchrY | Meaning |
|---|---|---|---|
| 46,XY (male) | ~1.0 | ~1.0 | Normal male |
| 46,XX (female) | ~2.0 | ~0.0 | Normal female |
| 47,XXY (Klinefelter) | ~2.0 | ~1.0 | Male with extra X |
| 47,XYY | ~1.0 | ~2.0 | Male with extra Y |
| 45,X (Turner) | ~1.0 | ~0.0 | Female with single X |

### Coverage Uniformity
- Flat depth across a chromosome = good sequencing
- Dips or spikes = possible CNVs, GC bias, or capture artifacts
- Systematic deviations across all chromosomes = library prep or sequencing problems

## Runtime
~5 seconds per sample.

## Notes
- Should be run as an early QC step after alignment (step 2). It only reads the `.bai` index, not the full BAM.
- Requires the BAM index (`.bai`) to exist alongside the BAM file.
- Works on any number of samples simultaneously — useful for batch QC.
- The HTML report opens in any browser, but it loads Chart.js and jQuery from public CDNs; offline, the page opens without its plots.
- For single-sample runs, the sex chromosome plot is still useful but the population-level clustering view is less informative.
