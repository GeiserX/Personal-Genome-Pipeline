# Step 16: Coverage QC and Sex Chromosome Verification

## What This Does
Ultra-fast whole-genome coverage profiling directly from the BAM index file. Infers sex chromosome copy number (CNchrX, CNchrY) and detects sex chromosome aneuploidies (XXY, XYY, X0). Produces per-chromosome depth uniformity plots.

## Why
Coverage QC catches alignment problems, sample swaps, and sequencing artifacts early — before spending hours on variant calling. Sex chromosome verification confirms sample identity and can reveal clinically significant aneuploidies like Klinefelter syndrome (XXY).

## Tool
- **goleft indexcov** (Brent Pedersen)

## Docker Image
```
quay.io/biocontainers/goleft:0.2.4--h9ee0642_1
```

## Command
```bash
mkdir -p ${GENOME_DIR}/${SAMPLE}/indexcov

docker run --rm \
  --cpus 2 --memory 4g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/goleft:0.2.4--h9ee0642_1 \
  goleft indexcov \
  --directory /genome/${SAMPLE}/indexcov \
  /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
```

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
- The HTML report is self-contained and can be opened in any browser.
- For single-sample runs, the sex chromosome plot is still useful but the population-level clustering view is less informative.
