# Step 10: Telomere Length Estimation

## What This Does
Estimates relative telomere content from WGS BAM files by quantifying telomeric repeat reads (TTAGGG/CCCTAA).

## Why
Telomere length correlates with cellular aging at a population level. Comparing telomere content between individuals of similar age, sequenced on the same platform, provides a rough relative comparison. However, telomere length alone is not established as a clinically important standalone risk marker for individuals — it provides only a rough estimate of aging rate and is influenced by many non-age factors (genetics, cell type, technical variables).

## Tool
- **TelomereHunter** (German Cancer Research Center)

## Docker Image
- `TELOMEREHUNTER_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

> The image is the Bioconda build of TelomereHunter 1.1.0. It runs as you, not as root, and with `--plotNone`: its plots fail (its R has no `dplyr`, and its PyPDF2 is written for Python 3), and nothing in the pipeline reads them. It replaced a patched, digest-pinned build from a personal Docker Hub account.

### Which TelomereHunter, and what changed

Step 10 runs TelomereHunter 1 (1.1.0). On the image test's input, the HG002 fixture BAM plus 1,100 planted unmapped telomeric reads (600 `TTAGGG`, 400 `CCCTAA`, 100 with one `TCAGGG`), three images were run side by side:

| Image | intratel_reads | tel_content | TCAGGG per intratelomeric read |
|---|---|---|---|
| the earlier digest-pinned build | 1100 | 5594.234887 | 0.0909 |
| Bioconda TelomereHunter 1.1.0 (this step) | 1100 | 5594.234887 | 0.0909 |
| Bioconda TelomereHunter2 1.0.12 | 1100 | 5594.234887 | 0.0909 |

So a value from an earlier run of this step stays comparable. TelomereHunter2 (GPL-3.0, Python 3, maintained) gives the same numbers on this input, but it is not switched in yet: without `-b` it skips the band classes instead of using hg19 bands, which the Nextflow warning and parameter help describe, and its summary has other columns in another order. Moving to it is a separate change.

## Command
```bash
export GENOME_DIR=/path/to/your/data
./scripts/10-telomere-hunter.sh your_sample
```

`ALIGN_DIR` (default `aligned`) picks the BAM, for example `ALIGN_DIR=aligned_bwamem2`. `THREADS` (default 4) caps the container's CPUs; TelomereHunter has no thread option of its own.

TelomereHunter sorts telomeric reads into intratelomeric, subtelomeric and junction classes by chromosome band. Without `-b` it uses its own hg19 bands (`telomerehunter --help`: "If no banding file is specified, the banding information of hg19 will be used"), which put the band ends at hg19 positions on a GRCh38 BAM. The script passes UCSC's GRCh38 bands, which `setup.sh` installs as `reference/cytoBand.hg38.txt` (chr1-22, X and Y; see [reference setup](00-reference-setup.md#small-pinned-data-files)). When they are not installed it runs without `-b` and says so. What the script runs:

```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm --network none --user "$(id -u):$(id -g)" \
  --cpus 4 --memory 4g \
  -v ${GENOME_DIR}:/genome \
  "${TELOMEREHUNTER_IMAGE}" \
  telomerehunter \
    -ibt /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    -o /genome/${SAMPLE}/telomere/${SAMPLE} \
    -p ${SAMPLE} \
    --plotNone \
    -b /genome/reference/cytoBand.hg38.txt

# Output: telomere content report in ${GENOME_DIR}/${SAMPLE}/telomere/${SAMPLE}/
```

TelomereHunter uses `-ibt` (input BAM tumor) for a single-sample analysis. No `--tumor_only` flag is needed — when only `-ibt` is provided (without `-ibc` for a matched control BAM), TelomereHunter runs in single-sample mode automatically.

## Key Metric
- **`tel_content`** — intratelomeric reads per million reads with 48-52% GC, the GC content of telomeric repeats (TelomereHunter 1.1.0 divides the intratelomeric read count by `total_reads_with_tel_gc`). It is a relative telomere content, not a telomere length; the report labels it "Telomere content (relative)"
- Higher values indicate longer/more abundant telomeres
- Compare between samples of known age for relative ranking

## Important Notes
- The container runs as you (`--user "$(id -u):$(id -g)"`), so the outputs are yours. A telomere directory an older version wrote as root has to be taken back first: `sudo chown -R "$(id -u):$(id -g)" "${GENOME_DIR}/${SAMPLE}/telomere"`.
- Short-read WGS (150bp reads) systematically underestimates true telomere length because reads cannot span long repetitive regions
- Results are useful as a **relative comparison** between samples, NOT as an absolute telomere length measurement
- Long-read sequencing (PacBio/ONT) provides more accurate telomere length if absolute values are needed
