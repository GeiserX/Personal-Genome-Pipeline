# Starting from a Vendor VCF

<!-- Every bash block on this page after the Setup block runs in CI, in order, as pasted
     (tests/e2e/vendor-vcf-intake-5-docs.sh). Keep them runnable; use ```text for anything else. -->

Many providers send a VCF, the file that lists where your DNA differs from the reference genome, and nothing else. You can start the pipeline from that file. This page says what runs without reads, what the pipeline checks before it starts, how to fix the three things vendor files most often get wrong, and how to take your name out of the file before you share any result.

## What runs on a VCF alone

With only `sample`, `vcf` and `vcf_index` in the [samplesheet](nextflow.md#samplesheet-format), these tools run: `pharmcat`, `cpic`, `clinvar`, `vep`, `vcfanno`, `slivar`, `clinical_filter`, `cpsr`, `roh`, `prs`, `ancestry`, `mito_haplogroup` and `html_report`.

Everything that reads the aligned reads has nothing to work on, so it does not run: coverage (`mosdepth`, and `multiqc`, which reads it), structural variants and copy-number changes (`manta`, `delly`, `cnvpytor`, `duphold`, `annotsv`, `survivor_merge`), HLA typing, repeat expansions, telomere length, mitochondrial heteroplasmy and the read-depth CYP2D6 calls (`cyrius`, `pypgx`). If your provider can send the BAM or FASTQ, get it. The [Vendor Compatibility Guide](vendor-guide.md) says what each provider delivers.

## What the pipeline checks first

Before any analysis, `VCF_PRECHECK` reads each VCF once. It stops the run, naming the sample and the fix, when:

- **the file holds more than one sample.** One samplesheet row is one person, and a joint-called family VCF has a genotype column per person, so the steps would mix them. The message gives the sample count, the first five names and the command that keeps one person's column: `bcftools view -s <name> -a -c 1 -Oz -o <name>.vcf.gz <file>`, then `bcftools index -t <name>.vcf.gz`. Give each person their own row.
- **no contig is named the chr way.** The pipeline needs GRCh38 names: `chr1` to `chr22`, `chrX`, `chrY`, `chrM`. A file named `1`, `2`, `MT` (Ensembl style) would give an empty mitochondrial haplogroup and a wrong ROH summary, both without an error. Step 2 below renames them.
- **the file is a gVCF and `pharmcat` is selected.** PharmCAT refuses a gVCF. A gVCF also lists the stretches where you match the reference ("reference blocks": ALT `<*>`, `<NON_REF>` or `.` with an `END`). PharmCAT also refuses any file whose name contains `.g.vcf` or `.genomic.vcf`, even a plain one, so the check stops on the name too. Step 3 below removes the blocks.
- **no record has FILTER=PASS.** See [FILTER=PASS required](nextflow.md#filterpass-required).

What the check cannot catch:

- **The genome build.** It reads names, not coordinates. A GRCh37 (hg19) file named `chr1` passes and gives wrong results everywhere. Check the build before anything else: [Genome Build](vendor-guide.md#genome-build-grch37-hg19-vs-grch38-hg38).
- **Unplaced scaffolds.** Records on contigs outside the 25 main ones (`GL000195.1`, `KI270706.1`) keep their names after the rename. The ClinVar screen leaves out records on contigs that ClinVar or the reference lacks and says how many. Nothing else reads them.

`run-all.sh` starts this pipeline, so it runs the check too. The numbered scripts need the same chr-named, variants-only file and do not run it: fix the file first.

## Setup

Every command below runs in the folder that holds your VCF, writes its new files there, and runs bcftools from the image the pipeline pins. Set these four lines first:

```bash
# Setup: change these four lines to your own paths.
PGP=/path/to/Personal-Genome-Pipeline        # this repository
REF=/path/to/GRCh38_no_alt_analysis_set.fasta   # with its .fai beside it
VCF=/path/to/your_vendor_file.vcf.gz         # with its .tbi beside it
LABEL=sample1                                # a neutral name for the samplesheet and the outputs
```

```bash
cd "$(dirname "$VCF")"
source "$PGP/versions.env"
bcftools() { docker run --rm -i -u "$(id -u):$(id -g)" -v "$PWD:$PWD" -w "$PWD" "$BCFTOOLS_IMAGE" bcftools "$@"; }
IN=$(basename "$VCF")
```

Each step reads `IN` and points it at the file it wrote, so skip a step your file does not need and the next one still works.

## 1. Look at the file

```bash
bcftools index -s "$IN" | cut -f1 | paste -sd' ' -
bcftools query -f '%ALT\n' "$IN" | grep -cE '^(<\*>|<NON_REF>|\.)$' || true
bcftools view --no-version -h "$IN" | grep -vE '^(##(fileformat|FILTER|INFO|FORMAT|ALT|contig)=|#CHROM)' || true
```

The first line lists the contigs that hold records. `chr1 chr2 ...` is what the pipeline needs, and `1 2 ...` needs step 2. The second counts reference-only records, and anything above 0 means a gVCF that needs step 3. The third prints the header lines step 4 removes. Those often hold your name, your sample id, a file name or the provider's command lines.

## 2. Rename the contigs (Ensembl names only)

```bash
for c in $(seq 1 22) X Y; do echo "$c chr$c"; done > chr_map.txt
echo "MT chrM" >> chr_map.txt
bcftools annotate --rename-chrs chr_map.txt -Oz -o "$LABEL.chr.vcf.gz" "$IN"
bcftools index -t "$LABEL.chr.vcf.gz"
IN="$LABEL.chr.vcf.gz"
```

This is the command the pipeline prints when it stops on contig names. The map covers the 25 main contigs. Scaffolds keep their Ensembl names, which no GRCh38 reference here has. The ClinVar screen leaves them out, and nothing else reads them, so dropping them is fine too.

## 3. Remove the reference blocks (gVCF only)

```bash
bcftools view -e 'INFO/END!="."' -Ou "$IN" | bcftools view --trim-unseen-allele -Oz -o "$LABEL.variants.vcf.gz"
bcftools index -t "$LABEL.variants.vcf.gz"
IN="$LABEL.variants.vcf.gz"
```

The first `view` drops every record with an `END`, which in a gVCF are the reference blocks. `--trim-unseen-allele` removes the `<*>` or `<NON_REF>` allele that gVCF callers add to every variant record. The new name has no `.g.vcf` in it, so PharmCAT accepts it.

This is a workaround, and it costs PharmCAT calls. Without the blocks, PharmCAT cannot tell "you match the reference here" from "this position was not covered", so a variants-only VCF leaves about half of its genes Unknown. A gVCF is the better PharmCAT input. Step 7 run as a script expands the blocks of the gVCF that step 3 writes from FASTQ; the Nextflow pipeline does not expand a vendor gVCF yet.

## 4. Take your name out of the file

```bash
bcftools view --no-version -h "$IN" | grep -E '^(##(fileformat|FILTER|INFO|FORMAT|ALT|contig)=|#CHROM)' > header.txt
echo "$LABEL" > names.txt
bcftools reheader -h header.txt -s names.txt -o "$LABEL.vcf.gz" "$IN"
bcftools index -t "$LABEL.vcf.gz"
```

This keeps only the header lines the tools need and renames the sample column to `LABEL`. It lists what to keep rather than what to delete, because providers name their command lines differently (`##bcftoolsCommand`, `##commandline`, `##cmdline`, `##GATKCommandLine`). `bcftools annotate -x` cannot remove these lines. Run it before the pipeline, so no output repeats the old name. [Before you share outputs](nextflow.md#before-you-share-outputs) lists what else in the outputs can identify you.

## 5. Run the pipeline

```bash
printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$LABEL" "$PWD/$LABEL.vcf.gz" "$PWD/$LABEL.vcf.gz.tbi" > samplesheet.csv
nextflow run "$PGP/main.nf" -profile docker --input samplesheet.csv --reference "$REF" \
    --tools pharmcat,cpic,roh,mito_haplogroup,html_report --outdir results
```

`results/sample1/sample1_report.html` (with your `LABEL`) shows PharmCAT, the CPIC lookup, runs of homozygosity and the mitochondrial haplogroup. Add `clinvar` with `--clinvar` and `--clinvar_index` pointing at `clinvar_pathogenic_chr.vcf.gz` from [Step 0](00-reference-setup.md) for the ClinVar screen; the other tools and their databases are in [Nextflow Execution](nextflow.md).
