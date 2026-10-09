# Common issues and FAQ

## Common Issues

| Problem | Cause | Fix |
|---|---|---|
| Container exits silently | Out of memory (OOM killed) | Increase Docker memory or reduce `--memory` flag. Check `docker logs <container>`. |
| "Permission denied" writing output | Container runs as non-root | Add `--user root` to `docker run` (already done in all scripts) |
| VEP cache download fails/times out | 26 GB download over unreliable connection | Use `wget -c` (supports resume). See [docs/13-vep-annotation.md](13-vep-annotation.md) |
| DeepVariant crashes on Mac | amd64 emulation + memory pressure | Reduce `--cpus` to 2 and `--memory` to 8g. Will be slow. |
| Wrong number of variants (too few) | Genome build mismatch | Ensure your BAM is aligned to GRCh38 (hg38), not hg19/GRCh37. Check with `samtools view -H your.bam \| grep SN:chr1` |
| 0-byte output files | Missing input or wrong path | Check that all input files exist. Run the script with `bash -x` for debug output. |
| "No such image" on `docker pull` | Image name/tag changed | Check the exact image name in the step's documentation. Biocontainer tags change frequently. |
| Very slow on macOS | Rosetta 2 emulation overhead | Expected. Consider running on a Linux machine or cloud instance for heavy steps. |

For detailed solutions, see:
- [docs/chip-data-guide.md](chip-data-guide.md) — using 23andMe/MyHeritage/AncestryDNA data with this pipeline
- [docs/troubleshooting.md](troubleshooting.md) — comprehensive troubleshooting guide organized by symptom
- [docs/lessons-learned.md](lessons-learned.md) — every failure encountered during development
- [docs/glossary.md](glossary.md) — alphabetical glossary of genomics terms

## FAQ

**Q: How much does WGS cost?**
$200-$1,000 depending on the vendor. Nebula/DNA Complete: $495 for 30X. Dante Labs: ~$300-600. Sequencing.com: $399. Novogene (research): ~$200-400. The pipeline itself is free.

**Q: I only have 23andMe/AncestryDNA data. Can I use this?**
Yes, partially. You can convert chip data to VCF and run pharmacogenomics (step 7), PRS (step 25), ClinVar screening (step 6), and ROH analysis (step 11). You cannot run alignment, variant calling, structural variants, repeat expansions, or ancestry analysis. See the **[chip data guide](chip-data-guide.md)** for conversion instructions, which steps work, and what to expect.

**Q: How long does the full pipeline take?**
On a 16-core/32GB desktop a default `run-all.sh` takes about 6-12 hours per sample, because the Nextflow pipeline it starts runs independent steps in parallel. Run it again after a failure or an interruption and it reuses every step that finished, so only the rest runs. [Hardware and storage requirements](hardware-requirements.md#runtime-per-step) has the time of each step, and the [pipeline overview](pipeline-overview.md#what-a-default-run-covers) lists which steps a default run includes and which are opt-in.

**Q: Do I need Java and Nextflow?**
For `run-all.sh`, yes: it starts the Nextflow pipeline, which needs Java 17 or later and Nextflow (26.04.7 is the release CI validates); without them it stops and prints the install line. Each step also runs as a script with Docker alone. [Full run](getting-started.md#full-run) has both.

**Q: Can I run this on a Raspberry Pi?**
No. Most bioinformatics Docker images are amd64 only, and a Pi doesn't have enough RAM. Minimum is a desktop/server with 16 GB RAM and an x86_64 CPU.

**Q: My data is aligned to hg19/GRCh37. What do I?**
Extract FASTQ from your BAM (`samtools fastq`) and re-align to GRCh38 using step 2. LiftOver is an alternative but introduces artifacts. Re-alignment is cleaner.

**Q: My BAM is GRCh38, but `validate-setup.sh` says it "was aligned to a different reference". Why?**
The pipeline aligns to the GRCh38 no-ALT analysis set (195 sequences). A BAM aligned to another GRCh38 file, such as one with ALT, HLA or decoy contigs (most vendor BAMs, and this pipeline's own BAMs from before the switch), has a different contig list, and reads at CYP2D6, the MHC and KIR were placed with mapping quality 0 where the reference holds two copies. Realign it: [Realigning after a reference change](realignment.md) has the commands and the checks.

**Q: I found a pathogenic variant. Should I be worried?**
Probably not. A typical genome shows 0-10 pathogenic/likely pathogenic ClinVar hits, almost all heterozygous (one copy) for recessive conditions. This means you're a **carrier**, not affected. Only worry if: (1) the variant is in a **dominant** gene, (2) you have **two** pathogenic variants in the same recessive gene, or (3) it is in a cancer predisposition gene (BRCA1/2, MLH1, etc.). See [interpreting-results.md](interpreting-results.md) for details.

**Q: What are VUS? Should I worry about them?**
VUS (Variants of Uncertain Significance) mean there is not enough evidence to classify the variant as pathogenic or benign. The majority will eventually be reclassified as benign. They are **not actionable** — do not change your medical care based on a VUS. Check back in 1-2 years with an updated ClinVar database.

**Q: How often should I re-run the analysis?**
ClinVar and other databases are updated monthly. Re-running the ClinVar screen (step 6, ~5 minutes) and CPSR (step 17, ~30 minutes) every 6-12 months with updated databases can catch newly classified variants. The compute-heavy steps (alignment, variant calling) do not need to be re-run unless you get new sequencing data or the pipeline changes its reference ([realignment](realignment.md)).

**Q: I ran the pipeline on two people (me and my partner). How do I compare?**
See [docs/multi-sample.md](multi-sample.md) for carrier cross-screening, pharmacogenomics comparison, and family analysis.

**Q: Is this clinically validated?**
No. This is a research/educational pipeline. It uses well-known open-source tools (DeepVariant, VEP, ClinVar, PharmCAT) but has not been through clinical validation. The tools themselves are research-grade and their results should not be treated as clinical diagnoses. Always discuss findings with a healthcare provider.

**Q: What about long-read sequencing (Nanopore, PacBio)?**
Supported since v0.3.0. See the [long-read guide](long-read-guide.md) for ONT and PacBio HiFi workflows using minimap2, Clair3, and Sniffles2. Most downstream VCF-based steps work as-is.

**Q: My vendor's download link expired. Can I still get my data?**
It depends on the vendor. Dante Labs deletes data after 30 days. Sequencing.com archives to cold storage (1-3 day retrieval). DNA Complete requires an active subscription. Novogene/BGI keep data for ~90 days. **Always download your raw data immediately.** See [docs/vendor-guide.md](vendor-guide.md) for vendor-specific details.

