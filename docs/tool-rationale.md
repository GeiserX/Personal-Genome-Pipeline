# Tool Selection Rationale

This document explains why each default tool was chosen and what alternatives are available. For pipeline steps where only one reasonable tool exists (e.g., PharmCAT for pharmacogenomics, haplogrep3 for mitochondrial haplogroups), there is no entry here -- the choice is self-evident.

For hands-on benchmarking instructions, see [benchmarking.md](benchmarking.md).

---

## 1. Alignment: minimap2 vs BWA-MEM2

| Property | minimap2 | BWA-MEM2 |
|---|---|---|
| **Script** | `scripts/02-alignment.sh` | `scripts/02a-alignment-bwamem2.sh` |
| **Output directory** | `aligned/` | `aligned_bwamem2/` |
| **Docker image** | `quay.io/biocontainers/minimap2:2.28` | `quay.io/biocontainers/bwa-mem2:2.2.1` |
| **Speed (30X WGS)** | ~1-2 hours | ~4-8 hours |
| **Index build time** | ~30 minutes | ~1 hour |
| **Index size** | ~7 GB (.mmi) | ~24 GB (multiple files) |
| **XS tag** | Not produced | Produced |
| **Germline SNP/indel accuracy** | Equivalent to BWA-MEM2 | Equivalent to minimap2 |
| **Somatic calling accuracy** | Slightly lower (no XS tag) | Slightly higher |

### Default: minimap2

minimap2 is the default because:

1. **Speed.** It is 2-4x faster than BWA-MEM2 for paired-end short reads, which matters on consumer hardware where alignment is the bottleneck step.
2. **Equivalent germline accuracy.** For single-sample germline variant calling (this pipeline's focus), minimap2-aligned BAMs produce indistinguishable SNP and indel call sets compared to BWA-MEM2. The [UMCCR comparison](https://umccr.org/blog/bwa-mem-vs-minimap2/) found no significant accuracy difference for germline calls.
3. **Simpler index.** A single `.mmi` file (~7 GB) vs multiple BWA-MEM2 index files (~24 GB total).

### When to use BWA-MEM2 instead

- **Strelka2 compatibility.** Strelka2's SNP calling relies on the XS (suboptimal alignment score) SAM tag to distinguish true variants from alignment artifacts. minimap2 does not produce XS tags. With minimap2 alignments, Strelka2's SNP precision drops measurably while its indel calling is unaffected. If you plan to run Strelka2, align with BWA-MEM2.
- **Somatic calling.** Somatic variant callers (Mutect2, Strelka2 somatic mode) benefit from the XS tag for artifact filtering. If you are doing tumor-normal analysis (not the default pipeline path), use BWA-MEM2.
- **Reproducibility with clinical pipelines.** Most GATK Best Practices workflows and clinical lab pipelines use BWA-MEM2. If you need your BAM to be directly comparable to a clinical lab's output, align with BWA-MEM2.

### Key tradeoffs

The alignment step is the single longest step in the pipeline when starting from FASTQ. Choosing minimap2 saves 2-6 hours of wall-clock time with no measurable loss in germline variant calling accuracy. The only functional downside is incompatibility with Strelka2's SNP mode, which is not part of the default pipeline.

**References:**
- [UMCCR: BWA-MEM2 vs minimap2 for short reads (2021)](https://umccr.org/blog/bwa-mem-vs-minimap2/)
- [Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics (2018)](https://doi.org/10.1093/bioinformatics/bty191)

---

## 2. SNP/Indel Calling: DeepVariant vs GATK HaplotypeCaller vs FreeBayes

| Property | DeepVariant | GATK HaplotypeCaller | FreeBayes |
|---|---|---|---|
| **Script** | `scripts/03-deepvariant.sh` | `scripts/03a-gatk-haplotypecaller.sh` | `scripts/03b-freebayes.sh` |
| **Output directory** | `vcf/` | `vcf_gatk/` | `vcf_freebayes/` |
| **Docker image** | `google/deepvariant:1.6.0` | `broadinstitute/gatk:4.6.1.0` | `quay.io/biocontainers/freebayes:1.3.6` |
| **Algorithm** | Deep learning (CNN) | Local haplotype assembly + PairHMM | Bayesian haplotype-based |
| **SNP F1 (30X)** | ~0.999 | ~0.998 | ~0.994 |
| **Indel F1 (30X)** | ~0.994 | ~0.983 | ~0.960 |
| **Runtime (30X, 8 cores)** | 2-4 hours (CPU) | 4-8 hours | 8-20+ hours |
| **GPU support** | Yes (significant speedup) | No | No |
| **Multi-threading** | Yes (`--num_shards`) | Yes (`--native-pair-hmm-threads`) | No (single-threaded) |
| **GVCF output** | Yes | Yes (script uses normal VCF; switch to `-ERC GVCF` for cohort workflows) | No |
| **Region restriction** | No (`INTERVALS` not supported) | Yes (`INTERVALS` env var) | Yes (`INTERVALS` env var) |
| **Extra reference files** | FASTA + FAI | FASTA + FAI + .dict | FASTA + FAI |

### Strelka2: a fourth small variant caller

Strelka2 (`scripts/03c-strelka2-germline.sh`, output to `vcf_strelka2/`) is a fourth option for SNV and indel calling. Despite sometimes being grouped with SV tools because it ships alongside Manta, Strelka2's germline mode is a **small variant caller** — it calls SNVs and indels up to ~49 bp, not structural variants. See [Kim et al. 2018](https://doi.org/10.1038/s41592-018-0051-x).

Key characteristics:
- **Fast:** ~1-2 hours on 30X WGS with 8 threads (fastest of the four callers)
- **Good accuracy:** Comparable to GATK for SNPs and indels
- **BWA-MEM2 recommended:** Strelka2's scoring model relies on XS (suboptimal alignment score) tags. minimap2 does not produce these, causing reduced SNP precision. Use BWA-MEM2 alignments for best results.
- **Complements Manta:** In Illumina's intended workflow, Manta detects SVs and Strelka2 detects small variants — they are complementary, not alternatives to each other

### Default: DeepVariant

DeepVariant is the default because:

1. **Highest accuracy.** It leads all callers in precision, recall, and F1 for both SNPs and indels on 30X Illumina WGS, consistently winning the PrecisionFDA Truth Challenges and GIAB benchmarks.
2. **Low false positive rate.** Its deep learning model was trained on real sequencing data and recognizes systematic artifacts (strand bias, mapping artifacts) that rule-based callers miss.
3. **Reasonable runtime.** 2-4 hours on 8 CPU cores, 1-2 hours with GPU acceleration. Faster than GATK, far faster than FreeBayes.
4. **No extra reference files.** Only needs the FASTA and FAI -- no sequence dictionary (.dict) required.

### GATK HaplotypeCaller: when to use it

- **Clinical reproducibility.** GATK HaplotypeCaller is the gold standard in clinical genetics labs and CLIA-certified pipelines. If you need to compare your calls with a clinical lab report, running GATK alongside DeepVariant helps identify discrepancies.
- **Cohort calling.** GATK's GVCF mode (`-ERC GVCF`) produces per-sample genomic VCFs that can be joint-genotyped across a cohort with `GenomicsDBImport` + `GenotypeGVCFs`. This is the standard workflow for family studies and population-scale projects.
- **Balanced performance.** GATK sits between DeepVariant (highest precision) and FreeBayes (highest sensitivity), offering a different precision/recall tradeoff that can complement either tool.

### FreeBayes: when to use it

- **Maximum sensitivity.** FreeBayes calls more variants than either DeepVariant or GATK. Some of these extra calls are real variants in difficult regions that the other callers missed. Many are false positives.
- **Second-opinion caller.** Running FreeBayes as a complement to DeepVariant is useful for research exploration. Variants called by both tools are very likely real. Variants unique to FreeBayes deserve manual inspection.
- **No GATK dependency.** FreeBayes does not require a sequence dictionary (.dict) file or Java.

### FreeBayes: caveats

- **Single-threaded.** FreeBayes does not support multi-threading. Full-genome calling on 30X data takes 8-20+ hours. Use the `INTERVALS` variable to restrict to a chromosome or region for faster results.
- **Requires post-filtering.** Raw FreeBayes output includes many low-quality calls. Always filter to PASS variants or apply quality thresholds (`QUAL > 20`, `DP > 10`) before downstream analysis.
- **Highest false positive rate.** Among the three callers, FreeBayes has the lowest precision. Its extra sensitivity comes at the cost of more false positives, particularly for indels.

### Key tradeoffs

| Priority | Recommended Caller |
|---|---|
| Highest accuracy (single sample) | DeepVariant |
| Clinical lab compatibility | GATK HaplotypeCaller |
| Maximum sensitivity (research) | FreeBayes (with quality filtering) |
| Fastest runtime | DeepVariant (especially with GPU) |
| Cohort / family joint calling | GATK HaplotypeCaller (GVCF mode) |
| Consensus approach (2+ callers agree) | DeepVariant + GATK (or all three) |

**References:**
- [PLOS ONE variant caller comparison (2024)](https://doi.org/10.1371/journal.pone.0339891)
- [Poplin et al. A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology (2018)](https://doi.org/10.1038/nbt.4235)
- [Van der Auwera & O'Connor. Genomics in the Cloud: Using Docker, GATK, and WDL in Terra. O'Reilly (2020)](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/)

---

## 3. Structural Variant Calling: Manta vs Delly vs CNVnator

The pipeline runs up to three SV callers and merges their output (step 22). Each uses a different detection strategy.

| Property | Manta | Delly | CNVnator |
|---|---|---|---|
| **Script** | `scripts/04-manta.sh` | `scripts/19-delly.sh` | `scripts/18-cnvnator.sh` |
| **Output directory** | `manta/` | `delly/` | `cnvnator/` |
| **Docker image** | `quay.io/biocontainers/manta:1.6.0` | `quay.io/biocontainers/delly:1.7.3` | `quay.io/biocontainers/cnvnator` |
| **Signal types** | Paired-end + split-read | Paired-end + split-read + read-depth | Read-depth only |
| **Best for** | DEL, DUP, INV, small indels | INV, BND (translocations) | Large CNVs (>1 kb) |
| **Runtime (30X)** | ~20-60 min | ~2-4 hours | ~2-4 hours |
| **Typical call count** | 3,000-5,000 SVs | 5,000-15,000 SVs | 500-2,000 CNVs |

### Default: Manta (with optional Delly + CNVnator for consensus)

Manta is the primary SV caller because:

1. **Speed.** 20-60 minutes vs 2-4 hours for Delly or CNVnator.
2. **Balanced accuracy.** Good detection of deletions, duplications, inversions, and insertions with reasonable false positive rates.
3. **Indel bonus.** Manta's `candidateSmallIndels.vcf.gz` captures indels in the 20-50 bp range that DeepVariant sometimes misses and that are below the size threshold of other SV callers.
4. **Well-maintained.** Illumina's Manta is widely used in clinical SV pipelines.

### Delly: strengths and role

- **Best for inversions and balanced translocations.** Delly combines all three signal types (paired-end, split-read, read-depth) and is the most accurate caller for inversions (INV) and breakend events (BND).
- **Higher sensitivity overall.** Delly calls more SVs than Manta (5,000-15,000 vs 3,000-5,000), catching events that Manta's more conservative filters miss.
- **Higher false positive rate.** The extra sensitivity comes with more false calls. This is why the pipeline uses multi-caller consensus (step 22) rather than trusting any single caller.

### CNVnator: strengths and role

- **Large CNV specialist.** CNVnator uses read-depth signal only, making it the best tool for large copy number variants (>1 kb) including deletions, duplications, and aneuploidies.
- **Complementary signal.** Because it uses a completely different detection method (depth binning) than Manta or Delly (paired-end and split-read), CNVnator's calls provide independent confirmation.
- **Limited SV types.** Only detects deletions and duplications. Does not call inversions, translocations, or insertions.

### Additional SV caller: TIDDIT

| Tool | Strength | Limitation | Status |
|---|---|---|---|
| **TIDDIT** | Large inversions, translocations; low memory usage | Lower sensitivity for small SVs; needs BWA index for assembly mode | Available (v0.2.0) |

**Note:** Strelka2 was previously listed here but is a **small variant caller** (SNVs + indels up to ~49 bp), not a structural variant caller. It has been reclassified under section 2 (SNP/Indel Calling) as `scripts/03c-strelka2-germline.sh`. See [Kim et al. 2018](https://doi.org/10.1038/s41592-018-0051-x).

### Key tradeoffs

| Priority | Recommended Approach |
|---|---|
| Quick SV scan | Manta only (~20-60 min) |
| Balanced analysis | Manta + Delly with consensus merge (~3-5 hours) |
| Maximum sensitivity | All three callers + consensus merge (~4-8 hours) |
| Large CNV focus | CNVnator alone or CNVnator + Manta |
| Inversions / translocations | Delly (strongest for these SV types) |

The consensus merge (step 22) keeps only SVs called by two or more callers, reducing false positives at the cost of losing some real single-caller-only events. For more thorough SV analysis, consider also running AnnotSV (step 5) on the consensus set.

---

## 4. When to Run Multiple Callers

### Research and exploration

Run all available callers and compare. This gives you the broadest view of your variant landscape and lets you identify caller-specific artifacts vs real variants. Use the benchmarking workflow described in [benchmarking.md](benchmarking.md) to measure concordance.

**Recommended configuration:**
- Alignment: minimap2 (speed)
- SNP/indel: DeepVariant + GATK HaplotypeCaller + FreeBayes
- SV: Manta + Delly + CNVnator with consensus merge (step 22)

**Time estimate:** ~12-24 hours total (callers can run in parallel after alignment).

### High-confidence analysis

Use a consensus of two or more callers. Variants called by multiple independent tools have lower false positive rates. This is the approach used by many production WGS pipelines.

**Recommended configuration:**
- Alignment: BWA-MEM2 (clinical compatibility, XS tags)
- SNP/indel: DeepVariant + GATK HaplotypeCaller (keep intersection)
- SV: Manta + Delly with consensus merge

**Time estimate:** ~10-16 hours total.

### Quick personal analysis

Stick with the defaults. The default tools were chosen for the best single-tool accuracy, and running one caller per step is sufficient for personal exploration.

**Recommended configuration:**
- Alignment: minimap2 (speed)
- SNP/indel: DeepVariant
- SV: Manta

**Time estimate:** ~4-8 hours total.

### Decision matrix

| Question | Answer | Action |
|---|---|---|
| Do I need to match a clinical lab? | Yes | Use BWA-MEM2 + GATK HC |
| Do I want maximum accuracy for a specific variant? | Yes | Run DeepVariant + GATK, keep shared calls |
| Am I exploring broadly and can tolerate false positives? | Yes | Run all three SNP callers |
| Do I plan to use Strelka2? | Yes | Align with BWA-MEM2 (XS tag required) |
| Do I have a GPU? | Yes | Run DeepVariant with GPU for ~3x speedup |
| Am I analyzing a family/cohort? | Yes | Use GATK HC in GVCF mode for joint genotyping |
| Am I only looking at structural variants? | Yes | Manta + Delly + CNVnator with consensus merge |
| Do I want results as fast as possible? | Yes | minimap2 + DeepVariant + Manta (defaults) |
