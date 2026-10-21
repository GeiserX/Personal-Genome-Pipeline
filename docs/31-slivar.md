# Step 31: Variant Prioritization (slivar)

Prioritizes clinically interesting variants using tiered filters and detects compound heterozygote candidates. Optionally annotates results with gnomAD gene constraint metrics.

## Why

A typical VEP-annotated VCF contains thousands of MODERATE/HIGH impact variants. Most are common population variants or benign polymorphisms. This step filters down to the variants most likely to be clinically relevant using a three-tier approach, then flags potential compound heterozygotes that could cause autosomal recessive disease.

slivar (by Brent Pedersen, author of vcfanno, mosdepth, duphold) is a streaming VCF filter that replaces the unmaintained GEMINI database approach. It uses JS expressions for flexible filtering without loading variants into a database.

## Prerequisites

- VEP-annotated VCF from step 13, or vcfanno-enriched VCF from step 30 (recommended). The step 30 output is used only when it is newer than the VEP output; an older one is ignored with a notice.
- Optional: gnomAD v4.1 gene constraint TSV at `${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv`
- Optional: the sample's BAM, `aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` picks another folder), to phase the compound-het pairs. Without it every pair is `unknown`.

## Docker Images

- `SLIVAR_IMAGE` (compound het detection)
- `BCFTOOLS_IMAGE` (variant filtering via split-vep)
- `PYPGX_IMAGE` (its pysam reads the BAM for the phase of each compound-het pair, as in step 06's read check)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Usage

```bash
export GENOME_DIR=/path/to/your/data
./scripts/31-slivar.sh your_sample
```

## Filter Tiers

"Rare" is the rule step 23 uses (`bin/csq_rank.awk`): VEP's `MAX_AF` (highest frequency in any 1000 Genomes or gnomAD exome/genome population) below 1%; without `MAX_AF`, the values of `gnomADe_AF` and `gnomADg_AF` below 1%. Exome frequency alone would call a variant that is common in genomes but absent from exomes rare. A variant with no frequency value at all is not rare: VEP found no gnomAD allele that matches it, which happens to common indels written differently from gnomAD. HIGH and MODERATE variants of unknown frequency go into their own tier, `freq_unknown`. With none of these fields the tiers are not filtered by frequency and the step prints a notice.

The tiers take the impact of the most severe consequence over every transcript (`bcftools +split-vep -s worst`). Step 23 ranks on the MANE Select transcript instead; the TSVs here show that consequence beside the worst one.

### Tier 1: rare_high
- PASS variants with HIGH VEP impact (stop-gain, frameshift, splice donor/acceptor)
- Rare

### Tier 2: rare_moderate_deleterious
- PASS variants with MODERATE VEP impact (missense, in-frame indel)
- Rare
- At least one deleterious predictor hit (if vcfanno annotations available):
  - CADD PHRED >= 20 (SNV and/or indel tags, whichever are present)
  - REVEL >= 0.5 (deliberately below ClinGen's PP3_Supporting threshold of 0.644, so this tier casts a wider net; a hit here is not PP3 evidence)
  - AlphaMissense "likely_pathogenic"
  - SpliceAI annotation present (presence check only — bcftools cannot parse the pipe-delimited delta scores; threshold filtering at >= 0.2 is done in step 23)
- Without vcfanno: all rare MODERATE variants included (same as step 23)
- Only annotation tags that exist in the VCF header are referenced — partial installs (e.g., CADD SNVs only) work correctly

### Tier 3: clinvar_pathogenic
- PASS variants with ClinVar pathogenic or likely_pathogenic (from VEP CLIN_SIG field)
- Excludes `conflicting_interpretations_of_pathogenicity` (the substring match `~"pathogenic"` would otherwise include these)
- No frequency filter (pathogenic variants can be common carriers)

### Tier 2b: freq_unknown
- PASS variants with HIGH or MODERATE VEP impact and no frequency value
- Not counted as rare; kept so the compound-het search still sees them

All tiers are merged and deduplicated into a single prioritized VCF.

## Compound Heterozygote Detection

slivar's `compound-hets` command groups heterozygous variants by gene from the prioritized VCF and reports pairs that could form compound heterozygotes (two different damaging variants in the same gene, one from each parent). It requires a PED file (`--ped`) describing sample relationships. For singleton samples (no trio), the `--allow-non-trios` flag is required.

The command outputs VCF to stdout with `INFO/slivar_comphet` annotations linking partner variants. Each VCF record represents a unique variant; the `slivar_comphet` field lists all its compound-het partners (format: `sample/GENE/PAIR_ID/chrom/pos/ref/alt`, comma-separated). A gene with N variants produces up to C(N,2) pairs but only N VCF records. The script counts unique pair IDs from this field and exports a human-readable TSV with columns: GENE, CHROM, POS, REF, ALT, IMPACT, Consequence, GT, PHASE, MAX_AF, gnomADe_AF, gnomADg_AF, FREQ_TIER, TRANSCRIPT, MANE_Consequence -- sorted by gene so that compound-het partners appear in consecutive rows. The frequency columns are VEP's (`.` for a field VEP did not write), `FREQ_TIER` is `rare`, `unknown`, `common` (a ClinVar record) or `none` (no frequency field), `TRANSCRIPT` is the set step 23 ranks on (`MANE`, `CANONICAL` or `ANY`), and `MANE_Consequence` is the consequence on the MANE Select transcript (`.` when the variant has none). The summary TSV ends with the same six columns, before the constraint columns. The Nextflow `SLIVAR_PRIORITIZE` module applies the same tiers and writes the same summary columns; it writes no per-variant compound-het TSV.

### Phase of each pair, from the reads

Two variants in one gene hit both copies of the gene only when they sit on different copies (in trans). On the same copy (in cis) the other copy is intact, which for a recessive gene usually means no disease. The genotypes of one sample do not tell the two apart, but the reads can: a read, or a read and its mate, that covers both positions shows both alt alleles together in cis, and one without the other in trans.

`bin/comphet_phase.py` does this for every pair slivar lists, with pysam on the sample's BAM. It counts the reads and read pairs that show an allele at both positions. A read pair counts once, with what its two mates show together; a pair whose mates disagree at a position is not counted. Reads that count are mapped, primary, not duplicates, not QC-failed, with mapping quality 20 or more. The verdict is `trans` or `cis` when at least 2 reads or read pairs cover both positions and at least 80% of them agree; otherwise `unknown`, with the reason. SNVs, MNVs and indels that share their first base with the reference are phased; a complex allele is not.

Two positions further apart than a read pair spans (about 500 bp for short reads) are never covered by one fragment, so the pair stays `unknown`. For those, a parent's genome (step 38, [family analysis](multi-sample.md#family-analysis-step-38-opt-in)) or long reads can phase the pair; this step does not run WhatsHap over whole genes.

The step writes `${SAMPLE}_compound_het_pairs.tsv`, one row per pair: GENE, CHROM, POS_1, REF_1, ALT_1, POS_2, REF_2, ALT_2, DISTANCE, PHASE, CIS_READS, TRANS_READS, PHASE_NOTE. CIS_READS and TRANS_READS are the reads or read pairs that show the two variants on the same copy and on different copies. The per-variant TSV gets a PHASE column after GT: `trans` when the variant is in trans with at least one partner, else `unknown` when one of its pairs is unknown, else `cis`. The HTML and text reports count the pairs by phase and list the first 20 trans ones by name.

In a Nextflow run the `COMPHET_PHASE` process does the same (same script, `PYPGX_IMAGE`) after `SLIVAR`, on its compound-het VCF and the sample's BAM; a sample without a BAM gets every pair `unknown`. It publishes `slivar/<id>_compound_het_pairs.tsv`. The Nextflow report still shows only the prioritized count for this step: `HTML_REPORT` has a place for the pairs table, but the annotation workflow does not pass it on yet.

**Important:** `cis` and `unknown` pairs are not two-hit findings. A `trans` pair is a candidate worth a look in a gene that fits your health history, and confirming it still needs a clinical test.

## Gene Constraint Enrichment

If `gnomad_v4.1_constraint.tsv` is available, the summary TSV is enriched with per-gene constraint metrics:

| Column | Description | Threshold |
|---|---|---|
| LOEUF | Loss-of-function observed/expected upper bound | < 0.35 = constrained |
| pLI | Probability of LoF intolerance | > 0.9 = constrained |
| mis_z | Missense Z-score (gnomAD v4.1 `mis.z_score`) | > 3.09 = constrained |
| CONSTRAINED | YES if LOEUF < 0.35 or pLI > 0.9; NO if LOEUF or pLI is known and neither says so; `.` if the gene has no constraint values | Flag column |

The values come from `bin/constraint_join.awk`, the loader step 23 and the Nextflow module run too: only canonical transcripts count, and of a gene's two canonical rows (Ensembl and RefSeq) the Ensembl one wins, so steps 23 and 31 report the same value for a gene. When rows carry gene symbols and not one matches the table, the step fails rather than write `.` in every row.

Variants in constrained genes are more likely to be pathogenic -- these genes are under strong purifying selection against damaging variants.

## Output

| File | Description |
|---|---|
| `slivar/${SAMPLE}_prioritized.vcf.gz` | All prioritized variants (merged, deduplicated) |
| `slivar/${SAMPLE}_slivar_summary.tsv` | Human-readable table with gene constraint |
| `slivar/${SAMPLE}_compound_hets.vcf.gz` | Compound het candidate variants (VCF) |
| `slivar/${SAMPLE}_compound_hets.tsv` | Compound het candidates (human-readable TSV, with PHASE) |
| `slivar/${SAMPLE}_compound_het_pairs.tsv` | Each compound-het pair: cis, trans or unknown, with the read counts |
| `slivar/${SAMPLE}_rare_high.vcf.gz` | Tier 1: HIGH impact |
| `slivar/${SAMPLE}_rare_moderate_del.vcf.gz` | Tier 2: MODERATE + deleterious |
| `slivar/${SAMPLE}_freq_unknown.vcf.gz` | Tier 2b: HIGH or MODERATE, frequency unknown (when the VEP output has a frequency field) |
| `slivar/${SAMPLE}_clinvar_path.vcf.gz` | Tier 3: ClinVar P/LP |

## Runtime

4 to 10 minutes, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)); almost all of it is the prioritising pass, where bcftools split-vep filters.

## Interpretation

A typical 30X WGS sample produces:
- **rare_high:** 50-200 variants (loss-of-function in rare alleles)
- **rare_moderate_deleterious:** 200-1,000 variants (depends on predictor availability)
- **clinvar_pathogenic:** 0-10 variants (most are heterozygous carriers)
- **compound het candidates:** 1,000-2,000 pairs across 100-200 genes (combinatorial: N variants in a gene = N*(N-1)/2 pairs). Most variants of a gene lie further apart than a read pair spans, so most pairs read `unknown`

Focus review on:
1. Variants in **constrained genes** (CONSTRAINED=YES)
2. **Homozygous** rare HIGH/MODERATE variants (potential recessive disease)
3. **Compound het pairs in trans** in known disease genes
4. Any **ClinVar pathogenic** variants, especially in dominant disease genes

See [docs/interpreting-results.md](interpreting-results.md) for pathogenicity score thresholds.
