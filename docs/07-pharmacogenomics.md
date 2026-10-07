# Step 7: Pharmacogenomics (PharmCAT)

## What This Does
Pharmacogenomic analysis — determines how you metabolize drugs based on your DNA using PharmCAT, a widely used research tool developed by CPIC/PharmGKB. Note: PharmCAT is a research tool, not a clinically validated diagnostic. Results should be confirmed by a certified laboratory before making prescribing decisions.

## Why
Identifies which drugs work well, which need dose adjustments, and which to avoid entirely. Covers 23 pharmacogenes affecting hundreds of medications.

## Tool
- **PharmCAT** v3.4.0 (Pharmacogenomics Clinical Annotation Tool, CPIC/PharmGKB)
- Upgraded from 2.15.5 to 3.2.0 in v0.3.0, and to 3.4.0 for the PharmVar and CPIC updates of 3.3.0 and 3.4.0. See `docs/lessons-learned.md` for migration notes (preprocessor rename, reporter flags, JSON property changes).

## Docker Image
- `PHARMCAT_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

# Step 0 (when step 3 wrote a gVCF): expand its reference blocks over
# PharmCAT's gene regions into a plain VCF, so a covered position where you
# match the reference is a 0/0 call and an uncovered one (./.) stays missing
docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  -v "${GENOME_DIR}:/genome" \
  "${PHARMCAT_IMAGE}" \
  sh -c 'bcftools convert --gvcf2vcf -f "$1" -R /pharmcat/pharmcat_regions.bed -Ou "$2" \
    | bcftools view --trim-alt-alleles -i "GT!=\"mis\"" -Oz -o "$3" --write-index=tbi' sh \
    "/genome/${REF_FASTA}" /data/${SAMPLE}.g.vcf.gz /data/${SAMPLE}.pgx_regions.vcf.gz

# Step 1: preprocess the VCF against the GRCh38 reference
# (-vcf /data/${SAMPLE}.vcf.gz when there is no gVCF)
docker run --rm \
  --cpus 2 --memory 4g \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  -v "${GENOME_DIR}:/genome" \
  "${PHARMCAT_IMAGE}" \
  python3 /pharmcat/pharmcat_vcf_preprocessor \
    -vcf /data/${SAMPLE}.pgx_regions.vcf.gz \
    -refFna "/genome/${REF_FASTA}" \
    -o /data/ \
    -bf ${SAMPLE}

# Step 2: rewrite backslashes in the header of PharmCAT's copy (see below)
docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  "${PHARMCAT_IMAGE}" \
  sh -c 'gzip -dc "$1" | awk "$3" > "$2"' sh \
    /data/${SAMPLE}.preprocessed.vcf.bgz /data/${SAMPLE}.pharmcat_input.vcf \
    '/^##/ { gsub(/\\"/, "\047"); gsub(/\\/, "/") } { print }'

# Step 3: run PharmCAT on the rewritten copy, with step 36's outside calls
# (HLA-A, HLA-B, an agreed CYP2D6) when that file is not empty; without it,
# leave out the second -v and -po
docker run --rm \
  --cpus 2 --memory 4g \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  -v ${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_outside_calls.tsv:/outside_calls.tsv:ro \
  "${PHARMCAT_IMAGE}" \
  java -jar /pharmcat/pharmcat.jar \
    -vcf /data/${SAMPLE}.pharmcat_input.vcf \
    -po /outside_calls.tsv \
    -o /data/ \
    -bf ${SAMPLE} \
    -reporterJson \
    -reporterHtml
```

### Input the preprocessor or PharmCAT refuses

- **gVCF.** PharmCAT refuses a gVCF, and decides by the file name too (`.g.vcf`, `.genomic.vcf`). Yet a gVCF is the better input: a variants-only VCF leaves about half of PharmCAT's genes Unknown, because PharmCAT cannot tell a reference call from a position that was not covered. So the script reads `vcf/${SAMPLE}.g.vcf.gz` when step 3 wrote one and expands its reference blocks into `${SAMPLE}.pgx_regions.vcf.gz` (step 0 above, deleted afterwards), which PharmCAT accepts; without a gVCF it reads `${SAMPLE}.vcf.gz`, and does not check whether that file is itself a gVCF (PharmCAT then stops on it). The Nextflow pipeline does the same with DeepVariant's gVCF, or with the gVCF a samplesheet row gives in its `gvcf` column (`run-all.sh` fills it from `vcf/${SAMPLE}.g.vcf.gz`); a gVCF given in the `vcf` column stops the run before any analysis when `pharmcat` is selected. For a vendor gVCF, remove the reference blocks and rename the file: [Starting from a Vendor VCF](vcf-first.md).
- **A backslash in a `##` header line.** PharmCAT up to 3.4.0 bundles vcf-parser 0.3.1, which stops with "Error parsing metadata: character to be escaped is missing" on one. The line is valid VCF; bcftools writes it for a soft filter with a quoted string (`bcftools filter -s LowDP -e 'FORMAT/DP<10 && GT!="0/0"'`). The script and the Nextflow module rewrite the header of PharmCAT's own copy (`${SAMPLE}.pharmcat_input.vcf`, deleted afterwards by the script): on `##` lines `\"` becomes `'` and any other `\` becomes `/`. PharmCAT's calls are the same with and without the rewrite. A newer PharmCAT is no fix yet: 3.4.0 still bundles vcf-parser 0.3.1 and fails the same way.

## Output
- HTML report with drug recommendations per gene, including HLA-A, HLA-B and CYP2D6 when [step 36](36-pgx-consensus.md) passed them as outside calls
- JSON report used by step 27 (`${SAMPLE}.report.json`)
- Preprocessed VCF (`${SAMPLE}.preprocessed.vcf.bgz`) generated as an intermediate
- `${SAMPLE}.missing_pgx_var.vcf`: the PGx positions absent from PharmCAT's input. With the gVCF these are the positions without coverage; with a variants-only VCF they include every position where you match the reference.
- Covers CYP2C19, CYP2D6, CYP2B6, CYP3A4/5, UGT1A1, DPYD, NAT2, TPMT, etc.
- Star allele calls with metabolizer status (Poor/Intermediate/Normal/Rapid/Ultra-rapid)

## Key Genes
| Gene | Drugs Affected | Example |
|---|---|---|
| CYP2C19 | SSRIs, PPIs, clopidogrel | One \*17 with one normal allele = rapid metabolizer → citalopram and escitalopram clear faster |
| CYP2D6 | 25% of all drugs, opioids, tamoxifen | Not called from a VCF; an outside call from step 36 when pypgx and Cyrius agree |
| UGT1A1 | Irinotecan, bilirubin clearance | Two \*28 alleles are associated with Gilbert's syndrome |
| DPYD | 5-FU, capecitabine (chemo) | Poor = lethal toxicity |
| NAT2 | Isoniazid, hydralazine | Slow acetylator = increased toxicity |

## Limitations
- **HLA-A, HLA-B and CYP2D6** are not called from a VCF: PharmCAT 3.4.0 reports them with no result (`callSource` `NONE`). [Step 36](36-pgx-consensus.md) gives PharmCAT T1K's HLA types (step 8) and a CYP2D6 call only when pypgx (step 32) and Cyrius (step 21, opt-in) agree on depth that passed its check; PharmCAT reads them with `-po` (`callSource` `OUTSIDE`). Run step 36 before this step, or run this step again after it. A disagreement, or one caller alone, leaves CYP2D6 without a result, on purpose.
- PharmCAT may disagree with lab reports on complex haplotypes (e.g., NAT2). Discrepancies can arise from different genome builds (hg19 vs hg38), different star allele definitions, or different variant calling pipelines. When a discrepancy matters clinically, compare both sets of raw variant calls and consult the PharmVar database for the current allele definitions — do not blindly trust either source.
- PharmCAT output structure changes across releases. If you upgrade PharmCAT, re-test step 27 (`27-cpic-lookup.sh`) because it parses the JSON output directly.

## Maintenance
- The pipeline is pinned to `PHARMCAT_IMAGE` for reproducibility. Upgraded from 2.15.5 in v0.3.0 (Apr 2026). Breaking changes in the 3.x series are documented in `docs/lessons-learned.md`.
- Treat **step 7 and step 27 as one upgrade unit**. If you bump PharmCAT, rerun both on a known sample and diff diplotypes, phenotypes, JSON structure, and CPIC recommendation text before merging.
- Recheck CPIC / ClinPGx guidance at least quarterly, or sooner if a drug-gene pair you expose in step 27 gets a meaningful update upstream.
