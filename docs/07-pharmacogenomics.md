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
# match the reference is a 0/0 call and an uncovered one (./.) stays missing;
# then bin/pgx_refblock_fill.py calls 0/0 the positions longer than one base
# (indels) that reference blocks with GQ >= 20 cover (see below)
docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  -v "${GENOME_DIR}:/genome" \
  -v "$PWD/bin:/pgp-bin:ro" \
  "${PHARMCAT_IMAGE}" \
  bash -euo pipefail -c 'records=$(mktemp)
    bcftools view -R /pharmcat/pharmcat_regions.bed -Ov -o "$records" "$2"
    bcftools convert --gvcf2vcf -f "$1" -R /pharmcat/pharmcat_regions.bed -Ou "$2" \
      | bcftools view --trim-alt-alleles -i "GT!=\"mis\"" -Ov \
      | python3 /pgp-bin/pgx_refblock_fill.py /pharmcat/pharmcat_positions.vcf.bgz "$records" \
      | bcftools view -Oz -o "$3" --write-index=tbi' _ \
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

# Step 3 (when step 20 wrote its chrM calls): add the MT-RNR1 call to a copy
# of step 36's outside calls (see below); without step 36's file, leave out
# --outside; without the BAM, leave out the samtools line and --depth
docker run --rm \
  -v "${GENOME_DIR}/${SAMPLE}:/s" \
  -v "$PWD/bin:/pgp-bin:ro" \
  "${PHARMCAT_IMAGE}" \
  bash -c 'samtools depth -a -Q 20 -q 20 -r chrM:648-1601 /s/aligned/$1_sorted.bam > /tmp/mtrnr1.depth &&
    python3 /pgp-bin/pgx_mito_merge.py --jar /pharmcat/pharmcat.jar \
      --mito /s/mito/$1_chrM_filtered.vcf.gz \
      --outside /s/pgx_consensus/$1_outside_calls.tsv \
      --depth /tmp/mtrnr1.depth \
      --out /s/vcf/$1.pharmcat_outside_calls.tsv' _ ${SAMPLE}

# Step 4: run PharmCAT on the rewritten copy, with the outside calls of step 3
# (or step 36's file when step 20 did not run) when that file is not empty;
# without it, leave out the second -v and -po. -w /tmp gives the jar a
# writable working directory for its own log (see below).
docker run --rm \
  --cpus 2 --memory 4g -w /tmp \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  -v ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.pharmcat_outside_calls.tsv:/outside_calls.tsv:ro \
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
- **Positions longer than one base.** The expansion writes one record per base, each with a one-base REF. 51 PharmCAT 3.4.0 positions have a REF or ALT longer than one base: the insertions and deletions of DPYD, UGT1A1 (\*28, the TA-repeat insertion at chr2:233760233), TPMT, CYP3A4, CYP3A5, CFTR, CYP2C9, NUDT15, RYR1, CYP2D6 and G6PD, and two two-base substitutions in RYR1. PharmCAT's preprocessor does not take a one-base record as a call at such a position: it sets the record aside (FILTER `PCATxINDEL`, or `PCATxREF` for the substitutions) and lists the position in `missing_pgx_var.vcf`. Where PharmCAT's REF is longer than one base (the deletions, UGT1A1 \*28, the substitutions), its matcher then discards the record ("REF in VCF (C) does not match expected reference (CAT)"), so the position has no call and the alleles defined there are uncallable, even where the gVCF shows you match the reference over the whole span. At an insertion with a one-base REF the matcher still reads the set-aside record as a reference call ("PharmCAT cannot validate this position"). So step 0 runs `bin/pgx_refblock_fill.py` on the expanded file. For each of these positions in PharmCAT's positions file (`/pharmcat/pharmcat_positions.vcf.bgz` in the image), it writes a 0/0 record with PharmCAT's own ID, REF and ALT (INFO flag `PGX_REFBLOCK`) in place of the one-base record, when every base of the REF span lies in reference blocks with GQ 20 or more (GQ is Phred-scaled: 20 means a 1 in 100 chance the reference call is wrong; genotype 0/0, or 0 on a haploid chrX; FILTER `PASS` or `.`) and no other record overlaps the span: no variant call, no `./.` block, no block below GQ 20. Otherwise the expanded file stays as it was there: the position stays missing, and at an insertion the matcher reads the one-base record as before. The step log names the positions it filled (`pgx_refblock_fill: N of 51 PharmCAT positions longer than one base called 0/0 from reference blocks with GQ >= 20: ...`). The GQ threshold applies to these positions only: a one-base position takes any 0/0 block, as before. The Nextflow module does the same.
- **MT-RNR1 from the chrM calls.** PharmCAT 3.4.0 has no MT-RNR1 definition for its matcher: it reads the gene only as an outside call, one allele (chrM is not a pair), so from a VCF alone MT-RNR1 is No Result on every sample. When [step 20](20-mtoolbox.md) wrote its Mutect2 chrM calls (`mito/${SAMPLE}_chrM_filtered.vcf.gz`), the script runs `bin/pgx_mito_merge.py` (step 3 above). It reads PharmCAT's own MT-RNR1 alleles from the jar (`Reference` and names such as `m.1555A>G`) and adds one MT-RNR1 line to a copy of step 36's outside calls (`${SAMPLE}.pharmcat_outside_calls.tsv`, deleted afterwards): a named single-base allele with FILTER `PASS` and an allele fraction (AF) of 0.95 or more, the level the pipeline calls homoplasmic; the one PharmCAT ranks highest (increased risk, then uncertain, then normal) when there are several, with the others named in the step log; or `Reference` when no named allele is called. Mutect2 writes variant records only, so a position without reads has no record either: `Reference` needs read evidence, 20 reads or more (`samtools depth`, MAPQ and base quality 20) at every position where PharmCAT names an MT-RNR1 allele, from the sample's BAM (`aligned/${SAMPLE}_sorted.bam`, the one step 20 read). Without the BAM, or with a position below 20 reads, no `Reference` is passed and the log names the positions; a homoplasmic named allele needs no depth, its record is the evidence. Mutect2 writes the genotype 0/1 even at AF 1.0, so the AF decides, not the genotype. No line is added, and the gene stays No Result with the reason in the step log, when a named allele is heteroplasmic (AF below 0.95), filtered, or has no AF; when a homoplasmic `PASS` record at a named position has an ALT PharmCAT does not name; when a `PASS` insertion or deletion overlaps a named position or lies within 10 bases of m.960 or m.961, where PharmCAT names deletions in a C run that the script does not map; and when the file has no `PASS` chrM record. Without step 20's file the step says so and MT-RNR1 stays No Result: run step 20 before this step. The PHARMCAT module does the same with MITO_VARIANTS' calls.
- **PharmCAT's own log file.** The jar's built-in `logback.xml` writes `pharmcat.log` into the working directory. In the image that is `/pharmcat`, which the calling user cannot write, so every run printed a logback status dump with a `pharmcat.log (Permission denied)` stack trace (the calls were not affected). The script runs the jar with `-w /tmp`, a writable directory that goes away with the container. A Nextflow task runs in its own writable work directory, so the module needs nothing: the log lands there.
- **A backslash in a `##` header line.** PharmCAT up to 3.4.0 bundles vcf-parser 0.3.1, which stops with "Error parsing metadata: character to be escaped is missing" on one. The line is valid VCF; bcftools writes it for a soft filter with a quoted string (`bcftools filter -s LowDP -e 'FORMAT/DP<10 && GT!="0/0"'`). The script and the Nextflow module rewrite the header of PharmCAT's own copy (`${SAMPLE}.pharmcat_input.vcf`, deleted afterwards by the script): on `##` lines `\"` becomes `'` and any other `\` becomes `/`. PharmCAT's calls are the same with and without the rewrite. A newer PharmCAT is no fix yet: 3.4.0 still bundles vcf-parser 0.3.1 and fails the same way.

## Output
- HTML report with drug recommendations per gene, including HLA-A, HLA-B and CYP2D6 when [step 36](36-pgx-consensus.md) passed them as outside calls, and MT-RNR1 (aminoglycosides) when step 20's chrM calls gave it a call
- JSON report used by step 27 (`${SAMPLE}.report.json`)
- Preprocessed VCF (`${SAMPLE}.preprocessed.vcf.bgz`) generated as an intermediate
- `${SAMPLE}.missing_pgx_var.vcf`: the PGx positions absent from PharmCAT's input. With the gVCF these are the positions without coverage, plus the positions longer than one base whose span is not wholly inside reference blocks with GQ 20 or more; with a variants-only VCF they include every position where you match the reference.
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
- **MT-RNR1** is an outside call from step 20's chrM calls (see above): PharmCAT 3.4.0 does not call it from a VCF. A heteroplasmic named allele, such as m.1555A>G at AF 0.4, is not passed, so the gene stays No Result and the step log names the allele and its AF: read it there.
- **CYP4F2 can come out No Result** when unphased heterozygous calls elsewhere in the gene, a weak one among them, leave no diplotype that matches. CPIC's warfarin guidance needs only the \*3 position (rs2108622). Step 27 then reads the genotype PharmCAT itself lists at that position in `report.json` (the gene's `variants`) and prints it in a section headed "Single-Position Fallback (not a PharmCAT call)", with the CPIC warfarin rule when you carry the variant. It does the same for VKORC1 at rs9923231 (-1639G>A). PharmCAT's own result for the gene stays No Result in the phenotypes table. A position with no genotype is reported as not callable. The report carries no GQ or read depth for the position: check it in your VCF before acting on it.
- PharmCAT may disagree with lab reports on complex haplotypes (e.g., NAT2). Discrepancies can arise from different genome builds (hg19 vs hg38), different star allele definitions, or different variant calling pipelines. When a discrepancy matters clinically, compare both sets of raw variant calls and consult the PharmVar database for the current allele definitions — do not blindly trust either source.
- PharmCAT output structure changes across releases. If you upgrade PharmCAT, re-test step 27 (`27-cpic-lookup.sh`) because it parses the JSON output directly.

## Maintenance
- The pipeline is pinned to `PHARMCAT_IMAGE` for reproducibility. Upgraded from 2.15.5 in v0.3.0 (Apr 2026). Breaking changes in the 3.x series are documented in `docs/lessons-learned.md`.
- Treat **step 7 and step 27 as one upgrade unit**. If you bump PharmCAT, rerun both on a known sample and diff diplotypes, phenotypes, JSON structure, and CPIC recommendation text before merging.
- Recheck CPIC / ClinPGx guidance at least quarterly, or sooner if a drug-gene pair you expose in step 27 gets a meaningful update upstream.
