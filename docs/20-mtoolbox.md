# Step 20: Mitochondrial Variant Calling and Heteroplasmy Detection

## What This Does
Calls mitochondrial DNA variants with heteroplasmy fractions — detecting variants present in only a fraction of your mitochondrial copies. Uses GATK Mutect2 in mitochondrial mode.

## Why
The main VCF (step 03) calls chrM as if it were a diploid nuclear contig. This step goes deeper, and step 12 then reads its calls for the haplogroup and the contamination check (run this step first):
- **Heteroplasmy detection**: Identifies variants present in only a fraction of mtDNA copies (clinically important for mitochondrial diseases)
- **Dedicated mitochondrial calling**: Mutect2's mitochondrial mode handles the unique properties of mtDNA (high copy number, circular genome, no recombination)
- **Somatic-grade sensitivity**: Detects variants at allele fractions as low as 1-3%

## Tool
- **GATK Mutect2** (Broad Institute) in `--mitochondria-mode`

> **Note:** This step was originally planned for MToolBox, but no working Docker image exists for MToolBox (see [lessons-learned.md](lessons-learned.md#mtoolbox-no-working-docker-image-exists)). GATK Mutect2 is the standard clinical alternative. The script is still called `scripts/20-mtoolbox.sh`; the name is historical.

## Docker Image
- `GATK_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
# Extract chrM reads
samtools view -b sorted.bam chrM > chrM.bam
samtools index chrM.bam

# Call variants in mitochondrial mode
gatk Mutect2 \
  -R reference.fasta \
  -I chrM.bam \
  -L chrM \
  --mitochondria-mode \
  --max-mnp-distance 0 \
  -O chrM_mutect2.vcf.gz

# Filter
gatk FilterMutectCalls \
  -R reference.fasta \
  -V chrM_mutect2.vcf.gz \
  --mitochondria-mode \
  -O chrM_mutect2_filtered.vcf.gz

# Mark possible NuMTs, given the median autosomal depth
gatk NuMTFilterTool \
  -R reference.fasta \
  -V chrM_mutect2_filtered.vcf.gz \
  --autosomal-coverage 30 \
  -O chrM_filtered.vcf.gz
```

The script runs these steps: `./scripts/20-mtoolbox.sh your_sample`. Then it annotates every call with ClinVar and gnomAD (below).

### ClinVar and gnomAD beside each call

The last stage (`bin/collect_summary.py mito-annotate`, in `PYTHON_IMAGE`) writes `${SAMPLE}_chrM_annotated.tsv`: one row per allele of the filtered VCF, with

- ClinVar's record of the same allele from `clinvar/clinvar_chr.vcf.gz`, the ClinVar file `setup.sh` installs: `clinvar_id`, `clinvar_significance` (CLNSIG), `clinvar_review` (CLNREVSTAT) and `clinvar_condition` (CLNDN);
- gnomAD v3.1's chrM sites of the same allele (`gnomad_mito/`, pinned in `versions.env`): `gnomad_af_hom` (the share of gnomAD's samples with the allele at 95% or more), `gnomad_af_het` (at 10% to 95%) and `gnomad_max_hl` (the highest level gnomAD saw).

Alleles are matched by position, REF and ALT after the bases they share are trimmed, so Mutect2's multi-allelic indels meet ClinVar's and gnomAD's single-allele records. A file that is not installed leaves its columns `.`, and the table's first lines say which. The reports list each PASS call that is heteroplasmic, or that ClinVar classes as anything but benign (a record with no classification, `not_provided`, does not count), with these columns, and count the ClinVar pathogenic or likely pathogenic ones.

Step 13's VEP annotates the main VCF, whose chrM records come from the nuclear caller: a heteroplasmic call only Mutect2 makes has no VEP record, so this step annotates its own calls. In the Nextflow pipeline `MITO_VARIANTS_ANNOTATE` does the same after `MITO_VARIANTS`, with the files of `--mito_clinvar` and `--gnomad_mito`; without them its table names each file as not given and leaves its columns `.`.

[MITOMAP](https://www.mitomap.org/) has the curated disease status of mtDNA variants. It stays a manual lookup: it answers scripted requests with a 403, so no step reads it, and none of its tables is shipped.

### NuMTs

NuMTs are copies of mitochondrial DNA in the nuclear genome. Their reads can map to chrM and look like low-level heteroplasmy. GATK's NuMTFilterTool marks an allele `possible_numt` when its depth is no more than such copies could give at the sample's autosomal depth. The script reads the median autosomal depth from step 16b's mosdepth output (`mosdepth/<sample>.mosdepth.summary.txt` and `.global.dist.txt`), so run step 16b first; `AUTOSOMAL_COVERAGE=30` sets it by hand. Without either, the filter runs at depth 0 and marks nothing, and the step says so. The Nextflow `MITO_VARIANTS` runs the same filter with the same depth, read from `MOSDEPTH`'s output when `mosdepth` is in `--tools` (a default tool); without it, at depth 0.

## Output
- `${SAMPLE}_chrM_mutect2.vcf.gz` — Raw mitochondrial variant calls
- `${SAMPLE}_chrM_mutect2_filtered.vcf.gz` — after FilterMutectCalls
- `${SAMPLE}_chrM_filtered.vcf.gz` — after FilterMutectCalls and NuMTFilterTool, with PASS or the reasons a call failed (`possible_numt` among them); the file the reports and step 12 (haplogroup, haplocheck) read
- Each variant includes an `AF` (allele fraction) field indicating heteroplasmy level
- `${SAMPLE}_chrM_annotated.tsv` — each allele of the filtered VCF with its FILTER, AF, ClinVar record and gnomAD v3.1 frequencies (above)

## Interpreting Heteroplasmy
| AF Level | Meaning |
|---|---|
| >0.95 | Homoplasmic — effectively fixed variant |
| 0.10-0.95 | Heteroplasmic — mixed population, clinically significant threshold varies |
| 0.03-0.10 | Low-level heteroplasmy — may be age-related somatic |
| <0.03 | Near detection limit |

Count only PASS calls: a call marked `possible_numt` or with another filter is not evidence of heteroplasmy. Calls near the two ends of the control region (about chrM:1-500 and 16,000-16,569) are less reliable: the reference is circular, so reads that span its end align poorly at both edges.

Length variants in the D-loop poly-C tracts (m.303-315 and m.16184-16193), one C more or less, are found in most people and often read as a mix of lengths. The reports list them as "length heteroplasmy in the D-loop poly-C tracts" and leave them out of the heteroplasmic count. A call counts as in a tract when the bases it changes (an indel: the bases after its anchor base) overlap one.

## Runtime
1 to 2 minutes per sample through the pipeline, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)).

## Notes
- `--mitochondria-mode` disables several filters inappropriate for mtDNA: no germline filtering, adjusted LOD thresholds, handles high copy number.
- `--max-mnp-distance 0` prevents merging nearby variants into multi-nucleotide polymorphisms.
- The GATK Docker image is large (~2.2 GB) but well-maintained and versioned.
- Disease annotation: ClinVar and gnomAD come with the step (`${SAMPLE}_chrM_annotated.tsv`); [MITOMAP](https://www.mitomap.org/) is a manual lookup in a browser, as it refuses scripted requests.
- Some mitochondrial diseases require heteroplasmy above a tissue-specific threshold (e.g., m.3243A>G MELAS requires >60% in blood).

## What is not done, and why

- **The shifted-reference pass for the control region.** Broad's mitochondria pipeline calls the region around the chrM start a second time on a reference rotated by 8 kb, so reads that span the circle's end align well. Here those calls stay less reliable (see above).
- **mtDNA-Server 2 (or mutserve) in place of this step.** Both were last released in December 2024 (mtDNA-Server 2.1.16, mutserve 2.0.3). The cheaper fix for the main source of false heteroplasmy, NuMT reads, is in place (NuMTFilterTool above), and haplocheck, the contamination check of the mtDNA-Server tools, runs in step 12 on these calls.

