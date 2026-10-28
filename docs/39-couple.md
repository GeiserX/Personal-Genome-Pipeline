# Step 39: Couple Carrier Screen

## What This Does

For two people who both ran the pipeline and are marked as partners in a pedigree, the step lists the recessive and X-linked genes where they carry a change, and says what that means for their children:

- **Shared risk.** Both partners carry a change in the same autosomal recessive gene. A child who inherits both changes can have the condition: 1 in 4 per pregnancy when both are heterozygous, 1 in 2 when one partner is homozygous.
- **X-linked findings.** A change on the X chromosome. From a female partner, each son has a 1 in 2 chance of inheriting it and can have the condition, and each daughter a 1 in 2 chance of being a carrier. From a male partner, every daughter is a carrier and no son inherits it.
- **One partner only.** A change in a recessive gene that only one partner carries. A child can be a carrier, not affected, unless the other partner carries a change the screen does not see.

Each finding comes with its ClinVar review stars, the condition, and a one-line explanation. A residual-risk list names what the screen does not see. Opt-in: it runs only when the pedigree marks two samples as partners.

This is a screen, not a diagnosis, and not medical advice. A result that matters for a pregnancy needs a clinical carrier test and a genetic counsellor, who can also say what a negative screen leaves out for your ancestry.

## Why

Before this step a couple had to compare two ClinVar tables by hand ([Multi-Sample Comparison](multi-sample.md#carrier-screening-for-partners)). That missed the pairs where one partner carries a small variant and the other a deletion, an SMN1 copy loss or a repeat expansion in the same gene, and it said nothing about X-linked genes.

## Tool

`bin/couple_report.py` (standard library Python). It reads, for each partner:

| Finding | From | Counted when |
|---|---|---|
| ClinVar pathogenic or likely pathogenic hit | step 06, `clinvar/<sample>_clinvar_hits.vcf` (required) | the genotype carries the alternate allele and ClinVar's `GENEINFO` names a gene of the carrier list |
| Deletion | step 19 (Delly), `delly/<sample>_sv.vcf.gz` | `SVTYPE=DEL`, `FILTER` PASS and a per-sample `FT` that is not a failed filter, the genotype carries it, it is 50 bp or longer, it removes at least one exon of a listed gene (GENCODE exons), when it is longer than 1,000 bp read depth shows the loss, and it does not match a common deletion (see below) |
| SMN1 copy number | step 35 (Parascopy), `paralogs/<sample>_smn_copy_number.tsv` | SMN1 exon 7 at 1 copy (a carrier) or 0, with filter PASS and quality 20 or more, by step 35's rule |
| FMR1 and FXN repeats | step 09 (ExpansionHunter), `expansion_hunter/<sample>_eh.vcf` | an FMR1 allele of 55 CGG or more (a premutation; more than 200 is a full mutation), an FXN allele of 66 GAA or more |

Short-read Delly can mark PASS a long deletion, even one of many megabases, that read depth contradicts, and such a call removes an exon of hundreds of listed genes at once. So a deletion longer than 1,000 bp counts only when Delly's own read-depth copy number (`FORMAT/RDCN`, the depth inside the deletion against its flanks) is below a normal region's copies: 1 on a male's chrX or chrY outside the pseudoautosomal regions, 2 everywhere else. When the pedigree gives no sex, every region is judged against 2. Delly infers the sex itself and scales RDCN to it, so a long artefact on an unsexed male's chrX or chrY reads 1 and still counts: give the sex in the pedigree and in the samplesheet. A heterozygous loss reads 1, a homozygous or hemizygous loss 0. A long deletion with RDCN at the normal count or more, or with no estimate (`.` or a negative value, no reads in the flanks), is set aside, and the table of what was assessed gives the count. 1,000 bp is Delly's own line between small and large SVs in its read counts (`indelsize` in Delly's source); a deletion that size or smaller rests on its split reads and read pairs alone. A deletion VCF without the RDCN field (not written by Delly 2.x) is read as before, without the depth check.

Two more rules keep common and small events out of the deletion findings:

- **Length.** A structural variant is 50 bp or longer. Delly also calls shorter deletions. Those are left to the small-variant caller and count here only through ClinVar, so a Delly deletion under 50 bp is set aside, and one the small-variant caller misses is missed.
- **Common deletions.** Some deletions are common in the population, and both partners can carry one at the same place, often on both copies (the H1/H2 haplotypes of 17q21.31, with their copy-number differences, are an example). A shared row from such an allele would report a risk for children that is not there. So a deletion is set aside when it matches a deletion of 1% or more frequency in AnnotSV's benign loss set: `annotsv_annotations/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/benign_Loss_SV_GRCh38.sorted.bed`, into which `setup.sh --annotsv-data` merges gnomAD-SV's common deletions ([step 05](05-annotsv.md#gnomad-sv-frequencies-and-mane-select-transcripts)). 1% is AnnotSV's own `-benignAF` default. A match is a 50% reciprocal overlap: each of the two covers at least half of the other, so breakpoints that differ by some hundred bases still match, while a short deletion inside a long common one, or a long one around a short common one, does not. A record without a frequency is not taken as common. The cost: a rare pathogenic deletion at the same place as a common one is set aside too.

The table of what was assessed counts, for each partner, the deletions over a listed exon set aside by each rule, and says whether the common-deletion filter ran. Without the benign set (no AnnotSV annotations installed) the length and depth rules still apply, common deletions are not filtered, and the report says so for each partner whose Delly calls were read.

Only the ClinVar screen is required. Each of the others is used when its file is there, and the report lists every one that was not.

### The carrier list

`assets/carrier_genes.tsv` holds 3,292 genes, each with its GRCh38 span, its inheritance and a condition name:

- the 113 genes of the ACMG tier 3 carrier screening recommendation (Gregg et al., Genet Med 2021;23(10):1793-1806): 97 autosomal recessive genes and 16 X-linked ones, with the condition names of its tables;
- marked as ACOG too, the genes of the conditions in ACOG Committee Opinions 690 and 691 (2017): cystic fibrosis, spinal muscular atrophy, the hemoglobinopathies, fragile X and the Ashkenazi Jewish panel. All 20 are also in the ACMG list;
- every other gene of an Orphanet disorder with autosomal recessive or X-linked inheritance and a disease-causing germline association Orphanet assessed (Orphadata of 2026-06-23, CC BY 4.0). These rows cover the recessive genes of a ClinVar hit outside the ACMG list. A gene on chrX outside the pseudoautosomal regions is X-linked, every other gene autosomal recessive.

The `severity` column names the alleles of a gene that are not a full loss of function: a hypomorph (BTD p.Asp444His, the GALT Duarte variants) or a pseudodeficiency allele (HEXA c.739C>T and c.745C>T). Such an allele is never a carrier finding on its own; the couple report lists it under what the screen does not see, and the sample report's own carrier list sets it aside.

A ClinVar hit in a gene outside the list (a dominant gene, for example) is not part of this screen; each partner's own ClinVar section lists it. An emptied or damaged list stops the step instead of reporting no shared gene.

## Docker Image

- `PYTHON_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Command

The partners are a pedigree row each, with each other in the `partner` column ([`assets/pedigree.example.tsv`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/assets/pedigree.example.tsv), the format of [step 38](multi-sample.md#the-pedigree)):

```text
sample	sex	father	mother	partner
partner_a	female	.	.	partner_b
partner_b	male	.	.	partner_a
```

Run step 06 for both (and steps 09, 19 and 35 for what they add), then:

```bash
GENOME_DIR=/path/to/data ./scripts/39-couple.sh my_family pedigree.tsv
```

Step 38 does not have to run first. The results go to `$GENOME_DIR/my_family/`, the family's folder; without the second argument the step reads `$GENOME_DIR/my_family/pedigree.tsv`. `CARRIER_GENES` gives another carrier list, `CARRIER_EXONS` another exon BED (default the GENCODE exons `setup.sh` installs, `reference/gencode.v<release>.basic.exons.bed`; without one, deletions are not used and the report says so), `CARRIER_BENIGN_SV` another benign loss BED (default AnnotSV's, above, under `$GENOME_DIR/annotsv_annotations/`; without one, common deletions are not filtered and the report says so). `run-all.sh` runs one sample, so it does not run this step.

In the Nextflow pipeline the screen runs inside the family workflow, as the process `COUPLE_CARRIER`, when `--pedigree` marks two samples of the samplesheet as partners and `clinvar` is in `--tools` (without it the run warns and leaves the screen out). It reads `CLINVAR_SCREEN`, and `DELLY`, `EXPANSION_HUNTER` and `PARASCOPY` when they run. Deletions need `--gencode_exons <reference/gencode.v<release>.basic.exons.bed>`. The common-deletion filter reads the benign loss set under `--annotsv_annotations` when that directory holds it, whether or not `annotsv` is in `--tools`.

## Output

| File | Contents |
|---|---|
| `<family>/couple_<a>_<b>.tsv` (Nextflow: `<outdir>/family/`) | one row per gene where either partner has a finding: `gene`, `inheritance` (AR or XL), `lists` (ACMG, ACOG, Orphanet), `condition`, `result` (`shared`, `x_linked` or `one_partner`), `offspring_risk`, `partner_a`, `findings_a`, `partner_b`, `findings_b`, `explanation`. A and B are sorted by name |
| `<family>/couple_<a>_<b>.md` | the same in words, a table of what was assessed for each partner, and the residual-risk list |
| `<sample>/couple/<sample>_couple_carrier.tsv` | the same table, for each partner's own reports: the HTML report (step 24) and the text report get a Couple Carrier Screen section with the shared and X-linked genes |

A finding reads, for example, `chr22:42070000 A>G het, Pathogenic, 2 stars, ClinVar 12345` or `deletion chr12:47855001-47858000 (3.0 kb) het, 3 exons of the gene (Delly)`.

## What the screen does not see

The report ends with this list, plus the inputs that were missing for each partner:

- changes ClinVar does not list as pathogenic or likely pathogenic: a new loss-of-function change, or a variant of uncertain significance;
- genes outside the carrier list;
- SMN1 silent carriers: two copies on one chromosome and none on the other read as two copies;
- changes in genes with a near-identical copy (HBA1 and HBA2, CYP21A2, GBA1, SMN1), which short reads place poorly, and gene conversions; a deletion there counts only when Delly calls it;
- repeat expansions other than FMR1 and FXN, and the exact size of a long FMR1 or FXN allele;
- deletions Delly does not call, long deletions read depth does not confirm, duplications, inversions and insertions;
- Delly deletions shorter than 50 bp, which count only when ClinVar lists them;
- with AnnotSV's benign set, a rare pathogenic deletion that matches a common one;
- low-level mosaicism, and whether two changes of one person in one gene sit on the same copy.

A screen that finds no shared gene lowers the chance of a child with one of these conditions; it does not remove it.

## Runtime

Not measured on whole genomes. The step reads only the partners' result files, never a BAM or the reference.

## Rebuilding the carrier list

The ACMG and ACOG rows are written by hand and kept; the Orphanet rows and every span are made again from the sources:

```bash
python3 bin/couple_report.py genes --table assets/carrier_genes.tsv \
  --orphanet-genes en_product6.xml --orphanet-inheritance en_product9_ages.xml \
  --gtf gencode.v50.basic.annotation.gtf.gz --par-bed assets/par_grch38.bed --out new_carrier_genes.tsv
```

with `en_product6.xml` and `en_product9_ages.xml` from https://www.orphadata.com/data/xml/ and the GENCODE GTF of the release in `GENCODE_RELEASE` (`scripts/lib/common.sh`). It stops when a hand-written gene has no GENCODE gene line (give the current HGNC symbol). The header of the new file names the releases.

## Testing

`tests/test_couple_report.py` checks the rules on synthetic step outputs at real coordinates, including a red control with the list emptied. The e2e case `tests/e2e/couple-carrier-1.sh` runs steps 06 and 39 on HG002 and a synthetic partner, with planted ClinVar records and a planted deletion, and compares the Nextflow pipeline's couple table and report with step 39's. `tests/e2e/couple-common-sv-1.sh` runs step 39 with a small benign loss set: a 30 bp deletion and a common one are set aside, a rare one counts, and without the set the common one counts and the report says the filter was off; a stub run of the pipeline stages the set into `COUPLE_CARRIER` when `--annotsv_annotations` holds it, and none when it does not.
