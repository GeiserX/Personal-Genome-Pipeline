# Step 42: Traits Panel

## What This Does

Reads your genotype at a small set of well-replicated variants and prints what each genotype means in plain words: lactase persistence (digesting milk as an adult), caffeine metabolism and coffee intake, alcohol flush, bitter taste, earwax type, the photic sneeze reflex, smelling asparagus in urine, muscle fibre type, eye, hair and skin pigmentation, red hair, hair thickness, folate metabolism, secretor status, vitamin D level, bilirubin (the Gilbert syndrome genotype), nicotine breakdown, two variants linked to a high Lp(a) level, and the APOE ε2/ε3/ε4 alleles. A rule table reads several sites together: the APOE pair and the TAS2R38 bitter-taste haplotypes (PAV, AVI), phased from the reads when that is needed. Each row names the primary paper it comes from. On by default.

## Why

These are the results most people find easiest to read, and the pipeline had no such section. SNPedia, the usual source, is CC BY-NC-SA, which a GPL-3.0 project cannot copy, so the table here is written from the primary papers instead: each row's reading follows the cited abstract or its open full text, and the alleles follow dbSNP's GRCh38 record. No text comes from SNPedia or any other database.

## Tool

- **bcftools** (`bcftools convert --gvcf2vcf`, `view`, `query`) reads the genotypes at the table's positions from the gVCF.
- **WhatsHap** (`whatshap phase`) phases the heterozygous sites of a rule that needs phase from the BAM's reads.
- **`bin/traits_report.py`** (Python standard library) reads each genotype against the table, and each rule against the rule table.

## Docker Image

- `BCFTOOLS_IMAGE`, `PYTHON_IMAGE` and `WHATSHAP_IMAGE`, pinned in `versions.env`; [Image versions](versions.md) lists the tags.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz` with its `.tbi`: the gVCF of [step 3](03-variant-calling.md). A gVCF holds reference blocks, so a position where you match the reference reads as ref/ref, and a position with no reads has no record and reads as "not covered". A variant-only VCF cannot tell those two apart, so the step does not read one. The gVCF must hold one sample: a gVCF with several sample columns stops the step with the count.
- `assets/traits.tsv`: the table, and `assets/genotype_rules.tsv`: the rules.
- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other steps), optional: WhatsHap reads it when a rule needs phase. Without it, or when WhatsHap fails (step 42 prints a warning and goes on, as `TRAITS_PHASE` does in Nextflow), such a rule is read only when at most one of its sites is heterozygous.

In the Nextflow pipeline (`traits` in `--tools`, part of the default list) the processes `TRAITS_SITES` and `TRAITS` run for every sample with a gVCF, and `TRAITS_PHASE` (WhatsHap) for every such sample with a BAM: one DeepVariant wrote (FASTQ, BAM or CRAM rows) or one the samplesheet gives in its `gvcf` and `gvcf_index` columns. A sample without one is skipped with a warning. From an existing VCF, `run-all.sh` passes a gVCF on only for a sample with a BAM and a VCF whose gVCF lies beside the VCF; any other sample that starts from a VCF gets step 42 listed as skipped, with the reason. A failed traits task is skipped and the run goes on, as for the other report-only tools.

## Command

```bash
./scripts/42-traits.sh your_name
```

## Output

| File | Contents |
|---|---|
| `traits/${SAMPLE}_traits.tsv` | one row per table row, then the APOE line: `rsid`, `gene`, `chrom`, `pos`, `trait`, `genotype` (reference allele first), `status`, `reading`, `ancestry`, `note`, `pmids` |
| `traits/${SAMPLE}_genotype_rules.tsv` | one line per rule: `rule`, `gene`, `trait`, `rsids`, `chrom`, `pos`, `genotypes` (one per site, reference allele first), `gq` (the GQ of each site's record, a reference block's own GQ for a homozygous reference site, `.` when the gVCF header has no FORMAT/GQ), `phase`, `result`, `status`, `reading`, `ancestry`, `note`, `pmids` |
| `traits/${SAMPLE}_traits_genotypes.tsv` | what bcftools read at the table's positions and the rule sites: chromosome, position, REF, ALT, the genotype as bases and its GQ |

`status` is `called`, `not_covered` (no genotype at the position, from no reads or a no-call: not read) or `other_allele` (an allele the table does not describe, or an indel starting at the position: not read). A record with a missing allele (a no-call or a half call such as `./T`) is dropped by the genotype filter, so its position reads as `not_covered`. A gVCF whose reference base differs from the table's GRCh38 base stops the step: the gVCF is on another build.

Both reports have a Traits section: every row with its genotype, reading, note and PubMed id, and a mark on the rows whose reading depends on ancestry. The APOE lines come last, after a warning, and in the HTML report they are folded away until you open them. The reports do not read the rules table yet: the TAS2R38 haplotype line is in `traits/${SAMPLE}_genotype_rules.tsv` and in the step's own output.

## The Table

`assets/traits.tsv` has one row per variant: GRCh38 chromosome and position, rsID, gene, the reference base (the GRCh38 base) and the allele the readings are about, the trait, one reading per genotype (ref/ref, ref/alt, alt/alt), whether ancestry changes the reading, a note, and the PubMed ids of the papers. Lines starting with `#` are comments.

| Trait | Variants | Papers (PubMed) |
|---|---|---|
| Lactase persistence | rs4988235 (-13910), rs182549 (-22018), rs145946881 (-14010), rs41380347 (-13915), rs41525747 (-13907) | 11788828, 17159977 |
| Alcohol flush, alcohol metabolism | rs671 (ALDH2*2), rs1229984 (ADH1B*2) | 2562960, 2014795 |
| Caffeine and coffee | rs762551 (CYP1A2), rs2472297 (CYP1A1-CYP1A2), rs4410790 (AHR) | 10233211, 21357676, 21490707 |
| Bitter taste (PTC, PROP) | rs713598, rs1726866, rs10246939 (TAS2R38) | 12595690, 15723792 |
| Earwax type | rs17822931 (ABCC11) | 16444273 |
| Photic sneeze reflex, asparagus smell | rs10427255, rs4481887 | 20585627, 27965198 |
| Muscle fibre type | rs1815739 (ACTN3 R577X) | 10192379, 12879365 |
| Eye colour | rs12913832 (HERC2), rs1800407 (OCA2) | 18252222, 20457092 |
| Pigmentation, freckles, red hair | rs16891982 (SLC45A2), rs1426654 (SLC24A5), rs12203592 (IRF4), rs1805007, rs1805008, rs1805009 (MC1R) | 15714523, 16357253, 18483556, 24267888, 9302268, 7581459 |
| Hair thickness | rs3827760 (EDAR) | 18065779 |
| Folate metabolism | rs1801133 (MTHFR C677T) | 7647779 |
| Secretor status | rs601338 (FUT2) | 7876235, 12692541 |
| Vitamin D level | rs2282679 in GC, the vitamin D binding protein gene | 20541252 |
| Bilirubin, Gilbert syndrome genotype | rs887829 (UGT1A1) | 19414484, 37456363 |
| Nicotine breakdown | rs1801272 (CYP2A6*2) | 2322567, 9827545 |
| Lp(a) level | rs10455872, rs3798220 (LPA); [step 46](46-lpa.md) adds the KIV-2 repeat count, and the report's Lp(a) section reads the three together | 20032323, 36036785, 29973585 |
| APOE ε2/ε3/ε4 | rs429358, rs7412 | 8346443, 9343467 |

Variants whose effect allele the papers we could read did not state were left out, among them SLC24A4 rs12896399 and TYR rs1393350 of the eye-colour panels.

### Adding a row

Read the primary paper first and write the reading from it, in your own words. Take the GRCh38 position, the reference base and the allele from dbSNP's record of the rsID, and the direction of the effect from the paper, checking the strand: an allele named on the gene's strand of a gene on the minus strand is the complement on GRCh38. Then run `python3 tests/test_traits_table.py --online`, which looks up every PubMed id and checks each rsID's GRCh38 position and alleles against dbSNP.

## Genotype Rules

Some results need several sites read together: APOE's ε alleles come from two sites, and the TAS2R38 taster and non-taster forms are haplotypes of three. Each such result is a rule in `assets/genotype_rules.tsv`, written as data, not code. A rule has four kinds of line:

| Line | Columns | What it says |
|---|---|---|
| `rule` | rule, gene, trait, needs_phase, ancestry, note, pmids | the rule and its sources |
| `site` | rule, rsid, chrom, pos, ref, alt | one site: GRCh38 position and plus-strand bases, as in dbSNP |
| `haplotype` | rule, name, alleles | a named haplotype: one base per site, in the order of the site lines |
| `combination` | rule, diplotypes, label, reading | what one pair of haplotypes reads as |

The two haplotypes a person's genotypes make are looked up by their alleles, and the pair gives the combination line. A heterozygous site does not say which of its alleles shares a copy with the alleles at the other sites, so:

- **`needs_phase no`** (APOE): the combination line is the one listing every pair the genotypes allow. Heterozygous at both APOE sites allows ε2/ε4 and ε1/ε3; that line reads ε2/ε4 and says it could also be the rare ε1/ε3.
- **`needs_phase yes`** (TAS2R38): with one heterozygous site or none, the two haplotypes follow from the genotypes. With two or more, WhatsHap phases them from the reads: a read, or the two reads of a pair, that covers two of the sites shows which alleles sit on one copy. When every heterozygous site of the rule lands in one phase block the rule is read (`phase` read_backed); otherwise it is `unphased`, not read. The three TAS2R38 sites lie within 741 bases, so a read pair can span them, but not always.

`status` of a rule line is `called`, `not_covered` (a site has no genotype), `other_allele` (an allele the rule does not describe), `other_haplotype` (a copy that is none of the rule's haplotypes; it names the alleles) or `unphased`. `phase` is `not_used` (needs_phase no), `not_needed`, `read_backed` or `unphased`.

The rules today:

| Rule | Sites | Haplotypes | Papers (PubMed) |
|---|---|---|---|
| APOE | rs429358, rs7412 | ε1, ε2, ε3, ε4 | 8346443, 9343467 |
| TAS2R38 | rs713598 (codon 49), rs1726866 (262), rs10246939 (296) | PAV (taster), AVI (non-taster) | 12595690, 15723792 |

Kim et al. (2003) found five TAS2R38 haplotypes worldwide. The rule reads the two common ones. A copy with AAV, AAI, PVI or a rarer one is `other_haplotype`, so the GRCh38 reference itself, which carries AAI (Ala49, Ala262, Ile296), is not read when every site matches it. The readings follow Bufe et al. (2005): PAV/PAV rated PTC and PROP the most bitter, PAV/AVI recognised them as bitter at nearly the same concentrations with a lower and more variable intensity, and AVI/AVI recognised PTC only at far higher concentrations.

### Checks on the rule table

- `bin/traits_report.py --check-table assets/traits.tsv --rules assets/genotype_rules.tsv` (also run by every step-42 run and by `tests/test_traits_report.py`): every haplotype allele is its site's ref or alt, the combination lines are exactly the pairs the genotypes can give (none missing, none that can never be read), and a site that is also a traits table row has the row's position and bases.
- `scripts/ci/check-rule-alleles.sh` checks every site against dbSNP: the rsID at that GRCh38 position, ref as dbSNP's REF and alt among its ALT alleles. A pair written the wrong way round, or as the minus-strand bases (TAS2R38 lies on the minus strand, where the papers name the codon bases), fails with the reason. The dbSNP records it compares with sit in the script, taken from dbSNP's GRCh38 VCF of the build `CONSUMER_DBSNP_BUILD` pins in `versions.env` (b157); `--online` reads them again over HTTPS and fails when they differ, and the e2e case `genotype-rules-1-phase` runs it. A C/G or A/T site is its own complement, so there a strand mistake keeps valid alleles: the check prints a note for each such site (rs713598 is one), and its haplotype alleles rest on the codon: Pro49 has C at the first base of codon 49 on the coding strand (NM_176817.5 c.145C, where Ala49 has G), so the PAV base on GRCh38 is its complement, G.

### Adding a rule

Write the rule, site, haplotype and combination lines from the primary papers. Take each site's position and bases from dbSNP's GRCh38 record (plus strand), and turn each haplotype's codon bases into plus-strand bases for a gene on the minus strand. A new site needs its dbSNP record in `scripts/ci/check-rule-alleles.sh`: run it with `--print-table` and add each printed line at the end of the script after `#: `. Then run `python3 tests/test_traits_report.py`.

## APOE

Many people choose not to learn their APOE result. Decide before you open the Traits section; the reports keep the APOE lines last, and the HTML report folds them away.

The APOE line is the APOE rule of `assets/genotype_rules.tsv`, written into the traits table after the two APOE rows; it is byte for byte the line the earlier hand code wrote (tests/test_traits_report.py compares the sha256 of every genotype pair). It combines rs429358 (Cys112Arg) and rs7412 (Arg158Cys) into the ε2, ε3 and ε4 alleles: Arg112 with Arg158 is ε4, Cys112 with Cys158 is ε2, Cys112 with Arg158 is ε3. The odds ratios it quotes for Alzheimer's disease, against ε3/ε3, are those of a meta-analysis of people of European descent from clinic and autopsy studies (Farrer et al., JAMA 1997): ε2/ε3 0.6, ε3/ε4 3.2, ε2/ε4 2.6, ε4/ε4 14.9. The ε4 effect was weaker in African American and Hispanic people and stronger in Japanese people. ε4 is a risk factor, not a diagnosis. Without phase a heterozygote at both sites reads as ε2/ε4, which could also be the rare ε1/ε3; a pair with the rare ε1 allele is not read.

[Step 25](25-prs.md) computes a polygenic score for late-onset Alzheimer's disease (PGS000334) whose two largest weights are the same two sites ([Before you look](interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you)). [Step 6](06-clinvar-screen.md) does not list APOE: it screens the ClinVar records classified as pathogenic or likely pathogenic, and on 2026-10-15 ClinVar classified rs429358 (variation 17864) as "conflicting classifications of pathogenicity; other; risk factor".

## Runtime

Seconds: bcftools reads the gVCF at the table's positions through its index (35 rows and 5 rule sites), and WhatsHap, when a rule needs it, reads the BAM over three sites within 741 bases.

## Notes

- A reading is what the cited studies found on average in the people they studied. It is not a diagnosis, and for most traits it says little about one person.
- Rows marked `ancestry` read differently in other populations; their note says how. Lactase persistence is the clearest case: -13910T explains it in Europeans, and the three African alleles have their own rows.
- Each traits table row is read alone. The genotype rules read several sites together: APOE, and the TAS2R38 taster haplotype with read-backed phase. Two different MC1R red-hair variants on the same or on different copies are still read one by one; their note says so.
- SNVs only; the table has no indel and no chrX, chrY or chrM row.
- Positions with no genotype (no reads, or a no-call) are listed as not read. On the e2e fixture, a small slice of the HG002 genome, only the table positions inside its slices are read, and the e2e case checks their genotypes against the GIAB truth.
- The licence: the readings, notes and table are this project's own text, GPL-3.0 like the rest of the repository.
