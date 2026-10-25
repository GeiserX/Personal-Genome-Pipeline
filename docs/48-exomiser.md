# Step 48: Symptom-Driven Variant Ranking (Exomiser)

## What This Does

Ranks the genes and variants of one person's genome by how well they could explain a list of symptoms. You give the symptoms as Human Phenotype Ontology (HPO) terms; [Exomiser](https://github.com/exomiser/Exomiser) scores each rare coding or splice variant by how damaging it looks, scores each gene by how closely its known diseases and animal models match the terms, and ranks the genes under each mode of inheritance. The step reports the top 10 ranks: gene, score, mode of inheritance, the variants behind the score, Exomiser's ACMG class for each, and the same in plain words.

Opt-in: it runs only when you give HPO terms.

**This is a ranking to bring to a clinical geneticist, not a diagnosis.** A high score means a rare variant, predicted to be damaging, in a gene linked to symptoms like the ones you gave. It does not show that the variant causes them. A low score or an empty list does not rule a genetic cause out.

## Why

The other steps report what is known to be pathogenic (ClinVar, CPSR) or what looks damaging (steps 23 and 31), with no idea why you are looking. When there is a question, such as a child with a set of symptoms or a condition that runs in the family, the useful order is the one that puts the genes that match those symptoms first. That is what Exomiser was built for: it combines the variant evidence with phenotype matching against human disease annotations and model organisms ([Smedley et al. 2015, Nature Protocols](https://doi.org/10.1038/nprot.2015.124); the cross-species matching, hiPHIVE, is [Robinson et al. 2014, Genome Research](https://doi.org/10.1101/gr.160325.113)). Its documentation reports that on 742 solved cases of the 100,000 Genomes Project the diagnosis was among the top 5 candidates in 94% of them ([Exomiser 15.1.0, ACMG Assignment](https://exomiser.readthedocs.io/en/15.1.0/acmg_assignment.html)). Those were patients with a rare disease and a clinician's choice of terms; a healthy person's genome with self-chosen terms is a different setting.

## Tool

- Exomiser 15.1.0 (AGPL-3.0), command `exomiser analyse`, with Exomiser's exome preset: PASS calls only; variants in UTRs, introns and between genes removed; a 2% ceiling on the highest population frequency (UK10K and gnomAD exomes and genomes); REVEL, MVP, AlphaMissense and SpliceAI for pathogenicity; the inheritance filter; then the OMIM and hiPHIVE prioritisers.
- `bin/exomiser_report.py` (standard library Python): checks the terms against the HPO release inside Exomiser's data, writes the analysis file, and reads Exomiser's variants table into the report.

The genome preset of Exomiser also scores non-coding variants, but it needs the separate REMM data and its regulatory features; this step does not install them.

## Docker Image

- `EXOMISER_IMAGE` (the `-bash` build of Exomiser 15.1.0), plus `PYTHON_IMAGE` for the terms and the report

The `-bash` build is used because the plain `15.1.0` tag is distroless: it has no shell, and Nextflow needs one in the task container. Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. The image is optional: `setup.sh --exomiser-data` pulls it.

## Data

Exomiser's data release `EXOMISER_DATA_VERSION` (2602, from 13 May 2026), installed by

```bash
./scripts/setup.sh --exomiser-data ${GENOME_DIR}
```

into `${GENOME_DIR}/reference/exomiser-2602/`, with `2602_hg38/` (variant frequencies, pathogenicity scores, transcripts) and `2602_phenotype/` (HPO, disease and model organism annotations). The two zip files are 37 GB together (23.6 GB and 13.5 GB) and unpack to 58 GB (32.7 GB and 25.7 GB). Setup downloads and unpacks one at a time and deletes each zip once it is unpacked, so it needs about 72 GB free while it runs. Exomiser publishes no checksums, so `versions.env` pins the sha256 of each zip as we measured it on 25 October 2026 (`EXOMISER_HG38_SHA256`, `EXOMISER_PHENOTYPE_SHA256`). A broken download resumes when setup runs again. The 2602 release notes list Exomiser 14.0.0 to 15.0.0 as compatible; the 15.1.0 release notes name no change to the data format.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.vcf.gz`: one sample, GRCh38
- the HPO terms, in `HPO` (bash) or `--hpo` (Nextflow)
- for the family mode, a step 38 folder (see below)

## Command

```bash
HPO=HP:0001250,HP:0001263 ./scripts/48-exomiser.sh your_name
```

With `run-all.sh`: set `HPO=...` and the run adds step 48 (it is reported as skipped when the data is not installed, or when `TOOLS` names `exomiser` and `HPO` is not set). In the Nextflow pipeline: `--hpo HP:0001250,HP:0001263 --exomiser_data ${GENOME_DIR}/reference/exomiser-2602` (processes `EXOMISER_PREPARE`, `EXOMISER` and `EXOMISER_REPORT`). `--hpo` turns the tool on; `exomiser` in `--tools` without `--hpo` stops the run. With more than one sample in the samplesheet, name the person the terms describe with `--exomiser_proband`.

Exomiser gets 16 GB of memory by default (`EXOMISER_MEM`, in GB; Java gets three quarters of it).

## How to pick the HPO terms

An HPO term is an id like `HP:0001250` (Seizure) for one clinical feature. Find them in the [HPO browser](https://hpo.jax.org): search for the feature in plain words, open the term, and copy its id. The ontology is described in [Gargano et al. 2024, Nucleic Acids Research](https://doi.org/10.1093/nar/gkad1005).

- **Use the most specific term that is true.** `HP:0011304` (Broad thumb) matches better than `HP:0001167` (Abnormality of the thumb). Exomiser gives credit for a near match, but a precise term separates genes better.
- **Give the features that stand out, usually 3 to 10.** A feature most people have, or one with an obvious other cause, adds noise rather than signal.
- **Leave out what is not there.** The step searches the terms you give; it does not read absent features.
- **Clinical notes help.** If a doctor has written the features down, the words in the notes are the ones to search the browser for. A clinical geneticist can choose the terms with you.

The step checks every term against the HPO release inside Exomiser's data. A term that is not HP: and seven digits, that the release does not know, or that it marks obsolete stops the step and says which (with the replacement the release names, if any). A term given by an old id that the release lists as an alternative is used under its current id, and the terms file keeps the id you gave. `${SAMPLE}_exomiser_terms.tsv` lists each term with its name, so you can check that the ids say what you meant.

## Family mode (bash only)

With the genotypes of relatives Exomiser can tell which variants fit the way the condition runs in the family: a variant shared by an unaffected parent is less likely to cause a dominant condition, and a de novo one is more likely. Run step 38 on the family first, then:

```bash
HPO=HP:0001156,HP:0001363 EXOMISER_FAMILY=my_family EXOMISER_UNAFFECTED=mother,father ./scripts/48-exomiser.sh child
```

`EXOMISER_FAMILY` is the step 38 folder (`family.vcf.gz` and `pedigree.tsv`); the sample on the command line is the person the terms describe, marked affected. `EXOMISER_AFFECTED` names other members who have the condition, `EXOMISER_UNAFFECTED` those who do not, and everyone else is unknown, which Exomiser does not filter on. The pedigree of step 38 has no affected column, so the step asks rather than guesses. The PED file it writes is `${SAMPLE}_exomiser.ped`. The Nextflow pipeline runs the single-sample mode only.

## Output

| File | Contents |
|---|---|
| `exomiser/${SAMPLE}_exomiser_top.tsv` | the top 10 ranks, one row per contributing variant: `rank`, `gene`, `moi` (AD, AR, XD, XR, MT), `exomiser_score`, `phenotype_score`, `variant_score`, `p_value`, `variant` (contig:position REF>ALT), `genotype`, `consequence`, `hgvs`, `acmg`, `acmg_evidence`, `disease`, `clinvar`, `max_frequency` |
| `exomiser/${SAMPLE}_exomiser.txt` | the same in words, with the terms searched and the note on what the ranking is |
| `exomiser/${SAMPLE}_exomiser_terms.tsv` | the terms searched: `hpo_id`, `label`, `given_as` |
| `exomiser/${SAMPLE}_exomiser.variants.tsv`, `.genes.tsv`, `.html` | Exomiser's own files, for the top 50 genes |

The scores run from 0 to 1. The phenotype score is how well the gene's known diseases and models match the terms; the variant score is how rare and damaging its variants look; the Exomiser score combines both with a model trained on solved cases of the 100,000 Genomes Project ([Exomiser 15.1.0, Interpreting the Results](https://exomiser.readthedocs.io/en/15.1.0/result_interpretation.html)). A gene can appear once per mode of inheritance. The p-value is Exomiser's estimate against random scores it draws afresh on every run ([CombinedScorePvalueCalculator](https://github.com/exomiser/Exomiser/blob/15.1.0/exomiser-core/src/main/java/org/monarchinitiative/exomiser/core/analysis/score/CombinedScorePvalueCalculator.java)), so it changes slightly between two runs of the same VCF and terms; the scores and the ranks do not.

The ACMG class is Exomiser's automatic, partial use of the ACMG/AMP variant classification standard ([Richards et al. 2015, Genetics in Medicine](https://doi.org/10.1038/gim.2015.30)): it applies the criteria its data can decide, and Exomiser's documentation says it is not a substitute for a clinical geneticist's classification ([ACMG Assignment](https://exomiser.readthedocs.io/en/15.1.0/acmg_assignment.html)).

Both reports show the top ranks as their "Symptom-Driven Ranking (Exomiser)" section, with the terms and the note.

## Runtime

Not measured on a whole genome here. In CI (a 4-core GitHub runner, `EXOMISER_MEM=6`), one run of step 48 on Exomiser's test release took 24 to 30 seconds, for the Pfeiffer exome alone and for its family of four alike (see Notes). The full hg38 release reads far more data, so expect longer.

## Notes

- What it does not see: structural variants, repeat expansions, mitochondrial heteroplasmy, and variants outside coding and splice regions. Steps 4, 9 and 20 cover the first three.
- The result depends on the terms. Different terms give a different ranking; try a second set if the first is unsure.
- The bash step runs Exomiser with no network: the terms and the results stay on your computer.
- CI runs Exomiser for real, on Exomiser's own test release (hg19, 2.8 GB) and the Pfeiffer syndrome exome from its repository, with the four terms its documentation gives for it (`EXOMISER_DATA` and `EXOMISER_ASSEMBLY=hg19` point the step at that release; only CI sets them). The full hg38 release is too large for a CI runner. The case checks that FGFR2 ranks first with the planted variant, that unrelated terms lower FGFR2's symptom match, that in the family run marking the sibling without the variant as affected removes FGFR2 from the dominant ranks, and that the Nextflow modules write the same files as the bash step.
