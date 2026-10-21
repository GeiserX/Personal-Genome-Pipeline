# Step 27: CPIC Drug-Gene Recommendation Lookup

## What This Does

Reads your PharmCAT report (step 7) and writes a plain-text list of the medications whose prescribing guidance depends on your result: for every gene where you are not a normal metabolizer, the drugs PharmCAT's own report matched to your phenotype, with the CPIC recommendation text. With pypgx output (step 32) it also compares the two callers and warns about a gene PharmCAT could not call while pypgx did. With the consensus table of [step 36](36-pgx-consensus.md) it lists the calls other tools gave PharmCAT (HLA-A and HLA-B from T1K, a CYP2D6 call pypgx and Cyrius agree on) and why a gene was held back. Opt-in: with a list of the medications you take, it also says when one of them is a strong or moderate CYP2D6 or CYP2C19 inhibitor, and for CYP2D6 gives the phenotype CPIC predicts on that drug (phenoconversion) beside your genotype phenotype.

## Why

PharmCAT produces detailed JSON and HTML reports, but finding the actionable parts takes time. This step distills them into one short file. The drug list comes from the PharmCAT report itself (its `drugs` section, maintained by PharmCAT with each data release), so a gene PharmCAT calls is never dropped because a table in this pipeline was missing it.

## Tool

- `bin/pgx_parse.py`, the one PharmCAT reader of the pipeline (Python standard library only). The Nextflow `CPIC_LOOKUP` module and the unit test `tests/test_cpic_parser.py` run the same file.

## Docker Image

- `PYTHON_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Input

- PharmCAT JSON report from step 7: the first `*.report.json` (or `*_pharmcat.json`) in `${GENOME_DIR}/${SAMPLE}/pharmcat/`, else in `${GENOME_DIR}/${SAMPLE}/vcf/`.
- Optional: `${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv` from step 32.
- Optional: `${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_pgx_consensus.tsv` from step 36.
- Optional, opt-in: `MEDICATIONS=FILE`, a plain text list of the drugs you take, one per line. Write the generic name first (`paroxetine 20 mg` is read as `paroxetine`; a brand name is not recognised). `#` starts a comment and blank lines are skipped. The file can sit anywhere; it is mounted read-only. A missing list, or a missing `assets/pgx_inhibitors.tsv`, stops the step. Nothing is ever inferred: without the list there is no phenoconversion. In the Nextflow pipeline the same list is `--medications FILE`. `run-all.sh` does not take one: run this step again with `MEDICATIONS` afterwards.

## Command

```bash
./scripts/27-cpic-lookup.sh your_name

# With a medication list (opt-in)
MEDICATIONS=/path/to/medications.txt ./scripts/27-cpic-lookup.sh your_name
```

## What the Script Does Internally

1. Locates the PharmCAT JSON report.
2. Runs `bin/pgx_parse.py cpic-report` in the Python image. It reads the gene calls (PharmCAT 3.x flat `genes` map, the 2.x map nested by source, and the older list are all read) and writes the phenotype table.
3. A gene counts as non-normal only when its phenotype names a changed metabolism or function (poor, intermediate, rapid or ultrarapid metabolizer; decreased, increased, poor or no function), a positive HLA test, a G6PD deficiency, malignant hyperthermia susceptibility, an increased MT-RNR1 risk, an ivacaftor-responsive CFTR, or a VKORC1 `-1639` genotype with the A allele (`-1639 GA` or `-1639 AA`, the warfarin dose-reduction case). A called gene whose phenotype is anything else (`n/a`, `no phenotype assigned`, the reference VKORC1 genotype `-1639 GG`, `Uncertain Susceptibility`, `Indeterminate`) is `unclassified`: it is listed in the section "Called Genes Without a Function Phenotype", is not counted, and gets no drug list. One exception: PharmCAT can leave a gene `Indeterminate` at gene level while the annotations of its drugs assign it a phenotype (for example a gene with unphased heterozygous positions plus positions missing from the VCF). When the drug annotations PharmCAT gives that gene (those for any of its possible diplotypes when there are any) all name the same phenotype, other than `Indeterminate`, `No Result` or `n/a`, the gene takes it, with the note `[gene level Indeterminate; phenotype from PharmCAT's drug annotations]` in the phenotype table and the text, and its status follows that phenotype. An annotation with no phenotype for the gene is not counted. Several possible diplotypes must also get the same drug guidance, because the drug list shows the first one's. With no annotation that names a phenotype, with annotations that disagree (one that says `Indeterminate` disagrees too), or with possible diplotypes whose guidance differs, it stays `Indeterminate` and `unclassified`. PharmCAT's own example report has 3 non-normal genes and 7 unclassified ones. For each gene with a non-normal phenotype it lists the drugs from the report's `drugs` section: CPIC's recommendation for the called diplotype per drug first, then the drugs DPWG or the FDA name. PharmCAT lists an annotation for every diplotype the sample may have, so only the ones for the called diplotype (or, without diplotype labels, its phenotype) are shown. When the report names no drug for the gene, it falls back to the gene's `relatedDrugs`, then to a small static table, and finally prints a line saying the gene is not in the drug table, so a gene is never skipped silently.
4. When PharmCAT lists more than one possible diplotype for a gene and their phenotypes differ (positions missing from the VCF leave it unable to choose), the gene is `ambiguous`: it is listed in its own section with the possible phenotypes and no drug guidance, never as its first diplotype.
5. With the pypgx summary it writes the comparison table and, for a gene PharmCAT reports as not called or ambiguous but pypgx called (CYP2D6 is the usual one), a warning in the recommendations naming every drug PharmCAT links to that gene. When pypgx's call is normal the line is a note instead.
6. With step 36's consensus table it writes a section "Calls From Other Tools": each gene passed to PharmCAT as an outside call, from which tool, and whether PharmCAT's report shows it as one (`callSource` `OUTSIDE`); and each gene held back, with the reason and what each caller said. For a held-back CYP2D6 it names the drugs CYP2D6 affects and says no guidance is given for them, and the pypgx-only warning of point 5 is not printed for it: one caller's call is what the consensus refused. Genes PharmCAT reports as outside calls carry `[outside call]` in the gene results. An outside call of CYP2D6, HLA-A or HLA-B in PharmCAT's report that the table does not pass on, with the same diplotype, comes from an earlier step 36 run that step 07 has not replaced: it gets no drug guidance, is listed as not called, and a warning says to rerun step 07. When the table is missing (a caller step ran again since) but PharmCAT's report has outside calls of those genes, the step passes the missing table, so none of them is confirmed. MT-RNR1's outside call comes from step 20's chrM calls through step 07, never from step 36: the table never lists it, so it is not checked against it and keeps its drug guidance.
7. With a medication list, each line is matched by its first word (leading commas and semicolons ignored; a line of only those is skipped) against `assets/pgx_inhibitors.tsv`: the strong and moderate inhibitors of CYP2D6 and CYP2C19 as the FDA's table of clinical CYP inhibitors classes them, each row with its source. For CYP2D6 the step applies CPIC's rule from its opioid guideline: "for patients on strong CYP2D6 inhibitors, the CYP2D6 activity score is adjusted to 0 and the predicted phenotype is a poor metabolizer", and for a moderate inhibitor "the activity score is multiplied by 0.5 and then converted to the predicted phenotype" (CPIC's bands: 0 poor, under 1.25 intermediate, 1.25 to 2.25 normal, above 2.25 ultrarapid). The activity score is the one PharmCAT gives the called diplotype. With several inhibitors of one gene the strongest counts. The result is written as `Poor Metabolizer (phenoconverted by paroxetine)` beside the genotype phenotype, never in its place, with a `WARNING:` line in the section "Medications and Phenoconversion", so the summary carries it. Nothing is converted for a gene that is not called, ambiguous or unclassified, or for a moderate inhibitor when the report gives no activity score: the line says why. CPIC has no activity score for CYP2C19, so a CYP2C19 inhibitor gets a note and no converted phenotype. A line that names no drug of the table is listed as "Not in the inhibitor table", never dropped. A list that names no drug converts nothing and leaves the phenotype table as it is without a list. An unreadable list or table is reported in that section; the rest of the report is written as usual.
8. A report that cannot be read, or that yields no gene, writes a "PARSING FAILED" report and the step exits 1. It never writes an all-clear report from a format it could not read.

The comparison used to be written by step 32. It moved here because `run-all.sh` starts steps 7 and 32 side by side, so step 32 could read a missing or previous-run PharmCAT report; step 27 runs after both.

## Output

| File | Contents |
|---|---|
| `cpic/${SAMPLE}_cpic_recommendations.txt` | Gene results, the medications for each non-normal gene, the called genes without a function phenotype, the calls from other tools (step 36), uncallable genes, and the pypgx warnings |
| `cpic/${SAMPLE}_phenotypes.tsv` | One row per gene: `Gene`, `Diplotype`, `Phenotype`, `Status` (`normal`, `non-normal`, `unclassified`, `ambiguous` or `not called`). With a medication list that names a drug, a fifth column `Phenoconverted_phenotype`: the converted phenotype, `not converted: <why>`, `not computed: ...` (CYP2C19), or `-` when no drug on the list inhibits the gene. `Phenotype` and `Status` stay the genotype's |
| `pypgx/${SAMPLE}_pharmcat_comparison.tsv` | PharmCAT and pypgx diplotypes side by side (only when step 32 ran) |

<figure markdown="span">
  ![The top of the step 27 text report: the sample, the PharmCAT version and the gene results table with each gene's diplotype and phenotype](images/demo-cpic-report.png){ loading=lazy }
  <figcaption>The top of the report for DEMO-001, an invented sample. The diplotypes and phenotypes are made up and are not anyone's result. HLA-A, HLA-B and CYP2D6 are marked as outside calls: step 36 passed them to PharmCAT from the BAM-based callers. Below this table the report lists the medications for each gene with a non-normal phenotype, the called genes without a function phenotype, the outside calls and the genes that could not be called.</figcaption>
</figure>

## Runtime

3 to 5 seconds, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)).

## Interpreting Results

### Gene Results Table

Lists every pharmacogene with its called diplotype and phenotype (placeholders, not a result):

```
Gene         Diplotype                      Phenotype
CYP2C19      *x/*y                          <phenotype>
CYP2D6       *x/*y                          <phenotype>
```

### Affected Medications

Only genes with a non-normal phenotype (point 3 above) appear here, each with the drugs and the CPIC recommendation PharmCAT matched to your result. A called gene whose phenotype names no changed function (`n/a`, `no phenotype assigned`, the reference VKORC1 genotype `-1639 GG`, `Uncertain Susceptibility`) is `unclassified`: it gets no drug list and is listed under "Called Genes Without a Function Phenotype". A variant allele in such a gene can still matter, so check PharmCAT's HTML report for it. Genes that could not be called are listed separately at the end. Their absence from the medications section does NOT mean normal function.

### Calls From Other Tools

When step 36 ran, this section lists HLA-A, HLA-B and CYP2D6: passed to PharmCAT (then their drug guidance is in the sections above, as for any gene PharmCAT called) or not, and why. `indeterminate` for CYP2D6 means the callers did not agree, only one ran, or the depth at CYP2D6 could not be trusted: no drug guidance is given for CYP2D6 then, and the drugs it affects are listed so you know what is not covered.

### Medications and Phenoconversion

Only with a medication list. A drug can make you metabolize like a different phenotype than your genes predict: on a strong CYP2D6 inhibitor a CYP2D6 normal metabolizer is predicted to be a poor metabolizer. This section names each inhibitor on your list, your genotype phenotype, the converted one with CPIC's rule, and the drugs the gene affects. The converted phenotype holds while you take the inhibitor. The drug guidance in "Affected Medications" stays the one for your genotype; read the CPIC guideline with the converted phenotype for the drugs listed here. Your list is only as complete as you made it, and the table holds strong and moderate inhibitors only: a weak inhibitor, a brand name or a drug missing from the table is not checked. This is an aid for your clinician, not advice: do not start, stop or change a medication because of it.

### PharmCAT and pypgx

When step 32 ran, this section names each gene PharmCAT could not call but pypgx did, with pypgx's call and the drugs it affects. Read those drugs with the pypgx call and [docs/32-pypgx.md](32-pypgx.md).

### What to do with the results

1. Check whether you take (or might be prescribed) any of the listed medications.
2. For any match, read the full CPIC guideline at [cpicpgx.org/guidelines](https://cpicpgx.org/guidelines/).
3. Share the report with your prescribing physician or pharmacist.

## Limitations

- The drug guidance is the one bundled with the pinned PharmCAT release; newer CPIC guidelines arrive with a PharmCAT update.
- The static fallback table is used only when the report names no drug for a gene.
- CYP2D6 has drug guidance here only when pypgx (step 32) and Cyrius (step 21, opt-in) agree and the depth check passed (step 36). Without Cyrius it is always held back.
- Phenoconversion covers CYP2D6 and CYP2C19 inhibitors in `assets/pgx_inhibitors.tsv` only (strong and moderate, per the FDA's table), and computes a converted phenotype for CYP2D6 only. Inducers and other genes are not covered.
- This is NOT medical advice. Always consult a healthcare professional before making medication changes.

## Links

- [CPIC Guidelines](https://cpicpgx.org/guidelines/)
- [PharmCAT](https://pharmcat.org/)
- [PharmGKB](https://www.pharmgkb.org/)
- [CPIC Gene-Drug Pairs Table](https://cpicpgx.org/genes-drugs/)
- [CPIC guideline for CYP2D6, OPRM1 and COMT and opioids](https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/) (the CYP2D6 inhibitor rule)
- [FDA: drugs that interact with CYP enzymes](https://www.fda.gov/drugs/drug-interactions-labeling/healthcare-professionals-fdas-examples-drugs-interact-cyp-enzymes-and-transporter-systems) (strong and moderate inhibitors)
