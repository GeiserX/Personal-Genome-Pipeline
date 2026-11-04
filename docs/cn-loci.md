# Copy Number from Depth (the locus catalog)

Several steps count copies of a gene from short-read depth: HBA1 and HBA2 in step 35, the RHD deletion in step 41 and the LPA KIV-2 repeat in step 46. Step 35 also counts three gene-dosage loci that no pinned caller models (HP, LCE3B/LCE3C and CCL3L1), CYP21A2 from depth and its pseudogene variant sites, and the interval rows: the recurrent CNV regions ClinGen curates (22q11.2, 1q21.1, 15q11.2, 16p11.2, 17p12 and the rest), three carrier deletion alleles and the Y chromosome AZF regions. They all use one engine, `bin/copy_number.py`, and one table, `bin/cn_loci.tsv`, with one row per locus. A new locus is a new row, not a new script.

The table sits in `bin/` because `bin/` is the one repository folder that every container mounts: the bash steps mount it as `/pgp-bin`, and Nextflow puts it on every task's PATH.

## The method

`samtools bedcov` sums the read depth over the row's segments. Take T as the depth over the target segments and C as the depth over the single-copy control segments. Then:

    estimate = scale x T / C - offset

The control segments carry two copies on a diploid genome, so `scale` turns the ratio into copies. For HBA, `scale` is 2 and T adds up the two genes, so four genes give 4. For the LPA KIV-2 block, `scale` is 14 (seven reference exon pairs on two chromosomes) and `offset` is 2 (one single-copy pair on each chromosome).

Before it calls, the engine checks, in this order:

1. **No data**: both bedcov tables are empty (the BAM has no reads on that contig): not assessed, with exit code 0.
2. **Broken control**: a control segment has no read while the other controls have at least the row's minimum depth. The engine stops with `BROKEN_CONTROL` and exit code 3. Either the row's coordinates are wrong, or the sample lacks that stretch. A number from such a run would mean nothing.
3. **Control depth** under `min_control_depth`: not assessed. With `control_depth` `flanks`, each flank must reach it.
4. **Flanks disagree** by more than `max_control_ratio` (`flanks` only): not assessed. A deletion that reaches one flank would skew the ratio.
5. **MAPQ 0 share** of the target depth under `mapq0_min`: not assessed. The reference holds every copy of a duplicated gene, so an aligner that reports every read gives reads that fit several copies mapping quality 0. With no such read, the aligner left them out (novoalign does by default) and the count would come out low. The check is skipped when a row without an offset rounds to no copy, since no read is left to judge by. A share over `mapq0_max` means the target is mostly misplaced reads: not assessed.
6. **Below the floor**: a row with an offset whose estimate is 0 or less: not assessed.

Then the row's `call` column decides:

- `round:D`: a whole number, when the estimate lies within D of one (HBA, D = 0.3);
- `cut:A,B`: the number of cut points the estimate reaches (RHD: 0 under 0.5, 1 under 1.5, else 2);
- `psv:D`: the module total of a gene and its pseudogene as `round:D` does, then the reads at the row's PSV sites ([the CYP21A2 row](#the-cyp21a2-row));
- `class:LOW,HIGH`: low, typical or high from the 95% interval (LPA KIV-2, 30 and 45). With `interval` `cv`, that interval comes from the spread of the per-segment depths: relative standard error sqrt(cv_target^2 / n_target + cv_control^2 / n_control), applied on the log scale. It covers sampling noise only.

`copy_number.py call --json` writes the result as one JSON block per locus: status, reason, estimate, copies, interval, class, the depths and the reading.

## The columns

| Column | Meaning |
|---|---|
| `locus` | the row's name, as `--locus` takes it |
| `kind` | `gene`, `interval` (a stretch between segmental duplications, controls outside them) or `organelle` |
| `chrom` | GRCh38 contig of every segment |
| `targets`, `controls`, `details` | `NAME:start-end[:label]`, separated by `;`, in BED coordinates (0-based start, end exclusive). Details get a copy number each but take no part in the call (the RHD exons); in a `psv` row they are the PSV sites, one base each, labelled `PSEUDO_POS GENE>PSEUDO text` |
| `target_depth` | `sum` (each target segment is one reference unit; their mean depths add up) or `pooled` (one length-weighted mean) |
| `control_depth` | `flanks` (the mean of the flank means, each flank checked) or `pooled` |
| `scale`, `offset` | as in the formula above |
| `bedcov` | the samtools bedcov options of each pass, `,` between passes, `.` for none. The second pass is the MAPQ >= 1 pass the MAPQ 0 share comes from. `tests/test_copy_number.py` checks that each step and module runs exactly these passes |
| `mapq0_min`, `mapq0_max` | the MAPQ 0 share guards, or `.` |
| `min_control_depth`, `max_control_ratio` | the control guards, or `.` |
| `call`, `interval` | as above |
| `readings` | `KEY=text`, separated by `\|`: the text for a count, a class, or `more` for any count above the keys |
| `input` | `bedcov` (this engine) or `cnvpytor-level` (the row records thresholds on CNVpytor's depth level, where 1.0 is two copies) |
| `owner` | the engine that calls the locus (`cn-depth`, this one) |
| `cross_check` | another engine that also reports the locus (`parascopy:RHCE`, `paraphase:hba`, `cnvpytor:cnvs` for step 18's CNVpytor calls) or `none` |
| `pmids` | the papers the method of the row comes from, or `.` |

## The loci

| Locus | Step | Targets against controls | Call |
|---|---|---|---|
| HBA | 35 | HBA2 and HBA1 gene bodies against two 40 kb flanks outside the `--SEA` deletion | whole number of alpha-globin genes ([step 35](35-paralogs.md)) |
| LPA_KIV2 | 46 | the 14 exons of the KIV-2 block against 15 single-copy LPA exons | repeat count with a 95% interval, low, typical or high ([step 46](46-lpa.md)) |
| RHD | 41 | the RHD gene against two flanks outside the common RHD deletion; the ten exons as details | 0, 1 or 2 copies ([step 41](41-bloodtype.md)) |
| HP | 35 | the two reference copies of the duplicated HP exons (exons 3-4 and 5-6) against flanks outside HP and HPR | copies of the duplicated exons: 2 for Hp1-1, 3 for Hp2-1, 4 for Hp2-2 |
| LCE3BC | 35 | LCE3C to LCE3B, inside the LCE3C_LCE3B deletion, against flanks outside it | 0, 1 or 2 copies: the deletion on both, one or no chromosome |
| CCL3L1 | 35 | CCL3L3, the one reference copy of CCL3L1/CCL3L3, against single-copy flanks outside the duplications | total copies of CCL3L1 and CCL3L3, with no reading |
| CYP21A2 | 35 | CYP21A2 and the matching stretch of CYP21A1P against flanks outside the RCCX duplication; six PSV sites as details | RCCX module total, CYP21A2 copies at each site, and a deletion or a possible conversion ([the CYP21A2 row](#the-cyp21a2-row)) |

### The gene-dosage rows

Step 35 (and the PARASCOPY module) runs these three rows one after the other and writes one line per locus to `<sample>_cn_dosage.tsv`: locus, status, copies, estimate, reason and reading. They give total copies only, never which copy or haplotype. A locus whose check fails is not assessed there, with reason `check_failed`, and the step goes on. Every row runs two bedcov passes, all reads and MAPQ >= 1, as the HBA row does; none of the three sets a MAPQ 0 guard, and the dosage table does not keep the MAPQ 0 share. The coordinates come from Ensembl 116, and every flank lies outside UCSC's segmental duplication track (`genomicSuperDups`), both checked on 2026-11-03.

- **HP**: the reference carries the Hp2 allele, the exon 3-4 duplication that Hp1 lacks (Boettger 2016, PMID 26901066): exons 3-4 and 5-6 in Ensembl's HP-201. The two copies are 98% alike, so a read can land on either. The row adds both up, as HBA does. 2 copies read as Hp1-1, 3 as Hp2-1, 4 as Hp2-2. A deletion of the whole gene (Hp-del) on one chromosome gives the same 2 as Hp1-1 when the other carries Hp2, and the reading says so.
- **LCE3BC**: the 32 kb LCE3C_LCE3B deletion (de Cid 2009, PMID 19169253) removes both genes, so the stretch from LCE3C to LCE3B lies inside it. 0 copies means the deletion on both chromosomes, a common genotype the paper links with psoriasis.
- **CCL3L1**: the no-ALT reference holds one copy of CCL3L1/CCL3L3 (CCL3L3), so every read of either lands there and the row counts both together. The count goes back to Gonzalez 2005 (PMID 15637236), but the disease links reported for it are disputed, so no reading is given. The test fails when a CCL3L1 reading names a disease.
- **DEFB4** has no row. The beta-defensin repeat that carries it (DEFB4B to DEFB103B) lies in the paralog copies of SPAG11A's Parascopy model, so Parascopy owns it under the rule below. The test plants a DEFB4 row there and checks that the owner check catches it.

The report does not show `<sample>_cn_dosage.tsv` yet; the table and the step log are where to read it.

### The CYP21A2 row

CYP21A2 and its pseudogene CYP21A1P lie in the duplicated RCCX module, about 98% identical (UCSC `genomicSuperDups` pairs chr6:31,980,532-32,013,273 with chr6:32,013,270-32,046,127). No pinned caller models CYP21A2 on short reads; Paraphase owns it on HiFi reads (step 35b) and is the row's `cross_check`, and C4 stays with Parascopy's C4A model, so both targets start where that model ends: CYP21A2 from chr6:32,038,951 to its last base (NM_000500.9) and the same stretch of CYP21A1P, chr6:32,006,221-32,008,909 (NR_040090.1), matched base by base. The flanks, chr6:31,945,001-31,975,000 and chr6:32,070,001-32,100,000, lie outside the duplication and every `genomicSuperDups` pair. Coordinates from UCSC's `ncbiRefSeq` and `genomicSuperDups` tracks, checked on 2026-11-04.

The call has two parts:

1. **The module total**: `scale` 2 times the two targets' depths added up over the flank depth, a whole number within 0.3, else not assessed (`between_counts`). Two modules per chromosome give 4.
2. **The PSV sites**: `copy_number.py bed --psv --locus CYP21A2` prints each site on CYP21A2 and on CYP21A1P, and `samtools mpileup -A -B -q 0 -Q 13 -d 100000` without a reference counts the bases there (every read whatever its mapping quality, bases of quality 13 or more; read starts, ends and indels are skipped). A read with CYP21A2's base counts for CYP21A2 at either position, since a converted copy's reads may land on the pseudogene. A site with 30 such reads or more gets `estimate = module total x CYP21A2 reads / (CYP21A2 + CYP21A1P reads)` and its ratio to the median of the sites.

The sites, each CYP21A2's base in the reference with CYP21A1P's base at the matching position (dbSNP positions, the CYP21A1P bases read from the reference, checked on 2026-11-04): In2G c.293-13A/C>G (rs6467, chr6:32,039,081; CYP21A2 carries A or C there, CYP21A1P G), I172N c.518T>A (rs6475), V237E c.713T>A (rs12530380, the exon 6 cluster), V281L c.844G>T (rs6471), Q318X c.955C>T (rs7755898) and R356W c.1069C>T (rs7769409).

`cyp21a2_copies` is the least count over the sites. The `finding` is `deletion` when the module total is under 4 and that least count under 2 (the 30 kb chimeric deletion leaves one module fewer and a chimera with CYP21A1P's sequence at its 5' sites), `conversion` at each site whose estimate lies 0.6 copies or more under the median's count, or under 2, `pseudogene_deletion` when the module total is under 4 with 2 copies everywhere, and `none` otherwise. At 30x a site has about 70 reads, and an estimate of 2 copies there has a standard deviation of about 0.25, so the 0.6 margin keeps a false conversion on a site under 1 in 100. Which chromosome carries a change is not claimed, and a two-copy result does not rule out two copies on one chromosome and none on the other. The row has no MAPQ 0 guard yet: an aligner that leaves out the reads that fit both genes (novoalign by default) would lower the module total and could read as a deletion, and no such BAM has been measured here to set a threshold. Use a BAM from minimap2 or BWA-MEM. Step 35 writes the table to `<sample>_cyp21a2_psv.tsv` ([step 35](35-paralogs.md#cyp21a2-from-depth-and-psv-sites)); the couple screen takes it with `--cyp21a2`.

### The interval rows

`scripts/ci/build-cn-intervals.py` writes these rows, everything in `bin/cn_loci.tsv` after its marker line; do not edit them by hand. Step 35 (and the PARASCOPY module) runs `copy_number.py bed --kind interval`, one `samtools bedcov` over every segment of every row, then `copy_number.py intervals`, which writes one line per row to `<sample>_cn_intervals.tsv`: locus, name, chrom, start, end, status, copies, expected, change (`deletion`, `duplication`, `no change`, or `.` when not called), estimate and reason. The step log lists the rows with a change. The rows give a total copy count only, never breakpoints.

- **The ClinGen recurrent regions**: every region row of ClinGen Dosage Sensitivity's recurrent CNV BED, version 2.1, GRCh38 (CC0, pinned by `CLINGEN_RECURRENT_CNV_URL` and its sha256 in `versions.env`). The file's breakpoint and LCR rows are the segmental duplications that bound the regions. They are cut out of the targets, as are GIAB v3.6 segmental duplications, assembly gaps, satellites and the regions of every pinned Parascopy model (CHRNA7 in 15q13.3 stays with Parascopy). Targets are cut to pieces of at most 500 kb, and their depth is pooled. Controls are two 100 kb windows, one on each side of the cluster of ClinGen rows the region belongs to, clear of all of that and of the targets of the hand rows below (so no row's control sits on another row's target); a side with no clear window within 3 Mb (a chromosome end, as next to 22q11.2 and 15q11.2) is left out. The BED spells its chrX rows `chrx`; the script normalises them, and the engine refuses any contig name that is not a GRCh38 chromosome. The classification of a CNV found here follows the ACMG and ClinGen technical standards (Riggs 2020, PMID 31690835); a call from this step is a lead for a clinical test, not a classification.
- **Two more ClinGen curations**: the 1p36 terminal region (ISCA-37434, haploinsufficiency score 3) and the Xq28 region that holds MECP2 (ISCA-46304, triplosensitivity score 3: the MECP2 duplication), from ClinGen's region curation list, read on 2026-11-04.
- **Carrier deletion alleles**, at their ClinVar GRCh38 spans: del(GJB6-D13S1830) (variation 5546; del Castillo 2002, PMID 11807148), CFTRdele2,3 (variation 66105) and a 621 bp HBB deletion that removes exon 3 (variation 3769040). The alpha-globin `--SEA` deletion is read by the HBA row: two alpha-globin genes, with the reading that names `--SEA`. The 619 bp HBB deletion and the Sicilian delta-beta deletion have no row: we could not confirm their GRCh38 spans in a primary source.
- **AZF** (Y chromosome microdeletions, Krausz 2024, PMID 37674303): AZFa as USP9Y and DDX3Y, AZFb as KDM5D and EIF1AY (single-copy genes the deletion removes), and AZFc as the four DAZ copies added up, against ZFY and NLGN4Y. A male has 1 copy of AZFa and AZFb and 4 DAZ copies; a complete AZFc deletion leaves 0, the partial gr/gr and b2/b3 deletions 2. Which copies are lost is not claimed. Gene spans from Ensembl, checked on 2026-11-04.

A row on chrX or chrY is read against the sample's own count there: its controls sit on the same chromosome, so a male's one X or Y gives 1 copy, and `expected` says so. The sex is the declared one: step 35's second argument, and in Nextflow the samplesheet's sex column. Step 16's inferred sex is not used, so both give the same table. A chrY row on a female is `not applicable`; without a sex, chrX and chrY rows are `not assessed` (`sex_unknown`). Every row needs 8x of control depth (4x on chrY), flanks within 1.5x of each other, and an estimate within 0.3 of a whole number; a control with no read makes that row `not assessed` (`broken_control`) and the others still run. A failed bedcov pass writes every row as `not assessed` (`check_failed`).

The couple screen (step 39) reads the table with `couple_report.py screen --cn-intervals ID=<id>_cn_intervals.tsv`: each carrier allele is present, absent or not assessable for each partner, and a male partner's AZF deletions are listed the same way (present means the deletion is there). With `--cnvpytor ID=<id>_cnvs.vcf.gz` (step 18), a CNVpytor call that disagrees with the depth call is printed beside it. Step 39 and the Nextflow couple screen do not pass these two inputs yet.

To change a row, change the script (its hand rows or a pinned download) and run `scripts/ci/build-cn-intervals.py`; `--check` fails when the catalog is not what the downloads give.

## One owner per locus

A gene gets a row only when no pinned caller already counts it. The pinned Parascopy data (version 1.7) models 167 genes, among them AMY1C, C4A, CFH, CFHR4, FCGR3A, CHRNA7, APOBEC3A and GYPA. Parascopy owns those, so the engine never gets a row for them. A model covers its main region and the paralog copies it counts with it: C4A's model also covers C4B. `assets/parascopy_models_v1.7.txt` lists every one of those regions, and `tests/test_copy_number.py` fails when a row's targets or controls overlap one. Naming the model in `cross_check` does not lift that. Planted rows on AMY1C and on C4B show the check catching them.

One row is exempt on purpose, and the test lists it by name:

- **RHD**: RHCE's model holds most of RHD as one of its paralog copies. Step 41 runs by default and Parascopy is opt-in, so the RhD deletion stays a depth row. `cross_check` names RHCE as a second opinion, not a second owner.

The test also fails when an exemption no longer matches an overlap, so the list cannot go stale.

pypgx owns the genes it calls from the BAM, the `bam` rows of `assets/pgx_genes.tsv` (UGT2B17 among them), so none of them gets a row here either; `tests/test_pypgx_genes.sh` fails when one does.

## Adding a locus

1. Check that no pinned caller owns the gene: Parascopy (`assets/parascopy_models_v1.7.txt`), pypgx, Paraphase.
2. Choose single-copy controls near the gene, outside any common deletion or duplication there and outside segmental duplications.
3. Add the row to `bin/cn_loci.tsv`, above its marker line, and cite the method in `pmids`. An interval row (everything after the marker) is never added by hand: add it to `scripts/ci/build-cn-intervals.py` and rerun the script.
4. Run it from a step: `copy_number.py bed --locus NAME`, `samtools bedcov` with the row's passes, then `copy_number.py call --locus NAME --bedcov ... --json ...`.
5. Add made-up bedcov tables for the row to `tests/test_copy_number.py`, at every count it can call and at each guard.

## How it is tested

- `tests/test_copy_number.py` section 9 runs the converter on a small ClinGen BED whose chrX rows are spelled `chrx` (the rows load on chrX; the same rows left as `chrx` are refused; an unknown contig stops the build). Made-up depths call a 22q11.2 A-D deletion at 1 copy and a 17p12 duplication at 3, and 2 for both unaltered. AZFc at half its DAZ depth on a male is a deletion; every AZF row on a female is not applicable. No interval row's control overlaps another row's target. The interval block of step 35, cut out and run on made-up bedcov tables, writes one line per row, and a failed bedcov writes every row as not assessed. `tests/test_couple_report.py` section 11 checks the three states and the CNVpytor cross-check.
- The e2e case `cn-intervals-1-catalog` rebuilds the interval rows from the pinned downloads and fails on any difference with the committed catalog. `cn-intervals-2-hg002` runs the interval blocks of step 35 and of the PARASCOPY module on HG002's real reads at the AZF genes, their controls and 22q11.2: a male with one copy of AZFa and AZFb, no change at 22q11.2 A-B, no call on a row off chrY and chr22 (they have no read; the chr22 rows share A-B's control window, so they are left out of that check), and the AZF rows not applicable when the same reads are read as female.

- `tests/test_copy_number.py` checks every row and refuses broken ones. For the three gene-dosage rows it calls made-up depths at 0 to 4 copies, cuts the dosage block out of step 35 and runs it on the made-up bedcov tables in `tests/fixtures/cn/` (HP 3, LCE3BC 1, CCL3L1 5 copies), checks that a row shifted by 1 kb no longer matches its table, and that the PARASCOPY module runs the same list. It runs the owner check with planted AMY1C, C4B and DEFB4 rows. It feeds made-up bedcov tables at 0, 1, 2, 3 and 4 copies and a MAPQ 0 share over the ceiling. It checks that a control with no read stops the run with `BROKEN_CONTROL`. It also checks that the depth tables of `tests/test_hba_copy_number.py` and `tests/test_lpa_kiv2.py` give, through the engine, the same tables the code before the engine wrote.
- The fake-docker case `cn-engine-calls` runs steps 35, 41 and 46 with made-up depths. Each step calls `copy_number.py` once for its regions, step 35 and 46 once more for the call, and none calls the old per-gene scripts. A control exon with no read fails step 46 with `BROKEN_CONTROL`.
- `tests/test_copy_number.py` section 10 checks the CYP21A2 row: its six sites, the pileup parser, and the fixtures `tests/fixtures/cn/CYP21A2_{normal,deletion,conversion}` (2 copies and no finding; 1 copy with the module total 3; an In2G ratio of 0.5 and a possible conversion), with the conversion's depths read with every site even as the red control (no finding). The e2e case `cyp21a2-short-reads-1-hg002` runs the step's and the module's blocks on HG002's reads, unmodified and with a planted conversion and a planted deletion.
- The e2e case `cn-engine-1-parity` runs steps 41 and 46 and the HBA block of step 35 on HG002's real reads at the three loci, at full depth and at an eighth of it. The code before the engine is fetched by its commit, 6056452, so the comparison still means something after this change is on main. The LPA and HBA tables and the RhD rows of step 41's table and evidence file must be byte for byte what it wrote.
