# Copy Number from Depth (the locus catalog)

Several steps count copies of a gene from short-read depth: HBA1 and HBA2 in step 35, the RHD deletion in step 41 and the LPA KIV-2 repeat in step 46. They all use one engine, `bin/copy_number.py`, and one table, `bin/cn_loci.tsv`, with one row per locus. A new locus is a new row, not a new script.

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
- `class:LOW,HIGH`: low, typical or high from the 95% interval (LPA KIV-2, 30 and 45). With `interval` `cv`, that interval comes from the spread of the per-segment depths: relative standard error sqrt(cv_target^2 / n_target + cv_control^2 / n_control), applied on the log scale. It covers sampling noise only.

`copy_number.py call --json` writes the result as one JSON block per locus: status, reason, estimate, copies, interval, class, the depths and the reading.

## The columns

| Column | Meaning |
|---|---|
| `locus` | the row's name, as `--locus` takes it |
| `kind` | `gene`, `interval` (a stretch between segmental duplications, controls outside them) or `organelle` |
| `chrom` | GRCh38 contig of every segment |
| `targets`, `controls`, `details` | `NAME:start-end[:label]`, separated by `;`, in BED coordinates (0-based start, end exclusive). Details get a copy number each but take no part in the call (the RHD exons) |
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
| `cross_check` | another engine that also reports the locus (`parascopy:RHCE`, `paraphase:hba`) or `none` |
| `pmids` | the papers the method of the row comes from, or `.` |

## The loci

| Locus | Step | Targets against controls | Call |
|---|---|---|---|
| HBA | 35 | HBA2 and HBA1 gene bodies against two 40 kb flanks outside the `--SEA` deletion | whole number of alpha-globin genes ([step 35](35-paralogs.md)) |
| LPA_KIV2 | 46 | the 14 exons of the KIV-2 block against 15 single-copy LPA exons | repeat count with a 95% interval, low, typical or high ([step 46](46-lpa.md)) |
| RHD | 41 | the RHD gene against two flanks outside the common RHD deletion; the ten exons as details | 0, 1 or 2 copies ([step 41](41-bloodtype.md)) |

## One owner per locus

A gene gets a row only when no pinned caller already counts it. The pinned Parascopy data (version 1.7) models 167 genes, among them AMY1C, C4A, CFH, CFHR4, FCGR3A, CHRNA7, APOBEC3A and GYPA. Parascopy owns those, so the engine never gets a row for them. A model covers its main region and the paralog copies it counts with it: C4A's model also covers C4B. `assets/parascopy_models_v1.7.txt` lists every one of those regions, and `tests/test_copy_number.py` fails when a row's targets or controls overlap one. Naming the model in `cross_check` does not lift that. Planted rows on AMY1C and on C4B show the check catching them.

One row is exempt on purpose, and the test lists it by name:

- **RHD**: RHCE's model holds most of RHD as one of its paralog copies. Step 41 runs by default and Parascopy is opt-in, so the RhD deletion stays a depth row. `cross_check` names RHCE as a second opinion, not a second owner.

The test also fails when an exemption no longer matches an overlap, so the list cannot go stale.

pypgx owns the genes it calls from the BAM, the `bam` rows of `assets/pgx_genes.tsv` (UGT2B17 among them), so none of them gets a row here either; `tests/test_pypgx_genes.sh` fails when one does.

## Adding a locus

1. Check that no pinned caller owns the gene: Parascopy (`assets/parascopy_models_v1.7.txt`), pypgx, Paraphase.
2. Choose single-copy controls near the gene, outside any common deletion or duplication there and outside segmental duplications.
3. Add the row to `bin/cn_loci.tsv` and cite the method in `pmids`.
4. Run it from a step: `copy_number.py bed --locus NAME`, `samtools bedcov` with the row's passes, then `copy_number.py call --locus NAME --bedcov ... --json ...`.
5. Add made-up bedcov tables for the row to `tests/test_copy_number.py`, at every count it can call and at each guard.

## How it is tested

- `tests/test_copy_number.py` checks every row and refuses broken ones. It runs the owner check with planted AMY1C and C4B rows. It feeds made-up bedcov tables at 0, 1, 2, 3 and 4 copies and a MAPQ 0 share over the ceiling. It checks that a control with no read stops the run with `BROKEN_CONTROL`. It also checks that the depth tables of `tests/test_hba_copy_number.py` and `tests/test_lpa_kiv2.py` give, through the engine, the same tables the code before the engine wrote.
- The fake-docker case `cn-engine-calls` runs steps 35, 41 and 46 with made-up depths. Each step calls `copy_number.py` once for its regions, step 35 and 46 once more for the call, and none calls the old per-gene scripts. A control exon with no read fails step 46 with `BROKEN_CONTROL`.
- The e2e case `cn-engine-1-parity` runs steps 41 and 46 and the HBA block of step 35 on HG002's real reads at the three loci, at full depth and at an eighth of it. The code before the engine is fetched by its commit, 6056452, so the comparison still means something after this change is on main. The LPA and HBA tables and the RhD rows of step 41's table and evidence file must be byte for byte what it wrote.
