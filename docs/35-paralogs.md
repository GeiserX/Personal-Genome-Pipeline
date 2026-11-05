# Step 35: Paralog Genes, SMN1 and SMN2 Copy Number (Parascopy), HBA1/HBA2 from Depth

> **OPT-IN:** a default run leaves this step out. Run it by name (`TOOLS=...,parascopy` for `run-all.sh`, `--tools ...,parascopy` for Nextflow) after installing its data once with `./scripts/setup.sh --parascopy-data ${GENOME_DIR}`.

> **PacBio HiFi reads?** [Step 35b](#step-35b-paralog-genes-from-pacbio-hifi-reads-paraphase) below runs Paraphase, which also covers HBA1/HBA2, CYP21A2, GBA1, PMS2, STRC, NCF1 and OTOA. It does not run on short reads.

## What This Does

Estimates how many copies of SMN1 and SMN2 the genome carries, from the short-read BAM, with [Parascopy](https://github.com/tprodanov/parascopy). For each stretch of the SMN1/SMN2 locus it reports the aggregate copy number of the two genes together (agCN) and, where the few sequence differences between them allow, the copy number of each (psCN), each with a quality.

The loci come from a catalog, `assets/paralog_loci.tsv`: one row per Parascopy model, and the step runs Parascopy once per row. The catalog has 37 rows: SMN1 and SMN2, and 36 more disease, carrier and blood-group genes with pseudogenes or paralogs, such as GBA1, PMS2, STRC, NCF1, C4A/C4B, the CFH cluster, RHCE and the amylase genes. See [The paralog catalog](#the-paralog-catalog) below.

It also counts the alpha-globin genes (HBA1 and HBA2 together) from read depth, for the alpha-thalassaemia deletions, because Parascopy has no model for them. See [HBA1 and HBA2 from depth](#hba1-and-hba2-from-depth) below.

## Why

Most spinal muscular atrophy (SMA) comes from the loss of both SMN1 copies, and most carriers have one SMN1 copy instead of two. SMN1 and SMN2 are near-identical copies on chr5, so short reads cannot be placed on one of them: their variant calls are unreliable, and a copy-number loss leaves no record in a VCF at all. No VCF-based step of this pipeline can see it. Parascopy models the reads of all copies together and the positions where the copies differ.

## Tool

- **Parascopy** 1.19.0 (Prodanov and Bansal, Nature Communications 2022), MIT licence

Not added, and why:

- **SMNCopyNumberCaller** (Illumina): archived since 2023, and under the PolyForm Strict licence (non-commercial use only), as Cyrius is.
- **Gauchian** (Illumina, GBA): the same licence, and no release since 2022.

## Docker Image

- `PARASCOPY_IMAGE` (optional: pulled the first time the step runs)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Data

`./scripts/setup.sh --parascopy-data ${GENOME_DIR}` installs, under `${GENOME_DIR}/reference/parascopy-1.7/` (`PARASCOPY_DATA_VERSION` in `versions.env`), Parascopy's precomputed data from [Zenodo record 15019940](https://zenodo.org/records/15019940) (CC-BY-4.0, about 50 MB, each archive checked against its md5):

- `homology_table/GRCh38.bed.gz`: the duplicated regions of GRCh38;
- `models_GRCh38_1KGP/<population>/`: model parameters estimated from 1000 Genomes samples of five populations (AFR, AMR, EAS, EUR, SAS).

Parascopy needs a BAM aligned to a reference without ALT contigs, which is the pipeline's default (see [realignment](realignment.md)).

## Command

```bash
./scripts/setup.sh --parascopy-data "$GENOME_DIR"   # once
./scripts/35-paralogs.sh your_name
```

| Variable | Default | Meaning |
|---|---|---|
| `PARASCOPY_POPULATION` | `EUR` | Whose model parameters to use: AFR, AMR, EAS, EUR or SAS. Pick the one closest to the sample's ancestry. |
| `PARALOG_LOCI` | `assets/paralog_loci.tsv` | The catalog of loci to run (bash step only; the Nextflow module reads the file in the repository). |
| `PARASCOPY_DEPTH_BED` | Parascopy's GRCh38 windows | Background windows (BED, every window the same size) for a BAM that covers only part of the genome. The depth is then not stratified by GC content (`--no-gc`): a few regions rarely span the GC range Parascopy's GC model needs, and it stops on them. The e2e test uses 100 bp windows over its chr20 slice. |

With Nextflow: `--tools ...,parascopy --parascopy_data ${GENOME_DIR}/reference/parascopy-1.7`, and `--parascopy_population`, `--parascopy_depth_bed` for the two variables.

## What the Script Does Internally

1. `parascopy depth`: background read depth of the sample, over Parascopy's GRCh38 windows (or `PARASCOPY_DEPTH_BED`).
2. For each row of the catalog: `samtools view -c` counts the reads over the row's region, then `parascopy cn-using models_GRCh38_1KGP/<population>/<model>` gives the copy number of the locus with those model parameters. A row whose model is not installed, whose region has no read (or whose reads `samtools` could not count), or on which Parascopy fails or writes no result file is not assessed: its table holds one row that says so, with the reason, and the other rows still run. Before any of this the step stops on a catalog row it cannot use: a malformed model, region or table name, a table name used twice, or the table names `hba` and `paraphase`, which other tables use.
3. Writes one row per region of each locus profile to `${SAMPLE}_<table>_copy_number.tsv` (`${SAMPLE}_smn_copy_number.tsv` for the SMN1 row).
4. `samtools bedcov` (in the Parascopy image) sums the depth over HBA2, HBA1 and two flanks, once over every read and once over reads of mapping quality 1 or more (`-Q 1`), and `bin/copy_number.py call --locus HBA` (the copy-number engine, on the HBA row of `bin/cn_loci.tsv`; see [Copy number from depth](cn-loci.md)) writes `${SAMPLE}_hba_copy_number.tsv`, laid out by `bin/hba_copy_number.py`. A failure here prints a warning and writes a table that says not assessed; the SMN result stands.

## Output

All output is written to `${GENOME_DIR}/${SAMPLE}/paralogs/`.

| File | Contents |
|---|---|
| `${SAMPLE}_smn_copy_number.tsv` | One row per region: `chrom`, `start`, `end`, `locus`, `agCN_filter`, `agCN`, `agCN_qual`, `psCN_filter`, `psCN`, `psCN_qual`, `homologous_regions`. A locus that was not assessed has one row over the catalog region (`chrom`, `start`, `end`), with `locus` (the model name), `agCN_filter` and `psCN_filter` (both `not assessed: <reason>`) and `.` in every other column |
| `${SAMPLE}_<table>_copy_number.tsv` | The same, for every other catalog row (`gba1`, `pms2_pms2cl`, `c4` and the rest of the table above) |
| `${SAMPLE}_hba_copy_number.tsv` | `metric`, `value`: `status` (`called` or `not assessed`), `alpha_globin_genes`, `estimate`, `hba2_depth`, `hba1_depth`, `flank_left_depth`, `flank_right_depth`, `mapq0_share` (the part of the depth over the two genes from reads of mapping quality 0), `message` |
| `${SAMPLE}_parascopy/<table>/` | Parascopy's own output for each catalog row (`smn/` for SMN1): `res.samples.bed.gz`, `res.paralog.bed.gz`, `psvs.vcf.gz` and the rest ([described upstream](https://github.com/tprodanov/parascopy/blob/main/docs/cn_output.md)) |

## Runtime

On the CI fixture with the one-row catalog the step took 8 seconds: one thread, the background depth over 4,000 windows of 100 bp on the chr20 slice, then the SMN1/SMN2 locus on the chr5 slice (case `pgx-4-parascopy`, e2e run 39544031036). The e2e case `paralog-rows-1` times the whole 37-row catalog on HG002 reads over every catalog locus on the fixture's chromosomes and fails past 30 minutes. A 30x genome has not been timed: there the background depth reads Parascopy's windows across the whole genome, which the fixture's slices do not have.

## Interpreting Results

- **SMN1** is at about chr5:70.92-70.95 Mb and **SMN2** at about chr5:70.05-70.08 Mb (GRCh38). A row whose region lies in SMN1 lists SMN2's region under `homologous_regions`.
- **agCN** is the number of copies of the two genes together; most people have 4 (two of each), and anything from 2 to 6 is common.
- **psCN** splits agCN between the copies: the first number is the copy number of the row's region, then one per region in `homologous_regions`. `?` means that copy could not be told apart.
- **Quality** is a Phred score: 20 means 99% likely right. Parascopy recommends 20 as the threshold. A filter other than `PASS` means the value may be wrong even with a high quality.
- One SMN1 copy (psCN 1 for SMN1 with quality 20 or more) suggests SMA carrier status. Zero copies of SMN1 with copies of SMN2 is the pattern of SMA itself. A two-copy result does not rule out carrier status: some carriers, more in some populations, have two SMN1 copies on one chromosome and none on the other (a "2+0" carrier), which a copy number cannot show.

**This is a research estimate, not a diagnostic test.** SMA carrier screening uses validated assays (MLPA or qPCR). Confirm any result with a clinical laboratory.

## What Is and Is Not Assessed

- **Assessed:** the copy number of SMN1 and SMN2, together and, where possible, each.
- **Assessed:** the number of alpha-globin genes, HBA1 and HBA2 together, from read depth (below).
- **Assessed:** the copies of the other 36 catalog loci, the gene and its paralogs together (agCN) (below).
- **Not assessed:** the "2+0" arrangement above; small variants inside SMN1 (a minority of SMA, and they need the reads placed on SMN1, which short reads cannot do reliably); which copy of a catalog locus other than SMN1 gained or lost one (its psCN is not read); the 130 Parascopy models with no catalog row; CYP21A2, for which Parascopy has no model; and which HBA gene or which chromosome lost a copy. The pinned data (version 1.7) holds 167 models per population, and none covers the alpha-globin locus on chr16 or CYP21A2. A PacBio HiFi BAM gets several of these genes, one by one, from step 35b.

## The Paralog Catalog

`assets/paralog_loci.tsv` lists the loci step 35 runs, one row per Parascopy model. Its columns:

| Column | Meaning |
|---|---|
| `model` | the model file in `models_GRCh38_1KGP/<population>/`, for example `SMN1.gz` |
| `genes` | the genes the model's region holds |
| `region` | the region of interest, GRCh38 `chrom:start-end` (1-based); a BAM with no read there gets the row not assessed, and Parascopy does not run on it |
| `reading` | what the copy number there means, per value |
| `pscn_reliable` | `yes` when the copy number of each copy (psCN) can be trusted there, `no` when only the copies together (agCN) can |
| `tag` | `carrier` (a recessive disease gene), `SF` (ACMG secondary finding), `disease` (a dominant or mosaic disease gene not on the SF list), `PGx` or `trait` |
| `cross_check` | another engine that calls the same locus: `paraphase` (step 35b, on HiFi reads), `cn-depth` (a row of `bin/cn_loci.tsv`), or `none` |
| `table` | the name in `${SAMPLE}_<table>_copy_number.tsv` and `${SAMPLE}_parascopy/<table>/` |
| `source` | where the reading comes from |

The SMN1 row (`SMN1.gz`, region chr5:70,951,946, SMN1 exon 7, table `smn`) is read for psCN, as above. The other 36 rows each cover their model's whole main region (from `assets/parascopy_models_v1.7.txt`), with `pscn_reliable` `no`: the reports read their agCN only.

| Model | Genes | Tag | Table |
|---|---|---|---|
| `GBA.gz` | GBA1, GBAP1 | carrier | `gba1` |
| `PMS2.gz` | PMS2, PMS2CL | SF | `pms2_pms2cl` |
| `STRC.gz` | STRC, STRCP1 | carrier | `strc` |
| `NCF1.gz` | NCF1, NCF1B, NCF1C | carrier | `ncf1` |
| `OTOA.gz` | OTOA, OTOAP1 | carrier | `otoa` |
| `NEB.gz` | NEB (the triplicate repeat) | carrier | `neb` |
| `CFH.gz` | CFH, CFHR1 | trait | `cfh` |
| `CFHR4.gz` | CFHR4, CFHR3 | trait | `cfhr4` |
| `C4A.gz` | C4A, C4B | trait | `c4` |
| `RHCE.gz` | RHCE (and the RHD sequence its model covers) | trait | `rhce` |
| `GYPA.gz` | GYPA, GYPB, GYPE | trait | `gypa` |
| `FCGR3A.gz` | FCGR3A, FCGR3B | trait | `fcgr3` |
| `FCGR2A.gz` | FCGR2A, FCGR2C | trait | `fcgr2a` |
| `FCGR1A.gz` | FCGR1A | trait | `fcgr1a` |
| `CD177.gz` | CD177 (HNA-2) | trait | `cd177` |
| `TPSB2.gz` | TPSB2, TPSAB1 | trait | `tpsb2` |
| `HSPA1A.gz` | HSPA1A, HSPA1B | trait | `hspa1` |
| `KLRC2.gz` | KLRC2 | trait | `klrc2` |
| `DEFA3.gz` | DEFA3, DEFA1, DEFA1B | trait | `defa1a3` |
| `CHRNA7.gz` | CHRNA7, CHRFAM7A | trait | `chrna7` |
| `HYDIN.gz` | HYDIN | carrier | `hydin` |
| `AMY1C.gz` | AMY1A, AMY1B, AMY1C | trait | `amy1` |
| `APOBEC3A.gz` | APOBEC3A, APOBEC3B | trait | `apobec3` |
| `ABCC6.gz` | ABCC6 | carrier | `abcc6` |
| `VWF.gz` | VWF | disease | `vwf` |
| `PROS1.gz` | PROS1 | disease | `pros1` |
| `SORD.gz` | SORD | carrier | `sord` |
| `GTF2I.gz` | GTF2I | disease | `gtf2i` |
| `GTF2H2C.gz` | GTF2H2C, GTF2H2 | trait | `gtf2h2` |
| `NOTCH2.gz` | NOTCH2 | disease | `notch2` |
| `KRT6C.gz` | KRT6C, KRT6A, KRT6B | disease | `krt6c` |
| `KRT86.gz` | KRT86 | disease | `krt86` |
| `PLG.gz` | PLG | carrier | `plg` |
| `PIK3CA.gz` | PIK3CA | disease | `pik3ca` |
| `CEL.gz` | CEL, CELP | disease | `cel` |
| `MUC5B.gz` | MUC5B | trait | `muc5b` |

Each row's `reading` says what its agCN means and its `source` the paper the reading comes from: Prodanov and Bansal 2022 (PMID 35680869) for the method on every row, Mandelker 2016 (PMID 27228465) for the clinical genes in segmental duplications, and per locus Sekar 2016 (C4, PMID 26814963), Cortese 2020 (SORD, PMID 32367058), Lyons 2016 (tryptase, PMID 27749843), Perry 2007 (amylase, PMID 17828263), Aitman 2006 (FCGR3B, PMID 16482158), Kidd 2007 (APOBEC3B deletion, PMID 17447845), Hughes 2006 and Zipfel 2007 (the CFHR3-CFHR1 deletion, PMIDs 16998489 and 17367211), Miyashita 2004 (KLRC2 deletion, PMID 14688071) and Fjeld 2015 (CEL-CELP, PMID 25774637).

One owner per locus. These loci are Parascopy's: no row of `bin/cn_loci.tsv` that the depth engine owns overlaps any region of their models. RHCE is the one exception, on purpose: its model covers most of RHD, so its agCN cross-checks the depth engine's RHD count (`cross_check` `cn-depth`), and RHD stays the depth engine's. The genes pypgx calls in BAM mode (UGT2B15, CYP2A6 through CYP2A7, SULT1A1 through SULT1A4) have no row. `assets/callability_genes.tsv` gives every catalog gene it lists the owner `parascopy` (its build script reads the catalog's `genes` column); a catalog gene its sources do not flag as hard to read, such as CFH or NEB, is not added there.

A new locus is a new row for a model already in the pinned data; it needs no code change. A rerun of the bash step removes the table of a row the catalog no longer has, so every catalog table in `paralogs/` comes from the last run. A row whose region has no read is not assessed and costs one `samtools view -c`, so a BAM of part of the genome runs Parascopy only where it has reads.

`assets/parascopy_models_v1.7.txt` lists every region of the 167 models of the pinned data (`models_GRCh38_1KGP_v1.7.tar.gz`, Zenodo record 15019940, archive sha256 `2e9be08b69936f13a49d52c7764f2414beac27a99b95e0c48b2e2ccafbb8c846`), the same for the five populations. Each gene there is the model file `<gene>.gz` of the `list.txt` the data installs. `tests/test_paralog_loci.py` refuses a catalog row whose model is not on that list: the gene GBA1, for example, has the model `GBA.gz`, not `GBA1.gz`.

## HBA1 and HBA2 from Depth

HBA1 and HBA2 are near-identical, so a short read from one may be placed on the other. But an aligner that reports every read puts each read from either on one of the two, with mapping quality 0 where it fits both equally. So the depth over both gene bodies together, against single-copy sequence on each side of the cluster, counts the genes: the HBA row of `bin/cn_loci.tsv`, run by `bin/copy_number.py` ([the locus catalog](cn-loci.md)), computes

    estimate = 2 x (depth over HBA2 + depth over HBA1) / mean flank depth

with every read counted whatever its mapping quality (duplicates, secondary and unmapped reads left out, as `samtools bedcov` does). The regions (GRCh38) are the MANE Select gene bodies of HBA2 (chr16:172,876-173,710) and HBA1 (chr16:176,680-177,522), and two 40 kb flanks, chr16:100,001-140,000 and chr16:190,001-230,000, outside the `--SEA` deletion and away from the HBZ and HBM genes.

| Genes | Reading |
|---|---|
| 4 | `aa/aa`, no deletion seen |
| 3 | one gene lost on one chromosome (`-a/aa`, as `-a3.7` or `-a4.2`): a silent carrier, marked unconfirmed (below) |
| 2 | alpha-thalassaemia trait: both genes lost on one chromosome (`--/aa`, as `--SEA`) or one on each (`-a/-a`); depth cannot tell these apart, and only the first puts a child at risk of Hb Bart's with a partner who carries the same |
| 1 | the pattern of HbH disease |
| 0 | the pattern of Hb Bart's hydrops fetalis |
| 5 or more | a triplication (`aaa/aa`) |

It says **not assessed**, with the reason, when either flank's mean depth is under 8x, when the two flanks differ by more than a factor of 1.5 (a deletion reaching a flank, or uneven coverage), when reads of mapping quality 0 make up under 2% of the depth over the two genes (unless the estimate rounds to 0, where no read is left to judge by), or when the estimate lies more than 0.3 from a whole number.

The mapping-quality check is there because the count needs the reads that fit both genes. The reference holds both genes, so such a read fits two places whatever the sample carries, and an aligner that reports it gives it mapping quality 0. minimap2 (step 02) and BWA-MEM do: on HG002's reads about a quarter of the depth over the two genes comes from such reads. novoalign, by default, leaves those reads unaligned. GIAB's novoalign BAM of HG002 has no read of mapping quality 0 over the genes, HBA2 sits at about half the depth of HBA1 and the flanks, and the count came out three (estimate 2.76) where HG002 has four. A BAM like that now gets not assessed, with the reason, instead of a low count.

Limits:

- It counts genes, not alleles: it cannot say which gene was lost, or whether two losses are on one chromosome or on two.
- The deletion must remove a gene body: a deletion that removes only regulatory sequence (the HS-40 enhancer) leaves the count at four.
- **One real genome checked.** HG002 has four genes (Paraphase on its HiFi reads), and this check counts four on its BWA-MEM reads (estimate 4.24 at 35x) and on the same reads realigned with minimap2 as step 02 does. No genome with a known deletion is in CI, so the accuracy on a real carrier has not been measured. Every count under four carries "unconfirmed: ... needs a clinical test before it means anything" in its message and in the report line. A count of four or more carries no such note.

## How It Is Tested

- `tests/test_paralog_loci.py` checks the model list (167 genes in `assets/parascopy_models_v1.7.txt`, each with one main region) and every catalog row (a model on the list, a GRCh38 region inside the model's main region, a known tag and cross-check, a table name used once, the SMN1 row), and that planted faults are refused: a model not on the list (`GBA1.gz`, and `CYP21A2.gz`, which the data has no model for), a bad region, a region outside its model, a repeated table name, a catalog without the SMN1 row. It checks one owner per locus: no catalog model overlaps a depth-engine row of `bin/cn_loci.tsv` but RHCE (without that exemption RHCE is refused, and a planted model over HBA is refused), and every catalog gene `assets/callability_genes.tsv` lists has the owner `parascopy` there (a planted `.` is refused). The catalog sets owners in that file but adds no rows to it: its rows come from GIAB's lists, gene body by gene body, and a model can cover only part of a gene. So 17 catalog genes are not there: seven pseudogenes (CELP, GBAP1, NCF1B, NCF1C, OTOAP1, PMS2CL, STRCP1), which it leaves out as not protein-coding, and ten genes to which none of its reasons applies (ABCC6, CFH, CFHR4, CHRNA7, DEFA3, KRT86, NEB, PIK3CA, PLG, VWF). This section reads them all. The test holds that list and fails when it changes. It also runs the step's own lines that write a not-assessed SMN table and checks that the report, the Paraphase side-by-side and the couple screen each read it as not assessed, with the reason, and that a catalog table without `agCN_filter` or `agCN_qual` reads as unreadable, not as not resolved.
- `tests/fake-docker/paralog-catalog-loop.sh` runs the step on a catalog of three rows against a fake Docker: three `parascopy cn-using` calls with the right model files, in order, one table each; another population's files with `PARASCOPY_POPULATION`; no call and a not-assessed table for a region without reads, a model that is not installed and a Parascopy run that fails, while the other rows run; a rerun with fewer rows removes the tables of the dropped rows and keeps the HBA and Paraphase tables; a bad row, or a region off chr1-22, X and Y or with its start past its end, stops the step before any container.
- The e2e case `paralog-catalog-1-gba` adds a GBA row (`GBA.gz`) to a copy of the catalog and runs the step on the fixture's slice plus GIAB's HG002 reads over GBA1 and GBAP1, sampled to about 30x: the GBA table must hold an integer agCN with its quality. On the fixture's slice alone, which has no read at GBA1, the same row must be not assessed. The case also checks that the committed models are the names in the installed `list.txt` of each population, that the SMN table is byte for byte what one `parascopy cn-using` of `SMN1.gz` writes, as the step did before the catalog, and that the loop of the Nextflow module, cut from the file, writes the same two tables in the same image.
- `tests/test_copy_number.py` checks the HBA row and that these depth tables give, through the engine, the same tables as the code before it. `tests/test_hba_copy_number.py` gives `bin/hba_copy_number.py` made-up depth tables: four genes at even depth, three with HBA2 at half depth, two, a triplication, an estimate between two whole numbers, low flanks, flanks that disagree, and no read of mapping quality 0 over the genes (the depths of the novoalign BAM, and a clean three), which must be not assessed.
- The e2e case `paralog-checks-and-doc-fixes-1-hba` writes synthetic 30x reads over chr16:95,001-235,000 for two samples, one with every second read over HBA2 left out (a planted one-gene loss), runs step 35 on each and checks four and three genes. A synthetic read whose bases occur more than once in that stretch gets mapping quality 0, as an aligner gives it. The HBA block of `modules/local/parascopy/main.nf`, cut from the file, then runs on the same BAMs in the same image and must write the same table. Last, the block runs on HG002's real Illumina reads at the cluster, three ways: GIAB's 60x novoalign BAM sampled to about 30x must be not assessed, naming the missing mapping quality 0 reads; Google's 35x NovaSeq BWA-MEM BAM must count four; and the same reads realigned with `minimap2 -a -x sr` to chr16:1-300,000 must count four. The truth is four, from Paraphase on HG002's HiFi reads. The run's notes record each estimate, the four depths and the share of mapping quality 0.
- The e2e case `paralog-rows-1` runs the whole catalog on the fixture's slice plus GIAB's HG002 reads over every catalog locus on the fixture's chromosomes, sampled to about 30x, timed (it fails past 30 minutes): every row's table must hold an integer agCN or say not assessed with the reason, and a header-only table fails that check. It then runs the step with a catalog that has no SMN1 row and checks that both reports still show the HBA count and the other loci. `tests/test_collect_summary.py` holds the same case on made-up tables, with the loci read in the catalog's order.
- The e2e case `hba-short-read-accuracy-1-nf-report` runs the Nextflow `HTML_REPORT` with and without an HBA table and checks the HBA line and the not-assessed list, then stub-runs every tool and checks that `main.nf` passes the table to each `HTML_REPORT` task that got the SMN table, against a copy of the pipeline with that join removed.

## In the Reports

Both reports (step 24 and the Nextflow `HTML_REPORT`) read `${SAMPLE}_smn_copy_number.tsv` into a section of its own, **SMN1/SMN2 Copy Number (Parascopy, short reads)**. It shows the copy number of SMN1 and of SMN2 at exon 7 (the row of the table, or its homologous region, that holds chr5:70,951,946 for SMN1 and chr5:70,076,526 for SMN2), with quality and filter; the two together (agCN) there; and a one-line reading: no SMN1 copy is the pattern of SMA, one copy that of a carrier, two copies is not a carrier by copy number with a 2+0 carrier not ruled out. A copy with a quality under 20 or a filter other than `PASS` is shown as not resolved and gets no reading, a table with no row over SMN1 exon 7 says that instead of a reading, and a table that says not assessed gives that reason instead. With step 35b's table too, the Paraphase section also sets the two tools' SMN1 side by side. When `${SAMPLE}_hba_copy_number.tsv` is there, the section adds one line, `HBA1 and HBA2 together (alpha-globin genes, from depth)`, with the count, its reading and the estimate, or the reason it was not assessed. Then one line per other catalog table, in the catalog's order (a table of a row the catalog does not name, from a `PARALOG_LOCI` of your own, last): the genes and their copies together (agCN). One number when every region of the locus profile has filter `PASS` and quality 20 or more and they agree; else each region with its value, quality and filter. The catalog's reading follows a call. A table that says not assessed gives the reason. In the HTML card these are a table: genes, copies together, reading. The HBA line and the other loci do not depend on the SMN table: a catalog without the SMN1 row still shows them, and SMN1 then says it has no table. The Nextflow `PARASCOPY` task writes the same tables to the same folder, and `HTML_REPORT` reads them the same way (`main.nf` passes it every catalog table, and both catalogs).

When neither step 35 nor step 35b left a result, both reports print `Paralog copy number: not checked (neither step 35 nor step 35b ran)`, naming a step that failed instead, so a missing check never reads as a normal result. Unless step 35's depth count made a call, or step 35b called both HBA1 and HBA2 (a run limited to other regions with `--paraphase_regions` has neither, and a `no call` row does not count), the "Not Assessed" list also names HBA1/HBA2 copy number.

## Notes

- The background depth comes from many windows across the genome. A BAM that holds only some regions (a targeted panel, or the e2e fixture) needs `PARASCOPY_DEPTH_BED` with windows where it has reads.
- Parascopy's population models are estimated from 1000 Genomes samples; a sample far from all five populations gets a less certain estimate.
- [Multi-sample](multi-sample.md) explains what a comparison of two samples can and cannot say here.

## Links

- [Parascopy](https://github.com/tprodanov/parascopy)
- [Prodanov and Bansal 2022](https://doi.org/10.1038/s41467-022-30930-3)
- [Parascopy precomputed data (Zenodo)](https://zenodo.org/records/15019940)

## Step 35b: Paralog Genes from PacBio HiFi Reads (Paraphase)

> **OPT-IN, PacBio HiFi only:** a default run leaves this step out. Run it by name (`TOOLS=...,paraphase` for `run-all.sh`, `--tools ...,paraphase` for Nextflow) on a BAM of PacBio HiFi reads. It has no data to install.

### What This Does

For a PacBio HiFi BAM, [Paraphase](https://github.com/PacificBiosciences/paraphase) gathers the reads of each family of near-identical genes, phases them into haplotypes (one per gene copy), and calls each gene's copy number and the variants each copy carries. It covers 165 GRCh38 regions (the top-level entries of `data/38/config.yaml` in Paraphase 4.1.0). The report shows the genes of `assets/paraphase_genes.tsv` one by one, each from the JSON field that catalog names:

| Gene | Why it matters | What the row says |
|---|---|---|
| SMN1, SMN2 | spinal muscular atrophy | copies of SMN1, of SMN2, and of SMN with exons 7-8 deleted |
| HBA1, HBA2 | alpha-thalassaemia | copies of each gene, the genotype (`aa/aa`, `-a/aa`, `--/aa`...), and the 3.7 kb and 4.2 kb deletions or duplications by name |
| CYP21A2 | congenital adrenal hyperplasia | copies of the RCCX module (each holds CYP21A2 or its pseudogene CYP21A1P), and Paraphase's annotation of each CYP21A2 allele |
| GBA1 | Gaucher disease, Parkinson's risk | copies of GBA1 apart from GBAP1 when Paraphase tells them apart, else of the family; the detail says which, and names a GBA1/GBAP1 fusion |
| PMS2 | Lynch syndrome | copies of PMS2, apart from PMS2CL |
| STRC | hearing loss | copies of STRC, apart from STRCP1 |
| NCF1 | chronic granulomatous disease | copies of NCF1, apart from NCF1B and NCF1C |
| OTOA | hearing loss | copies of OTOA apart from its pseudogene when Paraphase tells them apart, else of the family; the detail says which |
| F8 intron 22 | haemophilia A | copies of the intron 22 homology region, and the intron 22 inversion or h1/h2 deletion by name (`sv_called`) |
| OPN1LW, OPN1MW | red-green colour vision | copies of each (`opn1lw_cn`, `opn1mw_cn`) |
| IKBKG | incontinentia pigmenti | copies of IKBKG apart from IKBKGP1 (`gene_cn`), and the haplotypes with the 11.7 kb deletion |
| NEB triplicate repeat | nemaline myopathy | copies of the repeat, with the haplotypes Paraphase assigns to TRI1, TRI2 and TRI3 (`repeat_name`) |
| CYP11B1, CYP11B2 | congenital adrenal hyperplasia, glucocorticoid-remediable aldosteronism | copies of both, and a fusion by name (`fusions_called`) |
| TPSAB1, TPSB2 | hereditary alpha-tryptasemia (Lyons 2016, PMID 27749843) | copies of the tryptase genes together; extra copies need the alpha and beta copies told apart to mean alpha-tryptasemia |
| AMY1A, AMY1B, AMY1C | salivary amylase (Perry 2007, PMID 17828263) | copies of the three genes together |
| CFH cluster | atypical haemolytic uraemic syndrome, age-related macular degeneration | copies of the cluster (Paraphase's `CFHclust`, built from its CFH and CFHR3 regions), and its deletions and duplications by name, such as CFHR3-CFHR1 |

A row whose field of `region_specific_info` Paraphase did not write is **not assessed**, with the field it lacks. The two catalog fields outside it work differently: a `final_haplotypes` row with no gene1 or gene2 haplotypes gives the family's `total_cn`, as its catalog line says, and a row whose region has no `total_cn` says `no call`. HBA's `genotype` and RCCX's `annotated_alleles` are the two exceptions, because their copy number does not come from that field: without a genotype, HBA1 and HBA2 keep their counts from the haplotype names and say `no call`, as Paraphase could not resolve the alleles; without `annotated_alleles`, the RCCX row keeps the module's copies and lists no allele. A catalog row that reads either field on another region stops the table. Every other region gets one row with its copy number. A new region shown one by one is a new row of the catalog.

### Why Only PacBio HiFi

Paraphase is written for HiFi reads: a read of 10 kb or more spans the few differences between the copies, so the copies can be told apart. Short reads cannot do that, which is why step 35 uses Parascopy. Its licence (the PacBio Software License Agreement, from Paraphase 3.3.0 on) also allows it only on data from a PacBio instrument. So the step reads the BAM header first and stops unless every read group says `PL:PACBIO`. Step 02b with `PLATFORM=hifi` writes such a BAM, as PacBio's pbmm2 does.

On a short-read genome these genes stay where step 35 leaves them: SMN1 and SMN2 from Parascopy, HBA1 and HBA2 counted together from depth, the others not assessed.

### Tool

- **Paraphase** 4.1.0 (Chen et al., Nature Communications 2025), PacBio Software License Agreement: free to use, binary only, for data from a PacBio instrument

### Docker Image

- `PARAPHASE_IMAGE` (optional: pulled the first time the step runs)
- `PYTHON_IMAGE` for the read-group check and the two tables (`bin/paraphase_report.py`)

### Command

```bash
PLATFORM=hifi ./scripts/02b-alignment-longread.sh your_name   # writes aligned_longread/
./scripts/35b-paraphase.sh your_name
```

| Variable | Default | Meaning |
|---|---|---|
| `ALIGN_DIR` | `aligned_longread` | The folder of `<sample>_sorted.bam`. |
| `PARAPHASE_REGIONS` | all regions | Paraphase region names, comma separated, for example `smn1,hba,rccx,GBA,pms2,strc,ncf1,OTOA`. |

With Nextflow: `--tools ...,paraphase` on a samplesheet row whose BAM is PacBio HiFi, and `--paraphase_regions` for the variable. Its three processes are report-only: on a BAM that is not PacBio HiFi, `PARAPHASE_CHECK` fails, the rest of the run goes on, and `run-all.sh` marks step 35b failed.

### What the Script Does Internally

1. `bin/paraphase_report.py check-bam`: every read group of the BAM must say `PL:PACBIO`, or the step stops with the reason.
2. `paraphase -b <bam> -r <reference> -o ... -p <sample>`: the calls, as Paraphase's JSON, one VCF per region and a BAM of the phased reads.
3. `bin/paraphase_report.py table`: the copy-number table and the variant table.

### Output

All output is written to `${GENOME_DIR}/${SAMPLE}/paralogs/`.

| File | Contents |
|---|---|
| `${SAMPLE}_paraphase_copy_number.tsv` | One row per gene: `region`, `gene`, `copy_number`, `status` (`called`, `no call` or `not assessed`), `detail`, `total_cn` (the whole family), `region_depth`, `genome_depth` |
| `${SAMPLE}_paraphase_variants.tsv` | One row per variant a gene copy carries: `region`, `chrom`, `pos`, `ref`, `alt`, `filter`, `haplotypes` (the copies that carry it) and `support` (their ALT reads over their depth) |
| `${SAMPLE}_paraphase/` | Paraphase's own output: `${SAMPLE}.paraphase.json`, `${SAMPLE}_paraphase_vcfs/`, `${SAMPLE}.paraphase.bam` ([described upstream](https://github.com/PacificBiosciences/paraphase#interpreting-the-output)) |

### Interpreting Results

- **no call** means Paraphase saw one haplotype and the depth did not settle whether it is one copy or two identical ones.
- **HBA1 and HBA2:** most people have four alpha-globin genes (`aa/aa`). `-a` is an allele with one gene lost, `--` one with both lost. The 3.7 kb deletion (`-a3.7`, HBA2 and HBA1 fused into one gene) and the 4.2 kb deletion (`-a4.2`, HBA2 lost) are named, and so are the matching duplications. A `--` allele is a large deletion such as `--SEA`; its breakpoints lie outside the region Paraphase reads, so which large deletion it is cannot be told.
- **CYP21A2:** most people have four RCCX copies, two per chromosome. The allele annotation names a CYP21A2 deletion or duplication and the common pseudogene-derived variants on each allele.
- **SMN1 by both tools:** when step 35 (Parascopy) ran on the same BAM, the report shows both SMN1 values and says whether they agree. Parascopy's value counts only with filter `PASS` and quality 20 or more.

**This is a research estimate, not a diagnostic test.** Confirm a result, for example an HBA deletion or one SMN1 copy, with a clinical laboratory.

### What Is and Is Not Assessed

- **Assessed (PacBio HiFi only):** copy number and phased variants in Paraphase's 165 GRCh38 regions; the HBA 3.7 kb and 4.2 kb deletions and duplications by name.
- **Not assessed:** short-read and Oxford Nanopore BAMs (refused); which large HBA deletion removes both genes; the F8 intron 1 inversion, whose homology region is short (Paraphase leaves it to a standard HiFi SV caller, its docs/F8.md).

### Runtime

Not measured on a whole genome. Paraphase reads only the reads of its regions and of a few hundred depth windows, so the whole BAM is never read through.

### How It Is Tested

- `tests/test_paraphase_report.py` checks the region catalog: it holds the rows the test expects (a dropped row goes red); the old eight regions give the copy-number table byte for byte as the hard-coded table before the catalog wrote it (a changed JSON field goes red); a synthetic JSON of each new entry, shaped as Paraphase 4.1.0's source writes it, gives its rows; a region without the field its row reads is not assessed. It also checks that this page and `scripts/35b-paraphase.sh` state 165 regions, and the e2e case `paralog-rows-1` counts the top-level entries of `data/38/config.yaml` inside the pinned image and fails on another number.
- The image's smoke test runs Paraphase's own demo, the HG01175 SMN1/SMN2 HiFi reads shipped with the pinned release (their sha256 checked first), and checks the result its documentation publishes: two SMN1, one SMN2, one SMN with exons 7-8 deleted.
- The e2e case `paralogs-paraphase-1` reads GIAB's HG002 PacBio Revio BAM over HTTPS for the eight regions above and Paraphase's depth windows. It checks that a short-read BAM is refused, that each gene has its row with the copy number Paraphase calls on this slice (two each of SMN1, SMN2, PMS2, STRC, NCF1, HBA1, HBA2 and GBA1, four RCCX modules, no SMN exon 7-8 deletion), that GBA1 and OTOA count the gene alone rather than the family, and a planted 4.2 kb deletion: the 4,257 bp between the two homologous HBA boxes are cut out of the reads of one haplotype, and Paraphase must name the deletion and report one HBA2 copy fewer. A control relabels that result as swapped HBA1 and HBA2 regions in Paraphase's config would, and the HBA2 check must then fail. The same BAM then runs through Nextflow beside a short-read BAM, which fails the check while the run goes on. No per-sample Paraphase results for HG002 are published, so those copy numbers are Paraphase's own calls on this slice with the pinned image, not a published truth set: a change in the slice, the image or the parsing that moves one fails the case.

### Links

- [Paraphase](https://github.com/PacificBiosciences/paraphase) and its [region list](https://github.com/PacificBiosciences/paraphase/blob/main/docs/regions.md)
- [Chen et al. 2025](https://doi.org/10.1038/s41467-025-57505-2)
