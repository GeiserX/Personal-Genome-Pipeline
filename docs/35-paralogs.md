# Step 35: Paralog Genes, SMN1 and SMN2 Copy Number (Parascopy), HBA1/HBA2 from Depth

> **OPT-IN:** a default run leaves this step out. Run it by name (`TOOLS=...,parascopy` for `run-all.sh`, `--tools ...,parascopy` for Nextflow) after installing its data once with `./scripts/setup.sh --parascopy-data ${GENOME_DIR}`.

> **PacBio HiFi reads?** [Step 35b](#step-35b-paralog-genes-from-pacbio-hifi-reads-paraphase) below runs Paraphase, which also covers HBA1/HBA2, CYP21A2, GBA1, PMS2, STRC, NCF1 and OTOA. It does not run on short reads.

## What This Does

Estimates how many copies of SMN1 and SMN2 the genome carries, from the short-read BAM, with [Parascopy](https://github.com/tprodanov/parascopy). For each stretch of the SMN1/SMN2 locus it reports the aggregate copy number of the two genes together (agCN) and, where the few sequence differences between them allow, the copy number of each (psCN), each with a quality.

The loci come from a catalog, `assets/paralog_loci.tsv`: one row per Parascopy model, and the step runs Parascopy once per row. Today the catalog has one row, SMN1 and SMN2. See [The paralog catalog](#the-paralog-catalog) below.

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
2. For each row of the catalog: `samtools view -c` counts the reads over the row's region, then `parascopy cn-using models_GRCh38_1KGP/<population>/<model>` gives the copy number of the locus with those model parameters. A row whose model is not installed, whose region has no read (or whose reads `samtools` could not count), or on which Parascopy fails is not assessed: its table holds one row that says so, with the reason, and the other rows still run. Before any of this the step stops on a catalog row it cannot use: a malformed model, region or table name, a table name used twice, or the table names `hba` and `paraphase`, which other tables use.
3. Writes one row per region of each locus profile to `${SAMPLE}_<table>_copy_number.tsv` (`${SAMPLE}_smn_copy_number.tsv` for the SMN1 row).
4. `samtools bedcov` (in the Parascopy image) sums the depth over HBA2, HBA1 and two flanks, once over every read and once over reads of mapping quality 1 or more (`-Q 1`), and `bin/hba_copy_number.py check` writes `${SAMPLE}_hba_copy_number.tsv`. A failure here prints a warning and writes a table that says not assessed; the SMN result stands.

## Output

All output is written to `${GENOME_DIR}/${SAMPLE}/paralogs/`.

| File | Contents |
|---|---|
| `${SAMPLE}_smn_copy_number.tsv` | One row per region: `chrom`, `start`, `end`, `locus`, `agCN_filter`, `agCN`, `agCN_qual`, `psCN_filter`, `psCN`, `psCN_qual`, `homologous_regions`. A locus that was not assessed has one row over the catalog region (`chrom`, `start`, `end`), with `locus` (the model name), `agCN_filter` and `psCN_filter` (both `not assessed: <reason>`) and `.` in every other column |
| `${SAMPLE}_<table>_copy_number.tsv` | The same, for every other catalog row (none in the catalog today) |
| `${SAMPLE}_hba_copy_number.tsv` | `metric`, `value`: `status` (`called` or `not assessed`), `alpha_globin_genes`, `estimate`, `hba2_depth`, `hba1_depth`, `flank_left_depth`, `flank_right_depth`, `mapq0_share` (the part of the depth over the two genes from reads of mapping quality 0), `message` |
| `${SAMPLE}_parascopy/<table>/` | Parascopy's own output for each catalog row (`smn/` for SMN1): `res.samples.bed.gz`, `res.paralog.bed.gz`, `psvs.vcf.gz` and the rest ([described upstream](https://github.com/tprodanov/parascopy/blob/main/docs/cn_output.md)) |

## Runtime

On the CI fixture the step takes 8 seconds: one thread, the background depth over 4,000 windows of 100 bp on the chr20 slice, then the SMN1/SMN2 locus on the chr5 slice (case `pgx-4-parascopy`, e2e run 39544031036). A 30x genome has not been timed: there the background depth reads Parascopy's windows across the whole genome, which the fixture's slices do not have.

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
- **Not assessed:** the "2+0" arrangement above; small variants inside SMN1 (a minority of SMA, and they need the reads placed on SMN1, which short reads cannot do reliably); GBA/GBAP1, PMS2, NCF1, STRC and the other paralog genes Parascopy has models for, which have no catalog row yet; CYP21A2, for which Parascopy has no model; and which HBA gene or which chromosome lost a copy. The pinned data (version 1.7) holds 167 models per population, and none covers the alpha-globin locus on chr16 or CYP21A2. A PacBio HiFi BAM gets the genes above from step 35b.

## The Paralog Catalog

`assets/paralog_loci.tsv` lists the loci step 35 runs, one row per Parascopy model. Its columns:

| Column | Meaning |
|---|---|
| `model` | the model file in `models_GRCh38_1KGP/<population>/`, for example `SMN1.gz` |
| `genes` | the genes the model's region holds |
| `region` | the region of interest, GRCh38 `chrom:start-end` (1-based); a BAM with no read there gets the row not assessed, and Parascopy does not run on it |
| `reading` | what the copy number there means, per value |
| `pscn_reliable` | `yes` when the copy number of each copy (psCN) can be trusted there, `no` when only the copies together (agCN) can |
| `tag` | `carrier`, `SF` (secondary finding), `PGx` or `trait` |
| `cross_check` | another engine that calls the same locus (`paraphase` on HiFi reads), or `none` |
| `table` | the name in `${SAMPLE}_<table>_copy_number.tsv` and `${SAMPLE}_parascopy/<table>/` |
| `source` | where the reading comes from |

The SMN1 row (`SMN1.gz`, region chr5:70,951,946, SMN1 exon 7, table `smn`) is the only row today. A new locus is a new row for a model already in the pinned data; it needs no code change. The reports read the SMN table only, so a new row's table is written and published but not yet shown in either report.

`assets/parascopy_models_v1.7.txt` is the model list of the pinned data: the `list.txt` of `models_GRCh38_1KGP_v1.7.tar.gz` (Zenodo record 15019940, archive sha256 `2e9be08b69936f13a49d52c7764f2414beac27a99b95e0c48b2e2ccafbb8c846`), 167 models, the same for the five populations. `tests/test_paralog_loci.py` refuses a catalog row whose model is not on that list: the gene GBA1, for example, has the model `GBA.gz`, not `GBA1.gz`.

## HBA1 and HBA2 from Depth

HBA1 and HBA2 are near-identical, so a short read from one may be placed on the other. But an aligner that reports every read puts each read from either on one of the two, with mapping quality 0 where it fits both equally. So the depth over both gene bodies together, against single-copy sequence on each side of the cluster, counts the genes: `bin/hba_copy_number.py` computes

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

- `tests/test_paralog_loci.py` checks the model list (167 names, its sha256) and every catalog row (a model on the list, a GRCh38 region, a table name used once, the SMN1 row), and that planted faults are refused: a model not on the list, a bad region, a repeated table name, a catalog without the SMN1 row. It also runs the step's own lines that write a not-assessed SMN table and checks that the report, the Paraphase side-by-side and the couple screen each read it as not assessed, with the reason.
- `tests/fake-docker/paralog-catalog-loop.sh` runs the step on a catalog of three rows against a fake Docker: three `parascopy cn-using` calls with the right model files, in order, one table each; another population's files with `PARASCOPY_POPULATION`; no call and a not-assessed table for a region without reads, a model that is not installed and a Parascopy run that fails, while the other rows run; a bad row stops the step before any container.
- The e2e case `paralog-catalog-1-gba` adds a GBA row (`GBA.gz`) to a copy of the catalog and runs the step on the fixture's slice plus GIAB's HG002 reads over GBA1 and GBAP1, sampled to about 30x: the GBA table must hold an integer agCN with its quality. On the fixture's slice alone, which has no read at GBA1, the same row must be not assessed. The case also checks that the committed model list is the installed `list.txt` of each population, that the SMN table is byte for byte what one `parascopy cn-using` of `SMN1.gz` writes, as the step did before the catalog, and that the loop of the Nextflow module, cut from the file, writes the same two tables in the same image.
- `tests/test_hba_copy_number.py` gives `bin/hba_copy_number.py` made-up depth tables: four genes at even depth, three with HBA2 at half depth, two, a triplication, an estimate between two whole numbers, low flanks, flanks that disagree, and no read of mapping quality 0 over the genes (the depths of the novoalign BAM, and a clean three), which must be not assessed.
- The e2e case `paralog-checks-and-doc-fixes-1-hba` writes synthetic 30x reads over chr16:95,001-235,000 for two samples, one with every second read over HBA2 left out (a planted one-gene loss), runs step 35 on each and checks four and three genes. A synthetic read whose bases occur more than once in that stretch gets mapping quality 0, as an aligner gives it. The HBA block of `modules/local/parascopy/main.nf`, cut from the file, then runs on the same BAMs in the same image and must write the same table. Last, the block runs on HG002's real Illumina reads at the cluster, three ways: GIAB's 60x novoalign BAM sampled to about 30x must be not assessed, naming the missing mapping quality 0 reads; Google's 35x NovaSeq BWA-MEM BAM must count four; and the same reads realigned with `minimap2 -a -x sr` to chr16:1-300,000 must count four. The truth is four, from Paraphase on HG002's HiFi reads. The run's notes record each estimate, the four depths and the share of mapping quality 0.
- The e2e case `hba-short-read-accuracy-1-nf-report` runs the Nextflow `HTML_REPORT` with and without an HBA table and checks the HBA line and the not-assessed list, then stub-runs every tool and checks that `main.nf` passes the table to each `HTML_REPORT` task that got the SMN table, against a copy of the pipeline with that join removed.

## In the Reports

Both reports (step 24 and the Nextflow `HTML_REPORT`) read `${SAMPLE}_smn_copy_number.tsv` into a section of its own, **SMN1/SMN2 Copy Number (Parascopy, short reads)**. It shows the copy number of SMN1 and of SMN2 at exon 7 (the row of the table, or its homologous region, that holds chr5:70,951,946 for SMN1 and chr5:70,076,526 for SMN2), with quality and filter; the two together (agCN) there; and a one-line reading: no SMN1 copy is the pattern of SMA, one copy that of a carrier, two copies is not a carrier by copy number with a 2+0 carrier not ruled out. A copy with a quality under 20 or a filter other than `PASS` is shown as not resolved and gets no reading, a table with no row over SMN1 exon 7 says that instead of a reading, and a table that says not assessed gives that reason instead. With step 35b's table too, the Paraphase section also sets the two tools' SMN1 side by side. When `${SAMPLE}_hba_copy_number.tsv` is there, the section adds one line, `HBA1 and HBA2 together (alpha-globin genes, from depth)`, with the count, its reading and the estimate, or the reason it was not assessed. The Nextflow `PARASCOPY` task writes the same table to the same folder, and `HTML_REPORT` reads it the same way.

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

For a PacBio HiFi BAM, [Paraphase](https://github.com/PacificBiosciences/paraphase) gathers the reads of each family of near-identical genes, phases them into haplotypes (one per gene copy), and calls each gene's copy number and the variants each copy carries. It covers 162 GRCh38 regions. The report shows these genes one by one:

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

Every other region gets one row with its copy number.

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
| `${SAMPLE}_paraphase_copy_number.tsv` | One row per gene: `region`, `gene`, `copy_number`, `status` (`called` or `no call`), `detail`, `total_cn` (the whole family), `region_depth`, `genome_depth` |
| `${SAMPLE}_paraphase_variants.tsv` | One row per variant a gene copy carries: `region`, `chrom`, `pos`, `ref`, `alt`, `filter`, `haplotypes` (the copies that carry it) and `support` (their ALT reads over their depth) |
| `${SAMPLE}_paraphase/` | Paraphase's own output: `${SAMPLE}.paraphase.json`, `${SAMPLE}_paraphase_vcfs/`, `${SAMPLE}.paraphase.bam` ([described upstream](https://github.com/PacificBiosciences/paraphase#interpreting-the-output)) |

### Interpreting Results

- **no call** means Paraphase saw one haplotype and the depth did not settle whether it is one copy or two identical ones.
- **HBA1 and HBA2:** most people have four alpha-globin genes (`aa/aa`). `-a` is an allele with one gene lost, `--` one with both lost. The 3.7 kb deletion (`-a3.7`, HBA2 and HBA1 fused into one gene) and the 4.2 kb deletion (`-a4.2`, HBA2 lost) are named, and so are the matching duplications. A `--` allele is a large deletion such as `--SEA`; its breakpoints lie outside the region Paraphase reads, so which large deletion it is cannot be told.
- **CYP21A2:** most people have four RCCX copies, two per chromosome. The allele annotation names a CYP21A2 deletion or duplication and the common pseudogene-derived variants on each allele.
- **SMN1 by both tools:** when step 35 (Parascopy) ran on the same BAM, the report shows both SMN1 values and says whether they agree. Parascopy's value counts only with filter `PASS` and quality 20 or more.

**This is a research estimate, not a diagnostic test.** Confirm a result, for example an HBA deletion or one SMN1 copy, with a clinical laboratory.

### What Is and Is Not Assessed

- **Assessed (PacBio HiFi only):** copy number and phased variants in Paraphase's 162 GRCh38 regions; the HBA 3.7 kb and 4.2 kb deletions and duplications by name.
- **Not assessed:** short-read and Oxford Nanopore BAMs (refused); which large HBA deletion removes both genes.

### Runtime

Not measured on a whole genome. Paraphase reads only the reads of its regions and of a few hundred depth windows, so the whole BAM is never read through.

### How It Is Tested

- The image's smoke test runs Paraphase's own demo, the HG01175 SMN1/SMN2 HiFi reads shipped with the pinned release (their sha256 checked first), and checks the result its documentation publishes: two SMN1, one SMN2, one SMN with exons 7-8 deleted.
- The e2e case `paralogs-paraphase-1` reads GIAB's HG002 PacBio Revio BAM over HTTPS for the eight regions above and Paraphase's depth windows. It checks that a short-read BAM is refused, that each gene has its row with the copy number Paraphase calls on this slice (two each of SMN1, SMN2, PMS2, STRC, NCF1, HBA1, HBA2 and GBA1, four RCCX modules, no SMN exon 7-8 deletion), that GBA1 and OTOA count the gene alone rather than the family, and a planted 4.2 kb deletion: the 4,257 bp between the two homologous HBA boxes are cut out of the reads of one haplotype, and Paraphase must name the deletion and report one HBA2 copy fewer. A control relabels that result as swapped HBA1 and HBA2 regions in Paraphase's config would, and the HBA2 check must then fail. The same BAM then runs through Nextflow beside a short-read BAM, which fails the check while the run goes on. No per-sample Paraphase results for HG002 are published, so those copy numbers are Paraphase's own calls on this slice with the pinned image, not a published truth set: a change in the slice, the image or the parsing that moves one fails the case.

### Links

- [Paraphase](https://github.com/PacificBiosciences/paraphase) and its [region list](https://github.com/PacificBiosciences/paraphase/blob/main/docs/regions.md)
- [Chen et al. 2025](https://doi.org/10.1038/s41467-025-57505-2)
