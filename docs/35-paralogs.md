# Step 35: Paralog Genes, SMN1 and SMN2 Copy Number (Parascopy)

> **OPT-IN:** a default run leaves this step out. Run it by name (`TOOLS=...,parascopy` for `run-all.sh`, `--tools ...,parascopy` for Nextflow) after installing its data once with `./scripts/setup.sh --parascopy-data ${GENOME_DIR}`.

> **PacBio HiFi reads?** [Step 35b](#step-35b-paralog-genes-from-pacbio-hifi-reads-paraphase) below runs Paraphase, which also covers HBA1/HBA2, CYP21A2, GBA1, PMS2, STRC, NCF1 and OTOA. It does not run on short reads.

## What This Does

Estimates how many copies of SMN1 and SMN2 the genome carries, from the short-read BAM, with [Parascopy](https://github.com/tprodanov/parascopy). For each stretch of the SMN1/SMN2 locus it reports the aggregate copy number of the two genes together (agCN) and, where the few sequence differences between them allow, the copy number of each (psCN), each with a quality.

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
| `PARASCOPY_DEPTH_BED` | Parascopy's GRCh38 windows | Background windows (BED, every window the same size) for a BAM that covers only part of the genome. The depth is then not stratified by GC content (`--no-gc`): a few regions rarely span the GC range Parascopy's GC model needs, and it stops on them. The e2e test uses 100 bp windows over its chr20 slice. |

With Nextflow: `--tools ...,parascopy --parascopy_data ${GENOME_DIR}/reference/parascopy-1.7`, and `--parascopy_population`, `--parascopy_depth_bed` for the two variables.

## What the Script Does Internally

1. `parascopy depth`: background read depth of the sample, over Parascopy's GRCh38 windows (or `PARASCOPY_DEPTH_BED`).
2. `parascopy cn-using models_GRCh38_1KGP/<population>/SMN1.gz`: the copy number of the SMN1/SMN2 locus with those model parameters.
3. Writes one row per region of the copy-number profile to `${SAMPLE}_smn_copy_number.tsv`.

## Output

All output is written to `${GENOME_DIR}/${SAMPLE}/paralogs/`.

| File | Contents |
|---|---|
| `${SAMPLE}_smn_copy_number.tsv` | One row per region: `chrom`, `start`, `end`, `locus`, `agCN_filter`, `agCN`, `agCN_qual`, `psCN_filter`, `psCN`, `psCN_qual`, `homologous_regions` |
| `${SAMPLE}_parascopy/` | Parascopy's own output: `res.samples.bed.gz`, `res.paralog.bed.gz`, `psvs.vcf.gz` and the rest ([described upstream](https://github.com/tprodanov/parascopy/blob/main/docs/cn_output.md)) |

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
- **Not assessed:** the "2+0" arrangement above; small variants inside SMN1 (a minority of SMA, and they need the reads placed on SMN1, which short reads cannot do reliably); GBA/GBAP1, PMS2, NCF1, STRC and the other paralog genes Parascopy has models for, which are not wired in yet for short reads; and CYP21A2 and HBA1/HBA2, for which Parascopy has no model. The pinned data (version 1.7) holds 167 models per population, and none covers the alpha-globin locus on chr16 or CYP21A2. HBA1/HBA2 copy number on short reads, the 3.7 kb and 4.2 kb alpha-thalassaemia deletions first, needs another method, such as a depth check over the two deletion spans against their flanks; none is built yet. A PacBio HiFi BAM gets the genes above from step 35b.

## In the Reports

Both reports (step 24 and the Nextflow `HTML_REPORT`) read `${SAMPLE}_smn_copy_number.tsv` into a section of its own, **SMN1/SMN2 Copy Number (Parascopy, short reads)**. It shows the copy number of SMN1 and of SMN2 at exon 7 (the row of the table, or its homologous region, that holds chr5:70,951,946 for SMN1 and chr5:70,076,526 for SMN2), with quality and filter; the two together (agCN) there; and a one-line reading: no SMN1 copy is the pattern of SMA, one copy that of a carrier, two copies is not a carrier by copy number with a 2+0 carrier not ruled out. A copy with a quality under 20 or a filter other than `PASS` is shown as not resolved and gets no reading. With step 35b's table too, the Paraphase section also sets the two tools' SMN1 side by side.

When neither step 35 nor step 35b left a result, both reports print `Paralog copy number: not checked (neither step 35 nor step 35b ran)`, naming a step that failed instead, so a missing check never reads as a normal result. Unless step 35b's table has both HBA1 and HBA2 rows (a run limited to other regions with `--paraphase_regions` has neither), the "Not Assessed" list also names HBA1/HBA2 copy number from short reads.

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

On a short-read genome these genes stay where step 35 leaves them: SMN1 and SMN2 from Parascopy, the others not assessed.

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
