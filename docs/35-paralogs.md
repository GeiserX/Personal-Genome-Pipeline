# Step 35: Paralog Genes, SMN1 and SMN2 Copy Number (Parascopy)

> **OPT-IN:** a default run leaves this step out. Run it by name (`TOOLS=...,parascopy` for `run-all.sh`, `--tools ...,parascopy` for Nextflow) after installing its data once with `./scripts/setup.sh --parascopy-data ${GENOME_DIR}`.

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

A few minutes for the background depth over a 30x BAM, and under a minute for the SMN1/SMN2 locus.

## Interpreting Results

- **SMN1** is at about chr5:70.92-70.95 Mb and **SMN2** at about chr5:70.05-70.08 Mb (GRCh38). A row whose region lies in SMN1 lists SMN2's region under `homologous_regions`.
- **agCN** is the number of copies of the two genes together; most people have 4 (two of each), and anything from 2 to 6 is common.
- **psCN** splits agCN between the copies: the first number is the copy number of the row's region, then one per region in `homologous_regions`. `?` means that copy could not be told apart.
- **Quality** is a Phred score: 20 means 99% likely right. Parascopy recommends 20 as the threshold. A filter other than `PASS` means the value may be wrong even with a high quality.
- One SMN1 copy (psCN 1 for SMN1 with quality 20 or more) suggests SMA carrier status. Zero copies of SMN1 with copies of SMN2 is the pattern of SMA itself. A two-copy result does not rule out carrier status: some carriers, more in some populations, have two SMN1 copies on one chromosome and none on the other (a "2+0" carrier), which a copy number cannot show.

**This is a research estimate, not a diagnostic test.** SMA carrier screening uses validated assays (MLPA or qPCR). Confirm any result with a clinical laboratory.

## What Is and Is Not Assessed

- **Assessed:** the copy number of SMN1 and SMN2, together and, where possible, each.
- **Not assessed:** the "2+0" arrangement above; small variants inside SMN1 (a minority of SMA, and they need the reads placed on SMN1, which short reads cannot do reliably); GBA/GBAP1, PMS2, NCF1, STRC and the other paralog genes Parascopy has models for, and CYP21A2 and HBA1/HBA2, for which it has none. SMN1/SMN2 comes first; the others are not wired in yet.

## Notes

- The background depth comes from many windows across the genome. A BAM that holds only some regions (a targeted panel, or the e2e fixture) needs `PARASCOPY_DEPTH_BED` with windows where it has reads.
- Parascopy's population models are estimated from 1000 Genomes samples; a sample far from all five populations gets a less certain estimate.
- [Multi-sample](multi-sample.md) explains what a comparison of two samples can and cannot say here.

## Links

- [Parascopy](https://github.com/tprodanov/parascopy)
- [Prodanov and Bansal 2022](https://doi.org/10.1038/s41467-022-30930-3)
- [Parascopy precomputed data (Zenodo)](https://zenodo.org/records/15019940)
