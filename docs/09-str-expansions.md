# Step 9: Short Tandem Repeat (STR) Expansion Screening

## What This Does
Screens for pathogenic repeat expansions — a class of mutations invisible to both DeepVariant and Manta. In these diseases, a short DNA sequence (3-6 bases) gets repeated too many times.

## Why
STR expansions cause ~40 known neurological/neuromuscular diseases including Huntington's, Fragile X, Friedreich's ataxia, ALS/FTD, myotonic dystrophy, and multiple spinocerebellar ataxias.

## Tool
- **ExpansionHunter** v5.0.0 (Illumina): multithreaded, with a bundled GRCh38 variant catalog of 31 pathogenic loci
- **The catalog: STRchive's 82 disease loci.** STRchive ([Hiatt et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40140942/)) publishes every known disease repeat as an ExpansionHunter catalog, with the benign, intermediate and pathogenic sizes and the motifs of each locus in a second file. `setup.sh` downloads both files of `STRCHIVE_VERSION` (`versions.env`, checked against their sha256) and joins them with `bin/strchive_convert.py` into `reference/strchive-<version>/catalog.json`. A full setup installs it; `./scripts/setup.sh --strchive-data <genome_dir>` installs only it (about 0.5 MB). Without that file the step runs on the 31 bundled loci and says so in its log.

## Docker Image
- `EXPANSIONHUNTER_IMAGE`
- `PYTHON_IMAGE` (setup only) runs the converter

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

- Binary: `ExpansionHunter` (on PATH)
- Bundled GRCh38 catalog: `/usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json` (31 pathogenic loci)

## Which Catalog Runs

| `EH_CATALOG` | Catalog |
|---|---|
| unset | `reference/strchive-<version>/catalog.json` when `setup.sh` built it, else the 31 bundled loci |
| `bundled` | the 31 bundled loci |
| a file under `GENOME_DIR` | that file |

The step copies the catalog that ran to `expansion_hunter/<sample>_eh_catalog.json`; with the bundled catalog it writes none and removes an older copy. Step 9b grades with that copy when it carries Stranger's thresholds, as STRchive's does, so both steps read the same loci and thresholds; otherwise step 9b uses Stranger's bundled GRCh38 catalog ([step 9b](09b-stranger.md)). `run-all.sh` passes the same file to the pipeline as `--expansion_catalog`, and the pipeline's `EXPANSION_HUNTER` publishes it under the same name.

What the converter writes for each locus:

- **The name** is the gene in capitals (`HTT`, `FMR1`, `C9ORF72`, `RFC1`), the names the bundled catalog uses. Where one gene has two loci the name is STRchive's id (`EIEE1_ARX`, `PRTS_ARX`, `HFG_HOXA13-I` to `-III`). At a locus of several repeats (HTT's CAG and CCG, RFC1's AAAAG and AAGGG) the pathogenic repeat carries the locus name and the others `<name>_<motif>`.
- **Stranger's thresholds** (`NormalMax`, `PathologicMin`) come from STRchive's `benign_max` and `pathogenic_min`. Where STRchive gives no benign maximum (the FAME loci, SCA31, XDP and the GOLGA8A locus of aFTLD-U), `NormalMax` is the reference allele's unit count, so a longer allele is never called normal.
- **Seven loci get no thresholds**, so step 9b gives them no status, and the step log lists them as "Not graded": MUC1, CEL, TYMS and PLIN4, whose reference repeat is longer than 300 bases, so no short read spans it; VWA1 and POLG, where STRchive gives no pathogenic size; and MIR7-2, where disease comes from a contraction. ExpansionHunter still writes a record for each.
- **The motifs** of the eleven loci where the motif decides, not the size: RFC1, the FAME loci (SAMD12, STARD7, YEATS2, MARCHF6, RAPGEF2, TNRC6A, RAI1), SCA31 (BEAN1), SCA37 (DAB1) and aFTLD-U (GOLGA8A). Step 9b counts them in the reads ([step 9b](09b-stranger.md)).

## Key Disease Thresholds

The thresholds step 9b grades with on the STRchive v2.28.0 catalog, in repeat units: Stranger calls a count up to the normal maximum `normal`, one in the pre-mutation band `pre_mutation`, and one from the full-mutation minimum `full_mutation`. `tests/test_str_ranges.py` checks every row against the catalog. With `EH_CATALOG=bundled`, step 9b uses Stranger's bundled catalog, whose thresholds differ at some loci.

| Disease | Gene | Repeat Unit | Normal | Pre-mutation | Full mutation |
|---|---|---|---|---|---|
| Huntington's | HTT | CAG | ≤26 | 27-35 | ≥36 |
| Fragile X | FMR1 | CGG | ≤44 | 45-200 | ≥201 |
| Friedreich's Ataxia | FXN | GAA | ≤33 | 34-55 | ≥56 |
| ALS/FTD | C9ORF72 | GGCCCC | ≤23 | 24-30 | ≥31 |
| Myotonic Dystrophy 1 | DMPK | CTG | ≤34 | 35-49 | ≥50 |
| SCA1 | ATXN1 | CAG | ≤35 | 36-38 | ≥39 |
| SCA2 | ATXN2 | CAG | ≤28 | 29-34 | ≥35 |

Stranger's pre-mutation band is everything between the two thresholds; it is wider than a clinical intermediate range at some loci. At FMR1 it covers both the gray zone and the premutation below. `setup.sh` also writes `reference/strchive-<version>/str_ranges.tsv`, one row per locus with these two thresholds and STRchive's own intermediate range.

## FMR1 Clinical Zones

FMR1 (Fragile X) has four distinct clinical zones — the intermediate zone (45-54 repeats) is often omitted but clinically relevant:

| Zone | Repeats | Clinical Significance |
|---|---|---|
| Normal | <45 | No risk |
| Intermediate (gray zone) | 45-54 | Not affected, but repeats may expand in offspring. Genetic counseling recommended for carriers. |
| Premutation | 55-200 | Risk of FXTAS (tremor/ataxia, males >50), FXPOI (premature ovarian insufficiency). Offspring at risk of full expansion. |
| Full mutation | >200 | Fragile X syndrome (intellectual disability, behavioral features). Penetrance varies by sex and methylation. |

## Command

```bash
./scripts/09-expansion-hunter.sh your_name male
# or: ./scripts/09-expansion-hunter.sh your_name female
```

The second argument (`male`/`female`) is **required** — it affects X-linked loci (FMR1, AR): males have one allele, females have two.

## Notes
- Uses ExpansionHunter **v5.0.0** (`EXPANSIONHUNTER_IMAGE`)
- v5 CLI: `--reads`, `--reference`, `--variant-catalog` (JSON file), `--output-prefix` (auto-generates .vcf, .json)
- **Why STRchive, not a larger catalog.** STRchive lists the loci with a published disease link, each with the sources of its thresholds. Larger catalogs (tens of thousands of loci) are for research, have no thresholds, and would make step 9b's grades meaningless. Illumina has published no newer catalog than the bundled one (2021).
- **What the STRchive catalog cannot do.** It does not make long repeats readable: the seven loci above stay ungraded, and an expansion longer than a read is still sized only roughly. The motif count of step 9b reads the motif, never the size.
- Short-read WGS can reliably detect expansions up to ~150 repeats; very large expansions (>1000) are less accurate
- **How well each allele was sized.** The VCF's `SO` field says which reads support each allele: `SPANNING` (reads cross the whole repeat, so it is sized), `FLANKING` (reads reach into the repeat from one side only, so the size is not determined) or `INREPEAT` (reads made only of repeat, so the allele is longer than a read). `ADSP` counts the spanning reads and `REPCI` gives the interval. The report uses them: a `FLANKING` or `INREPEAT` allele, one with fewer than 2 spanning reads, or one whose interval reaches the pathogenic minimum is listed as not sized and never as normal ([how the report reads them](interpreting-results.md#str-expansions-step-9))
