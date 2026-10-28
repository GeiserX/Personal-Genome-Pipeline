# Step 9b: STR Clinical Annotation (Stranger)

## What This Does

Annotates the ExpansionHunter VCF (step 9) with clinical pathogenicity status for each repeat locus. Raw repeat counts become labelled calls: **normal**, **pre_mutation**, or **full_mutation**, with disease name, OMIM number, inheritance mode, and the specific repeat-size thresholds used for classification.

This turns step 9 output from a table of numbers into actionable clinical calls.

## Why

ExpansionHunter reports the number of repeats at each locus but applies no pathogenicity judgement. Stranger adds that judgement using the ClinGen/OMIM short tandem repeat database, so you can immediately see whether a repeat count falls in the normal, pre-mutation, or full-mutation range without manually cross-referencing disease thresholds.

## Tool

- **Stranger** v0.10.2 (Clinical Genomics Stockholm): annotates STR VCFs with pathogenicity labels from a repeat catalog

## Docker Image

- `STRANGER_IMAGE`
- `PYPGX_IMAGE` for the RFC1 motif count (it has pysam)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

- Binary: `stranger` (on PATH)
- Bundled repeat catalog: Stranger's GRCh38 catalog (`variant_catalog_grch38.json`, installed inside the container). Stranger's own default is its GRCh37 catalog; the script and the Nextflow module pass the GRCh38 one explicitly, because ExpansionHunter (step 9) calls the GRCh38 loci

## Which Catalog Grades

By default, the catalog step 9 ran on (`expansion_hunter/<sample>_eh_catalog.json`) when it carries Stranger's thresholds: that is the STRchive catalog `setup.sh` builds ([step 9](09-str-expansions.md#which-catalog-runs)). Its thresholds are STRchive's, so steps 9 and 9b read one list of loci, and every graded record carries `SourceDisplay=STRchive-<version>`. The log names it. Without it (step 9 ran on the bundled 31 loci, or on a catalog without thresholds) Stranger's bundled GRCh38 catalog grades, as before. `STRANGER_REPEATS=file` picks another catalog (TSV or JSON), `STRANGER_REPEATS=bundled` the bundled GRCh38 one. The Nextflow `STRANGER` module takes `--expansion_catalog` and makes the same choice.

The thresholds differ from the bundled catalog's at some loci. HTT, for example: STRchive's benign maximum is 26 and its pathogenic minimum 36, so 27 to 35 CAG is `pre_mutation`, where the bundled catalog calls up to 34 normal. Seven STRchive loci get no status at all (MUC1, CEL, TYMS, PLIN4, VWA1, POLG, MIR7-2; [why](09-str-expansions.md#which-catalog-runs)).

## Command

```bash
./scripts/09b-stranger.sh your_name
```

Step 09 must run first:

```bash
./scripts/09-expansion-hunter.sh your_name male   # or female
./scripts/09b-stranger.sh your_name
```

`STRANGER_REPEATS` picks the catalog (above). `ALIGN_DIR` (default `aligned`) names the BAM folder for the RFC1 motif count. The output is written to a temporary name and renamed when Stranger succeeds, so a failed run leaves no empty or partial VCF that a rerun would take as done.

## Output

| File | Description |
|---|---|
| `expansion_hunter/<sample>_eh_stranger.vcf` | Annotated VCF with `STR_STATUS` and disease metadata in INFO fields |
| `expansion_hunter/<sample>_rfc1_motif.tsv` | Only when Stranger flags RFC1: the repeat motif of the reads at the locus and a verdict (see the RFC1 note below) |
| `expansion_hunter/<sample>_str_motifs.tsv` | Only with the STRchive catalog, and only when Stranger flags one of its other motif loci: the same columns, one row per flagged locus (see the motif note below) |

Key INFO fields added by Stranger:

| Field | Values | Meaning |
|---|---|---|
| `STR_STATUS` | `normal`, `pre_mutation`, `full_mutation` | Pathogenicity call for this allele |
| `Disease` | e.g. `HD` | Disease associated with this locus |
| `InheritanceMode` | `AD`, `AR`, `XD`, `XR` | Inheritance mode |
| `STR_NORMAL_MAX` | integer | Upper bound of normal repeat range |
| `STR_PATHOLOGIC_MIN` | integer | Lower bound of clearly pathogenic range |
| `HGNCId`, `Source`, `SourceId`, `SourceDisplay` | | Gene and the source of the locus definition (`STRchive`, its locus id and release with the STRchive catalog) |

## Runtime

About a second, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)). Stranger is a pure-Python VCF annotator — it reads the VCF once and writes to stdout. The RFC1 motif count, when RFC1 is flagged, reads only the reads at that one locus; it has not been timed on a whole genome.

## Notes

- Requires step 09 (ExpansionHunter) to have run first. The script exits cleanly with an informational message if the EH VCF is absent.
- The STRchive catalog grades 75 of its 82 loci; Stranger's bundled GRCh38 catalog has thresholds for 29 of ExpansionHunter's 31 bundled loci (none for NIPA1 and NOTCH2NL).
- `STR_STATUS` is per-allele: a heterozygous locus may have one normal and one pre-mutation allele.
- Pre-mutation alleles at FMR1 (55-200 CGG) and ATXN1 carry carrier risk even without current disease. See `docs/interpreting-results.md` for guidance.
- Short-read WGS has limited ability to size very large expansions (>150 repeats) accurately; `full_mutation` calls at loci like FMR1 and C9orf72 should be confirmed with orthogonal methods.
- **RFC1 (CANVAS): the flag is resolved by the motif of the reads.** CANVAS requires the **AAGGG** motif (or ACAGG), **biallelic**, at **~400–2000+** repeats. The two catalogs size RFC1 differently. The bundled one sizes only the degenerate `AARRG` unit, and Stranger's bundled thresholds flag any allele above the reference's 11 units as `full_mutation`, so a modest expansion of the common benign AAAAG motif (e.g. 51/73) is flagged `full_mutation`. The STRchive catalog sizes two repeats, AAAAG then AAGGG, and Stranger grades only the AAGGG one: 12 to 399 units is `pre_mutation`, 400 or more `full_mutation`, and an expansion of AAAAG alone is not flagged. When Stranger flags RFC1, the script runs `bin/rfc1_motif.py` in `PYPGX_IMAGE` on the sample's BAM (`aligned/<sample>_sorted.bam`). It takes the bases of each read over the repeat and 50 bases either side with the soft-clipped bases next to them, and the whole sequence of unmapped reads placed at their mate's position within 1,000 bases of the repeat and of mates mapped elsewhere, finds the longest tandem run of one pentamer in each read (any rotation, either strand) and counts the read when the run is longer than the reference allele's 11 units, so reads of an unexpanded allele do not outvote those of an expanded one, and writes `<sample>_rfc1_motif.tsv` with the share of each motif and a verdict:

  | Verdict | When | The report says |
  |---|---|---|
  | `pathogenic_motif` | 3 or more reads, and at least 20% of them, carry AAGGG or ACAGG | confirm with repeat-primed PCR |
  | `benign_motif` | otherwise, AAAAG or AAAGG carries at least half of the reads | the flag is not supported |
  | `uncertain_motif` | otherwise, with 3 or more reads | the dominant motif is not a known benign one: confirm if symptoms fit |
  | `not_resolved` | fewer than 3 reads carry a run | the motif is not resolved |
  | `not_checked` | no BAM | confirm the flag with repeat-primed PCR |

  A rerun of the script keeps an existing table unless it is `not_checked`, so a BAM added later gets counted.

  In the Nextflow pipeline, `STRANGER_RFC1` runs the same count after `STRANGER`, in the same image, on the BAM of the samplesheet row (or the one the pipeline aligned), and writes `<sample>/expansion_hunter/<sample>_rfc1_motif.tsv` only when Stranger flags RFC1. `HTML_REPORT` gets the table, so the Nextflow report shows the same verdict as step 24's.

  The count reads the motif; it does not size a long expansion, and a heterozygous AAGGG allele is a carrier state, not CANVAS. CI checks it on made-up reads (`tests/test_rfc1_motif.py`, `tests/e2e/expansions-and-rfc1-1-motif.sh`); how it does on real carriers, whose reads the aligner may clip or misplace differently, has not been measured. Confirm any `pathogenic_motif` or `uncertain_motif` with motif-aware, flanking-PCR or repeat-primed-PCR testing, and only if cerebellar-ataxia/neuropathy/vestibular symptoms are present.

  With the STRchive catalog, the count takes RFC1's motifs and its two repeats (AAAAG and AAGGG, chr4:39348424-39348483) from the catalog: benign AAAAG and AAAGG, pathogenic AAGGG, ACAGG and CAGGG. STRchive lists AAAGG as pathogenic, and its own notes say why: expansions of 600 or more AAAGG units are seen in CANVAS, while 100 to 380 units are seen in unaffected people ([PMID 37450567](https://pubmed.ncbi.nlm.nih.gov/37450567/)). Short reads cannot size an allele that long, so the count keeps AAAGG benign, as [Cortese et al. 2019](https://pubmed.ncbi.nlm.nih.gov/30926972/) classed it. On the CI's made-up reads the verdicts are the same as with the built-in motifs.
- **The other loci where the motif decides.** At the FAME loci (familial adult myoclonic epilepsy: SAMD12, STARD7, MARCHF6, YEATS2, TNRC6A, RAPGEF2, RAI1; [Ishiura et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29507423/)), SCA31 (BEAN1), SCA37 (DAB1) and aFTLD-U (GOLGA8A), disease comes from an inserted motif (TTTCA at the FAME loci), while an expansion of the reference motif is benign. The STRchive catalog names both classes of motif for each. When Stranger flags one of them, step 9b runs the same count (`bin/rfc1_motif.py --flagged`) with that locus's motifs and writes one row per flagged locus to `<sample>_str_motifs.tsv`, with the verdicts of the RFC1 table. A recessive locus's note asks for a pathogenic expansion on both copies; these are dominant. A rerun redoes a `not_checked` table, as for RFC1. `STRANGER_RFC1` in Nextflow writes the same file. The HTML report does not read this table yet: open it beside the report. CI checks it on made-up SAMD12 reads (`tests/test_str_catalog.py`, `tests/e2e/str-catalog-1-strchive.sh`); how it does on real carriers has not been measured.
