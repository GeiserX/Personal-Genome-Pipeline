# Step 9b: STR Clinical Annotation (Stranger)

## What This Does

Annotates the ExpansionHunter VCF (step 9) with clinical pathogenicity status for each repeat locus. Raw repeat counts become labelled calls: **normal**, **pre_mutation**, or **full_mutation**, with disease name, OMIM number, inheritance mode, and the specific repeat-size thresholds used for classification.

This turns step 9 output from a table of numbers into actionable clinical calls.

## Why

ExpansionHunter reports the number of repeats at each locus but applies no pathogenicity judgement. Stranger adds that judgement using the ClinGen/OMIM short tandem repeat database, so you can immediately see whether a repeat count falls in the normal, pre-mutation, or full-mutation range without manually cross-referencing disease thresholds.

## Tool

- **Stranger** v0.10.2 (Clinical Genomics Stockholm) — annotates STR VCFs with pathogenicity labels from a curated repeat catalog covering HTT, FMR1, C9orf72, DMPK, ATXN*, RFC1, and ~40 other loci

## Docker Image

- `STRANGER_IMAGE`
- `PYPGX_IMAGE` for the RFC1 motif count (it has pysam)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

- Binary: `stranger` (on PATH)
- Bundled repeat catalog: Stranger's GRCh38 catalog (`variant_catalog_grch38.json`, installed inside the container). Stranger's own default is its GRCh37 catalog; the script and the Nextflow module pass the GRCh38 one explicitly, because ExpansionHunter (step 9) calls the GRCh38 loci

## Command

```bash
./scripts/09b-stranger.sh your_name
```

Step 09 must run first:

```bash
./scripts/09-expansion-hunter.sh your_name male   # or female
./scripts/09b-stranger.sh your_name
```

A custom repeat catalog (TSV or JSON) can be supplied via the `STRANGER_REPEATS` environment variable; the bundled GRCh38 catalog is used when it is not set. `ALIGN_DIR` (default `aligned`) names the BAM folder for the RFC1 motif count. The output is written to a temporary name and renamed when Stranger succeeds, so a failed run leaves no empty or partial VCF that a rerun would take as done.

## Output

| File | Description |
|---|---|
| `expansion_hunter/<sample>_eh_stranger.vcf` | Annotated VCF with `STR_STATUS` and disease metadata in INFO fields |
| `expansion_hunter/<sample>_rfc1_motif.tsv` | Only when Stranger flags RFC1: the repeat motif of the reads at the locus and a verdict (see the RFC1 note below) |

Key INFO fields added by Stranger:

| Field | Values | Meaning |
|---|---|---|
| `STR_STATUS` | `normal`, `pre_mutation`, `full_mutation` | Pathogenicity call for this allele |
| `Disease` | e.g. `HD` | Disease associated with this locus |
| `InheritanceMode` | `AD`, `AR`, `XD`, `XR` | Inheritance mode |
| `STR_NORMAL_MAX` | integer | Upper bound of normal repeat range |
| `STR_PATHOLOGIC_MIN` | integer | Lower bound of clearly pathogenic range |
| `HGNCId`, `Source`, `SourceId` | | Gene and the source of the locus definition |

## Runtime

About a second, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)). Stranger is a pure-Python VCF annotator — it reads the VCF once and writes to stdout. The RFC1 motif count, when RFC1 is flagged, reads only the reads at that one locus; it has not been timed on a whole genome.

## Notes

- Requires step 09 (ExpansionHunter) to have run first. The script exits cleanly with an informational message if the EH VCF is absent.
- The bundled catalog covers ~40 STR loci with established clinical thresholds. Custom catalogs can be used via `STRANGER_REPEATS=/path/to/catalog.tsv`.
- `STR_STATUS` is per-allele: a heterozygous locus may have one normal and one pre-mutation allele.
- Pre-mutation alleles at FMR1 (55-200 CGG) and ATXN1 carry carrier risk even without current disease. See `docs/interpreting-results.md` for guidance.
- Short-read WGS has limited ability to size very large expansions (>150 repeats) accurately; `full_mutation` calls at loci like FMR1 and C9orf72 should be confirmed with orthogonal methods.
- **RFC1 (CANVAS): the flag is resolved by the motif of the reads.** CANVAS requires the **AAGGG** motif (or ACAGG), **biallelic**, at **~400–2000+** repeats. ExpansionHunter sizes only the degenerate `AARRG` unit, and the catalog's `STR_PATHOLOGIC_MIN` for RFC1 is far below the clinical threshold, so a modest expansion of the common benign AAAAG motif (e.g. 51/73) is flagged `full_mutation`. When Stranger flags RFC1, the script runs `bin/rfc1_motif.py` in `PYPGX_IMAGE` on the sample's BAM (`aligned/<sample>_sorted.bam`). It takes the bases of each read over the repeat and 50 bases either side with the soft-clipped bases next to them, and the whole sequence of unmapped reads placed at their mate's position within 1,000 bases of the repeat and of mates mapped elsewhere, finds the longest tandem run of one pentamer in each read (any rotation, either strand) and counts the read when the run is longer than the reference allele's 11 units, so reads of an unexpanded allele do not outvote those of an expanded one, and writes `<sample>_rfc1_motif.tsv` with the share of each motif and a verdict:

  | Verdict | When | The report says |
  |---|---|---|
  | `pathogenic_motif` | 3 or more reads, and at least 20% of them, carry AAGGG or ACAGG | confirm with repeat-primed PCR |
  | `benign_motif` | otherwise, AAAAG or AAAGG carries at least half of the reads | the flag is not supported |
  | `uncertain_motif` | otherwise, with 3 or more reads | the dominant motif is not a known benign one: confirm if symptoms fit |
  | `not_resolved` | fewer than 3 reads carry a run | the motif is not resolved |
  | `not_checked` | no BAM, or the Nextflow pipeline (its STRANGER process gets no BAM) | confirm the flag with repeat-primed PCR |

  The count reads the motif; it does not size a long expansion, and a heterozygous AAGGG allele is a carrier state, not CANVAS. CI checks it on made-up reads (`tests/test_rfc1_motif.py`, `tests/e2e/expansions-and-rfc1-1-motif.sh`); how it does on real carriers, whose reads the aligner may clip or misplace differently, has not been measured. Confirm any `pathogenic_motif` or `uncertain_motif` with motif-aware, flanking-PCR or repeat-primed-PCR testing, and only if cerebellar-ataxia/neuropathy/vestibular symptoms are present.
