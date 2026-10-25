# Step 45: Y-STR Markers

## What This Does

Reads the repeat counts of a male sample's Y-chromosome short tandem repeats (Y-STRs) from the genome, for the 111 markers of FamilyTreeDNA's Y-111 test, and lists them in FamilyTreeDNA's order. The Y-12, Y-25, Y-37 and Y-67 tests are the first 12, 25, 37 and 67 of them. Opt-in.

## Why

Surname projects and Y-DNA matching compare men by these counts. A whole genome covers most of the markers, so a sequenced man can compare himself with a project's members without a separate test. Step 37 gives the haplogroup, the deep branch of the paternal line; the STR counts change faster and separate the recent branches of one surname.

## Tool

- **HipSTR** 0.7 (Willems et al., Nat Methods 2017, PMID 28436466; GPL-2.0), which genotypes STRs from aligned reads and models PCR stutter. It was validated on Y-STRs against capillary electrophoresis calls (Willems et al., Am J Hum Genet 2016, PMID 27126583).
- `bin/ystr_report.py` turns HipSTR's calls into counts in FamilyTreeDNA's order. It runs with the Python of the HipSTR image.

## Docker Image

- `HIPSTR_IMAGE` (Bioconda's HipSTR build)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## The marker table

`assets/ystr_panel.tsv` holds one row per marker, in FamilyTreeDNA's order. Its header names each source:

- **Order, names and panels**: FamilyTreeDNA's "Y-STR Results Guide" (https://help.familytreedna.com/hc/en-us/articles/4408063356303-Y-STR-Results-Guide), read on 2026-10-25 through the Internet Archive copy of 2026-04-13, because the live page refuses scripted requests.
- **Regions** (GRCh38, the motif length and the reference copies): the Y-STR set HipSTR's authors annotated, `human/all_annot_markers.hg38.fixed.bed` of HipSTR-references at a fixed commit (MIT licence). Some markers are several regions there (DYS389II, DYS448, DYS449, DYS450, DYS504, DYS510, DYS532, DYS552); the table keeps those parts as the source gives them.
- **Motif**: NIST STRBase's repeat structure, where it has one.

The 111 markers fall into four kinds:

| Kind | Markers | What the step does |
|---|---|---|
| single region | 76 | HipSTR genotypes it |
| several regions | 2: DYS449, DYS552 | HipSTR genotypes each part; the count is their sum |
| span | 6: DYS389II, DYS448, DYS450, DYS504, DYS510, DYS532 | HipSTR genotypes one region from the first part's start to the last part's end; the count is the parts' reference copies plus the change HipSTR finds. DYS389II (164 bp) and DYS448 (156 bp) are longer than 150 bp reads can span, so from such reads they are not callable |
| multi-copy | 16: DYS385a/b, DYS459a/b, DYS464a-d, YCAIIa/b, CDYa/b, DYF395S1a/b, DYS413a/b | not called. Their copies lie in duplicated, palindromic Y sequence, and short reads cannot tell one copy from another |
| no region | 11: DYS439, DYS447, DYS437, DYS710, DYS716, DYS441, Y-GGAAT-1B07, DYS712, DYS650, DYS715, DYS726 | not called: the marker set has no region for them |

Why the span markers are one region: HipSTR 0.7 genotypes every region on its own. When the parts of these six went to it one by one, on HG002's 60x reads it wrote no record for one part of each, though their neighbours had 25 to 70 reads, so none of the six got a count. Their parts lie 7 to 47 bp apart, or overlap another marker's region (DYS389II starts where DYS389I starts). DYS449 and DYS552, whose parts it did call, stay as parts.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps), and the reference with its `.fai`
- `${SAMPLE}/indexcov/indexcov-indexcov.ped` from step 16: the sex indexcov infers from the reads. The step runs only when it is male, as step 37 does. A female or undetermined sample is skipped with one line and exit 0.

In the Nextflow pipeline (`y_str` in `--tools`) the `Y_STR` process runs for the rows whose samplesheet sex is male; `INDEXCOV` has checked that sex against the reads.

## Command

```bash
./scripts/45-y-str.sh your_name
```

With `run-all.sh`: `TOOLS=...,y_str` (opt-in: a default run does not list it; with `y_str` named, a female sample is listed as skipped). There is no data to install.

Inside one container the step writes HipSTR's region file from the table, then runs:

```bash
HipSTR --bams <bam> --fasta <reference> --regions regions.bed --str-vcf <sample>_y_str.vcf.gz \
  --bam-samps <sample> --bam-libs <sample> --haploid-chrs chrY --def-stutter-model --min-reads 5 \
  --max-str-len 150 --output-filters --log <sample>_y_str.log
```

- `--haploid-chrs chrY`: one Y, so one allele per marker.
- `--def-stutter-model` and `--min-reads 5`: HipSTR learns a stutter model per locus from many samples and by default skips a locus with fewer than 100 reads; for one sample its README recommends the default stutter model and a lower minimum.
- `--max-str-len 150`: HipSTR skips regions over 100 bp by default; the table decides which repeats are too long (below).
- `--bam-samps` and `--bam-libs`: every read is the sample's, whatever read groups the BAM has.

## Output

| File | Contents |
|---|---|
| `y_str/${SAMPLE}_y_str.tsv` | one row per marker in FamilyTreeDNA's order: `order`, `panel` (the smallest test that has it), `marker`, `count`, `confidence`, `reads` (reads HipSTR used), `stutter_reads`, `flank_indel_reads`, `q` (HipSTR's posterior), `note` |
| `y_str/${SAMPLE}_y_str.csv` | `marker,value,confidence` in the same order, for a spreadsheet or a project's table |
| `y_str/${SAMPLE}_y_str.vcf.gz` | HipSTR's calls |
| `y_str/${SAMPLE}_y_str.log` | HipSTR's log: why it skipped a region |

How a count is made: HipSTR gives, per region, how many base pairs the call differs from the reference allele. The count is the region's reference copies (rounded to whole repeats, halves up) plus that difference divided by the motif length. Base pairs left over are written after a dot, the usual notation of a partial repeat (`11.2`: eleven repeats and two bases).

Confidence:

| Value | Meaning |
|---|---|
| `high` | HipSTR's own haploid filters pass (posterior at least 0.9, at most 15% of the reads with a stutter artifact and at most 15% with an indel in the flanks), and at least 10 reads were used |
| `low` | a call that fails one of those; the note says which |
| `no call` | HipSTR wrote no genotype for the marker or for one of its parts. When it wrote no record at all, the note quotes the reason from HipSTR's log, such as `Skipping locus with too few reads`, or `Aborting genotyping of the locus as the sequence downstream of the repeat is too repetitive for accurate genotyping` when the bases beside the repeat are themselves repetitive |
| `not callable` | multi-copy, no region, or a repeat longer than 150 bp reads can span with 10 bases on each side (DYS452, 155 bp in the reference, and the DYS389II and DYS448 spans) |

Both reports get a Y-STR section: the number of markers per confidence and every called marker with its count.

## Comparing with FamilyTreeDNA

- **A count can differ from FamilyTreeDNA's by a fixed number of repeats.** Labs do not all count the same repeat blocks of a marker, and the region edges of the marker set are not always the edges a lab counts. FamilyTreeDNA does not publish its conventions marker by marker, and we found no published table of the offsets between these regions and FamilyTreeDNA's counts (searched 2026-10-25), so the table's `offset` column is empty and no count is adjusted. A difference that is the same for every man of a project is a convention, not a mutation. The safe comparison is with a relative or project member who also has a genome, or with one known FamilyTreeDNA result of the same man.
- FamilyTreeDNA reports DYS389II as the sum total of DYS389I and DYS389II (its guide says so). The DYS389II span here starts where DYS389I starts, so its count would cover both too. At 164 bp it is longer than 150 bp reads can span, so from such reads DYS389II is not callable; DYS389I is.
- FamilyTreeDNA lists DYS385, DYS459, DYS464, YCAII, CDY, DYF395S1 and DYS413 as two or four values; those are the not callable rows here.

## Runtime

Seconds to a few minutes: HipSTR reads 86 short regions of chrY.

## Notes

- Short reads at 30x put about 15 reads on each Y region (one copy). A repeat near the read length gets fewer reads that span it, so its call is `low` or `no call` more often.
- HipSTR is built for Illumina short reads; long-read BAMs are not tested here.
- The e2e case reads HG002's reads at every region from GIAB's 60x BAM, checks each high call, the span markers included, against a count made from the reads themselves, keeps HipSTR's log with the case logs, and checks the step against the Nextflow process on the fixture's own BAM.
