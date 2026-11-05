# Step 6: ClinVar Pathogenic Variant Screening

## What This Does
Intersects your sample VCF against the ClinVar database of known pathogenic variants, identifying any positions where your genome carries a clinically reported disease variant.

## Why
ClinVar is the most widely used public database of clinically reported variants. This screen catches pathogenic SNPs and indels that have been submitted by clinical labs — carrier status, dominant disease risk, and pharmacogenomic flags. Note that ClinVar entries vary in evidence quality (see star ratings in [interpreting-results.md](interpreting-results.md)).

## Tool
- **bcftools norm**, **isec** and **annotate** — split and left-align both files, keep the sample's records whose allele is in ClinVar, and copy ClinVar's gene, significance and review status onto them
- **bin/clinvar_join_adjacent.awk** with **bcftools +fill-from-fasta** — build the combined allele of two neighbouring records, so a delins the caller wrote as two records is matched as one allele too
- **`bin/clip_check.py`** (pysam) — the cDNA read check: whether each hit's alt reads look like spliced cDNA rather than genomic DNA, and which hits sit at a known lab-construct site
- **`bin/cdna_scan.py`** (pysam) — the cDNA scan: the loci of the whole genome where reads jump from one exon to the next
- **`bin/clinvar_hits.awk`** — the hits by review stars, and the always-report list of known deep intronic sites (`mode=sites`)

## Docker Image
- `BCFTOOLS_IMAGE`
- `PYPGX_IMAGE`, for pysam (the cDNA read check and the cDNA scan; step 32 uses the same image)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- Sample VCF from DeepVariant (step 3)
- Reference FASTA and its `.fai`: `reference/GRCh38_no_alt_analysis_set.fasta` (used to left-align indels and to fill in the bases between the two parts of a combined allele)
- `clinvar_pathogenic_chr.vcf.gz` from reference setup (step 00) — chr-prefixed, filtered to Pathogenic/Likely_pathogenic only
- For the always-report list, optional: the sample's gVCF, `vcf/<sample>.g.vcf.gz` with its `.tbi` (step 3 writes it). Without it a site with no record in the VCF is "no call", since nothing tells a reference call from no coverage
- For the cDNA read check and the cDNA scan, optional: the sample's BAM, `aligned/<sample>_sorted.bam` with its `.bai` (`ALIGN_DIR` picks another folder), and `reference/gencode.v50.basic.exons.bed`, the exons of GENCODE's basic annotation that `setup.sh` installs. Without either, the step still screens and writes every hit as "not checked", with the reason

## Command
```bash
export GENOME_DIR=/path/to/data
./scripts/06-clinvar-screen.sh <sample_name>

# For long-read Clair3 output:
VCF_DIR=vcf_clair3 ./scripts/06-clinvar-screen.sh <sample_name>

# Another call floor (default GQ 30 and DP 10):
CLINVAR_MIN_GQ=20 CLINVAR_MIN_DP=8 ./scripts/06-clinvar-screen.sh <sample_name>

# Another always-report list (same columns as assets/deep_intronic_sites.tsv):
CLINVAR_ALWAYS_REPORT=/path/to/sites.tsv ./scripts/06-clinvar-screen.sh <sample_name>

# Without the genome-wide cDNA scan:
CDNA_SCAN=false ./scripts/06-clinvar-screen.sh <sample_name>
```

In the Nextflow pipeline the same floor is `ext.min_gq` and `ext.min_dp` of `CLINVAR_SCREEN`, set in a config: `withName: 'CLINVAR_SCREEN' { ext.min_gq = 20; ext.min_dp = 8 }`. The cDNA scan is skipped with `withName: 'CLINVAR_CDNA_CHECK' { ext.cdna_scan = false }`. The Nextflow pipeline always reads `assets/deep_intronic_sites.tsv`.

### What the Script Does

1. Lists the contigs each file holds records on (`bcftools index -s`) and keeps those the reference has too: ClinVar's contigs that are in the reference, and the sample's contigs among those. `bcftools norm` stops at the first record whose contig the reference lacks, such as a provider's unplaced scaffold or one of the `NT_` contigs the full ClinVar file has, and a record on a contig the other file lacks cannot match anyway. The step prints how many of the sample's records it leaves out and on which contigs. When no contig is shared it stops with an error instead of reporting zero hits, and says which file is named the other way (`1, 2, ...` against `chr1, chr2, ...`).
2. Builds `clinvar/clinvar_pathogenic_chr.norm.vcf.gz` beside the ClinVar file the first time, and again whenever the source file is newer: ClinVar's records on the reference's contigs, multiallelic records split (`bcftools norm -m -any`) and indels left-aligned against the reference. Every sample then reuses it. The Nextflow module reuses this file when it sits beside the `--clinvar` file and is newer; otherwise it normalises ClinVar inside the task.
3. Filters the sample VCF to PASS records. If the VCF has no PASS record at all (callers that leave FILTER as `.`), it uses `-f .,PASS` instead and prints a notice saying so. If no record is left, the step stops with an error.
4. Splits and left-aligns the sample the same way. `bcftools isec` matches only identical REF/ALT, so a pathogenic allele inside a multiallelic record (genotype `1/2`) is found only after this split. The split writes the other alt allele as missing, not as the reference (`bcftools norm --multi-overlaps .`): a `1/2` record gives `1/.` and `./1`, where the default would give `1/0` and `0/1`. Both copies of such a call carry an alt allele, so neither part may read as het. A poly-T tract is where this shows most: a call of 5 T on one copy and 9 T on the other, against a reference of 7, is one `1/2` record, and its 5 T allele is a ClinVar hit. The reports word such a hit "one copy, other copy not called ref (1/.)", never het; the ACMG and carrier lists count it as one copy of two. This step does not size the tract: a short-read call of a homopolymer length can be off by one or two bases, so check the reads.
5. Builds the combined alleles of neighbouring records. A caller can write one complex allele, a delins, as two records a few bases apart on the same haplotype, for example a deletion and an insertion. Matched one by one, the parts can hit two single-submitter ClinVar records and miss the better reviewed record ClinVar holds for the whole allele. `bcftools norm -m +any` joins only records at the same position, so it cannot build that allele. `bin/clinvar_join_adjacent.awk` pairs two records when the second starts at most 10 bases after the first ends (records that overlap are never paired) and their genotypes allow both ALT alleles on one haplotype: phased with the same genotype in the same phase set (`FORMAT/PS`), or unphased with the same genotype (`0/1` and `1/0` count as the same het call; an unphased call with a missing allele such as `./1` is never paired). A phased call with a missing allele, such as `1|.`, has its ALT on a known haplotype: it pairs with a record in the same phase set whose ALT sits on that haplotype (`1|0` or `1|.`), and the combined record keeps `.` on the haplotype where either part is missing. That `.` means the copy was not called there; when the part comes from a split multiallelic record, it holds the other alt allele. Two unphased het calls can also sit on different haplotypes, and the genotype cannot tell which, so the notes of a combined allele and its parts built from them end with "unphased, so the parts may sit on different haplotypes". The combined record runs from the first part's start to the second part's end; `bcftools +fill-from-fasta` takes the bases between the parts from the reference, and `bcftools norm` trims and left-aligns it the way ClinVar was. The gene is not checked here, because the sample VCF has no gene annotation yet: a combined allele only counts when ClinVar holds a record for it. The step prints how many combined alleles it built. Only pairs are built: an allele split into three or more records is matched by its pairs only.
6. Keeps the sample's records and combined alleles whose allele is in ClinVar (`bcftools isec -n=2 -w1`), copies ClinVar's `ID`, `GENEINFO`, `CLNSIG` and `CLNREVSTAT` onto them (`bcftools annotate --pair-logic exact`), and keeps only those whose genotype carries an ALT allele (`bcftools view -i 'GT~"[1-9]"'`, a non-zero allele index anywhere in the genotype). The match is by allele, so without this a `0/0` or `./.` record, the `0/0` half of a split multiallelic record and a reference-only row with ALT `.` would be listed as hits. A half call such as `./1` carries the ALT allele and is a hit; `GT="alt"` would drop it. The filter runs after the split, where it sees each allele's own genotype. No command line goes into the hits file's header (`--no-version`), so ClinVar's path on your machine stays out of it.
7. Marks the hits that belong together. A matched combined allele is a hit of its own and carries `INFO/PGP_JOINED_FROM`, its two parts as `CHROM:POS:REF:ALT`. Each of its parts that is a hit too keeps its row and gets `INFO/PGP_PART_OF`, the combined record it belongs to, so the reader sees one allele with the combined record's own review status, not two unrelated hits. Two hits in one gene that could be one allele by the rule of item 5 but matched no combined record each get `INFO/PGP_ONE_ALLELE`, the other hit: possibly one allele, check ClinVar for a combined record.
8. Writes a call floor other than the default into the hits file's header, for example `##PGP_CALL_FLOOR=GQ=20,DP=8` (below). At the default the file has no such line and reads as GQ 30 and DP 10.
9. Writes the always-report list (below): the sample's call at every site of `assets/deep_intronic_sites.tsv`.
10. Runs the cDNA read check on each hit and the cDNA scan of every exon (below).
11. Prints the hits grouped by ClinVar's review stars (`bin/clinvar_hits.awk`): 4 practice guideline, 3 expert panel, 2 multiple submitters with no conflict, 1 a single submitter or conflicting classifications, 0 no assertion criteria, each with its combined-allele note when it has one and "LOW-CONFIDENCE CALL (GQ .., DP ..): check the reads" when the call is below the floor, and the count of low-confidence calls. The screening file keeps every Pathogenic/Likely_pathogenic submission whatever its review status, so a zero-star hit counts as a hit; the stars say how much weight it deserves. The Nextflow module runs the same steps and prints the same counts; its `CLINVAR_CDNA_CHECK` process runs the cDNA read check and the cDNA scan.

A delins that the caller writes as two records, while ClinVar holds it as one, is matched exactly once, as the combined record of item 5: `tests/test_clinvar_join_adjacent.py` check 8 runs the real step on such a fixture and counts the hits.

### The always-report list

Some pathogenic alleles sit deep inside an intron, far from any exon. Filters built around exons and splice sites skip them, and a report that lists only hits cannot say whether such a site was read as the reference or not read at all ([Vaz-Drago 2017](https://pubmed.ncbi.nlm.nih.gov/28497172/)). `assets/deep_intronic_sites.tsv` lists the best known ones: CFTR c.3718-2477C>T and c.1680-886A>G, CEP290 c.2991+1655A>G, USH2A c.7595-2144A>G, ATM c.5763-1050A>G and GLA c.640-801G>A, each a ClinVar Pathogenic record with an expert-panel or multi-submitter review, at the GRCh38 position and alleles of ClinVar's canonical SPDI (checked 2026-11-04).

Every run writes `clinvar/${SAMPLE}_clinvar_deep_intronic.tsv`, one row per site, carrier or not. Its `call` is the genotype in the words the hits use (het, hom, ...) when the sample carries the site's allele; `ref` for a reference call, with the GQ of the `0/0` record or of the gVCF reference block over the site (`ref (low confidence)` when that GQ is below the call floor, `CLINVAR_MIN_GQ` or `ext.min_gq`, default 30); `other allele` when the sample carries another allele there; or `no call`, with the reason (no record in the VCF or the gVCF, a gVCF record that is not a reference call, or no gVCF at all). `detail` holds the GQ and DP, the other allele or the reason. The records read are the split PASS records of item 4 at the site's position (`bcftools query -T`) and the gVCF records over it (`bcftools query -R`). Both reports print every row under the ClinVar hits.

### Lab-construct sites and the cDNA scan

`assets/construct_variants.tsv` lists alleles that lab expression constructs and transgenic disease models carry: APP Swedish and London, PSEN1 M146L, A246E and L286V, MAPT P301L and P301S, SNCA A53T, SOD1 G93A and LRRK2 G2019S, at ClinVar's GRCh38 positions (checked 2026-11-04). Plasmid or cDNA of these genes in a library gives P/LP calls at exactly these sites ([Kim 2020](https://pubmed.ncbi.nlm.nih.gov/32814883/); [Kim 2016](https://pubmed.ncbi.nlm.nih.gov/27334474/)). A hit at one of them gets the site's name in the read-check table's `construct_site` column, and the step's log and both reports print it with its read counts, whatever the verdict ("KNOWN CONSTRUCT SITE ..."; a hit with no read counts, such as one with no BAM, gives the reason instead). The reports list the best-reviewed 50 hits, and a construct-site hit below them still gets its own line. A real carrier is possible at every site, LRRK2 G2019S above all, so this says to look at the reads, not that the call is wrong.

The read check sees only the reads at a hit. Contamination that carries no ClinVar allele, or whose junction reads miss the hit, needs a wider look: `bin/cdna_scan.py` reads the reads over every exon of the exon file. A read counts for a junction when one piece of it has a clip of at least 10 bases within 3 bases of an exon end and another piece, on the same contig and strand, 30 bases to 1 Mb further on and next along the read as aligned (starting within 3 read bases of where the first piece stops), a clip at an exon start (or a CIGAR `N` gap runs from an exon end to an exon start). A junction counts with at least 2 reads, junctions less than 20 kb apart form one locus, and a locus is listed with 2 junctions or more: one junction alone can be a chimeric read or a deletion. Duplicate and secondary alignments are skipped. Supplementary records are read too, so a read whose primary alignment lies outside every exon (on a vector sequence, say) still counts, and a read counts once per junction however many of its records show it. The table, `clinvar/${SAMPLE}_clinvar_cdna_scan.tsv`, starts with a `# cdna_scan:` line saying what was read, or why the scan did not run (no BAM, no exon file, `CDNA_SCAN=false`). A listed locus is cDNA or plasmid contamination, or a processed copy of the gene the person carries ([ten Berk de Boer 2023](https://pubmed.ncbi.nlm.nih.gov/37323659/)); the reads tell which. Both reports print the loci; the HTML ClinVar card counts them.

### Low-confidence calls

A P/LP call can rest on very few reads: a novel frameshift at depth 5 and genotype quality 3 inside a GC-rich repeat is called like any other. Low-GQ calls are also where caller versions disagree. A hit whose call has a genotype quality (`FORMAT/GQ`) below `CLINVAR_MIN_GQ` (default 30) or a depth (`FORMAT/DP`) below `CLINVAR_MIN_DP` (default 10) is a **low-confidence call**. It stays in the hits, marked "low-confidence call (GQ .., DP ..): check the reads" in the step's output and in both reports, and the reports count it apart: the headline P/LP count and the counts by stars leave it out. A call without GQ or DP is compared on the value it has, and a call with neither is never marked. A floor other than the default is the `##PGP_CALL_FLOOR` line of the hits file, so the reports read the floor the run used; a file without the line used the default. A floor of 0 marks nothing.

### Population frequency

Both reports put the gnomAD allele frequency beside each hit: VEP's `gnomADe_AF` and `gnomADg_AF` from step 13's `vep/${SAMPLE}_vep.vcf.gz`, read through its `.tbi` at the hit's position, REF and ALT, and the higher of the two is shown. A hit above 1% is marked "common: unlikely to be a high-penetrance finding". A hit VEP has no record for is "gnomAD AF unknown". Without step 13's output, or without its `.tbi`, every hit says unknown and the report says why. The match is by position and alleles as both files hold them: an indel that VEP's input holds in another representation (a multiallelic record step 6 split and trimmed, for example) finds no record and reads as unknown, never as rare.

### The cDNA read check

A sample can hold spliced cDNA or expression-vector DNA, a known lab contaminant. Its reads align to a gene's exons. A read that crosses an exon junction is clipped exactly at the exon edge, and the rest of the read aligns to the next exon as a supplementary alignment (an `SA` tag). Variants the construct carries are then called like any other, and the screen lists them as ClinVar hits. Genomic reads at the same site run on into the intron instead. Neither haplocheck (mtDNA only) nor VerifyBamID2 (human-to-human mixture only) can see this.

For each hit, `bin/clip_check.py` reads the BAM at the hit's position and splits the reads by the allele they show (a SNV, an MNV, or a simple insertion or deletion; other alleles are not checked). Duplicate and secondary alignments are left out, and a read counts once. A read counts as clipped at an exon edge when it has a clip of at least 10 bases (soft or hard) within 3 bases of an exon edge on the same side: a clip on its right end at an exon's end, on its left end at an exon's start. Its own alignment and each of its `SA` alignments are tested. The hit is **flagged**, "possible cDNA/transgene contamination: check the reads", when at least 2 reads carry the alt allele, at least half of them are clipped that way, and at most 20% of the ref reads are. Otherwise it is **clear**, or **not checked** with the reason: no BAM, no exon file, fewer than 2 alt reads, a contig the BAM lacks, or a complex allele. Only the reads at the hits are read, so the check takes seconds.

A flag is a reason to look at the reads (in IGV, for example), not a verdict. A real variant in a sample that also holds the construct is flagged too, and a construct whose junction reads all fall outside the hit is not seen. Mates on the next exon are not used.

In the Nextflow pipeline the process `CLINVAR_CDNA_CHECK` runs the same `bin/clip_check.py` on `CLINVAR_SCREEN`'s hits, in the same image, with the sample's BAM from the samplesheet (or the one the pipeline aligned) and the exon BED of `--gencode_exons` (`run-all.sh` passes the one `setup.sh` installs). It writes `<sample>/clinvar/<sample>_clinvar_cdna_check.tsv`, the same table, and `<sample>_clinvar_cdna_scan.tsv`. A sample without a BAM, or a run without `--gencode_exons`, gets every hit as not checked, with the reason. `CLINVAR_SCREEN` reads the sample's gVCF (DEEPVARIANT's, or the samplesheet's `gvcf` column) for the always-report list. `HTML_REPORT` gets the three tables too, so the Nextflow HTML report shows the same cDNA column and flagged count as step 24's.

## Runtime

The [runtime table](hardware-requirements.md#runtime-per-step) gives the measured time of this step. That measurement is older than item 5 above, which streams the split sample VCF once more through awk and bcftools, and older than the cDNA scan, which reads every read over the exons of the BAM; neither has been timed on a whole genome. `CDNA_SCAN=false` leaves the scan out.

## Output

| File | Description |
|---|---|
| `clinvar/${SAMPLE}_clinvar_hits.vcf` | **The hits: your records that carry a ClinVar Pathogenic/Likely_pathogenic allele, with ClinVar's ID, gene (`GENEINFO`), significance (`CLNSIG`) and review status (`CLNREVSTAT`). Your genotype, GQ and DP are in the sample column, and a `##PGP_CALL_FLOOR` header line holds the call floor of a run that changed it.** A combined allele of two neighbouring records is a record of its own, with `PGP_JOINED_FROM`; its parts carry `PGP_PART_OF`, and two hits that may be one allele carry `PGP_ONE_ALLELE` |
| `clinvar/${SAMPLE}_clinvar_hits.tsv` | The same hits, one row each: `chrom`, `pos`, `ref`, `alt`, `genotype`, `clinvar_id`, `geneinfo`, `clnsig`, `clnrevstat` |
| `clinvar/${SAMPLE}_clinvar_cdna_check.tsv` | The cDNA read check, one row per hit: `chrom`, `pos`, `ref`, `alt`, `alt_reads`, `alt_clipped`, `ref_reads`, `ref_clipped` (the reads with each allele, and how many of them are clipped at an exon edge), `cdna_check` (`flagged`, `clear` or `not_checked`), `cdna_note` (the flag text or the reason it was not checked) and `construct_site` (the gene and name of a known lab-construct site, empty otherwise) |
| `clinvar/${SAMPLE}_clinvar_cdna_scan.tsv` | The cDNA scan: a `# cdna_scan:` line, then one row per listed locus: `chrom`, `start` and `end` (the first junction's exon end and the last junction's exon start, 1-based), `junctions`, `reads` (distinct reads, so a read across two junctions counts once) and `junction_list` (`end-start:reads` per junction) |
| `clinvar/${SAMPLE}_clinvar_deep_intronic.tsv` | The always-report list, one row per site: the site's `chrom`, `pos`, `ref`, `alt`, `gene`, `hgvs` and `clinvar_id`, then `call` and `detail` |
| `clinvar/${SAMPLE}_pass.vcf.gz` | Filtered, split and left-aligned sample VCF (intermediate). It holds the records on the contigs shared with ClinVar and the reference only; versions before the shared-contig restriction kept every contig |
| `clinvar/clinvar_pathogenic_chr.norm.vcf.gz` (in `${GENOME_DIR}`) | Normalised ClinVar, shared by every sample |

Both reports (step 24 and `generate-report.sh`) read the hits file and show gene, genotype (het, hom, or "one copy, other copy not called ref" for the half of a split `1/2` call), significance, review status, stars and the combined-allele note ("combined allele of ...", "part of a combined allele: ...", both followed by "; unphased, so the parts may sit on different haplotypes" on an unphased het call, or "possibly one allele with ...: check ClinVar for a combined record") for each hit, the best-reviewed first, with the ClinVar file's release date, the gnomAD frequency and the low-confidence mark. The headline count leaves the low-confidence calls out and a line beside it counts them. When the cDNA check table is there, both print its result beside each hit (the HTML hits table has a "cDNA read check" column) and the ClinVar card counts the flagged hits. Older versions of this step wrote `clinvar/isec/0002.vcf`; that file holds only the sample's side of the intersection, with no gene or significance, and nothing reads it any more. Rerun step 6 to get the hits file.

`${SAMPLE}_clinvar_hits.tsv` is the table to read or load elsewhere. To list the hits from the VCF yourself:

```bash
source versions.env   # from the repository root
docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
  bcftools query -f '%CHROM:%POS %REF>%ALT [%GT] %INFO/GENEINFO %INFO/CLNSIG %INFO/CLNREVSTAT\n' \
  /genome/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf
```

## Interpreting Results

This step screens against **Pathogenic and Likely_pathogenic variants only** — benign and VUS entries are excluded at the database level (see step 00 reference setup). Every hit is an allele that ClinVar classifies as Pathogenic or Likely_pathogenic, and the hit shows your genotype for it.

| Scenario | Meaning | Action |
|---|---|---|
| Heterozygous + autosomal recessive | Healthy carrier | Note for family planning only |
| Homozygous + autosomal recessive | Possibly affected — requires clinical confirmation | Investigate — confirm with clinical evaluation and ClinVar review status |
| Any genotype + autosomal dominant | Possibly affected — requires clinical confirmation | Investigate — check penetrance, ClinVar review status, and phenotype |
| Compound het (two variants, same gene) | Potentially affected (recessive) | Check if variants are on different alleles (phasing) |

## Limitations

- The reports show ClinVar's review status (`CLNREVSTAT`), but a status such as "criteria provided, single submitter" still rests on one lab. Always check the full ClinVar entry before acting on any result.
- Matching is by allele after both files are split and left-aligned, in the bash script and in the Nextflow module alike, plus the combined allele of two neighbouring records (item 5). A complex variant can still be missed when the caller splits it into three or more records, when its parts are more than 10 bases apart, when their genotypes differ or they sit on different phase sets, or when ClinVar's record normalises to another representation. Two hits in one gene within 10 bases whose genotypes allow one allele, and that matched no combined record, are flagged: check ClinVar for a combined record by hand.
- A hit is only as good as the reads behind it. A flagged hit from the cDNA read check, a hit at a known lab-construct site, or a locus of the cDNA scan may come from contaminating DNA, not from your genome; check its reads before reading anything into it.
- The always-report list holds six sites. A "ref" there is a reference call at that one position, not a clean result for the gene.
- The screen sees only small variants in the VCF. Copy-number losses and gene deletions (for example SMN1 in spinal muscular atrophy) are invisible to it.
- A common allele with a P/LP record is often low-penetrance or misclassified. That holds for a heterozygous hit in a dominant gene as much as for a homozygous one, so read every hit with its gnomAD frequency, and look a hit up by hand when its frequency is unknown. The report's hit-count badge is never green for this reason.
- The call floor reads only GQ and DP. A call above it can still be wrong (in a repeat, a segmental duplication or a region of no-calls), and a call below it can be right; the mark says to look at the reads.
- Results are **research-grade**, not clinical diagnoses. Do not make medical decisions based solely on this output.

## Important Notes
- Most hits will be **heterozygous carriers of recessive conditions** — this is normal and expected
- A typical 30X WGS shows 0-10 pathogenic/likely pathogenic overlaps. The majority are benign carrier states for recessive conditions
- Focus review on: homozygous pathogenic, any autosomal dominant pathogenic, and compound heterozygous variants in the same gene
- The ClinVar pathogenic database must be chr-prefixed to match the BAM/VCF coordinate system (done in step 00)
- ClinVar is updated monthly — re-download periodically to catch newly classified variants
