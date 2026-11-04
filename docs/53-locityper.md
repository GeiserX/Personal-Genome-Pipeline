# Step 53: Hard Loci Second Opinion (Locityper)

## What This Does

Genotypes five of the hardest polymorphic regions of the genome with Locityper: HLA-DRB1, HLA-DQB1, the KIR cluster, CYP2D6 and FCGR. For each one it names the two haplotypes of a pangenome panel that best explain the sample's reads there, with a quality and a pass or fail filter. Opt-in.

It is a second opinion. The report shows it beside the call of the step that owns the gene, never in its place.

## Why

At these genes a wrong call from short reads is likely, and the pipeline has one caller or none at each: T1K (step 08) for HLA-DRB1, HLA-DQB1 and KIR, pypgx (step 32) and Cyrius (step 21) for CYP2D6, nothing for FCGR. With one caller a wrong call cannot be seen.

Locityper (Prodanov et al., Nat Genet 2025, PMID 41107551, doi 10.1038/s41588-025-02362-4) works another way. It does not call variants or alleles from alignments. It takes the haplotypes of many assembled genomes at a locus and finds the pair whose k-mers, read depth and read alignments fit the sample's reads best. A caller with a different method gives an independent check where the primary caller is weakest.

## Tool

- Locityper, the release `LOCITYPER_IMAGE` pins (1.7.14 today; MIT licence, https://github.com/tprodanov/locityper): `locityper preproc` and `locityper genotype` on the BAM
- `bin/locityper_report.py` (standard library Python, run with the image's Python) for the table

## Docker Image

- `LOCITYPER_IMAGE`: the Bioconda build, `quay.io/biocontainers/locityper`. It carries Locityper, jellyfish, samtools, minimap2, strobealign, tabix and Python. 1.7.15 exists; its only change is in `locityper align`, which this step does not run, so the pin stays at 1.7.14.

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. The line is marked optional, so `setup.sh` pulls it only with `--locityper-data`.

## Database

```bash
./scripts/setup.sh --locityper-data ${GENOME_DIR}
# Creates: ${GENOME_DIR}/reference/locityper-hprc-v2.1-mc-grch38/ (counts.jf, db/, loci.bed, DONE)
```

The panel is the HPRC release 2 Minigraph-Cactus pangenome against GRCh38, as a VCF (`LOCITYPER_PANEL_URL`, `LOCITYPER_PANEL_VERSION` in `versions.env`). It holds the two haplotypes of each of its 232 samples. The full VCF is 6.7 GB; `setup.sh` reads only the loci and 250 kb either side, over https, through the VCF's `.tbi`, whose sha256 `LOCITYPER_PANEL_TBI_SHA256` pins. Nothing of it is committed.

Then it follows Locityper's documentation (https://locityper.vercel.app, "Preparing targets"): jellyfish counts the 25-mers of `REF_FASTA`, `locityper target` cuts the panel's haplotypes at each locus, and `locityper augment` aligns them to each other, which the weighted distance needs. A new reference, image, panel or locus list builds again. `LOCITYPER_DATA` puts the database elsewhere and `LOCITYPER_LOCI` replaces the locus list.

## The loci

`assets/locityper_loci.bed` lists them. The regions are those of the Locityper authors' own GRCh38 database of 326 loci (Zenodo record 15879519, https://zenodo.org/records/15879519), except FCGR, where their FCGR3A and FCGR2B loci overlap and are joined.

| Locus | GRCh38 region | Primary caller here |
|---|---|---|
| FCGR | chr1:161,506,762-161,681,047 | none |
| HLA-DRB1 | chr6:32,484,348-32,593,538 (the authors' HLA-DRB5 locus, which holds DRB1) | T1K, step 08 |
| HLA-DQB1 | chr6:32,652,604-32,669,077 | T1K, step 08 |
| KIR | chr19:54,723,154-54,868,521 (the whole cluster) | T1K with `KIR=true`, step 08 |
| CYP2D6 | chr22:42,122,285-42,152,960 | pypgx (step 32) and Cyrius (step 21) |

We left C4, AMY1 and LPA out. They are not in the authors' database, and on the fixture's HG002 reads a C4 call picked a haplotype with one RCCX module where HG002 has two (QV 4.6 to the true haplotype), at a quality of 161. A wrong copy number with a high quality is what a second opinion must not give.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai`, aligned to `REF_FASTA`, the reference the database was built from. Locityper reads the BAM's contigs by their order in the reference.

## Command

```bash
./scripts/53-locityper.sh your_name
```

With `run-all.sh`: `TOOLS=...,locityper`. It is opt-in, so a run lists step 53 only when `TOOLS` names it, and it is skipped for a sample with no BAM or without the database. In the Nextflow pipeline: `locityper` in `--tools` with `--locityper_data` (process `LOCITYPER`, for each sample with a BAM).

`locityper preproc` measures the sample's read depth, insert size and error profile on a background region, by default Locityper's own for GRCh38 (chr17:72,062,001-76,562,000). `LOCITYPER_BG_REGION` (bash) or `--locityper_bg_region` (Nextflow) picks another, for a reference without that region.

## How the table is made

`bin/locityper_report.py` reads `gt/loci/<locus>/res.json.gz` and writes one row per locus. A locus Locityper did not finish reads "not called", with Locityper's exit code when it stopped.

The filter is the one the authors used for their HLA and KIR benchmark (Zenodo 15879519, README): a call passes with no Locityper warning, a weighted distance of 30 or less, and fewer than 1000 or fewer than 20% unexplained reads. The weighted distance is the edit distance from the call to the other likely genotypes, weighted by their probability: a high value means the reads fit very different genotypes about as well.

Locityper names haplotypes of panel people, such as `HG00733.1`, not alleles. The two callers are shown side by side and not compared: an allele name from T1K and a panel haplotype are not the same kind of answer.

## Output

```
${SAMPLE}/locityper/
├── ${SAMPLE}_locityper.tsv     # one row per locus: the two haplotypes, quality, reads, filter
├── ${SAMPLE}_locityper.log     # Locityper's log
└── work/
    ├── preproc/                # the sample's depth, insert size and error profile
    └── gt/loci/<locus>/res.json.gz
```

Columns: `locus`, `region`, `status`, `genotype`, `haplotype_1`, `haplotype_2`, `quality` (Phred), `total_reads`, `unexplained_reads`, `weighted_distance`, `warnings`, `filter`, `reason`.

Both reports carry a section "Locityper Second Opinion at Hard Loci" in the medicines chapter: per locus Locityper's call with its filter, and the primary caller's call.

## How it was checked

The image test (`tests/smoke/locityper.sh`, the `LOCITYPER_IMAGE` row of `tests/smoke/commands.tsv`) runs the database build and the step on the fixture's GIAB HG002 reads at HLA-DRB1, HLA-DQB1 and CYP2D6. The panel holds HG002's own two assembled haplotypes, which are the truth. Both called haplotypes at each locus must match them at QV 20 or more (one difference per 100 bases or fewer), aligned with `locityper align`.

The red control takes out of the panel HG002's haplotypes and every other haplotype within QV 20 of them, rebuilds the database and genotypes again; the same check must then fail at every locus. Leaving HG002 out with `--leave-out` is not enough: other people in the panel carry near copies of HG002's haplotypes, and the call stays close.

On the fixture, run under emulation before the image test existed: CYP2D6 and HLA-DQB1 came back as HG002's own two haplotypes, and HLA-DRB1 as one of them plus a haplotype one base in 98 kb from the other (QV 49.9). HLA-DQB1 and HLA-DRB1 had low qualities (7.3 and 3.1): near-identical haplotypes compete there.

## Runtime

The fixture's three loci took about 5 minutes to build and 8 minutes to augment on an Apple Silicon Mac under x86 emulation. We have not measured a whole genome yet.

## Notes

- A second opinion, not a diagnosis, and not used by PharmCAT or any other step.
- A fail does not make the primary call wrong, and a pass does not make it right: look at both and at the reads.
- The panel is 232 people. A haplotype not in it can only be called as its nearest panel haplotype.

## Sources

- Prodanov T, et al. Nat Genet 2025. PMID 41107551, doi 10.1038/s41588-025-02362-4 (read 2026-11-03)
- Locityper documentation, https://locityper.vercel.app (read 2026-11-03)
- Locityper authors' databases and benchmark filter, Zenodo record 15879519, https://zenodo.org/records/15879519 (read 2026-11-03)
- HPRC release 2 Minigraph-Cactus GRCh38 VCF, https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/ (read 2026-11-03)
