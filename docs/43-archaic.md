# Step 43: Archaic Ancestry (Neanderthal and Denisovan)

## What This Does

Estimates how much of the genome came from Neanderthals and from Denisovans, and lists the stretches where it did. IBDmix compares the genome with three high-coverage archaic genomes, the Altai and Vindija 33.19 Neanderthals and the Altai Denisovan, and calls each stretch of 50 kb or more that matches one of them better than modern allele frequencies explain. Opt-in: it runs once `setup.sh --archaic` has installed its table.

The report gives the Neanderthal and the Denisovan share, the number and total length of the segments, and the ten longest with the protein-coding genes inside them. A TSV holds every segment.

## Why

Most people outside Africa carry about 2% Neanderthal DNA, and people with ancestry from Oceania and parts of Asia carry Denisovan DNA too. Consumer DNA tests give a single percentage with no method behind it. This step uses a published method on the whole genome and says where each segment is.

## Tool

- [IBDmix](https://github.com/PrincetonUniversity/IBDmix) 1.0.1 (Chen et al. 2020, Cell 180:677, doi:10.1016/j.cell.2020.01.012)
- `bin/archaic_report.py` (standard library Python): IBDmix's genotype files, then the segments and the summary

IBDmix needs no unadmixed reference population. It scores each site from the modern allele frequency and the two genotypes, so it finds archaic segments in Africans too, where methods that use Africans as the unadmixed reference cannot.

### One genome and the stand-ins

IBDmix takes the modern allele frequency of each site from the modern genomes in its input. With one genome the frequency is that person's own genotype, and IBDmix drops every site below its minimum allele count, so nothing scores. The step therefore adds 100 stand-in people whose genotypes carry the 1000 Genomes phase 3 allele count of the chosen population at each site (`ARCHAIC_POP`, default `ALL`; `AFR`, `AMR`, `EAS`, `EUR` or `SAS`). IBDmix scores each modern column from the site's frequency and that column's genotype only, so your column gets the score it would get inside a population sample with that frequency. The stand-ins' own rows are dropped.

### The calls

- IBDmix with the options of Chen et al.: LOD 4 (`-d 4`), minor allele count 1 (`-m 1`), archaic error 0.01 (`-a`), modern error 0.002 (`-e`), `-c 2`
- Segments with LOD 4 or more and 50 kb or longer, the cutoffs of the paper
- **Neanderthal**: the union of the segments against Altai and against Vindija. **Denisovan**: the segments against the Denisovan that overlap no Neanderthal segment, because the Neanderthal and Denisovan genomes share much of their sequence.
- The share is the segment length over the analysed span: from the first to the last site of the table on each autosome the VCF has calls on.

## Docker Image

- `IBDMIX_IMAGE` for IBDmix, `PYTHON_IMAGE` for the genotype files and the report
- `setup.sh --archaic` also uses `BCFTOOLS_IMAGE` and `PICARD_IMAGE` (the liftover)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Data

`scripts/setup.sh --archaic <genome_dir>` builds `reference/archaic/snpad2016-1kgp3-v1/`, one `chrN.tsv.gz` per autosome on GRCh38: position, REF, ALT, the ALT count of each archaic genome (9 when it was not called there) and the 1000 Genomes ALT frequency of all and of each continental group. The rows are the 1000 Genomes biallelic SNVs where an archaic genome was called, plus every site where one of the three carries a non-reference allele. Those are the sites IBDmix's own `generate_gt` keeps from the same files.

| Source | What | Pin |
|---|---|---|
| MPI-EVA, https://cdna.eva.mpg.de/neandertal/Vindija/VCF/ | the snpAD genotypes of Altai, Vindija33.19 and Denisova, GRCh37, `chrN_mq25_mapab100.vcf.gz` | each file against the folder's `md5sums.txt`, whose sha256 `versions.env` pins |
| 1000 Genomes phase 3 | `ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz` (1.5 GB) | md5 from the release's `current.tree` |
| UCSC | `hg19ToHg38.over.chain.gz` | sha256 |

MPI-EVA publishes md5 lists, not sha256 sums, for 66 files of about 2.5 GB each. Pinning the sha256 of the three lists pins every file's md5 without downloading 165 GB to hash it.

The table is built on GRCh37 and lifted once to GRCh38 with Picard `LiftoverVcf`. A site that does not lift, or lifts to another chromosome, is dropped. The archaic files are downloaded one chromosome at a time (about 165 GB in all, under 13 GB at once) and deleted once the chromosome is built, so an interrupted install resumes at the chromosome it stopped in. `ARCHAIC_REGIONS` (a GRCh37 BED, chromosomes named 1 to 22) builds only those regions.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.vcf.gz`: one sample. A site with no record is read as the reference allele, a record whose FILTER is not PASS as unknown, and a site where you carry another ALT than the table's is left out.
- The GENCODE gene lines `setup.sh` installs, for the genes of each segment (without them the segments are listed without genes)

## Command

```bash
./scripts/setup.sh --archaic /path/to/genome_dir
./scripts/43-archaic.sh your_name
```

With `run-all.sh`: step 43 runs once the table is installed (`ARCHAIC_DIR=none` skips it, `ARCHAIC_POP` picks the frequencies). In the Nextflow pipeline: `archaic` in `--tools` with `--archaic_data <the table folder>` and `--archaic_pop`; `--hla_genes` gives the genes (processes `ARCHAIC_GENOTYPES`, `ARCHAIC_IBDMIX`, `ARCHAIC_REPORT`).

## Output

| File | Contents |
|---|---|
| `archaic/${SAMPLE}_archaic.tsv` | key and value: population, analysed span, Neanderthal and Denisovan segment count, total bp, percent and longest segment, `longest_1` to `longest_10` |
| `archaic/${SAMPLE}_archaic_segments.tsv` | every segment: `chrom`, `start`, `end`, `length_bp`, `ancestry` (Neanderthal or Denisovan), `matched` (the archaic genomes it matches), `lod`, `genes` |
| `archaic/${SAMPLE}_archaic.txt` | the same in plain text, with the caveat |
| `archaic/work/ibdmix/` | IBDmix's own output per chromosome and archaic genome (bash step) |

Both reports show the shares, the analysed span and the ten longest segments in their "Archaic Ancestry" section.

## Runtime

Not measured on a whole genome. IBDmix runs once per chromosome and archaic genome, `THREADS` chromosomes at a time. The install is bound by the download.

## Notes

- **An estimate, not a measurement.** The share moves with the cutoffs (LOD 4, 50 kb) and with the frequencies used.
- Part of what IBDmix finds in anyone is old ancestry shared with the archaic genomes, not gene flow from them. Chen et al. find Neanderthal sequence in every population they studied, Africans included.
- A segment counts when either of your two copies carries it. The share is of the analysed span, not of the whole genome.
- The table has autosomes only: chrX and chrY are not assessed.
- Sites the VCF has no record for are read as the reference allele. A VCF that holds only variants (no reference calls) cannot tell a reference call from a gap in coverage.
- On the e2e fixture the case builds the table for GRCh37 chr22:25-35 Mb. It checks HG002's GIAB v4.2.1 genotypes there, a copy of HG002 with one chromosome over chr22:28-29 Mb taken from Vindija (a Neanderthal segment must cover it), and a person built from the table who carries only the majority allele (no segment).

## References

- Chen L, Wolf AB, Fu W, Li L, Akey JM. Identifying and Interpreting Apparent Neanderthal Ancestry in African Individuals. Cell 180:677-687 (2020). doi:10.1016/j.cell.2020.01.012
- Prüfer K, et al. The complete genome sequence of a Neanderthal from the Altai Mountains. Nature 505:43-49 (2014). doi:10.1038/nature12886
- Prüfer K, et al. A high-coverage Neandertal genome from Vindija Cave in Croatia. Science 358:655-658 (2017). doi:10.1126/science.aao1887
- Meyer M, et al. A high-coverage genome sequence from an archaic Denisovan individual. Science 338:222-226 (2012). doi:10.1126/science.1224344
