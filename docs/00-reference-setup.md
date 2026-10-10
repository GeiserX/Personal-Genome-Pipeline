# Step 0: Reference Data Setup

One-time downloads required before running the pipeline. Each heading below gives the size of its download; [Hardware and storage requirements](hardware-requirements.md#shared-reference-data-one-time) adds them up (about 73 GB for a default run, 248 GB with the optional annotation databases).

> **Estimated time:** 1-3 hours depending on internet speed. The two VEP caches (26 GB for step 13, 24 GB for step 17) are the largest downloads of a default run.

## GRCh38 Reference Genome

The foundation for everything. All tools need this. The pipeline uses NCBI's **GRCh38 no-ALT analysis set**: chr1-22, X, Y and M, the unplaced and unlocalized scaffolds and the EBV genome, 195 sequences, with UCSC names (`chr1`, `chrM`). `setup.sh` downloads it, checks it and stores it as `reference/GRCh38_no_alt_analysis_set.fasta`:

```bash
./scripts/setup.sh ${GENOME_DIR}
```

By hand, the same steps:

```bash
export GENOME_DIR=/path/to/your/data
mkdir -p ${GENOME_DIR}/reference
cd ${GENOME_DIR}/reference
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

# Download (~0.9 GB compressed, ~3.2 GB unpacked) and its index
wget -c ${NCBI}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget -c ${NCBI}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
wget ${NCBI}/md5checksums.txt

# Verify both against the md5 NCBI lists for them; each line must say OK
grep -E ' \./GCA_000001405\.15_GRCh38_no_alt_analysis_set\.fna\.(gz|fai)$' md5checksums.txt | md5sum -c -

# Unpack under the name the scripts read
gzip -dc GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GRCh38_no_alt_analysis_set.fasta
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai GRCh38_no_alt_analysis_set.fasta.fai
wc -l GRCh38_no_alt_analysis_set.fasta.fai   # 195
```

### Why the no-ALT analysis set

GRCh38 also has ALT contigs: second copies of regions that vary a lot between people, such as the MHC (the HLA genes), KIR and the CYP2D6 locus. A reference with them (the Broad `hg38` FASTA, 3,366 sequences with ALT, HLA and decoy contigs) only helps an aligner that is run ALT-aware, and none of the aligners here is. With them, a read that matches the primary copy and an ALT copy equally well gets mapping quality 0, and callers ignore it: depth thins at exactly those loci, and depth-based callers can report a deletion that is not there. The [measured loss](realignment.md#how-much-depth-alt-contigs-cost) is in the realignment page.

The no-ALT analysis set has each of those regions once. It is the reference DeepVariant's case studies use, the one the GIAB v4.2.1 benchmark is defined on, and the one PharmCAT normalises against. chr1-22, X, Y and M are the same sequence as in the Broad file, so ClinVar, the VEP cache and the score files stay valid.

Two other references were considered and are not the default:

- **The no-ALT set plus the hs38d1 decoys** (`GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz` in the same NCBI directory, 2,580 sequences). The decoys catch reads from sequence missing in the primary assembly, which removes some false calls for callers less robust than DeepVariant. None of the tools here expects them, and its 2,385 extra contigs are more that CNVpytor and Delly have to be kept away from. It works through `REF_FASTA` (below).
- **GIAB's GRCh38 file that also masks false duplications.** It is the natural next step if a benchmark run shows a gain.

GSTT1 lies on an ALT contig in GRCh38, so pypgx (step 32) cannot call it on the default reference; the step says so and calls the other genes.

Every sample aligned to another reference has to be aligned again. [Realigning after a reference change](realignment.md) says what to redo, what to keep and how to check a BAM.

### The reference path on every page

The step scripts read the reference from `${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta`. The commands on the other pages of these docs write that path as `${REF_FASTA}`, relative to `GENOME_DIR`: on the host it is `${GENOME_DIR}/${REF_FASTA}`, and inside a container that mounts `GENOME_DIR` at `/genome` it is `/genome/${REF_FASTA}`. Set it once in the shell where you paste those commands:

```bash
export REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta
```

The scripts read `REF_FASTA` too, as a path relative to `GENOME_DIR` or as an absolute path inside it, so the same variable points both at another reference. `setup.sh` then also needs `REF_FASTA_URL` (a `.gz` URL is unpacked) and `REF_FASTA_MD5` (an md5, the URL of a checksum file that lists the download, or empty), and `REF_FAI_MD5` the same way for the `.fai` published beside it (empty builds the index with samtools). For the decoy variant:

```bash
export REF_FASTA=reference/GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
export REF_FASTA_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
./scripts/setup.sh ${GENOME_DIR}   # checks both files against NCBI's md5checksums.txt
```

`validate-setup.sh` fails when the reference has ALT or HLA contigs. To keep such a reference on purpose, set `ALLOW_ALT_REFERENCE=true` as well; the check then warns instead.

### Why GRCh38?

This pipeline uses **GRCh38** (also called hg38) exclusively. It's the current standard genome build with:
- Corrected mitochondrial sequence (rCRS)
- Better representation of centromeres and telomeres
- `chr` prefix naming (chr1, chr2, ..., chrX, chrY, chrM)

If your data is on **GRCh37/hg19**, extract FASTQ from BAM and re-align. See [vendor-guide.md](vendor-guide.md#genome-build-grch37-hg19-vs-grch38-hg38).

## ClinVar Database

Updated monthly by NCBI. Contains known pathogenic/benign variant classifications. `setup.sh` downloads it, checks it against NCBI's published md5, builds the two files the steps read from it and records its release date:

| File | What it is |
|---|---|
| `clinvar/clinvar.vcf.gz` (+ `.tbi`) | NCBI's file as published, chromosomes named `1, 2 ... MT` |
| `clinvar/clinvar_chr.vcf.gz` (+ `.tbi`) | the same records with `chr1, chr2 ... chrM`, the names the BAMs and VCFs use |
| `clinvar/clinvar_pathogenic_chr.vcf.gz` (+ `.tbi`) | the Pathogenic and Likely_pathogenic records, read by step 6 |
| `clinvar/RELEASE` | the release date, from the file's `##fileDate` line |

To move to the current release, run:

```bash
./scripts/setup.sh --refresh clinvar ${GENOME_DIR}
```

It downloads the new file under `clinvar/.refresh/`, checks the md5, builds both derived files there, and only then replaces all six files and `RELEASE`. Step 6's normalised copy (`clinvar_pathogenic_chr.norm.vcf.gz`) is removed, so step 6 rebuilds it from the new release. A failed download or build leaves the installed release as it was. If moving the new files in fails part way, the refresh exits with an error and leaves no `RELEASE`, so `validate-setup.sh` reports the date as unknown until a refresh completes. Do not re-download with a plain `wget`: it writes `clinvar.vcf.gz.1` beside the old file, and the derived files stay built from the old one.

`validate-setup.sh` prints the release date and warns when it is more than 35 days old.

> **Tip:** Refresh ClinVar monthly for the latest classifications. ClinVar adds ~1000 new pathogenic variants per month. Rerun step 6 afterwards.

## Small pinned data files

`setup.sh` also installs these small files, each from a fixed commit or release and checked before it is stored (`./scripts/setup.sh --sample-qc-data ${GENOME_DIR}` installs only the last two). A failed download does not stop setup; the next run tries again, and `validate-setup.sh` lists what is missing.

| File | Source | Used by | Without it |
|---|---|---|---|
| `reference/delly_human.hg38.excl.tsv` | Delly's GRCh38 exclude map at a pinned commit (sha256 checked) | step 19 (`delly sr -x`) | Delly runs without it: slower, and with calls in centromeres, telomeres and the unplaced scaffolds |
| `reference/cytoBand.hg38.txt` | UCSC's GRCh38 chromosome bands, chr1-22, X and Y (sha256 checked) | step 10 (`telomerehunter -b`) | TelomereHunter falls back to its hg19 bands |
| `hla/IPD-IMGT-HLA_<release>/hla.dat` | IPD-IMGT/HLA release `HLA_DB_RELEASE` (3.65.0) from the IMGTHLA repository (md5 checked) | step 8 | step 8 is skipped |
| `reference/gencode.v50.basic.genes.gtf` | the gene lines of GENCODE 50's basic annotation (md5 checked) | step 8 (gene positions for T1K) | step 8 is skipped |
| `reference/somalier/sites.hg38.vcf.gz` | somalier's GRCh38 sites, 17,766 SNPs (sha256 checked) | step 33 (`somalier extract`), Nextflow `--somalier_sites` | step 33 stops and says to install it |
| `reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.{UD,mu,bed}` | VerifyBamID2's 1000 Genomes panel of 100,000 markers, from the VerifyBamID v2.0.3 release (sha256 checked) | step 33 (`verifybamid2`), Nextflow `--verifybamid2_panel` | step 33 stops and says to install it |

`setup.sh` also writes the reference's sequence dictionary (`GRCh38_no_alt_analysis_set.dict`), which GATK, Picard and `chip-to-vcf.sh` need.

## AnnotSV Annotations (~5 GB)

Required for step 5. `setup.sh` downloads them (about 5 GB, about 20 GB unpacked) into `${GENOME_DIR}/annotsv_annotations/`; its source server is slow, so this can take an hour or two, and an interrupted download resumes when `setup.sh` runs again.

## VEP Cache (~26 GB)

Ensembl Variant Effect Predictor annotation database. Required for step 13 and the pipeline's VEP. Without it, a run skips VEP, vcfanno, the clinical filter and slivar. Install it before the first run:

```bash
./scripts/setup.sh --vep-cache ${GENOME_DIR}
```

It downloads the release-116 cache, checks it against Ensembl's `CHECKSUMS` file, unpacks it into a temporary folder and moves it to `${GENOME_DIR}/vep_cache/homo_sapiens/116_GRCh38/` only when it is complete. The tarball (~26 GB) and the unpacked cache (~30 GB) are both on disk until the install ends, then the tarball is deleted. An interrupted download resumes when you run the command again.

By hand, the same steps:

```bash
mkdir -p ${GENOME_DIR}/vep_cache/tmp
cd ${GENOME_DIR}/vep_cache/tmp

# Download (manual wget is more reliable than VEP INSTALL.pl)
# The -c flag enables resume if the download is interrupted
wget -c https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz

# Extract to parent directory (~30 GB extracted)
cd ${GENOME_DIR}/vep_cache
tar xzf tmp/homo_sapiens_vep_116_GRCh38.tar.gz
# Creates: ${GENOME_DIR}/vep_cache/homo_sapiens/116_GRCh38/

# Optional: delete the tarball to save 26 GB
# rm tmp/homo_sapiens_vep_116_GRCh38.tar.gz
```

> **Warning:** The VEP `INSTALL.pl` script downloads to a temporary directory that may lack write permissions inside Docker. Always download manually with `wget -c`. See [lessons-learned.md](lessons-learned.md) for details.

## PCGR/CPSR Ref Data Bundle (~7 GB)

Required for step 17 (CPSR cancer predisposition screening). Includes ClinVar, gnomAD, CancerMine, and other databases. PCGR 2.x uses a separate, smaller ref data bundle — VEP cache is mounted independently.

> **Important:** PCGR 2.3.2 bundles VEP 115, which requires a **release-115** cache — different from the release-116 cache used by step 13 above. See the next section for the VEP 115 download. The bundle date (`PCGR_DATA_BUNDLE`) and the cache release (`PCGR_VEP_CACHE_RELEASE`) are in `versions.env` and move with `PCGR_IMAGE`.

```bash
mkdir -p ${GENOME_DIR}/pcgr_data
cd ${GENOME_DIR}/pcgr_data

# Download (~7 GB)
wget -c https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20260620.grch38.tgz

# Extract and organize into version-stamped directory
tar xzf pcgr_ref_data.20260620.grch38.tgz
mkdir -p 20260620 && mv data/ 20260620/
# Creates: ${GENOME_DIR}/pcgr_data/20260620/data/

# Optional: delete the tarball to save 7 GB
# rm pcgr_ref_data.20260620.grch38.tgz
```

## VEP 115 Cache for CPSR (~24 GB)

PCGR 2.3.2 (step 17) bundles VEP 115 internally, which needs the **release-115** cache. This is separate from the release-116 cache used by step 13. Both coexist in the same `vep_cache/` directory under different subdirectories (`116_GRCh38/` and `115_GRCh38/`). A `113_GRCh38/` cache left from PCGR 2.2.5 is no longer read and can be deleted.

```bash
mkdir -p ${GENOME_DIR}/vep_cache/tmp
cd ${GENOME_DIR}/vep_cache/tmp

# Download VEP 115 cache (~24 GB)
wget -c https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz

# Extract alongside the existing release-116 cache
cd ${GENOME_DIR}/vep_cache
tar xzf tmp/homo_sapiens_vep_115_GRCh38.tar.gz
# Creates: ${GENOME_DIR}/vep_cache/homo_sapiens/115_GRCh38/
```

> If you only run step 13 (VEP annotation) and skip step 17 (CPSR), you only need the release-116 cache. If you only run step 17, you only need release-115.

## T1K HLA Reference (Optional)

Only needed for step 8 (HLA typing). Step 8 builds its T1K index itself, the first time it runs (a few minutes), from the `hla.dat` and GENCODE gene lines `setup.sh` installs (see "Small pinned data files" above). The index goes to `t1k_idx/t1k-<T1K version>_imgt-<release>_gencode-<release>/`, so a new T1K image or a new database release builds a new index instead of reusing the old one. Step 8 writes the release it typed against to `hla_t1k/database_release.txt`.

To type against another IPD-IMGT/HLA release, set `HLA_DB_RELEASE` (for example `HLA_DB_RELEASE=3.64.0`) for both `setup.sh` and step 8.

The coordinate file takes each HLA gene's GRCh38 position from the GENCODE annotation (`t1k-build.pl -g`, as T1K's README describes). Built from the FASTA or its `.fai` instead, every gene gets `-1 -1` coordinates and T1K extracts no reads; step 8 stops when a typed gene has no coordinates.

> **Note:** HLA typing from WGS is challenging. For clinical HLA typing, dedicated lab assays are more reliable. See [lessons-learned.md](lessons-learned.md#t1k-coordinate-file-with-wrong-values).

## CNVpytor GC/Mask Resources (Optional, for Step 18)

Only needed for step 18 (CNVpytor CNV calling). The 1.3.2 biocontainer ships **without** the reference GC/mask files and its built-in `-download` is broken, so download the pinned v1.3.2 resources once and mount them at run time.

```bash
mkdir -p ${GENOME_DIR}/reference/cnvpytor
BASE=https://github.com/abyzovlab/CNVpytor/raw/v1.3.2/cnvpytor/data

# gc_hg38/mask_hg38 are the files actually used; the other genomes' files only
# need to exist to satisfy CNVpytor's (over-eager) global resource check.
for f in gc_hg38.pytor mask_hg38.pytor gc_hg19.pytor mask_hg19.pytor \
         gc_chm13v2.0.pytor gc_chm13v1.1.pytor gc_kn99.pytor; do
  curl -fsSL -o ${GENOME_DIR}/reference/cnvpytor/$f "$BASE/$f"
done
```

Pinning to the `v1.3.2` git tag (not `master`) keeps runs reproducible; the files total ~90 MB. `scripts/18-cnvpytor.sh` bind-mounts this directory onto the container's package data path, so no network access is needed at run time.

## Somatic Calling Resources (Optional, for Step 29)

Only needed if you plan to run step 29 (somatic variant calling with Mutect2 tumor-only mode). These resources significantly reduce false positives. See [29-mutect2-somatic.md](29-mutect2-somatic.md) for details.

```bash
mkdir -p ${GENOME_DIR}/somatic

# gnomAD AF-only VCF (~3 GB) — germline allele frequencies for filtering
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
  -O ${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi \
  -O ${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz.tbi

# Panel of Normals (~17 MB) — recurrent technical artifacts from 1000 Genomes
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz \
  -O ${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi \
  -O ${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz.tbi
```

## Annotation Databases (Optional, for Steps 30-31)

Deep pathogenicity scoring and variant prioritization. These databases power vcfanno (step 30) and slivar (step 31). All are optional — step 30 gracefully skips any missing track.

> **Total download:** ~175 GB, the sum of the sections below (SpliceAI ~91 GB, CADD ~83 GB). If disk space is tight, start with REVEL + AlphaMissense (~1.2 GB combined) — they provide the highest value per byte for missense variant interpretation.

### CADD v1.7 Pre-scored (SNVs + Indels) — ~83 GB

The most widely used deleteriousness metric. Scores all possible SNVs genome-wide plus gnomAD indels.

> **License:** Free for non-commercial/academic use. Commercial use requires a license from the University of Washington. The pipeline does not redistribute CADD data — users download directly from the source.

```bash
mkdir -p ${GENOME_DIR}/annotations
cd ${GENOME_DIR}/annotations

# Whole-genome SNV scores (~81.5 GB, pre-indexed)
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi

# gnomAD v4.0 indel scores (~1.2 GB, pre-indexed)
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi
```

> **Note:** CADD TSVs use chromosome names without `chr` prefix (1, 2, 3...). vcfanno handles this via column mapping in the TOML config — no manual renaming needed.

### SpliceAI Pre-scored — ~91 GB

Deep learning splice-site variant prediction. Catches pathogenic intronic variants that VEP's rule-based splice prediction misses.

The pipeline uses the **masked** score files. SpliceAI's authors recommend the masked scores for variant interpretation and the raw scores for alternative-splicing research. In the masked files, a gain at an annotated splice site and a loss at a site that is not an annotated splice site are set to 0, because neither changes how a variant is interpreted.

> **License:** The precomputed SpliceAI scores are free for academic and not-for-profit use only; any other use requires a commercial license from Illumina. They are not Apache 2.0. The pipeline does not redistribute them.

```bash
cd ${GENOME_DIR}/annotations

# SNV splice scores (~27 GB)
wget -c https://download.molgeniscloud.org/downloads/vip/resources/GRCh38/spliceai_scores.masked.snv.hg38.vcf.gz
wget -c https://download.molgeniscloud.org/downloads/vip/resources/GRCh38/spliceai_scores.masked.snv.hg38.vcf.gz.tbi

# Indel splice scores (~64 GB)
wget -c https://download.molgeniscloud.org/downloads/vip/resources/GRCh38/spliceai_scores.masked.indel.hg38.vcf.gz
wget -c https://download.molgeniscloud.org/downloads/vip/resources/GRCh38/spliceai_scores.masked.indel.hg38.vcf.gz.tbi
```

Step 30 and `validate-setup.sh` find the masked files under these names. They also accept the raw files (`spliceai_scores.raw.*`) if you already have them; when both sets are present, step 30 uses the raw ones.

> **Alternative source:** Illumina BaseSpace at `https://basespace.illumina.com/s/otSPW8hnhaZR` (requires free account).

### REVEL v1.3 — ~0.6 GB

Ensemble pathogenicity scoring for missense variants, combining 13 individual tools. Recommended by ClinGen for missense variant classification.

```bash
source versions.env   # from the repository root
cd ${GENOME_DIR}/annotations

# Download and prepare for vcfanno
wget -c https://zenodo.org/record/7072866/files/revel-v1.3_all_chromosomes.zip
unzip revel-v1.3_all_chromosomes.zip

# Convert to tabix-indexed TSV for GRCh38
# Extract GRCh38 columns, add chr prefix, sort, bgzip, index.
# The GATK image is used because it ships bgzip and tabix; the bcftools image does not.
# The source file has nine comma-separated columns:
#   chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
# The table keeps five of them, so it gets its own five-name header.
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "${GATK_IMAGE}" \
  bash -c '
    set -euo pipefail
    cd /genome/annotations
    printf "#chr\tpos\tref\talt\tREVEL\n" > revel_grch38.tsv
    # Skip the header, keep rows with a GRCh38 position (col 3), add the chr prefix, sort
    tail -n+2 revel_with_transcript_ids | \
      awk -F"," -v OFS="\t" "\$3 != \".\" {print \"chr\"\$1, \$3, \$4, \$5, \$8}" | \
      sort -T /genome/annotations -k1,1V -k2,2n >> revel_grch38.tsv
    bgzip -f revel_grch38.tsv
    tabix -f -s 1 -b 2 -e 2 revel_grch38.tsv.gz
    # Clean up intermediates (only reached if every command above succeeded)
    rm -f revel_with_transcript_ids revel-v1.3_all_chromosomes.zip
  '

# Check: five columns under a five-name header
gzip -dc revel_grch38.tsv.gz | head -2
```

> **Thresholds:** ClinGen's calibration (Pejaver et al. 2022) gives REVEL >= 0.644 as PP3_Supporting, >= 0.773 as PP3_Moderate and >= 0.932 as PP3_Strong for missense variants. There is no PP3 Very Strong level. See [interpreting-results.md](interpreting-results.md#revel-rare-exome-variant-ensemble-learner) for the BP4 (benign) levels.

### AlphaMissense — ~613 MB

DeepMind's protein-structure-informed missense classifier. Complements REVEL with structural context. Covers all ~71 million possible human missense variants.

> **License:** CC BY-NC-SA 4.0 (non-commercial, share-alike). The pipeline does not redistribute these scores.

```bash
source versions.env   # from the repository root
cd ${GENOME_DIR}/annotations

# Download pre-scored GRCh38 predictions
wget -c https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz

# Index for vcfanno (skip header lines starting with #). The GATK image has
# tabix; the bcftools image does not.
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "${GATK_IMAGE}" \
  tabix -s 1 -b 2 -e 2 -S 1 /genome/annotations/AlphaMissense_hg38.tsv.gz
```

> **Thresholds:** am_pathogenicity < 0.34 = likely benign, > 0.564 = likely pathogenic. These are AlphaMissense's own class boundaries, not ACMG evidence levels.

### gnomAD v4.1 Constraint Metrics — ~91 MB

Per-gene loss-of-function intolerance scores. Essential for interpreting novel variants in constrained genes.

```bash
cd ${GENOME_DIR}/annotations

# Gene-level constraint table
wget -c -O gnomad_v4.1_constraint.tsv \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
```

> **Key metrics:** LOEUF (oe_lof_upper) < 0.35 = constrained for LoF; pLI >= 0.9 = intolerant to LoF; missense Z > 3.09 = constrained for missense.

## Docker Images — Pre-Pull All

`setup.sh` pulls every image in `versions.env` except the lines marked `# optional`, and `scripts/setup.sh --pull-only` pulls the same list without the rest of setup. A failed pull is tried three times (`FETCH_TRIES` and `FETCH_WAIT` change that), and Docker's own message is printed: `toomanyrequests` is Docker Hub's rate limit, `manifest unknown` a tag that does not exist. `--parascopy-data` and `--yleaf-data` pull their step's optional image along with its data. To pull every image in advance, the optional ones included, from the repository root in bash:

```bash
source versions.env
for var in $(grep -oE '^[A-Z0-9_]+_IMAGE=' versions.env | tr -d '='); do
  docker pull "${!var}"
done
```

## Alternative Callers (Optional, for Benchmarking)

Only needed if you plan to run alternative variant callers. See [benchmarking.md](benchmarking.md).

The alternative aligners and callers and hap.py (`benchmark-variants.sh`) are the lines marked `# optional` in `versions.env`. `setup.sh` skips them, and each one is pulled the first time its script runs. [Image versions](versions.md) lists them under "Alternative aligners and callers (optional)".

### GATK Sequence Dictionary

GATK HaplotypeCaller requires a `.dict` file alongside the reference FASTA. `setup.sh` creates it; to create it by hand:

```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see "The reference path on every page" above
docker run --rm --user root \
  -v ${GENOME_DIR}:/genome \
  "${GATK_IMAGE}" \
  gatk CreateSequenceDictionary \
    -R /genome/${REF_FASTA}
```

### BWA-MEM2 Index

BWA-MEM2 requires its own index files (different from minimap2's `.mmi`). Building them needs a lot of memory: the bwa-mem2 README states 28 GB per Gbp of reference, which is **about 90 GB of RAM** for the 3.2 Gbp GRCh38 FASTA. The finished index is about 10 GB on disk, and aligning with it needs about 10 GB of RAM, so only the one-time build is the problem.

The practical route is to build the index once on a machine (or a rented cloud instance) with at least 96 GB of RAM, then copy the five index files next to the FASTA on your own machine. With less memory the build is killed (exit code 137).

```bash
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see "The reference path on every page" above
source versions.env   # from the repository root: BWAMEM2_IMAGE, the image step 02a uses
docker run --rm --user root \
  --cpus 8 --memory 96g \
  -v ${GENOME_DIR}:/genome \
  "${BWAMEM2_IMAGE}" \
  bwa-mem2 index "/genome/${REF_FASTA}"
# Creates: .0123, .amb, .ann, .bwt.2bit.64, .pac alongside the FASTA
```

GRIDSS (step 4b) needs the classic BWA index instead (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`); the two are not interchangeable. See [04b-gridss.md](04b-gridss.md) for how to build it.

### GIAB Truth Set (for hap.py Benchmarking)

Download a GIAB truth set for benchmarking variant callers. **HG002** (Ashkenazi Jewish male) is preferred because its truth set covers more difficult genomic regions. HG001 (NA12878) is an alternative used in the [quick test](quick-test.md).

**Important:** Truth set benchmarking only works when the query VCF comes from the **same biological sample** as the truth set. You must sequence HG002 (or HG001) DNA, not your own sample. See [benchmarking.md](benchmarking.md) for details.

```bash
mkdir -p ${GENOME_DIR}/giab
cd ${GENOME_DIR}/giab

# HG002 truth set (recommended, GRCh38 v4.2.1)
# Pinned to the NISTv4.2.1 directory: GIAB's latest/ directory moves to each new
# release (it now holds v5.0q), so latest/ URLs for v4.2.1 files return 404.
wget -c https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -c https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget -c https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# Alternative: HG001/NA12878 (used in quick-test.md)
# wget -c https://giab.s3.amazonaws.com/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
# wget -c https://giab.s3.amazonaws.com/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
# wget -c https://giab.s3.amazonaws.com/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

**Total Docker image size:** ~10-15 GB (compressed, after layer deduplication). Alternative tools add ~3-5 GB.

## Disk Space Summary

[Hardware and storage requirements](hardware-requirements.md#shared-reference-data-one-time) has the one table of every download above, with the totals.

> **Tip:** If disk space is tight, you can skip the VEP cache (step 13) and PCGR bundle (step 17) initially. The core pipeline (steps 2-3-6-7) only needs the reference FASTA and ClinVar (~3.5 GB unpacked).

## Verifying Your Setup

After all downloads, verify everything is in place:

```bash
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see "The reference path on every page" above
echo "Checking reference setup..."
[ -f "${GENOME_DIR}/${REF_FASTA}" ] && echo "  GRCh38 FASTA: OK" || echo "  GRCh38 FASTA: MISSING"
[ -f "${GENOME_DIR}/${REF_FASTA}.fai" ] && echo "  FASTA index: OK" || echo "  FASTA index: MISSING"
[ -f "${GENOME_DIR}/clinvar/clinvar_chr.vcf.gz" ] && echo "  ClinVar (chr): OK" || echo "  ClinVar: MISSING"
[ -d "${GENOME_DIR}/vep_cache/homo_sapiens/116_GRCh38" ] && echo "  VEP cache: OK" || echo "  VEP cache: MISSING"
[ -d "${GENOME_DIR}/pcgr_data/20260620/data" ] && echo "  PCGR data: OK" || echo "  PCGR data: MISSING"
[ -d "${GENOME_DIR}/vep_cache/homo_sapiens/115_GRCh38" ] && echo "  VEP 115 cache (CPSR): OK" || echo "  VEP 115 cache (CPSR): MISSING"

echo "Annotation databases (optional, for steps 30-31):"
for DB_PAIR in \
  "whole_genome_SNVs.tsv.gz:CADD SNVs" \
  "gnomad.genomes.r4.0.indel.tsv.gz:CADD indels" \
  "spliceai_scores.masked.snv.hg38.vcf.gz:SpliceAI SNVs" \
  "spliceai_scores.masked.indel.hg38.vcf.gz:SpliceAI indels" \
  "revel_grch38.tsv.gz:REVEL" \
  "AlphaMissense_hg38.tsv.gz:AlphaMissense" \
  "gnomad_v4.1_constraint.tsv:gnomAD constraint"; do
  DB_FILE="${DB_PAIR%%:*}"
  DB_NAME="${DB_PAIR#*:}"
  if [ -f "${GENOME_DIR}/annotations/${DB_FILE}" ]; then
    # Check for .tbi index (not needed for plain TSV files)
    if [[ "$DB_FILE" == *.vcf.gz ]] || [[ "$DB_FILE" == *.tsv.gz && "$DB_FILE" != *constraint* ]]; then
      [ -f "${GENOME_DIR}/annotations/${DB_FILE}.tbi" ] && echo "  ${DB_NAME}: OK (+index)" || echo "  ${DB_NAME}: OK (WARNING: .tbi index missing)"
    else
      echo "  ${DB_NAME}: OK"
    fi
  else
    echo "  ${DB_NAME}: not downloaded"
  fi
done
echo "Done."
```
