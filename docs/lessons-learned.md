# Lessons Learned

Every failure encountered during pipeline development (Mar 2026), documented so they don't happen again.

## Docker Image Issues

### AnnotSV: Official image doesn't exist
- **Failed:** `bioinfochrustrasbourg/annotsv:3.4.4` — no such image on Docker Hub
- **Fix:** Use `getwilds/annotsv:latest` instead (Fred Hutch maintained)

### SnpEff/SnpSift: Combined package
- **Failed:** `quay.io/biocontainers/snpsift:5.2--hdfd78af_1` — no such manifest
- **Fix:** Use `quay.io/biocontainers/snpeff:5.2--hdfd78af_1` — SnpSift is bundled inside the SnpEff package

### SnpEff Database: Azure blob storage outage
- **Failed:** `snpEff download GRCh38.105` and all database names — Azure blob storage returned 0-byte files for ALL URLs including `GRCh38.105`, `GRCh38.mane.1.2.ensembl`
- **Root cause:** SnpEff databases are hosted on Azure blob storage which was down
- **Fix:** Pivot to Ensembl VEP as alternative functional annotation tool

### ExpansionHunter: Different image formats
- **Failed:** `quay.io/biocontainers/expansionhunter:5.0.0--hd03093a_1` — manifest not found. Also `mgibio/expansionhunter:latest` — not found
- **Fix:** Use `weisburd/expansionhunter:latest` — binary at `/ExpansionHunter/bin/ExpansionHunter`, variant catalogs at `/pathogenic_repeats/GRCh38/`
- **Note:** This is ExpansionHunter v2.5.5 which uses `--repeat-specs` (directory), not `--variant-catalog` (single JSON)

### ExpansionHunter: Missing required --log parameter
- **Failed:** Container exits immediately with "the option '--log' is required but missing"
- **Fix:** Always include `--log /output/sample_eh.log` in the command

### StellarPGx: Empty Docker image
- **Failed:** `twesigomwedavid/stellarpgx:latest` — image exists but contains no StellarPGx binaries
- **Status:** UNSOLVED. No working Docker image found for StellarPGx as of Mar 2026.

## Tool-Specific Issues

### TelomereHunter: Permission denied
- **Failed:** `OSError: [Errno 13] Permission denied: '/output/sergio'` when writing output (sample name in error will vary)
- **Fix:** Add `--user root` flag to `docker run`

### TelomereHunter: pip install on host fails
- **Failed:** `pip install telomerehunter` on the host gives `UnicodeDecodeError` — Python environment issues
- **Fix:** Use Docker image `lgalarno/telomerehunter:latest` instead of native install

### HLA-LA: Graph not serialized (3 failures)
- **Failed attempt 1:** `zlskidmore/hla-la:latest` — graph at `src/additionalReferences/PRG_MHC_GRCh38_withIMGT/` exists but is NOT serialized. HLA-LA exits silently.
- **Failed attempt 2:** Copied graph files to `graphs/PRG_MHC_GRCh38_withIMGT/` — "graph not complete"
- **Failed attempt 3:** Ran `--action prepareGraph` in detached container — container exited during prep, graph still not serialized
- **Root cause:** HLA-LA requires a pre-serialized graph (~40GB RAM to prepare, takes hours). No standard Docker image includes it.
- **Fix:** Use `jiachenzdocker/hla-la:latest` (27.5GB image with pre-built graph) OR switch to T1K

### HLA-LA: Binary crash with pre-built graph image
- **Failed:** `jiachenzdocker/hla-la:latest` — read extraction succeeds but `HLA-LA` C++ binary crashes during graph alignment even with 32GB RAM and 8 threads. Error: "HLA-LA execution not successful."
- **Root cause:** Likely an incompatibility between the pre-built binary and the BAM data format, or an unmet memory requirement (the graph deserialization may need >32GB)
- **Status:** UNSOLVED. HLA-LA from WGS BAMs is unreliable in Docker. Alternative: use Sanitas clinical HLA typing results, or use arcas-hla or T1K with partial coordinates

### T1K: Coordinate file with wrong values
- **Failed:** `t1k-build.pl -d hla.dat -g reference.fasta.fai` produced coordinate file with `chr19 -1 -1 +` for all HLA genes
- **Root cause:** The `-g` parameter expects the actual reference FASTA (3.1GB), not the FAI index (158KB). The `AddGeneCoord.pl` script needs to align HLA sequences against the genome to find coordinates.
- **Fix:** Use `-g Homo_sapiens_assembly38.fasta` (the full FASTA, not the .fai)

### T1K: BAM extraction produces 0-byte FASTQ
- **Failed:** `bam-extractor` runs but produces empty `_candidate_1.fq` and `_candidate_2.fq`
- **Root cause:** Coordinate file had `-1 -1` coordinates (see above), so no genomic region was extracted
- **Fix:** Regenerate coordinate file with proper reference FASTA

### Cyrius CYP2D6: Inconclusive on short-read WGS
- **Result:** Cyrius returned `None` for CYP2D6 star alleles
- **Root cause:** CYP2D6 has extensive pseudogene homology (CYP2D7, CYP2D8) making short-read WGS unreliable
- **Mitigation:** Rely on Sanitas lab calls for CYP2D6; consider long-read sequencing in future

## bcftools/htslib Issues

### bgzip/tabix not in bcftools image PATH
- **Failed:** `staphb/bcftools:1.21` does not include `bgzip` or `tabix` in `$PATH`
- **Fix:** Use `bcftools view -Oz -o output.vcf.gz` (native bgzip output) and `bcftools index -t` (native tabix index) instead of piping to `bgzip`/`tabix`
- **Alternative:** Use `quay.io/biocontainers/samtools:1.21` which includes all htslib tools

### MIS VCF conversion: 0-byte output files
- **Failed:** All 22 chr files in `mis_ready/` were 0 bytes
- **Root cause:** Script used `bgzip -c > output.vcf.gz` which failed silently because bgzip wasn't available
- **Fix:** See above — use `bcftools view -Oz -o`

## VEP Cache Issues

### VEP INSTALL.pl: Permission denied on temp directory
- **Failed:** `Cannot open Local file /opt/vep/.vep/tmp/homo_sapiens_vep_112_GRCh38.tar.gz`
- **Fix:** Run with `--user root` and pre-create the temp dir: `mkdir -p /opt/vep/.vep/tmp && chmod 777 /opt/vep/.vep/tmp`

### VEP INSTALL.pl: Silent download failure
- **Failed:** Container exits after "downloading..." with no cache files extracted
- **Root cause:** The 17GB download timed out or was interrupted; INSTALL.pl doesn't retry
- **Fix:** Download manually with `wget -c` (supports resume), then extract with `tar xzf`, then run VEP with `--cache --dir_cache /path/to/cache`

## General Docker Tips

### Always use resource limits
```bash
docker run --cpus 4 --memory 8g ...
```
Without limits, tools like DeepVariant or minimap2 will consume ALL available RAM and crash the host.

### Build for amd64 from macOS
```bash
docker build --platform linux/amd64 ...
```
macOS is arm64; most Linux servers are amd64. Images built on Mac without `--platform` won't run on amd64 servers.

### Use --rm for one-shot containers
Always use `--rm` for analysis containers to avoid accumulating stopped containers. Use `-d` (detached) for long-running jobs.

### Always use --user root for write access
Most bioinformatics containers run as non-root users. If writing to bind-mounted volumes, add `--user root` to avoid permission issues.

## CI / Workflow Issues

### ShellCheck warnings still fail the GitHub Action
- **Observed:** `ludeeus/action-shellcheck@master` exits non-zero even when configured with `severity: warning`
- **Impact:** Once the check is required on `main`, "warning-only" findings still block merges
- **Fix:** Clear ShellCheck warnings before enabling required status checks, or explicitly relax the workflow instead of assuming warnings are advisory only

### Protect `main` only after CI is green
- **Observed:** Required status checks become a trap if you enable branch protection while the default branch or active PR branch is still red
- **Fix:** Get `Lint` and `Smoke Tests` green first, then enable required checks, block force-pushes, and block deletion

## MToolBox Issues

### MToolBox: No working Docker image exists
- **Failed:** `robertopreste/mtoolbox:latest` — "repository does not exist or may require docker login"
- **Also checked:** No image on quay.io/biocontainers, ghcr.io, or Docker Hub
- **Root cause:** MToolBox was never officially containerized. GitHub issue #107 (Mar 2022) confirms: "Not at the moment."
- **Fix:** Use GATK Mutect2 in mitochondrial mode instead (`broadinstitute/gatk:latest`). Mutect2 handles mitochondrial heteroplasmy detection natively and is well-maintained.

## CNVnator Issues

### CNVnator: Biocontainer tag with wrong build hash
- **Failed:** `quay.io/biocontainers/cnvnator:0.4.1--py312hc02a2a2_7` — manifest not found
- **Fix:** Use `quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11`. Biocontainer hashes encode the conda build hash and change between builds. Always verify at quay.io/repository/biocontainers/cnvnator.

## Delly Issues

### Delly: Biocontainer tag doesn't exist
- **Failed:** `quay.io/biocontainers/delly:1.2.9--ha41ced6_0` — manifest not found
- **Fix:** Use `quay.io/biocontainers/delly:1.7.3--hd6466ae_0` (latest as of Mar 2026). Biocontainer tags are version-specific and change frequently — always verify at quay.io/repository/biocontainers/delly.

### Delly: SV annotation phase takes 2-3 hours
- **Observed:** Delly's "SV annotation" step runs for 2-3 hours at 100% CPU on a 30X WGS genome. No new log output appears during this time, which can look like the process is stuck.
- **This is normal.** Delly genotypes every candidate SV site against the reference, which is CPU-intensive. Total runtime for 30X WGS: ~3-4 hours.
- **Tip:** Use `docker stats` to confirm the container is still using CPU. If CPU is at 0%, the process may actually be stuck.

### Delly: Output is BCF format, not VCF
- **Gotcha:** Delly writes BCF (binary VCF), not VCF. The output file has a `.bcf` extension.
- **Fix:** Convert with `bcftools view input.bcf -Oz -o output.vcf.gz` and index with `bcftools index -t output.vcf.gz`. The pipeline script handles this automatically.

### CNVnator: ROOT file appears empty (266 bytes) during tree extraction
- **Observed:** During the `-tree` step, the `.root` file stays at 266 bytes (just the ROOT header) until the entire BAM is parsed.
- **This is normal.** For a 30X WGS (~80GB BAM), the tree step takes ~5-10 minutes. The ROOT file grows to ~900MB-1.2GB only at the very end when the tree is flushed to disk.
- **If the container exits and the file is still 266 bytes:** Check if a corrupt `.root` file from a previous failed run is blocking it. Delete and retry.

## CPSR/PCGR Issues

### CPSR: --pcgr_dir path confusion
- **Failed:** `cpsr --pcgr_dir /genome/pcgr_data/data` → "Data directory (/genome/pcgr_data/data/data) does not exist"
- **Root cause:** CPSR internally appends `/data` to whatever `--pcgr_dir` you pass. If you point to the `data/` directory inside the extracted bundle, it looks for `data/data/`.
- **Fix:** Point `--pcgr_dir` to the **parent** of the `data/` directory: `--pcgr_dir /genome/pcgr_data` (not `/genome/pcgr_data/data`)

### CPSR: Docker image is inside PCGR
- **Failed:** `sigven/cpsr:2.0.0` does not exist on Docker Hub
- **Fix:** Use `sigven/pcgr:1.4.1` which bundles both `pcgr` and `cpsr` binaries at `/usr/local/bin/`

## Michigan Imputation Server Notes

### Minimum 20 samples per job
- MIS is designed for genotyping array data from cohort studies
- Single-sample WGS submissions may be rejected (20-sample minimum)
- For individual WGS: imputation adds minimal value (you already have 90%+ variant coverage)
- Main benefit for WGS would be **phasing**, not imputation

### Registration required
- Must create account at imputationserver.sph.umich.edu
- API tokens expire after 30 days
- Results auto-deleted after 7 days

### TOPMed panel is best for Europeans
- TOPMed Freeze 8 (r2): 132K samples, 705M variants, GRCh38 native
- Uses `chr` prefix (which GRCh38 BAMs already have)
- HRC r1.1 (32K samples) is European-centric but hg19 only

## Chip Data Conversion (Genotyping Arrays → VCF)

### plink silently corrupts single-sample homozygous ALT genotypes
- **What failed:** `plink --23file` (1.9) to import + `plink2 --ref-from-fa force` to fix REF/ALT
- **Why:** For single-sample data, ALL homozygous positions are monomorphic. plink's `.bim` stores only one allele for these. `--ref-from-fa` cannot create a proper ALT because there's no second allele slot. Homozygous ALT genotypes silently become homozygous REF.
- **Verified:** rs9939609 (FTO), genotype=AA, REF=T. plink: `REF=A, ALT=., GT=0/0` (WRONG). bcftools: `REF=T, ALT=A, GT=1/1` (CORRECT). ~66K positions (11%) corrupted.
- **Fix:** Use `bcftools convert --tsv2vcf -f <reference.fa>`. Single command, no intermediate binary format.

### plink 1.9 --23file quirks
- `--allow-extra-chr` cannot be used with `--23file`
- Female samples with Y calls (MyHeritage GSA PAR region) error with sex=F
- Sex inference defaults to male unless explicitly set

### MyHeritage CSV must be converted to TSV
- Quoted CSV with `"RSID","CHROMOSOME","POSITION","RESULT"` columns
- Strip `##` comments, header, quotes; convert commas to tabs

### bcftools hg19 VCF needs chr prefix before liftover
- hg19 reference uses numeric chromosomes; chain file expects chr prefix
- `bcftools annotate --rename-chrs` between conversion and liftover

### PharmCAT chip vs WGS results (MyHeritage GSA, verified 2026-03-31)
- **Correct:** CYP2B6 (*1/*6), CYP4F2 (*1/*6), DPYD (*5/*5), NUDT15 (*1/*2)
- **Missed:** CYP2C19 (25 missing), VKORC1 (1 missing)
- **Wrong:** CYP3A5 *1/*1 (should be *3/*3, 4 missing positions)
- Total: 888 missing PGx positions from the GSA chip

### ROH and PRS need special flags for chip data
- ROH: `-G30` required (no FORMAT/PL in chip VCF)
- PRS: `no-mean-imputation` required (single sample lacks allele frequencies)
- PRS matching: chip ~12% of large scoring files vs WGS ~28%
