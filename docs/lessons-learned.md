# Lessons Learned

Every failure encountered during pipeline development (Mar 2026), documented so they don't happen again.

## Docker Image Issues

### AnnotSV: Official image doesn't exist
- **Failed:** `bioinfochrustrasbourg/annotsv:3.4.4` — no such image on Docker Hub
- **Fix:** Use `getwilds/annotsv:3.4.4` instead (Fred Hutch maintained)

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
- **Failed:** `OSError: [Errno 13] Permission denied: '/output/<sample>'` when writing output (the sample name in the error path will vary)
- **Fix:** Add `--user root` flag to `docker run`

### TelomereHunter: pip install on host fails
- **Failed:** `pip install telomerehunter` on the host gives `UnicodeDecodeError` — Python environment issues
- **Fix:** Use the Docker image (`lgalarno/telomerehunter`, digest-pinned in `versions.env`) instead of native install

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

### CPSR: --pcgr_dir path confusion (PCGR 1.x, historical)
- **Failed:** `cpsr --pcgr_dir /genome/pcgr_data/data` → "Data directory (/genome/pcgr_data/data/data) does not exist"
- **Root cause:** CPSR 1.x internally appends `/data` to whatever `--pcgr_dir` you pass. If you point to the `data/` directory inside the extracted bundle, it looks for `data/data/`.
- **Fix (1.x):** Point `--pcgr_dir` to the **parent** of the `data/` directory: `--pcgr_dir /genome/pcgr_data` (not `/genome/pcgr_data/data`)
- **Superseded by PCGR 2.x** — the `--pcgr_dir` flag no longer exists. See migration notes below.

### CPSR: Docker image is inside PCGR
- **Failed:** `sigven/cpsr:2.0.0` does not exist on Docker Hub
- **Fix:** Use `sigven/pcgr:2.2.5` which bundles both `pcgr` and `cpsr` binaries at `/usr/local/bin/`

### PCGR 2.x Migration (1.4.1 to 2.2.5)
- **CLI completely changed:** The `--pcgr_dir` flag is gone. Replaced by `--refdata_dir` (for the ref data bundle) and `--vep_dir` (for the VEP cache). These are separate mount points inside the container.
- **Docker volumes changed:** PCGR 1.x used a single `-v ${GENOME_DIR}:/genome` mount. PCGR 2.x requires four separate mounts: VEP cache (`:/mnt/.vep`), ref data bundle (`:/mnt/bundle`), input VCFs (`:/mnt/inputs`), and outputs (`:/mnt/outputs`).
- **Data bundle is smaller and different:** The old monolithic ~21 GB bundle (`pcgr.databundle.grch38.20220203.tgz`) that included VEP cache is replaced by a smaller ~5 GB ref data bundle (`pcgr_ref_data.20250314.grch38.tgz`). VEP cache is now provided separately (reuse the same cache from step 13).
- **Bundle extraction requires extra step:** After `tar xzf`, the extracted `data/` directory must be moved into a version-stamped directory: `mkdir -p 20250314 && mv data/ 20250314/`. The `--refdata_dir` mount points to this version directory.
- **Docker tag:** `sigven/pcgr:1.4.1` → `sigven/pcgr:2.2.5`. The image still bundles both `pcgr` and `cpsr` binaries.
- **Old data bundle is incompatible:** If you have the 1.x bundle, you must download the 2.x bundle fresh. The directory structure and expected paths are completely different.

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

## Alternative Variant Caller Issues

### FreeBayes 1.3.7: SIGILL crash (exit code 132)
- **Failed:** `quay.io/biocontainers/freebayes:1.3.7--h1870644_0` — `freebayes --version` works, but actual variant calling triggers `SIGILL` (illegal instruction, exit code 132)
- **Tested on:** Intel i5-14500
- **Root cause:** Likely a build-time CPU optimization mismatch in the 1.3.7 biocontainer binary
- **Fix:** Use `quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2` which works correctly

### FreeBayes: Memory grows to ~13 GB on full genome
- **Observed:** FreeBayes memory usage grows unpredictably during full-genome runs: 463MB at 30 min, 6.4GB at 60 min, 12.7GB at 90 min, then stabilizes ~12GB
- **Original limit:** `--memory 16g` was too tight — would have OOM-killed at 80% usage
- **Fix:** Use `--memory 32g` for full-genome runs. Peak observed was 12.8GB but growth is non-linear and region-dependent.

### FreeBayes: Single-threaded, no parallelism
- **Observed:** FreeBayes has no `-t` or `--threads` flag. Full 30X WGS takes 8-12 hours.
- **Workaround:** Use `--region chr22` (or `INTERVALS=chr22`) for quick testing (~20-40 min)
- **For production:** Consider GNU parallel with per-chromosome regions, then merge VCFs

### GATK HaplotypeCaller: bcftools index fails on existing .tbi
- **Failed:** `bcftools index -t` fails with "index file exists" after GATK already creates its own `.tbi`
- **Fix:** Use `bcftools index -ft` (with `-f` force flag) to overwrite the GATK-generated index

### GATK HaplotypeCaller: Requires .dict file
- **Failed:** GATK HaplotypeCaller fails if `Homo_sapiens_assembly38.dict` is missing
- **Fix:** Generate once with `gatk CreateSequenceDictionary -R /genome/reference/Homo_sapiens_assembly38.fasta`

### bcftools isec: -R vs -r for region strings
- **Failed:** `bcftools isec -R chr22` treats `-R` (uppercase) as a BED file path, fails with "file not found"
- **Fix:** Use `-r chr22` (lowercase) for region strings. `-R` expects a file.

### TIDDIT >=3.9: Requires BWA index for local assembly
- **Failed:** `tiddit --sv` exits with "The reference must be indexed using bwa index; run bwa index, or skip local assembly (--skip_assembly)"
- **Root cause:** TIDDIT 3.9+ uses local assembly for breakpoint refinement, which requires BWA index files alongside the reference
- **Fix:** Use `--skip_assembly` when using minimap2 alignments (no BWA index available). If using BWA-MEM2 alignment, the index files are compatible.

### TIDDIT: Image tag 3.7.0 doesn't exist on quay.io
- **Failed:** `quay.io/biocontainers/tiddit:3.7.0--py312h24f4cff_1` — manifest unknown
- **Fix:** Use `quay.io/biocontainers/tiddit:3.9.5--py312h6e8b409_0`. Always verify tags at quay.io/repository/biocontainers/tiddit.

### Strelka2: --callRegions needs bgzipped + tabixed BED
- **Failed:** `--callRegions reference.fasta.fai` → "Can't find expected call-regions bed index file"
- **Fix:** Create a proper bgzipped BED file with tabix index. Use GATK container for bgzip/tabix (not in bcftools or samtools staphb images).

### bgzip/tabix not in staphb/samtools or staphb/bcftools images
- **Observed:** Neither `staphb/samtools:1.20` nor `staphb/bcftools:1.21` include `bgzip` or `tabix` in PATH
- **Fix:** Use `broadinstitute/gatk:4.6.1.0` which has both at `/usr/bin/bgzip` and `/usr/bin/tabix`. Or use `bcftools view -Oz` as a bgzip alternative.

### FreeBayes chr22 variant count (3x more than DeepVariant)
- **Observed:** FreeBayes calls ~247K variants on chr22 vs DeepVariant ~91K and GATK ~69K
- **Interpretation:** The ~200K FreeBayes-unique variants are mostly false positives. FreeBayes maximizes sensitivity at the cost of precision.
- **Recommendation:** Always quality-filter FreeBayes output with `bcftools filter` or `vcffilter` before use.

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

## v0.3.0 Tool Additions (Apr 2026)

### ExpansionHunter v5.0.0: Completely different CLI from v2.5.5
- **Old (v2.5.5):** `--bam`, `--ref-fasta`, `--repeat-specs` (directory), `--vcf`, `--json`, `--log` (all required)
- **New (v5.0.0):** `--reads`, `--reference`, `--variant-catalog` (single JSON), `--output-prefix`, `--threads`
- **The `--log` flag is gone** in v5.0.0. Do NOT pass it or the command will fail.
- **Biocontainer image:** `quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5`. Binary is `ExpansionHunter` (on PATH, not at `/ExpansionHunter/bin/`).
- **Bundled catalog:** `/usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json` (31 loci). No need to download separately.

### GRIDSS: Requires BWA index (not minimap2)
- **Failed:** GRIDSS exits with "BWA index not found" when using default minimap2 alignment
- **Fix:** Either (a) align with BWA-MEM2 first (`02a-alignment-bwamem2.sh`), or (b) generate BWA index files separately. The `04b-gridss.sh` script validates this and prints instructions.
- **Note:** GRIDSS outputs ALL SVs as BND (breakend) notation. Standard DEL/DUP/INV types require post-processing conversion for SURVIVOR merge compatibility.

### GRIDSS: 32 GB memory requirement
- **Observed:** GRIDSS assembly-based SV calling needs ~28 GB JVM heap for 30X WGS
- **Fix:** Container runs with `--memory 32g` and `-Xmx28g` JVM argument. Will fail silently on machines with < 32 GB RAM.

### GRIDSS: ENCODE blacklist download
- **Observed:** GRIDSS benefits from an ENCODE blacklist to suppress known artifact regions
- **Fix:** Script auto-downloads `ENCFF356LFX.bed.gz` (hg38 blacklist) from ENCODE on first run. If download fails (offline), GRIDSS runs without it (lower precision but still functional).

### fastp: Maximum 16 threads despite --workers flag
- **Observed:** fastp accepts `-w` (workers) up to 16. Values above 16 are clamped to 16. The `-w` flag controls I/O worker threads; actual adapter detection is single-threaded.
- **Recommendation:** Use `-w` matching `THREADS` up to 16. For most WGS runs, `-w 4` is sufficient.

### fastp: BGI/MGI adapter auto-detection
- **Observed:** fastp's `--detect_adapter_for_pe` works for Illumina, BGI, and MGI adapters without specifying adapter sequences. BGI/MGI adapters are compiled into fastp's `knownadapters.h`.
- **No action needed:** The `--detect_adapter_for_pe` flag handles all common sequencing platforms.

### mosdepth: --fast-mode skips per-base output
- **Observed:** `--fast-mode` uses a simpler, faster counting method and does NOT write the per-base `.per-base.bed.gz` file. This saves ~2 GB of output and cuts runtime by ~40%.
- **Fix:** Always use `--fast-mode` unless per-base resolution is specifically needed.

### MultiQC: Auto-discovers fastp JSON by content, not filename
- **Observed:** MultiQC identifies fastp output by looking for `"before_filtering": {` in JSON files, not by filename pattern. Files must end in `.json`.
- **Tip:** Ensure fastp's `-j` output uses `.json` extension.

### Octopus: No issues observed
- Docker image `dancooke/octopus:0.7.4` works out of the box for germline calling
- Supports `--threads` for parallelism (unlike FreeBayes)
- Typical memory: 8-12 GB for 30X WGS (much less than FreeBayes peak of 13 GB)

## PharmCAT 3.x Migration (2.15.5 to 3.2.0)

### Preprocessor script renamed (no .py extension)
- **Old (2.15.5):** `python3 /pharmcat/pharmcat_vcf_preprocessor.py`
- **New (3.2.0):** `python3 /pharmcat/pharmcat_vcf_preprocessor`
- **Impact:** Step 7 preprocessor command must drop the `.py` suffix or the container exits with "No such file"

### Reporter flags: must be explicit for both formats
- **Old (2.15.5):** `-reporterJson` produced JSON; HTML was always generated by default
- **New (3.2.0):** If ANY format flag is specified, ONLY those formats are produced. To get both HTML and JSON, you must pass `-reporterJson -reporterHtml`
- **Impact:** Step 7 now passes both flags explicitly. Without `-reporterHtml`, the HTML report (used for manual review) would silently stop being generated.

### JSON property rename: wildtypeAllele to referenceAllele
- **Old (2.15.5):** `wildtypeAllele` property in gene result JSON objects
- **New (3.2.0):** Renamed to `referenceAllele`
- **Impact:** Step 27's JSON parser does not use this property directly, so no code change was needed. Any downstream scripts or notebooks that parse `wildtypeAllele` must be updated.

### New features in 3.2.0
- **NAT2 calling:** PharmCAT 3.x includes improved NAT2 acetylator status calling
- **BCF support:** Preprocessor now accepts BCF input files directly (no conversion needed)
- **Single-gene calling:** New `-g` flag allows running PharmCAT on a single gene (useful for targeted re-analysis)

## Nextflow version compatibility (2026-06)

### The pipeline runs cleanly on Nextflow 25.10.4; 24.x and 26.x currently fail at parse time
- **Observed:** A full run requires **Nextflow 25.10.4** (the validated version). Other versions fail before any process executes:
  - **26.04.4** — the strict config parser rejects top-level `def`/variable declarations in `nextflow.config` ("Variable declarations cannot be mixed with config statements"), and then the `def check_max(...)` function in `conf/base.config` ("Unexpected input: '('").
  - **24.04.4** — the DSL2 module parser flags the optional annotation inputs in `modules/local/vcfanno/main.nf` as "Variable already defined in the process scope" (`cadd_snv`/`cadd_indel`/`spliceai_*`/`revel`/`alphamissense`, referenced inside the `def has_nochr`/`def has_chr` expressions). 25.10.4 tolerates this; 24.04.4 does not.
- **Fix status:** `nextflow.config` is strict-parser-clean — the execution-report timestamp is inlined into each report path (no top-level `def`; see #30/#31), which also preserves per-run report history. Full NF-26 support is still pending: migrating `conf/base.config`'s `check_max()` → `process.resourceLimits` and refactoring the vcfanno input scope (tracked in `docs/sota-update-2026-06.md`). Pin `NXF_VER=25.10.4` to run.
- **Tip:** `NXF_VER=25.10.4 nextflow run main.nf ...`. The `manifest.nextflowVersion` floor is raised to `25.10.0` so the known-broken 24.x is rejected up front; 26.x is gated by comment until the migration lands.

### CYP2D6 structural alleles: pypgx resolves *5 deletions where Cyrius and PharmCAT return "No Result"
- **Observed:** For a homozygous CYP2D6 whole-gene deletion (*5/*5), **Cyrius can return `None/None`** (Total_CN null — its copy-number consensus cannot resolve the locus) and **PharmCAT reports `Unknown/Unknown — No Result`** (it does not call the structural *5 from a plain VCF), while **pypgx (BAM-based, SV-aware) resolves `*5/*5 — Poor Metabolizer` with `SV_detected: Yes`.**
- **Impact:** Updates the older "rely on lab calls" note above — for CYP2D6 deletion/duplication alleles, pypgx on the BAM is the authoritative caller. Do **not** read a Cyrius `None/None` as "no deletion." Keep all three callers (PharmCAT star alleles, Cyrius, pypgx) and reconcile; pypgx wins for CNV/SV-driven star alleles (a Poor Metabolizer has no functional CYP2D6 → major impact on CYP2D6-cleared drugs such as codeine/tramadol/tamoxifen).

## CNVpytor migration (2026-07)

### CNVpytor 1.3.2 biocontainer ships without GC/mask data and its downloader is broken
- **Failed:** `cnvpytor -his` aborts with `Some reference genome resource files are missing. Run 'cnvpytor -download'` — the `quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0` image's `cnvpytor/data/` dir contains only an empty `readme.txt`.
- **Failed:** `cnvpytor -download` itself crashes in 1.3.2 (`AttributeError: 'PosixPath' object has no attribute 'split'` in `genome.py`), so it cannot self-heal.
- **Fix:** Pre-download the pinned **v1.3.2** GC/mask files and bind-mount them onto the container's package data dir (`/usr/local/lib/python3.12/site-packages/cnvpytor/data`). Runs then work fully offline (verified with `--network none`). See `docs/00-reference-setup.md`.

### CNVpytor's resource check requires every genome's files to exist, not just hg38
- **Observed:** `genome.py:check_resources()` iterates every bundled reference genome (hg19, hg38, chm13v2.0, chm13v1.1, kn99) and `os.path.exists()`-checks each `gc_file`/`mask_file`. Mounting only `gc_hg38.pytor`+`mask_hg38.pytor` still fails the check.
- **Fix:** Provide all seven files. Only `gc_hg38.pytor`/`mask_hg38.pytor` are actually read for an hg38 BAM; the others just need to exist (the check is existence-only). All seven total ~90 MB.

### CNVpytor container has no bcftools/bgzip
- **Observed:** Unlike the old CNVnator biocontainer, `cnvpytor:1.3.2` bundles no bcftools/bgzip/tabix/samtools.
- **Fix:** The bash step runs a separate `staphb/bcftools` container for VCF normalization; the Nextflow module splits into two processes (`CNVPYTOR` calls → `CNVPYTOR_VCF` reheader/sort/index). `cnvpytor -view` emits a proper VCFv4.2 (SVTYPE/END/SVLEN, ALT DEL/DUP/LOH, GT/CN) but only carries `##contig` lines for processed chromosomes — reheader from the reference `.fai` before merging with Manta/Delly in step 22.

## Real-data lessons from full-WGS runs (2026-07)

### CNVpytor chokes on GRCh38 ALT/HLA/decoy contigs — restrict `-rd` to canonical chromosomes
- **Failed:** On a real full-reference GRCh38 BAM (hundreds of ALT/HLA/decoy contigs), `cnvpytor -rd` (no `-chrom`) crashes/stalls during read-depth import and produces **no calls** (only a tiny stub `.pytor`). The bundled GC-correction data covers only the main chromosomes. A chr20-only validation did **not** surface this.
- **Fix:** Pass `-chrom chr1 … chr22 chrX chrY` to the `-rd` step (`scripts/18-cnvpytor.sh` + `modules/local/cnvpytor`). Validated on two real 30× genomes (≈2200 and ≈2400 canonical CNV calls, ~99% concordant with the prior CNVnator counts). Lesson: validate depth-based callers on a **full** BAM, not a single chromosome.

### pypgx 0.27.0 biocontainer ships pandas 3.0.3 → every gene fails
- **Failed:** `pypgx:0.27.0--pyh106432d_0` fails all genes at runtime with `pandas.errors.LossySetitemError` / `TypeError: Invalid value` — its code assigns floats into int columns, which pandas ≥2.1 (the image bundles **3.0.3**) rejects. The container smoke test only runs `pypgx --version`, so CI cannot see it.
- **Status:** Needs a pypgx build with a compatible pandas (<2.1) or an upstream fix; tracked follow-up. Until then step 32 is non-functional and CYP2D6/others must come from PharmCAT/Cyrius.

### Stranger over-flags RFC1 (CANVAS) from short reads — do not read it as a diagnosis
- **Observed:** Stranger can report RFC1 `STR_STATUS=full_mutation` for a modest expansion (e.g. 51/73 of the degenerate `AARRG` motif). CANVAS requires the **AAGGG** motif specifically, **biallelic**, at **~400–2000+** repeats — short-read ExpansionHunter cannot resolve AAGGG vs the benign AAAAG, and the catalog's `STR_PATHOLOGIC_MIN=12` is not the clinical threshold.
- **Interpretation:** Treat an RFC1 flag as **uninterpretable from short-read WGS** — confirm with motif-aware/flanking-PCR testing only if clinically indicated. (Documented in `docs/09b-stranger.md`.)
