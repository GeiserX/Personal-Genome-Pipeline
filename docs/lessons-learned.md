# Lessons Learned

Every failure encountered during pipeline development (Mar 2026), documented so they don't happen again.

## Docker Image Issues

### AnnotSV: Official image doesn't exist
- **Failed:** `bioinfochrustrasbourg/annotsv:3.4.4` — no such image on Docker Hub
- **Fix (then):** Used `getwilds/annotsv:3.4.4` (Fred Hutch maintained). That image is retired here; the pipeline now pins the Bioconda build `quay.io/biocontainers/annotsv` (see `ANNOTSV_IMAGE` in `versions.env`).

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
- **Status:** UNSOLVED. HLA-LA from WGS BAMs is unreliable in Docker. Alternative: a clinical lab result, if you have one, or arcas-hla or T1K with partial coordinates

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
- **Mitigation:** Use a clinical lab result, if you have one, for CYP2D6; consider long-read sequencing in future

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
- **Fix:** Use GATK Mutect2 in mitochondrial mode instead (`broadinstitute/gatk:4.6.2.0`, the version the scripts pin). Mutect2 handles mitochondrial heteroplasmy detection natively and is well-maintained.

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
- **Fix:** Use `--skip_assembly` when no classic BWA index is available. TIDDIT's assembly step calls classic `bwa`, which needs the classic index files (`.amb .ann .bwt .pac .sa`). BWA-MEM2's index (`.bwt.2bit.64`, `.0123`) is a different format and does not work, even after aligning with BWA-MEM2.

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
- **Verified:** on one homozygous-alternate chip position (chip genotype AA where the reference base is T): plink wrote `REF=A, ALT=., GT=0/0` (WRONG); bcftools wrote `REF=T, ALT=A, GT=1/1` (CORRECT). ~66K positions (11%) corrupted.
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
- **Correct:** CYP2B6, CYP4F2, DPYD and NUDT15 matched the WGS diplotypes
- **Missed:** CYP2C19 (25 missing), VKORC1 (1 missing)
- **Wrong:** CYP3A5 came out as the reference diplotype because 4 of its defining positions are not on the chip; WGS called a different diplotype
- **Fix:** treat a chip-only call as unresolved when any of its defining positions are missing from the chip, and confirm it from WGS before relying on it
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
- **Fix:** Build the classic BWA index (`bwa index` on the reference, about an hour). GRIDSS runs classic `bwa` internally, so BWA-MEM2's index files (`.bwt.2bit.64`) do not help, and aligning with BWA-MEM2 first does not provide the index it needs. The `04b-gridss.sh` script checks for the classic files and prints the command.
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

### The pipeline runs on Nextflow 25.10.x and 26.04.x; 24.x fails at parse time
- **Observed:** A full run requires **Nextflow 25.10.8** (the validated version, `NEXTFLOW_VERSION` in `versions.env`; 25.10.4 when this was found). Other versions fail before any process executes:
  - **26.04.4** — the strict config parser rejects top-level `def`/variable declarations in `nextflow.config` ("Variable declarations cannot be mixed with config statements"), and then the `def check_max(...)` function in `conf/base.config` ("Unexpected input: '('").
  - **24.04.4** — the DSL2 module parser flags the optional annotation inputs in `modules/local/vcfanno/main.nf` as "Variable already defined in the process scope" (`cadd_snv`/`cadd_indel`/`spliceai_*`/`revel`/`alphamissense`, referenced inside the `def has_nochr`/`def has_chr` expressions). 25.10.4 tolerates this; 24.04.4 does not.
- **Fix status:** `nextflow.config` is strict-parser-clean — the execution-report timestamp is inlined into each report path (no top-level `def`; see #30/#31), which also preserves per-run report history. `conf/base.config` has since replaced `check_max()` with `process.resourceLimits`, which clears the second 26.x parse error above. The last 26.x error was `VCF_PRECHECK` setting `FILTER_STATUS` and `FILTER_COUNTS` without declaring them as `env()` outputs; it declares them now. CI (`nextflow.yml`) runs `nextflow lint` and the all-tools stub run on both the pinned 25.10.8 and the newest 26.04.x release; the real-data runs (`e2e.yml`) use 25.10.8. The vcfanno input scope is not refactored, so 24.x still fails.
- **Tip:** `NXF_VER=25.10.8 nextflow run main.nf ...` runs the release the real-data tests use; `run-all.sh` sets `NXF_VER` from `versions.env` unless you set it yourself. The `manifest.nextflowVersion` floor is `!>=25.10.4`, so 24.x is rejected up front; there is no cap, and 26.04.x passes the lint and stub checks.

### CYP2D6 copy number: no single caller settles it (withdrawn lesson)
- **Withdrawn:** an earlier version of this entry told readers to trust pypgx on the BAM over Cyrius and PharmCAT for CYP2D6 deletion and duplication alleles, because pypgx made a deletion call where Cyrius returned `None/None` and PharmCAT `No Result`. That conclusion does not hold.
- **Why:** pypgx calls copy number from read depth. On a reference with ALT contigs and an aligner that is not run ALT-aware, reads from the CYP2D locus can split between the primary copy and an ALT copy, depth on the primary drops, and a depth-based caller can report a whole-gene deletion that is not there.
- **Rule:** on a reference with ALT contigs, compare CYP2D6 depth with its flanks before trusting any copy-number call from pypgx or Cyrius, and report CYP2D6 only when two callers agree. A Cyrius `None/None` means "no call", not "no deletion" and not "deletion". If you have a clinical lab result for CYP2D6, it outranks all of these.

### A reference with ALT contigs and an aligner that is not run ALT-aware thin the depth at paralogous loci
- **Failed:** `setup.sh` installed the Broad `Homo_sapiens_assembly38.fasta` (3,366 sequences, with ALT, HLA and decoy contigs), and step 02 aligned with plain `minimap2 -x sr`. A read that matches a primary locus and its ALT copy equally well got MAPQ 0, and callers skip MAPQ 0 reads, so depth thinned at CYP2D6, the MHC and KIR. Measured on the HG002 fixture (`ALT depth A/B` workflow): depth at MAPQ >= 1 was 0.0x at CYP2D6 and HLA-A where the no-ALT analysis set gives 16.8x and 26.9x, and 85% of the MHC reads were aligned to ALT or HLA contigs. See [realignment](realignment.md#how-much-depth-alt-contigs-cost).
- **Root cause:** ALT contigs only help an aligner that is run ALT-aware (BWA-MEM with the `.alt` file and its post-processing, or DRAGEN's graph reference); the pipeline runs none. The extra contigs also forced per-contig workarounds in CNVpytor and Delly.
- **Fix:** the default reference is NCBI's GRCh38 no-ALT analysis set, `reference/GRCh38_no_alt_analysis_set.fasta` (195 sequences), under a new name so an old file is never read by accident. `validate-setup.sh` fails on a reference with `_alt` or `HLA-` contigs (unless `ALLOW_ALT_REFERENCE=true`) and on a BAM whose `@SQ` names and lengths differ from the `.fai`, naming the first difference.
- **Rule:** a change of reference means realigning every sample, from FASTQ or from the old BAM ([realignment](realignment.md)). Never run a step on a BAM whose header does not match the reference; `validate-setup.sh` checks it.

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

### pypgx image must match the bundle version — 0.27.0 mismatched 0.26.0 and shipped pandas 3.0.3
- **Failed:** `pypgx:0.27.0--pyh106432d_0` failed **all** genes at runtime with `pandas.errors.LossySetitemError` / `TypeError: Invalid value` — pypgx assigns floats into int columns, which pandas ≥2.1 (the 0.27.0 image bundles **3.0.3**) rejects. It was also mismatched against the pipeline's **0.26.0** `pypgx-bundle`. The smoke test only runs `pypgx --version`, so CI could not catch it.
- **Fix (v0.8.2):** pin the image to **`pypgx:0.26.0--pyh7e72e81_0`** to match the bundle (pandas 2.3.3, no error). Verified end-to-end on a real 30× genome: **23/23 genes succeed**. Lesson: the pypgx image version MUST track the `pypgx-bundle` branch, and bumping one without the other (and without a real gene-call test) silently breaks step 32.

### Stranger over-flags RFC1 (CANVAS) from short reads — do not read it as a diagnosis
- **Observed:** Stranger can report RFC1 `STR_STATUS=full_mutation` for a modest expansion (e.g. 51/73 of the degenerate `AARRG` motif). CANVAS requires the **AAGGG** motif specifically, **biallelic**, at **~400–2000+** repeats — short-read ExpansionHunter cannot resolve AAGGG vs the benign AAAAG, and the catalog's `STR_PATHOLOGIC_MIN=12` is not the clinical threshold.
- **Interpretation:** Treat an RFC1 flag as **uninterpretable from short-read WGS** — confirm with motif-aware/flanking-PCR testing only if clinically indicated. (Documented in `docs/09b-stranger.md`.)

## Vendor VCF intake (2026-10)

### PharmCAT's Java step stops on a backslash in a valid header line
- **Failed:** `pharmcat.jar` exits with `Error parsing metadata: character to be escaped is missing` when a `##FILTER` or `##INFO` line holds a backslash. bcftools writes such a line itself for a soft filter with a quoted string (`bcftools filter -s LowDP -e 'FORMAT/DP<10 && GT!="0/0"'`), and the VCF 4.3 spec asks for the escape. The preprocessor copies the line through; the bundled vcf-parser 0.3.1 refuses it. PharmCAT 3.4.0 bundles the same parser.
- **Fix:** step 07 and the PHARMCAT module rewrite the header of PharmCAT's own copy (`\"` to `'`, any other `\` to `/`, on `##` lines) before the Java step. The calls are the same as without the line. Remove the rewrite once a PharmCAT release bundles a newer vcf-parser; bumping PharmCAT is not the fix.

### Ensembl contig names fail silently
- **Failed:** with contigs named `1`, `MT`, `bcftools view -r chrM` returns no record with exit 0, so haplogrep3 writes a header-only file, and the ROH summary's `chrX|chrY` filter lets `X` through. Only the ClinVar screen failed loudly, and only when selected.
- **Fix:** `VCF_PRECHECK` stops the run before any analysis when no contig holding records is chr-named, and prints the `bcftools annotate --rename-chrs` command. Lesson: a step that selects a region by name must not read zero records as a clean result.

### `bcftools norm -f` stops at the first contig the reference lacks
- **Failed:** exit 255 (`The sequence "NT_113889.1" was not found`) on the full chr-renamed ClinVar file, and on a sample record on a scaffold the reference lacks (a half-done rename, or the no-alt reference).
- **Fix:** the ClinVar screen reads both files only on the contigs they share with each other and the reference, with a targets file whose end is a constant (`chrom 1 2147483647`): a header without contig lengths would otherwise give an empty region and zero records.

### `workflow.onComplete` saw `workflow` as null
- **Failed:** `Cannot get property 'success' on null object` on every run under Nextflow 25.10: the handler runs with the script binding's variable map as its delegate, and a map answers null for a name it lacks.
- **Fix:** the handler reads local variables (`run_info`, `run_log`, `outdir`) set in the workflow body; closures resolve those where they are written.
