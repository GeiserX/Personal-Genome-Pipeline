# Troubleshooting Guide

Definitive reference for diagnosing and fixing problems with the genomics pipeline. Organized by symptom so you can Ctrl+F your error message and find the fix.

**Before troubleshooting:** Run `docker info` to confirm Docker is running, and check `docker stats` for resource usage. Most failures are either Docker resource limits or wrong file paths inside containers.

---

## Table of Contents

- [Docker Issues](#docker-issues)
- [Input Data Issues](#input-data-issues)
- [Per-Step Troubleshooting](#per-step-troubleshooting)
- [Performance Issues](#performance-issues)
- [Output Issues](#output-issues)
- [Getting Help](#getting-help)

---

## Docker Issues

### Container exits immediately with no output

**Symptom:** `docker run` returns instantly. No error message. Exit code 137.

**Cause:** The container was OOM-killed (Out Of Memory). Docker silently kills containers that exceed their `--memory` limit. Exit code 137 = SIGKILL from the kernel OOM killer.

**How to confirm:**
```bash
# Check the last container's exit status
docker ps -a --latest --format "{{.Status}}"
# Output like "Exited (137)" confirms OOM

# Check system logs for OOM events
dmesg | grep -i "oom\|killed" | tail -10

# On Docker Desktop (Mac/Windows)
# Check Docker Desktop > Troubleshoot > Logs
```

**Fix:**
1. Increase the `--memory` flag in the script that failed
2. Reduce parallelism (run fewer steps simultaneously)
3. Increase Docker Desktop memory allocation (see [Docker Desktop not enough memory](#docker-desktop-not-enough-memory))
4. For DeepVariant, reduce `--num_shards` (each shard needs ~2-4 GB)

---

### Container exits with code 1 but no error message

**Symptom:** Script fails, exit code 1, but the container printed nothing useful.

**Cause:** Most bioinformatics tools print errors to stderr, which may not be captured depending on how you ran the script.

**Fix:**
```bash
# Re-run with full debug output
bash -x scripts/03-deepvariant.sh your_name 2>&1 | tee debug.log

# Or run the docker command manually with -it for interactive output
docker run --rm -it \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  google/deepvariant:1.6.0 \
  /bin/bash
# Then run the command inside the container to see the error
```

---

### "No such image" / "manifest for ... not found"

**Symptom:** `Error response from daemon: manifest for quay.io/biocontainers/toolname:tag not found`

**Cause:** Biocontainer image tags change frequently. The hash suffix (e.g., `--hd03093a_1`) encodes the conda build environment and breaks between releases.

**Fix:**
1. Check the current tag at the container registry:
   ```bash
   # For quay.io images
   # Visit: https://quay.io/repository/biocontainers/TOOLNAME?tab=tags

   # For Docker Hub images
   docker search toolname
   ```
2. Known problematic images and their correct tags (as of March 2026):

   | Tool | Wrong Tag | Correct Tag |
   |---|---|---|
   | CNVnator | `cnvnator:0.4.1--py312hc02a2a2_7` | `cnvnator:0.4.1--py312h99c8fb2_11` |
   | Delly | `delly:1.2.9--ha41ced6_0` | `delly:1.7.3--hd6466ae_0` |
   | ExpansionHunter | `expansionhunter:5.0.0--hd03093a_1` | `weisburd/expansionhunter:latest` |
   | AnnotSV | `bioinfochrustrasbourg/annotsv:3.4.4` | `getwilds/annotsv:latest` |
   | CPSR | `sigven/cpsr:2.0.0` | `sigven/pcgr:1.4.1` (bundles both) |
   | SnpSift | `quay.io/biocontainers/snpsift:5.2--hdfd78af_1` | `quay.io/biocontainers/snpeff:5.2--hdfd78af_1` (bundled) |
   | MToolBox | `robertopreste/mtoolbox:latest` | Does not exist. Use `broadinstitute/gatk:4.6.1.0` instead |

3. If an image disappears entirely, check if the tool has an official Docker image on GitHub Container Registry (`ghcr.io`), Docker Hub, or the tool's documentation.

---

### "Permission denied" writing to mounted volumes

**Symptom:** `OSError: [Errno 13] Permission denied: '/genome/sample/output'` or similar.

**Cause:** Most bioinformatics Docker images run as a non-root user (e.g., UID 1000). When writing to bind-mounted host directories, the container user may not have write access.

**Fix:** Add `--user root` to the `docker run` command. All pipeline scripts already include this flag, but if you are running commands manually:

```bash
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  your-image:tag \
  your-command
```

**Note:** Files created with `--user root` will be owned by root on the host. To fix ownership afterward:
```bash
sudo chown -R $(whoami) ${GENOME_DIR}/${SAMPLE}/
```

---

### Docker Desktop not enough memory (Mac/Windows)

**Symptom:** Containers fail with exit code 137 even though your machine has enough RAM.

**Cause:** Docker Desktop runs inside a virtual machine with a fixed memory allocation. Default is often 2-4 GB, which is far too little for genomics.

**Fix for macOS:**
1. Open Docker Desktop > Settings > Resources
2. Set Memory to at least **16 GB** (32 GB recommended)
3. Set Swap to at least **4 GB**
4. Click Apply & Restart

**Fix for Windows (WSL2):**
Create or edit `%UserProfile%\.wslconfig`:
```ini
[wsl2]
memory=24GB
swap=8GB
processors=8
```
Then restart WSL: `wsl --shutdown` and reopen your terminal.

**Verification:**
```bash
docker info | grep "Total Memory"
# Should show close to what you configured
```

---

### Containers stuck / running forever

**Symptom:** A container has been running for much longer than expected. No output. Not using CPU.

**How to diagnose:**
```bash
# Check what's running
docker ps

# Check resource usage — is the container actually doing work?
docker stats --no-stream

# Check container logs
docker logs <container_id>

# If CPU is at 0% for a long time, the container may be stuck
```

**Common causes:**
1. **Waiting for stdin:** Some tools expect interactive input. Fix: add `-i` flag or ensure all parameters are passed on the command line.
2. **I/O blocked on slow storage:** NFS/SMB/USB mounts can stall Docker I/O. Check `docker stats` — if CPU is 0% and memory is stable, the container is likely waiting for I/O.
3. **Java heap space:** Tools like PharmCAT or GATK may be stuck in garbage collection. Fix: increase `--memory`.

**Fix:**
```bash
# Kill the stuck container
docker kill <container_id>

# Or if you started with run-all.sh, kill all pipeline containers
docker ps --format "{{.ID}} {{.Image}}" | grep -E "deepvariant|samtools|bcftools|manta|vep|pcgr|pharmcat" | awk '{print $1}' | xargs docker kill
```

---

### "No space left on device" mid-analysis

**Symptom:** Container crashes with "No space left on device", "ENOSPC", or "write failed: No space left".

**Cause:** Docker needs free space in TWO places:
1. **Your data directory** (`GENOME_DIR`) for output files
2. **Docker's storage driver** (usually `/var/lib/docker`) for container layers and temp files

**Diagnosis:**
```bash
# Check host filesystem
df -h ${GENOME_DIR}
df -h /var/lib/docker   # Linux
# Docker Desktop stores data differently — check Docker Desktop settings

# Check Docker disk usage
docker system df
```

**Fix:**
```bash
# Clean up Docker (removes unused containers, images, build cache)
docker system prune -a
# WARNING: This removes ALL unused images. You'll need to re-pull them.

# Safer: only remove stopped containers and dangling images
docker container prune
docker image prune

# Delete intermediate files from completed steps
rm -f ${GENOME_DIR}/${SAMPLE}/cnvnator/*.root       # 5-15 GB each
rm -f ${GENOME_DIR}/${SAMPLE}/delly/*.bcf            # After VCF conversion
```

**Space requirements per sample:**

| Phase | Space Needed | Cumulative |
|---|---|---|
| FASTQ input | 60-90 GB | 90 GB |
| BAM (step 2) | 80-120 GB | 210 GB |
| VCF + analyses | 10-30 GB | 240 GB |
| VEP output | 2-5 GB | 245 GB |
| CNVnator ROOT file | 5-15 GB | 260 GB |
| **Recommended free** | | **500 GB** |

---

## Input Data Issues

### Wrong genome build (hg19 vs hg38)

**Symptom:** Variant calling produces far fewer variants than expected (e.g., 100K instead of 4-5 million). Or ClinVar screen finds 0 hits.

**How to detect:**
```bash
# Check BAM header for chromosome naming and lengths
docker run --rm -v "${GENOME_DIR}:/genome" staphb/samtools:1.20 \
  samtools view -H /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam | grep "^@SQ" | head -3

# GRCh38 (hg38): SN:chr1  LN:248956422
# GRCh37 (hg19): SN:1     LN:249250621   (or SN:chr1 LN:249250621)

# For VCF files
docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
  bcftools view -h /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz | grep "^##contig" | head -3
```

**Key differences:**

| Feature | GRCh37 / hg19 | GRCh38 / hg38 |
|---|---|---|
| chr1 length | 249,250,621 | 248,956,422 |
| Chromosome prefix | Often no `chr` | Always `chr` |
| Mitochondria name | `MT` (or `chrM`) | `chrM` |
| ALT contigs | No | Yes |

**Fix:** Extract FASTQ from BAM and re-align to GRCh38:
```bash
docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
  bash -c "samtools sort -n /genome/${SAMPLE}/aligned/old_hg19.bam | \
           samtools fastq -1 /genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz \
                          -2 /genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz -"

# Then run alignment (step 2) which uses GRCh38
./scripts/02-alignment.sh $SAMPLE
```

Do **not** use LiftOver on BAM files. Re-alignment from FASTQ is cleaner and avoids coordinate translation artifacts.

---

### CRAM files instead of BAM

**Symptom:** Pipeline scripts expect `.bam` files but you have `.cram` files.

**Fix:** Convert CRAM to BAM (requires the reference genome used for encoding, which is usually GRCh38):
```bash
docker run --rm --user root \
  -v ${GENOME_DIR}:/genome \
  staphb/samtools:1.20 \
  samtools view -b \
    -T /genome/reference/Homo_sapiens_assembly38.fasta \
    -o /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    /genome/${SAMPLE}/aligned/${SAMPLE}.cram

docker run --rm --user root \
  -v ${GENOME_DIR}:/genome \
  staphb/samtools:1.20 \
  samtools index /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
```

**If you get an error about mismatched references:** The CRAM was encoded against a different reference build. You need the exact FASTA used during encoding, or extract FASTQ and re-align:
```bash
docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
  bash -c "samtools sort -n /genome/${SAMPLE}/aligned/${SAMPLE}.cram | \
           samtools fastq -1 /genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz \
                          -2 /genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz -"
```

---

### Paired FASTQ naming issues

**Symptom:** Step 2 (alignment) fails with "ERROR: File not found" for R1 or R2.

**Cause:** The pipeline expects FASTQ files named exactly:
```
${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz
${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz
```

Vendors use different naming conventions:
```
# Illumina standard
Sample_S1_L001_R1_001.fastq.gz / Sample_S1_L001_R2_001.fastq.gz

# BGI/DNBSEQ
V350012345_L01_1.fq.gz / V350012345_L01_2.fq.gz

# Other
sample.1.fastq.gz / sample.2.fastq.gz
```

**Fix:** Create symlinks or rename to match the expected pattern:
```bash
cd ${GENOME_DIR}/${SAMPLE}/fastq/

# If you have multiple lanes, concatenate first:
cat *_R1_*.fastq.gz > ${SAMPLE}_R1.fastq.gz
cat *_R2_*.fastq.gz > ${SAMPLE}_R2.fastq.gz

# Or create symlinks for single files:
ln -s original_name_R1.fastq.gz ${SAMPLE}_R1.fastq.gz
ln -s original_name_R2.fastq.gz ${SAMPLE}_R2.fastq.gz
```

**Warning:** If concatenating multi-lane FASTQ, ensure R1 and R2 files are concatenated in the same lane order. Mismatched pairs will produce corrupt alignments.

---

### Corrupt or truncated downloads

**Symptom:** Tools crash with "unexpected end of file", "truncated file", or "not in gzip format".

**How to detect:**
```bash
# Test gzip integrity (for FASTQ.gz, VCF.gz)
gzip -t ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz
# If corrupt: "unexpected end of file" or "invalid compressed data"

# Test BAM integrity
docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
  samtools quickcheck /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
# No output = OK. Error message = corrupt.

# Check VCF can be read
docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
  bcftools view -h /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz > /dev/null
```

**Fix:** Re-download the file. For large files, use `wget -c` (supports resume):
```bash
wget -c https://url/to/your/data.fastq.gz
```

For reference data downloads that fail repeatedly, see the per-step sections for [VEP cache](#vep-cache-download-failures) and PCGR data bundle.

---

### BAM not sorted or not indexed

**Symptom:** Tools fail with "file is not sorted", "index file not found", or "[E::hts_idx_load3] Could not load index".

**How to detect:**
```bash
# Check if BAM is sorted
docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
  samtools view -H /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam | grep "^@HD"
# Should show: SO:coordinate

# Check if index exists
ls -la ${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai
```

**Fix:**
```bash
# Sort the BAM (if not already sorted)
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  staphb/samtools:1.20 \
  samtools sort -@ 4 \
    -o /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    /genome/${SAMPLE}/aligned/${SAMPLE}_unsorted.bam

# Create index (always needed)
docker run --rm --user root \
  -v ${GENOME_DIR}:/genome \
  staphb/samtools:1.20 \
  samtools index /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
```

**Note:** The BAM index file must have the exact same base name as the BAM. If your BAM is `sample_sorted.bam`, the index must be `sample_sorted.bam.bai` (not `sample_sorted.bai`).

---

## Per-Step Troubleshooting

### Step 2: Alignment (minimap2 + samtools)

**Problem: "Cannot find FASTQ" but files exist**
The script looks for `${SAMPLE}_R1.fastq.gz` and `${SAMPLE}_R2.fastq.gz`. See [Paired FASTQ naming issues](#paired-fastq-naming-issues).

**Problem: Alignment produces 0-byte BAM**
Usually means minimap2 failed silently. Run with `bash -x` to see the actual error. Common causes:
- Wrong reference genome path
- FASTQ files are corrupt (run `gzip -t` on both)
- Docker `--memory` too low (minimap2 needs ~6-10 GB for the GRCh38 index)

**Problem: minimap2 index build takes forever**
The `.mmi` index build is a one-time step (~30 minutes). If it seems stuck, check that the reference FASTA is not corrupt and the output path is writable.

---

### Step 3: DeepVariant crashes on Mac (amd64 emulation)

**Symptom:** DeepVariant container exits with code 137 or crashes with memory errors on Apple Silicon Mac.

**Cause:** DeepVariant is an amd64-only image. On Apple Silicon Macs, Docker runs it through Rosetta 2 emulation, which adds ~2x memory overhead and 3-5x CPU overhead. The default `--memory 32g` and `--cpus 8` settings are too aggressive for emulation.

**Fix:**
```bash
# Edit scripts/03-deepvariant.sh or run manually with reduced resources:
docker run --rm \
  --cpus 2 --memory 12g \
  -v "${GENOME_DIR}:/genome" \
  google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="/genome/reference/Homo_sapiens_assembly38.fasta" \
    --reads="/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    --output_vcf="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --num_shards=2
```

**Expected runtime on Mac:** 8-16 hours (vs 2-4 hours on native Linux amd64).

**Alternative:** If DeepVariant repeatedly crashes on your Mac, run ONLY step 3 on a cloud Linux instance (e.g., Hetzner CCX33 for ~$0.18/hr) and download the VCF. All other steps can proceed on Mac.

---

### Step 4: Manta produces no variants

**Symptom:** `diploidSV.vcf.gz` is empty or has only header lines.

**Common causes:**
1. BAM has too few reads (low-coverage WGS below ~10X)
2. BAM was not indexed (Manta requires `.bai`)
3. Reference FASTA mismatch — chromosome names must match between BAM and FASTA

**Diagnosis:**
```bash
# Check BAM read count
docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
  samtools flagstat /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
# "mapped" count should be 500M+ for 30X WGS

# Check chromosome naming consistency
docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
  samtools view -H /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam | grep "^@SQ" | head -1
# Must use "chr" prefix (SN:chr1) to match GRCh38 reference
```

---

### Step 5: AnnotSV produces no output

**Symptom:** AnnotSV runs but the output TSV is empty or missing.

**Common causes:**
1. **Input VCF path wrong:** AnnotSV expects the Manta `diploidSV.vcf.gz`. The script checks both `manta/` and `manta2/` directories.
2. **No PASS variants:** If Manta flagged all SVs as filtered (no PASS), AnnotSV may produce empty output.
3. **Image version mismatch:** Some AnnotSV versions require specific internal database versions.

**Fix:**
```bash
# Verify the Manta VCF has PASS variants
docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
  bcftools view -f PASS /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz | grep -c -v "^#"
# Should be > 0. Typical: 5,000-9,000 for 30X WGS.
```

If the Manta VCF is valid but AnnotSV still produces nothing, try running with verbose output:
```bash
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  getwilds/annotsv:latest \
  AnnotSV \
    -SVinputFile "/genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz" \
    -outputFile "/genome/${SAMPLE}/annotsv/${SAMPLE}_test.tsv" \
    -genomeBuild GRCh38 \
    -annotationMode both 2>&1 | tee annotsv_debug.log
```

---

### Step 7: PharmCAT VCF format requirements

**Symptom:** PharmCAT exits with "Input VCF does not meet requirements" or produces empty report with no gene calls.

**Cause:** PharmCAT is strict about VCF format:
- Must be block-gzipped (`.vcf.gz`, not plain `.vcf`)
- Must have a tabix index (`.vcf.gz.tbi`)
- Must be aligned to GRCh38
- Must contain a GT (genotype) FORMAT field
- Multi-allelic sites must be decomposed

**Fix:**
```bash
# Verify your VCF has the GT field
docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
  bcftools query -f '[%GT]\n' /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz | head -1
# Should show something like "0/1" or "1/1"

# If your VCF is plain text, compress and index it:
docker run --rm --user root -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
  bash -c "bcftools view /genome/${SAMPLE}/vcf/${SAMPLE}.vcf -Oz \
    -o /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"

# If multi-allelic sites are an issue, normalize:
docker run --rm --user root -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
  bcftools norm -m -both \
    -f /genome/reference/Homo_sapiens_assembly38.fasta \
    /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    -Oz -o /genome/${SAMPLE}/vcf/${SAMPLE}_norm.vcf.gz
```

**Note:** PharmCAT includes a VCF preprocessor. If direct input fails, try the preprocessor first:
```bash
docker run --rm \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  -v "${GENOME_DIR}/reference:/ref" \
  pgkb/pharmcat:2.15.5 \
  java -cp /pharmcat/pharmcat.jar org.pharmgkb.pharmcat.VcfPreprocessor \
    -vcf "/data/${SAMPLE}.vcf.gz" \
    -refFasta /ref/Homo_sapiens_assembly38.fasta
```

---

### Step 9: ExpansionHunter fails immediately

**Symptom:** Container exits with "the option '--log' is required but missing".

**Cause:** The `weisburd/expansionhunter:latest` image (v2.5.5) requires the `--log` parameter, which is not optional.

**Fix:** Ensure the command includes `--log`:
```bash
--log "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.log"
```
All pipeline scripts already include this. If running manually, do not omit it.

**Symptom:** ExpansionHunter uses `--variant-catalog` but the image expects `--repeat-specs`.

**Cause:** Version mismatch. The `weisburd/expansionhunter` image uses v2.5.5 which expects `--repeat-specs` pointing to a directory, not `--variant-catalog` pointing to a single JSON file (v5.x syntax).

**Fix:** Use the correct syntax for the image:
```bash
--repeat-specs /pathogenic_repeats/GRCh38/    # v2.5.5 (weisburd image)
# NOT: --variant-catalog /path/to/catalog.json  # v5.x syntax
```

---

### Step 13: VEP cache download failures

**Symptom:** VEP cache download times out, produces a partial file, or extraction fails.

**Cause:** The VEP cache is ~22-26 GB, hosted on Ensembl FTP servers that can be slow. The built-in `INSTALL.pl` downloader does not support resume and may fail silently.

**Fix — manual download with resume:**
```bash
mkdir -p ${GENOME_DIR}/vep_cache/tmp
cd ${GENOME_DIR}/vep_cache/tmp

# wget -c resumes interrupted downloads
wget -c https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz

# Verify download size (~22 GB)
ls -lh homo_sapiens_vep_112_GRCh38.tar.gz

# Extract (takes 10-20 minutes, expands to ~30 GB)
cd ${GENOME_DIR}/vep_cache
tar xzf tmp/homo_sapiens_vep_112_GRCh38.tar.gz

# Verify extraction
ls ${GENOME_DIR}/vep_cache/homo_sapiens/112_GRCh38/
# Should contain info.txt, variation_set_*.gz, and many other files
```

**Symptom:** VEP INSTALL.pl fails with "Cannot open Local file" permission error.

**Cause:** The VEP container runs as a non-root user who cannot write to `/opt/vep/.vep/tmp/`.

**Fix:** Use `--user root` when running VEP and pre-create the temp directory. Or (recommended) download the cache manually as shown above.

---

### Step 17: CPSR `--pcgr_dir` path confusion

**Symptom:** `Data directory (/genome/pcgr_data/data/data) does not exist`

**Cause:** CPSR internally appends `/data` to whatever path you provide as `--pcgr_dir`. If you point it to the `data/` subdirectory, it looks for `data/data/`.

**Fix:** Point `--pcgr_dir` to the **parent** of the `data/` directory:
```bash
# CORRECT:
--pcgr_dir /genome/pcgr_data
# This makes CPSR look for: /genome/pcgr_data/data/grch38/ (exists)

# WRONG:
--pcgr_dir /genome/pcgr_data/data
# This makes CPSR look for: /genome/pcgr_data/data/data/ (does NOT exist)
```

**Additional CPSR issue — wrong Docker image:**
CPSR does not have its own Docker image. It is bundled inside the PCGR image:
```bash
# CORRECT:
docker run ... sigven/pcgr:1.4.1 cpsr ...

# WRONG (image does not exist):
docker run ... sigven/cpsr:2.0.0 ...
```

---

### Step 18: CNVnator empty ROOT file

**Symptom:** CNVnator steps 2-5 fail because the ROOT file from step 1 is empty or corrupt. Or the final output has 0 CNV calls.

**Common causes:**
1. **BAM path wrong inside container.** The path must be the containerized path (`/genome/...`), not the host path.
2. **ROOT file not writable.** Add `--user root`.
3. **BAM chromosome naming mismatch.** CNVnator needs chromosomes matching the reference FASTA.

**Diagnosis:**
```bash
# Check ROOT file size (should be several GB for 30X WGS)
ls -lh ${GENOME_DIR}/${SAMPLE}/cnvnator/${SAMPLE}.root
# If < 1 MB, step 1 (-tree) failed silently.

# Verify the image tag exists before running
docker pull quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11
```

**Fix:** Re-run step 1 with verbose output and confirm the BAM path resolves correctly inside the container.

---

### Step 19: Delly long runtime / high memory

**Symptom:** Delly has been running for 8+ hours on a 30X genome.

**Cause:** Delly examines every read pair across the genome for split-read and paired-end SV evidence. On high-coverage samples or with slow storage, this takes longer.

**Expected runtimes:**

| Coverage | CPU Cores | Expected Runtime |
|---|---|---|
| 30X | 4 | 2-4 hours |
| 30X | 8 | 1.5-3 hours |
| 50X+ | 4 | 4-8 hours |

**If Delly exceeds 12 hours:**
1. Check `docker stats` to verify it is still using CPU (not stuck)
2. Confirm storage is not a bottleneck (NFS/SMB mounts are much slower)
3. Consider running only on specific chromosomes to reduce scope:
   ```bash
   # Run Delly on chr1-chr22 only (skip ALT contigs)
   docker run --rm --user root \
     --cpus 4 --memory 8g \
     -v "${GENOME_DIR}:/genome" \
     quay.io/biocontainers/delly:1.7.3--hd6466ae_0 \
     delly call \
       -g /genome/reference/Homo_sapiens_assembly38.fasta \
       -o /genome/${SAMPLE}/delly/${SAMPLE}_sv.bcf \
       /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
   # Delly processes all standard chromosomes by default; slowness usually
   # comes from ALT contigs. You can exclude them with the -x (exclude) flag
   # using a BED file of regions to skip.
   ```

---

### Step 20: GATK Mutect2 mitochondrial mode

**Problem: "Cannot create sequence dictionary"**
The script tries to create `Homo_sapiens_assembly38.dict` if it does not exist. If you get permission errors, create it manually:
```bash
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  broadinstitute/gatk:4.6.1.0 \
  gatk CreateSequenceDictionary \
    -R /genome/reference/Homo_sapiens_assembly38.fasta \
    -O /genome/reference/Homo_sapiens_assembly38.dict
```

**Problem: 0 mitochondrial variants called**
Possible causes:
- BAM has no chrM reads (check with `samtools idxstats | grep chrM`)
- Mitochondrial chromosome name mismatch (`MT` vs `chrM`). GRCh38 uses `chrM`.

---

### HLA Typing (Step 8): Known difficulties

HLA typing from WGS data is unreliable in Docker. The two main tools have unresolved issues:

**HLA-LA:** Requires a pre-serialized 40 GB graph. Most Docker images do not include it. The one image that does (`jiachenzdocker/hla-la:latest`) crashes during graph alignment. **Status: UNSOLVED.**

**T1K:** Works partially, but the coordinate file build step requires the full reference FASTA (not the `.fai` index). With correctly built coordinates, T1K can call some HLA alleles but may have ~50% unmapped alleles.

**Recommendation:** For clinical HLA typing, rely on dedicated lab assays. Short-read WGS is not ideal for this highly polymorphic region.

---

## Performance Issues

### Everything is slow on macOS Apple Silicon

**Cause:** All bioinformatics Docker images are amd64/x86_64. On Apple Silicon Macs (M1-M4), Docker Desktop runs them through Rosetta 2 emulation. This adds:
- **2-5x CPU overhead** (emulation is compute-expensive)
- **~2x memory overhead** (emulation runtime consumes extra RAM)
- **I/O overhead** from the Docker Desktop VM

**Expected slowdown by step:**

| Step | Native Linux | Mac (Rosetta 2) | Slowdown |
|---|---|---|---|
| DeepVariant | 2-4 hr | 8-16 hr | 3-5x |
| minimap2 alignment | 1-2 hr | 3-6 hr | 3x |
| VEP annotation | 2-4 hr | 4-8 hr | 2x |
| Manta | 20 min | 1-2 hr | 3-4x |
| CNVnator | 2-4 hr | 6-12 hr | 3x |
| bcftools steps | 1-5 min | 2-10 min | 2x |

**Mitigation strategies:**
1. **Reduce parallelism.** Run one heavy step at a time instead of many in parallel.
2. **Reduce resource limits.** Use `--cpus 2 --memory 8g` instead of `--cpus 8 --memory 32g` to avoid emulation thrashing.
3. **Offload the heaviest steps.** Run steps 2, 3, 18, 19 on a remote Linux machine and bring back the outputs. Everything after step 3 needs only the VCF or BAM.
4. **Use a cloud instance.** A Hetzner CCX33 (8 vCPU, 32 GB, ~$0.18/hr) will run the full pipeline in 6-10 hours for about $1-2.

---

### WSL2 file I/O is extremely slow

**Symptom:** Steps that read large files (alignment, DeepVariant, VEP) take 10-50x longer than expected.

**Cause:** If your data is on a Windows drive (`/mnt/c/`, `/mnt/d/`), WSL2 accesses it through the 9P filesystem protocol, which is catastrophically slow for large random-access files.

**Fix:** Move ALL genomics data to the Linux filesystem:
```bash
# Inside WSL2
mkdir -p ~/genome_data
# Move or copy data from /mnt/c/ to ~/genome_data/
cp -r /mnt/c/Users/you/genome_data/* ~/genome_data/
export GENOME_DIR=~/genome_data
```

Data on the Linux filesystem (`~/`, `/home/`, `/tmp/`) uses ext4 natively and performs at full speed.

---

### Monitoring container resources

Use `docker stats` to see real-time resource usage:

```bash
# Live dashboard (updates every second)
docker stats

# One-time snapshot
docker stats --no-stream

# Output columns:
# CONTAINER  CPU%  MEM USAGE/LIMIT  NET I/O  BLOCK I/O
```

**What to look for:**
- **CPU at 0%** for a long time: Container may be stuck or waiting for I/O
- **MEM USAGE near LIMIT**: Container is about to be OOM-killed. Increase `--memory`.
- **NET I/O > 0**: Container is downloading something (may be slow due to network)
- **BLOCK I/O very high**: I/O-bound step. Faster storage would help.

---

### When to reduce `--cpus` or `--memory`

**Reduce `--cpus` when:**
- Running multiple steps in parallel on a machine with limited cores
- Running on macOS with Rosetta 2 emulation (more shards = more overhead)
- Your machine has fewer than 8 cores

**Reduce `--memory` when:**
- Docker Desktop memory allocation is less than the script's `--memory` flag
- Running multiple containers simultaneously
- System becomes unresponsive during analysis

**Safe minimum values per step:**

| Step | Min --cpus | Min --memory | Notes |
|---|---|---|---|
| 2 (minimap2) | 2 | 8g | Slower but works |
| 3 (DeepVariant) | 2 | 8g | Set `--num_shards=2` to match |
| 4 (Manta) | 2 | 4g | |
| 6 (ClinVar) | 1 | 1g | Very light |
| 7 (PharmCAT) | 1 | 2g | |
| 9 (ExpansionHunter) | 2 | 2g | |
| 13 (VEP) | 2 | 4g | Set `--fork 2` to match |
| 17 (CPSR) | 2 | 4g | |
| 18 (CNVnator) | 2 | 4g | |
| 19 (Delly) | 2 | 4g | |

---

### Steps that can run in parallel vs must be sequential

```
SEQUENTIAL (must run in order):
  Step 2 (alignment) ──> Step 3 (variant calling)
  Step 4 (Manta) ──> Step 5 (AnnotSV)
  Step 4 (Manta) ──> Step 15 (duphold)

PARALLEL after Step 3 completes (all independent):
  ┌─ Step 4 (Manta)         ← needs BAM
  ├─ Step 6 (ClinVar)       ← needs VCF
  ├─ Step 7 (PharmCAT)      ← needs VCF
  ├─ Step 9 (ExpansionHunter) ← needs BAM
  ├─ Step 10 (TelomereHunter) ← needs BAM
  ├─ Step 11 (ROH)          ← needs VCF
  ├─ Step 12 (haplogrep3)   ← needs VCF
  ├─ Step 13 (VEP)          ← needs VCF
  ├─ Step 14 (imputation)   ← needs VCF
  ├─ Step 16 (indexcov)      ← needs BAM index only
  ├─ Step 17 (CPSR)         ← needs VCF
  ├─ Step 18 (CNVnator)     ← needs BAM
  ├─ Step 19 (Delly)        ← needs BAM
  └─ Step 20 (Mutect2 mito) ← needs BAM

RESOURCE CAUTION:
  Steps 3, 18, 19 are the most CPU/RAM intensive.
  Do not run more than 2 of these simultaneously on a 32 GB machine.
```

---

## Output Issues

### 0-byte output files

**Symptom:** Output file exists but is 0 bytes.

**Most common cause:** The input path was wrong inside the Docker container. Docker mounts the host path to `/genome`, so all paths inside the container must start with `/genome/`, not the host path.

**Diagnosis checklist:**
1. **Check that the input file exists at the expected path:**
   ```bash
   ls -la ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz
   ```
2. **Run the script with debug output:**
   ```bash
   bash -x scripts/06-clinvar-screen.sh your_name
   ```
3. **Inspect Docker mount mapping.** The script mounts `${GENOME_DIR}:/genome`. Inside the container, `${GENOME_DIR}/sample/vcf/sample.vcf.gz` becomes `/genome/sample/vcf/sample.vcf.gz`.
4. **Check for bgzip/tabix path issues.** The `staphb/bcftools:1.21` image does not include `bgzip` or `tabix` in `$PATH`. Use `bcftools view -Oz -o` instead of piping to `bgzip`. See [lessons-learned.md](lessons-learned.md#bgzip-tabix-not-in-bcftools-image-path).

**Other causes:**
- Tool crashed silently before writing output (check exit code and stderr)
- Disk full (check `df -h`)
- Permission denied in output directory (add `--user root`)

---

### Wrong number of variants (too few)

**Symptom:** DeepVariant VCF has far fewer variants than expected.

**Expected variant counts for 30X WGS (GRCh38):**

| Metric | Expected Range | Notes |
|---|---|---|
| Total variants (SNP + indel) | 4.5-5.5 million | Total rows in VCF |
| PASS variants | 4.0-4.8 million | After DeepVariant quality filtering |
| SNPs | 3.8-4.5 million | ~85% of total |
| Indels | 600K-900K | ~15% of total |
| Heterozygous | 2.5-3.2 million | ~55-60% of total |
| Homozygous ALT | 1.8-2.3 million | ~40-45% of total |

**If total variants < 3 million:**
1. **Genome build mismatch.** Your BAM may be hg19, but the reference is hg38. See [Wrong genome build](#wrong-genome-build-hg19-vs-hg38).
2. **Low coverage.** Check average depth:
   ```bash
   docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
     samtools depth -a /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam | \
     awk '{sum+=$3; n++} END {print "Average depth:", sum/n}'
   ```
   30X WGS should show average depth of 25-35.
3. **BAM is mostly unmapped.** Check alignment rate:
   ```bash
   docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
     samtools flagstat /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
   ```
   Mapped percentage should be > 98%.

**If total variants > 7 million:**
1. Likely includes many false positives. Filter to PASS only:
   ```bash
   docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
     bcftools view -f PASS /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz | \
     bcftools stats | grep "number of records"
   ```
2. Multi-allelic sites may be inflating the count. Normalize:
   ```bash
   docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
     bcftools norm -m -both \
       -f /genome/reference/Homo_sapiens_assembly38.fasta \
       /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz | bcftools stats | grep "number of records"
   ```

---

### Wrong number of variants (too many ClinVar hits)

**Symptom:** ClinVar screen returns hundreds of "pathogenic" hits, which seems too many.

**Cause:** The ClinVar pathogenic filter includes `CLNSIG~"Pathogenic"`, which matches variants where pathogenicity depends on context. Some variants are "Pathogenic" for one condition but "Benign" for another (compound `CLNSIG` fields like `Pathogenic/Likely_benign`).

**Expected ClinVar hit counts:**

| Category | Expected Count |
|---|---|
| Total ClinVar matches (any significance) | 50-150 |
| Pathogenic + Likely pathogenic | 0-10 |
| True actionable findings | 0-3 |

**If you see more than ~15 pathogenic hits:** The ClinVar VCF may not be filtered correctly, or chromosome naming may be mismatched (causing position-only matches without allele verification). Re-run the ClinVar screen from [step 0 setup](00-reference-setup.md).

---

### How to verify outputs are correct

Quick sanity checks for each major output:

```bash
SAMPLE=your_name

# VCF: check variant count and type distribution
docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
  bcftools stats /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz | grep "^SN"

# BAM: check alignment rate and depth
docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.20 \
  samtools flagstat /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam

# Manta SV count (expect 5,000-9,000 total)
docker run --rm -v ${GENOME_DIR}:/genome staphb/bcftools:1.21 \
  bcftools view /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz | grep -c -v "^#"

# PharmCAT: check the HTML report exists and has content
ls -la ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.report.html
# Should be > 50 KB

# CPSR: check the HTML report exists
ls -la ${GENOME_DIR}/${SAMPLE}/cpsr/${SAMPLE}.cpsr.grch38.html
# Should be > 100 KB

# VEP: check annotated VCF exists and has annotations
head -50 ${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf | grep "CSQ="
# Should see consequence annotations
```

---

### Expected output sizes per step

| Step | Output File(s) | Expected Size |
|---|---|---|
| 2 | `*_sorted.bam` + `.bai` | 80-120 GB + 5-8 MB |
| 3 | `*.vcf.gz` + `.tbi` | 80-200 MB + 1-2 MB |
| 4 | `diploidSV.vcf.gz` | 1-5 MB |
| 5 | `*_sv_annotated.tsv` | 25-35 MB |
| 6 | `isec/0002.vcf` (shared hits) | < 1 MB |
| 7 | `*.report.html` | 50-200 KB |
| 9 | `*_eh.json` + `*_eh.vcf` | < 1 MB each |
| 10 | `*_summary.tsv` + plots | 50-200 MB total |
| 11 | `*_roh.txt` | 1-5 MB |
| 12 | `*_haplogroup.txt` | < 1 KB |
| 13 | `*_vep.vcf` (uncompressed) | 2-5 GB |
| 17 | `*.cpsr.grch38.html` + TSV | 50-200 MB total |
| 18 | `*_cnvs.txt` + `.root` | Text: < 1 MB, ROOT: 5-15 GB |
| 19 | `*_sv.vcf.gz` | 5-20 MB |
| 20 | `*_chrM_filtered.vcf.gz` | < 1 MB |

If any output is significantly smaller than expected (especially 0 bytes), see [0-byte output files](#0-byte-output-files).

---

## Getting Help

### Before opening an issue

1. **Check this guide first.** Ctrl+F your error message.
2. **Check [docs/lessons-learned.md](lessons-learned.md)** for known Docker image issues and tool quirks.
3. **Run with debug output:** `bash -x scripts/XX-stepname.sh your_name 2>&1 | tee debug.log`
4. **Collect diagnostic info:**
   ```bash
   docker version
   docker info | grep -E "Total Memory|CPUs|Storage Driver"
   uname -a
   df -h ${GENOME_DIR}
   ```

### Opening a GitHub issue

File an issue at: **[github.com/GeiserX/genomics-pipeline/issues](https://github.com/GeiserX/genomics-pipeline/issues)**

Include:
- Which step failed (step number and script name)
- Full error message (copy-paste, not screenshot)
- Your platform (Linux distro, macOS version, WSL2 version)
- Docker version (`docker version`)
- Available RAM and disk space
- Input data type (FASTQ, BAM, or VCF) and approximate size
- Whether you are using native Linux or Docker Desktop

### Useful diagnostic commands

```bash
# System info
uname -a                           # OS and architecture
docker version                     # Docker version
docker info | head -30             # Docker config summary
docker system df                   # Docker disk usage

# Container debugging
docker ps -a --latest              # Last container status
docker logs <container_id>         # Container output
docker inspect <container_id>      # Full container metadata

# File validation
samtools quickcheck file.bam       # BAM integrity (via Docker)
gzip -t file.vcf.gz               # Gzip integrity
bcftools view -h file.vcf.gz      # VCF header readability

# Resource monitoring
docker stats --no-stream           # Current container resource usage
df -h                              # Disk space
free -h                            # RAM (Linux only)
```
