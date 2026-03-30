# Step 22: Structural Variant Consensus Merge

> **EXPERIMENTAL:** This step uses a heuristic position-binning approach that may over-count calls from the same caller. Results should be treated as a rough intersection, not a true consensus merge. For production use, consider SURVIVOR or Jasmine with proper multi-sample VCF merging.

## What This Does

Performs a rough intersection of structural variant (SV) calls from multiple independent callers — Manta (step 4), Delly (step 19), and CNVnator (step 18). SVs are binned by chromosome, position (1 kb windows), and SV type; bins with calls from two or more callers are retained. This is an approximation, not a true breakpoint-aware merge like SURVIVOR or Jasmine would produce.

## Why

Individual SV callers each have distinct biases and false-positive profiles:

- **Manta**: Fast, sensitive for smaller SVs and indels (paired-end + split-read)
- **Delly**: Strongest for inversions and balanced translocations (paired-end + split-read + depth)
- **CNVnator**: Best for large CNVs (read-depth only)

Taking the intersection across callers reduces false positives. An SV seen by two independent algorithms using different signal types is more likely to be real. Note that dedicated SV comparison tools (SURVIVOR, Jasmine) use breakpoint distance, size similarity, and strand matching for more accurate merging than the position-binning heuristic used here.

## Tool

- **bcftools** (for merging and overlap detection)

The script uses a breakpoint-binning approach with bcftools rather than SURVIVOR, since SURVIVOR Docker image availability is unreliable. SVs are grouped by chromosome, binned position (1 kb windows), and SV type. Bins with calls from 2+ callers are kept.

## Docker Image

```
staphb/bcftools:1.21
```

## Input

At least two of the following (the script auto-detects which are available):

| Caller | Expected path |
|---|---|
| Manta (step 4) | `${GENOME_DIR}/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz` |
| Delly (step 19) | `${GENOME_DIR}/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz` |
| CNVnator (step 18) | `${GENOME_DIR}/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.vcf.gz` or `_cnvs.txt` |

If CNVnator output is in TXT format (its native output), the script automatically converts it to VCF before merging.

## Command

```bash
./scripts/22-survivor-merge.sh your_name
```

## What the Script Does Internally

1. Scans for available SV VCFs from Manta, Delly, and CNVnator
2. If CNVnator output is only in TXT format, converts it to VCF (adding proper headers, SV type, and END coordinates)
3. Requires at least 2 callers to proceed (exits with an error otherwise)
4. Extracts PASS variants from each caller and bins them by `chromosome + position/1000 + SVTYPE`
5. Keeps bins where 2+ callers contributed a call (consensus SVs)
6. Writes a sorted, compressed, and indexed consensus VCF
7. Reports the total count of consensus SVs

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_sv_consensus.vcf.gz` | Consensus SVs called by 2+ callers |
| `${SAMPLE}_sv_consensus.vcf.gz.tbi` | Tabix index |
| `sv_files.txt` | List of input VCFs used |
| `consensus_raw.txt` | Intermediate merged records |

All output is written to `${GENOME_DIR}/${SAMPLE}/sv_merged/`.

## Runtime

~5-15 minutes (mostly I/O reading the input VCFs).

## Interpreting Results

A typical 30X WGS genome produces:

- **Manta**: 3,000-5,000 SVs
- **Delly**: 5,000-15,000 SVs
- **CNVnator**: 500-2,000 CNVs

After consensus filtering, expect **200-1,000 multi-caller SVs**. These have lower false-positive rates than single-caller calls, though the 1 kb binning heuristic is less precise than dedicated tools like SURVIVOR or Jasmine.

SV types in the output:
- **DEL** -- Deletion (missing segment)
- **DUP** -- Duplication (extra copy of a segment)
- **INV** -- Inversion (segment flipped in orientation)
- **BND** -- Breakend / Translocation (segment moved to another chromosome)
- **INS** -- Insertion

### Quick inspection

```bash
# Count consensus SVs by type
docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
  bcftools query -f '%INFO/SVTYPE\n' \
    /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz | sort | uniq -c | sort -rn
```

## Limitations

- The 1 kb breakpoint-binning approach is an approximation. True SURVIVOR merge uses more sophisticated overlap criteria (breakpoint distance, SV type matching, strand, size similarity). Some near-boundary SVs may be missed or incorrectly grouped.
- CNVnator-to-VCF conversion produces minimal VCF records (no genotype, no quality scores). These SVs carry less metadata than Manta/Delly calls.
- Single-caller SVs are discarded even if they are real. If you suspect a specific SV, check the individual caller outputs directly.
- BND (translocation) breakpoints from different callers may not bin together well due to how breakends are represented.

## Notes

- Run this step only after completing at least two of: step 4 (Manta), step 19 (Delly), step 18 (CNVnator).
- All three callers are independent of each other and can run in parallel after alignment.
- The consensus VCF can be annotated with VEP or loaded into IGV for visual inspection.
- For clinical-grade SV analysis, consider also running AnnotSV on the consensus set.

## Links

- [SURVIVOR (original tool)](https://github.com/fritzsedlazeck/SURVIVOR) -- the gold-standard SV merge tool, used as conceptual basis
- [Manta](https://github.com/Illumina/manta)
- [Delly](https://github.com/dellytools/delly)
- [CNVnator](https://github.com/abyzovlab/CNVnator)
