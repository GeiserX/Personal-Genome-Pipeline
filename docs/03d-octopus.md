# Step 3d: Octopus Variant Caller (Alternative)

Haplotype-aware Bayesian variant caller. Optional 5th caller for benchmarking or potential FreeBayes replacement.

---

## What It Does

Octopus jointly calls SNPs and indels using a haplotype-based model that considers multiple candidates simultaneously. Unlike GATK's fixed-ploidy model, Octopus uses a variable-ploidy approach that handles complex regions more accurately.

## Why

- **Better than FreeBayes in benchmarks** — higher precision with comparable sensitivity, and significantly faster
- **Lower memory than FreeBayes** — typically < 8 GB vs FreeBayes' 13+ GB peaks
- **Multithreaded** — unlike FreeBayes which is single-threaded
- **Haplotype-aware** — considers multiple variant candidates in a window, reducing false calls from alignment artifacts

## Tool

**Octopus** v0.7.4 — Bayesian haplotype-based mutation calling.

- Paper: Cooke et al., Nature Biotechnology 2021 (doi:10.1038/s41587-021-00861-3)
- Source: [github.com/luntergroup/octopus](https://github.com/luntergroup/octopus)

## Docker Image

```
dancooke/octopus:0.7.4
```

## Command

```bash
export GENOME_DIR=/path/to/data
./scripts/03d-octopus.sh <sample_name>

# Test on one chromosome first:
INTERVALS=chr22 ./scripts/03d-octopus.sh <sample_name>
```

## Output

| File | Location | Description |
|---|---|---|
| VCF | `vcf_octopus/<sample>.vcf.gz` | Germline variant calls |

## Runtime

| Dataset | Threads | Time | Memory |
|---|---|---|---|
| 30X WGS | 8 | 2-4 hours | ~8-12 GB |
| chr22 only | 4 | ~5-10 min | < 4 GB |

## Notes

- Outputs to `vcf_octopus/` to maintain isolation from the default DeepVariant calls
- Supports `ALIGN_DIR` variable for BWA-MEM2 alignments: `ALIGN_DIR=aligned_bwamem2 ./scripts/03d-octopus.sh sample`
- Supports `INTERVALS` for region-restricted testing: `INTERVALS=chr22 ./scripts/03d-octopus.sh sample`
- The benchmark script (`benchmark-variants.sh`) auto-discovers `vcf_octopus/` for pairwise concordance analysis
- Octopus is a good candidate to replace FreeBayes (which is slow, single-threaded, and memory-heavy) as the default alternative caller
