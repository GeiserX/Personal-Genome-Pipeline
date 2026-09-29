# Prerequisites and platform notes

## Prerequisites

### Hardware Requirements

| Resource | Minimum | Recommended | Notes |
|---|---|---|---|
| **CPU** | 4 cores | 16+ cores | DeepVariant scales linearly with cores |
| **RAM** | 16 GB | 32 GB | Some steps need 8-16 GB; pipeline limits each container |
| **Disk** | 500 GB free | 1 TB+ | See [detailed breakdown](hardware-requirements.md) |
| **Internet** | Broadband | 100+ Mbps | ~70-75 GB core downloads + ~104 GB optional annotation databases |
| **OS** | Linux (amd64) | Ubuntu 22.04+ | macOS/ARM works but slower (see below) |

> **Disk space is the #1 surprise.** A single 30X WGS sample produces 60-90 GB of FASTQ, 30-80 GB of BAM, plus reference genomes and databases. See [docs/hardware-requirements.md](hardware-requirements.md) for the full breakdown.

### Software

| Software | Version | Install |
|---|---|---|
| Docker | 20.10+ | [docs.docker.com/get-docker](https://docs.docker.com/get-docker/) |
| bash | 4.0+ | Pre-installed on Linux; macOS ships 3.2 — install via `brew install bash` |
| wget or curl | Any | For downloading references |
| python3 *(optional)* | 3.6+ | Used by long-read alignment (02b) for symlink resolution. Falls back to `readlink -f` on GNU/Linux if absent |

That's it. Every analysis tool runs inside Docker -- no conda environments, no Python version conflicts, no compilation.

### Reference Data (One-Time Downloads)

| Resource | Size | Required For |
|---|---|---|
| GRCh38 reference FASTA + index | ~3.5 GB | All steps |
| ClinVar database | ~200 MB | Step 6 (ClinVar screen) |
| VEP cache | ~26 GB | Step 13 (VEP annotation) |
| PCGR/CPSR data bundle + VEP 113 cache | ~31 GB | Step 17 (cancer predisposition) |
| Docker images (all steps) | ~10-15 GB | All steps |
| Annotation databases (CADD, SpliceAI, REVEL, AlphaMissense) | ~104 GB | Steps 30-31 (optional) |
| **Total one-time setup (core)** | **~70-75 GB** | |
| **Total with annotation enrichment** | **~175 GB** | |

See [docs/00-reference-setup.md](00-reference-setup.md) for download instructions.

## Platform Notes

### Linux (Recommended)
Best performance. Docker runs natively. All pipeline images are linux/amd64. No issues.

### macOS (Intel)
Works fine. Docker Desktop runs a Linux VM, so there's a ~10-20% I/O overhead on file operations. Set Docker Desktop memory to at least 16 GB (Preferences > Resources).

### macOS (Apple Silicon / M1-M4)
Works but **slower**. All bioinformatics Docker images are amd64 and run under Rosetta 2 emulation (2-5x performance penalty). DeepVariant and BWA-MEM2 are the most affected. Set Docker Desktop to use Rosetta 2 for amd64 emulation (enabled by default on newer versions).

### Windows (WSL2)
Works. Install Docker Desktop with WSL2 backend. **Critical:** Keep all genomics data on the Linux filesystem (`~/data/`, not `/mnt/c/`). Accessing Windows drives from WSL2 is 10-50x slower due to the 9P protocol. Set WSL2 memory in `%UserProfile%\.wslconfig`:
```ini
[wsl2]
memory=24GB
swap=8GB
```

### Unraid / NAS Servers
Works great for long-running analyses. Use `--cpus` and `--memory` Docker flags (already set in all scripts) to avoid starving other services. Consider running in detached mode (`-d` flag) for multi-hour steps.

