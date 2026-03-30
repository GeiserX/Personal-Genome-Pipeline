# Step 1: ORA to FASTQ Conversion

## What This Does
Converts Illumina's proprietary ORA-compressed sequencing files into standard FASTQ format. ORA is Illumina's lossless compression format (~3x smaller than gzip FASTQ).

## Why
Raw sequencing data from Illumina DRAGEN comes in ORA format. All downstream tools expect FASTQ.

## Tool
- **orad** (Illumina ORA decompression tool)
- Not available as Docker image — must be installed natively

## Prerequisites
- `orad` binary (download from Illumina)
- Sufficient disk space: ORA→FASTQ expands ~3x (e.g., 30GB ORA → 90GB FASTQ)

## Command
```bash
# Convert paired-end ORA files to FASTQ
orad --ora-reference /path/to/oradata/reference \
  sample_R1.ora sample_R2.ora \
  --output-directory /path/to/output/

# Output: sample_R1.fastq.gz, sample_R2.fastq.gz
```

## Notes
- ORA reference files must match the sequencing run (provided alongside ORA files)
- If you receive FASTQ.gz directly (e.g., from a resequencing service), skip this step
