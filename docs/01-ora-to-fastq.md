# Step 1: ORA to FASTQ Conversion

## What This Does
Converts Illumina's proprietary ORA-compressed sequencing files into standard FASTQ format. ORA is Illumina's lossless compression format (~5x smaller than gzipped FASTQ).

## Why
Raw sequencing data from Illumina DRAGEN comes in ORA format. All downstream tools expect FASTQ.

## Tool
- **orad** (Illumina ORA decompression tool)
- Not available as Docker image — must be installed natively

## Prerequisites
- `orad` binary (download from Illumina)
- Sufficient disk space: ORA→FASTQ.gz expands ~5x (e.g., a 30X genome is 15-20 GB of ORA and 60-90 GB of FASTQ.gz)

## Command

The script takes three arguments and decompresses one ORA file per call. Set `ORAD` if `orad` is not at `/opt/orad/bin/orad`.

```bash
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_name

# Arguments: <sample> <ora_reference_dir> <ora_file>; run once for R1 and once for R2
./scripts/01-ora-to-fastq.sh $SAMPLE /path/to/oradata /path/to/${SAMPLE}_S1_L001_R1_001.fastq.ora
./scripts/01-ora-to-fastq.sh $SAMPLE /path/to/oradata /path/to/${SAMPLE}_S1_L001_R2_001.fastq.ora
```

The output lands in `${GENOME_DIR}/${SAMPLE}/fastq/` and keeps the ORA file's name (`..._R1_001.fastq.gz`). Step 1b and step 2 read `${SAMPLE}_R1.fastq.gz` and `${SAMPLE}_R2.fastq.gz`, so rename the two files:

```bash
cd ${GENOME_DIR}/${SAMPLE}/fastq
mv ${SAMPLE}_S1_L001_R1_001.fastq.gz ${SAMPLE}_R1.fastq.gz
mv ${SAMPLE}_S1_L001_R2_001.fastq.gz ${SAMPLE}_R2.fastq.gz
```

If the sample was split across several lanes (`L001`, `L002`, ...), decompress every file, then join each read direction in lane order. Concatenated gzip files are valid gzip:

```bash
cd ${GENOME_DIR}/${SAMPLE}/fastq
cat ${SAMPLE}_S1_L00*_R1_001.fastq.gz > ${SAMPLE}_R1.fastq.gz
cat ${SAMPLE}_S1_L00*_R2_001.fastq.gz > ${SAMPLE}_R2.fastq.gz
```

## Notes
- ORA reference files must match the sequencing run (provided alongside ORA files)
- If you receive FASTQ.gz directly (e.g., from a resequencing service), skip this step
