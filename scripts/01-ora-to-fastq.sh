#!/usr/bin/env bash
# ORA to FASTQ — Convert Illumina ORA-compressed files to standard FASTQ
# Input: .ora file(s) from Illumina DRAGEN sequencing
# Output: .fastq.gz in $GENOME_DIR/<sample>/fastq/
# NOTE: orad is a native binary (not Docker), typically at /opt/orad/bin/orad
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <ora_reference_dir> <ora_file>}
ORA_REF=${2:?Usage: $0 <sample_name> <ora_reference_dir> <ora_file>}
ORA_FILE=${3:?Usage: $0 <sample_name> <ora_reference_dir> <ora_file>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
ORAD=${ORAD:-/opt/orad/bin/orad}
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/fastq"

echo "=== ORA to FASTQ: ${SAMPLE} ==="
echo "ORA file: ${ORA_FILE}"
echo "ORA reference: ${ORA_REF}"
echo "Output: ${OUTPUT_DIR}/"

# Validate inputs
if [ ! -f "$ORA_FILE" ]; then
  echo "ERROR: ORA file not found: ${ORA_FILE}" >&2
  exit 1
fi

if [ ! -d "$ORA_REF" ]; then
  echo "ERROR: ORA reference directory not found: ${ORA_REF}" >&2
  exit 1
fi

if ! command -v "$ORAD" &>/dev/null; then
  echo "ERROR: orad binary not found at ${ORAD}" >&2
  echo "Download from Illumina and set ORAD=/path/to/orad" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

"$ORAD" \
  --ora-reference "$ORA_REF" \
  --output-directory "$OUTPUT_DIR" \
  "$ORA_FILE"

echo "=== ORA to FASTQ complete ==="
echo "FASTQ files: ${OUTPUT_DIR}/"
ls -lh "$OUTPUT_DIR"/*.fastq.gz 2>/dev/null || echo "WARNING: No .fastq.gz files found in output"
