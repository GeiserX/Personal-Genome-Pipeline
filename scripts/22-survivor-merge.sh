#!/usr/bin/env bash
# 22-survivor-merge.sh — [EXPERIMENTAL] Merge SV calls from multiple callers
# Usage: ./scripts/22-survivor-merge.sh <sample_name>
#
# Finds structural variants called by 2+ callers at overlapping positions.
# Uses a simplified bcftools-based approach (no SURVIVOR Docker image).
#
# EXPERIMENTAL: This uses a heuristic position-binning approach that may
# over-count calls from the same caller. Results should be treated as a
# rough intersection, not a true consensus merge. For production use,
# consider SURVIVOR or Jasmine with proper multi-sample VCF merging.
#
# Requires: Output from steps 4 (Manta), 19 (Delly), and/or 18 (CNVnator)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

OUTDIR="${GENOME_DIR}/${SAMPLE}/sv_merged"
mkdir -p "$OUTDIR"

echo "============================================"
echo "  Step 22: SV Consensus Merge (SURVIVOR)"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# Collect available SV VCFs
SV_FILES=()
CALLERS=""

# Manta
MANTA_VCF="${GENOME_DIR}/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz"
if [ -f "$MANTA_VCF" ]; then
  SV_FILES+=("/genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz")
  CALLERS="${CALLERS}Manta "
  echo "  [OK] Manta SVs found"
else
  echo "  [--] Manta SVs not found (run step 4 first)"
fi

# Delly
DELLY_VCF="${GENOME_DIR}/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz"
if [ -f "$DELLY_VCF" ]; then
  SV_FILES+=("/genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz")
  CALLERS="${CALLERS}Delly "
  echo "  [OK] Delly SVs found"
else
  echo "  [--] Delly SVs not found (run step 19 first)"
fi

# CNVnator (convert to VCF format if only TXT exists)
CNVNATOR_TXT="${GENOME_DIR}/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.txt"
CNVNATOR_VCF="${GENOME_DIR}/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.vcf.gz"
if [ -f "$CNVNATOR_VCF" ]; then
  SV_FILES+=("/genome/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.vcf.gz")
  CALLERS="${CALLERS}CNVnator "
  echo "  [OK] CNVnator CNVs found (VCF)"
elif [ -f "$CNVNATOR_TXT" ]; then
  echo "  [OK] CNVnator CNVs found (TXT — converting to VCF)..."
  # Convert CNVnator TXT to simple VCF (with contig headers for bcftools)
  {
    echo "##fileformat=VCFv4.2"
    echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    echo "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">"
    # Add contig headers from reference .fai (required by bcftools sort)
    REF_FAI="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta.fai"
    if [ -f "$REF_FAI" ]; then
      awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "$REF_FAI"
    fi
    echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
    awk '{
      split($2,a,":");
      split(a[2],b,"-");
      svtype="DEL";
      alt="<DEL>";
      svlen=-(b[2]-b[1]);
      if($1=="duplication") { svtype="DUP"; alt="<DUP>"; svlen=b[2]-b[1]; }
      printf "%s\t%s\t.\tN\t%s\t.\tPASS\tSVTYPE=%s;END=%s;SVLEN=%d\n",
        a[1], b[1], alt, svtype, b[2], svlen;
    }' "$CNVNATOR_TXT"
  } > "${GENOME_DIR}/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.vcf"

  docker run --rm --user root \
    -v "${GENOME_DIR}:/genome" \
    staphb/bcftools:1.21 \
    bash -c "bcftools sort /genome/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.vcf -Oz \
      -o /genome/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.vcf.gz"

  SV_FILES+=("/genome/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.vcf.gz")
  CALLERS="${CALLERS}CNVnator "
else
  echo "  [--] CNVnator CNVs not found (run step 18 first)"
fi

echo ""

if [ ${#SV_FILES[@]} -lt 2 ]; then
  echo "ERROR: Need at least 2 SV callers for consensus merge."
  echo "  Found: ${CALLERS:-none}"
  echo "  Run steps 4 (Manta) and 19 (Delly) first."
  exit 1
fi

echo "Merging SVs from: ${CALLERS}"
echo "  Callers found: ${#SV_FILES[@]}"
echo ""

# Create file list for SURVIVOR
echo "[1/3] Preparing SV file list..."
FILE_LIST="/genome/${SAMPLE}/sv_merged/sv_files.txt"
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "
    > ${FILE_LIST}
    for f in ${SV_FILES[*]}; do
      echo \$f >> ${FILE_LIST}
    done
    cat ${FILE_LIST}
  "

# Strategy: Extract BED-like positions from each caller separately, then find
# positions seen in 2+ callers using per-caller tagging.
echo "[2/3] Finding consensus SVs (breakpoints within 1kb, 2+ callers)..."

# Step A: Extract SV positions per caller as "caller\tchr\tbin\tsvtype" for counting
docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "
    CALLER_IDX=0
    for VCF_FILE in ${SV_FILES[*]}; do
      CALLER_IDX=\$((CALLER_IDX + 1))
      (bcftools view -f PASS \"\$VCF_FILE\" 2>/dev/null || bcftools view \"\$VCF_FILE\" 2>/dev/null) | \
        grep -v '^#' | \
        awk -F'\t' -v caller=\$CALLER_IDX '{
          chrom=\$1; pos=\$2; info=\$8;
          end=pos;
          if(match(info, /END=[0-9]+/)) end=substr(info, RSTART+4, RLENGTH-4);
          svtype=\"UNK\";
          if(match(info, /SVTYPE=[A-Z]+/)) svtype=substr(info, RSTART+7, RLENGTH-7);
          bin=int(pos/1000);
          # Output: chrom, pos, end, svtype, caller_id, 8-col VCF line
          printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\tN\t<%s>\t.\tPASS\tSVTYPE=%s;END=%s\n\",
            chrom, bin, svtype, caller, pos, chrom, pos, svtype, svtype, end;
        }'
    done > /genome/${SAMPLE}/sv_merged/all_sv_tagged.tsv

    # Step B: Find bins seen by 2+ distinct callers, emit one clean 8-col VCF record each
    # Uses pipe-delimited caller string (mawk-compatible — no nested arrays)
    awk -F'\t' '{
      key=\$1\"_\"\$2\"_\"\$3;
      caller=\$4;
      if(!(key in seen)) {
        seen[key]=caller;
        line[key]=\$6\"\t\"\$7\"\t.\tN\t\"\$10\"\t.\tPASS\t\"\$13;
      } else if(index(seen[key], caller) == 0) {
        seen[key]=seen[key]\"|\"caller;
      }
    } END {
      for(k in seen) {
        n=split(seen[k], a, \"|\");
        if(n >= 2) print line[k];
      }
    }' /genome/${SAMPLE}/sv_merged/all_sv_tagged.tsv | \
      sort -k1,1V -k2,2n > /genome/${SAMPLE}/sv_merged/consensus_raw.txt

    # Step C: Build a valid VCF with contig headers
    {
      echo '##fileformat=VCFv4.2'
      echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">'
      echo '##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">'
      if [ -f /genome/reference/Homo_sapiens_assembly38.fasta.fai ]; then
        awk '{printf \"##contig=<ID=%s,length=%s>\\n\", \$1, \$2}' /genome/reference/Homo_sapiens_assembly38.fasta.fai
      fi
      printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
      cat /genome/${SAMPLE}/sv_merged/consensus_raw.txt
    } > /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf

    # Compress and index
    bcftools sort /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf -Oz \
      -o /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz
    bcftools index -t /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz
  "

echo "[3/3] Counting results..."

CONSENSUS_COUNT=0
if [ -f "${OUTDIR}/${SAMPLE}_sv_consensus.vcf.gz" ]; then
  CONSENSUS_COUNT=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -H "/genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz" 2>/dev/null | wc -l || echo 0)
fi

# Also count per-caller for comparison
echo ""
echo "============================================"
echo "  SV Consensus Merge complete: ${SAMPLE}"
echo "  Callers: ${CALLERS}"
echo "  Consensus SVs (2+ callers): ${CONSENSUS_COUNT}"
echo ""
echo "  Output: ${OUTDIR}/${SAMPLE}_sv_consensus.vcf.gz"
echo "============================================"
echo ""
echo "Multi-caller SVs (2+ callers) have lower false-positive rates than single-caller calls."
echo "Note: 1kb position binning is approximate — see docs for limitations."
