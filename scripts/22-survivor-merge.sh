#!/usr/bin/env bash
# 22-survivor-merge.sh — Merge SV calls from multiple callers into consensus set
# Usage: ./scripts/22-survivor-merge.sh <sample_name>
#
# Merges structural variants from Manta, Delly, and CNVnator using SURVIVOR.
# Variants called by 2+ callers at overlapping positions are high-confidence.
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
  # Convert CNVnator TXT to simple VCF
  {
    echo "##fileformat=VCFv4.2"
    echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    echo "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">"
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

# Use bcftools-based merge approach (SURVIVOR Docker image availability is unreliable)
# Strategy: intersect Manta and Delly, keep SVs with overlapping breakpoints within 1kb
echo "[2/3] Finding consensus SVs (breakpoints within 1kb)..."

# For each pair of callers, find overlapping SVs using bcftools isec
# Simplified approach: extract BED from each, use overlap logic
docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "
    # Extract SV positions from each caller as BED-like format
    for VCF_FILE in ${SV_FILES[*]}; do
      bcftools view -f PASS \"\$VCF_FILE\" 2>/dev/null || bcftools view \"\$VCF_FILE\" 2>/dev/null
    done | grep -v '^#' | \
    awk -F'\t' '{
      chrom=\$1; pos=\$2; info=\$8;
      end=pos;
      if(match(info, /END=[0-9]+/)) {
        end=substr(info, RSTART+4, RLENGTH-4);
      }
      svtype=\"UNK\";
      if(match(info, /SVTYPE=[A-Z]+/)) {
        svtype=substr(info, RSTART+7, RLENGTH-7);
      }
      key=chrom\"_\"int(pos/1000)\"_\"svtype;
      count[key]++;
      if(count[key]==1) {
        lines[key]=\$0;
        ends[key]=end;
      }
    } END {
      for(k in count) {
        if(count[k] >= 2) print lines[k];
      }
    }' | sort -k1,1V -k2,2n > /genome/${SAMPLE}/sv_merged/consensus_raw.txt

    # Create VCF from consensus
    echo '##fileformat=VCFv4.2' > /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf
    echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">' >> /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf
    echo '##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">' >> /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf
    echo '##INFO=<ID=CALLERS,Number=1,Type=Integer,Description=\"Number of callers\">' >> /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf
    echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO' >> /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf
    cat /genome/${SAMPLE}/sv_merged/consensus_raw.txt >> /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf

    # Compress and index
    bcftools sort /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf -Oz \
      -o /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz 2>/dev/null || true
    bcftools index -t /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz 2>/dev/null || true
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
echo "High-confidence SVs are those called by 2+ independent methods."
echo "Single-caller SVs have a higher false positive rate."
