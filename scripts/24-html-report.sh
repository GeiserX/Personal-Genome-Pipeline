#!/usr/bin/env bash
# 24-html-report.sh — Generate a self-contained HTML report of all pipeline results
# Usage: ./scripts/24-html-report.sh <sample_name>
#
# Scans all output directories and produces a single HTML file that can be
# opened in any browser. No internet connection needed to view it.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTPUT="${SAMPLE_DIR}/${SAMPLE}_report.html"

echo "============================================"
echo "  Step 24: HTML Report Generator"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${OUTPUT}"
echo "============================================"
echo ""

# Collect data from all pipeline outputs
echo "[1/2] Scanning pipeline outputs..."

# --- Variant counts ---
VCF_TOTAL="N/A"
VCF_PASS="N/A"
VCF_SNPS="N/A"
VCF_INDELS="N/A"
if [ -f "${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz" ]; then
  STATS=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools stats "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null | grep "^SN" || true)
  VCF_TOTAL=$(echo "$STATS" | grep "number of records:" | awk '{print $NF}' || echo "N/A")
  VCF_SNPS=$(echo "$STATS" | grep "number of SNPs:" | awk '{print $NF}' || echo "N/A")
  VCF_INDELS=$(echo "$STATS" | grep "number of indels:" | awk '{print $NF}' || echo "N/A")
  VCF_PASS=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -f PASS -H "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
fi

# --- ClinVar ---
CLINVAR_HITS="N/A"
CLINVAR_DETAILS=""
if [ -d "${SAMPLE_DIR}/clinvar/isec" ]; then
  ISEC_FILE="${SAMPLE_DIR}/clinvar/isec/0002.vcf"
  if [ -f "$ISEC_FILE" ]; then
    CLINVAR_HITS=$(grep -c -v "^#" "$ISEC_FILE" 2>/dev/null || echo "0")
    CLINVAR_DETAILS=$(grep -v "^#" "$ISEC_FILE" 2>/dev/null | head -20 | \
      awk -F'\t' '{
        gene="."; clnsig=".";
        if(match($8,/GENEINFO=[^;]+/)) gene=substr($8,RSTART+9,RLENGTH-9);
        if(match($8,/CLNSIG=[^;]+/)) clnsig=substr($8,RSTART+7,RLENGTH-7);
        printf "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n",$1,$2,$4,$5,gene"|"clnsig;
      }' || true)
  fi
fi

# --- PharmCAT ---
# Step 7 writes to vcf/ (PharmCAT outputs alongside the VCF)
PHARMCAT_STATUS="Not run"
if find "${SAMPLE_DIR}/vcf" -maxdepth 1 -name "*.report.html" 2>/dev/null | grep -q .; then
  PHARMCAT_STATUS="Complete"
elif find "${SAMPLE_DIR}/pharmcat" -maxdepth 1 -name "*.report.html" 2>/dev/null | grep -q .; then
  PHARMCAT_STATUS="Complete"
fi

# --- ExpansionHunter ---
EH_STATUS="Not run"
EH_DETAILS=""
EH_FILE=$(find "${SAMPLE_DIR}/expansion_hunter" -maxdepth 1 -name "*_eh.vcf" 2>/dev/null | head -1)
if [ -n "$EH_FILE" ] && [ -f "$EH_FILE" ]; then
  EH_STATUS="Complete"
  # Check key loci — REPCN (repeat copy number) is FORMAT field 3 (GT:SO:REPCN:...)
  for LOCUS in HTT FMR1 C9ORF72 ATXN1 DMPK; do
    REPEAT=$(grep -w "$LOCUS" "$EH_FILE" 2>/dev/null | head -1 | \
      awk -F'\t' '{split($10,a,":"); print a[3]}' || echo "N/A")
    EH_DETAILS="${EH_DETAILS}<tr><td>${LOCUS}</td><td>${REPEAT:-N/A}</td></tr>"
  done
fi

# --- Manta SVs ---
MANTA_TOTAL="N/A"
MANTA_PASS="N/A"
if [ -f "${SAMPLE_DIR}/manta/results/variants/diploidSV.vcf.gz" ]; then
  MANTA_TOTAL=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -H "/genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
  MANTA_PASS=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -f PASS -H "/genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
fi

# --- CPSR ---
CPSR_STATUS="Not run"
if find "${SAMPLE_DIR}/cpsr" -maxdepth 1 -name "*.cpsr.grch38.html" 2>/dev/null | grep -q .; then
  CPSR_STATUS="Complete"
fi

# --- ROH ---
ROH_TOTAL_MB="N/A"
ROH_MAX_MB="N/A"
if [ -f "${SAMPLE_DIR}/vcf/${SAMPLE}_roh.txt" ]; then
  ROH_TOTAL_MB=$(awk '$1=="RG" {sum+=$6} END {printf "%.1f", sum/1000000}' \
    "${SAMPLE_DIR}/vcf/${SAMPLE}_roh.txt" 2>/dev/null || echo "N/A")
  ROH_MAX_MB=$(awk '$1=="RG" {if($6>max) max=$6} END {printf "%.1f", max/1000000}' \
    "${SAMPLE_DIR}/vcf/${SAMPLE}_roh.txt" 2>/dev/null || echo "N/A")
fi

# --- Haplogroup ---
# Step 12 writes to mito/ (not haplogrep/)
HAPLOGROUP="N/A"
if [ -f "${SAMPLE_DIR}/mito/${SAMPLE}_haplogroup.txt" ]; then
  HAPLOGROUP=$(awk -F'\t' 'NR==2 {gsub(/"/, "", $2); print $2}' "${SAMPLE_DIR}/mito/${SAMPLE}_haplogroup.txt" 2>/dev/null || echo "N/A")
fi

# --- Telomere ---
# Step 10 writes to telomere/${SAMPLE}/${SAMPLE}/${SAMPLE}_summary.tsv
# tel_content (GC-corrected telomeric reads per million) is column 11
TELOMERE="N/A"
TEL_FILE="${SAMPLE_DIR}/telomere/${SAMPLE}/${SAMPLE}/${SAMPLE}_summary.tsv"
if [ -f "$TEL_FILE" ]; then
  TELOMERE=$(awk -F'\t' 'NR==2 {print $11}' "$TEL_FILE" 2>/dev/null || echo "N/A")
fi

# --- Clinical filter ---
CLINICAL_TOTAL="N/A"
CLINICAL_HIGH="N/A"
if [ -f "${SAMPLE_DIR}/clinical/${SAMPLE}_clinical.vcf.gz" ]; then
  CLINICAL_TOTAL=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
fi
if [ -f "${SAMPLE_DIR}/clinical/${SAMPLE}_high_impact.vcf.gz" ]; then
  CLINICAL_HIGH=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
fi

# --- Mitochondrial ---
# Step 20 writes to mito/ (not mtoolbox/). Heteroplasmy is detected via AF field, not GT.
MITO_PASS="N/A"
MITO_HETERO="N/A"
MITO_FILE="${SAMPLE_DIR}/mito/${SAMPLE}_chrM_filtered.vcf.gz"
if [ -f "$MITO_FILE" ]; then
  MITO_PASS=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -f PASS -H "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
  # Heteroplasmic = AF < 0.95 (not homoplasmic)
  MITO_HETERO=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools query -f '[%AF]\n' -i 'FILTER="PASS"' "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz" 2>/dev/null | \
    awk '{if($1+0 < 0.95) c++} END {print c+0}' || echo "N/A")
fi

# --- CNVnator ---
CNV_COUNT="N/A"
if [ -f "${SAMPLE_DIR}/cnvnator/${SAMPLE}_cnvs.txt" ]; then
  CNV_COUNT=$(wc -l < "${SAMPLE_DIR}/cnvnator/${SAMPLE}_cnvs.txt" 2>/dev/null || echo "N/A")
fi

# --- Delly ---
DELLY_PASS="N/A"
DELLY_FILE="${SAMPLE_DIR}/delly/${SAMPLE}_sv.vcf.gz"
if [ -f "$DELLY_FILE" ]; then
  DELLY_PASS=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -f PASS -H "/genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
fi

echo "[2/2] Generating HTML report..."

# Generate HTML
cat > "$OUTPUT" << 'HTMLHEAD'
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Personal Genome Pipeline Report</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         background: #f5f5f5; color: #333; line-height: 1.6; }
  .container { max-width: 1100px; margin: 0 auto; padding: 20px; }
  .header { background: linear-gradient(135deg, #1a5276 0%, #2e86c1 100%);
            color: white; padding: 30px; border-radius: 12px; margin-bottom: 24px; }
  .header h1 { font-size: 28px; margin-bottom: 8px; }
  .header .meta { opacity: 0.85; font-size: 14px; }
  .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 20px; margin-bottom: 24px; }
  .card { background: white; border-radius: 10px; padding: 24px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
  .card h2 { font-size: 18px; color: #1a5276; margin-bottom: 16px;
             padding-bottom: 8px; border-bottom: 2px solid #eee; }
  .stat { display: flex; justify-content: space-between; padding: 8px 0;
          border-bottom: 1px solid #f0f0f0; }
  .stat:last-child { border-bottom: none; }
  .stat .label { color: #666; }
  .stat .value { font-weight: 600; }
  .badge { display: inline-block; padding: 2px 10px; border-radius: 12px;
           font-size: 13px; font-weight: 600; }
  .badge-green { background: #d5f5e3; color: #196f3d; }
  .badge-yellow { background: #fef9e7; color: #7d6608; }
  .badge-red { background: #fadbd8; color: #922b21; }
  .badge-blue { background: #d6eaf8; color: #1a5276; }
  .badge-gray { background: #eee; color: #666; }
  table { width: 100%; border-collapse: collapse; font-size: 14px; }
  th, td { padding: 8px 12px; text-align: left; border-bottom: 1px solid #eee; }
  th { background: #f8f9fa; font-weight: 600; color: #555; }
  .full-width { grid-column: 1 / -1; }
  .disclaimer { background: #fff3cd; border: 1px solid #ffc107; border-radius: 8px;
                padding: 16px; margin-top: 24px; font-size: 14px; }
  .footer { text-align: center; color: #999; font-size: 13px; margin-top: 24px; padding: 16px; }
  @media (max-width: 700px) { .grid { grid-template-columns: 1fr; } }
</style>
</head>
<body>
<div class="container">
HTMLHEAD

# Header
cat >> "$OUTPUT" << EOF
<div class="header">
  <h1>Genomic Analysis Report</h1>
  <div class="meta">
    Sample: <strong>${SAMPLE}</strong> &nbsp;|&nbsp;
    Generated: $(date '+%Y-%m-%d %H:%M') &nbsp;|&nbsp;
    Pipeline: <a href="https://github.com/GeiserX/personal-genome-pipeline" style="color:#aed6f1">personal-genome-pipeline</a>
  </div>
</div>

<div class="grid">
EOF

# Card: Variant Summary
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>Variant Calling</h2>
    <div class="stat"><span class="label">Total variants</span><span class="value">${VCF_TOTAL}</span></div>
    <div class="stat"><span class="label">PASS variants</span><span class="value">${VCF_PASS}</span></div>
    <div class="stat"><span class="label">SNPs</span><span class="value">${VCF_SNPS}</span></div>
    <div class="stat"><span class="label">Indels</span><span class="value">${VCF_INDELS}</span></div>
  </div>
EOF

# Card: ClinVar
CLINVAR_BADGE="badge-green"
if [ "$CLINVAR_HITS" != "N/A" ] && [ "$CLINVAR_HITS" -gt 5 ] 2>/dev/null; then
  CLINVAR_BADGE="badge-yellow"
fi
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>ClinVar Screening</h2>
    <div class="stat"><span class="label">ClinVar matches</span>
      <span class="value"><span class="badge ${CLINVAR_BADGE}">${CLINVAR_HITS}</span></span></div>
    <div class="stat"><span class="label">Status</span>
      <span class="value">$([ "$CLINVAR_HITS" = "N/A" ] && echo '<span class="badge badge-gray">Not run</span>' || echo '<span class="badge badge-green">Complete</span>')</span></div>
  </div>
EOF

# Card: Pharmacogenomics
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>Pharmacogenomics</h2>
    <div class="stat"><span class="label">PharmCAT report</span>
      <span class="value"><span class="badge $([ "$PHARMCAT_STATUS" = "Complete" ] && echo "badge-green" || echo "badge-gray")">${PHARMCAT_STATUS}</span></span></div>
    <div class="stat"><span class="label">Tip</span>
      <span class="value" style="font-weight:normal;font-size:13px">Open the HTML report in your browser for full drug-gene details</span></div>
  </div>
EOF

# Card: Structural Variants
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>Structural Variants</h2>
    <div class="stat"><span class="label">Manta SVs (total)</span><span class="value">${MANTA_TOTAL}</span></div>
    <div class="stat"><span class="label">Manta SVs (PASS)</span><span class="value">${MANTA_PASS}</span></div>
    <div class="stat"><span class="label">Delly SVs (PASS)</span><span class="value">${DELLY_PASS}</span></div>
    <div class="stat"><span class="label">CNVnator CNVs</span><span class="value">${CNV_COUNT}</span></div>
  </div>
EOF

# Card: Cancer Predisposition
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>Cancer Predisposition</h2>
    <div class="stat"><span class="label">CPSR report</span>
      <span class="value"><span class="badge $([ "$CPSR_STATUS" = "Complete" ] && echo "badge-green" || echo "badge-gray")">${CPSR_STATUS}</span></span></div>
    <div class="stat"><span class="label">Tip</span>
      <span class="value" style="font-weight:normal;font-size:13px">Open the CPSR HTML report for tier classification details</span></div>
  </div>
EOF

# Card: Repeat Expansions
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>Repeat Expansions</h2>
    <div class="stat"><span class="label">ExpansionHunter</span>
      <span class="value"><span class="badge $([ "$EH_STATUS" = "Complete" ] && echo "badge-green" || echo "badge-gray")">${EH_STATUS}</span></span></div>
EOF
if [ -n "$EH_DETAILS" ]; then
  cat >> "$OUTPUT" << EOF
    <table><tr><th>Locus</th><th>Repeat Count</th></tr>${EH_DETAILS}</table>
EOF
fi
echo "  </div>" >> "$OUTPUT"

# Card: Ancestry & Identity
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>Ancestry &amp; Identity</h2>
    <div class="stat"><span class="label">Mitochondrial haplogroup</span><span class="value">${HAPLOGROUP}</span></div>
    <div class="stat"><span class="label">ROH total</span><span class="value">${ROH_TOTAL_MB} MB</span></div>
    <div class="stat"><span class="label">ROH largest segment</span><span class="value">${ROH_MAX_MB} MB</span></div>
    <div class="stat"><span class="label">Telomere content</span><span class="value">${TELOMERE}</span></div>
  </div>
EOF

# Card: Mitochondrial
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>Mitochondrial Analysis</h2>
    <div class="stat"><span class="label">chrM variants (PASS)</span><span class="value">${MITO_PASS}</span></div>
    <div class="stat"><span class="label">Heteroplasmic</span><span class="value">${MITO_HETERO}</span></div>
  </div>
EOF

# Card: Clinical Filter
cat >> "$OUTPUT" << EOF
  <div class="card">
    <h2>Clinical Variant Filter</h2>
    <div class="stat"><span class="label">Total interesting variants</span><span class="value">${CLINICAL_TOTAL}</span></div>
    <div class="stat"><span class="label">HIGH impact (LoF)</span><span class="value">${CLINICAL_HIGH}</span></div>
  </div>
EOF

# ClinVar detail table (if hits exist)
if [ -n "$CLINVAR_DETAILS" ]; then
  cat >> "$OUTPUT" << EOF
  <div class="card full-width">
    <h2>ClinVar Hits (Top 20)</h2>
    <table>
      <tr><th>Chr</th><th>Position</th><th>Ref</th><th>Alt</th><th>Gene | Significance</th></tr>
      ${CLINVAR_DETAILS}
    </table>
  </div>
EOF
fi

# Close grid, add disclaimer and footer
cat >> "$OUTPUT" << 'HTMLFOOT'
</div>

<div class="disclaimer">
  <strong>Disclaimer:</strong> This report is for educational and research purposes only.
  It is not a clinical diagnosis. Always discuss genomic findings with a qualified healthcare
  professional before making any medical decisions. Variants of Uncertain Significance (VUS)
  are not clinically actionable.
</div>

<div class="footer">
  Generated by <a href="https://github.com/GeiserX/personal-genome-pipeline">personal-genome-pipeline</a>
  — 100% local analysis, no data uploaded
</div>

</div>
</body>
</html>
HTMLFOOT

REPORT_SIZE=$(wc -c < "$OUTPUT" 2>/dev/null || echo 0)
REPORT_KB=$(( REPORT_SIZE / 1024 ))

echo ""
echo "============================================"
echo "  HTML report generated: ${OUTPUT}"
echo "  Size: ${REPORT_KB} KB"
echo "============================================"
echo ""
echo "Open in your browser:"
echo "  open ${OUTPUT}          # macOS"
echo "  xdg-open ${OUTPUT}     # Linux"
echo "  start ${OUTPUT}         # Windows (WSL)"
