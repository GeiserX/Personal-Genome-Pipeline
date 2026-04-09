/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HTML_REPORT — Summary HTML report of key pipeline results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generates a self-contained HTML file summarising ClinVar hits, pharmacogenomics,
    cancer predisposition (CPSR), clinical filtering, and slivar variant prioritization.
    Each input section is optional — when a path is empty the section is omitted.
    NOTE: BAM-analysis outputs (HLA, ExpansionHunter, telomere, Cyrius, coverage) and
    SV outputs are not included — see their individual output directories for results.

    Equivalent to: scripts/24-html-report.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HTML_REPORT {
    tag "$meta.id"
    label 'process_low'

    container 'staphb/bcftools:1.21'

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(clinvar_dir), path(pharmcat_html), path(clinical_vcf), path(cpsr_html), path(slivar_vcf)

    output:
    tuple val(meta), path("${meta.id}_report.html"), emit: html_report
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_clinvar  = clinvar_dir   && !clinvar_dir.name.startsWith('EMPTY')   ? true : false
    def has_pharmcat = pharmcat_html && !pharmcat_html.name.startsWith('EMPTY') ? true : false
    def has_clinical = clinical_vcf  && !clinical_vcf.name.startsWith('EMPTY')  ? true : false
    def has_cpsr     = cpsr_html     && !cpsr_html.name.startsWith('EMPTY')     ? true : false
    def has_slivar   = slivar_vcf    && !slivar_vcf.name.startsWith('EMPTY')    ? true : false
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # --- Collect ClinVar hits ---
    CLINVAR_HITS="N/A"
    CLINVAR_ROWS=""
    if [ "${has_clinvar}" = "true" ] && [ -d "${clinvar_dir}" ]; then
        ISEC_FILE="${clinvar_dir}/0002.vcf"
        if [ -f "\$ISEC_FILE" ]; then
            CLINVAR_HITS=\$(grep -c -v '^#' "\$ISEC_FILE" 2>/dev/null || echo "0")
            CLINVAR_ROWS=\$(grep -v '^#' "\$ISEC_FILE" 2>/dev/null | head -20 | \\
                awk -F'\\t' '{
                    gene="."; clnsig=".";
                    if(match(\$8,/GENEINFO=[^;]+/)) gene=substr(\$8,RSTART+9,RLENGTH-9);
                    if(match(\$8,/CLNSIG=[^;]+/)) clnsig=substr(\$8,RSTART+7,RLENGTH-7);
                    # Escape all fields rendered into HTML
                    gsub(/&/,"\\&amp;",\$1); gsub(/</,"\\&lt;",\$1); gsub(/>/,"\\&gt;",\$1);
                    gsub(/&/,"\\&amp;",\$4); gsub(/</,"\\&lt;",\$4); gsub(/>/,"\\&gt;",\$4);
                    gsub(/&/,"\\&amp;",\$5); gsub(/</,"\\&lt;",\$5); gsub(/>/,"\\&gt;",\$5);
                    gsub(/&/,"\\&amp;",gene); gsub(/</,"\\&lt;",gene); gsub(/>/,"\\&gt;",gene); gsub(/"/,"\\&quot;",gene);
                    gsub(/&/,"\\&amp;",clnsig); gsub(/</,"\\&lt;",clnsig); gsub(/>/,"\\&gt;",clnsig); gsub(/"/,"\\&quot;",clnsig);
                    printf "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\\n",\$1,\$2,\$4,\$5,gene"|"clnsig;
                }' || true)
        fi
    fi

    # --- PharmCAT status ---
    PHARMCAT_STATUS="Not run"
    if [ "${has_pharmcat}" = "true" ] && [ -f "${pharmcat_html}" ]; then
        PHARMCAT_STATUS="Complete"
    fi

    # --- Clinical filter counts ---
    CLINICAL_TOTAL="N/A"
    if [ "${has_clinical}" = "true" ] && [ -f "${clinical_vcf}" ]; then
        CLINICAL_TOTAL=\$(bcftools view -H "${clinical_vcf}" 2>/dev/null | wc -l || echo "N/A")
    fi

    # --- CPSR status ---
    CPSR_STATUS="Not run"
    if [ "${has_cpsr}" = "true" ] && [ -f "${cpsr_html}" ]; then
        CPSR_STATUS="Complete"
    fi

    # --- Slivar status ---
    SLIVAR_STATUS="Not run"
    SLIVAR_COUNT="N/A"
    if [ "${has_slivar}" = "true" ] && [ -f "${slivar_vcf}" ]; then
        SLIVAR_STATUS="Complete"
        SLIVAR_COUNT=\$(bcftools view -H "${slivar_vcf}" 2>/dev/null | wc -l || echo "N/A")
    fi

    # --- ClinVar badge colour ---
    CLINVAR_BADGE="badge-green"
    if [ "\$CLINVAR_HITS" != "N/A" ] && [ "\$CLINVAR_HITS" -gt 5 ] 2>/dev/null; then
        CLINVAR_BADGE="badge-yellow"
    fi

    # --- Generate HTML ---
    cat > ${meta.id}_report.html << 'HTMLHEAD'
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

    cat >> ${meta.id}_report.html << EOF
<div class="header">
  <h1>Genomic Analysis Report</h1>
  <div class="meta">
    Sample: <strong>${meta.id}</strong> &nbsp;|&nbsp;
    Generated: \$(date '+%Y-%m-%d %H:%M') &nbsp;|&nbsp;
    Pipeline: <a href="https://github.com/GeiserX/Personal-Genome-Pipeline" style="color:#aed6f1">Personal-Genome-Pipeline</a>
  </div>
</div>
<div class="grid">
EOF

    # Card: ClinVar
    cat >> ${meta.id}_report.html << EOF
  <div class="card">
    <h2>ClinVar Screening</h2>
    <div class="stat"><span class="label">ClinVar matches</span>
      <span class="value"><span class="badge \${CLINVAR_BADGE}">\${CLINVAR_HITS}</span></span></div>
    <div class="stat"><span class="label">Status</span>
      <span class="value">\$([ "\$CLINVAR_HITS" = "N/A" ] && echo '<span class="badge badge-gray">Not run</span>' || echo '<span class="badge badge-green">Complete</span>')</span></div>
  </div>
EOF

    # Card: Pharmacogenomics
    cat >> ${meta.id}_report.html << EOF
  <div class="card">
    <h2>Pharmacogenomics</h2>
    <div class="stat"><span class="label">PharmCAT report</span>
      <span class="value"><span class="badge \$([ "\$PHARMCAT_STATUS" = "Complete" ] && echo "badge-green" || echo "badge-gray")">\${PHARMCAT_STATUS}</span></span></div>
    <div class="stat"><span class="label">Tip</span>
      <span class="value" style="font-weight:normal;font-size:13px">Open the PharmCAT HTML report for full drug-gene details</span></div>
  </div>
EOF

    # Card: Cancer Predisposition
    cat >> ${meta.id}_report.html << EOF
  <div class="card">
    <h2>Cancer Predisposition</h2>
    <div class="stat"><span class="label">CPSR report</span>
      <span class="value"><span class="badge \$([ "\$CPSR_STATUS" = "Complete" ] && echo "badge-green" || echo "badge-gray")">\${CPSR_STATUS}</span></span></div>
    <div class="stat"><span class="label">Tip</span>
      <span class="value" style="font-weight:normal;font-size:13px">Open the CPSR HTML report for tier classification details</span></div>
  </div>
EOF

    # Card: Clinical Filter
    cat >> ${meta.id}_report.html << EOF
  <div class="card">
    <h2>Clinical Variant Filter</h2>
    <div class="stat"><span class="label">Total interesting variants</span><span class="value">\${CLINICAL_TOTAL}</span></div>
  </div>
EOF

    # Card: Slivar
    cat >> ${meta.id}_report.html << EOF
  <div class="card">
    <h2>Variant Prioritization (Slivar)</h2>
    <div class="stat"><span class="label">Status</span>
      <span class="value"><span class="badge \$([ "\$SLIVAR_STATUS" = "Complete" ] && echo "badge-green" || echo "badge-gray")">\${SLIVAR_STATUS}</span></span></div>
    <div class="stat"><span class="label">Prioritized variants</span><span class="value">\${SLIVAR_COUNT}</span></div>
    <div class="stat"><span class="label">Tip</span>
      <span class="value" style="font-weight:normal;font-size:13px">Rare HIGH/MODERATE + deleterious + ClinVar pathogenic tiers</span></div>
  </div>
EOF

    # ClinVar detail table
    if [ -n "\$CLINVAR_ROWS" ]; then
        cat >> ${meta.id}_report.html << EOF
  <div class="card full-width">
    <h2>ClinVar Hits (Top 20)</h2>
    <table>
      <tr><th>Chr</th><th>Position</th><th>Ref</th><th>Alt</th><th>Gene | Significance</th></tr>
      \${CLINVAR_ROWS}
    </table>
  </div>
EOF
    fi

    # Close grid, disclaimer, footer
    cat >> ${meta.id}_report.html << 'HTMLFOOT'
</div>

<div class="disclaimer">
  <strong>Disclaimer:</strong> This report is for educational and research purposes only.
  It is not a clinical diagnosis. Always discuss genomic findings with a qualified healthcare
  professional before making any medical decisions. Variants of Uncertain Significance (VUS)
  are not clinically actionable.
</div>

<div class="footer">
  Generated by <a href="https://github.com/GeiserX/Personal-Genome-Pipeline">Personal-Genome-Pipeline</a>
  — 100% local analysis, no data uploaded
</div>

</div>
</body>
</html>
HTMLFOOT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
