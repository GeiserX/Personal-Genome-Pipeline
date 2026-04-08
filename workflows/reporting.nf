/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORTING — HTML Report Generation & MultiQC Aggregation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Produces two kinds of report:
    1. Per-sample consolidated HTML report (variant summary, ClinVar, PGx, CPSR, slivar)
    2. Cross-sample MultiQC dashboard aggregating QC outputs (fastp, mosdepth, etc.)

    Both modules are gated on params.tools containing their tool name.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HTML_REPORT } from '../modules/local/html_report/main'
include { MULTIQC     } from '../modules/local/multiqc/main'

workflow REPORTING {

    take:
    ch_report_inputs  // channel: [meta, clinvar_dir, pharmcat_html, clinical_vcf, cpsr_html, slivar_vcf]
    ch_multiqc_files  // channel: flat collection of QC files (fastp, mosdepth, samtools, etc.)

    main:
    ch_versions     = Channel.empty()
    ch_html_reports = Channel.empty()
    ch_multiqc_html = Channel.empty()

    //
    // MODULE 1: Per-sample consolidated HTML report
    //
    if (params.tools && params.tools.split(',').collect{ it.trim() }.contains('html_report')) {
        HTML_REPORT(ch_report_inputs)
        ch_html_reports = HTML_REPORT.out.html_report
        ch_versions     = ch_versions.mix(HTML_REPORT.out.versions)
    }

    //
    // MODULE 2: Cross-sample MultiQC aggregation
    //
    if (params.tools && params.tools.split(',').collect{ it.trim() }.contains('multiqc')) {
        MULTIQC(ch_multiqc_files.collect())
        ch_multiqc_html = MULTIQC.out.report
        ch_versions     = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    html_reports = ch_html_reports  // channel: [meta, html]
    multiqc_html = ch_multiqc_html  // channel: path(html)
    versions     = ch_versions      // channel: path(versions.yml)
}
