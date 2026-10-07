/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORTING — HTML Report Generation & MultiQC Aggregation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Produces two kinds of report:
    1. Per-sample consolidated HTML report, rendered by bin/render_report.py
       as the bash report is (QC, sample identity and contamination, variant
       counts, ClinVar, PGx, CPIC, CPSR, clinical filter, slivar, ROH, mito
       haplogroup)
    2. Cross-sample MultiQC dashboard aggregating QC outputs (fastp, mosdepth, etc.)

    Both modules are gated on params.tools containing their tool name.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HTML_REPORT } from '../modules/local/html_report/main'
include { MULTIQC     } from '../modules/local/multiqc/main'

workflow REPORTING {

    take:
    ch_report_inputs  // channel: [meta, vcf, [output files of the steps that ran for the sample]]
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
        // Guard: only run MultiQC when QC inputs exist. VCF-only runs produce
        // no mosdepth summaries and MultiQC would fail with "No analysis
        // results found" and create no output files. The skip is logged, so
        // a run with multiqc selected never ends without a word about it.
        ch_multiqc_files
            .collect()
            .ifEmpty([])
            .map { files ->
                if (files.isEmpty()) {
                    log.warn "multiqc skipped: it reads mosdepth's coverage summaries, and this run produced none " +
                             "(mosdepth needs a BAM in the samplesheet and 'mosdepth' in --tools)."
                }
                files
            }
            .filter { files -> !files.isEmpty() }
            .set { ch_multiqc_gated }

        MULTIQC(ch_multiqc_gated)
        ch_multiqc_html = MULTIQC.out.report
        ch_versions     = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    html_reports = ch_html_reports  // channel: [meta, html]
    multiqc_html = ch_multiqc_html  // channel: path(html)
    versions     = ch_versions      // channel: path(versions.yml)
}
