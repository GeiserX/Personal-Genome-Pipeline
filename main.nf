#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Personal Genome Pipeline — Post-processing & Clinical Interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Accepts VCF + BAM from any upstream caller (e.g. nf-core/sarek) and runs
    pharmacogenomics, ClinVar screening, PRS, ancestry, and clinical reporting.

    https://github.com/GeiserX/Personal-Genome-Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.input) {
    error "Please provide a samplesheet with --input <samplesheet.csv>"
}

if (!params.reference) {
    error "Please provide a reference FASTA with --reference <path/to/GRCh38.fasta>"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PGX } from './workflows/pgx'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // Parse samplesheet CSV
    // Expected columns: sample,vcf,vcf_index,bam,bam_index
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.vcf) {
                error "Samplesheet must have 'sample' and 'vcf' columns. Got: ${row.keySet()}"
            }
            def meta = [id: row.sample]
            def vcf = file(row.vcf, checkIfExists: true)
            def vcf_index = file(row.vcf_index, checkIfExists: true)
            def bam = row.bam ? file(row.bam, checkIfExists: true) : []
            def bam_index = row.bam_index ? file(row.bam_index, checkIfExists: true) : []
            [meta, vcf, vcf_index, bam, bam_index]
        }
        .set { ch_input }

    // Reference genome
    ch_reference = Channel.value(file(params.reference, checkIfExists: true))

    // ClinVar database (optional)
    ch_clinvar = params.clinvar
        ? Channel.value(file(params.clinvar, checkIfExists: true))
        : Channel.value([])

    ch_clinvar_index = params.clinvar_index
        ? Channel.value(file(params.clinvar_index, checkIfExists: true))
        : Channel.value([])

    // Branch: VCF feeds both PharmCAT and ClinVar in parallel
    ch_vcf = ch_input.map { meta, vcf, vcf_index, bam, bam_index ->
        [meta, vcf, vcf_index]
    }

    // Run pharmacogenomics + ClinVar workflow
    PGX(ch_vcf, ch_reference, ch_clinvar, ch_clinvar_index)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (workflow.success) {
        log.info ""
        log.info "Pipeline completed successfully!"
        log.info "Results: ${params.outdir}"
        log.info ""
    } else {
        log.error "Pipeline failed. Check .nextflow.log for details."
    }
}
