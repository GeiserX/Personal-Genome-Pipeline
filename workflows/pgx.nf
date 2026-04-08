/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PGX — Pharmacogenomics & ClinVar Screening Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Demonstrates the core channel-branching pattern:
    VCF input feeds BOTH PharmCAT and ClinVar screen in parallel.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHARMCAT_PREPROCESS } from '../modules/local/pharmcat/main'
include { PHARMCAT            } from '../modules/local/pharmcat/main'
include { CLINVAR_SCREEN      } from '../modules/local/clinvar_screen/main'

workflow PGX {

    take:
    ch_vcf            // channel: [meta, vcf, vcf_index]
    ch_reference      // channel: val(path) — reference FASTA
    ch_clinvar        // channel: val(path) — ClinVar VCF or []
    ch_clinvar_index  // channel: val(path) — ClinVar VCF index or []

    main:
    ch_versions = Channel.empty()

    //
    // BRANCH 1: PharmCAT pharmacogenomics
    // Runs on every sample unconditionally
    //
    PHARMCAT_PREPROCESS(ch_vcf, ch_reference)
    PHARMCAT(PHARMCAT_PREPROCESS.out.preprocessed_vcf)
    ch_versions = ch_versions.mix(PHARMCAT.out.versions)

    //
    // BRANCH 2: ClinVar pathogenic screen
    // Runs IN PARALLEL with PharmCAT — only if --clinvar is provided
    //
    if (params.clinvar) {
        CLINVAR_SCREEN(
            ch_vcf,
            ch_clinvar.first(),
            ch_clinvar_index.first()
        )
        ch_versions = ch_versions.mix(CLINVAR_SCREEN.out.versions)
    }

    emit:
    pharmcat_html = PHARMCAT.out.html_report
    pharmcat_json = PHARMCAT.out.json_report
    versions      = ch_versions
}
