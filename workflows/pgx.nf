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
include { PYPGX               } from '../modules/local/pypgx/main'
include { CPIC_LOOKUP          } from '../modules/local/cpic_lookup/main'

workflow PGX {

    take:
    ch_vcf            // channel: [meta, vcf, vcf_index]
    ch_reference      // channel: val(path) — reference FASTA
    ch_reference_fai  // channel: val(path) — reference FASTA index
    ch_clinvar        // channel: val(path) — ClinVar VCF or []
    ch_clinvar_index  // channel: val(path) — ClinVar VCF index or []
    ch_bam            // channel: [meta, bam, bai]
    ch_pypgx_bundle   // channel: val(path) — pypgx-bundle directory

    main:
    ch_versions    = Channel.empty()
    ch_clinvar_dir = Channel.empty()

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
            ch_clinvar,
            ch_clinvar_index,
            ch_reference
        )
        ch_clinvar_dir = CLINVAR_SCREEN.out.isec_dir
        ch_versions    = ch_versions.mix(CLINVAR_SCREEN.out.versions)
    }

    //
    // BRANCH 3: PyPGx star allele calling with SV detection
    // BAM-based analysis for CYP2D6/CYP2A6/GSTM1/GSTT1 + VCF-based for ~19 genes
    //
    ch_pypgx_results  = Channel.empty()
    ch_pypgx_summary  = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('pypgx')) {
        PYPGX(
            ch_bam,
            ch_reference,
            ch_reference_fai,
            ch_pypgx_bundle
        )
        ch_pypgx_results = PYPGX.out.results
        ch_pypgx_summary = PYPGX.out.summary
        ch_versions      = ch_versions.mix(PYPGX.out.versions)
    }

    //
    // BRANCH 4: CPIC drug-gene recommendation lookup
    // Parses PharmCAT JSON output for actionable prescribing guidance
    //
    ch_cpic_recommendations = Channel.empty()
    ch_cpic_phenotypes      = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('cpic')) {
        CPIC_LOOKUP(PHARMCAT.out.json_report)
        ch_cpic_recommendations = CPIC_LOOKUP.out.recommendations
        ch_cpic_phenotypes      = CPIC_LOOKUP.out.phenotypes
        ch_versions             = ch_versions.mix(CPIC_LOOKUP.out.versions)
    }

    emit:
    clinvar_dir          = ch_clinvar_dir
    pharmcat_html        = PHARMCAT.out.html_report
    pharmcat_json        = PHARMCAT.out.json_report
    pypgx_results        = ch_pypgx_results
    pypgx_summary        = ch_pypgx_summary
    cpic_recommendations = ch_cpic_recommendations
    cpic_phenotypes      = ch_cpic_phenotypes
    versions             = ch_versions
}
