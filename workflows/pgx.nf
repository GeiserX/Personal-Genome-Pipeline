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
    ch_gvcf           // channel: [meta, gvcf, gvcf_index] — DEEPVARIANT's, for the samples it called

    main:
    ch_versions    = Channel.empty()
    ch_clinvar_dir = Channel.empty()

    //
    // BRANCH 1: PharmCAT pharmacogenomics
    //
    ch_pharmcat_html = Channel.empty()
    ch_pharmcat_json = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('pharmcat')) {
        // A sample with a gVCF has its reference blocks expanded over
        // PharmCAT's regions, so a PGx position where it matches the
        // reference is a 0/0 call instead of missing; [] for the others.
        ch_pharmcat_input = ch_vcf
            .map { meta, vcf, idx -> [meta.id, meta, vcf, idx] }
            .join(ch_gvcf.map { meta, gvcf, gidx -> [meta.id, gvcf, gidx] }, remainder: true)
            .filter { row -> row[1] != null }
            .map { row -> [row[1], row[2], row[3], row[4] ?: [], row[5] ?: []] }
        PHARMCAT_PREPROCESS(ch_pharmcat_input, ch_reference, ch_reference_fai)
        PHARMCAT(PHARMCAT_PREPROCESS.out.preprocessed_vcf)
        ch_pharmcat_html = PHARMCAT.out.html_report
        ch_pharmcat_json = PHARMCAT.out.json_report
        ch_versions = ch_versions.mix(PHARMCAT_PREPROCESS.out.versions, PHARMCAT.out.versions)
    }

    //
    // BRANCH 2: ClinVar pathogenic screen
    // Requires both --clinvar database AND 'clinvar' in --tools
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('clinvar') && params.clinvar) {
        CLINVAR_SCREEN(
            ch_vcf,
            ch_clinvar,
            ch_clinvar_index,
            ch_reference,
            ch_reference_fai
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
        // Join BAM and VCF channels by sample ID so pypgx gets both per sample
        ch_bam_vcf = ch_bam
            .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
            .join(ch_vcf.map { meta, vcf, idx -> [meta.id, vcf, idx] })
            .map { id, meta, bam, bai, vcf, idx -> tuple(meta, bam, bai, vcf, idx) }

        PYPGX(
            ch_bam_vcf,
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
    // Reads PharmCAT JSON output for actionable prescribing guidance. With
    // pypgx, its summary joins by sample, so the report can say when pypgx
    // called a gene PharmCAT could not (CYP2D6 from read depth).
    //
    ch_cpic_recommendations = Channel.empty()
    ch_cpic_phenotypes      = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('cpic') &&
        params.tools.split(',').collect{it.trim()}.contains('pharmcat')) {
        if (params.tools.split(',').collect{it.trim()}.contains('pypgx')) {
            // remainder: a sample whose PYPGX task failed still gets its CPIC report
            ch_cpic_input = ch_pharmcat_json
                .map { meta, json -> [meta.id, meta, json] }
                .join(ch_pypgx_summary.map { meta, tsv -> [meta.id, tsv] }, remainder: true)
                .filter { row -> row[1] != null && row[2] != null }
                .map { row -> tuple(row[1], row[2], row[3] ?: []) }
        } else {
            ch_cpic_input = ch_pharmcat_json.map { meta, json -> tuple(meta, json, []) }
        }
        CPIC_LOOKUP(ch_cpic_input)
        ch_cpic_recommendations = CPIC_LOOKUP.out.recommendations
        ch_cpic_phenotypes      = CPIC_LOOKUP.out.phenotypes
        ch_versions             = ch_versions.mix(CPIC_LOOKUP.out.versions)
    }

    emit:
    clinvar_dir          = ch_clinvar_dir
    pharmcat_html        = ch_pharmcat_html
    pharmcat_json        = ch_pharmcat_json
    pypgx_results        = ch_pypgx_results
    pypgx_summary        = ch_pypgx_summary
    cpic_recommendations = ch_cpic_recommendations
    cpic_phenotypes      = ch_cpic_phenotypes
    versions             = ch_versions
}
