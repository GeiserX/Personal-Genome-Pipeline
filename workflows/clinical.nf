/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINICAL — Clinical Screening & Population Genetics Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs cancer predisposition (CPSR), runs of homozygosity, polygenic risk scores,
    ancestry PCA, and mitochondrial haplogroup classification in parallel from a
    single input VCF.

    Each module is gated on params.tools containing the tool name.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CPSR               } from '../modules/local/cpsr/main'
include { ROH                } from '../modules/local/roh/main'
include { PRS                } from '../modules/local/prs/main'
include { ANCESTRY           } from '../modules/local/ancestry/main'
include { MITO_EXTRACT_CHRM  } from '../modules/local/mito_haplogroup/main'
include { MITO_HAPLOGROUP    } from '../modules/local/mito_haplogroup/main'

workflow CLINICAL {

    take:
    ch_vcf              // channel: [meta, vcf, vcf_index]
    ch_pcgr_data        // channel: path — PCGR 2.x reference data bundle
    ch_vep_cache_cpsr   // channel: path — VEP 113 cache for CPSR
    ch_pgs_scoring      // channel: path — PGS Catalog scoring files directory
    ch_ancestry_ref     // channel: path — ancestry reference panel

    main:
    ch_versions = Channel.empty()

    // Initialise output channels with empty defaults
    ch_cpsr_html          = Channel.empty()
    ch_roh_regions        = Channel.empty()
    ch_prs_scores         = Channel.empty()
    ch_ancestry_results   = Channel.empty()
    ch_haplogroup         = Channel.empty()

    //
    // MODULE 1: CPSR — Cancer predisposition screening
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('cpsr')) {
        CPSR(ch_vcf, ch_pcgr_data, ch_vep_cache_cpsr)
        ch_cpsr_html = CPSR.out.html_report
        ch_versions  = ch_versions.mix(CPSR.out.versions)
    }

    //
    // MODULE 2: ROH — Runs of homozygosity
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('roh')) {
        ROH(ch_vcf)
        ch_roh_regions = ROH.out.roh_regions
        ch_versions    = ch_versions.mix(ROH.out.versions)
    }

    //
    // MODULE 3: PRS — Polygenic risk scores
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('prs')) {
        PRS(ch_vcf, ch_pgs_scoring)
        ch_prs_scores = PRS.out.scores
        ch_versions   = ch_versions.mix(PRS.out.versions)
    }

    //
    // MODULE 4: ANCESTRY — Population PCA
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('ancestry')) {
        ANCESTRY(ch_vcf, ch_ancestry_ref)
        ch_ancestry_results = ANCESTRY.out.ancestry_tsv
        ch_versions         = ch_versions.mix(ANCESTRY.out.versions)
    }

    //
    // MODULES 5+6: MITO — chrM extraction then haplogroup classification (sequential chain)
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('mito_haplogroup')) {
        MITO_EXTRACT_CHRM(ch_vcf)
        MITO_HAPLOGROUP(MITO_EXTRACT_CHRM.out.chrm_vcf)
        ch_haplogroup = MITO_HAPLOGROUP.out.haplogroup
        ch_versions   = ch_versions.mix(
            MITO_EXTRACT_CHRM.out.versions,
            MITO_HAPLOGROUP.out.versions
        )
    }

    emit:
    cpsr_html         = ch_cpsr_html
    roh_regions       = ch_roh_regions
    prs_scores        = ch_prs_scores
    ancestry_results  = ch_ancestry_results
    haplogroup        = ch_haplogroup
    versions          = ch_versions
}
