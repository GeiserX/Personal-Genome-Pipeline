/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINICAL — Clinical Screening & Population Genetics Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs cancer predisposition (CPSR), runs of homozygosity, polygenic scores
    (pgsc_calc; with an ancestry panel also the sample's projection onto it,
    step 26), and the mitochondrial haplogroup (from MITO_VARIANTS' Mutect2
    calls when they exist, with the haplocheck contamination check) in
    parallel from a single input VCF.

    Each module is gated on params.tools containing the tool name.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CPSR               } from '../modules/local/cpsr/main'
include { ROH                } from '../modules/local/roh/main'
include { PRS_PREPARE        } from '../modules/local/prs/main'
include { PRS_SCORE_SITES    } from '../modules/local/prs/main'
include { PRS                } from '../modules/local/prs/main'
include { PRS_SUMMARY        } from '../modules/local/prs/main'
include { MITO_EXTRACT_CHRM  } from '../modules/local/mito_haplogroup/main'
include { MITO_HAPLOGROUP    } from '../modules/local/mito_haplogroup/main'
include { HAPLOCHECK         } from '../modules/local/mito_haplogroup/main'

workflow CLINICAL {

    take:
    ch_vcf              // channel: [meta, vcf, vcf_index]
    ch_pcgr_data        // channel: path — PCGR 2.x reference data bundle
    ch_vep_cache_cpsr   // channel: path — VEP cache for CPSR (PCGR_VEP_CACHE_RELEASE, 115)
    ch_pgs_scoring      // channel: path — PGS Catalog scoring files directory
    ch_ancestry_ref     // channel: path — pgsc_calc's ancestry panel (.tar.zst), or []
    ch_ancestry_sites   // channel: path — the panel's GRCh38 SNVs setup.sh writes beside it, or []
    ch_pgs_labels       // channel: path — assets/pgs_scores.tsv (the catalog's trait of each id)
    ch_gvcf             // channel: [meta, gvcf, gvcf_index] — DEEPVARIANT's, for the samples it called
    ch_reference        // channel: val(path) — reference FASTA (gVCF expansion)
    ch_reference_fai    // channel: val(path) — reference .fai
    ch_mito_vcf         // channel: [meta, chrM VCF] — MITO_VARIANTS' Mutect2 calls, for the samples it ran for

    main:
    ch_versions = Channel.empty()

    // Initialise output channels with empty defaults
    ch_cpsr_html          = Channel.empty()
    ch_roh_regions        = Channel.empty()
    ch_prs_scores         = Channel.empty()
    ch_ancestry_results   = Channel.empty()
    ch_haplogroup         = Channel.empty()
    ch_haplocheck         = Channel.empty()

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
    // MODULE 3: PRS — Polygenic scores (pgsc_calc), and with --ancestry_ref
    // the sample's projection onto the panel (step 26's ancestry table)
    //
    // Runs only with --pgs_scoring: without scoring files there is nothing
    // to score.
    def tools_here = params.tools ? params.tools.split(',').collect{it.trim()} : []
    if (tools_here.contains('ancestry') && !params.ancestry_ref) {
        log.warn "ancestry skipped: --ancestry_ref (pgsc_calc's panel, setup.sh --ancestry-panel) is not set."
    } else if (tools_here.contains('ancestry') && !(tools_here.contains('prs') && params.pgs_scoring)) {
        log.warn "ancestry skipped: it is computed by pgsc_calc in prs; add prs to --tools and set --pgs_scoring."
    }
    if (tools_here.contains('prs')) {
        if (!params.pgs_scoring) {
            log.warn "prs skipped: --pgs_scoring is not set."
        } else {
            def panel_name = params.ancestry_ref ? file(params.ancestry_ref).name.replaceFirst(/\.tar\.zst$/, '') : ''
            // With a gVCF the score (and panel) positions are genotyped from
            // it, so a site where the sample matches the reference counts as 0/0.
            ch_prs_input = ch_vcf
                .map { meta, vcf, idx -> [meta.id, meta, vcf, idx] }
                .join(ch_gvcf.map { meta, gvcf, gidx -> [meta.id, gvcf, gidx] }, remainder: true)
                .filter { row -> row[1] != null }
                .branch { row ->
                    gvcf: row[4] != null
                    vcf:  true
                }
            // A sample whose samplesheet gives the sex is scored on chrX too
            // (pgsc_calc's plink2 needs the sex to read chrX); one without
            // keeps the scores' autosomal rows only. PRS_PREPARE formats the
            // scores once per set the run needs.
            ch_prs_rows = ch_prs_input.gvcf.map { row -> [row[1], row[4], row[5], 'gvcf'] }
                .mix(ch_prs_input.vcf.map { row -> [row[1], row[2], row[3], 'vcf'] })
                .map { meta, f, idx, kind -> [meta.sex ? 'with_x' : 'autosomes', meta, f, idx, kind] }
            PRS_PREPARE(ch_prs_rows.map { row -> row[0] }.unique(), ch_pgs_scoring, ch_pgs_labels)
            ch_prs_sets = ch_prs_rows.combine(PRS_PREPARE.out.scores, by: 0)
            ch_sample_pgs = ch_prs_sets.map { _set, meta, _f, _idx, _kind, pgs, _alleles -> [meta.id, pgs] }
            // Both are cut to the score (and panel) positions: a small file
            // for pgsc_calc to convert.
            PRS_SCORE_SITES(
                ch_prs_sets.map { _set, meta, f, idx, kind, _pgs, alleles -> [meta, f, idx, kind, alleles] },
                ch_ancestry_sites,
                ch_reference,
                ch_reference_fai
            )
            PRS(
                PRS_SCORE_SITES.out.vcf
                    .map { meta, vcf, kind -> [meta.id, meta, vcf, kind] }
                    .join(ch_sample_pgs)
                    .map { _id, meta, vcf, kind, pgs -> [meta, vcf, kind, pgs] },
                ch_ancestry_ref
            )
            PRS_SUMMARY(
                PRS.out.results
                    .map { meta, res, kind -> [meta.id, meta, res, kind] }
                    .join(ch_sample_pgs)
                    .map { _id, meta, res, kind, pgs -> [meta, res, kind, pgs] },
                panel_name
            )
            ch_prs_scores       = PRS_SUMMARY.out.summary
            ch_ancestry_results = PRS_SUMMARY.out.ancestry
            ch_versions = ch_versions.mix(PRS_PREPARE.out.versions, PRS_SCORE_SITES.out.versions,
                                          PRS.out.versions, PRS_SUMMARY.out.versions)
        }
    }

    //
    // MITO: haplogroup from the Mutect2 chrM calls when MITO_VARIANTS ran for
    // the sample (then also the haplocheck contamination check), else from
    // the chrM records of its VCF
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('mito_haplogroup')) {
        ch_mito_in = ch_vcf
            .map { meta, vcf, idx -> [meta.id, meta, vcf, idx] }
            .join(ch_mito_vcf.map { meta, m -> [meta.id, m] }, remainder: true)
            .filter { row -> row[1] != null }
            .branch { row ->
                mutect2: row[4] != null
                vcf:     true
            }
        MITO_EXTRACT_CHRM(
            ch_mito_in.mutect2.map { row -> [row[1], row[4], [], 'mutect2'] }
                .mix(ch_mito_in.vcf.map { row -> [row[1], row[2], row[3], 'vcf'] })
        )
        MITO_HAPLOGROUP(MITO_EXTRACT_CHRM.out.chrm_vcf)
        // haplocheck needs Mutect2's allele fractions: the samples it called only
        HAPLOCHECK(
            MITO_EXTRACT_CHRM.out.chrm_vcf
                .map { meta, vcf, idx -> [meta.id, meta, vcf, idx] }
                .join(ch_mito_vcf.map { meta, m -> [meta.id, true] })
                .map { _id, meta, vcf, idx, _m -> [meta, vcf, idx] }
        )
        ch_haplogroup = MITO_HAPLOGROUP.out.haplogroup
        ch_haplocheck = HAPLOCHECK.out.report
        ch_versions   = ch_versions.mix(
            MITO_EXTRACT_CHRM.out.versions,
            MITO_HAPLOGROUP.out.versions,
            HAPLOCHECK.out.versions
        )
    }

    emit:
    cpsr_html         = ch_cpsr_html
    roh_regions       = ch_roh_regions
    prs_scores        = ch_prs_scores
    ancestry_results  = ch_ancestry_results
    haplogroup        = ch_haplogroup
    haplocheck        = ch_haplocheck     // [meta, <id>_haplocheck.txt]
    versions          = ch_versions
}
