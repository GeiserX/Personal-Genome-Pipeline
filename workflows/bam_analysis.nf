/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAM_ANALYSIS — Parallel BAM-based analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs HLA typing, repeat expansion detection, telomere length estimation,
    coverage statistics, mitochondrial variant calling, and CYP2D6 star alleles
    ALL in parallel from a single BAM input.

    Each module is gated on params.tools containing the tool name.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HLA_TYPING       } from '../modules/local/hla_typing/main'
include { EXPANSION_HUNTER } from '../modules/local/expansion_hunter/main'
include { TELOMERE_HUNTER  } from '../modules/local/telomere_hunter/main'
include { MOSDEPTH         } from '../modules/local/mosdepth/main'
include { MITO_VARIANTS    } from '../modules/local/mito_variants/main'
include { CYRIUS           } from '../modules/local/cyrius/main'

workflow BAM_ANALYSIS {

    take:
    ch_bam               // channel: [meta, bam, bai]
    ch_reference         // channel: val(path) — reference FASTA
    ch_reference_fai     // channel: val(path) — reference .fai index
    ch_reference_dict    // channel: val(path) — reference .dict
    ch_expansion_catalog // channel: val(path) — ExpansionHunter variant catalog JSON
    ch_hla_dat           // channel: val(path) — Pre-downloaded IPD-IMGT/HLA hla.dat

    main:
    ch_versions = Channel.empty()

    // Initialise output channels with empty defaults
    ch_hla_alleles      = Channel.empty()
    ch_expansion_vcf    = Channel.empty()
    ch_telomere_results = Channel.empty()
    ch_coverage         = Channel.empty()
    ch_mito_vcf         = Channel.empty()
    ch_cyrius_results   = Channel.empty()

    //
    // MODULE 1: HLA Typing (T1K)
    // Gates on: params.tools contains 'hla_typing'
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('hla_typing')) {
        HLA_TYPING(ch_bam, ch_reference, ch_hla_dat)
        ch_hla_alleles = HLA_TYPING.out.hla_alleles
        ch_versions    = ch_versions.mix(HLA_TYPING.out.versions)
    }

    //
    // MODULE 2: ExpansionHunter (STR expansions)
    // Gates on: params.tools contains 'expansion_hunter'
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('expansion_hunter')) {
        EXPANSION_HUNTER(
            ch_bam,
            ch_reference,
            ch_reference_fai,
            ch_expansion_catalog        )
        ch_expansion_vcf = EXPANSION_HUNTER.out.vcf
        ch_versions      = ch_versions.mix(EXPANSION_HUNTER.out.versions)
    }

    //
    // MODULE 3: TelomereHunter (telomere length)
    // Gates on: params.tools contains 'telomere_hunter'
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('telomere_hunter')) {
        TELOMERE_HUNTER(ch_bam)
        ch_telomere_results = TELOMERE_HUNTER.out.telomere_results
        ch_versions         = ch_versions.mix(TELOMERE_HUNTER.out.versions)
    }

    //
    // MODULE 4: mosdepth (coverage statistics)
    // Gates on: params.tools contains 'mosdepth'
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('mosdepth')) {
        MOSDEPTH(ch_bam)
        ch_coverage = MOSDEPTH.out.summary
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    }

    //
    // MODULE 5: Mitochondrial variant calling (GATK Mutect2)
    // Gates on: params.tools contains 'mito_variants'
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('mito_variants')) {
        MITO_VARIANTS(
            ch_bam,
            ch_reference,
            ch_reference_fai,
            ch_reference_dict        )
        ch_mito_vcf = MITO_VARIANTS.out.mito_vcf
        ch_versions = ch_versions.mix(MITO_VARIANTS.out.versions)
    }

    //
    // MODULE 6: Cyrius (CYP2D6 star alleles)
    // Gates on: params.tools contains 'cyrius'
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('cyrius')) {
        CYRIUS(ch_bam)
        ch_cyrius_results = CYRIUS.out.cyp2d6_results
        ch_versions       = ch_versions.mix(CYRIUS.out.versions)
    }

    emit:
    hla_alleles      = ch_hla_alleles
    expansion_vcf    = ch_expansion_vcf
    telomere_results = ch_telomere_results
    coverage         = ch_coverage
    mito_vcf         = ch_mito_vcf
    cyrius_results   = ch_cyrius_results
    versions         = ch_versions
}
