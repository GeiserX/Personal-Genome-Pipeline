/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATION — Variant annotation, enrichment, prioritization, and clinical filtering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sequential pipeline:
      VEP → VCFANNO → SLIVAR    (sequential dependency)
                    → CLINICAL_FILTER  (branches from VCFANNO output)

    Each module is gated on params.tools containing the tool name.
    Modules that are skipped pass their input through to dependents.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VEP             } from '../modules/local/vep/main'
include { VCFANNO         } from '../modules/local/vcfanno/main'
include { SLIVAR          } from '../modules/local/slivar/main'
include { CLINICAL_FILTER } from '../modules/local/clinical_filter/main'

workflow ANNOTATION {

    take:
    ch_vcf                // channel: [meta, vcf, vcf_index]
    ch_reference          // channel: path — reference FASTA
    ch_vep_cache          // channel: path — VEP cache directory or []
    ch_cadd_snv           // channel: path — CADD SNV file or []
    ch_cadd_snv_index     // channel: path — CADD SNV index or []
    ch_cadd_indel         // channel: path — CADD indel file or []
    ch_cadd_indel_index   // channel: path — CADD indel index or []
    ch_spliceai_snv       // channel: path — SpliceAI SNV VCF or []
    ch_spliceai_snv_index // channel: path — SpliceAI SNV index or []
    ch_spliceai_indel     // channel: path — SpliceAI indel VCF or []
    ch_spliceai_indel_index // channel: path — SpliceAI indel index or []
    ch_revel              // channel: path — REVEL file or []
    ch_revel_index        // channel: path — REVEL index or []
    ch_alphamissense      // channel: path — AlphaMissense file or []
    ch_alphamissense_index // channel: path — AlphaMissense index or []
    ch_gnomad_constraint  // channel: path — gnomAD constraint TSV or []
    ch_slivar_bin         // channel: path — pre-built slivar static binary

    main:
    ch_versions = Channel.empty()

    // Track the VCF channel as it flows through the pipeline.
    // Each step either transforms it or passes it through.
    ch_current_vcf = ch_vcf

    //
    // STEP 1: VEP — Ensembl Variant Effect Predictor
    // Adds consequence, SIFT, PolyPhen, gnomAD AF, ClinVar, regulatory annotations
    //
    ch_vep_vcf = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('vep')) {
        VEP(
            ch_current_vcf,
            ch_reference,
            ch_vep_cache        )
        ch_versions = ch_versions.mix(VEP.out.versions)
        ch_vep_vcf  = VEP.out.vcf

        // Combine VCF + index for downstream
        ch_current_vcf = VEP.out.vcf
            .join(VEP.out.vcf_index)
            .map { meta, vcf, idx -> [meta, vcf, idx] }
    }

    //
    // STEP 2: VCFANNO — Enrich with CADD, SpliceAI, REVEL, AlphaMissense
    // Requires VEP output (or raw VCF if VEP is skipped)
    //
    ch_enriched_vcf = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('vcfanno')) {
        VCFANNO(
            ch_current_vcf,
            ch_cadd_snv,
            ch_cadd_snv_index,
            ch_cadd_indel,
            ch_cadd_indel_index,
            ch_spliceai_snv,
            ch_spliceai_snv_index,
            ch_spliceai_indel,
            ch_spliceai_indel_index,
            ch_revel,
            ch_revel_index,
            ch_alphamissense,
            ch_alphamissense_index        )
        ch_versions     = ch_versions.mix(VCFANNO.out.versions)
        ch_enriched_vcf = VCFANNO.out.vcf

        // Update current VCF for downstream
        ch_current_vcf = VCFANNO.out.vcf
            .join(VCFANNO.out.vcf_index)
            .map { meta, vcf, idx -> [meta, vcf, idx] }
    }

    //
    // STEP 3a: SLIVAR — Variant prioritization + compound het detection
    // Sequential from VCFANNO (or VEP, or raw VCF)
    //
    ch_slivar_vcf = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('slivar')) {
        if (!params.tools.split(',').collect{it.trim()}.contains('vep')) {
            error "slivar requires VEP-annotated input (IMPACT/CSQ fields). Add 'vep' to --tools or remove 'slivar'."
        }
        SLIVAR(
            ch_current_vcf,
            ch_gnomad_constraint,
            ch_slivar_bin        )
        ch_versions  = ch_versions.mix(SLIVAR.out.versions)
        ch_slivar_vcf = SLIVAR.out.vcf
    }

    //
    // STEP 3b: CLINICAL_FILTER — Extract clinically relevant variants
    // Branches from VCFANNO output (parallel with SLIVAR)
    //
    ch_clinical_vcf = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('clinical_filter')) {
        if (!params.tools.split(',').collect{it.trim()}.contains('vep')) {
            error "clinical_filter requires VEP-annotated input (IMPACT field). Add 'vep' to --tools or remove 'clinical_filter'."
        }
        CLINICAL_FILTER(
            ch_current_vcf
        )
        ch_versions     = ch_versions.mix(CLINICAL_FILTER.out.versions)
        ch_clinical_vcf = CLINICAL_FILTER.out.vcf
    }

    emit:
    vep_vcf      = ch_vep_vcf           // channel: [meta, vcf]
    enriched_vcf = ch_enriched_vcf      // channel: [meta, vcf]
    slivar_vcf   = ch_slivar_vcf        // channel: [meta, vcf]
    clinical_vcf = ch_clinical_vcf      // channel: [meta, vcf]
    versions     = ch_versions          // channel: versions.yml
}
