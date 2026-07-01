/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SV — Structural Variant Calling, Filtering & Annotation Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs three SV callers in parallel (Manta, Delly, CNVpytor), annotates
    Manta output with duphold depth metrics, classifies SVs via AnnotSV,
    and merges consensus calls from all callers.

    DAG:
      BAM ──┬── MANTA ──── DUPHOLD ──── ANNOTSV
            ├── DELLY
            └── CNVPYTOR
                         └── SURVIVOR_MERGE (collects all SV VCFs)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MANTA          } from '../modules/local/manta/main'
include { DELLY          } from '../modules/local/delly/main'
include { CNVPYTOR; CNVPYTOR_VCF } from '../modules/local/cnvpytor/main'
include { DUPHOLD        } from '../modules/local/duphold/main'
include { ANNOTSV        } from '../modules/local/annotsv/main'
include { SURVIVOR_MERGE } from '../modules/local/survivor_merge/main'

workflow SV {

    take:
    ch_bam            // channel: [meta, bam, bai]
    ch_reference      // channel: val(path) -- reference FASTA
    ch_reference_fai  // channel: val(path) -- reference FASTA index

    main:
    ch_versions = Channel.empty()

    // ── Parallel SV callers ─────────────────────────────────────────────

    //
    // MANTA: paired-end + split-read SV caller
    //
    ch_manta_vcf = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('manta')) {
        MANTA(ch_bam, ch_reference, ch_reference_fai)
        ch_manta_vcf = MANTA.out.diploid_sv
        ch_versions  = ch_versions.mix(MANTA.out.versions)
    }

    //
    // DELLY: paired-end + split-read + read-depth SV caller
    //
    ch_delly_vcf = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('delly')) {
        DELLY(ch_bam, ch_reference, ch_reference_fai)
        ch_delly_vcf = DELLY.out.sv_vcf
        ch_versions  = ch_versions.mix(DELLY.out.versions)
    }

    //
    // CNVPYTOR: depth-based CNV caller (orthogonal to paired-end methods)
    //
    ch_cnvpytor_calls = Channel.empty()
    ch_cnvpytor_vcf   = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('cnvpytor')) {
        // value channel so the pinned resources broadcast to every sample (not just the first)
        ch_cnvpytor_res   = Channel.value(file(params.cnvpytor_resources, checkIfExists: true))
        CNVPYTOR(ch_bam, ch_cnvpytor_res)
        ch_cnvpytor_calls = CNVPYTOR.out.cnv_calls
        CNVPYTOR_VCF(CNVPYTOR.out.raw_vcf, ch_reference_fai)
        ch_cnvpytor_vcf   = CNVPYTOR_VCF.out.cnv_vcf
        ch_versions       = ch_versions.mix(CNVPYTOR.out.versions, CNVPYTOR_VCF.out.versions)
    }

    // ── Duphold annotation (Manta -> DUPHOLD) ───────────────────────────

    //
    // DUPHOLD: depth-based SV quality annotation on Manta output
    //
    ch_duphold_vcf = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('duphold') && !params.tools.split(',').collect{it.trim()}.contains('manta')) {
        error "Tool 'duphold' requires 'manta' output — add 'manta' to --tools or remove 'duphold'."
    }
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('duphold')) {
        // Combine Manta SV VCF with BAM for duphold input
        ch_duphold_input = ch_manta_vcf
            .join(ch_bam)
            .map { meta, sv_vcf, bam, bai ->
                [meta, sv_vcf, bam, bai]
            }

        DUPHOLD(ch_duphold_input, ch_reference, ch_reference_fai)
        ch_duphold_vcf = DUPHOLD.out.annotated_vcf
        ch_versions    = ch_versions.mix(DUPHOLD.out.versions)
    }

    // ── AnnotSV classification (DUPHOLD -> ANNOTSV) ─────────────────────

    //
    // ANNOTSV: ACMG pathogenicity classification of duphold-filtered SVs
    //
    ch_annotsv_tsv = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('annotsv') && !params.tools.split(',').collect{it.trim()}.contains('duphold')) {
        error "Tool 'annotsv' requires 'duphold' (and 'manta') output — add both to --tools or remove 'annotsv'."
    }
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('annotsv')) {
        ANNOTSV(ch_duphold_vcf)
        ch_annotsv_tsv = ANNOTSV.out.annotated_tsv
        ch_versions    = ch_versions.mix(ANNOTSV.out.versions)
    }

    // ── Consensus merge (all SV VCFs -> SURVIVOR_MERGE) ─────────────────

    //
    // SURVIVOR_MERGE: bcftools-based heuristic merge of 2+ callers
    //
    ch_merged_sv = Channel.empty()
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('survivor_merge')) {
        // Collect all available SV VCFs by meta.id
        ch_all_sv_vcfs = Channel.empty()
            .mix(
                ch_manta_vcf,
                ch_delly_vcf,
                ch_cnvpytor_vcf
            )
            .groupTuple()
            .map { meta, vcfs ->
                [meta, vcfs.flatten()]
            }

        SURVIVOR_MERGE(ch_all_sv_vcfs, ch_reference_fai)
        ch_merged_sv = SURVIVOR_MERGE.out.merged_vcf
        ch_versions  = ch_versions.mix(SURVIVOR_MERGE.out.versions)
    }

    emit:
    manta_vcf      = ch_manta_vcf
    delly_vcf      = ch_delly_vcf
    cnvpytor_calls = ch_cnvpytor_calls
    duphold_vcf    = ch_duphold_vcf
    annotsv_tsv    = ch_annotsv_tsv
    merged_sv      = ch_merged_sv
    versions       = ch_versions
}
