#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Personal Genome Pipeline — Post-processing & Clinical Interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Accepts VCF + BAM from any upstream caller (e.g. nf-core/sarek) and runs
    pharmacogenomics, variant annotation, clinical screening, BAM analysis,
    structural variant calling, and consolidated reporting.

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

// ─── Fail-fast: warn when enabled tools lack required databases ────────
def tools_list = params.tools ? params.tools.split(',').collect{it.trim()} : []

def db_requirements = [
    ['vep',              'vep_cache',      '--vep_cache'],
    ['cpsr',             'pcgr_data',      '--pcgr_data'],
    ['cpsr',             'vep_cache_cpsr', '--vep_cache_cpsr'],
    ['expansion_hunter', 'expansion_catalog', '--expansion_catalog'],
]

db_requirements.each { tool, param_name, flag ->
    if (tools_list.contains(tool) && !params[param_name]) {
        log.warn "Tool '${tool}' is enabled but ${flag} is not set — it will likely fail or produce empty output. " +
                 "Either provide ${flag} or remove '${tool}' from --tools."
    }
}

// ClinVar: paired inputs required together
if (params.clinvar && !params.clinvar_index) {
    error "When --clinvar is provided, --clinvar_index must also be provided."
}
if (!params.clinvar && params.clinvar_index) {
    error "When --clinvar_index is provided, --clinvar must also be provided."
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PGX          } from './workflows/pgx'
include { ANNOTATION   } from './workflows/annotation'
include { CLINICAL     } from './workflows/clinical'
include { BAM_ANALYSIS } from './workflows/bam_analysis'
include { SV           } from './workflows/sv'
include { REPORTING    } from './workflows/reporting'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // ─── Parse samplesheet ──────────────────────────────────────────────
    // Expected columns: sample,vcf,vcf_index,bam,bam_index
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.vcf || !row.vcf_index) {
                error "Samplesheet must have 'sample', 'vcf', and 'vcf_index' columns. Got: ${row.keySet()}"
            }
            // Sanitize sample ID — used in shell commands, file paths, and HTML output
            if (!(row.sample ==~ /^[a-zA-Z0-9._-]+$/)) {
                error "Sample name '${row.sample}' contains invalid characters. Use only a-z, A-Z, 0-9, '.', '_', '-'"
            }
            def meta = [id: row.sample]
            def vcf = file(row.vcf, checkIfExists: true)
            def vcf_index = file(row.vcf_index, checkIfExists: true)
            def bam = row.bam ? file(row.bam, checkIfExists: true) : []
            def bam_index = row.bam_index ? file(row.bam_index, checkIfExists: true) : []
            [meta, vcf, vcf_index, bam, bam_index]
        }
        .set { ch_input }

    // ─── Branch input channels ──────────────────────────────────────────
    ch_vcf = ch_input.map { meta, vcf, vcf_index, bam, bam_index ->
        [meta, vcf, vcf_index]
    }

    ch_bam = ch_input
        .filter { meta, vcf, vcf_index, bam, bam_index -> bam }
        .map { meta, vcf, vcf_index, bam, bam_index -> [meta, bam, bam_index] }

    // ─── Reference genome ───────────────────────────────────────────────
    ch_reference      = Channel.value(file(params.reference, checkIfExists: true))
    ch_reference_fai  = Channel.value(file("${params.reference}.fai"))
    ch_reference_dict = Channel.value(
        file(params.reference.replaceAll(/\.(fasta|fa)$/, '.dict'))
    )

    // ─── Optional reference databases ───────────────────────────────────
    // Empty list [] = "no file" — standard Nextflow pattern for optional path inputs.
    // Processes check truthiness (e.g., `if (myfile)`) to skip absent databases.
    def empty = file("${projectDir}/assets/stub/EMPTY")

    // ClinVar
    ch_clinvar       = params.clinvar       ? Channel.value(file(params.clinvar, checkIfExists: true))       : Channel.value([])
    ch_clinvar_index = params.clinvar_index  ? Channel.value(file(params.clinvar_index, checkIfExists: true)) : Channel.value([])

    // VEP caches
    ch_vep_cache      = Channel.value(params.vep_cache      ? file(params.vep_cache, checkIfExists: true)      : [])
    ch_vep_cache_cpsr = Channel.value(params.vep_cache_cpsr ? file(params.vep_cache_cpsr, checkIfExists: true) : [])

    // PCGR/CPSR data bundle
    ch_pcgr_data = Channel.value(params.pcgr_data ? file(params.pcgr_data, checkIfExists: true) : [])

    // PyPGx bundle
    ch_pypgx_bundle = Channel.value(params.pypgx_bundle ? file(params.pypgx_bundle, checkIfExists: true) : [])

    // Annotation score databases (CADD, SpliceAI, REVEL, AlphaMissense)
    ch_cadd_snv             = Channel.value(params.cadd_snv             ? file(params.cadd_snv, checkIfExists: true)             : [])
    ch_cadd_snv_index       = Channel.value(params.cadd_snv_index       ? file(params.cadd_snv_index, checkIfExists: true)       : [])
    ch_cadd_indel           = Channel.value(params.cadd_indel           ? file(params.cadd_indel, checkIfExists: true)           : [])
    ch_cadd_indel_index     = Channel.value(params.cadd_indel_index     ? file(params.cadd_indel_index, checkIfExists: true)     : [])
    ch_spliceai_snv         = Channel.value(params.spliceai_snv         ? file(params.spliceai_snv, checkIfExists: true)         : [])
    ch_spliceai_snv_index   = Channel.value(params.spliceai_snv_index   ? file(params.spliceai_snv_index, checkIfExists: true)   : [])
    ch_spliceai_indel       = Channel.value(params.spliceai_indel       ? file(params.spliceai_indel, checkIfExists: true)       : [])
    ch_spliceai_indel_index = Channel.value(params.spliceai_indel_index ? file(params.spliceai_indel_index, checkIfExists: true) : [])
    ch_revel                = Channel.value(params.revel                ? file(params.revel, checkIfExists: true)                : [])
    ch_revel_index          = Channel.value(params.revel_index          ? file(params.revel_index, checkIfExists: true)          : [])
    ch_alphamissense        = Channel.value(params.alphamissense        ? file(params.alphamissense, checkIfExists: true)        : [])
    ch_alphamissense_index  = Channel.value(params.alphamissense_index  ? file(params.alphamissense_index, checkIfExists: true)  : [])
    ch_gnomad_constraint    = Channel.value(params.gnomad_constraint    ? file(params.gnomad_constraint, checkIfExists: true)    : [])

    // PGS scoring & ancestry reference
    ch_pgs_scoring  = Channel.value(params.pgs_scoring  ? file(params.pgs_scoring, checkIfExists: true)  : [])
    ch_ancestry_ref = Channel.value(params.ancestry_ref ? file(params.ancestry_ref, checkIfExists: true) : [])

    // ExpansionHunter variant catalog
    ch_expansion_catalog = Channel.value(params.expansion_catalog ? file(params.expansion_catalog, checkIfExists: true) : [])

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 1: PGX — Pharmacogenomics & ClinVar screening
    // ═══════════════════════════════════════════════════════════════════
    PGX(
        ch_vcf,
        ch_reference,
        ch_reference_fai,
        ch_clinvar,
        ch_clinvar_index,
        ch_bam,
        ch_pypgx_bundle
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 2: ANNOTATION — VEP → vcfanno → slivar / clinical_filter
    // ═══════════════════════════════════════════════════════════════════
    ANNOTATION(
        ch_vcf,
        ch_reference,
        ch_vep_cache,
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
        ch_alphamissense_index,
        ch_gnomad_constraint
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 3: CLINICAL — CPSR, ROH, PRS, ancestry, mito haplogroup
    // ═══════════════════════════════════════════════════════════════════
    CLINICAL(
        ch_vcf,
        ch_pcgr_data,
        ch_vep_cache_cpsr,
        ch_pgs_scoring,
        ch_ancestry_ref
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 4: BAM_ANALYSIS — HLA, STR, telomere, coverage, mito, CYP2D6
    // ═══════════════════════════════════════════════════════════════════
    BAM_ANALYSIS(
        ch_bam,
        ch_reference,
        ch_reference_fai,
        ch_reference_dict,
        ch_expansion_catalog
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 5: SV — Structural variant calling & annotation (opt-in)
    // ═══════════════════════════════════════════════════════════════════
    SV(
        ch_bam,
        ch_reference,
        ch_reference_fai
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 6: REPORTING — HTML report & MultiQC
    // ═══════════════════════════════════════════════════════════════════

    // Build per-sample report inputs by joining available outputs.
    // Uses remainder:true so samples without a given output get null → EMPTY.
    ch_report_inputs = ch_vcf
        .map { meta, vcf, idx -> [meta.id, meta] }
        .join(PGX.out.clinvar_dir.map             { meta, f -> [meta.id, f] }, remainder: true)
        .join(PGX.out.pharmcat_html.map           { meta, f -> [meta.id, f] }, remainder: true)
        .join(ANNOTATION.out.clinical_vcf.map     { meta, f -> [meta.id, f] }, remainder: true)
        .join(CLINICAL.out.cpsr_html.map          { meta, f -> [meta.id, f] }, remainder: true)
        .join(ANNOTATION.out.slivar_vcf.map       { meta, f -> [meta.id, f] }, remainder: true)
        .map { items ->
            def meta     = items[1]
            def clinvar  = items[2] ?: empty
            def pharmcat = items[3] ?: empty
            def clinical = items[4] ?: empty
            def cpsr     = items[5] ?: empty
            def slivar   = items[6] ?: empty
            [meta, clinvar, pharmcat, clinical, cpsr, slivar]
        }

    // Collect QC files for MultiQC (versions + mosdepth summaries)
    ch_multiqc_files = Channel.empty()
        .mix(
            BAM_ANALYSIS.out.coverage.map { meta, f -> f },
            PGX.out.versions,
            ANNOTATION.out.versions,
            CLINICAL.out.versions,
            BAM_ANALYSIS.out.versions,
            SV.out.versions
        )

    REPORTING(
        ch_report_inputs,
        ch_multiqc_files
    )
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
