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

    // ─── Validate inputs ────────────────────────────────────────────────
    if (!params.input) {
        error "Please provide a samplesheet with --input <samplesheet.csv>"
    }

    if (!params.reference) {
        error "Please provide a reference FASTA with --reference <path/to/GRCh38.fasta>"
    }

    def tools_list = params.tools ? params.tools.split(',').collect { it.trim() }.findAll { it } : []

    // Every name the workflows gate on. A name outside this list would be
    // ignored silently (e.g. 'clinvar_screen' instead of 'clinvar').
    def known_tools = [
        'pharmcat', 'cpic', 'clinvar', 'pypgx',
        'vep', 'vcfanno', 'slivar', 'clinical_filter',
        'cpsr', 'roh', 'prs', 'ancestry', 'mito_haplogroup',
        'hla_typing', 'expansion_hunter', 'stranger', 'telomere_hunter', 'mosdepth', 'mito_variants', 'cyrius',
        'manta', 'delly', 'cnvpytor', 'duphold', 'annotsv', 'survivor_merge',
        'html_report', 'multiqc',
    ]
    def unknown_tools = tools_list.findAll { !known_tools.contains(it) }
    if (unknown_tools) {
        error "unknown tool${unknown_tools.size() > 1 ? 's' : ''} ${unknown_tools.join(', ')} in --tools. " +
              "Known tools: ${known_tools.join(', ')}."
    }

    // Fail-fast: stop when enabled tools lack required databases

    def db_requirements = [
        ['vep',              'vep_cache',         '--vep_cache'],
        ['cpsr',             'pcgr_data',         '--pcgr_data'],
        ['cpsr',             'vep_cache_cpsr',    '--vep_cache_cpsr'],
        ['expansion_hunter', 'expansion_catalog', '--expansion_catalog'],
        ['hla_typing',       'hla_dat',           '--hla_dat'],
        ['slivar',           'slivar_bin',        '--slivar_bin'],
        ['clinvar',          'clinvar',           '--clinvar'],
        ['clinvar',          'clinvar_index',     '--clinvar_index'],
        ['pypgx',            'pypgx_bundle',      '--pypgx_bundle'],
        ['annotsv',          'annotsv_annotations', '--annotsv_annotations'],
        ['cnvpytor',         'cnvpytor_resources', '--cnvpytor_resources'],
    ]

    db_requirements.each { tool, param_name, flag ->
        if (tools_list.contains(tool) && !params[param_name]) {
            error "Tool '${tool}' is enabled in --tools but ${flag} is not set. " +
                  "Either provide ${flag} or remove '${tool}' from --tools."
        }
    }

    // cpic requires pharmcat (it parses PharmCAT JSON output)
    if (tools_list.contains('cpic') && !tools_list.contains('pharmcat')) {
        error "Tool 'cpic' is enabled in --tools but 'pharmcat' is not. " +
              "CPIC lookup requires PharmCAT JSON output — add 'pharmcat' to --tools or remove 'cpic'."
    }

    // stranger requires expansion_hunter (it annotates EH STR VCF output)
    if (tools_list.contains('stranger') && !tools_list.contains('expansion_hunter')) {
        error "Tool 'stranger' is enabled in --tools but 'expansion_hunter' is not. " +
              "Stranger annotates ExpansionHunter VCF output — add 'expansion_hunter' to --tools or remove 'stranger'."
    }

    // survivor_merge keeps calls seen by two or more callers; with one caller
    // it would write a header-only consensus and report success.
    def sv_callers = ['manta', 'delly', 'cnvpytor'].findAll { tools_list.contains(it) }
    if (tools_list.contains('survivor_merge') && sv_callers.size() < 2) {
        error "Tool 'survivor_merge' needs at least two SV callers in --tools (manta, delly, cnvpytor); " +
              "got ${sv_callers ? sv_callers.join(', ') : 'none'}. Add another caller or remove 'survivor_merge'."
    }

    // ClinVar: paired inputs required together
    if (params.clinvar && !params.clinvar_index) {
        error "When --clinvar is provided, --clinvar_index must also be provided."
    }
    if (!params.clinvar && params.clinvar_index) {
        error "When --clinvar_index is provided, --clinvar must also be provided."
    }

    // ─── Parse samplesheet ──────────────────────────────────────────────
    // Expected columns: sample,vcf,vcf_index,bam,bam_index[,sex]
    // Rows are read and checked here, before any task starts, so a bad row
    // stops the run at once.
    def samplesheet_rows = file(params.input, checkIfExists: true).splitCsv(header: true, strip: true)
    def seen_samples = [] as Set
    samplesheet_rows.each { row ->
        if (!row.sample || !row.vcf || !row.vcf_index) {
            error "Samplesheet must have 'sample', 'vcf', and 'vcf_index' columns. Got: ${row.keySet()}"
        }
        // Sanitize sample ID — used in shell commands, file paths, and HTML output
        if (!(row.sample ==~ /^[a-zA-Z0-9._-]+$/)) {
            error "Sample name '${row.sample}' contains invalid characters. Use only a-z, A-Z, 0-9, '.', '_', '-'"
        }
        // Sample ids name the output directory and key every per-sample join
        if (!seen_samples.add(row.sample)) {
            error "Sample '${row.sample}' appears more than once in ${params.input}. Each sample needs exactly one row."
        }
        // Validate BAM/BAI are provided together
        if (row.bam && !row.bam_index) {
            error "Sample '${row.sample}': 'bam' provided without 'bam_index'. Both are required together."
        }
        if (!row.bam && row.bam_index) {
            error "Sample '${row.sample}': 'bam_index' provided without 'bam'. Both are required together."
        }
        // Sex sets the chrX ploidy for ExpansionHunter (its default is female)
        if (row.sex && !(row.sex.toLowerCase() in ['male', 'female'])) {
            error "Sample '${row.sample}': sex '${row.sex}' is not recognised. Use 'male' or 'female'."
        }
        if (tools_list.contains('expansion_hunter') && !row.sex) {
            error "Sample '${row.sample}': expansion_hunter needs the sample's sex to genotype chrX loci. " +
                  "Add a 'sex' column (male or female) to the samplesheet, or remove 'expansion_hunter' from --tools."
        }
    }

    Channel
        .fromList(samplesheet_rows)
        .map { row ->
            def meta = [id: row.sample, sex: row.sex ? row.sex.toLowerCase() : null]
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
    ch_reference_fai  = Channel.value(file("${params.reference}.fai", checkIfExists: true))
    ch_reference_dict = Channel.value(
        file(params.reference.replaceAll(/\.(fasta|fa|fna)$/, '.dict'), checkIfExists: true)
    )

    // ─── Optional reference databases ───────────────────────────────────
    // Empty list [] = "no file" — standard Nextflow pattern for optional path inputs.
    // Processes check truthiness (e.g., `if (myfile)`) to skip absent databases.
    // Per-slot sentinel files avoid Nextflow staging collisions when multiple
    // optional report inputs are all absent in the same process invocation.
    def empty_clinvar  = file("${projectDir}/assets/stub/EMPTY_CLINVAR")
    def empty_pharmcat = file("${projectDir}/assets/stub/EMPTY_PHARMCAT")
    def empty_clinical = file("${projectDir}/assets/stub/EMPTY_CLINICAL")
    def empty_cpsr     = file("${projectDir}/assets/stub/EMPTY_CPSR")
    def empty_slivar   = file("${projectDir}/assets/stub/EMPTY_SLIVAR")

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

    // HLA reference database (IPD-IMGT/HLA hla.dat)
    ch_hla_dat = Channel.value(params.hla_dat ? file(params.hla_dat, checkIfExists: true) : [])

    // Slivar static binary
    ch_slivar_bin = Channel.value(params.slivar_bin ? file(params.slivar_bin, checkIfExists: true) : [])

    // AnnotSV annotation directory (the biocontainer ships no annotation data)
    ch_annotsv_annotations = Channel.value(params.annotsv_annotations ? file(params.annotsv_annotations, checkIfExists: true) : [])

    // Delly exclude map (regions skipped by delly call -x)
    ch_delly_exclude = Channel.value(params.delly_exclude ? file(params.delly_exclude, checkIfExists: true) : [])

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
        ch_gnomad_constraint,
        ch_slivar_bin
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
        ch_expansion_catalog,
        ch_hla_dat
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 5: SV — Structural variant calling & annotation (opt-in)
    // ═══════════════════════════════════════════════════════════════════
    SV(
        ch_bam,
        ch_reference,
        ch_reference_fai,
        ch_delly_exclude,
        ch_annotsv_annotations
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
            def clinvar  = items[2] ?: empty_clinvar
            def pharmcat = items[3] ?: empty_pharmcat
            def clinical = items[4] ?: empty_clinical
            def cpsr     = items[5] ?: empty_cpsr
            def slivar   = items[6] ?: empty_slivar
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

    // ─── Completion handler ─────────────────────────────────────────────
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
}
