/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PGX — Pharmacogenomics & ClinVar Screening Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The VCF feeds PharmCAT and the ClinVar screen in parallel. The BAM feeds
    the CYP2D6 callers (pypgx, and Cyrius when opted in), each after the
    CYP2D6 depth check. PGX_CONSENSUS turns their calls and T1K's HLA types
    into PharmCAT's outside calls: HLA-A and HLA-B from T1K, CYP2D6 only when
    pypgx and Cyrius agree and the depth check passed. PharmCAT waits for
    them, and the CPIC lookup lists what was passed and what was held back.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHARMCAT_PREPROCESS } from '../modules/local/pharmcat/main'
include { PHARMCAT            } from '../modules/local/pharmcat/main'
include { CLINVAR_SCREEN      } from '../modules/local/clinvar_screen/main'
include { CYP2D6_DEPTH        } from '../modules/local/pgx_consensus/main'
include { PGX_CONSENSUS       } from '../modules/local/pgx_consensus/main'
include { PYPGX               } from '../modules/local/pypgx/main'
include { CYRIUS              } from '../modules/local/cyrius/main'
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
    ch_gvcf           // channel: [meta, gvcf, gvcf_index] — DEEPVARIANT's, or the samplesheet's gvcf column
    ch_hla_alleles    // channel: [meta, <id>_hla_genotype.tsv] — HLA_TYPING's
    ch_cyrius_install // channel: val(path) — setup.sh --cyrius install directory or []

    main:
    ch_versions    = Channel.empty()
    ch_clinvar_dir = Channel.empty()
    def tools = params.tools ? params.tools.split(',').collect { it.trim() } : []

    //
    // BRANCH 1: ClinVar pathogenic screen
    // Requires both --clinvar database AND 'clinvar' in --tools
    //
    if (tools.contains('clinvar') && params.clinvar) {
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
    // BRANCH 2: CYP2D6 from the BAM, after the depth check
    // pypgx: BAM-based CYP2D6/CYP2A6/GSTM1/GSTT1 + VCF-based for ~19 genes.
    // Cyrius (opt-in): the second CYP2D6 caller.
    //
    ch_depth = Channel.empty()
    if (tools.contains('pypgx') || tools.contains('cyrius')) {
        CYP2D6_DEPTH(ch_bam)
        ch_depth    = CYP2D6_DEPTH.out.regions.map { meta, q0, q1 -> [meta.id, q0, q1] }
        ch_versions = ch_versions.mix(CYP2D6_DEPTH.out.versions)
    }

    ch_pypgx_results  = Channel.empty()
    ch_pypgx_summary  = Channel.empty()
    ch_depth_check    = Channel.empty()
    if (tools.contains('pypgx')) {
        // Join BAM, VCF and the depth by sample ID so pypgx gets all per sample
        ch_bam_vcf = ch_bam
            .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
            .join(ch_vcf.map { meta, vcf, idx -> [meta.id, vcf, idx] })
            .join(ch_depth)
            .map { id, meta, bam, bai, vcf, idx, q0, q1 -> tuple(meta, bam, bai, vcf, idx, q0, q1) }

        PYPGX(
            ch_bam_vcf,
            ch_reference,
            ch_reference_fai,
            ch_pypgx_bundle
        )
        ch_pypgx_results = PYPGX.out.results
        ch_pypgx_summary = PYPGX.out.summary
        ch_depth_check   = PYPGX.out.depth_check
        ch_versions      = ch_versions.mix(PYPGX.out.versions)
    }

    ch_cyrius_results = Channel.empty()
    if (tools.contains('cyrius')) {
        CYRIUS(
            ch_bam.map { meta, bam, bai -> [meta.id, meta, bam, bai] }
                .join(ch_depth)
                .map { id, meta, bam, bai, q0, q1 -> tuple(meta, bam, bai, q0, q1) },
            ch_cyrius_install,
            Channel.value(file("${projectDir}/scripts/cyrius-constraints.txt", checkIfExists: true))
        )
        ch_cyrius_results = CYRIUS.out.cyp2d6_results
        ch_versions       = ch_versions.mix(CYRIUS.out.versions)
        if (!tools.contains('pypgx')) {
            ch_depth_check = CYRIUS.out.depth_check
        }
    }

    //
    // BRANCH 3: the outside calls PharmCAT reads (HLA from T1K, an agreed
    // CYP2D6), for each sample with a BAM. remainder: a step that did not run
    // for the sample, or whose task failed, gives [] and the table says so.
    //
    ch_outside_calls = Channel.empty()
    ch_consensus     = Channel.empty()
    if (tools.contains('pharmcat') &&
        (tools.contains('hla_typing') || tools.contains('pypgx') || tools.contains('cyrius'))) {
        ch_consensus_in = ch_bam
            .map { meta, bam, bai -> [meta.id, meta] }
            .join(ch_hla_alleles.map    { meta, f -> [meta.id, f] }, remainder: true)
            .join(ch_pypgx_summary.map  { meta, f -> [meta.id, f] }, remainder: true)
            .join(ch_cyrius_results.map { meta, f -> [meta.id, f] }, remainder: true)
            .join(ch_depth_check.map    { meta, f -> [meta.id, f] }, remainder: true)
            .filter { row -> row[1] != null }
            .map { row -> [row[1], row[2] ?: [], row[3] ?: [], row[4] ?: [], row[5] ?: []] }
        PGX_CONSENSUS(ch_consensus_in)
        ch_outside_calls = PGX_CONSENSUS.out.calls
        ch_consensus     = PGX_CONSENSUS.out.consensus
        ch_versions      = ch_versions.mix(PGX_CONSENSUS.out.versions)
    }

    //
    // BRANCH 4: PharmCAT pharmacogenomics
    //
    ch_pharmcat_html = Channel.empty()
    ch_pharmcat_json = Channel.empty()
    if (tools.contains('pharmcat')) {
        // A sample with a gVCF has its reference blocks expanded over
        // PharmCAT's regions, so a PGx position where it matches the
        // reference is a 0/0 call instead of missing; [] for the others.
        ch_pharmcat_input = ch_vcf
            .map { meta, vcf, idx -> [meta.id, meta, vcf, idx] }
            .join(ch_gvcf.map { meta, gvcf, gidx -> [meta.id, gvcf, gidx] }, remainder: true)
            .filter { row -> row[1] != null }
            .map { row -> [row[1], row[2], row[3], row[4] ?: [], row[5] ?: []] }
        PHARMCAT_PREPROCESS(ch_pharmcat_input, ch_reference, ch_reference_fai)
        // PharmCAT waits for the sample's outside calls; [] without them
        PHARMCAT(
            PHARMCAT_PREPROCESS.out.preprocessed_vcf
                .map { meta, vcf -> [meta.id, meta, vcf] }
                .join(ch_outside_calls.map { meta, f -> [meta.id, f] }, remainder: true)
                .filter { row -> row[1] != null }
                .map { row -> [row[1], row[2], row[3] ?: []] }
        )
        ch_pharmcat_html = PHARMCAT.out.html_report
        ch_pharmcat_json = PHARMCAT.out.json_report
        ch_versions = ch_versions.mix(PHARMCAT_PREPROCESS.out.versions, PHARMCAT.out.versions)
    }

    //
    // BRANCH 5: CPIC drug-gene recommendation lookup
    // Reads PharmCAT JSON output for actionable prescribing guidance. With
    // pypgx, its summary joins by sample, so the report can say when pypgx
    // called a gene PharmCAT could not; with PGX_CONSENSUS, its table says
    // which calls came from other tools and which were held back.
    //
    ch_cpic_recommendations = Channel.empty()
    ch_cpic_phenotypes      = Channel.empty()
    if (tools.contains('cpic') && tools.contains('pharmcat')) {
        // remainder: a sample whose PYPGX or PGX_CONSENSUS task failed still
        // gets its CPIC report
        ch_cpic_input = ch_pharmcat_json
            .map { meta, json -> [meta.id, meta, json] }
            .join(ch_pypgx_summary.map { meta, tsv -> [meta.id, tsv] }, remainder: true)
            .join(ch_consensus.map { meta, tsv -> [meta.id, tsv] }, remainder: true)
            .filter { row -> row[1] != null && row[2] != null }
            .map { row -> tuple(row[1], row[2], row[3] ?: [], row[4] ?: []) }
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
    cyrius_results       = ch_cyrius_results
    pgx_consensus        = ch_consensus
    cpic_recommendations = ch_cpic_recommendations
    cpic_phenotypes      = ch_cpic_phenotypes
    versions             = ch_versions
}
