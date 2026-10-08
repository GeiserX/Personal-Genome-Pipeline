#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Personal Genome Pipeline — whole genome, from reads to report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Each samplesheet row starts from FASTQ (trimmed, aligned and called
    here), from a BAM or CRAM (called here), or from a VCF with an optional
    BAM or CRAM from any other caller (e.g. nf-core/sarek). From there it runs
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

include { UPSTREAM     } from './workflows/upstream'
include { PGX          } from './workflows/pgx'
include { ANNOTATION   } from './workflows/annotation'
include { CLINICAL     } from './workflows/clinical'
include { BAM_ANALYSIS } from './workflows/bam_analysis'
include { SV           } from './workflows/sv'
include { REPORTING    } from './workflows/reporting'
include { VCF_PRECHECK } from './modules/local/vcf_precheck/main'
include { CRAM_TO_BAM  } from './modules/local/cram_archive/main'
include { CRAM_ARCHIVE } from './modules/local/cram_archive/main'

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
    if (!(params.sex_check in ['fail', 'warn'])) {
        error "--sex_check must be 'fail' or 'warn', got '${params.sex_check}'."
    }
    if (!(params.freemix_warn instanceof Number) || params.freemix_warn < 0 || params.freemix_warn >= 1) {
        error "--freemix_warn must be a fraction from 0 up to (not including) 1 (FREEMIX above it is reported " +
              "as possible contamination), got '${params.freemix_warn}'."
    }
    if (params.reference.endsWith('.gz')) {
        error "--reference ${params.reference} is compressed. The tools here need a plain FASTA with a .fai index: " +
              "decompress it (gunzip, or bgzip -d) and run 'samtools faidx' on the result."
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
        'sample_qc', 'cram_archive', 'parascopy', 'y_haplogroup', 'html_report', 'multiqc',
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
        ['hla_typing',       'hla_genes',         '--hla_genes'],
        ['clinvar',          'clinvar',           '--clinvar'],
        ['clinvar',          'clinvar_index',     '--clinvar_index'],
        ['pypgx',            'pypgx_bundle',      '--pypgx_bundle'],
        ['annotsv',          'annotsv_annotations', '--annotsv_annotations'],
        ['cnvpytor',         'cnvpytor_resources', '--cnvpytor_resources'],
        ['sample_qc',        'somalier_sites',     '--somalier_sites'],
        ['sample_qc',        'verifybamid2_panel', '--verifybamid2_panel'],
        ['cyrius',           'cyrius_install',     '--cyrius_install'],
        ['parascopy',        'parascopy_data',     '--parascopy_data'],
        ['y_haplogroup',     'yleaf_data',         '--yleaf_data'],
    ]

    db_requirements.each { tool, param_name, flag ->
        if (tools_list.contains(tool) && !params[param_name]) {
            error "Tool '${tool}' is enabled in --tools but ${flag} is not set. " +
                  "Either provide ${flag} or remove '${tool}' from --tools."
        }
    }

    // --slivar_bin is accepted for one release so an old command line still
    // starts; slivar now runs from the pinned image (SLIVAR_IMAGE).
    if (params.slivar_bin) {
        log.warn "--slivar_bin is deprecated and ignored: slivar now runs from the pinned image " +
                 "(SLIVAR_IMAGE in versions.env). Remove the option; the next release drops it."
    }

    // KIR is a second T1K pass of hla_typing, against IPD-KIR
    if (params.kir && !tools_list.contains('hla_typing')) {
        error "--kir types the KIR genes in a second T1K pass of hla_typing: add 'hla_typing' to --tools."
    }
    if (params.kir && !params.kir_dat) {
        error "--kir needs --kir_dat, IPD-KIR's kir.dat (scripts/setup.sh --kir-data installs it)."
    }
    if (tools_list.contains('parascopy') && !(params.parascopy_population in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'])) {
        error "--parascopy_population must be one of AFR, AMR, EAS, EUR, SAS (the 1000 Genomes models), " +
              "got '${params.parascopy_population}'."
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
    // Columns: sample, then one starting point per row:
    //   fastq_1,fastq_2  reads: trimmed, aligned and called here
    //   bam,bam_index    a BAM without a VCF: called here
    //   cram,crai        a CRAM instead of the BAM (read with --reference,
    //                    which must be the FASTA it was written with)
    //   vcf,vcf_index    a VCF from any caller, with bam,bam_index or
    //                    cram,crai optional, and gvcf,gvcf_index optional:
    //                    the gVCF of the same calls (step 03 writes one), read
    //                    by PharmCAT and PRS as DEEPVARIANT's would be
    // and sex (male or female), required on every row that is called.
    // Rows are read and checked here, before any task starts, so a bad row
    // stops the run at once.
    def samplesheet_rows = file(params.input, checkIfExists: true).splitCsv(header: true, strip: true)
    def seen_samples = [] as Set
    samplesheet_rows.each { row ->
        if (!row.sample) {
            error "Samplesheet row without a 'sample' value. Columns found: ${row.keySet()}"
        }
        // Sanitize sample ID — used in shell commands, file paths, and HTML output
        if (!(row.sample ==~ /^[a-zA-Z0-9._-]+$/)) {
            error "Sample name '${row.sample}' contains invalid characters. Use only a-z, A-Z, 0-9, '.', '_', '-'"
        }
        // Sample ids name the output directory and key every per-sample join
        if (!seen_samples.add(row.sample)) {
            error "Sample '${row.sample}' appears more than once in ${params.input}. Each sample needs exactly one row."
        }
        // Pairs go together
        [['fastq_1', 'fastq_2'], ['bam', 'bam_index'], ['cram', 'crai'], ['vcf', 'vcf_index'], ['gvcf', 'gvcf_index']].each { a, b ->
            if (row[a] && !row[b]) {
                error "Sample '${row.sample}': '${a}' provided without '${b}'. Both are required together."
            }
            if (!row[a] && row[b]) {
                error "Sample '${row.sample}': '${b}' provided without '${a}'. Both are required together."
            }
        }
        if (row.fastq_1 && (row.bam || row.cram || row.vcf)) {
            error "Sample '${row.sample}': a row starts from FASTQ or from a BAM, CRAM or VCF, not both. " +
                  "Remove the fastq columns to use the BAM, CRAM or VCF, or the other columns to align the reads."
        }
        if (row.gvcf && !row.vcf) {
            error "Sample '${row.sample}': a gvcf goes with the vcf it belongs to. A row without a VCF is " +
                  "called here, and DEEPVARIANT writes its own gVCF: remove the gvcf columns, or add the vcf."
        }
        if (row.bam && row.cram) {
            error "Sample '${row.sample}': the row has both a BAM and a CRAM. Keep one of them."
        }
        if (row.cram && !(row.cram ==~ /.*\.cram$/)) {
            error "Sample '${row.sample}': the cram column names ${row.cram}, which does not end in .cram. " +
                  "A BAM goes in the bam column."
        }
        if (!row.fastq_1 && !row.bam && !row.cram && !row.vcf) {
            error "Sample '${row.sample}': the row has no input. Give fastq_1 and fastq_2, or bam and bam_index, " +
                  "or cram and crai, or vcf and vcf_index (with an optional bam and bam_index, or cram and crai)."
        }
        // Sex sets the chrX/chrY ploidy of DeepVariant and of ExpansionHunter
        // (whose default is female), and INDEXCOV checks it against the BAM
        if (row.sex && !(row.sex.toLowerCase() in ['male', 'female'])) {
            error "Sample '${row.sample}': sex '${row.sex}' is not recognised. Use 'male' or 'female'."
        }
        if ((row.fastq_1 || ((row.bam || row.cram) && !row.vcf)) && !row.sex) {
            error "Sample '${row.sample}': DeepVariant calls this row and needs the sample's sex: a male sample " +
                  "is called haploid on chrX and chrY outside the pseudoautosomal regions. Add a 'sex' column " +
                  "(male or female)."
        }
        if (tools_list.contains('expansion_hunter') && (row.bam || row.cram) && !row.sex) {
            error "Sample '${row.sample}': expansion_hunter needs the sample's sex to genotype chrX loci. " +
                  "Add a 'sex' column (male or female) to the samplesheet, or remove 'expansion_hunter' from --tools."
        }
    }

    ch_rows = Channel
        .fromList(samplesheet_rows)
        .map { row -> [[id: row.sample, sex: row.sex ? row.sex.toLowerCase() : null], row] }

    ch_fastq = ch_rows
        .filter { meta, row -> row.fastq_1 }
        .map { meta, row -> [meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)] }
    ch_bam_call = ch_rows
        .filter { meta, row -> row.bam && !row.vcf }
        .map { meta, row -> [meta, file(row.bam, checkIfExists: true), file(row.bam_index, checkIfExists: true)] }
    ch_bam_given = ch_rows
        .filter { meta, row -> row.bam && row.vcf }
        .map { meta, row -> [meta, file(row.bam, checkIfExists: true), file(row.bam_index, checkIfExists: true)] }
    ch_vcf_given = ch_rows
        .filter { meta, row -> row.vcf }
        .map { meta, row -> [meta, file(row.vcf, checkIfExists: true), file(row.vcf_index, checkIfExists: true)] }
    ch_gvcf_given = ch_rows
        .filter { meta, row -> row.gvcf }
        .map { meta, row -> [meta, file(row.gvcf, checkIfExists: true), file(row.gvcf_index, checkIfExists: true)] }

    // ─── Reference genome ───────────────────────────────────────────────
    ch_reference      = Channel.value(file(params.reference, checkIfExists: true))
    ch_reference_fai  = Channel.value(file("${params.reference}.fai", checkIfExists: true))
    ch_par_bed        = Channel.value(file("${projectDir}/assets/par_grch38.bed", checkIfExists: true))

    // ─── CRAM rows ──────────────────────────────────────────────────────
    // CRAM_TO_BAM writes each CRAM out as a BAM in the work directory, checked
    // against the CRAM (samtools flagstat), and the row goes on as a BAM row:
    // every BAM step then reads it as it reads any BAM. The decoding needs the
    // reference the CRAM was written with; samtools stops on another one.
    def cram_ids = samplesheet_rows.findAll { row -> row.cram }.collect { row -> row.sample } as Set
    def cram_vcf_ids = samplesheet_rows.findAll { row -> row.cram && row.vcf }.collect { row -> row.sample } as Set
    ch_cram = ch_rows
        .filter { meta, row -> row.cram }
        .map { meta, row -> [meta, file(row.cram, checkIfExists: true), file(row.crai, checkIfExists: true)] }
    CRAM_TO_BAM(ch_cram, ch_reference, ch_reference_fai)
    ch_bam_call  = ch_bam_call.mix(CRAM_TO_BAM.out.bam.filter { meta, bam, bai -> !cram_vcf_ids.contains(meta.id) })
    ch_bam_given = ch_bam_given.mix(CRAM_TO_BAM.out.bam.filter { meta, bam, bai -> cram_vcf_ids.contains(meta.id) })

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 0: UPSTREAM — FASTQ to BAM, sex check, BAM to VCF and gVCF
    // ═══════════════════════════════════════════════════════════════════
    UPSTREAM(
        ch_fastq,
        ch_bam_call,
        ch_bam_given,
        ch_reference,
        ch_reference_fai,
        ch_par_bed
    )

    // Every sample's VCF: the one given, or the one DeepVariant called; and
    // its gVCF, when the samplesheet gives one or DeepVariant wrote one
    ch_vcf_input = ch_vcf_given.mix(UPSTREAM.out.vcf)
    ch_gvcf      = ch_gvcf_given.mix(UPSTREAM.out.gvcf)

    // ─── Input check ────────────────────────────────────────────────────
    // VCF_PRECHECK reads each VCF once before any analysis. Two problems stop
    // the run here, with the fix in the message:
    //   - no contig is chr-named (1, MT): the mito haplogroup comes out empty
    //     and chrX leaks into the ROH summary, both with exit 0;
    //   - a gVCF (or a file named like one) with pharmcat selected: PharmCAT
    //     refuses it.
    // Downstream filters keep FILTER=PASS only. A VCF with no PASS record
    // (FILTER '.' everywhere) would give zero hits in every step, so stop,
    // or with --allow_unfiltered use a copy where '.' reads as PASS.
    VCF_PRECHECK(ch_vcf_input)

    ch_vcf_checked = ch_vcf_input
        .map { meta, vcf, idx -> [meta.id, meta, vcf, idx] }
        .join(VCF_PRECHECK.out.status.map { meta, status, counts, contig_style, contigs_seen, gvcf ->
            [meta.id, status, counts, contig_style, contigs_seen, gvcf]
        })
        .map { id, meta, vcf, idx, status, counts, contig_style, contigs_seen, gvcf ->
            if (contig_style == 'other') {
                error "Sample '${id}': no contig in ${vcf.name} is named the chr way (it has ${contigs_seen}...). " +
                      "The pipeline needs GRCh38 contigs named chr1 to chr22, chrX, chrY and chrM; without them the " +
                      "mito haplogroup comes out empty and the ROH summary is wrong. Rename them, then point the " +
                      "samplesheet at the new file:\n" +
                      '    for c in $(seq 1 22) X Y; do echo "$c chr$c"; done > chr_map.txt\n' +
                      '    echo "MT chrM" >> chr_map.txt\n' +
                      "    bcftools annotate --rename-chrs chr_map.txt -Oz -o ${id}.chr.vcf.gz ${vcf.name}\n" +
                      "    bcftools index -t ${id}.chr.vcf.gz\n" +
                      "Contigs outside the map (unplaced scaffolds) keep their names, and the ClinVar screen leaves " +
                      "them out. A GRCh37 file needs more than a rename: see docs/vcf-first.md."
            }
            if (gvcf == 'blocks' && tools_list.contains('pharmcat')) {
                error "Sample '${id}': ${vcf.name} is a gVCF (it has reference blocks: ALT <*>, <NON_REF> or '.' " +
                      "with INFO/END, or a ##GVCFBlock header line), and PharmCAT refuses a gVCF. Remove the " +
                      "reference blocks and give the file a name without .g.vcf or .genomic.vcf (docs/vcf-first.md " +
                      "has the commands), or run without pharmcat and cpic."
            }
            if (gvcf == 'name' && tools_list.contains('pharmcat')) {
                error "Sample '${id}': ${vcf.name} has no reference blocks, but its name contains .g.vcf or " +
                      ".genomic.vcf, and PharmCAT refuses any file named so. Rename the file and its index, for " +
                      "example to ${id}.vcf.gz and ${id}.vcf.gz.tbi, and point the samplesheet at them."
            }
            [id, meta, vcf, idx, status, counts]
        }
        .branch { id, meta, vcf, idx, status, counts ->
            pass:       status == 'pass'
            unfiltered: status == 'unfiltered'
            no_pass:    true
        }

    ch_vcf = ch_vcf_checked.pass
        .map { id, meta, vcf, idx, status, counts -> [meta, vcf, idx] }
        .mix(
            ch_vcf_checked.unfiltered
                .map { id, meta, vcf, idx, status, counts -> [id, meta, vcf.name, counts] }
                .join(VCF_PRECHECK.out.relaxed.map { meta, vcf, idx -> [meta.id, vcf, idx] })
                .map { id, meta, name, counts, vcf, idx ->
                    log.warn "Sample '${id}': no PASS record in ${name} (${counts}); --allow_unfiltered is set, " +
                             "so records with FILTER '.' are used as PASS."
                    [meta, vcf, idx]
                },
            ch_vcf_checked.no_pass
                .map { id, meta, vcf, idx, status, counts ->
                    def remedy = params.allow_unfiltered
                        ? "Filter the VCF with your caller's recommended filters; --allow_unfiltered cannot " +
                          "help here, because no record has FILTER '.'."
                        : "Filter the VCF with your caller's recommended filters, or rerun with " +
                          "--allow_unfiltered to treat FILTER '.' as PASS."
                    error "Sample '${id}': no record in ${vcf.name} has FILTER=PASS (${counts}). " +
                          "Every PASS-only step would report zero hits. ${remedy}"
                }
        )

    // Every BAM, given or aligned here, once INDEXCOV has checked its sex
    ch_bam = UPSTREAM.out.bam

    // ─── Sequence dictionary ────────────────────────────────────────────
    // Only MITO_VARIANTS (GATK) reads the sequence dictionary
    ch_reference_dict = Channel.value([])
    if (tools_list.contains('mito_variants')) {
        if (!(params.reference ==~ /.*\.(fasta|fa|fna)$/)) {
            error "mito_variants needs the reference's .dict next to it, found by replacing the FASTA extension; " +
                  "--reference ${params.reference} does not end in .fasta, .fa or .fna."
        }
        ch_reference_dict = Channel.value(
            file(params.reference.replaceAll(/\.(fasta|fa|fna)$/, '.dict'), checkIfExists: true)
        )
    }

    // ─── Optional reference databases ───────────────────────────────────
    // Empty list [] = "no file" — standard Nextflow pattern for optional path inputs.
    // Processes check truthiness (e.g., `if (myfile)`) to skip absent databases.

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
    // The panel's GRCh38 SNVs, which setup.sh --ancestry-panel writes beside
    // it: PRS_SCORE_SITES genotypes them from the gVCF, so the projection
    // sees the sites where the sample matches the reference too.
    def ancestry_sites = params.ancestry_ref ? file(params.ancestry_ref.toString().replaceFirst(/\.tar\.zst$/, '') + '_GRCh38_sites.tsv') : null
    if (ancestry_sites && !ancestry_sites.exists() && !workflow.stubRun) {
        error "ERROR: --ancestry_ref ${params.ancestry_ref} has no site list beside it (${ancestry_sites}); run scripts/setup.sh --ancestry-panel <genome_dir> again."
    }
    ch_ancestry_sites = Channel.value(ancestry_sites && ancestry_sites.exists() ? ancestry_sites : [])
    // The catalog's trait label of each PGS id (scripts/ci/check-pgs-labels.sh checks them).
    ch_pgs_labels = Channel.value(file("${projectDir}/assets/pgs_scores.tsv", checkIfExists: true))
    if (params.pgsc_calc && !file("${params.pgsc_calc}/main.nf").exists()) {
        error "ERROR: --pgsc_calc ${params.pgsc_calc} holds no main.nf; it is the pgsc_calc checkout setup.sh makes (tools/pgsc_calc-<release>)."
    }

    // ExpansionHunter variant catalog
    ch_expansion_catalog = Channel.value(params.expansion_catalog ? file(params.expansion_catalog, checkIfExists: true) : [])

    // HLA reference database (IPD-IMGT/HLA hla.dat) and the gene annotation
    // T1K takes the genes' GRCh38 coordinates from
    ch_hla_dat   = Channel.value(params.hla_dat   ? file(params.hla_dat, checkIfExists: true)   : [])
    ch_hla_genes = Channel.value(params.hla_genes ? file(params.hla_genes, checkIfExists: true) : [])

    // AnnotSV annotation directory (the biocontainer ships no annotation data)
    ch_annotsv_annotations = Channel.value(params.annotsv_annotations ? file(params.annotsv_annotations, checkIfExists: true) : [])

    // Delly exclude map (regions skipped by delly sr -x)
    ch_delly_exclude = Channel.value(params.delly_exclude ? file(params.delly_exclude, checkIfExists: true) : [])

    // Manta call regions (configManta.py --callRegions): a bgzipped BED with its .tbi beside it
    ch_manta_call_regions       = Channel.value(params.manta_call_regions ? file(params.manta_call_regions, checkIfExists: true) : [])
    ch_manta_call_regions_index = Channel.value(params.manta_call_regions ? file("${params.manta_call_regions}.tbi", checkIfExists: true) : [])

    // UCSC GRCh38 chromosome bands for TelomereHunter (-b)
    ch_cytoband = Channel.value(params.cytoband ? file(params.cytoband, checkIfExists: true) : [])

    // Sample identity and contamination (sample_qc): somalier's sites VCF and
    // the folder of VerifyBamID2's marker panel (.UD, .mu, .bed)
    ch_somalier_sites     = Channel.value(params.somalier_sites     ? file(params.somalier_sites, checkIfExists: true)     : [])
    ch_verifybamid2_panel = Channel.value(params.verifybamid2_panel ? file(params.verifybamid2_panel, checkIfExists: true) : [])

    // Opt-in steps: Cyrius's install (setup.sh --cyrius), IPD-KIR's kir.dat
    // (--kir), Parascopy's homology table and models, and its background
    // windows for a BAM that covers part of the genome
    ch_cyrius_install = Channel.value(params.cyrius_install ? file(params.cyrius_install, checkIfExists: true) : [])
    ch_kir_dat        = Channel.value(params.kir_dat ? file(params.kir_dat, checkIfExists: true) : [])
    ch_parascopy_data = Channel.value(params.parascopy_data ? file(params.parascopy_data, checkIfExists: true) : [])
    ch_parascopy_bed  = Channel.value(params.parascopy_depth_bed ? file(params.parascopy_depth_bed, checkIfExists: true) : [])
    ch_yleaf_data     = Channel.value(params.yleaf_data ? file(params.yleaf_data, checkIfExists: true) : [])

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 1: BAM_ANALYSIS — HLA (and KIR), STR, telomere, coverage,
    // mito, SMN1/SMN2, sample QC, Y haplogroup. First, because PGX reads the HLA types.
    // ═══════════════════════════════════════════════════════════════════
    BAM_ANALYSIS(
        ch_bam,
        ch_reference,
        ch_reference_fai,
        ch_reference_dict,
        ch_expansion_catalog,
        ch_hla_dat,
        ch_hla_genes,
        ch_cytoband,
        ch_somalier_sites,
        ch_verifybamid2_panel,
        ch_kir_dat,
        ch_parascopy_data,
        ch_parascopy_bed,
        ch_yleaf_data
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 2: PGX — Pharmacogenomics & ClinVar screening; PharmCAT gets
    // the HLA types and an agreed CYP2D6 call as outside calls
    // ═══════════════════════════════════════════════════════════════════
    PGX(
        ch_vcf,
        ch_reference,
        ch_reference_fai,
        ch_clinvar,
        ch_clinvar_index,
        ch_bam,
        ch_pypgx_bundle,
        ch_gvcf,
        BAM_ANALYSIS.out.hla_alleles,
        ch_cyrius_install
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 3: ANNOTATION — VEP → vcfanno → slivar / clinical_filter
    // ═══════════════════════════════════════════════════════════════════
    ANNOTATION(
        ch_vcf,
        ch_reference,
        ch_reference_fai,
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
        ch_clinvar,
        ch_clinvar_index
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 4: CLINICAL — CPSR, ROH, PRS and ancestry (pgsc_calc), mito
    // haplogroup (from the Mutect2 chrM calls of BAM_ANALYSIS when they exist)
    // ═══════════════════════════════════════════════════════════════════
    CLINICAL(
        ch_vcf,
        ch_pcgr_data,
        ch_vep_cache_cpsr,
        ch_pgs_scoring,
        ch_ancestry_ref,
        ch_ancestry_sites,
        ch_pgs_labels,
        ch_gvcf,
        ch_reference,
        ch_reference_fai,
        BAM_ANALYSIS.out.mito_vcf
    )

    // ─── CRAM archive (opt-in: cram_archive) ────────────────────────────
    // The CRAM of every BAM, checked against it, published beside it; a row
    // given as CRAM already has one. The BAM is never deleted here.
    ch_cram_versions = Channel.empty()
    if (tools_list.contains('cram_archive')) {
        CRAM_ARCHIVE(ch_bam.filter { meta, bam, bai -> !cram_ids.contains(meta.id) }, ch_reference, ch_reference_fai)
        ch_cram_versions = CRAM_ARCHIVE.out.versions
    }

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 5: SV — Structural variant calling & annotation (opt-in)
    // ═══════════════════════════════════════════════════════════════════
    SV(
        ch_bam,
        ch_reference,
        ch_reference_fai,
        ch_delly_exclude,
        ch_annotsv_annotations,
        ch_manta_call_regions,
        ch_manta_call_regions_index
    )

    // ═══════════════════════════════════════════════════════════════════
    // WORKFLOW 6: REPORTING — HTML report & MultiQC
    // ═══════════════════════════════════════════════════════════════════

    // Per-sample report inputs: the sample's VCF, and the output files of the
    // steps that ran for it, as one list (HTML_REPORT links each where
    // bin/collect_summary.py reads it). remainder: true keeps a sample a step
    // did not run for; that step's slot is null and drops out of the list. The
    // join also makes the report wait for every selected step it shows.
    ch_report_inputs = ch_vcf
        .map { meta, vcf, idx -> [meta.id, meta, vcf] }
        .join(PGX.out.clinvar_dir.map             { meta, f -> [meta.id, f] }, remainder: true)
        .join(PGX.out.pharmcat_json.map           { meta, f -> [meta.id, f] }, remainder: true)
        .join(PGX.out.pharmcat_html.map           { meta, f -> [meta.id, f] }, remainder: true)
        .join(PGX.out.cpic_phenotypes.map         { meta, f -> [meta.id, f] }, remainder: true)
        .join(PGX.out.cpic_recommendations.map    { meta, f -> [meta.id, f] }, remainder: true)
        .join(ANNOTATION.out.clinical_vcf.map     { meta, f -> [meta.id, f] }, remainder: true)
        .join(ANNOTATION.out.slivar_vcf.map       { meta, f -> [meta.id, f] }, remainder: true)
        .join(CLINICAL.out.cpsr_html.map          { meta, f -> [meta.id, f] }, remainder: true)
        .join(CLINICAL.out.roh_regions.map        { meta, f -> [meta.id, f] }, remainder: true)
        .join(CLINICAL.out.haplogroup.map         { meta, f -> [meta.id, f] }, remainder: true)
        .join(CLINICAL.out.haplocheck.map         { meta, f -> [meta.id, f] }, remainder: true)
        .join(BAM_ANALYSIS.out.y_haplogroup.map   { meta, f -> [meta.id, f] }, remainder: true)
        .join(CLINICAL.out.prs_scores.map         { meta, f -> [meta.id, f] }, remainder: true)
        .join(CLINICAL.out.ancestry_results.map   { meta, f -> [meta.id, f] }, remainder: true)
        .join(BAM_ANALYSIS.out.coverage.map       { meta, f -> [meta.id, f] }, remainder: true)
        .join(BAM_ANALYSIS.out.sample_qc.map      { meta, f -> [meta.id, f] }, remainder: true)
        .filter { items -> items[1] != null }
        .map { items -> [items[1], items[2], items[3..-1].findAll { f -> f != null }] }

    // QC files for MultiQC (mosdepth summaries)
    ch_multiqc_files = BAM_ANALYSIS.out.coverage.map { meta, f -> f }

    REPORTING(
        ch_report_inputs,
        ch_multiqc_files
    )

    // ─── Software versions ──────────────────────────────────────────────
    // Every task writes a versions.yml whose first line is the tag of the
    // image it ran in; one file collects them, one block per process.
    Channel.empty()
        .mix(
            UPSTREAM.out.versions,
            VCF_PRECHECK.out.versions,
            PGX.out.versions,
            ANNOTATION.out.versions,
            CLINICAL.out.versions,
            BAM_ANALYSIS.out.versions,
            SV.out.versions,
            REPORTING.out.versions,
            CRAM_TO_BAM.out.versions,
            ch_cram_versions
        )
        .map { f -> f.text }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)

    // ─── Completion handler ─────────────────────────────────────────────
    // Nextflow runs the handler with the script binding's variable map as its
    // delegate, and a map answers null for a name it lacks: `workflow` inside
    // the handler was null ("Cannot get property 'success' on null object" on
    // every run under 25.10), so the completion message never printed. Local
    // variables are resolved where the closure is written, so the handler
    // reads only these.
    def run_info = workflow
    def run_log  = log
    def outdir   = params.outdir
    run_info.onComplete {
        if (run_info.success) {
            run_log.info ""
            run_log.info "Pipeline completed successfully!"
            run_log.info "Results: ${outdir}"
            run_log.info ""
        } else {
            run_log.error "Pipeline failed. Check .nextflow.log for details."
        }
    }
}
