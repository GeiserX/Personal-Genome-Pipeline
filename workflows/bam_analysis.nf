/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAM_ANALYSIS — Parallel BAM-based analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs HLA typing (with an opt-in KIR pass), repeat expansion detection,
    telomere length estimation, coverage statistics, mitochondrial variant
    calling, SMN1/SMN2 copy number (opt-in) and the sample identity and
    contamination check, and the Y haplogroup (opt-in) ALL in parallel
    from a single BAM input. The CYP2D6
    callers run in the PGX workflow, which also reads the HLA types.

    Each module is gated on params.tools containing the tool name.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { T1K_BUILD        } from '../modules/local/t1k_build/main'
include { HLA_TYPING       } from '../modules/local/hla_typing/main'
include { KIR_BUILD        } from '../modules/local/hla_typing/main'
include { KIR_TYPING       } from '../modules/local/hla_typing/main'
include { EXPANSION_HUNTER } from '../modules/local/expansion_hunter/main'
include { STRANGER         } from '../modules/local/stranger/main'
include { TELOMERE_HUNTER  } from '../modules/local/telomere_hunter/main'
include { MOSDEPTH         } from '../modules/local/mosdepth/main'
include { MITO_VARIANTS    } from '../modules/local/mito_variants/main'
include { PARASCOPY        } from '../modules/local/parascopy/main'
include { SOMALIER         } from '../modules/local/somalier/main'
include { SOMALIER_RELATE  } from '../modules/local/somalier/main'
include { SAMPLE_QC        } from '../modules/local/somalier/main'
include { VERIFYBAMID2     } from '../modules/local/verifybamid2/main'
include { Y_POSITIONS; Y_PILEUP; Y_HAPLOGROUP } from '../modules/local/y_haplogroup/main'

workflow BAM_ANALYSIS {

    take:
    ch_bam               // channel: [meta, bam, bai]
    ch_reference         // channel: val(path) — reference FASTA
    ch_reference_fai     // channel: val(path) — reference .fai index
    ch_reference_dict    // channel: val(path) — reference .dict
    ch_expansion_catalog // channel: val(path) — ExpansionHunter variant catalog JSON
    ch_hla_dat           // channel: val(path) — Pre-downloaded IPD-IMGT/HLA hla.dat
    ch_hla_genes         // channel: val(path) — gene annotation (GTF) T1K takes coordinates from
    ch_cytoband          // channel: val(path) — UCSC GRCh38 chromosome bands or []
    ch_somalier_sites    // channel: val(path) — somalier sites VCF or []
    ch_verifybamid2_panel // channel: val(path) — folder of VerifyBamID2's .UD/.mu/.bed panel or []
    ch_kir_dat           // channel: val(path) — IPD-KIR kir.dat (--kir) or []
    ch_parascopy_data    // channel: val(path) — Parascopy's homology table and models or []
    ch_parascopy_bed     // channel: val(path) — Parascopy background windows or []
    ch_yleaf_data        // channel: val(path) — Yleaf's data folder (--yleaf_data) or []

    main:
    ch_versions = Channel.empty()

    // Initialise output channels with empty defaults
    ch_hla_alleles      = Channel.empty()
    ch_expansion_vcf    = Channel.empty()
    ch_stranger_vcf     = Channel.empty()
    ch_telomere_results = Channel.empty()
    ch_coverage         = Channel.empty()
    ch_mito_vcf         = Channel.empty()
    ch_kir_genotype     = Channel.empty()
    ch_smn_copy_number  = Channel.empty()
    ch_sample_qc        = Channel.empty()
    ch_y_haplogroup     = Channel.empty()

    //
    // MODULE 1: HLA Typing (T1K)
    // Gates on: params.tools contains 'hla_typing'
    // T1K_BUILD reads only value channels, so it runs once and its index
    // serves every sample's HLA_TYPING task.
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('hla_typing')) {
        T1K_BUILD(ch_hla_dat, ch_hla_genes)
        HLA_TYPING(ch_bam, T1K_BUILD.out.seq_fa, T1K_BUILD.out.coord_fa)
        ch_hla_alleles = HLA_TYPING.out.hla_alleles
        ch_versions    = ch_versions.mix(T1K_BUILD.out.versions, HLA_TYPING.out.versions)

        // Opt-in: --kir, a second T1K pass over the KIR genes (main.nf
        // requires --kir_dat with it)
        if (params.kir) {
            KIR_BUILD(ch_kir_dat, ch_hla_genes)
            KIR_TYPING(ch_bam, KIR_BUILD.out.seq_fa, KIR_BUILD.out.coord_fa, KIR_BUILD.out.release)
            ch_kir_genotype = KIR_TYPING.out.kir_genotype
            ch_versions     = ch_versions.mix(KIR_BUILD.out.versions, KIR_TYPING.out.versions)
        }
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
    // MODULE 2b: Stranger (STR clinical annotation)
    // Gates on: params.tools contains 'stranger'
    // Requires expansion_hunter to also be in params.tools (enforced in main.nf)
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('stranger')) {
        STRANGER(ch_expansion_vcf)
        ch_stranger_vcf = STRANGER.out.vcf
        ch_versions     = ch_versions.mix(STRANGER.out.versions)
    }

    //
    // MODULE 3: TelomereHunter (telomere length)
    // Gates on: params.tools contains 'telomere_hunter'
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('telomere_hunter')) {
        if (!params.cytoband) {
            log.warn "telomere_hunter: --cytoband is not set, so TelomereHunter classifies reads by its own hg19 " +
                     "chromosome bands on GRCh38 positions. Pass UCSC's GRCh38 bands (scripts/setup.sh installs " +
                     "reference/cytoBand.hg38.txt)."
        }
        TELOMERE_HUNTER(ch_bam, ch_cytoband)
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
    // MODULE 6: Parascopy (SMN1/SMN2 copy number)
    // Gates on: params.tools contains 'parascopy' (opt-in; main.nf requires
    // --parascopy_data with it)
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('parascopy')) {
        PARASCOPY(ch_bam, ch_reference, ch_reference_fai, ch_parascopy_data, ch_parascopy_bed)
        ch_smn_copy_number = PARASCOPY.out.copy_number
        ch_versions        = ch_versions.mix(PARASCOPY.out.versions)
    }

    //
    // MODULE 7: Sample identity and contamination (somalier, VerifyBamID2)
    // Gates on: params.tools contains 'sample_qc' (main.nf requires
    // --somalier_sites and --verifybamid2_panel with it)
    // SOMALIER_RELATE runs once over every sample, so two rows that are the
    // same person are reported. A sex that disagrees with the samplesheet
    // stops the run as INDEXCOV's check does (--sex_check warn logs it);
    // contamination only warns.
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('sample_qc')) {
        SOMALIER(ch_bam, ch_reference, ch_reference_fai, ch_somalier_sites)
        // Each sample's samplesheet sex as somalier's pedigree (-9 where it
        // gives none), so somalier's own table shows the declared sex too.
        ch_somalier_ped = SOMALIER.out.extract
            .map { meta, f, id_file ->
                def code = meta.sex == 'male' ? '1' : (meta.sex == 'female' ? '2' : '-9')
                "${meta.id}\t${meta.id}\t-9\t-9\t${code}\t-9".toString()
            }
            .collectFile(name: 'samples.ped', newLine: true, sort: true)
        SOMALIER_RELATE(SOMALIER.out.extract.map { meta, f, id_file -> f }.collect(), ch_somalier_sites, ch_somalier_ped)
        VERIFYBAMID2(ch_bam, ch_reference, ch_reference_fai, ch_verifybamid2_panel)
        ch_qc_in = SOMALIER.out.extract
            .map { meta, f, id_file -> [meta.id, meta, id_file] }
            .join(VERIFYBAMID2.out.selfsm.map { meta, selfsm, marker_check -> [meta.id, selfsm, marker_check] })
            .map { _id, meta, id_file, selfsm, marker_check -> [meta, id_file, selfsm, marker_check] }
        SAMPLE_QC(ch_qc_in, SOMALIER_RELATE.out.samples, SOMALIER_RELATE.out.pairs)
        ch_sample_qc = SAMPLE_QC.out.table.map { meta, tsv ->
            def qc = [:]
            tsv.text.readLines().drop(1).each { line ->
                def kv = line.split('\t', 2)
                qc[kv[0]] = kv.size() > 1 ? kv[1] : ''
            }
            if (qc.sex_check == 'mismatch') {
                def msg = "Sample '${meta.id}': the samplesheet says sex ${meta.sex}, but somalier infers " +
                          "${qc.inferred_sex} from the reads (chrX sites ${qc.x_sites}: ${qc.x_het} heterozygous, " +
                          "${qc.x_hom_alt} homozygous ALT; chrY depth ratio ${qc.y_depth_ratio}). Either the sample " +
                          "is not the one you think, the declared sex is wrong, or the sample has a sex-chromosome " +
                          "aneuploidy. DeepVariant's chrX/chrY ploidy and ExpansionHunter take the declared sex."
                if (params.sex_check == 'warn') {
                    log.warn "${msg} --sex_check warn is set: going on with ${meta.sex}."
                } else {
                    error "${msg} Correct the sex column, or rerun with --sex_check warn to go on with the declared sex."
                }
            } else {
                log.info "Sample '${meta.id}': somalier infers ${qc.inferred_sex} from the reads; " +
                         "declared ${meta.sex ?: 'nothing'} (${qc.sex_check_reason})."
            }
            if (qc.contamination == 'warn') {
                log.warn "Sample '${meta.id}': VerifyBamID2 estimates FREEMIX ${qc.freemix}, above --freemix_warn " +
                         "${qc.freemix_warn_above}: about that share of the reads may come from another person. " +
                         "Calls, above all heterozygous ones, are less reliable; see docs/33-sample-qc.md."
            }
            if (qc.same_person_as) {
                log.warn "Sample '${meta.id}': somalier finds the same person in ${qc.same_person_as}: a duplicate " +
                         "row or a sample swap."
            }
            [meta, tsv]
        }
        ch_versions = ch_versions.mix(SOMALIER.out.versions, SOMALIER_RELATE.out.versions,
                                      VERIFYBAMID2.out.versions, SAMPLE_QC.out.versions)
    }

    //
    // MODULE 8: Y haplogroup (Yleaf), opt-in: 'y_haplogroup' in --tools
    // Male samples only. The sex is the samplesheet's, which INDEXCOV has
    // checked against the reads before any BAM reaches this workflow (a
    // mismatch stops the run unless --sex_check warn).
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('y_haplogroup')) {
        ch_y_bam = ch_bam.filter { meta, bam, bai ->
            if (meta.sex == 'male') {
                return true
            }
            log.info "Sample '${meta.id}': y_haplogroup skipped: the sample is ${meta.sex ?: 'of no declared sex'}, " +
                     "and a Y haplogroup needs a male sample."
            return false
        }
        Y_POSITIONS(ch_yleaf_data)
        Y_PILEUP(ch_y_bam, Y_POSITIONS.out.positions.first())
        Y_HAPLOGROUP(Y_PILEUP.out.pileup, ch_reference, ch_yleaf_data)
        ch_y_haplogroup = Y_HAPLOGROUP.out.haplogroup
        ch_versions     = ch_versions.mix(Y_POSITIONS.out.versions, Y_PILEUP.out.versions, Y_HAPLOGROUP.out.versions)
    }

    emit:
    hla_alleles      = ch_hla_alleles
    expansion_vcf    = ch_expansion_vcf
    stranger_vcf     = ch_stranger_vcf
    telomere_results = ch_telomere_results
    coverage         = ch_coverage
    mito_vcf         = ch_mito_vcf
    kir_genotype     = ch_kir_genotype    // [meta, <id>_kir_genotype.tsv, database_release.txt]
    smn_copy_number  = ch_smn_copy_number // [meta, <id>_smn_copy_number.tsv]
    sample_qc        = ch_sample_qc       // [meta, <id>_sample_qc.tsv]
    y_haplogroup     = ch_y_haplogroup    // [meta, <id>_y_haplogroup.txt]
    versions         = ch_versions
}
