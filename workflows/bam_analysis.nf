/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAM_ANALYSIS — Parallel BAM-based analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs HLA typing (with an opt-in KIR pass), repeat expansion detection,
    telomere length estimation, coverage statistics, mitochondrial variant
    calling, SMN1/SMN2 copy number (opt-in), paralog genes from PacBio HiFi
    reads (opt-in) and the sample identity and
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
include { HLA_QC_READS     } from '../modules/local/hla_typing/main'
include { HLA_QC           } from '../modules/local/hla_typing/main'
include { EXPANSION_HUNTER } from '../modules/local/expansion_hunter/main'
include { STRANGER         } from '../modules/local/stranger/main'
include { STRANGER_RFC1    } from '../modules/local/stranger/main'
include { TELOMERE_HUNTER  } from '../modules/local/telomere_hunter/main'
include { MOSDEPTH         } from '../modules/local/mosdepth/main'
include { MITO_VARIANTS    } from '../modules/local/mito_variants/main'
include { MITO_VARIANTS_ANNOTATE } from '../modules/local/mito_variants/main'
include { PARASCOPY        } from '../modules/local/parascopy/main'
include { PARAPHASE_CHECK; PARAPHASE; PARAPHASE_REPORT } from '../modules/local/paraphase/main'
include { SOMALIER         } from '../modules/local/somalier/main'
include { SOMALIER_RELATE  } from '../modules/local/somalier/main'
include { SAMPLE_QC        } from '../modules/local/somalier/main'
include { VERIFYBAMID2     } from '../modules/local/verifybamid2/main'
include { Y_POSITIONS; Y_PILEUP; Y_HAPLOGROUP } from '../modules/local/y_haplogroup/main'
include { Y_STR } from '../modules/local/y_str/main'

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
    ch_vcf               // channel: [meta, vcf, tbi] — the sample's VCF, for SAMPLE_QC's chrX allele fractions
    ch_mito_clinvar      // channel: val(path) — whole ClinVar VCF with chr names (--mito_clinvar) or []
    ch_mito_clinvar_index // channel: val(path) — its .tbi or []
    ch_gnomad_mito       // channel: val(path) — gnomAD v3.1 chrM sites (--gnomad_mito) or []

    main:
    ch_versions = Channel.empty()

    // Initialise output channels with empty defaults
    ch_hla_alleles      = Channel.empty()
    ch_expansion_vcf    = Channel.empty()
    ch_stranger_vcf     = Channel.empty()
    ch_rfc1_motif       = Channel.empty()   // [meta, <id>_rfc1_motif.tsv] when Stranger flags RFC1
    ch_telomere_results = Channel.empty()
    ch_coverage         = Channel.empty()
    ch_mito_vcf         = Channel.empty()
    ch_mito_annotated   = Channel.empty()   // [meta, <id>_chrM_annotated.tsv]
    ch_kir_genotype     = Channel.empty()
    ch_hla_optitype     = Channel.empty()   // [meta, <id>_optitype_result.tsv] with --hla_qc
    ch_hla_concordance  = Channel.empty()   // [meta, <id>_hla_concordance.tsv] when OptiType typed reads
    ch_smn_copy_number  = Channel.empty()
    ch_hba_copy_number  = Channel.empty()
    ch_paraphase_cn     = Channel.empty()
    ch_paraphase_vars   = Channel.empty()
    ch_sample_qc        = Channel.empty()
    ch_y_haplogroup     = Channel.empty()
    ch_y_str            = Channel.empty()   // [meta, <id>_y_str.tsv] (y_str, male samples)
    ch_y_str_csv        = Channel.empty()   // [meta, <id>_y_str.csv]
    ch_mosdepth_depth   = Channel.empty()   // [meta, summary, global_dist] when mosdepth runs
    ch_mosdepth_regions = Channel.empty()   // [meta, regions.bed.gz] when mosdepth runs

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

        // Opt-in: --hla_qc, OptiType as a second caller of HLA-A, -B and -C,
        // on the reads of those genes; PGX_CONSENSUS passes HLA-A and HLA-B
        // to PharmCAT only when it agrees with T1K at two fields
        if (params.hla_qc) {
            HLA_QC_READS(ch_bam, T1K_BUILD.out.coord_fa)
            HLA_QC(
                HLA_QC_READS.out.reads.map { meta, r1, r2 -> [meta.id, meta, r1, r2] }
                    .join(HLA_TYPING.out.hla_alleles.map { meta, f -> [meta.id, f] })
                    .map { id, meta, r1, r2, f -> tuple(meta, r1, r2, f) }
            )
            ch_hla_optitype    = HLA_QC.out.result
            ch_hla_concordance = HLA_QC.out.concordance
            ch_versions        = ch_versions.mix(HLA_QC_READS.out.versions, HLA_QC.out.versions)
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
    // STRANGER grades with ExpansionHunter's catalog when it carries
    // Stranger's thresholds (the STRchive one), else with Stranger's GRCh38 one.
    // STRANGER_RFC1 counts the motif of the reads at RFC1 from the sample's
    // BAM when Stranger flags RFC1 (step 9b); it writes nothing otherwise.
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('stranger')) {
        STRANGER(ch_expansion_vcf, ch_expansion_catalog)
        STRANGER_RFC1(
            STRANGER.out.vcf.map { meta, vcf -> [meta.id, meta, vcf] }
                .join(ch_bam.map { meta, bam, bai -> [meta.id, bam, bai] })
                .map { id, meta, vcf, bam, bai -> tuple(meta, vcf, bam, bai) },
            ch_expansion_catalog
        )
        ch_stranger_vcf = STRANGER.out.vcf
        ch_rfc1_motif   = STRANGER_RFC1.out.motif
        ch_versions     = ch_versions.mix(STRANGER.out.versions, STRANGER_RFC1.out.versions)
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
        ch_coverage       = MOSDEPTH.out.summary
        ch_mosdepth_depth = MOSDEPTH.out.summary.join(MOSDEPTH.out.global_dist)
        ch_mosdepth_regions = MOSDEPTH.out.regions_bed
        ch_versions       = ch_versions.mix(MOSDEPTH.out.versions)
    }

    //
    // MODULE 5: Mitochondrial variant calling (GATK Mutect2, then NuMTFilterTool)
    // Gates on: params.tools contains 'mito_variants'
    // With mosdepth, each BAM waits for its mosdepth summary and distribution:
    // NuMTFilterTool needs the median autosomal depth, as step 20 reads it from
    // step 16b. Without mosdepth the filter runs at depth 0, as step 20 does.
    // failOnMismatch: a BAM with no mosdepth row stops the run instead of
    // dropping the sample from MITO_VARIANTS.
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('mito_variants')) {
        ch_mito_in = params.tools.split(',').collect{it.trim()}.contains('mosdepth')
            ? ch_bam.join(ch_mosdepth_depth, failOnMismatch: true)
            : ch_bam.map { meta, bam, bai -> [meta, bam, bai, [], []] }
        MITO_VARIANTS(
            ch_mito_in,
            ch_reference,
            ch_reference_fai,
            ch_reference_dict        )
        // Step 20's last stage: ClinVar and gnomAD for each call
        MITO_VARIANTS_ANNOTATE(MITO_VARIANTS.out.mito_vcf, ch_mito_clinvar, ch_mito_clinvar_index, ch_gnomad_mito)
        ch_mito_vcf       = MITO_VARIANTS.out.mito_vcf
        ch_mito_annotated = MITO_VARIANTS_ANNOTATE.out.annotated
        ch_versions       = ch_versions.mix(MITO_VARIANTS.out.versions, MITO_VARIANTS_ANNOTATE.out.versions)
    }

    //
    // MODULE 6: Parascopy (SMN1/SMN2 copy number)
    // Gates on: params.tools contains 'parascopy' (opt-in; main.nf requires
    // --parascopy_data with it)
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('parascopy')) {
        PARASCOPY(ch_bam, ch_reference, ch_reference_fai, ch_parascopy_data, ch_parascopy_bed)
        ch_smn_copy_number = PARASCOPY.out.copy_number
        ch_hba_copy_number = PARASCOPY.out.hba
        ch_versions        = ch_versions.mix(PARASCOPY.out.versions)
    }

    //
    // MODULE 6b: Paraphase (paralog genes from PacBio HiFi reads)
    // Gates on: params.tools contains 'paraphase' (opt-in). PARAPHASE_CHECK
    // passes a BAM on only when every read group says PL:PACBIO; the three
    // processes are report-only (nextflow.config), so another BAM fails the
    // check and the run goes on.
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('paraphase')) {
        PARAPHASE_CHECK(ch_bam)
        ch_paraphase_bam = ch_bam.join(PARAPHASE_CHECK.out.ok).map { meta, bam, bai, _ok -> [meta, bam, bai] }
        PARAPHASE(ch_paraphase_bam, ch_reference, ch_reference_fai)
        PARAPHASE_REPORT(PARAPHASE.out.results)
        ch_paraphase_cn   = PARAPHASE_REPORT.out.copy_number
        ch_paraphase_vars = PARAPHASE_REPORT.out.variants
        ch_versions     = ch_versions.mix(PARAPHASE_CHECK.out.versions, PARAPHASE.out.versions,
                                          PARAPHASE_REPORT.out.versions)
    }

    //
    // MODULE 7: Sample identity and contamination (somalier, VerifyBamID2)
    // Gates on: params.tools contains 'sample_qc' (main.nf requires
    // --somalier_sites and --verifybamid2_panel with it)
    // SOMALIER_RELATE runs once over every sample, so two rows that are the
    // same person are reported. A sex that disagrees with the samplesheet
    // stops the run as INDEXCOV's check does (--sex_check warn logs it);
    // contamination only warns. SAMPLE_QC also gets mosdepth's windows and
    // the VCF, when they exist, for the mosaic loss of chrY or chrX.
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
            .first()   // a value channel, so every SAMPLE_QC task gets the tables
        SOMALIER_RELATE(SOMALIER.out.extract.map { meta, f, id_file -> f }.collect(), ch_somalier_sites, ch_somalier_ped)
        VERIFYBAMID2(ch_bam, ch_reference, ch_reference_fai, ch_verifybamid2_panel)
        ch_qc_in = SOMALIER.out.extract
            .map { meta, f, id_file -> [meta.id, meta, id_file] }
            .join(VERIFYBAMID2.out.selfsm.map { meta, selfsm, marker_check -> [meta.id, selfsm, marker_check] })
            .join(ch_mosdepth_regions.map { meta, f -> [meta.id, f] }, remainder: true)
            .join(ch_vcf.map { meta, vcf, _tbi -> [meta.id, vcf] }, remainder: true)
            .filter { items -> items[1] != null }   // a VCF-only row has no BAM
            .map { _id, meta, id_file, selfsm, marker_check, regions, vcf ->
                [meta, id_file, selfsm, marker_check, regions ?: [], vcf ?: []] }
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
                          "aneuploidy. DeepVariant's chrX/chrY ploidy, ExpansionHunter and Delly take the declared sex."
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
        Y_PILEUP(ch_y_bam, Y_POSITIONS.out.positions.first(), Y_POSITIONS.out.tree_positions.first())
        Y_HAPLOGROUP(Y_PILEUP.out.pileup, ch_reference, ch_yleaf_data)
        ch_y_haplogroup = Y_HAPLOGROUP.out.haplogroup
        ch_versions     = ch_versions.mix(Y_POSITIONS.out.versions, Y_PILEUP.out.versions, Y_HAPLOGROUP.out.versions)
    }

    //
    // MODULE 9: Y-STR repeat counts (HipSTR), opt-in: 'y_str' in --tools
    // Male samples only, on the sex INDEXCOV has checked, as for MODULE 8.
    //
    if (params.tools && params.tools.split(',').collect{it.trim()}.contains('y_str')) {
        ch_ystr_bam = ch_bam.filter { meta, bam, bai ->
            if (meta.sex == 'male') {
                return true
            }
            log.info "Sample '${meta.id}': y_str skipped: the sample is ${meta.sex ?: 'of no declared sex'}, " +
                     "and Y-STRs need a male sample."
            return false
        }
        Y_STR(ch_ystr_bam, ch_reference, ch_reference_fai,
              Channel.value(file("${projectDir}/assets/ystr_panel.tsv", checkIfExists: true)))
        ch_y_str     = Y_STR.out.table
        ch_y_str_csv = Y_STR.out.csv
        ch_versions  = ch_versions.mix(Y_STR.out.versions)
    }

    emit:
    hla_alleles      = ch_hla_alleles
    expansion_vcf    = ch_expansion_vcf
    stranger_vcf     = ch_stranger_vcf
    rfc1_motif       = ch_rfc1_motif      // [meta, <id>_rfc1_motif.tsv]
    telomere_results = ch_telomere_results
    coverage         = ch_coverage
    mito_vcf         = ch_mito_vcf
    mito_annotated   = ch_mito_annotated  // [meta, <id>_chrM_annotated.tsv]
    kir_genotype     = ch_kir_genotype    // [meta, <id>_kir_genotype.tsv, database_release.txt]
    hla_optitype     = ch_hla_optitype    // [meta, <id>_optitype_result.tsv] (--hla_qc)
    hla_concordance  = ch_hla_concordance // [meta, <id>_hla_concordance.tsv] (--hla_qc, reads typed)
    smn_copy_number  = ch_smn_copy_number // [meta, <id>_smn_copy_number.tsv]
    hba_copy_number  = ch_hba_copy_number // [meta, <id>_hba_copy_number.tsv]
    paraphase_cn     = ch_paraphase_cn    // [meta, <id>_paraphase_copy_number.tsv]
    paraphase_vars   = ch_paraphase_vars  // [meta, <id>_paraphase_variants.tsv]
    sample_qc        = ch_sample_qc       // [meta, <id>_sample_qc.tsv]
    y_haplogroup     = ch_y_haplogroup    // [meta, <id>_y_haplogroup.txt]
    y_str            = ch_y_str           // [meta, <id>_y_str.tsv]
    y_str_csv        = ch_y_str_csv       // [meta, <id>_y_str.csv]
    versions         = ch_versions
}
