/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UPSTREAM — From reads to a checked BAM, a VCF and a gVCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTQ rows:   FASTP (unless --skip_trim) -> ALIGN_MINIMAP2 -> ALIGN_MARKDUP
    BAM rows without a VCF join them for calling.

    Every BAM, also the one of a VCF+BAM row, then goes through INDEXCOV,
    and nothing downstream sees a BAM before its sex check passed: a sample
    whose declared sex disagrees with the sex indexcov infers stops the run
    here (--sex_check warn logs it and goes on), before DeepVariant or
    ExpansionHunter use the declared sex.

    DEEPVARIANT calls the FASTQ and BAM-only rows; a VCF+BAM row keeps its VCF.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP          } from '../modules/local/fastp/main'
include { MINIMAP2_INDEX } from '../modules/local/align_minimap2/main'
include { ALIGN_MINIMAP2 } from '../modules/local/align_minimap2/main'
include { ALIGN_MARKDUP  } from '../modules/local/align_minimap2/main'
include { INDEXCOV       } from '../modules/local/indexcov/main'
include { DEEPVARIANT    } from '../modules/local/deepvariant/main'

workflow UPSTREAM {

    take:
    ch_fastq          // channel: [meta, fastq_1, fastq_2]  rows to align and call
    ch_bam_call       // channel: [meta, bam, bai]          BAM rows without a VCF: call them
    ch_bam_given      // channel: [meta, bam, bai]          the BAM of a VCF+BAM row
    ch_reference      // channel: val(path) — reference FASTA
    ch_reference_fai  // channel: val(path) — reference .fai
    ch_par_bed        // channel: val(path) — assets/par_grch38.bed

    main:
    ch_versions = Channel.empty()

    //
    // Read QC and trimming
    //
    ch_reads = ch_fastq
    if (!params.skip_trim) {
        FASTP(ch_fastq)
        ch_reads    = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions)
    }

    //
    // The minimap2 sr index: --minimap2_index, or <reference base>.sr.mmi
    // beside the FASTA (the name step 02 gives it), or else built once in
    // this run, and only when a FASTQ row needs it.
    //
    def ref_base = params.reference.replaceAll(/\.gz$/, '').replaceAll(/\.[^.\/]+$/, '')
    def mmi = params.minimap2_index ? file(params.minimap2_index, checkIfExists: true) : file("${ref_base}.sr.mmi")
    ch_index = Channel.value(mmi)
    if (!mmi.exists()) {
        MINIMAP2_INDEX(ch_fastq.first().combine(ch_reference).map { row -> row[-1] })
        ch_index    = MINIMAP2_INDEX.out.index.first()
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }

    //
    // Alignment, duplicate marking
    //
    ALIGN_MINIMAP2(ch_reads, ch_index)
    ALIGN_MARKDUP(ALIGN_MINIMAP2.out.sam)
    ch_versions = ch_versions.mix(ALIGN_MINIMAP2.out.versions, ALIGN_MARKDUP.out.versions)

    ch_to_call = ALIGN_MARKDUP.out.bam.mix(ch_bam_call)
    ch_all_bam = ch_to_call.mix(ch_bam_given)

    //
    // Sex check: every BAM waits for it
    //
    INDEXCOV(ch_all_bam)
    ch_versions = ch_versions.mix(INDEXCOV.out.versions)

    ch_sex_checked = INDEXCOV.out.sex.map { meta, tsv ->
        def fields   = tsv.text.readLines()[1].split('\t')
        def inferred = fields[0]
        def cn       = "CNchrX=${fields[1]}, CNchrY=${fields[2]}"
        if (meta.sex && inferred != meta.sex) {
            def msg = "Sample '${meta.id}': the samplesheet says sex ${meta.sex}, but indexcov infers ${inferred} " +
                      "from the BAM index (${cn}). Either the sample is not the one you think, the declared sex " +
                      "is wrong, or the sample has a sex-chromosome aneuploidy. DeepVariant's chrX/chrY ploidy " +
                      "and ExpansionHunter take the declared sex."
            if (params.sex_check == 'warn') {
                log.warn "${msg} --sex_check warn is set: going on with ${meta.sex}."
            } else {
                error "${msg} Correct the sex column, or rerun with --sex_check warn to go on with the declared sex."
            }
        } else {
            log.info "Sample '${meta.id}': indexcov infers ${inferred} (${cn}); declared ${meta.sex ?: 'nothing'}."
        }
        [meta.id, inferred]
    }

    ch_bam = ch_all_bam
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
        .join(ch_sex_checked)
        .map { _id, meta, bam, bai, _inferred -> [meta, bam, bai] }

    //
    // Small variant calling, on the checked BAMs of the rows that need it
    //
    ch_call = ch_bam
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
        .join(ch_to_call.map { meta, _bam, _bai -> [meta.id, true] })
        .map { _id, meta, bam, bai, _call -> [meta, bam, bai] }

    DEEPVARIANT(ch_call, ch_reference, ch_reference_fai, ch_par_bed)
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

    emit:
    bam        = ch_bam                  // [meta, bam, bai]: every BAM, after the sex check
    vcf        = DEEPVARIANT.out.vcf     // [meta, vcf, tbi]
    gvcf       = DEEPVARIANT.out.gvcf    // [meta, g.vcf.gz, tbi]
    versions   = ch_versions
}
