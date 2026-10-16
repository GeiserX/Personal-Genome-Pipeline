/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FAMILY — Related samples together (opt-in: --pedigree)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    main.nf runs this workflow only when --pedigree names two or more
    samples of the samplesheet; a run without it is untouched. Every member
    needs a gVCF: the one DEEPVARIANT writes, or the samplesheet's gvcf
    column.

      FAMILY_PREPARE  check the pedigree and each gVCF's sample name
      GLNEXUS         joint call of the gVCFs (GLnexus, DeepVariantWGS)
      GLNEXUS_VCF     the joint VCF, family.vcf.gz
      FAMILY_RELATE   somalier extract on each gVCF, relate over the family
      FAMILY_CHECK    declared against measured relatedness: stops the run
                      on a pair that disagrees (a sample swap)
      WHATSHAP_PED    phase the children of complete trios by inheritance
      FAMILY_REPORT   allele origin per child, de novo candidates (both
                      parents sequenced), compound heterozygotes, the report

    Everything is published under <outdir>/family/. docs/multi-sample.md
    describes the pedigree, the outputs and the limits.

    Equivalent to: scripts/38-family.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FAMILY_PREPARE; GLNEXUS; GLNEXUS_VCF } from '../modules/local/glnexus/main'
include { FAMILY_RELATE; FAMILY_CHECK         } from '../modules/local/somalier/main'
include { WHATSHAP_PED; FAMILY_REPORT         } from '../modules/local/whatshap_ped/main'

workflow FAMILY {

    take:
    ch_gvcf           // channel: [meta, g.vcf.gz, index] of each member
    ch_candidates     // channel: [meta, clinical.vcf.gz] of the members clinical_filter ran for
    ch_pedigree       // channel: val(path) — the --pedigree file
    ch_sites          // channel: val(path) — somalier's sites VCF
    ch_reference      // channel: val(path) — reference FASTA
    ch_reference_fai  // channel: val(path) — reference .fai
    ch_par_bed        // channel: val(path) — assets/par_grch38.bed

    main:
    // One list of the members, sorted by id, so every task sees them in
    // the same order and a rerun hits the cache
    ch_members = ch_gvcf.toSortedList { a, b -> a[0].id <=> b[0].id }
    ch_ids     = ch_members.map { rows -> rows.collect { r -> r[0].id } }
    ch_sexes   = ch_members.map { rows -> rows.collect { r -> r[0].sex ?: '' } }
    ch_files   = ch_members.map { rows -> rows.collect { r -> r[1] } }
    ch_indexes = ch_members.map { rows -> rows.collect { r -> r[2] } }

    FAMILY_PREPARE(ch_pedigree, ch_ids, ch_sexes, ch_files)
    GLNEXUS(ch_ids, ch_files, FAMILY_PREPARE.out.members)
    GLNEXUS_VCF(GLNEXUS.out.bcf)
    FAMILY_RELATE(ch_ids, ch_files, ch_indexes, ch_sites, ch_reference, ch_reference_fai, FAMILY_PREPARE.out.ped)
    FAMILY_CHECK(FAMILY_PREPARE.out.members, FAMILY_RELATE.out.pairs)
    WHATSHAP_PED(GLNEXUS_VCF.out.vcf, FAMILY_PREPARE.out.trios, FAMILY_CHECK.out.relationships)

    ch_cands   = ch_candidates.toSortedList { a, b -> a[0].id <=> b[0].id }
    ch_phased  = WHATSHAP_PED.out.vcf.map { vcf, tbi -> vcf }.ifEmpty([])
    FAMILY_REPORT(
        FAMILY_PREPARE.out.members,
        FAMILY_CHECK.out.relationships,
        GLNEXUS_VCF.out.vcf,
        ch_phased,
        ch_cands.map { rows -> rows.collect { r -> r[0].id } },
        ch_cands.map { rows -> rows.collect { r -> r[1] } },
        ch_par_bed
    )

    ch_versions = Channel.empty().mix(
        FAMILY_PREPARE.out.versions, GLNEXUS.out.versions, GLNEXUS_VCF.out.versions,
        FAMILY_RELATE.out.versions, FAMILY_CHECK.out.versions, WHATSHAP_PED.out.versions,
        FAMILY_REPORT.out.versions
    )

    emit:
    vcf           = GLNEXUS_VCF.out.vcf              // [family.vcf.gz, tbi]
    relationships = FAMILY_CHECK.out.relationships   // relationships.tsv
    report        = FAMILY_REPORT.out.report         // family_report.md
    versions      = ch_versions
}
