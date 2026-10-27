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
                      parents sequenced), compound heterozygotes, the
                      chromosome map of each child, the report
      FAMILY_IBD_*    the DNA each pair shares by descent (Beagle against
                      step 26b's panel, then hap-ibd): only with
                      --ancestry_fine_panel
      COUPLE_CARRIER  the carrier screen of two partners (step 39): only
                      when the pedigree marks two members as partners and
                      main.nf hands their ClinVar hits over

    Everything is published under <outdir>/family/, and each partner's copy
    of the couple table under <outdir>/<id>/couple/. docs/multi-sample.md
    describes the pedigree, the outputs and the limits; docs/39-couple.md the
    couple screen.

    Equivalent to: scripts/38-family.sh and scripts/39-couple.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FAMILY_PREPARE; GLNEXUS; GLNEXUS_VCF } from '../modules/local/glnexus/main'
include { FAMILY_RELATE; FAMILY_CHECK         } from '../modules/local/somalier/main'
include { WHATSHAP_PED; FAMILY_REPORT         } from '../modules/local/whatshap_ped/main'
include { COUPLE_CARRIER                      } from '../modules/local/couple/main'
include { FAMILY_IBD_SITES; FAMILY_IBD_PHASE; FAMILY_IBD_SEGMENTS; FAMILY_IBD_REPORT } from '../modules/local/ibd/main'

workflow FAMILY {

    take:
    ch_gvcf           // channel: [meta, g.vcf.gz, index] of each member
    ch_candidates     // channel: [meta, clinical.vcf.gz] of the members clinical_filter ran for
    ch_pedigree       // channel: val(path) — the --pedigree file
    ch_sites          // channel: val(path) — somalier's sites VCF
    ch_reference      // channel: val(path) — reference FASTA
    ch_reference_fai  // channel: val(path) — reference .fai
    ch_par_bed        // channel: val(path) — assets/par_grch38.bed
    ch_couple_outputs // channel: path — the partners' ClinVar hits, Delly VCF, ExpansionHunter VCF and SMN1
                      //   table, named <id>_*; empty when the pedigree marks no partners
    ch_carrier_genes  // channel: val(path) — assets/carrier_genes.tsv
    ch_exons          // channel: val(path) — GENCODE exons BED (--gencode_exons), or []
    ch_ibd_panel      // channel: val(path) — step 26b's panel folder (--ancestry_fine_panel); empty without it,
                      //   and then no shared-segment task runs

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

    // The shared segments: the joint call phased against the panel, then
    // hap-ibd, after the relationship check.
    FAMILY_IBD_SITES(GLNEXUS_VCF.out.vcf, ch_ibd_panel, FAMILY_CHECK.out.relationships)
    FAMILY_IBD_PHASE(FAMILY_IBD_SITES.out.sites, ch_ibd_panel)
    FAMILY_IBD_SEGMENTS(FAMILY_IBD_PHASE.out.phased, FAMILY_IBD_SITES.out.sites, ch_ibd_panel)
    FAMILY_IBD_REPORT(ch_pedigree, FAMILY_IBD_SITES.out.sites, FAMILY_IBD_PHASE.out.phased,
                      FAMILY_IBD_SEGMENTS.out.segments, ch_ibd_panel)

    // The couple screen reads each partner's own step outputs, not the joint
    // call: one task for every couple of the pedigree, once those are all
    // there. No output (no partners) means no task.
    COUPLE_CARRIER(
        ch_pedigree, ch_ids, ch_sexes, ch_couple_outputs.collect(sort: true),
        ch_carrier_genes, ch_exons, ch_par_bed
    )

    ch_versions = Channel.empty().mix(
        FAMILY_PREPARE.out.versions, GLNEXUS.out.versions, GLNEXUS_VCF.out.versions,
        FAMILY_RELATE.out.versions, FAMILY_CHECK.out.versions, WHATSHAP_PED.out.versions,
        FAMILY_REPORT.out.versions, COUPLE_CARRIER.out.versions,
        FAMILY_IBD_SITES.out.versions, FAMILY_IBD_PHASE.out.versions, FAMILY_IBD_SEGMENTS.out.versions,
        FAMILY_IBD_REPORT.out.versions
    )

    emit:
    vcf           = GLNEXUS_VCF.out.vcf              // [family.vcf.gz, tbi]
    relationships = FAMILY_CHECK.out.relationships   // relationships.tsv
    report        = FAMILY_REPORT.out.report         // family_report.md
    couple_tables = COUPLE_CARRIER.out.partner_tables  // <id>_couple_carrier.tsv of each partner (a list per task)
    versions      = ch_versions
}
