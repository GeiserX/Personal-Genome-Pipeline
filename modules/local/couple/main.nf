/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COUPLE_CARRIER — Carrier screen of two partners
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    For every two samples the pedigree marks as partners, bin/couple_report.py
    screen (the rule scripts/39-couple.sh runs too) collects each partner's
    findings in the genes of assets/carrier_genes.tsv: ClinVar pathogenic and
    likely pathogenic hits (CLINVAR_SCREEN), Delly deletions over an exon
    of 50 bp or more (with --gencode_exons; one matching a common deletion of
    AnnotSV's benign loss set is set aside, the set read from
    --annotsv_annotations when it is there), SMN1 copy number (PARASCOPY) and
    the FMR1 and FXN repeats (EXPANSION_HUNTER). It writes couple_<a>_<b>.tsv and .md, published
    under the family's folder, <outdir>/<family>/ (--family_name, default family/), and one copy of the table per partner,
    <id>_couple_carrier.tsv, published under <outdir>/<id>/couple/ for that
    partner's report.

    The step outputs come as one list; each file is matched to its sample by
    its name (<id>_clinvar_hits.vcf, <id>_sv.vcf.gz, <id>_eh.vcf,
    <id>_smn_copy_number.tsv).

    Equivalent to: scripts/39-couple.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COUPLE_CARRIER {
    label 'process_single'

    publishDir "${params.outdir}", mode: params.publish_dir_mode, saveAs: { fn ->
        fn == 'versions.yml' ? null :
        fn.endsWith('_couple_carrier.tsv') ? "${fn.substring(0, fn.length() - 19)}/couple/${fn}" : "${params.family_name ?: 'family'}/${fn}"
    }

    input:
    path(pedigree, stageAs: 'pedigree.tsv')  // the family's own rows (FAMILY_PREPARE's pedigree.tsv)
    val(ids)             // the family members' samplesheet ids, sorted
    val(sexes)           // each id's samplesheet sex ('' when the row gives none), in the same order
    path(outputs, stageAs: 'in/*')
    path(genes, stageAs: 'carrier_genes.tsv')
    path(exons)          // GENCODE exons BED, or [] without --gencode_exons
    path(benign_sv)      // AnnotSV's benign_Loss_SV_GRCh38.sorted.bed, or [] without --annotsv_annotations
    path(par_bed)

    output:
    path "couple_*.tsv",          emit: table
    path "couple_*.md",           emit: report
    path "*_couple_carrier.tsv",  emit: partner_tables
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def files = outputs instanceof List ? outputs : [outputs]
    def kinds = [clinvar: '_clinvar_hits.vcf', delly: '_sv.vcf.gz', eh: '_eh.vcf', smn: '_smn_copy_number.tsv']
    def args = files.collect { f ->
        // A staged input's name is its stage path (in/<file>): keep the file name only
        def n = f.name.tokenize('/')[-1]
        def k = kinds.find { kind, suffix -> n.endsWith(suffix) && ids.contains(n.substring(0, n.length() - suffix.length())) }
        k ? "--${k.key} ${n.substring(0, n.length() - k.value.length())}=in/${n}" : ''
    }.findAll { a -> a }.sort().join(' ')
    def sx = [ids, sexes].transpose().findAll { p -> p[1] }.collect { p -> "--declared-sex ${p[0]}=${p[1]}" }.join(' ')
    def ex = exons ? "--exons ${exons}" : ''
    def bn = benign_sv ? "--benign-sv ${benign_sv}" : ''
    """
    couple_report.py screen --pedigree pedigree.tsv --genes carrier_genes.tsv --par-bed ${par_bed} \\
        ${args} ${ex} ${bn} ${sx} --out-dir .

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    def pair = ids.take(2)
    """
    touch couple_${pair.join('_')}.tsv couple_${pair.join('_')}.md
    ${pair.collect { i -> "touch ${i}_couple_carrier.tsv" }.join('\n    ')}
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
