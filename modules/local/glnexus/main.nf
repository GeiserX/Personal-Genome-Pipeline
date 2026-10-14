/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FAMILY_PREPARE, GLNEXUS, GLNEXUS_VCF — The family's joint VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FAMILY_PREPARE reads the pedigree (--pedigree) with bin/family_report.py
    prepare, the rule scripts/38-family.sh runs too: it stops the run on a
    pedigree problem, or on a gVCF whose sample column is not the
    samplesheet's id, before the joint call starts. It writes members.tsv,
    somalier's pedigree (family.ped) and WhatsHap's trios (trios.ped).

    GLNEXUS joint-calls the members' DeepVariant gVCFs with GLnexus's
    DeepVariantWGS preset into one BCF, and GLNEXUS_VCF writes it out as
    family.vcf.gz with its index. Each gVCF is linked as <id>.g.vcf.gz, so
    the sample order of the joint VCF is the ids' order.

    Equivalent to: scripts/38-family.sh (its first three stages)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FAMILY_PREPARE {
    label 'process_single'

    publishDir { "${params.outdir}/family" }, mode: params.publish_dir_mode

    input:
    path(pedigree, stageAs: 'pedigree.tsv')
    val(ids)             // the members' samplesheet ids, sorted
    val(sexes)           // each id's samplesheet sex ('' when the row gives none), in the same order
    path(gvcfs, stageAs: 'g?/*')

    output:
    path "members.tsv",  emit: members
    path "family.ped",   emit: ped
    path "trios.ped",    emit: trios
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gl = gvcfs instanceof List ? gvcfs : [gvcfs]
    def args = [ids, gl].transpose().collect { p -> "--gvcf ${p[0]}=${p[1]}" }.join(' ') + ' ' +
        [ids, sexes].transpose().findAll { p -> p[1] }.collect { p -> "--declared-sex ${p[0]}=${p[1]}" }.join(' ')
    """
    family_report.py prepare --pedigree pedigree.tsv ${args} --out-dir .

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'sample\\tsex\\tfather\\tmother\\tpartner\\tin_pedigree\\tsequenced\\n' > members.tsv
    touch family.ped trios.ped
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process GLNEXUS {
    label 'process_medium'

    input:
    val(ids)
    path(gvcfs, stageAs: 'g?/*')
    path(members)        // FAMILY_PREPARE's: the joint call waits for the pedigree check

    output:
    path "family.bcf",   emit: bcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gl = gvcfs instanceof List ? gvcfs : [gvcfs]
    def links = [ids, gl].transpose().collect { p -> "ln -s \"\$(readlink -f ${p[1]})\" gvcf/${p[0]}.g.vcf.gz" }.join('\n    ')
    def mem = Math.max(1, (task.memory.toGiga() as int) - 1)
    """
    mkdir gvcf
    ${links}
    # GLnexus refuses a database folder that exists: a fresh one each task
    glnexus_cli --config DeepVariantWGS --dir glnexus_db --threads ${task.cpus} --mem-gbytes ${mem} \\
        gvcf/*.g.vcf.gz > family.bcf
    rm -rf glnexus_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch family.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process GLNEXUS_VCF {
    label 'process_low'

    publishDir { "${params.outdir}/family" }, mode: params.publish_dir_mode

    input:
    path(bcf)

    output:
    tuple path("family.vcf.gz"), path("family.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view --threads ${task.cpus} -Oz -o family.vcf.gz ${bcf}
    bcftools index -t family.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > family.vcf.gz
    touch family.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
