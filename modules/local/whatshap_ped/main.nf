/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WHATSHAP_PED, FAMILY_REPORT — Who passed on which allele
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WHATSHAP_PED phases the family's joint VCF with whatshap phase --ped over
    the trios FAMILY_PREPARE found (children whose two parents are
    sequenced), by inheritance alone: no reads are given. A site where the
    child and both parents are heterozygous stays unphased, as no genotype
    tells it. With no trio it writes nothing.

    FAMILY_REPORT runs bin/family_report.py report (the rule
    scripts/38-family.sh runs too) on the phased VCF, or on the joint VCF
    when there is no trio: per child with a sequenced parent, the allele from
    each parent at every site where the genotypes tell; de novo candidates
    when both parents are sequenced; compound heterozygotes in trans or in
    cis among the child's clinical-filter candidates (step 23's VCF, when
    clinical_filter ran); and family_report.md.

    Equivalent to: scripts/38-family.sh (its last two stages)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process WHATSHAP_PED {
    label 'process_medium'

    publishDir { "${params.outdir}/family" }, mode: params.publish_dir_mode

    input:
    tuple path(vcf), path(vcf_index)
    path(trios)
    path(relationships)  // FAMILY_CHECK's: phasing waits for the relationship check

    output:
    tuple path("family.phased.vcf.gz"), path("family.phased.vcf.gz.tbi"), emit: vcf, optional: true
    path "versions.yml",                                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [ -s ${trios} ]; then
        whatshap phase --ped ${trios} --use-ped-samples -o family.phased.vcf.gz ${vcf}
        python -c "import pysam; pysam.tabix_index('family.phased.vcf.gz', preset='vcf', force=True)"
    else
        echo "No child has both parents sequenced: nothing to phase by inheritance."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    if [ -s ${trios} ]; then
        echo | gzip > family.phased.vcf.gz
        touch family.phased.vcf.gz.tbi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process FAMILY_REPORT {
    label 'process_low'

    publishDir { "${params.outdir}/family" }, mode: params.publish_dir_mode

    input:
    path(members)
    path(relationships)
    tuple path(vcf), path(vcf_index)
    path(phased)         // WHATSHAP_PED's VCF, or [] with no trio
    val(candidate_ids)
    path(candidates, stageAs: 'c?/*')
    path(par_bed)

    output:
    path "*_inheritance.tsv.gz", emit: inheritance, optional: true
    path "de_novo.tsv",          emit: de_novo
    path "compound_het.tsv",     emit: compound_het
    path "family_report.md",     emit: report
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cl = candidates instanceof List ? candidates : (candidates ? [candidates] : [])
    def cand = [candidate_ids, cl].transpose().collect { p -> "--candidates ${p[0]}=${p[1]}" }.join(' ')
    def ph = phased ? "--phased ${phased}" : ''
    """
    family_report.py report --members ${members} --relationships ${relationships} \\
        --vcf ${vcf} ${ph} ${cand} --par-bed ${par_bed} --out-dir .

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch de_novo.tsv compound_het.tsv family_report.md
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
