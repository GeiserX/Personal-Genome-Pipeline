/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRAITS — The traits panel of step 42, from the gVCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRAITS_SITES  genotypes the positions of assets/traits.tsv from the gVCF:
                  each reference block over a position becomes a record with
                  the reference base, and a position without reads has no
                  line (it is "not covered", never the reference)
    TRAITS        <id>_traits.tsv: the reading of each genotype from the
                  table, and the APOE line (bin/traits_report.py)

    Runs for the samples with a gVCF only (DEEPVARIANT's, or one the
    samplesheet gives); main.nf logs the samples it skips.

    Equivalent to: scripts/42-traits.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TRAITS_SITES {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/traits" }, mode: params.publish_dir_mode, pattern: '*_traits_genotypes.tsv'

    input:
    tuple val(meta), path(gvcf), path(gvcf_index)
    path(table)
    path(reference)
    path(reference_fai)  // staged beside the FASTA

    output:
    tuple val(meta), path("${meta.id}_traits_genotypes.tsv"), emit: genotypes
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The same commands as scripts/42-traits.sh
    """
    grep -v '^#' ${table} | awk -F'\\t' 'NR > 1 { print \$1 "\\t" \$2 }' | LC_ALL=C sort -k1,1 -k2,2n -u > sites.tsv
    bcftools convert --gvcf2vcf -f ${reference} -R sites.tsv -Ou ${gvcf} \\
        | bcftools view -T sites.tsv --trim-alt-alleles -i 'GT!="mis"' -Ou \\
        | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%TGT]\\n' > ${meta.id}_traits_genotypes.tsv
    rm -f sites.tsv

    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_traits_genotypes.tsv
    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process TRAITS {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/traits" }, mode: params.publish_dir_mode, pattern: '*_traits.tsv'

    input:
    tuple val(meta), path(genotypes)
    path(table)

    output:
    tuple val(meta), path("${meta.id}_traits.tsv"), emit: traits
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    traits_report.py --table ${table} --genotypes ${genotypes} --sample ${meta.id} --out ${meta.id}_traits.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'rsid\\tgene\\tchrom\\tpos\\ttrait\\tgenotype\\tstatus\\treading\\tancestry\\tnote\\tpmids\\n' > ${meta.id}_traits.tsv
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
