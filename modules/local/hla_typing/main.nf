/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HLA_TYPING — HLA allele typing from WGS BAM using T1K
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Types HLA-A, B, C (Class I) and DRB1, DQB1, DPB1 (Class II)
    at 4-digit resolution using IPD-IMGT/HLA database against GRCh38.

    The index (allele sequences and their GRCh38 coordinates) comes from
    T1K_BUILD, built once per run for every sample; this process types one
    BAM against it.

    Equivalent to: step 2 of scripts/08-hla-typing.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HLA_TYPING {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/hla" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(seq_fa)    // T1K_BUILD: allele sequences
    path(coord_fa)  // T1K_BUILD: their GRCh38 coordinates

    output:
    tuple val(meta), path("*_hla_genotype.tsv"), emit: hla_alleles
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run-t1k \\
        -b ${bam} \\
        -f ${seq_fa} \\
        -c ${coord_fa} \\
        --preset hla-wgs \\
        -t ${task.cpus} \\
        --od ./ \\
        -o ${prefix}_hla

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_hla_genotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
