/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CYP2D6_DEPTH and PGX_CONSENSUS — what the BAM-based callers tell PharmCAT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CYP2D6_DEPTH runs mosdepth over CYP2D6 and two flanks, once for all reads
    and once for MAPQ >= 1, before pypgx or Cyrius calls copy number from
    depth. PYPGX and CYRIUS judge the two files with
    bin/cyp2d6_depth_check.py, which also refuses regions other than its own,
    so the copy of them below cannot drift unnoticed.

    PGX_CONSENSUS runs bin/pgx_outside_calls.py: PharmCAT's outside-call file
    (HLA-A and HLA-B from T1K; CYP2D6 only when pypgx and Cyrius agree and the
    depth check passed) and the consensus table that says what each caller
    said. PHARMCAT reads the file with -po.

    Equivalent to: the depth check of scripts/32-pypgx.sh and
    scripts/21-cyrius.sh, and scripts/36-pgx-consensus.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CYP2D6_DEPTH {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.cyp2d6.q0.regions.bed.gz"), path("${meta.id}.cyp2d6.q1.regions.bed.gz"), emit: regions
    path "versions.yml",                                                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The regions of `cyp2d6_depth_check.py bed` (GRCh38, 0-based).
    """
    printf 'chr22\\t42050000\\t42100000\\tflank\\nchr22\\t42123192\\t42132032\\tCYP2D6\\nchr22\\t42200000\\t42250000\\tflank\\n' \\
        > cyp2d6_regions.bed
    for q in 0 1; do
        mosdepth -n -c chr22 -t ${task.cpus} -Q "\$q" -b cyp2d6_regions.bed "${meta.id}.cyp2d6.q\$q" ${bam}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'chr22\\t42050000\\t42100000\\tflank\\t30\\nchr22\\t42123192\\t42132032\\tCYP2D6\\t20\\nchr22\\t42200000\\t42250000\\tflank\\t30\\n' \\
        | gzip -c > ${meta.id}.cyp2d6.q0.regions.bed.gz
    cp ${meta.id}.cyp2d6.q0.regions.bed.gz ${meta.id}.cyp2d6.q1.regions.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process PGX_CONSENSUS {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/pgx_consensus" }, mode: params.publish_dir_mode

    input:
    // each file is [] when its step did not run for the sample
    tuple val(meta), path(hla_genotype), path(pypgx_summary), path(cyrius_tsv), path(depth_check)

    output:
    tuple val(meta), path("${meta.id}_outside_calls.tsv"), emit: calls
    tuple val(meta), path("${meta.id}_pgx_consensus.tsv"), emit: consensus
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = []
    if (hla_genotype)  { args << "--hla ${hla_genotype}" }
    if (pypgx_summary) { args << "--pypgx ${pypgx_summary}" }
    if (cyrius_tsv)    { args << "--cyrius ${cyrius_tsv}" }
    if (depth_check)   { args << "--depth-check ${depth_check}" }
    """
    pgx_outside_calls.py \\
        --calls ${meta.id}_outside_calls.tsv \\
        --consensus ${meta.id}_pgx_consensus.tsv \\
        ${args.join(' ')}

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_outside_calls.tsv
    printf 'Gene\\tResult\\tOutside_call\\tReason\\tEvidence\\n' > ${meta.id}_pgx_consensus.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
