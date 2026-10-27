/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LPA — the LPA KIV-2 repeat count estimated from read depth (Lp(a))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('lpa' in --tools). Three processes, as in scripts/46-lpa.sh:
      LPA_REGIONS  bin/lpa_kiv2.py bed, once: the KIV-2 and single-copy LPA
                   exons (the regions live only in that script)
      LPA_DEPTH    samtools bedcov -j over them, all MAPQ and MAPQ >= 1, for
                   each sample with a BAM; empty tables when the BAM has no chr6
      LPA          bin/lpa_kiv2.py estimate: the count, its 95% interval and
                   the message

    Equivalent to: scripts/46-lpa.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process LPA_REGIONS {
    label 'process_single'

    output:
    path "lpa_regions.bed", emit: bed
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    lpa_kiv2.py bed > lpa_regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch lpa_regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process LPA_DEPTH {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/lpa" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(regions)   // LPA_REGIONS' BED

    output:
    tuple val(meta), path("${meta.id}_lpa_kiv2_depth.tsv"), path("${meta.id}_lpa_kiv2_depth_q1.tsv"), emit: depth
    path "versions.yml",                                                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // bedcov stops on a contig the BAM has no @SQ line for: without chr6 both
    // tables stay empty. The same commands as scripts/46-lpa.sh.
    """
    if samtools view -H ${bam} | awk -F '\\t' '\$1 == "@SQ" {for (i = 2; i <= NF; i++) if (\$i == "SN:chr6") f = 1} END {exit !f}'; then
        samtools bedcov -j ${regions} ${bam} > ${meta.id}_lpa_kiv2_depth.tsv
        samtools bedcov -j -Q 1 ${regions} ${bam} > ${meta.id}_lpa_kiv2_depth_q1.tsv
    else
        : > ${meta.id}_lpa_kiv2_depth.tsv
        : > ${meta.id}_lpa_kiv2_depth_q1.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_lpa_kiv2_depth.tsv ${meta.id}_lpa_kiv2_depth_q1.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process LPA {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/lpa" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(depth), path(depth_q1)

    output:
    tuple val(meta), path("${meta.id}_lpa_kiv2.tsv"), emit: table
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    lpa_kiv2.py estimate --bedcov ${depth} --bedcov-mapq1 ${depth_q1} --out ${meta.id}_lpa_kiv2.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\nmessage\\tstub\\n' > ${meta.id}_lpa_kiv2.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
