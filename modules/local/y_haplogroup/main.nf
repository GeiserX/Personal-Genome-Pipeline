/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Y_HAPLOGROUP — Y-chromosome haplogroup (paternal line) with Yleaf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('y_haplogroup' in --tools). Runs on male samples only:
    workflows/bam_analysis.nf passes the BAMs whose sex is male, the sex
    INDEXCOV has checked against the reads.

    The Yleaf image has no samtools, which Yleaf calls for a BAM, and no
    marker tables (--yleaf_data, setup.sh --yleaf-data), so the work is
    three processes, as in scripts/37-y-haplogroup.sh (bin/yleaf_run.py
    says why):
      Y_POSITIONS   Yleaf's GRCh38 marker positions (Yleaf image, once)
      Y_PILEUP      samtools idxstats and the pileup at those positions
      Y_HAPLOGROUP  Yleaf on that pileup, offline; Hg NA means too few Y
                    markers had reads to call a haplogroup

    Equivalent to: scripts/37-y-haplogroup.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process Y_POSITIONS {
    label 'process_single'

    input:
    path(yleaf_data)   // --yleaf_data: reference/yleaf-<release>/data (setup.sh --yleaf-data)

    output:
    path "y_positions.txt", emit: positions
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    yleaf_run.py positions --data ${yleaf_data} y_positions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'chrY\\t2781480\\n' > y_positions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process Y_PILEUP {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)
    path(positions)

    output:
    tuple val(meta), path(bam), path(bai), path("${meta.id}_idxstats.txt"), path("${meta.id}_y_pileup.txt"), emit: pileup
    path "versions.yml",                                                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // -AQ20q1: the flags Yleaf gives samtools mpileup at its default quality 20
    """
    samtools idxstats ${bam} > ${meta.id}_idxstats.txt
    samtools mpileup -l ${positions} -AQ20q1 ${bam} > ${meta.id}_y_pileup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_idxstats.txt ${meta.id}_y_pileup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process Y_HAPLOGROUP {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/y_haplogroup" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(idxstats), path(pileup)
    path(reference)
    path(yleaf_data)

    output:
    tuple val(meta), path("${meta.id}_y_haplogroup.txt"), emit: haplogroup
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    yleaf_run.py predict --data ${yleaf_data} --bam ${bam} --reference "\$(readlink -f ${reference})" \\
        --idxstats ${idxstats} --pileup ${pileup} --out yleaf
    cp yleaf/hg_prediction.hg ${meta.id}_y_haplogroup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'Sample_name\\tHg\\tHg_marker\\tTotal_reads\\tValid_markers\\tQC-score\\tQC-1\\tQC-2\\tQC-3\\n${meta.id}\\tNA\\t\\t0\\t0\\t0\\t0\\t0\\t0\\n' > ${meta.id}_y_haplogroup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
