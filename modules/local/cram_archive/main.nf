/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CRAM_ARCHIVE, CRAM_TO_BAM — Alignments as CRAM, about half the size of a BAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CRAM_ARCHIVE writes <id>_sorted.cram (+ .crai) beside the published BAM
    and fails the task unless the CRAM passes `samtools quickcheck` and its
    `samtools flagstat` equals the BAM's line for line. It never deletes the
    BAM: delete it yourself once the CRAM is there, or use
    scripts/34-cram-archive.sh --delete-bam. A CRAM can only be read with the
    reference it was written with: keep that FASTA.

    CRAM_TO_BAM reads a samplesheet row given as cram,crai: it writes the BAM
    every BAM step of this pipeline reads, in the work directory only, and
    checks it against the CRAM the same way.

    Equivalent to: scripts/34-cram-archive.sh (and its --restore)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CRAM_ARCHIVE {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/aligned" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_sorted.cram"), path("${meta.id}_sorted.cram.crai"), emit: cram
    tuple val(meta), path("${meta.id}_sorted.cram.flagstat"),                            emit: flagstat
    path "versions.yml",                                                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -@ ${task.cpus} -C --reference ${reference} -o ${meta.id}_sorted.cram ${bam}
    samtools index -@ ${task.cpus} ${meta.id}_sorted.cram
    samtools quickcheck -v ${meta.id}_sorted.cram
    samtools flagstat -@ ${task.cpus} ${bam} > bam.flagstat
    samtools flagstat -@ ${task.cpus} --input-fmt-option reference=${reference} ${meta.id}_sorted.cram \\
        > ${meta.id}_sorted.cram.flagstat
    if ! diff bam.flagstat ${meta.id}_sorted.cram.flagstat >&2; then
        echo "ERROR: the CRAM of ${meta.id} does not hold the same reads as its BAM (flagstat above: < BAM, > CRAM)" >&2
        exit 1
    fi
    echo "${meta.id}: BAM \$(du -k ${bam} | cut -f1) kB, CRAM \$(du -k ${meta.id}_sorted.cram | cut -f1) kB, flagstat equal"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sorted.cram ${meta.id}_sorted.cram.crai ${meta.id}_sorted.cram.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process CRAM_TO_BAM {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_from_cram.bam"), path("${meta.id}_from_cram.bam.bai"), emit: bam
    path "versions.yml",                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -@ ${task.cpus} -b --reference ${reference} -o ${meta.id}_from_cram.bam ${cram}
    samtools index -@ ${task.cpus} ${meta.id}_from_cram.bam
    samtools quickcheck -v ${meta.id}_from_cram.bam
    samtools flagstat -@ ${task.cpus} --input-fmt-option reference=${reference} ${cram} > cram.flagstat
    samtools flagstat -@ ${task.cpus} ${meta.id}_from_cram.bam > bam.flagstat
    if ! diff cram.flagstat bam.flagstat >&2; then
        echo "ERROR: the BAM written from ${cram} does not hold the same reads (flagstat above: < CRAM, > BAM)" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_from_cram.bam ${meta.id}_from_cram.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
