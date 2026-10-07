/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGN_MINIMAP2 — FASTQ to a sorted, duplicate-marked, indexed BAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Three processes, the commands of scripts/02-alignment.sh:

    MINIMAP2_INDEX  minimap2 -x sr -d: the index with the preset the reads are
                    mapped with. main.nf uses <reference base>.sr.mmi beside
                    the FASTA when it is there (the file step 02 builds) or
                    --minimap2_index, and runs this only when neither exists.
    ALIGN_MINIMAP2  minimap2 -a -x sr with the read group step 02 writes
                    (ID, SM and LB the sample, PL ILLUMINA): GATK rejects reads
                    without one and callers take the sample name from SM.
    ALIGN_MARKDUP   samtools fixmate -m | sort | markdup, then index and
                    quickcheck, so PCR and optical duplicates carry the 0x400
                    flag.

    The script pipes minimap2 into samtools across two containers. A task
    runs in one image, and no image in versions.env holds both tools, so
    ALIGN_MINIMAP2 writes the SAM compressed (gzip -1) and ALIGN_MARKDUP reads
    it: the same commands and the same BAM, plus one temporary file in the
    work directory (docs/nextflow.md, Bash vs Nextflow parity).

    Equivalent to: scripts/02-alignment.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MINIMAP2_INDEX {
    tag "${reference.baseName}"
    label 'process_high'

    input:
    path(reference)

    output:
    path("${reference.baseName}.sr.mmi"), emit: index
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    minimap2 -x sr -t ${task.cpus} -d ${reference.baseName}.sr.mmi.tmp ${reference}
    mv ${reference.baseName}.sr.mmi.tmp ${reference.baseName}.sr.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${reference.baseName}.sr.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process ALIGN_MINIMAP2 {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path(index)

    output:
    tuple val(meta), path("${meta.id}.sam.gz"), emit: sam
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The read group as step 02 writes it: the script gets a literal \t,
    // which minimap2 turns into a tab.
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}"
    """
    minimap2 -t ${task.cpus} -a -x sr \\
        -R '${rg}' \\
        ${index} \\
        ${fastq_1} \\
        ${fastq_2} \\
        | gzip -1 > ${meta.id}.sam.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.sam.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process ALIGN_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/aligned" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}_sorted.bam"), path("${meta.id}_sorted.bam.bai"), emit: bam
    path "versions.yml",                                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // samtools sort -m is per thread: the label's memory covers cpus x 1 GB.
    """
    mkdir -p sort_tmp
    gzip -dc ${sam} \\
        | samtools fixmate -u -m - - \\
        | samtools sort -u -@ ${task.cpus} -m 1G -T sort_tmp/sort - \\
        | samtools markdup -@ ${task.cpus} -T sort_tmp/markdup - ${meta.id}_sorted.bam
    rm -rf sort_tmp
    samtools index -@ ${task.cpus} ${meta.id}_sorted.bam
    samtools quickcheck -v ${meta.id}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sorted.bam ${meta.id}_sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
