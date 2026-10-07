/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTP — Read QC and adapter trimming before alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Adapter detection by pair overlap, Q20 sliding-window trimming from both
    ends, polyG tail removal and a 36 bp length floor, with the flags of
    scripts/01b-fastp-qc.sh. The trimmed reads go to fastq_trimmed/, under the
    names the script writes, so they never share a name with the staged input
    (fastp would write through the input's symlink).

    Skipped with --skip_trim (the script's SKIP_TRIM=true).

    Equivalent to: scripts/01b-fastp-qc.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTP {
    tag "$meta.id"
    label 'process_medium'

    // The reports only: the trimmed reads are as large as the input and go
    // straight to alignment, so they stay in the work directory.
    publishDir { "${params.outdir}/${meta.id}" }, mode: params.publish_dir_mode, pattern: 'fastq_trimmed/*_fastp.*'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("fastq_trimmed/${meta.id}_R1.fastq.gz"), path("fastq_trimmed/${meta.id}_R2.fastq.gz"), emit: reads
    tuple val(meta), path("fastq_trimmed/${meta.id}_fastp.json"),  emit: json
    tuple val(meta), path("fastq_trimmed/${meta.id}_fastp.html"),  emit: html
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p fastq_trimmed
    fastp \\
        -i ${fastq_1} \\
        -I ${fastq_2} \\
        -o fastq_trimmed/${meta.id}_R1.fastq.gz \\
        -O fastq_trimmed/${meta.id}_R2.fastq.gz \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 \\
        --cut_front \\
        --cut_tail \\
        --cut_mean_quality 20 \\
        --length_required 36 \\
        -g \\
        -R ${meta.id} \\
        -j fastq_trimmed/${meta.id}_fastp.json \\
        -h fastq_trimmed/${meta.id}_fastp.html \\
        -w ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p fastq_trimmed
    echo | gzip > fastq_trimmed/${meta.id}_R1.fastq.gz
    echo | gzip > fastq_trimmed/${meta.id}_R2.fastq.gz
    touch fastq_trimmed/${meta.id}_fastp.json fastq_trimmed/${meta.id}_fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
