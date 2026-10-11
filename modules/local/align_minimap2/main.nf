/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGN_MINIMAP2 — FASTQ to a sorted, duplicate-marked, indexed BAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two processes, the commands of scripts/02-alignment.sh:

    MINIMAP2_INDEX  minimap2 -x sr -d: the index with the preset the reads are
                    mapped with. main.nf uses <reference base>.sr.mmi beside
                    the FASTA when it is there (the file step 02 builds) or
                    --minimap2_index, and runs this only when neither exists.
    ALIGN_MINIMAP2  minimap2 -a -x sr with the read group step 02 writes
                    (ID, SM and LB the sample, PL ILLUMINA: GATK rejects reads
                    without one and callers take the sample name from SM),
                    piped into samtools fixmate -m | sort | markdup, then
                    quickcheck, so PCR and optical duplicates carry the 0x400
                    flag. markdup writes the .bai while it writes the BAM
                    (--write-index), so no second pass reads the BAM.

    The script pipes minimap2 into samtools across two containers. A task
    runs in one image, so ALIGN_MINIMAP2 runs in MINIMAP2_SAMTOOLS_IMAGE, a
    BioContainers image that holds both, at the releases of MINIMAP2_IMAGE
    and SAMTOOLS_IMAGE (scripts/ci/check-images.sh fails when they part).
    No SAM reaches the disk, and the task writes the minimap2 and samtools
    versions it ran into versions.yml.

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

    publishDir { "${params.outdir}/${meta.id}/aligned" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path(index)

    output:
    tuple val(meta), path("${meta.id}_sorted.bam"), path("${meta.id}_sorted.bam.bai"), emit: bam
    path "versions.yml",                                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The read group as step 02 writes it: the script gets a literal \t,
    // which minimap2 turns into a tab.
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}"
    // The memory is the process_high label's (conf/base.config: 32 GB, twice
    // that on the retry after an out-of-memory exit, capped by --max_memory).
    // minimap2 and the samtools stages run at the same time in it, where step
    // 02 gives them two containers (32 GB, and THREADS + 4 GB for samtools:
    // 12 GB at 8 threads). samtools sort -m is per thread, so sort gets a
    // quarter of the task's memory: 1 GB per thread at 32 GB and 8 CPUs, the
    // -m 1G of step 02. The rest is minimap2's (about
    // 20 GB planned on GRCh38, 10 GB measured on the 1.8 Gb test reference,
    // docs/02-alignment.md) and markdup's, so the first attempt runs near
    // its 32 GB.
    def mem_mb  = task.memory ? task.memory.toMega() : 8192
    def sort_mb = Math.max(256, (mem_mb / 4 / task.cpus) as int)
    // The ##idx##<name>.bai suffix makes markdup's index a .bai, not a .csi.
    // markdup finishes the index just before the BAM's last block; the touch
    // makes the index the newer file, as step 02 does. test -s fails the task
    // if markdup wrote no index, so touch never creates an empty one.
    //
    // With pipefail the pipe's status is that of its last failing stage, so
    // a kill for memory hides behind the stage after it: minimap2 killed
    // mid-record exits 137, fixmate then exits 1 on the cut record, and exit
    // 1 is not retried (conf/base.config retries 137 and the other signal
    // exits). The || block exits with the last signal status in the pipe
    // (above 128), else 1, so the retry with twice the memory still happens.
    // The last one, not the first: a killed markdup leaves the stages before
    // it to die of SIGPIPE (141, not retried). A retry runs in a new work
    // directory, so the block also removes the failed attempt's sort spills.
    """
    mkdir -p sort_tmp
    minimap2 -t ${task.cpus} -a -x sr \\
        -R '${rg}' \\
        ${index} \\
        ${fastq_1} \\
        ${fastq_2} \\
        | samtools fixmate -@ ${task.cpus} -u -m - - \\
        | samtools sort -u -@ ${task.cpus} -m ${sort_mb}M -T sort_tmp/sort - \\
        | samtools markdup -@ ${task.cpus} --write-index -T sort_tmp/markdup - ${meta.id}_sorted.bam##idx##${meta.id}_sorted.bam.bai \\
        || {
            pipe_status=("\${PIPESTATUS[@]}")
            rc=1
            for s in "\${pipe_status[@]}"; do
                if [ "\$s" -gt 128 ]; then rc=\$s; fi
            done
            rm -rf sort_tmp
            echo "ERROR: the alignment pipe failed, exit statuses \${pipe_status[*]}" >&2
            exit "\$rc"
        }
    rm -rf sort_tmp
    samtools quickcheck -v ${meta.id}_sorted.bam
    test -s ${meta.id}_sorted.bam.bai
    touch ${meta.id}_sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2-samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sorted.bam ${meta.id}_sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2-samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
