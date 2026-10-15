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
                    Then the read count, the rule of step 02
                    (scripts/lib/common.sh, read_count_check): the BAM's
                    primary records (samtools view -c -F 0x900) must equal
                    the FASTQ records of R1 plus R2, counted while minimap2
                    reads them, because minimap2 writes every read once,
                    mapped or not. A difference fails the task (exit 1, not
                    retried). Both counts go to qc/<id>_read_count.tsv.
                    ext.read_count_check = 'warn' (a -c config: withName
                    'ALIGN_MINIMAP2') keeps the BAM and the table says
                    mismatch, as READ_COUNT_CHECK=warn does for step 02.

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

    publishDir { "${params.outdir}/${meta.id}/aligned" }, mode: params.publish_dir_mode,
        saveAs: { f -> f.endsWith('_read_count.tsv') ? null : f }
    publishDir { "${params.outdir}/${meta.id}/qc" }, mode: params.publish_dir_mode, pattern: '*_read_count.tsv'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path(index)

    output:
    tuple val(meta), path("${meta.id}_sorted.bam"), path("${meta.id}_sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_read_count.tsv"),                               emit: read_count
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
    def read_count_check = task.ext.read_count_check ?: 'fail'
    if (!(read_count_check in ['fail', 'warn'])) {
        error "ALIGN_MINIMAP2: ext.read_count_check must be 'fail' or 'warn', got '${read_count_check}'"
    }
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
    //
    // The read counters run beside minimap2, one gzip per FASTQ, each in a
    // subshell with pipefail so that wait returns a gzip error too. awk
    // prints the records (lines / 4) or exits 3 when the lines are not a
    // multiple of 4, as FASTQ_RECORDS_AWK in scripts/lib/common.sh does.
    def count_awk = 'END { if (NR % 4) { printf "%.0f lines, not a multiple of 4: a record is cut short\\n", NR > "/dev/stderr"; exit 3 } printf "%.0f\\n", NR / 4 }'
    """
    mkdir -p sort_tmp
    ( set -o pipefail; gzip -cdf ${fastq_1} | awk '${count_awk}' > r1.count ) &
    r1_pid=\$!
    ( set -o pipefail; gzip -cdf ${fastq_2} | awk '${count_awk}' > r2.count ) &
    r2_pid=\$!
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

    wait \$r1_pid && wait \$r2_pid \\
        || { echo "ERROR: could not count the records of ${fastq_1} and ${fastq_2} (see the message above)" >&2; exit 1; }
    n1=\$(cat r1.count)
    n2=\$(cat r2.count)
    reads_in=\$((n1 + n2))
    primary=\$(samtools view -c -F 0x900 -@ ${task.cpus} ${meta.id}_sorted.bam)
    missing=\$((reads_in - primary))
    src="${fastq_1} \$n1 records + ${fastq_2} \$n2 records"
    verdict=ok
    [ "\$missing" -eq 0 ] || verdict=mismatch
    printf 'key\\tvalue\\ninput_reads\\t%s\\ninput_source\\t%s\\nprimary_records\\t%s\\nmissing_reads\\t%s\\nread_count_check\\t%s\\n' \\
        "\$reads_in" "\$src" "\$primary" "\$missing" "\$verdict" > ${meta.id}_read_count.tsv
    if [ "\$verdict" = ok ]; then
        echo "Read count: the BAM holds \$primary primary records, the \$reads_in reads that went in (\$src)."
    else
        level=ERROR
        [ "${read_count_check}" = fail ] || level=WARNING
        if [ "\$missing" -gt 0 ]; then
            pct=\$(awk -v m="\$missing" -v n="\$reads_in" 'BEGIN { printf "%.1f", 100 * m / n }')
            echo "\$level: the BAM holds \$primary primary records, but \$reads_in reads went in (\$src): \$missing reads (\$pct%) are missing." >&2
        else
            echo "\$level: the BAM holds \$primary primary records, more than the \$reads_in reads that went in (\$src): it was not made from these reads." >&2
        fi
        echo "  The aligner writes every read it is given, mapped or not, once, so the two counts must be equal." >&2
        echo "  Depth, copy number (CYP2D6), CNV, telomere content and every call that depends on coverage would be off." >&2
        echo "  An R1 and R2 with different numbers of records is one cause: minimap2 and bwa-mem2 stop at the end" >&2
        echo "  of the shorter file with only a warning. See docs/troubleshooting.md, \\"The BAM lost reads\\"." >&2
        if [ "${read_count_check}" = warn ]; then
            echo "  ext.read_count_check = 'warn': going on." >&2
        else
            echo "  For an input you know is odd, ext.read_count_check = 'warn' on ALIGN_MINIMAP2 goes on anyway." >&2
            exit 1
        fi
    fi

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
    printf 'key\\tvalue\\n' > ${meta.id}_read_count.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2-samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
