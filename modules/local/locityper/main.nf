/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCITYPER — A second opinion at hard polymorphic genes (step 53)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('locityper' in --tools, with --locityper_data: setup.sh
    --locityper-data), for each sample with a BAM. One task, as
    scripts/53-locityper.sh: locityper preproc measures the sample's read
    depth, insert size and error profile on its background region
    (--locityper_bg_region, or Locityper's own for GRCh38), locityper
    genotype picks two pangenome haplotypes at each locus of the database,
    and bin/locityper_report.py, with the image's Python, writes one row per
    locus. Locityper's exit code goes to the table; a locus it did not finish
    reads "not called".

    Equivalent to: scripts/53-locityper.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process LOCITYPER {
    tag "$meta.id"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/locityper" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(data)      // reference/locityper-<version>/: counts.jf, db/, loci.bed

    output:
    tuple val(meta), path("${meta.id}_locityper.tsv"), emit: table
    tuple val(meta), path("${meta.id}_locityper.log"), path("work"), emit: run
    path "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bg = params.locityper_bg_region ? "-b ${params.locityper_bg_region}" : ''
    """
    mkdir -p work
    rc=0
    { locityper preproc -a ${bam} -r ${fasta} -j ${data}/counts.jf -o work/preproc -@ ${task.cpus} ${bg} \\
        && locityper genotype -a ${bam} -d ${data}/db -p work/preproc -o work/gt -@ ${task.cpus}; } \\
        > ${meta.id}_locityper.log 2>&1 || rc=\$?
    [ "\$rc" -eq 0 ] || { echo "WARNING: Locityper exited \$rc; last lines of its log:" >&2; tail -n 20 ${meta.id}_locityper.log >&2; }
    locityper_report.py --gt-dir work/gt --loci ${data}/loci.bed --exit-code \$rc --out ${meta.id}_locityper.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        locityper: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p work
    touch ${meta.id}_locityper.log
    printf 'locus\\tregion\\tstatus\\tgenotype\\thaplotype_1\\thaplotype_2\\tquality\\ttotal_reads\\tunexplained_reads\\tweighted_distance\\twarnings\\tfilter\\treason\\n' > ${meta.id}_locityper.tsv
    printf 'CYP2D6\\tchr22:42122285-42152960\\tnot called\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\tstub\\n' >> ${meta.id}_locityper.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        locityper: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
