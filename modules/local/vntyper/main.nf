/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VNTYPER — MUC1 VNTR frameshifts (ADTKD-MUC1) with VNtyper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('vntyper' in --tools), for each sample with a BAM. One task, as
    scripts/52-vntyper.sh: VNtyper runs on the BAM (fast mode, GRCh38, from
    /opt/vntyper, where its reference paths point), then bin/vntyper_report.py,
    with the image's Python, writes the table: positive, possible, negative or
    not assessed, with the sensitivity at the sample's depth over the VNTR.
    VNtyper's exit code goes to the table; a run that stops reads "not
    assessed". The image's entrypoint takes only its own command words, so
    nextflow.config's docker profile replaces it for this process.

    Equivalent to: scripts/52-vntyper.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VNTYPER {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/vntyper" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_vntyper.tsv"), emit: table
    tuple val(meta), path("${meta.id}_vntyper.log"), path("run"), emit: run
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # VNtyper refuses to write its log inside the folder of its input BAM (both
    # the link's folder and the file's), and Nextflow stages the BAM in the task
    # folder, where run/ goes: the BAM is linked from a folder of its own.
    W=\$PWD
    mkdir -p in
    ln -s "\$W/${bam}" in/${bam}
    ln -s "\$W/${bai}" in/${bai}
    rc=0
    ( cd /opt/vntyper && vntyper pipeline --bam "\$W/in/${bam}" -o "\$W/run" --sample-name ${meta.id} --threads ${task.cpus} \\
        --fast-mode --reference-assembly hg38 ) > ${meta.id}_vntyper.log 2>&1 || rc=\$?
    [ "\$rc" -eq 0 ] || { echo "WARNING: VNtyper exited \$rc; last lines of its log:" >&2; tail -n 20 ${meta.id}_vntyper.log >&2; }
    mkdir -p run
    vntyper_report.py --run-dir run --exit-code \$rc --log ${meta.id}_vntyper.log --out ${meta.id}_vntyper.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vntyper: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p run
    touch ${meta.id}_vntyper.log
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\nmessage\\tstub\\n' > ${meta.id}_vntyper.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vntyper: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
