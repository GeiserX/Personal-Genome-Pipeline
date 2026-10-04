/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VERIFYBAMID2 — How much of the sample is someone else's DNA?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Estimates FREEMIX, the share of reads from another person, from allele
    balance at the 100,000 markers of the 1000 Genomes panel that
    --verifybamid2_panel holds (a folder with the .UD, .mu and .bed files
    scripts/setup.sh installs). VerifyBamID2 refuses to estimate when fewer
    than 1,000 markers have reads (a targeted or sliced BAM); the task then
    runs it again with --DisableSanityCheck and writes "skipped" to
    <id>.marker_check, so the report says on how many markers FREEMIX rests.
    SAMPLE_QC compares FREEMIX with --freemix_warn; contamination never stops
    the run.

    Equivalent to: the VerifyBamID2 part of scripts/33-sample-qc.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VERIFYBAMID2 {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/qc/verifybamid2" }, mode: params.publish_dir_mode,
        pattern: "*.{selfSM,Ancestry,log}"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(panel)

    output:
    tuple val(meta), path("${meta.id}.selfSM"), path("${meta.id}.marker_check"), emit: selfsm
    path "${meta.id}.log",                                         emit: log
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    UD=\$(find -L ${panel} -maxdepth 1 -name '*.UD' | head -n 1)
    if [ -z "\$UD" ]; then
        echo "ERROR: --verifybamid2_panel ${panel} holds no .UD file (scripts/setup.sh installs the panel)" >&2
        exit 1
    fi
    PREFIX=\${UD%.UD}
    vb2() {
        verifybamid2 --SVDPrefix "\$PREFIX" --Reference ${reference} --BamFile ${bam} \\
            --Output ${meta.id} --NumThread ${task.cpus} "\$@" > ${meta.id}.log 2>&1
    }
    MARKER_CHECK=passed
    if ! vb2; then
        if grep -q 'Insufficient Available markers' ${meta.id}.log; then
            echo "Fewer than 1,000 panel markers have reads: running again with --DisableSanityCheck"
            MARKER_CHECK=skipped
            vb2 --DisableSanityCheck || { tail -n 20 ${meta.id}.log >&2; exit 1; }
        else
            tail -n 20 ${meta.id}.log >&2
            exit 1
        fi
    fi
    echo "\$MARKER_CHECK" > ${meta.id}.marker_check
    cat ${meta.id}.selfSM

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    echo passed > ${meta.id}.marker_check
    printf '#SEQ_ID\\tRG\\tCHIP_ID\\t#SNPS\\t#READS\\tAVG_DP\\tFREEMIX\\n%s\\tNA\\tNA\\t100000\\t1\\t30\\t0.001\\n' ${meta.id} > ${meta.id}.selfSM
    touch ${meta.id}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
