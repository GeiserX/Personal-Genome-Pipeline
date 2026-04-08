/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MANTA — Structural variant calling (deletions, duplications, inversions, translocations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two-step process:
    1. Configure Manta (configManta.py)
    2. Run workflow (runWorkflow.py)

    Equivalent to: scripts/04-manta.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MANTA {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_2'

    publishDir "${params.outdir}/${meta.id}/manta", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("results/variants/diploidSV.vcf.gz"),    emit: diploid_sv
    tuple val(meta), path("results/variants/candidateSV.vcf.gz"),  emit: candidate_sv
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Step 1: Configure Manta
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${reference} \\
        --runDir manta_run

    # Step 2: Run Manta workflow
    manta_run/runWorkflow.py -j ${task.cpus}

    # Move results to expected output location
    mv manta_run/results .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version 2>&1 | sed 's/.*version //' || echo '1.6.0')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p results/variants
    touch results/variants/diploidSV.vcf.gz
    touch results/variants/candidateSV.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: 1.6.0
    END_VERSIONS
    """
}
