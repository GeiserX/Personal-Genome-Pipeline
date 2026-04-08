/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MITO_VARIANTS — Mitochondrial variant calling with heteroplasmy detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uses GATK Mutect2 in mitochondrial mode to call variants on chrM,
    including low-frequency heteroplasmic variants (AF < 0.95).

    Four-step process:
    1. Extract chrM reads from BAM
    2. Ensure sequence dictionary exists
    3. Run Mutect2 --mitochondria-mode
    4. Filter variants with FilterMutectCalls --mitochondria-mode

    Equivalent to: scripts/20-mtoolbox.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MITO_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/${meta.id}/mito", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(reference_dict)

    output:
    tuple val(meta), path("*_chrM_filtered.vcf.gz"),     emit: mito_vcf
    tuple val(meta), path("*_chrM_filtered.vcf.gz.tbi"), emit: mito_vcf_index, optional: true
    tuple val(meta), path("*_chrM_mutect2.vcf.gz.stats"), emit: stats
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Step 1: Extract chrM reads
    gatk PrintReads \\
        -I ${bam} \\
        -L chrM \\
        -O ${prefix}_chrM.bam

    # Step 2: Run Mutect2 in mitochondrial mode
    gatk Mutect2 \\
        -R ${reference} \\
        -I ${prefix}_chrM.bam \\
        -L chrM \\
        --mitochondria-mode \\
        --max-mnp-distance 0 \\
        -O ${prefix}_chrM_mutect2.vcf.gz

    # Step 3: Filter variants
    gatk FilterMutectCalls \\
        -R ${reference} \\
        -V ${prefix}_chrM_mutect2.vcf.gz \\
        --mitochondria-mode \\
        -O ${prefix}_chrM_filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep 'GATK' | sed 's/.*v//' | sed 's/).*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chrM_filtered.vcf.gz
    touch ${prefix}_chrM_mutect2.vcf.gz.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: 4.6.1.0
    END_VERSIONS
    """
}
