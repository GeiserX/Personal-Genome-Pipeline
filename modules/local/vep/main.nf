/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VEP — Ensembl Variant Effect Predictor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Full functional annotation: consequence, SIFT, PolyPhen, gnomAD AF, ClinVar,
    regulatory features, etc. Uses offline cache for reproducibility.

    Equivalent to: scripts/13-vep-annotation.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VEP {
    tag "$meta.id"
    label 'process_high'

    container 'ensemblorg/ensembl-vep:release_112.0'

    publishDir "${params.outdir}/${meta.id}/vep", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(reference)
    path(vep_cache)

    output:
    tuple val(meta), path("*_vep.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_vep.vcf.gz.tbi"), emit: vcf_index
    tuple val(meta), path("*_vep.vcf_summary.html"), emit: stats, optional: true
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${meta.id}_vep.vcf \\
        --vcf \\
        --cache \\
        --dir_cache ${vep_cache} \\
        --offline \\
        --assembly GRCh38 \\
        --everything \\
        --af_gnomade \\
        --force_overwrite \\
        --fork ${task.cpus} \\
        --fasta ${reference} \\
        ${args}

    bgzip -c ${meta.id}_vep.vcf > ${meta.id}_vep.vcf.gz
    tabix -p vcf ${meta.id}_vep.vcf.gz
    rm -f ${meta.id}_vep.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$(vep --help 2>&1 | grep 'ensembl-vep' | sed 's/.*: //' || echo 'release_112.0')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_vep.vcf.gz
    touch ${meta.id}_vep.vcf.gz.tbi
    touch ${meta.id}_vep.vcf_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: release_112.0
    END_VERSIONS
    """
}
