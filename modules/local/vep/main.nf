/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VEP — Ensembl Variant Effect Predictor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Full functional annotation: consequence, SIFT, PolyPhen, gnomAD AF, ClinVar,
    regulatory features, etc. Uses offline cache for reproducibility. With
    --clinvar, the same ClinVar file the screen reads is added with --custom:
    its CLNSIG, CLNREVSTAT and CLNDN as ClinVar_CLNSIG, ClinVar_CLNREVSTAT and
    ClinVar_CLNDN, so the clinical filter's ClinVar tier follows a refresh of
    that file instead of the cache release.

    VEP reads the PASS records only ('.' too, for a caller that writes no
    FILTER): the clinical filter and slivar, after it, keep PASS records alone. The VEP image has
    no bcftools, so awk selects them on the FILTER column, the records
    `bcftools view -f PASS,.` keeps in scripts/13-vep-annotation.sh.

    Equivalent to: scripts/13-vep-annotation.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VEP {
    tag "$meta.id"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/vep" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(reference)
    path(reference_fai)  // staged beside the FASTA, so no task builds its own
    path(vep_cache)
    path(clinvar)        // [] without --clinvar
    path(clinvar_index)  // staged beside it

    output:
    tuple val(meta), path("*_vep.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_vep.vcf.gz.tbi"), emit: vcf_index
    tuple val(meta), path("*_vep.vcf_summary.html"), emit: stats, optional: true
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom = clinvar ? "--custom file=${clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN" : ''
    """
    bgzip -dc ${vcf} | awk -F '\\t' '/^#/ || \$7 == "." || (";" \$7 ";") ~ /;PASS;/' | bgzip -c > ${meta.id}.pass.vcf.gz

    vep \\
        --input_file ${meta.id}.pass.vcf.gz \\
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
        ${custom} \\
        ${args}

    bgzip -c ${meta.id}_vep.vcf > ${meta.id}_vep.vcf.gz
    tabix -p vcf ${meta.id}_vep.vcf.gz
    rm -f ${meta.id}_vep.vcf ${meta.id}.pass.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_vep.vcf.gz
    touch ${meta.id}_vep.vcf.gz.tbi
    touch ${meta.id}_vep.vcf_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
