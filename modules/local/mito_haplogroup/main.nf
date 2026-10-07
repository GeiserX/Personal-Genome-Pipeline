/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mitochondrial Haplogroup — Determine maternal lineage from mtDNA variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two-step process:
    1. Extract chrM variants from full-genome VCF (bcftools)
    2. Classify haplogroup with haplogrep3

    Equivalent to: scripts/12-mito-haplogroup.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MITO_EXTRACT_CHRM {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*_chrM.vcf.gz"), path("*_chrM.vcf.gz.tbi"), emit: chrm_vcf
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view -r chrM ${vcf} -Oz -o ${meta.id}_chrM.vcf.gz
    bcftools index -t ${meta.id}_chrM.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_chrM.vcf.gz
    touch ${meta.id}_chrM.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process MITO_HAPLOGROUP {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/mito" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(chrm_vcf), path(chrm_vcf_index)

    output:
    tuple val(meta), path("*_haplogroup.txt"), emit: haplogroup
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # haplogrep3 reads haplogrep3.yaml and its trees from the working
    # directory. In the task directory it finds neither and fetches the
    # config from GitHub, which fails without network (and was an unpinned
    # download before). So it runs in the directory of its binary, where the
    # image keeps them, as tests/smoke/haplogrep3.sh does.
    WD=\$PWD
    cd "\$(dirname "\$(readlink -f "\$(command -v haplogrep3)")")"
    haplogrep3 classify \\
        --tree phylotree-fu-rcrs@1.2 \\
        --input "\$WD/${chrm_vcf}" \\
        --output "\$WD/${meta.id}_haplogroup.txt" \\
        --extend-report
    REPORTED=\$( { haplogrep3 --version 2>&1 || true; } | awk '!v && match(\$0, /[0-9]+\\.[0-9]+(\\.[0-9]+)*/) { v = substr(\$0, RSTART, RLENGTH) } END { print (v != "" ? v : "unknown") }')
    cd "\$WD"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep3: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
        haplogrep3_reported: \${REPORTED}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_haplogroup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep3: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
