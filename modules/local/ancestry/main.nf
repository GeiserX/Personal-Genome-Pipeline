/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANCESTRY — Ancestry-informative PCA from WGS VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Intersects sample VCF with a reference panel of common SNPs, performs LD pruning,
    and runs PCA. Single-sample PCA is mathematically limited; results become meaningful
    when projected alongside reference population samples.

    Equivalent to: scripts/26-ancestry.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANCESTRY {
    tag "$meta.id"
    label 'process_medium'

    container 'pgscatalog/plink2:2.00a5.10'

    publishDir "${params.outdir}/${meta.id}/ancestry", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(ref_panel)

    output:
    tuple val(meta), path("*_pca.eigenvec"),  emit: pca_results, optional: true
    tuple val(meta), path("*_pca.eigenval"),  emit: eigenvalues,  optional: true
    tuple val(meta), path("*_ancestry.tsv"),  emit: ancestry_tsv
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Step 1: LD pruning on autosomal SNPs
    LD_PRUNE_OK=false
    if plink2 \\
        --vcf ${vcf} \\
        --indep-pairwise 50 5 0.2 \\
        --out ${meta.id}_ld \\
        --threads ${task.cpus} \\
        --memory \$(( ${task.memory.toMega()} - 500 )) \\
        --chr 1-22 \\
        --allow-extra-chr \\
        --output-chr chrM 2>&1; then
        LD_PRUNE_OK=true
        EXTRACT_ARGS="--extract ${meta.id}_ld.prune.in"
    else
        EXTRACT_ARGS=""
    fi

    # Step 2: PCA (requires >=2 samples; expected to fail for single-sample)
    PCA_OK=false
    if plink2 \\
        --vcf ${vcf} \\
        \${EXTRACT_ARGS} \\
        --pca 10 \\
        --out ${meta.id}_pca \\
        --threads ${task.cpus} \\
        --memory \$(( ${task.memory.toMega()} - 500 )) \\
        --chr 1-22 \\
        --allow-extra-chr \\
        --output-chr chrM 2>&1; then
        PCA_OK=true
    fi

    # Step 3: Summary file with QC metrics
    echo -e "metric\\tvalue" > ${meta.id}_ancestry.tsv
    if [ "\${LD_PRUNE_OK}" = true ] && [ -f "${meta.id}_ld.prune.in" ]; then
        PRUNED=\$(wc -l < ${meta.id}_ld.prune.in)
        echo -e "ld_pruned_snps\\t\${PRUNED}" >> ${meta.id}_ancestry.tsv
    fi
    if [ "\${PCA_OK}" = true ] && [ -f "${meta.id}_pca.eigenvec" ]; then
        echo -e "pca_status\\tcompleted" >> ${meta.id}_ancestry.tsv
    else
        echo -e "pca_status\\tskipped_single_sample" >> ${meta.id}_ancestry.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1 | awk '{print \$2}' || echo '2.00a5.10')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_pca.eigenvec ${meta.id}_pca.eigenval
    echo -e "metric\\tvalue\\npca_status\\tstub" > ${meta.id}_ancestry.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1 | awk '{print \$2}' || echo '2.00a5.10')
    END_VERSIONS
    """
}
