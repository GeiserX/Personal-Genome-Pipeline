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
    def has_ref = ref_panel ? true : false
    """
    # Step 1: Intersect with reference panel if provided
    # Uses plink2 variant ID matching (equivalent to bcftools isec in bash script)
    EXTRACT_REF=""
    SHARED_COUNT=0
    if [ "${has_ref}" = "true" ] && [ -f "${ref_panel}" ]; then
        echo "Extracting variant IDs from reference panel..."
        plink2 \\
            --vcf ${ref_panel} \\
            --set-all-var-ids '@:#' \\
            --output-chr chrM \\
            --write-snplist \\
            --out ref_vars \\
            --chr 1-22 \\
            --allow-extra-chr \\
            --threads ${task.cpus} \\
            --memory \$(( ${task.memory.toMega()} / 4 )) 2>&1 || true

        if [ -f "ref_vars.snplist" ] && [ -s "ref_vars.snplist" ]; then
            SHARED_COUNT=\$(wc -l < ref_vars.snplist)
            echo "  Reference panel variants: \${SHARED_COUNT}"
            EXTRACT_REF="--extract ref_vars.snplist"
        else
            echo "WARNING: Could not parse reference panel. Proceeding without intersection."
        fi
    else
        echo "WARNING: No reference panel provided. Skipping intersection step."
        echo "  PCA results will be limited without population reference."
    fi

    # Step 2: Convert VCF to plink2 format (assigns position-based IDs)
    # Must convert BEFORE extraction — plink2 applies --extract before
    # --set-all-var-ids, so VCF with '.' IDs won't match the snplist.
    plink2 \\
        --vcf ${vcf} \\
        --set-all-var-ids '@:#' \\
        --new-id-max-allele-len 100 \\
        --make-pgen \\
        --out ${meta.id}_converted \\
        --threads ${task.cpus} \\
        --memory \$(( ${task.memory.toMega()} - 500 )) \\
        --chr 1-22 \\
        --allow-extra-chr \\
        --output-chr chrM 2>&1

    # Step 2b: Extract reference panel intersection (IDs now match)
    if [ -n "\${EXTRACT_REF}" ]; then
        plink2 \\
            --pfile ${meta.id}_converted \\
            \${EXTRACT_REF} \\
            --make-pgen \\
            --out ${meta.id} \\
            --threads ${task.cpus} \\
            --memory \$(( ${task.memory.toMega()} - 500 )) 2>&1
    else
        plink2 \\
            --pfile ${meta.id}_converted \\
            --make-pgen \\
            --out ${meta.id} \\
            --threads ${task.cpus} \\
            --memory \$(( ${task.memory.toMega()} - 500 )) 2>&1
    fi

    # Step 3: LD pruning on autosomal SNPs
    LD_PRUNE_OK=false
    if plink2 \\
        --pfile ${meta.id} \\
        --indep-pairwise 50 5 0.2 \\
        --out ${meta.id}_ld \\
        --threads ${task.cpus} \\
        --memory \$(( ${task.memory.toMega()} - 500 )) 2>&1; then
        LD_PRUNE_OK=true
        EXTRACT_ARGS="--extract ${meta.id}_ld.prune.in"
    else
        EXTRACT_ARGS=""
    fi

    # Step 4: PCA (requires >=2 samples; expected to fail for single-sample)
    PCA_OK=false
    if plink2 \\
        --pfile ${meta.id} \\
        \${EXTRACT_ARGS} \\
        --pca 10 \\
        --out ${meta.id}_pca \\
        --threads ${task.cpus} \\
        --memory \$(( ${task.memory.toMega()} - 500 )) 2>&1; then
        PCA_OK=true
    fi

    # Count variants actually loaded into plink2 (post-extraction)
    LOADED_COUNT=0
    if [ -f "${meta.id}.pvar" ]; then
        LOADED_COUNT=\$(grep -c -v '^#' ${meta.id}.pvar 2>/dev/null || echo "0")
    fi

    # Step 5: Summary file with QC metrics
    echo -e "metric\\tvalue" > ${meta.id}_ancestry.tsv
    if [ "\${LOADED_COUNT}" -gt 0 ]; then
        if [ "${has_ref}" = "true" ] && [ -n "\${EXTRACT_REF}" ]; then
            echo -e "ref_panel_shared_snps\\t\${LOADED_COUNT}" >> ${meta.id}_ancestry.tsv
        else
            echo -e "autosomal_snps\\t\${LOADED_COUNT}" >> ${meta.id}_ancestry.tsv
        fi
    fi
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
