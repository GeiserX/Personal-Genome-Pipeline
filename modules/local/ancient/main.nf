/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANCIENT — Step 44: where the genome falls among ancient genomes (AADR, smartpca)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANCIENT_SITES   genotypes the panel's GRCh38 sites from the gVCF: each
                    reference block over a site becomes a record with the
                    reference base, and a site without reads has no line
    ANCIENT_GENO    the sample as one more column of the panel's present-day
                    EIGENSTRAT genotypes (bin/ancient_report.py person-geno)
    ANCIENT_PCA     smartpca: axes from the present-day West Eurasians, the
                    sample projected (lsqproject)
    ANCIENT_REPORT  the nearest ancient groups, the summary and the plot, on
                    the panel's axes (bin/ancient_report.py report)

    The panel is the folder setup.sh --ancient-panel writes
    (reference/ancient/aadr_<release>). Runs for the samples with a gVCF only;
    main.nf logs the samples it skips.

    Equivalent to: scripts/44-ancient.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANCIENT_SITES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(gvcf), path(gvcf_index)
    path(panel)
    path(reference)
    path(reference_fai)  // staged beside the FASTA

    output:
    tuple val(meta), path("${meta.id}_ancient_genotypes.tsv"), emit: genotypes
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The same commands as scripts/44-ancient.sh
    """
    n=\$(bcftools query -l ${gvcf} | wc -l | tr -d ' ')
    if [ "\$n" -ne 1 ]; then
        echo "ERROR: Sample '${meta.id}': ${gvcf} holds \$n samples; step 44 reads a one-sample gVCF (keep one with bcftools view -s <sample>)." >&2
        exit 1
    fi
    awk -F'\\t' -v OFS='\\t' 'NR > 1 {print \$2, \$3}' ${panel}/sites_grch38.tsv > sites.tsv
    bcftools convert --gvcf2vcf -f ${reference} -R sites.tsv -Ou ${gvcf} \\
        | bcftools view -T sites.tsv --trim-alt-alleles -i 'GT!="mis"' -Ou \\
        | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%TGT]\\n' > ${meta.id}_ancient_genotypes.tsv
    rm -f sites.tsv

    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_ancient_genotypes.tsv
    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process ANCIENT_GENO {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(genotypes)
    path(panel)

    output:
    tuple val(meta), path("merged.geno"), path("merged.ind"), path("stats.tsv"), emit: merged
    path "versions.yml",                                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ancient_report.py person-geno --sites ${panel}/sites_grch38.tsv --genotypes ${genotypes} \\
        --snp ${panel}/moderns.snp --moderns-geno ${panel}/moderns.geno --out merged.geno --stats stats.tsv
    { cat ${panel}/moderns.ind; printf '%s\\tU\\tPerson\\n' '${meta.id}'; } > merged.ind

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch merged.geno merged.ind stats.tsv
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process ANCIENT_PCA {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(geno), path(ind), path(stats)
    path(panel)

    output:
    tuple val(meta), path("run.evec"), path(stats), emit: evec
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The parameters of scripts/44-ancient.sh and setup.sh's ancient_smartpca_par
    """
    printf '%s\\n' "genotypename: ${geno}" "snpname: ${panel}/moderns.snp" "indivname: ${ind}" \\
        "poplistname: ${panel}/poplist.txt" "evecoutname: run.evec" "evaloutname: run.eval" "lsqproject: YES" \\
        "numoutevec: 10" "numoutlieriter: 0" "hashcheck: NO" "numthreads: ${task.cpus}" > run.par
    smartpca -p run.par > smartpca.log 2>&1 || { tail -n 20 smartpca.log >&2; exit 1; }

    printf '"%s":\\n    eigensoft: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch run.evec
    printf '"%s":\\n    eigensoft: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process ANCIENT_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/ancient" }, mode: params.publish_dir_mode, pattern: "${meta.id}_ancient*"

    input:
    tuple val(meta), path(evec), path(stats)
    path(panel)

    output:
    tuple val(meta), path("${meta.id}_ancient.tsv"),      emit: table
    tuple val(meta), path("${meta.id}_ancient_info.tsv"), emit: info
    tuple val(meta), path("${meta.id}_ancient_pca.svg"),  emit: svg
    tuple val(meta), path("${meta.id}_ancient.txt"),      emit: text
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    release=\$(awk -F'\\t' '\$1 == "release" {print \$2}' ${panel}/DONE)
    ancient_report.py report --panel-evec ${panel}/panel.evec --run-evec ${evec} --ancients ${panel}/ancients.tsv \\
        --sample ${meta.id} --release "\$release" --stats ${stats} \\
        --out-tsv ${meta.id}_ancient.tsv --out-info ${meta.id}_ancient_info.tsv \\
        --out-txt ${meta.id}_ancient.txt --out-svg ${meta.id}_ancient_pca.svg

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'rank\\tgroup\\tn\\tdistance\\tpc1\\tpc2\\tdate_bp\\tdates\\tcountry\\tlocality\\n' > ${meta.id}_ancient.tsv
    printf 'key\\tvalue\\n' > ${meta.id}_ancient_info.tsv
    touch ${meta.id}_ancient_pca.svg ${meta.id}_ancient.txt
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
