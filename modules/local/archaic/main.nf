/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ARCHAIC — Neanderthal and Denisovan segments with IBDmix (step 43)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('archaic' in --tools, with --archaic_data). Three processes, as in
    scripts/43-archaic.sh:
      ARCHAIC_GENOTYPES  bin/archaic_report.py genotypes: IBDmix's genotype file of
                         each chromosome (the archaic genomes, the sample's VCF and
                         the 1000 Genomes stand-ins of --archaic_pop)
      ARCHAIC_IBDMIX     IBDmix against each archaic genome, LOD 4, the options of
                         Chen et al. 2020
      ARCHAIC_REPORT     bin/archaic_report.py report: the segments of 50 kb or
                         more, the share, the ten longest with their genes

    Equivalent to: scripts/43-archaic.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ARCHAIC_GENOTYPES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)
    path(archaic_data)   // the folder setup.sh --archaic writes

    output:
    tuple val(meta), path("gt"), emit: genotypes
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    archaic_report.py genotypes --table-dir ${archaic_data} --vcf ${vcf} --pop ${params.archaic_pop} --out-dir gt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir gt
    printf 'chrom\\tfirst\\tlast\\tsites\\tsites_with_your_alt\\tsites_other_allele\\nchr22\\t1\\t2\\t1\\t0\\t0\\n' > gt/spans.tsv
    printf 'SAMPLE\\n' > gt/samples.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process ARCHAIC_IBDMIX {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(gt)

    output:
    tuple val(meta), path("ibdmix"), emit: segments
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // One chromosome after another, the three archaic genomes of each; the
    // same options as scripts/43-archaic.sh.
    """
    mkdir ibdmix
    cp ${gt}/spans.tsv ibdmix/spans.tsv
    for f in ${gt}/chr*.gt.gz; do
        c=\$(basename "\$f" .gt.gz)
        gzip -dc "\$f" > "\$c.gt"
        for a in Altai Vindija Denisova; do
            ibdmix -g "\$c.gt" -s ${gt}/samples.txt -n "\$a" -o "ibdmix/\$c.\$a.txt" \\
                -d 4 -m 1 -a 0.01 -e 0.002 -c 2
        done
        rm -f "\$c.gt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ibdmix: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir ibdmix
    cp ${gt}/spans.tsv ibdmix/spans.tsv
    for a in Altai Vindija Denisova; do printf 'ID\\tchrom\\tstart\\tend\\tslod\\n' > "ibdmix/chr22.\$a.txt"; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ibdmix: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process ARCHAIC_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/archaic" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(ibdmix)
    path(genes)   // GENCODE's gene lines (--hla_genes), or [] for none

    output:
    tuple val(meta), path("${meta.id}_archaic.tsv"),          emit: summary
    tuple val(meta), path("${meta.id}_archaic_segments.tsv"), emit: segments
    tuple val(meta), path("${meta.id}_archaic.txt"),          emit: text
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genes_arg = genes ? "--genes ${genes}" : ''
    """
    archaic_report.py report --sample ${meta.id} --ibdmix-dir ${ibdmix} --spans ${ibdmix}/spans.tsv \\
        --pop ${params.archaic_pop} ${genes_arg} --out-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'key\\tvalue\\nsample\\t${meta.id}\\n' > ${meta.id}_archaic.tsv
    printf 'chrom\\tstart\\tend\\tlength_bp\\tancestry\\tmatched\\tlod\\tgenes\\n' > ${meta.id}_archaic_segments.tsv
    echo stub > ${meta.id}_archaic.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
