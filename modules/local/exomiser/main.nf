/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXOMISER — Rank the variants that could explain a person's symptoms (HPO terms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in: runs when --hpo names the symptoms of one sample (--exomiser_proband,
    or the only sample of the samplesheet). Three processes, as in
    scripts/48-exomiser.sh:
      EXOMISER_PREPARE  bin/exomiser_report.py prepare: checks the terms against
                        the HPO release of the phenotype data and writes the
                        analysis (Exomiser's exome preset)
      EXOMISER          Exomiser on the sample's VCF, offline, with its data
                        folder (--exomiser_data, setup.sh --exomiser-data)
      EXOMISER_REPORT   bin/exomiser_report.py report: the top ten ranks with
                        their variants, and the same in words

    The family mode of step 48 (EXOMISER_FAMILY) is bash-only.

    Equivalent to: scripts/48-exomiser.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EXOMISER_PREPARE {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/exomiser" }, mode: params.publish_dir_mode, pattern: '*_terms.tsv'

    input:
    tuple val(meta), path(vcf)
    val(hpo)              // --hpo: the terms, comma separated
    path(exomiser_data)   // --exomiser_data: reference/exomiser-<release> (setup.sh --exomiser-data)

    output:
    tuple val(meta), path(vcf), path("${meta.id}_analysis.yml"), emit: analysis
    tuple val(meta), path("${meta.id}_exomiser_terms.tsv"),      emit: terms
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The HPO ontology of the phenotype data: hp.obo, or phenix/hp.obo in the
    // releases that keep it there. Without one the terms are checked for their
    // form only, as in step 48.
    """
    obo=""
    for o in ${exomiser_data}/*_phenotype/hp.obo ${exomiser_data}/*_phenotype/phenix/hp.obo; do
        if [ -s "\$o" ]; then obo="\$o"; break; fi
    done
    exomiser_report.py prepare --hpo '${hpo}' \${obo:+--obo "\$obo"} --sample ${meta.id} --vcf "${vcf}" \\
        --assembly ${params.exomiser_assembly} --out-dir results --analysis ${meta.id}_analysis.yml \\
        --terms ${meta.id}_exomiser_terms.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf '{}\\n' > ${meta.id}_analysis.yml
    printf 'hpo_id\\tlabel\\tgiven_as\\nHP:0001250\\t\\tHP:0001250\\n' > ${meta.id}_exomiser_terms.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process EXOMISER {
    tag "$meta.id"
    cpus 4
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    // The image's entrypoint is bash, which would read Nextflow's command as a
    // script to run. A process's containerOptions replaces the docker
    // profile's '--network none', so it is repeated here: Exomiser reads only
    // its data folder. Singularity runs the command without the entrypoint.
    containerOptions { workflow.containerEngine == 'singularity' ? '' : '--entrypoint "" --network none' }

    publishDir { "${params.outdir}/${meta.id}/exomiser" }, mode: params.publish_dir_mode,
        saveAs: { f -> f.startsWith('results/') ? f.substring(8) : null }

    input:
    tuple val(meta), path(vcf), path(analysis)
    path(exomiser_data)

    output:
    tuple val(meta), path("results/${meta.id}_exomiser.variants.tsv"), emit: variants
    path "results/${meta.id}_exomiser.*",                              emit: results
    path "versions.yml",                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exomiser reads its data folder and versions from these variables (its
    // installation docs); the version is the one <version>_phenotype folder.
    def xmx = Math.max(1, (task.memory.toGiga() * 3 / 4) as int)
    def asm = params.exomiser_assembly.toUpperCase()
    """
    set -- ${exomiser_data}/*_phenotype
    if [ \$# -ne 1 ] || [ ! -d "\$1" ]; then
        echo "ERROR: ${exomiser_data} must hold one <version>_phenotype folder (setup.sh --exomiser-data)" >&2
        exit 1
    fi
    version=\$(basename "\$1" _phenotype)
    mkdir -p results
    JAVA_OPTS="-Xmx${xmx}g" EXOMISER_DATA_DIRECTORY="\$(readlink -f ${exomiser_data})" \\
        EXOMISER_${asm}_DATA_VERSION="\$version" EXOMISER_PHENOTYPE_DATA_VERSION="\$version" \\
        exomiser analyse --analysis ${analysis}
    test -f results/${meta.id}_exomiser.variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        exomiser: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p results
    printf '#RANK\\tGENE_SYMBOL\\n' > results/${meta.id}_exomiser.variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        exomiser: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process EXOMISER_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/exomiser" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(variants), path(terms)

    output:
    tuple val(meta), path("${meta.id}_exomiser_top.tsv"), emit: top
    tuple val(meta), path("${meta.id}_exomiser.txt"),     emit: text
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    exomiser_report.py report --sample ${meta.id} --variants ${variants} --terms ${terms} \\
        --out ${meta.id}_exomiser_top.tsv --text ${meta.id}_exomiser.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'rank\\tgene\\n' > ${meta.id}_exomiser_top.tsv
    printf 'Symptoms searched: HP:0001250.\\n' > ${meta.id}_exomiser.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
