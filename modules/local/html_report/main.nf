/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HTML_REPORT — The sample's HTML report, rendered by bin/render_report.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The same code scripts/24-html-report.sh runs: bin/collect_summary.py
    reads the step outputs into summary.json (the bash step's name for it),
    and render_report.py writes the report from it, so a number on the
    Nextflow report is the number the bash report shows for the same file.

    The inputs are the outputs of the steps that ran for this sample, as one
    list (main.nf joins them, so the report waits for each). They are linked
    into the folder layout collect_summary.py reads (clinvar/, pharmcat/,
    cpic/, roh/, mito/, ...); a step that did not run has no file, and its
    card says "Not run". The called or given VCF comes separately, because
    its name is the user's.

    Equivalent to: scripts/24-html-report.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HTML_REPORT {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf, stageAs: 'input/sample.vcf.gz'), path(outputs, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}_report.html"), emit: html_report
    tuple val(meta), path("summary.json"),             emit: summary
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def id = meta.id
    """
    # Where collect_summary.py looks for each output (bin/collect_summary.py)
    S=layout/${id}
    mkdir -p "\$S"
    place() { mkdir -p "\$S/\$1" && ln -s "\$(readlink -f "\$2")" "\$S/\$1/\$(basename "\$2")"; }
    place vcf input/sample.vcf.gz && mv "\$S/vcf/sample.vcf.gz" "\$S/vcf/${id}.vcf.gz"
    for f in input/*; do
        case "\$(basename "\$f")" in
            sample.vcf.gz) ;;
            ${id}_clinvar_hits.vcf)                       place clinvar "\$f" ;;
            ${id}.report.json|${id}.report.html)          place pharmcat "\$f" ;;
            ${id}_phenotypes.tsv|${id}_cpic_recommendations.txt) place cpic "\$f" ;;
            ${id}_clinical.vcf.gz)                        place clinical "\$f" ;;
            ${id}.cpsr.grch38.html)                       place cpsr "\$f" ;;
            ${id}_prioritized.vcf.gz)                     place slivar "\$f" ;;
            ${id}_roh.txt)                                place roh "\$f" ;;
            ${id}_haplogroup.txt)                         place mito "\$f" ;;
            ${id}.mosdepth.summary.txt)                   place coverage "\$f" ;;
            ${id}_sample_qc.tsv)                          place qc "\$f" ;;
            *) echo "ERROR: HTML_REPORT got an input it has no place for: \$(basename "\$f")" >&2; exit 1 ;;
        esac
    done
    find layout | sort

    render_report.py \\
        --sample ${id} \\
        --sample-dir "\$S" \\
        --json summary.json \\
        --declared-sex "${meta.sex ?: ''}" \\
        -o ${id}_report.html

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_report.html summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
