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
    cpic/, roh/, mito/, telomere/, hla_t1k/, hla_optitype/, kir/, pypgx/, annotsv/, traits/, chip/, couple/, ...); a step that did
    not run, or whose task failed and was ignored, has no file, and its card
    says "Not run". The called or given VCF comes separately, because its
    name is the user's. indexcov's folder becomes indexcov/ itself (the
    sex check reads indexcov/indexcov-indexcov.ped), and CPIC_LOOKUP's
    PharmCAT comparison goes beside the pypgx summary, where
    collect_summary.py looks for it. assets/acmg_sf.tsv (the ACMG SF rules)
    and assets/carrier_genes.tsv (the carrier list) come in the same list
    and go to assets/, which --assets names. The report's charts also go
    to ${id}_charts.json beside it, for the live page (docs/24-html-report.md).

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
    tuple val(meta), path("${meta.id}_charts.json"),   emit: charts
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
            ${id}_clinvar_hits.vcf|${id}_clinvar_cdna_check.tsv) place clinvar "\$f" ;;
            ${id}.report.json|${id}.report.html)          place pharmcat "\$f" ;;
            ${id}_phenotypes.tsv|${id}_cpic_recommendations.txt) place cpic "\$f" ;;
            ${id}_clinical.vcf.gz)                        place clinical "\$f" ;;
            ${id}.cpsr.grch38.html)                       place cpsr "\$f" ;;
            ${id}_prioritized.vcf.gz|${id}_compound_het_pairs.tsv) place slivar "\$f" ;;
            ${id}_roh.txt)                                place roh "\$f" ;;
            ${id}_haplogroup.txt|${id}_haplocheck.txt)    place mito "\$f" ;;
            ${id}_y_haplogroup.txt)                       place y_haplogroup "\$f" ;;
            ${id}_y_str.tsv)                              place y_str "\$f" ;;
            ${id}_paraphase_copy_number.tsv|${id}_paraphase_variants.tsv|${id}_smn_copy_number.tsv|${id}_hba_copy_number.tsv) place paralogs "\$f" ;;
            ${id}.mosdepth.summary.txt)                   place coverage "\$f" ;;
            ${id}_sample_qc.tsv)                          place qc "\$f" ;;
            ${id}_prs_summary.tsv)                        place prs "\$f" ;;
            ${id}_ancestry.tsv)                           place ancestry "\$f" ;;
            ${id}_ancestry_fine.tsv|${id}_painting.png)   place ancestry_fine "\$f" ;;
            ${id})                                        place telomere "\$f" ;;
            ${id}_chrM_filtered.vcf.gz|${id}_chrM_annotated.tsv) place mito "\$f" ;;
            ${id}_hla_genotype.tsv)                       place hla_t1k "\$f" ;;
            ${id}_optitype_result.tsv|${id}_hla_concordance.tsv) place hla_optitype "\$f" ;;
            ${id}_kir_genotype.tsv)                       place kir "\$f" ;;
            ${id}_eh.vcf|${id}_eh_stranger.vcf|${id}_rfc1_motif.tsv) place expansion_hunter "\$f" ;;
            diploidSV.vcf.gz)                             place manta "\$f" ;;
            ${id}_sv.vcf.gz)                              place delly "\$f" ;;
            ${id}_cnvs.txt)                               place cnvpytor "\$f" ;;
            ${id}_sv_consensus.vcf.gz)                    place sv_merged "\$f" ;;
            ${id}_sv_class4_5_checked.tsv)                place annotsv "\$f" ;;
            ${id}_pypgx_summary.tsv|${id}_pharmcat_comparison.tsv|${id}_cyp2d6_depth_check.tsv) place pypgx "\$f" ;;
            indexcov)                                     ln -s "\$(readlink -f "\$f")" "\$S/indexcov" ;;
            ${id}_cyp2d6.tsv)                             place cyrius "\$f" ;;
            ${id}_pgx_consensus.tsv)                      place pgx_consensus "\$f" ;;
            ${id}_consumer_export.tsv)                    place consumer_export "\$f" ;;
            ${id}_traits.tsv)                             place traits "\$f" ;;
            ${id}_bloodtype.tsv)                          place bloodtype "\$f" ;;
            ${id}_archaic.tsv)                            place archaic "\$f" ;;
            ${id}_ancient.tsv|${id}_ancient_info.tsv|${id}_ancient_pca.svg) place ancient "\$f" ;;
            ${id}_lpa_kiv2.tsv)                           place lpa "\$f" ;;
            ${id}_vntyper.tsv)                            place vntyper "\$f" ;;
            ${id}_viral.tsv)                              place viral "\$f" ;;
            ${id}_exomiser_top.tsv|${id}_exomiser_terms.tsv) place exomiser "\$f" ;;
            ${id}_chip.tsv|${id}_chip_depth.tsv)          place chip "\$f" ;;
            ${id}_mocha.tsv|${id}_mocha_summary.tsv)      place mosaic "\$f" ;;
            ${id}_couple_carrier.tsv)                     place couple "\$f" ;;
            acmg_sf.tsv|carrier_genes.tsv)                mkdir -p assets && ln -s "\$(readlink -f "\$f")" assets/ ;;
            *) echo "ERROR: HTML_REPORT got an input it has no place for: \$(basename "\$f")" >&2; exit 1 ;;
        esac
    done
    find layout | sort

    render_report.py \\
        --assets assets \\
        --sample ${id} \\
        --sample-dir "\$S" \\
        --json summary.json \\
        --declared-sex "${meta.sex ?: ''}" \\
        -o ${id}_report.html \\
        --charts-json ${id}_charts.json

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_report.html summary.json ${meta.id}_charts.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
