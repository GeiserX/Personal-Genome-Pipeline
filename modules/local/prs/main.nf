/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRS — Polygenic scores with pgsc_calc, and the ancestry of step 26
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRS_PREPARE      writes every file of --pgs_scoring as the custom GRCh38
                     scoring file pgsc_calc reads (bin/collect_summary.py
                     prs-format), labelled from assets/pgs_scores.tsv
    PRS_SCORE_SITES  genotypes the score positions (and, with a panel, the
                     panel's) from the gVCF, so a site where the sample
                     matches the reference is a real 0/0; without a gVCF it
                     cuts the variant-only VCF to those positions
    PRS              runs pgsc_calc (PGSC_CALC_VERSION) on the host: it is a
                     Nextflow pipeline of its own and starts its own
                     containers, with the images versions.env pins and no
                     network; with --ancestry_ref it projects the sample onto
                     the panel and adjusts each score for ancestry
    PRS_SUMMARY      <id>_prs_summary.tsv (bin/collect_summary.py prs-table),
                     and with the panel the ancestry table of step 26

    Equivalent to: scripts/25-prs.sh (and scripts/26-ancestry.sh)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PRS_PREPARE {
    label 'process_single'

    input:
    path(scoring_dir)
    path(labels)

    output:
    path("pgs"),               emit: scores
    path("score_alleles.tsv"), emit: alleles
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    collect_summary.py prs-format --scores ${scoring_dir} --labels ${labels} --out pgs --alleles score_alleles.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    mkdir -p pgs
    touch pgs/PGS000000.txt.gz score_alleles.tsv
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

// gvcf2vcf expands each reference block over a position into a 0/0 record
// with the reference base, -T keeps the positions, --trim-alt-alleles drops
// <*>, and a no-call is dropped, so a position without coverage stays
// missing. A 0/0 record has ALT '.', which no allele matches: the awk sets
// ALT to the position's first candidate (a score allele, the panel's ALT)
// that is not the reference. The same commands as step 25.
process PRS_SCORE_SITES {
    tag "$meta.id"
    label 'process_low'

    input:
    // kind gvcf: the sample's gVCF, its reference blocks expanded at the
    // positions; kind vcf: its variant-only VCF, cut to the positions.
    tuple val(meta), path(gvcf), path(gvcf_index), val(input_kind)
    path(score_alleles)
    path(panel_sites)    // [] without an ancestry panel
    path(reference)
    path(reference_fai)  // staged beside the FASTA

    output:
    tuple val(meta), path("${meta.id}_score_sites.vcf.gz"), val(input_kind), emit: vcf
    path "versions.yml",                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def panel = panel_sites ? "cut -f1,2,4 ${panel_sites};" : ''
    """
    { cat ${score_alleles}; ${panel} } | LC_ALL=C sort -u -k1,1 -k2,2n -k3,3 > alleles.tsv
    cut -f1,2 alleles.tsv | uniq > sites.tsv
    if [ ! -s sites.tsv ]; then
        echo "ERROR: no score position to genotype" >&2
        exit 1
    fi

    if [ "${input_kind}" = vcf ]; then
        bcftools view -T sites.tsv -Ov ${gvcf} | gzip -c > ${meta.id}_score_sites.vcf.gz
    else
    bcftools convert --gvcf2vcf -f ${reference} -R sites.tsv -Ou ${gvcf} \\
        | bcftools view -T sites.tsv --trim-alt-alleles -i 'GT!="mis"' -Ov \\
        | awk -F'\\t' -v OFS='\\t' -v alleles=alleles.tsv '
            BEGIN { while ((getline l < alleles) > 0) { split(l, f, "\\t"); k = f[1] ":" f[2]; ea[k] = ea[k] " " f[3] } }
            /^#/ { print; next }
            \$5 == "." {
                n = split(ea[\$1 ":" \$2], c, " ")
                for (i = 1; i <= n; i++) if (c[i] != \$4) { \$5 = c[i]; break }
            }
            { print }' | gzip -c > ${meta.id}_score_sites.vcf.gz
    fi
    rm -f alleles.tsv sites.tsv

    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_score_sites.vcf.gz
    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

// Runs on the host (scripts/ci/gen-containers-config.sh gives it no
// container): pgsc_calc is a Nextflow pipeline and starts its own containers,
// each with --network none. ext.pipeline_version and ext.pipeline_images
// come from versions.env through conf/containers.config. With --pgsc_calc
// (the checkout setup.sh makes) it runs offline; without, Nextflow fetches
// pgscatalog/pgsc_calc at that release and its nf-schema plugin first.
process PRS {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/prs" }, mode: params.publish_dir_mode,
        saveAs: { f -> f == 'versions.yml' ? null : f }

    input:
    // input_kind: gvcf when vcf is PRS_SCORE_SITES' output, vcf for the
    // sample's variants-only VCF.
    tuple val(meta), path(vcf, stageAs: 'target/*'), val(input_kind)
    path(scores, stageAs: 'pgs')
    path(panel)          // [] without an ancestry panel

    output:
    tuple val(meta), path("pgsc_calc"), val(input_kind), emit: results
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def release  = task.ext.pipeline_version
    def pipeline = params.pgsc_calc ? "${params.pgsc_calc}/main.nf" : "pgscatalog/pgsc_calc -r ${release}"
    def ancestry = panel ? "--run_ancestry \$(readlink -f ${panel})" : ''
    def labels = task.ext.pipeline_images.tokenize(';').collect { kv ->
        def i = kv.indexOf('=')
        def image = kv.substring(i + 1)
        // pgsc_calc sets docker.registry to quay.io: name Docker Hub images in full.
        def first = image.tokenize('/')[0]
        if (!(first.contains('.') || first.contains(':') || first == 'localhost')) { image = "docker.io/${image}" }
        "\"    withLabel: '${kv.substring(0, i)}' { ext.docker = '${image}'; ext.docker_version = '' }\""
    }.join(' ')
    def engine  = workflow.containerEngine ?: 'docker'
    // With the checkout setup.sh makes (it installs the nf-schema plugin
    // pgsc_calc needs too) nothing is fetched; without, Nextflow pulls both.
    def offline = params.pgsc_calc ? 'export NXF_OFFLINE=true' : ''
    """
    VCF=\$(readlink -f target/*)
    case "\$VCF" in
        *.vcf.gz) PREFIX=\${VCF%.vcf.gz} ;;
        *) echo "ERROR: PRS needs a .vcf.gz, got \$VCF" >&2; exit 1 ;;
    esac
    printf 'sampleset,path_prefix,chrom,format\\nsample,%s,,vcf\\n' "\$PREFIX" > samplesheet.csv
    # The images of versions.env for pgsc_calc's process labels, and no network for its containers.
    printf '%s\\n' 'process {' ${labels} '}' "docker.runOptions = '-u \$(id -u):\$(id -g) --network none'" > images.config
    cat images.config
    ${offline}
    rc=0
    mkdir -p pgsc_calc
    if [ "\$(gzip -dc "\$VCF" | grep -vc '^#' || true)" -eq 0 ]; then
        # Not one score position is in the input: nothing for pgsc_calc to match.
        echo "None of the score positions is in this input; no score." | tee pgsc_calc/ZERO_MATCHES
        echo "pgsc_calc not run: no score position in the input" > pgsc_calc.log
        touch pgsc_calc.nextflow.log
    else
    nextflow -log pgsc_calc.nextflow.log run ${pipeline} \\
        -profile ${engine} -c images.config -work-dir work -ansi-log false \\
        --input samplesheet.csv --target_build GRCh38 \\
        --scorefile "\$(readlink -f pgs)/*.txt.gz" \\
        ${ancestry} \\
        --outdir pgsc_calc \\
        --max_cpus ${task.cpus} --max_memory '${task.memory.toGiga()}.GB' > pgsc_calc.log 2>&1 || rc=\$?
    fi
    cat pgsc_calc.log
    mv pgsc_calc.log pgsc_calc.nextflow.log pgsc_calc/
    if [ "\$rc" -ne 0 ]; then
        if grep -qE 'ZeroMatchesError|No match candidates found for any scoring files|All scores fail to meet match threshold' pgsc_calc/pgsc_calc*.log; then
            echo "None of the scores matched enough of its variants in this input; no score." > pgsc_calc/ZERO_MATCHES
        else
            echo "ERROR: pgsc_calc failed (exit \$rc); see pgsc_calc/pgsc_calc.log" >&2
            exit 1
        fi
    fi
    # Its run reports carry the time of the run in their names.
    rm -rf work .nextflow pgsc_calc/pipeline_info

    printf '"%s":\\n    pgsc_calc: %s\\n' "${task.process}" "${release}" > versions.yml
    """

    stub:
    """
    mkdir -p pgsc_calc
    printf '"%s":\\n    pgsc_calc: %s\\n' "${task.process}" "${task.ext.pipeline_version}" > versions.yml
    """
}

process PRS_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/prs" }, mode: params.publish_dir_mode, pattern: '*_prs_summary.tsv'
    publishDir { "${params.outdir}/${meta.id}/ancestry" }, mode: params.publish_dir_mode, pattern: '*_ancestry.tsv'

    input:
    tuple val(meta), path(results), val(input_kind)
    path(scores, stageAs: 'pgs')
    val(panel_name)      // '' without an ancestry panel

    output:
    tuple val(meta), path("${meta.id}_prs_summary.tsv"),             emit: summary
    tuple val(meta), path("${meta.id}_ancestry.tsv"), optional: true, emit: ancestry
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def anc = panel_name ? "--panel '${panel_name}' --ancestry-out ${meta.id}_ancestry.tsv" : ''
    """
    ZERO=""
    if [ -f ${results}/ZERO_MATCHES ]; then ZERO=--zero-matches; fi
    collect_summary.py prs-table --sample ${meta.id} --results ${results} --sampleset sample \\
        --scores pgs --input-kind ${input_kind} \$ZERO ${anc} --out ${meta.id}_prs_summary.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'Condition\\tPGS_ID\\tScore_SUM\\tVariants_Matched\\tVariants_Total\\tMatched_Pct\\tPercentile\\tAncestry_Group\\tInput\\n' > ${meta.id}_prs_summary.tsv
    if [ -n "${panel_name}" ]; then printf 'key\\tvalue\\nsample\\t%s\\n' ${meta.id} > ${meta.id}_ancestry.tsv; fi
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
