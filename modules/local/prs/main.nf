/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRS — Polygenic scores with pgsc_calc, and the ancestry of step 26
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRS_PREPARE      writes every file of --pgs_scoring as the custom GRCh38
                     scoring file pgsc_calc reads (bin/collect_summary.py
                     prs-format), labelled from assets/pgs_scores.tsv: once
                     with the chrX rows (set with_x, for samples whose sex the
                     samplesheet gives) and once without (set autosomes), as
                     the samples of the run need
    PRS_SCORE_SITES  genotypes the score positions (and, with a panel, the
                     panel's) from the gVCF, so a site where the sample
                     matches the reference is a real 0/0; without a gVCF it
                     cuts the variant-only VCF to those positions
    PRS              runs pgsc_calc (PGSC_CALC_VERSION) on the host: it is a
                     Nextflow pipeline of its own and starts its own
                     containers, with the images versions.env pins and no
                     network; with --ancestry_ref it projects the sample onto
                     the panel and adjusts each score for ancestry. With the
                     samplesheet's sex its plink2 gets the sex (--update-sex),
                     which it needs to read chrX
    PRS_SUMMARY      <id>_prs_summary.tsv (bin/collect_summary.py prs-table),
                     and with the panel the ancestry table of step 26

    Equivalent to: scripts/25-prs.sh (and scripts/26-ancestry.sh)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PRS_PREPARE {
    tag "$score_set"
    label 'process_single'

    input:
    val(score_set)       // with_x: keep the chrX rows; autosomes: drop them
    path(scoring_dir)
    path(labels)

    output:
    tuple val(score_set), path("pgs"), path("score_alleles.tsv"), emit: scores
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def keep_x = score_set == 'with_x' ? '--keep-x' : ''
    """
    collect_summary.py prs-format --scores ${scoring_dir} --labels ${labels} --out pgs --alleles score_alleles.tsv ${keep_x}

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
    // score_alleles: PRS_PREPARE's, of the sample's set.
    tuple val(meta), path(gvcf), path(gvcf_index), val(input_kind), path(score_alleles)
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
// (the copy setup.sh makes) it runs offline; without, the task fetches the
// release's archive (checked against ext.pipeline_sha256) and Nextflow its
// nf-schema plugin first.
process PRS {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/prs" }, mode: params.publish_dir_mode,
        saveAs: { f -> f == 'versions.yml' ? null : f }

    input:
    // input_kind: gvcf when vcf is PRS_SCORE_SITES' output, vcf for the
    // sample's variants-only VCF. scores: PRS_PREPARE's, of the sample's set.
    tuple val(meta), path(vcf, stageAs: 'target/*'), val(input_kind), path(scores, stageAs: 'pgs')
    path(panel)          // [] without an ancestry panel

    output:
    tuple val(meta), path("pgsc_calc"), val(input_kind), emit: results
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def release  = task.ext.pipeline_version
    // Absolute: the task runs in its own work folder, where a relative
    // --pgsc_calc would name nothing.
    def pipeline = params.pgsc_calc ? "${file(params.pgsc_calc).toAbsolutePath()}/main.nf" : 'pgsc_calc-src/main.nf'
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
    // The sex for pgsc_calc's plink2 (--update-sex: the VCF's sample name, 1
    // male or 2 female), added to the ext.args its conf/modules.config gives
    // PLINK2_VCF at this release; the folder is mounted into its containers.
    def sex_code = meta.sex == 'male' ? '1' : (meta.sex == 'female' ? '2' : '')
    // With the copy setup.sh makes (it installs the nf-schema plugin pgsc_calc
    // needs too) nothing is fetched. Without, the task fetches GitHub's
    // archive of the release (checked against PGSC_CALC_SHA256, as setup.sh
    // and step 25 do; not `nextflow run owner/repo`, whose GitHub API calls
    // shared runners exhaust) and Nextflow fetches the plugin.
    def offline = params.pgsc_calc ? 'export NXF_OFFLINE=true' :
        "curl -fsSL --retry 3 -o pgsc_calc.tar.gz https://github.com/PGScatalog/pgsc_calc/archive/refs/tags/${release}.tar.gz\n" +
        "    if command -v sha256sum >/dev/null; then SUM=sha256sum; else SUM='shasum -a 256'; fi\n" +
        "    echo '${task.ext.pipeline_sha256}  pgsc_calc.tar.gz' | \$SUM -c -\n" +
        "    mkdir pgsc_calc-src && tar -xzf pgsc_calc.tar.gz -C pgsc_calc-src --strip-components 1 && rm -f pgsc_calc.tar.gz"
    """
    VCF=\$(readlink -f target/*)
    case "\$VCF" in
        *.vcf.gz) PREFIX=\${VCF%.vcf.gz} ;;
        *) echo "ERROR: PRS needs a .vcf.gz, got \$VCF" >&2; exit 1 ;;
    esac
    printf 'sampleset,path_prefix,chrom,format\\nsample,%s,,vcf\\n' "\$PREFIX" > samplesheet.csv
    SEX_LINE=''
    MOUNT=''
    if [ -n "${sex_code}" ]; then
        mkdir -p sex
        VCF_ID=\$(gzip -dc "\$VCF" | awk -F'\\t' '/^#CHROM/ {id = \$10} END {print id}')
        printf '#IID\\tSEX\\n%s\\t%s\\n' "\$VCF_ID" ${sex_code} > sex/sex.tsv
        SEX_LINE="    withName: 'PLINK2_VCF' { ext.args = '--new-id-max-allele-len 100 missing --update-sex \\"\$PWD/sex/sex.tsv\\"' }"
        MOUNT=" -v \\"\$PWD/sex:\$PWD/sex:ro\\""
    fi
    # The images of versions.env for pgsc_calc's process labels, and no network for its containers.
    printf '%s\\n' 'process {' ${labels} "\$SEX_LINE" '}' "docker.runOptions = '-u \$(id -u):\$(id -g) --network none\$MOUNT'" > images.config
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
    nextflow -log pgsc_calc.nextflow.log run "${pipeline}" \\
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
        # As step 25: every score under the minimum overlap, or no score variant at all.
        if grep -q 'All scores fail to meet match threshold' pgsc_calc/pgsc_calc*.log; then
            echo "Every score matched under pgsc_calc's minimum overlap; no score." > pgsc_calc/BELOW_THRESHOLD
        elif grep -qE 'ZeroMatchesError|No match candidates found for any scoring files' pgsc_calc/pgsc_calc*.log; then
            echo "None of the score variants is in this input; no score." > pgsc_calc/ZERO_MATCHES
        else
            echo "ERROR: pgsc_calc failed (exit \$rc); see pgsc_calc/pgsc_calc.log" >&2
            exit 1
        fi
    fi
    # Its run reports carry the time of the run in their names.
    rm -rf work .nextflow pgsc_calc/pipeline_info pgsc_calc-src

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
    tuple val(meta), path(results), val(input_kind), path(scores, stageAs: 'pgs')
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
    FLAGS=()
    if [ -f ${results}/ZERO_MATCHES ]; then FLAGS=(--zero-matches); fi
    if [ -f ${results}/BELOW_THRESHOLD ]; then FLAGS=(--below-threshold ${results}/pgsc_calc.log); fi
    collect_summary.py prs-table --sample ${meta.id} --results ${results} --sampleset sample \\
        --scores pgs --input-kind ${input_kind} \${FLAGS[@]+"\${FLAGS[@]}"} ${anc} --out ${meta.id}_prs_summary.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'Condition\\tPGS_ID\\tScore_SUM\\tVariants_Matched\\tVariants_Total\\tMatched_Pct\\tPercentile\\tAncestry_Group\\tInput\\n' > ${meta.id}_prs_summary.tsv
    if [ -n "${panel_name}" ]; then printf 'key\\tvalue\\nsample\\t%s\\n' ${meta.id} > ${meta.id}_ancestry.tsv; fi
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
