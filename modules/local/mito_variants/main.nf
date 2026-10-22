/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MITO_VARIANTS — Mitochondrial variant calling with heteroplasmy detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uses GATK Mutect2 in mitochondrial mode to call variants on chrM,
    including low-frequency heteroplasmic variants (AF < 0.95).

    Four steps:
    1. Extract chrM reads from the BAM (GATK PrintReads; step 20 uses
       samtools view for the same reads)
    2. Run Mutect2 --mitochondria-mode
    3. Filter variants with FilterMutectCalls --mitochondria-mode
    4. Mark possible NuMTs (nuclear copies of chrM) with NuMTFilterTool at
       the median autosomal depth, read from MOSDEPTH's summary and
       global distribution with the awk of scripts/20-mtoolbox.sh. Without
       mosdepth in --tools the depth is 0 and the filter marks nothing, as
       step 20 does without step 16b's output.

    MITO_VARIANTS_ANNOTATE, after it, annotates each call with its ClinVar
    record and gnomAD v3.1's chrM sites (bin/collect_summary.py
    mito-annotate), as step 20's last stage does.

    Equivalent to: scripts/20-mtoolbox.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MITO_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/mito" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(mosdepth_summary), path(mosdepth_dist)  // [] [] without mosdepth
    path(reference)
    path(reference_fai)
    path(reference_dict)

    output:
    tuple val(meta), path("*_chrM_filtered.vcf.gz"),     emit: mito_vcf
    tuple val(meta), path("*_chrM_filtered.vcf.gz.tbi"), emit: mito_vcf_index, optional: true
    tuple val(meta), path("*_chrM_mutect2.vcf.gz.stats"), emit: stats
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // With mosdepth in --tools the workflow joins its output in; a join that
    // lost it would run the NuMT filter at depth 0 and mark nothing, silently.
    if ((params.tools ?: '').toString().split(',').collect { it.trim() }.contains('mosdepth') && !mosdepth_summary) {
        error "MITO_VARIANTS got no mosdepth output for '${meta.id}' although mosdepth is in --tools: the NuMT filter would run at depth 0"
    }
    def depth_files = mosdepth_summary ? "${mosdepth_summary} ${mosdepth_dist}" : ''
    """
    # Step 1: Extract chrM reads
    gatk PrintReads \\
        -I ${bam} \\
        -L chrM \\
        -O ${prefix}_chrM.bam

    # Step 2: Run Mutect2 in mitochondrial mode
    gatk Mutect2 \\
        -R ${reference} \\
        -I ${prefix}_chrM.bam \\
        -L chrM \\
        --mitochondria-mode \\
        --max-mnp-distance 0 \\
        -O ${prefix}_chrM_mutect2.vcf.gz

    # Step 3: Filter variants
    gatk FilterMutectCalls \\
        -R ${reference} \\
        -V ${prefix}_chrM_mutect2.vcf.gz \\
        --mitochondria-mode \\
        -O ${prefix}_chrM_mutect2_filtered.vcf.gz

    # Step 4: Mark possible NuMTs. The median autosomal depth is the depth at
    # which half of the chr1-22 bases are covered at least that deep, from
    # mosdepth's per-chromosome distribution (scripts/20-mtoolbox.sh, step 5).
    DEPTH_FILES="${depth_files}"
    if [ -n "\${DEPTH_FILES}" ]; then
        AUTOSOMAL_COVERAGE=\$(awk -F'\\t' '
            FNR == NR { if (\$1 ~ /^chr[0-9]+\$/) len[\$1] = \$2; next }
            (\$1 in len) { at_least[\$2] += len[\$1] * \$3 }
            END {
              for (c in len) total += len[c]
              best = 0
              if (total > 0) for (k in at_least) if (at_least[k] / total >= 0.5 && k + 0 > best) best = k + 0
              print best
            }' \${DEPTH_FILES})
        echo "Median autosomal coverage: \${AUTOSOMAL_COVERAGE} (from \${DEPTH_FILES})"
    else
        AUTOSOMAL_COVERAGE=0
        echo "mosdepth is not in --tools: NuMTFilterTool runs at depth 0 and marks nothing"
    fi
    gatk NuMTFilterTool \\
        -R ${reference} \\
        -V ${prefix}_chrM_mutect2_filtered.vcf.gz \\
        --autosomal-coverage \${AUTOSOMAL_COVERAGE} \\
        -O ${prefix}_chrM_filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // With mosdepth in --tools the workflow joins its output in; a join that
    // lost it would run the NuMT filter at depth 0 and mark nothing, silently.
    if ((params.tools ?: '').toString().split(',').collect { it.trim() }.contains('mosdepth') && !mosdepth_summary) {
        error "MITO_VARIANTS got no mosdepth output for '${meta.id}' although mosdepth is in --tools: the NuMT filter would run at depth 0"
    }
    """
    touch ${prefix}_chrM_filtered.vcf.gz
    touch ${prefix}_chrM_mutect2.vcf.gz.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}


/*
    MITO_VARIANTS_ANNOTATE: each ALT of the filtered chrM calls with its
    ClinVar record (CLNSIG, CLNREVSTAT, CLNDN of the same allele, from
    --mito_clinvar, the whole ClinVar VCF with chr names and its .tbi beside
    it) and gnomAD v3.1's chrM sites (AF_hom, AF_het, max_hl, from
    --gnomad_mito, the file GNOMAD_MITO_URL in versions.env names), into
    <id>_chrM_annotated.tsv. A file not given is named in the table's header
    and its columns are '.', as step 20 does without it.

    Equivalent to: stage 6 of scripts/20-mtoolbox.sh
*/
process MITO_VARIANTS_ANNOTATE {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/mito" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)
    path(clinvar)        // [] without --mito_clinvar
    path(clinvar_index)  // staged beside it
    path(gnomad)         // [] without --gnomad_mito

    output:
    tuple val(meta), path("${meta.id}_chrM_annotated.tsv"), emit: annotated
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cv = clinvar ? "--clinvar ${clinvar}" : ''
    def gn = gnomad ? "--gnomad ${gnomad}" : ''
    """
    collect_summary.py mito-annotate --vcf ${vcf} ${cv} ${gn} --out ${meta.id}_chrM_annotated.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'chrom\\tpos\\tref\\talt\\tfilter\\taf\\tclinvar_id\\tclinvar_significance\\tclinvar_review\\tclinvar_condition\\tgnomad_af_hom\\tgnomad_af_het\\tgnomad_max_hl\\n' > ${meta.id}_chrM_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
