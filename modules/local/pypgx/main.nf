/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PyPGx — Comprehensive pharmacogenomic star allele calling with SV detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Three-phase process from BAM:
    1. Prepare depth-of-coverage and control statistics for SV genes
    2. Call the bam genes of assets/pgx_genes.tsv (the 13 genes pypgx 0.26.0
       has a copy-number model for) with SV detection
    3. Call its vcf genes from the VCF

    CYP2D6 copy number comes from read depth: the depth CYP2D6_DEPTH measured
    is judged first (bin/cyp2d6_depth_check.py), and when the reads there are
    multi-mapped the summary's CYP2D6 row says Indeterminate. Beside the
    verdict, the check gives the MAPQ >= 20 and gene-body (without REP6)
    shares it counts from the BAM.

    The summary is bin/pgx_outside_calls.py pypgx-summary, as in the bash
    step: the defining variants of each call looked up in the VCF, a call
    resting on one below the GQ floor (task.ext.gq_floor, default 20) marked
    'low-confidence: check the reads', and a gene the table marks
    callable_without_alt=no (GSTT1) 'not callable on this reference' on a BAM
    without its ALT contig.

    Equivalent to: scripts/32-pypgx.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PYPGX {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/pypgx" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(vcf_index), path(depth_q0), path(depth_q1)
    path(reference)
    path(reference_fai)
    path(pypgx_bundle)

    output:
    tuple val(meta), path("${meta.id}_pypgx_results"), emit: results
    tuple val(meta), path("${meta.id}_pypgx_summary.tsv"), emit: summary
    tuple val(meta), path("${meta.id}_cyp2d6_depth_check.tsv"), emit: depth_check
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The gene table scripts/32-pypgx.sh reads: gene, mode (bam|vcf),
    // callable_without_alt (yes|no), alt_contig, second_caller,
    // pharmcat_outside_call. The lists land in this script, so a changed
    // table changes the task hash.
    def pgx_rows = file("${projectDir}/assets/pgx_genes.tsv", checkIfExists: true).readLines()
        .findAll { line -> line && !line.startsWith('#') && !line.startsWith('gene\t') }
        .collect { line -> line.split('\t') as List }
    def bad_rows = pgx_rows.findAll { r -> r.size() != 6 || !(r[0] ==~ /[A-Za-z0-9-]+/) || !(r[1] in ['bam', 'vcf']) || !(r[2] in ['yes', 'no']) || (r[2] == 'no' && !(r[3] ==~ /[A-Za-z0-9_.-]+/)) || (r[2] == 'no' && r[3] == '.') }
    if (pgx_rows.isEmpty()) {
        error "assets/pgx_genes.tsv lists no gene"
    }
    if (bad_rows || pgx_rows.collect { r -> r[0] }.unique().size() != pgx_rows.size()) {
        error "assets/pgx_genes.tsv has malformed rows or a gene twice: ${bad_rows.collect { r -> r.join(' ') }}"
    }
    def bam_genes = pgx_rows.findAll { r -> r[1] == 'bam' }.collect { r -> r[0] }.join(' ')
    def vcf_genes = pgx_rows.findAll { r -> r[1] == 'vcf' }.collect { r -> r[0] }.join(' ')
    // "GENE:CONTIG" for each bam gene that needs an ALT contig
    def alt_needs = pgx_rows.findAll { r -> r[1] == 'bam' && r[2] == 'no' }.collect { r -> "${r[0]}:${r[3]}" }.join(' ')
    def gq_floor = task.ext.gq_floor != null ? task.ext.gq_floor : 20
    """
    OUTBASE="${meta.id}_pypgx_results"
    mkdir -p "\$OUTBASE"

    # CYP2D6 depth check (CYP2D6_DEPTH measured it), with the MAPQ >= 20 and
    # gene-body shares counted here from the BAM, as scripts/32-pypgx.sh does
    MEASURE=""
    if cyp2d6_depth_check.py measure --bam ${bam} --out cyp2d6_mapq_shares.tsv && [ -s cyp2d6_mapq_shares.tsv ]; then
        MEASURE="--measure cyp2d6_mapq_shares.tsv"
    else
        echo "WARNING: the MAPQ >= 20 and gene-body depth at CYP2D6 was not measured; the check says 'not measured'"
    fi
    cyp2d6_depth_check.py check --all ${depth_q0} --mapq1 ${depth_q1} \$MEASURE --out ${meta.id}_cyp2d6_depth_check.tsv
    DEPTH_STATUS=\$(awk -F'\\t' '\$1 == "status" {print \$2}' ${meta.id}_cyp2d6_depth_check.tsv)

    # pypgx reads its bundle from ~/pypgx-bundle. HOME is the task directory,
    # so nothing is written inside the image (a read-only Singularity image
    # works) and the bundle is linked there unless it is already staged
    # under that name.
    export HOME="\$PWD"
    if [ -d "${pypgx_bundle}" ] && [ "${pypgx_bundle}" != "pypgx-bundle" ]; then
        ln -sfn "\$PWD/${pypgx_bundle}" "\$HOME/pypgx-bundle"
    fi

    DOC="\$OUTBASE/depth_of_coverage.zip"
    CTRL="\$OUTBASE/control_statistics.zip"
    FAILED=""
    SUCCEEDED=0
    BAM_GENES="${bam_genes}"

    # A gene the table marks callable_without_alt=no lies on an ALT contig of
    # GRCh38 (GSTT1 on chr22_KI270879v1_alt). A BAM aligned to a reference
    # without ALT contigs (the default) has no such contig, and depth
    # preparation then fails for every SV gene. Leave such a gene out in that
    # case: the summary says it is not callable on this reference, as
    # scripts/32-pypgx.sh does.
    NOT_CALLABLE=""
    for NEED in ${alt_needs}; do
      G=\${NEED%%:*} C=\${NEED#*:}
      rc=0
      python3 -c "import pysam, sys; sys.exit(0 if sys.argv[2] in pysam.AlignmentFile(sys.argv[1]).references else 3)" ${bam} "\$C" || rc=\$?
      if [ "\$rc" = 3 ]; then
        echo "NOTICE: the BAM has no \$C contig (reference without ALT contigs); \$G is not callable on this reference and is skipped"
        BAM_GENES=\$(echo "\$BAM_GENES" | tr ' ' '\\n' | grep -vxF "\$G" | paste -sd ' ' - || true)
        NOT_CALLABLE="\$NOT_CALLABLE --not-callable \$G"
      elif [ "\$rc" != 0 ]; then
        echo "ERROR: could not read the contigs of ${bam} (exit \$rc)" >&2
        exit 1
      fi
    done

    # Phase 1: Prepare depth of coverage for SV genes
    echo "--- Preparing depth of coverage for SV genes: \${BAM_GENES} ---"
    DOC_OK=true
    if [ -z "\$BAM_GENES" ]; then
      DOC_OK=false
    elif ! pypgx prepare-depth-of-coverage \\
      "\$DOC" ${bam} --assembly GRCh38 --genes \$BAM_GENES 2>&1; then
      echo "ERROR: prepare-depth-of-coverage failed"
      for GENE in \$BAM_GENES; do FAILED="\${FAILED} \${GENE}"; done
      DOC_OK=false
    fi

    # Phase 2: Compute control statistics (VDR)
    if [ "\$DOC_OK" = true ]; then
      echo "--- Computing control statistics (VDR) ---"
      if ! pypgx compute-control-statistics \\
        VDR "\$CTRL" ${bam} --assembly GRCh38 2>&1; then
        echo "WARNING: compute-control-statistics failed"
        CTRL=""
      fi
    fi

    # Phase 3a: BAM-based genes with SV detection
    if [ "\$DOC_OK" = true ]; then
      for GENE in \$BAM_GENES; do
        echo "--- Calling \${GENE} (BAM + VCF) ---"
        EXTRA=""
        [ -f "\$CTRL" ] && EXTRA="--control-statistics \$CTRL"
        pypgx run-ngs-pipeline "\$GENE" "\$OUTBASE/\${GENE}" \\
          --variants ${vcf} \\
          --depth-of-coverage "\$DOC" \\
          --assembly GRCh38 \\
          --force \\
          \$EXTRA 2>&1 \\
          && SUCCEEDED=\$((SUCCEEDED + 1)) \\
          || { echo "WARNING: \${GENE} failed"; FAILED="\${FAILED} \${GENE}"; }
      done
    fi

    # Phase 3b: VCF-based genes
    for GENE in ${vcf_genes}; do
      echo "--- Calling \${GENE} (VCF-based) ---"
      pypgx run-ngs-pipeline "\$GENE" "\$OUTBASE/\${GENE}" \\
        --variants ${vcf} \\
        --assembly GRCh38 \\
        --force 2>&1 \\
        && SUCCEEDED=\$((SUCCEEDED + 1)) \\
        || { echo "WARNING: \${GENE} failed"; FAILED="\${FAILED} \${GENE}"; }
    done

    echo "pypgx pipeline: \${SUCCEEDED} genes succeeded"
    [ -n "\$FAILED" ] && echo "Failed genes:\${FAILED}"
    [ "\$SUCCEEDED" -gt 0 ] || exit 1

    # Consolidate results into the summary TSV, as scripts/32-pypgx.sh does
    pgx_outside_calls.py pypgx-summary --outbase "\$OUTBASE" --bam-genes "${bam_genes}" --vcf-genes "${vcf_genes}" \\
        --vcf ${vcf} --gq-floor ${gq_floor} \$NOT_CALLABLE --out ${meta.id}_pypgx_summary.tsv

    # A CYP2D6 call from multi-mapped depth is not a call: the row says so
    # (the call stays in CYP2D6/results.zip), as scripts/32-pypgx.sh does.
    if [ "\$DEPTH_STATUS" != ok ]; then
        awk -F'\\t' -v OFS='\\t' '\$1 == "CYP2D6" {\$2 = "Indeterminate"; \$3 = "Indeterminate (CYP2D6 depth check)"; \$4 = "N/A (CYP2D6 depth check)"; \$6 = "no call (CYP2D6 depth check)"} {print}' \\
            ${meta.id}_pypgx_summary.tsv > summary.tmp
        mv summary.tmp ${meta.id}_pypgx_summary.tsv
        echo "WARNING: \$(awk -F'\\t' '\$1 == "message" {print \$2}' ${meta.id}_cyp2d6_depth_check.tsv)"
    fi

    # printf, not a here-document: it writes the same file whether or not
    # Nextflow strips this script's indent.
    printf '"%s":\\n    pypgx: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    mkdir -p ${meta.id}_pypgx_results
    printf 'Gene\\tDiplotype\\tPhenotype\\tCNV_call\\tSource\\tStatus\\tDefining_variants\\tMin_GQ\\tWeakest_variant\\tDefining_filter\\tQuality\\n' > ${meta.id}_pypgx_summary.tsv
    printf 'metric\\tvalue\\nstatus\\tok\\nmessage\\tstub\\n' > ${meta.id}_cyp2d6_depth_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
