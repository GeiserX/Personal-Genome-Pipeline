/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PyPGx — Comprehensive pharmacogenomic star allele calling with SV detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Three-phase process from BAM:
    1. Prepare depth-of-coverage and control statistics for SV genes
    2. Call BAM-based genes (CYP2D6, CYP2A6, GSTM1, GSTT1) with SV detection
    3. Call VCF-based genes (~19 additional pharmacogenes)

    CYP2D6 copy number comes from read depth: the depth CYP2D6_DEPTH measured
    is judged first (bin/cyp2d6_depth_check.py), and when the reads there are
    multi-mapped the summary's CYP2D6 row says Indeterminate. Beside the
    verdict, the check gives the MAPQ >= 20 and gene-body (without REP6)
    shares it counts from the BAM.

    The summary is bin/pgx_outside_calls.py pypgx-summary, as in the bash
    step: the defining variants of each call looked up in the VCF, a call
    resting on one below the GQ floor (task.ext.gq_floor, default 20) marked
    'low-confidence: check the reads', and GSTT1 'not callable on this
    reference' on a BAM without its ALT contig.

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
    def bam_genes = "CYP2D6 CYP2A6 GSTM1 GSTT1"
    def vcf_genes = "CYP1A2 CYP2B6 CYP2C9 CYP2C19 CYP3A4 CYP3A5 CYP4F2 DPYD TPMT NUDT15 UGT1A1 SLCO1B1 VKORC1 NAT2 COMT MTHFR ABCB1 G6PD IFNL3"
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

    # GSTT1 lies on chr22_KI270879v1_alt in GRCh38. A BAM aligned to a
    # reference without ALT contigs (the default) has no such contig, and
    # depth preparation then fails for every SV gene. Leave GSTT1 out in that
    # case: the summary says it is not callable on this reference, as
    # scripts/32-pypgx.sh does.
    NOT_CALLABLE=""
    rc=0
    python3 -c "import pysam, sys; sys.exit(0 if 'chr22_KI270879v1_alt' in pysam.AlignmentFile(sys.argv[1]).references else 3)" ${bam} || rc=\$?
    if [ "\$rc" = 3 ]; then
      echo "NOTICE: the BAM has no chr22_KI270879v1_alt contig (reference without ALT contigs); GSTT1 is not callable on this reference and is skipped"
      BAM_GENES="CYP2D6 CYP2A6 GSTM1"
      NOT_CALLABLE="--not-callable GSTT1"
    elif [ "\$rc" != 0 ]; then
      echo "ERROR: could not read the contigs of ${bam} (exit \$rc)" >&2
      exit 1
    fi

    # Phase 1: Prepare depth of coverage for SV genes
    echo "--- Preparing depth of coverage for SV genes: \${BAM_GENES} ---"
    DOC_OK=true
    if ! pypgx prepare-depth-of-coverage \\
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
