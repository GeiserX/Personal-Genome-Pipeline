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
    multi-mapped the summary's CYP2D6 row says Indeterminate.

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
    """
    OUTBASE="${meta.id}_pypgx_results"
    mkdir -p "\$OUTBASE"

    # CYP2D6 depth check (CYP2D6_DEPTH measured it)
    cyp2d6_depth_check.py check --all ${depth_q0} --mapq1 ${depth_q1} --out ${meta.id}_cyp2d6_depth_check.tsv
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
    # case and say so, as scripts/32-pypgx.sh does.
    if ! python3 -c "import pysam, sys; sys.exit(0 if 'chr22_KI270879v1_alt' in pysam.AlignmentFile(sys.argv[1]).references else 1)" ${bam}; then
      echo "NOTICE: the BAM has no chr22_KI270879v1_alt contig (reference without ALT contigs); GSTT1 cannot be called from depth and is skipped"
      BAM_GENES=\$(echo "\$BAM_GENES" | tr " " "\\n" | grep -vx GSTT1 | tr "\\n" " ")
      FAILED="\${FAILED} GSTT1"
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

    # Consolidate results into summary TSV
    python3 -c "
import os, sys, subprocess, csv

outbase = '\$OUTBASE'
bam_genes = '${bam_genes}'.split()
vcf_genes = '${vcf_genes}'.split()
all_genes = bam_genes + vcf_genes

summary_path = '${meta.id}_pypgx_summary.tsv'
rows = []

for gene in all_genes:
    results_zip = f'{outbase}/{gene}/results.zip'
    if not os.path.isfile(results_zip):
        rows.append([gene, 'FAILED', 'N/A', 'N/A', 'BAM' if gene in bam_genes else 'VCF'])
        continue

    diplotype = 'N/A'
    phenotype = 'N/A'
    try:
        out = subprocess.run(
            ['pypgx', 'print-data', results_zip],
            capture_output=True, text=True
        )
        if out.returncode == 0:
            lines = out.stdout.rstrip().split('\\n')
            if len(lines) >= 2:
                headers = lines[0].split('\\t')
                values = lines[1].split('\\t')
                if 'Genotype' in headers:
                    idx = headers.index('Genotype')
                    if idx < len(values):
                        diplotype = values[idx]
                if 'Phenotype' in headers:
                    idx = headers.index('Phenotype')
                    if idx < len(values):
                        phenotype = values[idx]
    except Exception as e:
        print(f'WARNING: Error extracting {gene}: {e}', file=sys.stderr)

    source = 'BAM' if gene in bam_genes else 'VCF'
    if gene in bam_genes:
        sv_detected = 'Yes' if any(x in (diplotype or '') for x in ['DEL', 'DUP', 'x2', 'x3', '*5']) else 'No'
    else:
        sv_detected = 'N/A'
    rows.append([gene, diplotype, phenotype, sv_detected, source])

with open(summary_path, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\\t')
    w.writerow(['Gene', 'Diplotype', 'Phenotype', 'SV_detected', 'Source'])
    w.writerows(rows)

called = sum(1 for r in rows if r[1] != 'FAILED')
print(f'Summary: {called}/{len(rows)} genes called')
"

    # A CYP2D6 call from multi-mapped depth is not a call: the row says so
    # (the call stays in CYP2D6/results.zip), as scripts/32-pypgx.sh does.
    if [ "\$DEPTH_STATUS" != ok ]; then
        awk -F'\\t' -v OFS='\\t' '\$1 == "CYP2D6" {\$2 = "Indeterminate"; \$3 = "Indeterminate (CYP2D6 depth check)"} {print}' \\
            ${meta.id}_pypgx_summary.tsv > summary.tmp
        mv summary.tmp ${meta.id}_pypgx_summary.tsv
        echo "WARNING: \$(awk -F'\\t' '\$1 == "message" {print \$2}' ${meta.id}_cyp2d6_depth_check.tsv)"
    fi

    # printf, not a here-document: the column-0 lines above stop Nextflow
    # from stripping this script's indent, and an indented END_VERSIONS
    # would not end a <<- here-document.
    printf '"%s":\\n    pypgx: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    mkdir -p ${meta.id}_pypgx_results
    printf 'Gene\\tDiplotype\\tPhenotype\\tSV_detected\\tSource\\n' > ${meta.id}_pypgx_summary.tsv
    printf 'metric\\tvalue\\nstatus\\tok\\nmessage\\tstub\\n' > ${meta.id}_cyp2d6_depth_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
