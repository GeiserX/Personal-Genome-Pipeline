/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CYRIUS — CYP2D6 star allele calling from WGS BAM (opt-in: --tools cyrius)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CYP2D6 is the hardest pharmacogene to call because of its pseudogene (CYP2D7)
    and complex structural variants (deletions, duplications, hybrids).
    Cyrius uses depth-based analysis specifically designed for CYP2D6.

    Cyrius 1.1.1 has had no release since 2021 and is under the PolyForm
    Strict licence 1.0.0 (non-commercial use only), so it is opt-in. It runs
    from the directory `scripts/setup.sh --cyrius` installs (--cyrius_install),
    hash-locked by scripts/cyrius-constraints.txt, with no network, as every
    other task. It is the second CYP2D6 caller PGX_CONSENSUS needs before a
    CYP2D6 call reaches PharmCAT.

    The depth CYP2D6_DEPTH measured is judged first
    (bin/cyp2d6_depth_check.py): when the reads there are multi-mapped, the
    call's Filter becomes CYP2D6_depth_unreliable and PGX_CONSENSUS does not
    pass it on.

    Equivalent to: scripts/21-cyrius.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CYRIUS {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/cyrius" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(depth_q0), path(depth_q1)
    path(cyrius_install)  // setup.sh --cyrius: GENOME_DIR/tools/cyrius-<version>
    path(cyrius_lock)     // scripts/cyrius-constraints.txt, which the install must come from

    output:
    tuple val(meta), path("*_cyp2d6.tsv"),                       emit: cyp2d6_results
    tuple val(meta), path("${meta.id}_cyp2d6_depth_check.tsv"),  emit: depth_check
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # The stamp setup.sh --cyrius writes: this image and this lock file.
    want="python=${task.container} lock=\$(sha256sum ${cyrius_lock} | cut -d' ' -f1)"
    if [ "\$(cat ${cyrius_install}/INSTALLED 2>/dev/null)" != "\$want" ]; then
        echo "ERROR: ${cyrius_install} is not a Cyrius install for this image and scripts/cyrius-constraints.txt" >&2
        echo "  (want '\$want'). Run scripts/setup.sh --cyrius <genome_dir> again." >&2
        exit 1
    fi
    cyp2d6_depth_check.py check --all ${depth_q0} --mapq1 ${depth_q1} --out ${meta.id}_cyp2d6_depth_check.tsv

    echo "${bam}" > manifest.txt
    PYTHONPATH="\$PWD/${cyrius_install}" python3 -m cyrius \\
        --manifest manifest.txt \\
        --genome 38 \\
        --prefix ${prefix}_cyp2d6 \\
        --outDir ./ \\
        --threads ${task.cpus}

    # Keep Cyrius's genotype for the record; the Filter says it cannot be used.
    if [ "\$(awk -F'\\t' '\$1 == "status" {print \$2}' ${meta.id}_cyp2d6_depth_check.tsv)" != ok ]; then
        awk -F'\\t' -v OFS='\\t' 'NR > 1 {\$3 = "CYP2D6_depth_unreliable"} {print}' ${prefix}_cyp2d6.tsv > cyp2d6.tmp
        mv cyp2d6.tmp ${prefix}_cyp2d6.tsv
        echo "WARNING: \$(awk -F'\\t' '\$1 == "message" {print \$2}' ${meta.id}_cyp2d6_depth_check.tsv)"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Sample\\tGenotype\\tFilter\\n${prefix}\\t*1/*1\\tPASS\\n' > ${prefix}_cyp2d6.tsv
    printf 'metric\\tvalue\\nstatus\\tok\\nmessage\\tstub\\n' > ${meta.id}_cyp2d6_depth_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
