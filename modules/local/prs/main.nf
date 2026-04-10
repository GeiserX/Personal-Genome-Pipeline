/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRS — Polygenic Risk Scores via plink2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calculates polygenic risk scores from VCF using PGS Catalog scoring files.
    Converts VCF to plink2 binary format, then runs --score for each scoring file
    found in the scoring directory.

    NOTE: Raw PRS from a single sample are NOT directly interpretable without a
    population reference distribution. Treat as exploratory, not clinical.

    Equivalent to: scripts/25-prs.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PRS {
    tag "$meta.id"
    label 'process_medium'

    container 'pgscatalog/plink2:2.00a5.10'

    publishDir "${params.outdir}/${meta.id}/prs", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(scoring_dir)

    output:
    tuple val(meta), path("*.sscore"),           emit: scores, optional: true
    tuple val(meta), path("*_prs_summary.tsv"),  emit: summary
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Step 1: Convert VCF to plink2 binary format
    plink2 \\
        --vcf ${vcf} \\
        --make-pgen \\
        --out ${meta.id} \\
        --threads ${task.cpus} \\
        --memory \$(( ${task.memory.toMega()} - 500 )) \\
        --set-all-var-ids '@:#' \\
        --new-id-max-allele-len 100 \\
        --chr 1-22 \\
        --allow-extra-chr \\
        --output-chr chrM

    # Step 2: Score each PGS file in scoring directory
    echo -e "Condition\\tPGS_ID\\tScore\\tVariants_Used\\tVariants_Total" > ${meta.id}_prs_summary.tsv

    for SCORE_FILE in ${scoring_dir}/*.txt.gz ${scoring_dir}/*.txt; do
        [ -f "\${SCORE_FILE}" ] || continue
        PGS_ID=\$(basename "\${SCORE_FILE}" | sed 's/\\(_hmPOS_GRCh38\\)\\?.txt\\(.gz\\)\\?\$//')

        # Format scoring file: extract chr:pos, effect_allele, effect_weight
        FORMATTED="\${PGS_ID}_formatted.tsv"
        (zcat "\${SCORE_FILE}" 2>/dev/null || cat "\${SCORE_FILE}") | grep -v "^#" | \\
            awk -F'\\t' 'NR==1 {
                for(i=1;i<=NF;i++) {
                    if(\$i=="chr_name") chr_col=i;
                    if(\$i=="chr_position") pos_col=i;
                    if(\$i=="effect_allele") ea_col=i;
                    if(\$i=="effect_weight") ew_col=i;
                    if(\$i=="hm_chr") chr_col=i;
                    if(\$i=="hm_pos") pos_col=i;
                }
                next
            }
            chr_col && pos_col && ea_col && ew_col {
                chr=\$chr_col; pos=\$pos_col; ea=\$ea_col; ew=\$ew_col;
                if(chr!="" && pos!="" && ea!="" && ew!="") {
                    if(chr !~ /^chr/) chr="chr"chr;
                    key=chr":"pos"\\t"ea;
                    if(!(key in seen)) { seen[key]=1; printf "%s:%s\\t%s\\t%s\\n", chr, pos, ea, ew; }
                }
            }' > "\${FORMATTED}" 2>/dev/null || true

        TOTAL_VARS=\$(wc -l < "\${FORMATTED}" 2>/dev/null || echo 0)
        [ "\${TOTAL_VARS}" -eq 0 ] && continue

        # Run plink2 --score
        plink2 \\
            --pfile ${meta.id} \\
            --score "\${FORMATTED}" 1 2 3 \\
                ignore-dup-ids \\
                no-mean-imputation \\
            --out "\${PGS_ID}" \\
            --threads ${task.cpus} \\
            --memory \$(( ${task.memory.toMega()} - 500 )) \\
            --allow-extra-chr 2>/dev/null || true

        if [ -f "\${PGS_ID}.sscore" ]; then
            SCORE=\$(awk 'NR==2 {print \$NF}' "\${PGS_ID}.sscore" 2>/dev/null || echo "N/A")
            USED_VARS=\$(awk 'NR==2 {print \$(NF-1)}' "\${PGS_ID}.sscore" 2>/dev/null || echo "N/A")
            echo -e "\${PGS_ID}\\t\${PGS_ID}\\t\${SCORE}\\t\${USED_VARS}\\t\${TOTAL_VARS}" >> ${meta.id}_prs_summary.tsv
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1 | awk '{print \$2}' || echo '2.00a5.10')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_prs_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1 | awk '{print \$2}' || echo '2.00a5.10')
    END_VERSIONS
    """
}
