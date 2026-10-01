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

    publishDir { "${params.outdir}/${meta.id}/prs" }, mode: params.publish_dir_mode

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

    # Step 2: Score each PGS file in scoring directory.
    # Only GRCh38-harmonised PGS Catalog files are accepted (#HmPOS_build=GRCh38);
    # an author-reported file is often GRCh37 and would score the wrong positions.
    echo -e "Condition\\tPGS_ID\\tScore_SUM\\tVariants_Matched\\tVariants_Total" > ${meta.id}_prs_summary.tsv

    for SCORE_FILE in ${scoring_dir}/*.txt.gz ${scoring_dir}/*.txt; do
        [ -f "\${SCORE_FILE}" ] || continue
        PGS_ID=\$(basename "\${SCORE_FILE}" | sed 's/\\(_hmPOS_GRCh38\\)\\?.txt\\(.gz\\)\\?\$//')

        BUILD=\$( { gzip -cdf "\${SCORE_FILE}" 2>/dev/null || true; } | awk -F= '/^#HmPOS_build=/ {b=\$2} !/^#/ {exit} END {print b}')
        if [ "\${BUILD}" != "GRCh38" ]; then
            echo "ERROR: \${SCORE_FILE} has #HmPOS_build='\${BUILD}', expected GRCh38 (use the PGS Catalog Harmonized/<id>_hmPOS_GRCh38 file)" >&2
            exit 1
        fi

        # Format scoring file: GRCh38 hm_chr:hm_pos, effect_allele, effect_weight
        FORMATTED="\${PGS_ID}_formatted.tsv"
        gzip -cdf "\${SCORE_FILE}" | \\
            awk -F'\\t' '/^#/ {next}
            !hdr {
                for(i=1;i<=NF;i++) {
                    if(\$i=="hm_chr") chr_col=i;
                    if(\$i=="hm_pos") pos_col=i;
                    if(\$i=="effect_allele") ea_col=i;
                    if(\$i=="effect_weight") ew_col=i;
                }
                hdr=1
                next
            }
            chr_col && pos_col && ea_col && ew_col {
                chr=\$chr_col; pos=\$pos_col; ea=\$ea_col; ew=\$ew_col;
                if(chr!="" && pos!="" && ea!="" && ew!="") {
                    if(chr !~ /^chr/) chr="chr"chr;
                    key=chr":"pos"\\t"ea;
                    if(!(key in seen)) { seen[key]=1; printf "%s:%s\\t%s\\t%s\\n", chr, pos, ea, ew; }
                }
            }' > "\${FORMATTED}"

        TOTAL_VARS=\$(wc -l < "\${FORMATTED}" | tr -d ' ')
        if [ "\${TOTAL_VARS}" -eq 0 ]; then
            echo "ERROR: no GRCh38 hm_chr/hm_pos/effect_allele/effect_weight rows in \${SCORE_FILE}" >&2
            exit 1
        fi

        # Run plink2 --score. cols=+scoresums adds SCORE1_SUM (the plain weighted sum).
        rm -f "\${PGS_ID}.sscore" "\${PGS_ID}.log"
        if ! plink2 \\
            --pfile ${meta.id} \\
            --score "\${FORMATTED}" 1 2 3 \\
                ignore-dup-ids \\
                no-mean-imputation \\
                cols=+scoresums \\
            --out "\${PGS_ID}" \\
            --threads ${task.cpus} \\
            --memory \$(( ${task.memory.toMega()} - 500 )) \\
            --allow-extra-chr; then
            if grep -q 'No valid variants' "\${PGS_ID}.log" 2>/dev/null; then
                echo -e "\${PGS_ID}\\t\${PGS_ID}\\tNA\\t0\\t\${TOTAL_VARS}" >> ${meta.id}_prs_summary.tsv
                continue
            fi
            echo "ERROR: plink2 --score failed for \${PGS_ID}" >&2
            exit 1
        fi

        # Columns by name; ALLELE_CT/2 = variants matched (two alleles per autosomal site)
        ROW=\$(awk -F'\\t' 'NR==1 { for(i=1;i<=NF;i++) col[\$i]=i; next }
            NR==2 && ("SCORE1_SUM" in col) && ("ALLELE_CT" in col) { print \$col["SCORE1_SUM"] "\\t" \$col["ALLELE_CT"]/2 }' "\${PGS_ID}.sscore")
        if [ -z "\${ROW}" ]; then
            echo "ERROR: \${PGS_ID}.sscore has no SCORE1_SUM or ALLELE_CT column" >&2
            exit 1
        fi
        echo -e "\${PGS_ID}\\t\${PGS_ID}\\t\${ROW}\\t\${TOTAL_VARS}" >> ${meta.id}_prs_summary.tsv
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
