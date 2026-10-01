/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SLIVAR — Variant prioritization and compound heterozygote detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prioritizes clinically interesting variants using tiered filters (rare HIGH,
    rare MODERATE + deleterious predictors, ClinVar pathogenic) and detects
    compound heterozygote candidates. Optionally annotates results with gnomAD
    gene constraint metrics (LOEUF, pLI).

    Equivalent to: scripts/31-slivar.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SLIVAR {
    tag "$meta.id"
    label 'process_medium'

    container 'staphb/bcftools:1.21'

    publishDir { "${params.outdir}/${meta.id}/slivar" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(gnomad_constraint)
    path(slivar_bin)

    output:
    tuple val(meta), path("*_prioritized.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_prioritized.vcf.gz.tbi"), emit: vcf_index
    tuple val(meta), path("*_compound_hets.vcf.gz"),   emit: compound_het_vcf
    tuple val(meta), path("*_slivar_summary.tsv"),     emit: summary_tsv
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_constraint = gnomad_constraint ? true : false
    """
    # --- Make the staged slivar binary executable ---
    chmod +x ${slivar_bin}

    # --- Generate PED file for single sample ---
    SAMPLE_NAME=\$(bcftools query -l ${vcf} | head -1)
    echo -e "\${SAMPLE_NAME}\\t\${SAMPLE_NAME}\\t0\\t0\\t0\\t-9" > ${meta.id}.ped

    # --- Which CSQ fields does the input carry? ---
    # Tested explicitly (as CLINICAL_FILTER does) instead of falling back with
    # `A | B || C | D`, which hid real failures behind the fallback branch.
    bcftools +split-vep -l ${vcf} > csq_fields.txt
    HAS_GNOMAD=0
    HAS_CLINSIG=0
    awk '\$2 == "gnomADe_AF" { f = 1 } END { exit !f }' csq_fields.txt && HAS_GNOMAD=1
    awk '\$2 == "CLIN_SIG" { f = 1 } END { exit !f }' csq_fields.txt && HAS_CLINSIG=1
    if [ "\${HAS_GNOMAD}" -eq 1 ]; then
        RARE_COLS="IMPACT,gnomADe_AF"
        RARE_AND=' && (gnomADe_AF<0.01 || gnomADe_AF=".")'
    else
        echo "WARNING: gnomADe_AF not in the CSQ fields; rare tiers are not filtered by frequency." >&2
        RARE_COLS="IMPACT"
        RARE_AND=""
    fi

    # --- Filter 1: rare_high (PASS + HIGH impact + gnomAD AF < 1%) ---
    bcftools view -f PASS ${vcf} | \\
        bcftools +split-vep - -c "\${RARE_COLS}" -s worst \\
            -i "IMPACT=\\"HIGH\\"\${RARE_AND}" \\
            -Oz -o ${meta.id}_rare_high.vcf.gz
    bcftools index -t ${meta.id}_rare_high.vcf.gz

    # --- Filter 2: rare_moderate_deleterious (two-pass: CSQ then INFO predictors) ---
    # Pass 1: extract rare MODERATE via split-vep (CSQ fields)
    bcftools view -f PASS ${vcf} | \\
        bcftools +split-vep - -c "\${RARE_COLS}" -s worst \\
            -i "IMPACT=\\"MODERATE\\"\${RARE_AND}" \\
            -Oz -o ${meta.id}_rare_moderate_all.vcf.gz
    bcftools index -t ${meta.id}_rare_moderate_all.vcf.gz

    # Pass 2: gate on deleteriousness predictors in INFO fields (added by vcfanno)
    # Build predictor expression from available INFO fields in the VCF header
    PREDICTOR_PARTS=""
    INFO_HEADER=\$(bcftools view -h ${meta.id}_rare_moderate_all.vcf.gz | grep '^##INFO' || true)

    echo "\$INFO_HEADER" | grep -q 'ID=CADD_PHRED,' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/CADD_PHRED>=20"
    echo "\$INFO_HEADER" | grep -q 'ID=CADD_PHRED_indel,' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/CADD_PHRED_indel>=20"
    echo "\$INFO_HEADER" | grep -q 'ID=REVEL' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/REVEL>=0.5"
    echo "\$INFO_HEADER" | grep -q 'ID=AM_class' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/AM_class=\\"likely_pathogenic\\""
    echo "\$INFO_HEADER" | grep -q 'ID=SpliceAI,' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/SpliceAI!=\\".\\""
    echo "\$INFO_HEADER" | grep -q 'ID=SpliceAI_indel,' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/SpliceAI_indel!=\\".\\""

    if [ -n "\${PREDICTOR_PARTS}" ]; then
        echo "  Filtering MODERATE variants with predictors: \${PREDICTOR_PARTS}"
        bcftools view -i "\${PREDICTOR_PARTS}" \\
            ${meta.id}_rare_moderate_all.vcf.gz \\
            -Oz -o ${meta.id}_rare_moderate_del.vcf.gz
    else
        echo "  WARNING: No vcfanno predictor annotations found — including all rare MODERATE variants."
        echo "  Run vcfanno for CADD/REVEL/AlphaMissense/SpliceAI filtering."
        cp ${meta.id}_rare_moderate_all.vcf.gz ${meta.id}_rare_moderate_del.vcf.gz
    fi
    bcftools index -t ${meta.id}_rare_moderate_del.vcf.gz

    # --- Filter 3: clinvar_pathogenic ---
    MERGE_FILES="${meta.id}_rare_high.vcf.gz ${meta.id}_rare_moderate_del.vcf.gz"
    if [ "\${HAS_CLINSIG}" -eq 1 ]; then
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c CLIN_SIG \\
                -i 'CLIN_SIG~"pathogenic" && CLIN_SIG!~"conflicting"' \\
                -Oz -o ${meta.id}_clinvar_path.vcf.gz
        bcftools index -t ${meta.id}_clinvar_path.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_clinvar_path.vcf.gz"
    else
        echo "ClinVar tier skipped (CLIN_SIG not in the CSQ fields)"
    fi

    # --- Merge tiers into prioritized VCF ---
    bcftools concat -a -D \${MERGE_FILES} | \\
        bcftools sort -Oz -o ${meta.id}_prioritized.vcf.gz
    bcftools index -t ${meta.id}_prioritized.vcf.gz

    # --- Compound heterozygote detection ---
    # slivar writes to a file first so its own exit status fails the task: a
    # crash must not read as "no compound hets".
    ./${slivar_bin} compound-hets \\
        --allow-non-trios \\
        --vcf ${meta.id}_prioritized.vcf.gz \\
        --ped ${meta.id}.ped \\
        > ${meta.id}_compound_hets.vcf \\
        2> ${meta.id}_compound_hets.log
    bcftools view ${meta.id}_compound_hets.vcf -Oz -o ${meta.id}_compound_hets.vcf.gz
    rm -f ${meta.id}_compound_hets.vcf

    # --- Generate summary TSV with optional gnomAD constraint enrichment ---
    bcftools +split-vep \\
        ${meta.id}_prioritized.vcf.gz \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%IMPACT\\t%SYMBOL\\t%Consequence\\t%Existing_variation[\\t%GT]\\n' \\
        -s worst -d > ${meta.id}_variants_raw.tsv

    if [ "${has_constraint}" = "true" ]; then
        # Join with gnomAD v4.1 constraint metrics, keyed on gene symbol. Only
        # canonical transcripts count; v4.1 lists an Ensembl and a RefSeq
        # canonical row per gene, and the Ensembl (ENST) one wins.
        awk -F'\\t' -v OFS='\\t' '
            function print_header() {
                print "CHROM", "POS", "REF", "ALT", "IMPACT", "SYMBOL", "Consequence", "Existing_variation", "GT", "LOEUF", "pLI", "mis_z", "CONSTRAINED"
            }
            NR == FNR {
                if (FNR == 1) {
                    for (i = 1; i <= NF; i++) col[\$i] = i
                    split("gene canonical transcript lof.oe_ci.upper lof.pLI mis.z_score", need, " ")
                    for (k in need) if (!(need[k] in col)) {
                        print "ERROR: column " need[k] " missing from the constraint table" > "/dev/stderr"
                        bad = 1; exit 3
                    }
                    next
                }
                if (\$col["canonical"] != "true") next
                g = \$col["gene"]
                ens = (\$col["transcript"] ~ /^ENST/)
                if ((g in val) && (src[g] || !ens)) next
                l = \$col["lof.oe_ci.upper"]; p = \$col["lof.pLI"]; m = \$col["mis.z_score"]
                if (l == "NA" || l == "") l = "."
                if (p == "NA" || p == "") p = "."
                if (m == "NA" || m == "") m = "."
                val[g] = l OFS p OFS m
                loeuf[g] = l; pli[g] = p; src[g] = ens
                next
            }
            !hdr { print_header(); hdr = 1 }
            /^#/ || NF < 6 { next }
            {
                rows++
                g = \$6
                if (g != "." && g != "") with_gene++
                if (g in val) {
                    matched[g] = 1
                    c = "NO"
                    if ((loeuf[g] != "." && loeuf[g] + 0 < 0.35) || (pli[g] != "." && pli[g] + 0 > 0.9)) c = "YES"
                    print \$0, val[g], c
                } else {
                    print \$0, ".", ".", ".", "NO"
                }
            }
            END {
                if (bad) exit 3
                if (!hdr) print_header()
                n = 0
                for (g in matched) n++
                printf "gnomAD constraint: %d of %d variant rows carry a gene symbol; %d distinct genes matched\\n", with_gene, rows, n > "/dev/stderr"
                if (with_gene > 0 && n == 0) {
                    print "ERROR: --gnomad_constraint is set but no gene in the variants matched it; check that the file is the gnomAD v4.1 constraint table" > "/dev/stderr"
                    exit 4
                }
            }
        ' ${gnomad_constraint} ${meta.id}_variants_raw.tsv > ${meta.id}_slivar_summary.tsv
    else
        {
            echo -e "CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT"
            cat ${meta.id}_variants_raw.tsv
        } > ${meta.id}_slivar_summary.tsv
    fi

    # Clean up intermediate files
    rm -f ${meta.id}_variants_raw.tsv ${meta.id}_rare_high.vcf.gz* ${meta.id}_rare_moderate_del.vcf.gz* \\
          ${meta.id}_rare_moderate_all.vcf.gz* ${meta.id}_clinvar_path.vcf.gz* csq_fields.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(./${slivar_bin} 2>&1 | grep -oP '[0-9]+\\.[0-9.]+' | head -1)
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_prioritized.vcf.gz
    touch ${meta.id}_prioritized.vcf.gz.tbi
    touch ${meta.id}_compound_hets.vcf.gz
    printf 'CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT\\n' > ${meta.id}_slivar_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: 0.3.4
        bcftools: 1.21
    END_VERSIONS
    """
}
