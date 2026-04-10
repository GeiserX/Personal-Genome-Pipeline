/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINICAL_FILTER — Extract clinically interesting variants from annotated VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Produces a small VCF of variants that are:
      - Rare (gnomAD AF < 1%) AND functionally impactful (HIGH/MODERATE VEP impact)
      - OR ClinVar pathogenic/likely pathogenic
      - OR high CADD score (>= 20) for non-coding variants
      - OR high SpliceAI delta score (>= 0.2) for cryptic splice variants
      - OR high REVEL/AlphaMissense for missense variants

    Uses bcftools +split-vep for VEP CSQ fields and bcftools view -i for
    INFO-level annotations from vcfanno. Python constraint enrichment is deferred.

    Equivalent to: scripts/23-clinical-filter.sh (bcftools portion only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINICAL_FILTER {
    tag "$meta.id"
    label 'process_low'

    container 'staphb/bcftools:1.21'

    publishDir "${params.outdir}/${meta.id}/clinical", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*_clinical.vcf.gz"),         emit: vcf
    tuple val(meta), path("*_clinical.vcf.gz.tbi"),     emit: vcf_index
    tuple val(meta), path("*_clinical_summary.tsv"),    emit: summary_tsv
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # --- Detect available annotation fields ---
    HAS_GNOMAD=0
    HAS_CLINVAR=0
    HAS_CADD=0
    HAS_CADD_INDEL=0
    HAS_SPLICEAI=0
    HAS_SPLICEAI_INDEL=0
    HAS_REVEL=0
    HAS_AM=0

    VEP_FIELDS=\$(bcftools +split-vep -l ${vcf} 2>/dev/null || echo "")
    if ! echo "\${VEP_FIELDS}" | grep -q 'IMPACT'; then
        echo "ERROR: CLINICAL_FILTER requires a VEP-annotated VCF with CSQ/IMPACT field." >&2
        echo "Enable 'vep' in --tools before 'clinical_filter'." >&2
        exit 1
    fi
    echo "\${VEP_FIELDS}" | grep -q 'gnomADe_AF' && HAS_GNOMAD=1
    echo "\${VEP_FIELDS}" | grep -q 'CLIN_SIG' && HAS_CLINVAR=1

    VCF_HEADER=\$(bcftools view -h ${vcf} 2>/dev/null || echo "")
    echo "\${VCF_HEADER}" | grep -q 'ID=CADD_PHRED,' && HAS_CADD=1
    echo "\${VCF_HEADER}" | grep -q 'ID=CADD_PHRED_indel,' && HAS_CADD_INDEL=1
    echo "\${VCF_HEADER}" | grep -q 'ID=SpliceAI,' && HAS_SPLICEAI=1
    echo "\${VCF_HEADER}" | grep -q 'ID=SpliceAI_indel,' && HAS_SPLICEAI_INDEL=1
    echo "\${VCF_HEADER}" | grep -q 'ID=REVEL' && HAS_REVEL=1
    echo "\${VCF_HEADER}" | grep -q 'ID=AM_pathogenicity' && HAS_AM=1

    # --- Filter 1: HIGH impact variants ---
    bcftools view -f PASS ${vcf} | \\
        bcftools +split-vep - -c IMPACT -s worst -i 'IMPACT="HIGH"' \\
            -Oz -o ${meta.id}_high_impact.vcf.gz
    bcftools index -t ${meta.id}_high_impact.vcf.gz

    # --- Filter 2: Rare MODERATE impact variants ---
    if [ "\${HAS_GNOMAD}" -eq 1 ]; then
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c IMPACT,gnomADe_AF -s worst \\
                -i 'IMPACT="MODERATE" && (gnomADe_AF<0.01 || gnomADe_AF=".")' \\
                -Oz -o ${meta.id}_rare_moderate.vcf.gz
    else
        # Without gnomAD AF, skip rarity filter — emit header-only VCF to avoid thousands of unfiltered MODERATE variants
        echo "WARN: gnomADe_AF not found — skipping rare MODERATE tier" >&2
        bcftools view -h ${vcf} | bgzip -c > ${meta.id}_rare_moderate.vcf.gz
    fi
    bcftools index -t ${meta.id}_rare_moderate.vcf.gz

    # --- Filter 3: ClinVar pathogenic/likely pathogenic ---
    CLINVAR_FILE=""
    if [ "\${HAS_CLINVAR}" -eq 1 ]; then
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c CLIN_SIG \\
                -i 'CLIN_SIG~"pathogenic" && CLIN_SIG!~"conflicting"' \\
                -Oz -o ${meta.id}_clinvar_pathogenic.vcf.gz
        bcftools index -t ${meta.id}_clinvar_pathogenic.vcf.gz
        CLINVAR_FILE="${meta.id}_clinvar_pathogenic.vcf.gz"
    fi

    # --- Filter 4: High CADD non-coding variants ---
    CADD_FILE=""
    CADD_EXPR=""
    [ "\${HAS_CADD}" -eq 1 ] && CADD_EXPR="INFO/CADD_PHRED>=20"
    if [ "\${HAS_CADD_INDEL}" -eq 1 ]; then
        [ -n "\${CADD_EXPR}" ] && CADD_EXPR="\${CADD_EXPR} || INFO/CADD_PHRED_indel>=20" || CADD_EXPR="INFO/CADD_PHRED_indel>=20"
    fi
    if [ -n "\${CADD_EXPR}" ]; then
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c IMPACT -s worst \\
                -i "IMPACT!=\\"HIGH\\" && IMPACT!=\\"MODERATE\\" && (\${CADD_EXPR})" \\
                -Oz -o ${meta.id}_cadd_high.vcf.gz
        bcftools index -t ${meta.id}_cadd_high.vcf.gz
        CADD_FILE="${meta.id}_cadd_high.vcf.gz"
    fi

    # --- Filter 5: SpliceAI cryptic splice variants ---
    SPLICEAI_FILE=""
    SPLICEAI_PREFILTER=""
    [ "\${HAS_SPLICEAI}" -eq 1 ] && SPLICEAI_PREFILTER='INFO/SpliceAI!="."'
    if [ "\${HAS_SPLICEAI_INDEL}" -eq 1 ]; then
        [ -n "\${SPLICEAI_PREFILTER}" ] && SPLICEAI_PREFILTER="\${SPLICEAI_PREFILTER} || "'INFO/SpliceAI_indel!="."' || SPLICEAI_PREFILTER='INFO/SpliceAI_indel!="."'
    fi
    if [ -n "\${SPLICEAI_PREFILTER}" ]; then
        bcftools view -f PASS -i "\${SPLICEAI_PREFILTER}" ${vcf} | \\
            awk -F'\\t' 'BEGIN{OFS="\\t"} /^#/{print;next} {
                dominated=0
                n=split(\$8, info_arr, ";")
                for(i=1;i<=n;i++){
                    if(info_arr[i] ~ /^SpliceAI=/){
                        sub(/^SpliceAI=/,"",info_arr[i])
                        split(info_arr[i],sp,"|")
                        for(j=3;j<=6;j++) if(sp[j]+0>=0.2) dominated=1
                    }
                    if(info_arr[i] ~ /^SpliceAI_indel=/){
                        sub(/^SpliceAI_indel=/,"",info_arr[i])
                        split(info_arr[i],sp,"|")
                        for(j=3;j<=6;j++) if(sp[j]+0>=0.2) dominated=1
                    }
                }
                if(dominated) print
            }' | bgzip -c > ${meta.id}_spliceai_high.vcf.gz
        tabix -p vcf ${meta.id}_spliceai_high.vcf.gz
        SPLICEAI_FILE="${meta.id}_spliceai_high.vcf.gz"
    fi

    # --- Filter 6: High-confidence deleterious missense (REVEL/AlphaMissense) ---
    MISSENSE_FILE=""
    MISSENSE_FILTER=""
    [ "\${HAS_REVEL}" -eq 1 ] && MISSENSE_FILTER="INFO/REVEL>=0.644"
    if [ "\${HAS_AM}" -eq 1 ]; then
        [ -n "\${MISSENSE_FILTER}" ] && MISSENSE_FILTER="\${MISSENSE_FILTER} || INFO/AM_pathogenicity>=0.564" || MISSENSE_FILTER="INFO/AM_pathogenicity>=0.564"
    fi
    if [ -n "\${MISSENSE_FILTER}" ]; then
        bcftools view -f PASS -i "\${MISSENSE_FILTER}" ${vcf} \\
            -Oz -o ${meta.id}_missense_deleterious.vcf.gz
        bcftools index -t ${meta.id}_missense_deleterious.vcf.gz
        MISSENSE_FILE="${meta.id}_missense_deleterious.vcf.gz"
    fi

    # --- Merge all tiers into combined clinical VCF ---
    MERGE_FILES="${meta.id}_high_impact.vcf.gz ${meta.id}_rare_moderate.vcf.gz"
    [ -n "\${CLINVAR_FILE}" ] && MERGE_FILES="\${MERGE_FILES} \${CLINVAR_FILE}"
    [ -n "\${CADD_FILE}" ] && MERGE_FILES="\${MERGE_FILES} \${CADD_FILE}"
    [ -n "\${SPLICEAI_FILE}" ] && MERGE_FILES="\${MERGE_FILES} \${SPLICEAI_FILE}"
    [ -n "\${MISSENSE_FILE}" ] && MERGE_FILES="\${MERGE_FILES} \${MISSENSE_FILE}"

    bcftools concat -a -D \${MERGE_FILES} | \\
        bcftools sort -Oz -o ${meta.id}_clinical.vcf.gz
    bcftools index -t ${meta.id}_clinical.vcf.gz

    # --- Generate summary TSV ---
    {
        echo -e "CHROM\\tPOS\\tREF\\tALT\\tGT\\tIMPACT\\tGENE\\tCADD_PHRED\\tREVEL\\tAM_CLASS\\tCSQ_EXCERPT"
        bcftools view -H ${meta.id}_clinical.vcf.gz | \\
        awk -F'\\t' '{
            gt="."
            split(\$9, fmt, ":")
            split(\$10, vals, ":")
            for(i in fmt) if(fmt[i]=="GT") gt=vals[i]
            impact="."
            if(\$8 ~ /HIGH/) impact="HIGH"
            else if(\$8 ~ /MODERATE/) impact="MODERATE"
            else if(\$8 ~ /pathogenic/) impact="CLINVAR"
            else if(\$8 ~ /CADD_PHRED/) impact="CADD"
            else if(\$8 ~ /SpliceAI/) impact="SPLICEAI"
            gene="."
            if(match(\$8, /SYMBOL=[^;|]+/)) gene=substr(\$8, RSTART+7, RLENGTH-7)
            cadd="."
            if(match(\$8, /CADD_PHRED=[^;]+/)) cadd=substr(\$8, RSTART+11, RLENGTH-11)
            revel="."
            if(match(\$8, /REVEL=[^;]+/)) revel=substr(\$8, RSTART+6, RLENGTH-6)
            am="."
            if(match(\$8, /AM_class=[^;]+/)) am=substr(\$8, RSTART+9, RLENGTH-9)
            csq="."
            match(\$8, /CSQ=[^;]+/)
            if(RSTART>0) csq=substr(\$8, RSTART, RLENGTH>150?150:RLENGTH)
            print \$1"\\t"\$2"\\t"\$4"\\t"\$5"\\t"gt"\\t"impact"\\t"gene"\\t"cadd"\\t"revel"\\t"am"\\t"csq
        }' 2>/dev/null || true
    } > ${meta.id}_clinical_summary.tsv

    # Clean up intermediate tier files
    rm -f ${meta.id}_high_impact.vcf.gz* ${meta.id}_rare_moderate.vcf.gz*
    rm -f ${meta.id}_clinvar_pathogenic.vcf.gz* ${meta.id}_cadd_high.vcf.gz*
    rm -f ${meta.id}_spliceai_high.vcf.gz* ${meta.id}_missense_deleterious.vcf.gz*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_clinical.vcf.gz
    touch ${meta.id}_clinical.vcf.gz.tbi
    printf 'CHROM\\tPOS\\tREF\\tALT\\tGT\\tIMPACT\\tGENE\\tCADD_PHRED\\tREVEL\\tAM_CLASS\\tCSQ_EXCERPT\\n' > ${meta.id}_clinical_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
