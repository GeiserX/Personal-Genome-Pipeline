/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINICAL_FILTER — Extract clinically interesting variants from annotated VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Produces a small VCF of PASS variants that are:
      - Rare (MAX_AF < 1%, else gnomADe_AF and gnomADg_AF) AND HIGH/MODERATE VEP impact
      - OR ClinVar pathogenic/likely pathogenic (VEP's CLIN_SIG; any frequency)
      - OR rare with a high CADD score (>= 20) outside HIGH/MODERATE
      - OR rare with a high SpliceAI delta score (>= 0.2), any gene of the value
      - OR rare with REVEL >= 0.644 or AlphaMissense >= 0.564
    With no frequency field no tier is filtered by frequency (as the script).

    Uses bcftools +split-vep for VEP CSQ fields and bcftools view -i for
    INFO-level annotations from vcfanno. Gene constraint columns are added by
    the script only (this process gets no constraint table).

    Equivalent to: scripts/23-clinical-filter.sh (bcftools portion only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINICAL_FILTER {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/clinical" }, mode: params.publish_dir_mode

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
    bcftools +split-vep -l ${vcf} | cut -f2 > csq_fields.txt
    for f in IMPACT SYMBOL Consequence; do
        if ! grep -qx "\${f}" csq_fields.txt; then
            echo "ERROR: CLINICAL_FILTER requires a VEP-annotated VCF with the CSQ \${f} field." >&2
            echo "Enable 'vep' in --tools before 'clinical_filter'." >&2
            exit 1
        fi
    done
    bcftools view -h ${vcf} > header.txt
    has_info() { grep -q "^##INFO=<ID=\$1," header.txt; }

    # --- Rarity: MAX_AF, else gnomADe_AF and gnomADg_AF (as scripts/23) ---
    FREQ_COLS=""
    FREQ_EXPR=""
    if grep -qx MAX_AF csq_fields.txt; then
        FREQ_COLS="MAX_AF:Float"
        FREQ_EXPR='(MAX_AF<0.01 || MAX_AF=".")'
    else
        for f in gnomADe_AF gnomADg_AF; do
            grep -qx "\${f}" csq_fields.txt || continue
            FREQ_COLS="\${FREQ_COLS:+\${FREQ_COLS},}\${f}:Float"
            FREQ_EXPR="\${FREQ_EXPR:+\${FREQ_EXPR} && }(\${f}<0.01 || \${f}=\\".\\")"
        done
    fi

    # --- PASS records every non-ClinVar tier starts from, rare when a frequency exists ---
    if [ -n "\${FREQ_EXPR}" ]; then
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c "\${FREQ_COLS}" -s worst -i "\${FREQ_EXPR}" -Oz -o rare_pass.vcf.gz
    else
        # Same as scripts/23: keep every MODERATE variant and say so.
        echo "NOTICE: no MAX_AF, gnomADe_AF or gnomADg_AF in the CSQ fields; no tier is filtered by frequency and every MODERATE variant is kept" >&2
        bcftools view -f PASS ${vcf} -Oz -o rare_pass.vcf.gz
    fi
    bcftools index -t rare_pass.vcf.gz

    # --- Filter 1: HIGH impact, rare ---
    bcftools +split-vep rare_pass.vcf.gz -c IMPACT -s worst -i 'IMPACT="HIGH"' \\
        -Oz -o ${meta.id}_high_impact.vcf.gz
    bcftools index -t ${meta.id}_high_impact.vcf.gz

    # --- Filter 2: MODERATE impact, rare ---
    bcftools +split-vep rare_pass.vcf.gz -c IMPACT -s worst -i 'IMPACT="MODERATE"' \\
        -Oz -o ${meta.id}_rare_moderate.vcf.gz
    bcftools index -t ${meta.id}_rare_moderate.vcf.gz
    MERGE_FILES="${meta.id}_high_impact.vcf.gz ${meta.id}_rare_moderate.vcf.gz"

    # --- Filter 3: ClinVar pathogenic/likely pathogenic, any frequency ---
    # VEP's CLIN_SIG comes from its cache release, not from --clinvar.
    if grep -qx CLIN_SIG csq_fields.txt; then
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c CLIN_SIG \\
                -i 'CLIN_SIG~"pathogenic" && CLIN_SIG!~"conflicting"' \\
                -Oz -o ${meta.id}_clinvar_pathogenic.vcf.gz
        bcftools index -t ${meta.id}_clinvar_pathogenic.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_clinvar_pathogenic.vcf.gz"
    fi

    # --- Filter 4: high CADD outside HIGH/MODERATE, rare ---
    CADD_EXPR=""
    has_info CADD_PHRED && CADD_EXPR="INFO/CADD_PHRED>=20"
    has_info CADD_PHRED_indel && CADD_EXPR="\${CADD_EXPR:+\${CADD_EXPR} || }INFO/CADD_PHRED_indel>=20"
    if [ -n "\${CADD_EXPR}" ]; then
        bcftools +split-vep rare_pass.vcf.gz -c IMPACT -s worst \\
            -i "IMPACT!=\\"HIGH\\" && IMPACT!=\\"MODERATE\\" && (\${CADD_EXPR})" \\
            -Oz -o ${meta.id}_cadd_high.vcf.gz
        bcftools index -t ${meta.id}_cadd_high.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_cadd_high.vcf.gz"
    fi

    # --- Filter 5: SpliceAI cryptic splice, rare; every gene of a ','-joined value ---
    SPLICEAI_PREFILTER=""
    has_info SpliceAI && SPLICEAI_PREFILTER='INFO/SpliceAI!="."'
    has_info SpliceAI_indel && SPLICEAI_PREFILTER="\${SPLICEAI_PREFILTER:+\${SPLICEAI_PREFILTER} || }INFO/SpliceAI_indel!=\\".\\""
    if [ -n "\${SPLICEAI_PREFILTER}" ]; then
        bcftools view -i "\${SPLICEAI_PREFILTER}" rare_pass.vcf.gz | \\
            awk -F'\\t' 'BEGIN{OFS="\\t"} /^#/{print;next} {
                hit=0
                n=split(\$8, kv, ";")
                for(i=1;i<=n;i++){
                    if(kv[i] !~ /^SpliceAI(_indel)?=/) continue
                    v=kv[i]; sub(/^[^=]*=/,"",v)
                    na=split(v, genes, ",")
                    for(a=1;a<=na;a++){
                        split(genes[a], sp, "|")
                        for(j=3;j<=6;j++) if(sp[j]!="" && sp[j]!="." && sp[j]+0>=0.2) hit=1
                    }
                }
                if(hit) print
            }' | bcftools view -Oz -o ${meta.id}_spliceai_high.vcf.gz
        bcftools index -t ${meta.id}_spliceai_high.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_spliceai_high.vcf.gz"
    fi

    # --- Filter 6: REVEL >= 0.644 (ClinGen PP3 Supporting) or AlphaMissense >= 0.564
    #     (its likely_pathogenic class boundary), rare ---
    MISSENSE_FILTER=""
    has_info REVEL && MISSENSE_FILTER="INFO/REVEL>=0.644"
    has_info AM_pathogenicity && MISSENSE_FILTER="\${MISSENSE_FILTER:+\${MISSENSE_FILTER} || }INFO/AM_pathogenicity>=0.564"
    if [ -n "\${MISSENSE_FILTER}" ]; then
        bcftools view -i "\${MISSENSE_FILTER}" rare_pass.vcf.gz -Oz -o ${meta.id}_missense_deleterious.vcf.gz
        bcftools index -t ${meta.id}_missense_deleterious.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_missense_deleterious.vcf.gz"
    fi

    # --- Merge all tiers into combined clinical VCF ---
    bcftools concat -a -D \${MERGE_FILES} | \\
        bcftools sort -Oz -o ${meta.id}_clinical.vcf.gz
    bcftools index -t ${meta.id}_clinical.vcf.gz

    # --- Summary TSV: gene, impact and consequence of the worst consequence ---
    COL_FREQ='.'; grep -qx MAX_AF csq_fields.txt && COL_FREQ='%MAX_AF'
    COL_CADD='.'; has_info CADD_PHRED && COL_CADD='%INFO/CADD_PHRED'
    COL_REVEL='.'; has_info REVEL && COL_REVEL='%INFO/REVEL'
    COL_AM='.'; has_info AM_class && COL_AM='%INFO/AM_class'
    {
        printf 'CHROM\\tPOS\\tREF\\tALT\\tGT\\tIMPACT\\tGENE\\tConsequence\\tMAX_AF\\tCADD_PHRED\\tREVEL\\tAM_CLASS\\n'
        bcftools +split-vep ${meta.id}_clinical.vcf.gz -s worst \\
            -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\t%IMPACT\\t%SYMBOL\\t%Consequence\\t\${COL_FREQ}\\t\${COL_CADD}\\t\${COL_REVEL}\\t\${COL_AM}\\n" \\
            | awk -F'\\t' 'BEGIN {OFS = "\\t"} \$7 == "" {\$7 = "."} {print}'
    } > ${meta.id}_clinical_summary.tsv

    # Clean up intermediate tier files
    rm -f rare_pass.vcf.gz* csq_fields.txt header.txt
    rm -f ${meta.id}_high_impact.vcf.gz* ${meta.id}_rare_moderate.vcf.gz*
    rm -f ${meta.id}_clinvar_pathogenic.vcf.gz* ${meta.id}_cadd_high.vcf.gz*
    rm -f ${meta.id}_spliceai_high.vcf.gz* ${meta.id}_missense_deleterious.vcf.gz*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_clinical.vcf.gz
    touch ${meta.id}_clinical.vcf.gz.tbi
    printf 'CHROM\\tPOS\\tREF\\tALT\\tGT\\tIMPACT\\tGENE\\tConsequence\\tMAX_AF\\tCADD_PHRED\\tREVEL\\tAM_CLASS\\n' > ${meta.id}_clinical_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
