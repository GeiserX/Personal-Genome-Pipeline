/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINICAL_FILTER — Extract clinically interesting variants from annotated VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Produces a small VCF of PASS variants that are:
      - Rare (MAX_AF < 1%, else gnomADe_AF and gnomADg_AF) AND HIGH/MODERATE VEP impact
      - OR HIGH/MODERATE with no known frequency: their own tier
        (*_freq_unknown.vcf.gz), never counted as rare
      - OR ClinVar pathogenic/likely pathogenic, any frequency: ClinVar_CLNSIG
        (the --clinvar file through VEP --custom), else the cache's CLIN_SIG
      - OR rare or of unknown frequency with a high CADD score (>= 20) outside HIGH/MODERATE
      - OR rare or of unknown frequency with a high SpliceAI delta score (>= 0.2), any gene of the value
      - OR rare or of unknown frequency with REVEL >= 0.644 or AlphaMissense >= 0.564
    With no frequency field no tier is filtered by frequency (as the script).

    bin/csq_rank.awk (on the task PATH, the program scripts/23 runs) ranks each
    record on its MANE Select transcript, else the canonical one, else any,
    gives it a frequency tier and checks its genotype call against the floors
    ext.min_gq (30), ext.min_dp (10), ext.het_ab_min and ext.het_ab_max
    (0.2, 0.8) and ext.hom_ab_min (0.8), as PGP_* INFO tags. The tiers select
    on those tags; bcftools view -i reads the INFO-level scores from vcfanno.
    Gene constraint columns are added by the script only (this process gets
    no constraint table).

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
    def min_gq = task.ext.min_gq != null ? task.ext.min_gq : 30
    def min_dp = task.ext.min_dp != null ? task.ext.min_dp : 10
    def het_lo = task.ext.het_ab_min != null ? task.ext.het_ab_min : 0.2
    def het_hi = task.ext.het_ab_max != null ? task.ext.het_ab_max : 0.8
    def hom_lo = task.ext.hom_ab_min != null ? task.ext.hom_ab_min : 0.8
    if (!(min_gq.toString() ==~ /[0-9]+/) || !(min_dp.toString() ==~ /[0-9]+/)) {
        error "CLINICAL_FILTER: ext.min_gq and ext.min_dp must be whole numbers, got '${min_gq}' and '${min_dp}'"
    }
    [het_lo, het_hi, hom_lo].each { ab ->
        if (!(ab.toString() ==~ /0(\.[0-9]+)?|1(\.0+)?/)) {
            error "CLINICAL_FILTER: ext.het_ab_min, ext.het_ab_max and ext.hom_ab_min must be fractions from 0 to 1, got '${ab}'"
        }
    }
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
    has_format() { grep -q "^##FORMAT=<ID=\$1," header.txt; }
    # rank KEEP: bin/csq_rank.awk on stdin (as scripts/23); KEEP=noncommon
    # drops the records with a known frequency of 1% or more.
    rank() {
        LC_ALL=C awk -v keep="\$1" -v min_gq=${min_gq} -v min_dp=${min_dp} \\
            -v het_lo=${het_lo} -v het_hi=${het_hi} -v hom_lo=${hom_lo} \\
            -f "\$(command -v csq_rank.awk)"
    }

    # --- Rarity: MAX_AF, else gnomADe_AF and gnomADg_AF (as scripts/23) ---
    HAS_FREQ=0
    for f in MAX_AF gnomADe_AF gnomADg_AF; do
        if grep -qx "\${f}" csq_fields.txt; then HAS_FREQ=1; fi
    done
    if [ "\${HAS_FREQ}" -eq 0 ]; then
        # Same as scripts/23: keep every MODERATE variant and say so.
        echo "NOTICE: no MAX_AF, gnomADe_AF or gnomADg_AF in the CSQ fields; no tier is filtered by frequency and every MODERATE variant is kept" >&2
    fi

    # --- PASS records every non-ClinVar tier starts from: ranked, the common ones left out ---
    bcftools view -f PASS ${vcf} | rank noncommon | bcftools view -Oz -o rare_pass.vcf.gz
    bcftools index -t rare_pass.vcf.gz
    # "rare" for the HIGH and MODERATE tiers: a known frequency below 1%, or
    # any record when the CSQ has no frequency field at all.
    RARE_FREQ='(INFO/PGP_FREQ="rare" || INFO/PGP_FREQ="none")'

    # --- Filter 1: HIGH impact on the ranked transcript, rare ---
    bcftools view -i "INFO/PGP_IMPACT=\\"HIGH\\" && \${RARE_FREQ}" rare_pass.vcf.gz \\
        -Oz -o ${meta.id}_high_impact.vcf.gz
    bcftools index -t ${meta.id}_high_impact.vcf.gz

    # --- Filter 2: MODERATE impact on the ranked transcript, rare ---
    bcftools view -i "INFO/PGP_IMPACT=\\"MODERATE\\" && \${RARE_FREQ}" rare_pass.vcf.gz \\
        -Oz -o ${meta.id}_rare_moderate.vcf.gz
    bcftools index -t ${meta.id}_rare_moderate.vcf.gz
    MERGE_FILES="${meta.id}_high_impact.vcf.gz ${meta.id}_rare_moderate.vcf.gz"

    # --- Filter 2b: HIGH or MODERATE with no known frequency: their own tier, not rare ---
    if [ "\${HAS_FREQ}" -eq 1 ]; then
        bcftools view -i '(INFO/PGP_IMPACT="HIGH" || INFO/PGP_IMPACT="MODERATE") && INFO/PGP_FREQ="unknown"' \\
            rare_pass.vcf.gz -Oz -o ${meta.id}_freq_unknown.vcf.gz
        bcftools index -t ${meta.id}_freq_unknown.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_freq_unknown.vcf.gz"
    fi

    # --- Filter 3: ClinVar pathogenic/likely pathogenic, any frequency ---
    # ClinVar_CLNSIG is the --clinvar file VEP annotated with --custom, so it
    # follows a ClinVar refresh; CLIN_SIG comes from the cache release and is
    # the fallback.
    if grep -qx ClinVar_CLNSIG csq_fields.txt; then
        echo "ClinVar tier: ClinVar_CLNSIG (the --clinvar file, through VEP --custom)" >&2
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c ClinVar_CLNSIG \\
                -i 'ClinVar_CLNSIG~"athogenic" && ClinVar_CLNSIG!~"onflicting"' \\
            | rank all | bcftools view -Oz -o ${meta.id}_clinvar_pathogenic.vcf.gz
        bcftools index -t ${meta.id}_clinvar_pathogenic.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_clinvar_pathogenic.vcf.gz"
    elif grep -qx CLIN_SIG csq_fields.txt; then
        echo "ClinVar tier: VEP's cached CLIN_SIG (set --clinvar for the current ClinVar file)" >&2
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c CLIN_SIG \\
                -i 'CLIN_SIG~"pathogenic" && CLIN_SIG!~"conflicting"' \\
            | rank all | bcftools view -Oz -o ${meta.id}_clinvar_pathogenic.vcf.gz
        bcftools index -t ${meta.id}_clinvar_pathogenic.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_clinvar_pathogenic.vcf.gz"
    fi

    # --- Filter 4: high CADD outside HIGH/MODERATE, rare ---
    CADD_EXPR=""
    has_info CADD_PHRED && CADD_EXPR="INFO/CADD_PHRED>=20"
    has_info CADD_PHRED_indel && CADD_EXPR="\${CADD_EXPR:+\${CADD_EXPR} || }INFO/CADD_PHRED_indel>=20"
    if [ -n "\${CADD_EXPR}" ]; then
        bcftools view rare_pass.vcf.gz \\
            -i "INFO/PGP_IMPACT!=\\"HIGH\\" && INFO/PGP_IMPACT!=\\"MODERATE\\" && (\${CADD_EXPR})" \\
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

    # --- Summary TSV (as scripts/23): IMPACT, GENE and Consequence of the
    #     ranked transcript, the worst consequence beside them, the frequency
    #     tier and the call check of bin/csq_rank.awk ---
    COL_FREQ='.'; grep -qx MAX_AF csq_fields.txt && COL_FREQ='%MAX_AF'
    COL_CADD='.'; has_info CADD_PHRED && COL_CADD='%INFO/CADD_PHRED'
    COL_REVEL='.'; has_info REVEL && COL_REVEL='%INFO/REVEL'
    COL_AM='.'; has_info AM_class && COL_AM='%INFO/AM_class'
    COL_GQ='.'; has_format GQ && COL_GQ='%GQ'
    COL_DP='.'; has_format DP && COL_DP='%DP'
    COL_AD='.'; has_format AD && COL_AD='%AD'
    {
        printf 'CHROM\\tPOS\\tREF\\tALT\\tGT\\tIMPACT\\tGENE\\tConsequence\\tMAX_AF\\tCADD_PHRED\\tREVEL\\tAM_CLASS\\tTRANSCRIPT\\tWORST_IMPACT\\tWORST_Consequence\\tFREQ_TIER\\tGQ\\tDP\\tAD\\tCALL_CHECK\\tCALL_WHY\\n'
        bcftools +split-vep ${meta.id}_clinical.vcf.gz -s worst \\
            -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\t%INFO/PGP_IMPACT\\t%INFO/PGP_GENE\\t%INFO/PGP_CSQ\\t\${COL_FREQ}\\t\${COL_CADD}\\t\${COL_REVEL}\\t\${COL_AM}\\t%INFO/PGP_TX\\t%IMPACT\\t%Consequence\\t%INFO/PGP_FREQ[\\t\${COL_GQ}\\t\${COL_DP}\\t\${COL_AD}]\\t%INFO/PGP_CALL\\t%INFO/PGP_CALL_WHY\\n" \\
            | awk -F'\\t' 'BEGIN {OFS = "\\t"} \$7 == "" {\$7 = "."} {print}'
    } > ${meta.id}_clinical_summary.tsv

    # Clean up intermediate tier files
    rm -f rare_pass.vcf.gz* csq_fields.txt header.txt
    rm -f ${meta.id}_high_impact.vcf.gz* ${meta.id}_rare_moderate.vcf.gz* ${meta.id}_freq_unknown.vcf.gz*
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
    printf 'CHROM\\tPOS\\tREF\\tALT\\tGT\\tIMPACT\\tGENE\\tConsequence\\tMAX_AF\\tCADD_PHRED\\tREVEL\\tAM_CLASS\\tTRANSCRIPT\\tWORST_IMPACT\\tWORST_Consequence\\tFREQ_TIER\\tGQ\\tDP\\tAD\\tCALL_CHECK\\tCALL_WHY\\n' > ${meta.id}_clinical_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
