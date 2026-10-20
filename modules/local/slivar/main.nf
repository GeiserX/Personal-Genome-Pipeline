/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SLIVAR — Variant prioritization and compound heterozygote detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prioritizes clinically interesting variants using tiered filters (rare HIGH,
    rare MODERATE + deleterious predictors, ClinVar pathogenic) and detects
    compound heterozygote candidates. Optionally annotates results with gnomAD
    gene constraint metrics (LOEUF, pLI).

    Rarity is the rule of scripts/23 and scripts/31: bin/csq_rank.awk (on the
    task PATH) gives each record PGP_FREQ, and a HIGH or MODERATE record with
    no known frequency goes into its own freq_unknown tier instead of the
    rare ones; the prioritized VCF still includes it. The summary TSV ends
    with MAX_AF, gnomADe_AF, gnomADg_AF, FREQ_TIER, TRANSCRIPT and
    MANE_Consequence, as the script's. This module writes no per-variant
    compound-het TSV (the script does).

    Three processes, each in its own pinned image:
      SLIVAR_PRIORITIZE  bcftools image: tiers, prioritized VCF, PED, summary TSV
      SLIVAR             slivar image:   compound-hets on the prioritized VCF
      COMPHET_PHASE      pypgx image:    each compound-het pair cis, trans or
                                         unknown from the sample's BAM, the
                                         pairs table scripts/31 writes

    Equivalent to: scripts/31-slivar.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SLIVAR_PRIORITIZE {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/slivar" }, mode: params.publish_dir_mode,
        pattern: "*_{prioritized.vcf.gz,prioritized.vcf.gz.tbi,slivar_summary.tsv}"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(gnomad_constraint)

    output:
    tuple val(meta), path("*_prioritized.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_prioritized.vcf.gz.tbi"), emit: vcf_index
    tuple val(meta), path("${meta.id}.ped"),           emit: ped
    tuple val(meta), path("*_slivar_summary.tsv"),     emit: summary_tsv
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_constraint = gnomad_constraint ? true : false
    """
    # --- Generate PED file for single sample (SLIVAR reads it) ---
    SAMPLE_NAME=\$(bcftools query -l ${vcf} | awk 'NR == 1')
    echo -e "\${SAMPLE_NAME}\\t\${SAMPLE_NAME}\\t0\\t0\\t0\\t-9" > ${meta.id}.ped

    # --- Which CSQ fields does the input carry? ---
    # Tested explicitly (as CLINICAL_FILTER does) instead of falling back with
    # `A | B || C | D`, which hid real failures behind the fallback branch.
    bcftools +split-vep -l ${vcf} > csq_fields.txt
    HAS_CLINSIG=0
    awk '\$2 == "CLIN_SIG" { f = 1 } END { exit !f }' csq_fields.txt && HAS_CLINSIG=1
    # Rarity, as scripts/31-slivar.sh: VEP's MAX_AF (highest frequency in any
    # 1000 Genomes or gnomAD exome/genome population), else gnomADe_AF and
    # gnomADg_AF, applied by bin/csq_rank.awk (PGP_FREQ). Exome frequency
    # alone calls a variant common in genomes but absent from exomes rare.
    has_csq() { awk -v f="\$1" '\$2 == f { x = 1 } END { exit !x }' csq_fields.txt; }
    rank() { LC_ALL=C awk -v keep="\$1" -f "\$(command -v csq_rank.awk)"; }
    HAS_FREQ=0
    for f in MAX_AF gnomADe_AF gnomADg_AF; do
        if has_csq "\${f}"; then HAS_FREQ=1; fi
    done
    if [ "\${HAS_FREQ}" -eq 0 ]; then
        echo "WARNING: no MAX_AF, gnomADe_AF or gnomADg_AF in the CSQ fields; rare tiers are not filtered by frequency." >&2
    fi
    # PASS records, ranked, the common ones left out. Rare: a known frequency
    # below 1%, or any record when the CSQ has no frequency field at all.
    bcftools view -f PASS ${vcf} | rank noncommon | bcftools view -Oz -o ranked.vcf.gz
    RARE_FREQ='(INFO/PGP_FREQ="rare" || INFO/PGP_FREQ="none")'

    # --- Filter 1: rare_high (PASS + HIGH impact + rare) ---
    bcftools +split-vep ranked.vcf.gz -c IMPACT -s worst -i 'IMPACT="HIGH"' | \\
        bcftools view -i "\${RARE_FREQ}" -Oz -o ${meta.id}_rare_high.vcf.gz
    bcftools index -t ${meta.id}_rare_high.vcf.gz

    # --- Filter 2: rare_moderate_deleterious (two-pass: CSQ then INFO predictors) ---
    # Pass 1: extract rare MODERATE via split-vep (CSQ fields)
    bcftools +split-vep ranked.vcf.gz -c IMPACT -s worst -i 'IMPACT="MODERATE"' | \\
        bcftools view -i "\${RARE_FREQ}" -Oz -o ${meta.id}_rare_moderate_all.vcf.gz
    bcftools index -t ${meta.id}_rare_moderate_all.vcf.gz

    # Pass 2: gate on deleteriousness predictors in INFO fields (added by vcfanno)
    # Build predictor expression from available INFO fields in the VCF header
    PREDICTOR_PARTS=""
    INFO_HEADER=\$(bcftools view -h ${meta.id}_rare_moderate_all.vcf.gz | grep '^##INFO' || true)

    grep -q 'ID=CADD_PHRED,' <<< "\$INFO_HEADER" && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/CADD_PHRED>=20"
    grep -q 'ID=CADD_PHRED_indel,' <<< "\$INFO_HEADER" && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/CADD_PHRED_indel>=20"
    grep -q 'ID=REVEL' <<< "\$INFO_HEADER" && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/REVEL>=0.5"
    grep -q 'ID=AM_class' <<< "\$INFO_HEADER" && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/AM_class=\\"likely_pathogenic\\""
    grep -q 'ID=SpliceAI,' <<< "\$INFO_HEADER" && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/SpliceAI!=\\".\\""
    grep -q 'ID=SpliceAI_indel,' <<< "\$INFO_HEADER" && \\
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

    MERGE_FILES="${meta.id}_rare_high.vcf.gz ${meta.id}_rare_moderate_del.vcf.gz"

    # --- Filter 2b: freq_unknown (PASS + HIGH or MODERATE + no known frequency; not rare) ---
    if [ "\${HAS_FREQ}" -eq 1 ]; then
        bcftools +split-vep ranked.vcf.gz -c IMPACT -s worst -i 'IMPACT="HIGH" || IMPACT="MODERATE"' | \\
            bcftools view -i 'INFO/PGP_FREQ="unknown"' -Oz -o ${meta.id}_freq_unknown.vcf.gz
        bcftools index -t ${meta.id}_freq_unknown.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_freq_unknown.vcf.gz"
    fi

    # --- Filter 3: clinvar_pathogenic (ranked too, for the summary columns) ---
    if [ "\${HAS_CLINSIG}" -eq 1 ]; then
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c CLIN_SIG \\
                -i 'CLIN_SIG~"pathogenic" && CLIN_SIG!~"conflicting"' \\
            | rank all | bcftools view -Oz -o ${meta.id}_clinvar_path.vcf.gz
        bcftools index -t ${meta.id}_clinvar_path.vcf.gz
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_clinvar_path.vcf.gz"
    else
        echo "ClinVar tier skipped (CLIN_SIG not in the CSQ fields)"
    fi

    # --- Merge tiers into prioritized VCF ---
    bcftools concat -a -D \${MERGE_FILES} | \\
        bcftools sort -Oz -o ${meta.id}_prioritized.vcf.gz
    bcftools index -t ${meta.id}_prioritized.vcf.gz

    # --- Generate summary TSV with optional gnomAD constraint enrichment ---
    # The last six columns as scripts/31: the frequency fields VEP wrote ('.'
    # for one it did not), the tier, the ranked transcript set and the MANE
    # Select consequence ('.' without a MANE transcript).
    EXTRA_FMT=""
    for f in MAX_AF gnomADe_AF gnomADg_AF; do
        if has_csq "\${f}"; then EXTRA_FMT="\${EXTRA_FMT}\\t%\${f}"; else EXTRA_FMT="\${EXTRA_FMT}\\t."; fi
    done
    {
        printf 'CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT\\tMAX_AF\\tgnomADe_AF\\tgnomADg_AF\\tFREQ_TIER\\tTRANSCRIPT\\tMANE_Consequence\\n'
        bcftools +split-vep \\
            ${meta.id}_prioritized.vcf.gz \\
            -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%IMPACT\\t%SYMBOL\\t%Consequence\\t%Existing_variation[\\t%GT]\${EXTRA_FMT}\\t%INFO/PGP_FREQ\\t%INFO/PGP_TX\\t%INFO/PGP_CSQ\\n" \\
            -s worst -d \\
            | awk -F'\\t' 'BEGIN {OFS = "\\t"} \$(NF - 1) != "MANE" {\$NF = "."} {print}'
    } > ${meta.id}_variants_raw.tsv

    if [ "${has_constraint}" = "true" ]; then
        # bin/constraint_join.awk (on the task PATH), the loader scripts/23 and
        # scripts/31 run: canonical rows only, the Ensembl row over the RefSeq
        # one, mis.z_score; it exits non-zero when rows carry gene symbols and
        # not one matches the table.
        awk -f "\$(command -v constraint_join.awk)" gene_col=SYMBOL constrained=1 \\
            ${gnomad_constraint} ${meta.id}_variants_raw.tsv > ${meta.id}_slivar_summary.tsv
    else
        mv ${meta.id}_variants_raw.tsv ${meta.id}_slivar_summary.tsv
    fi

    # Clean up intermediate files
    rm -f ${meta.id}_variants_raw.tsv ${meta.id}_rare_high.vcf.gz* ${meta.id}_rare_moderate_del.vcf.gz* \\
          ${meta.id}_rare_moderate_all.vcf.gz* ${meta.id}_clinvar_path.vcf.gz* ${meta.id}_freq_unknown.vcf.gz* \\
          ranked.vcf.gz csq_fields.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_prioritized.vcf.gz
    touch ${meta.id}_prioritized.vcf.gz.tbi
    touch ${meta.id}.ped
    printf 'CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT\\tMAX_AF\\tgnomADe_AF\\tgnomADg_AF\\tFREQ_TIER\\tTRANSCRIPT\\tMANE_Consequence\\n' > ${meta.id}_slivar_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process SLIVAR {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/slivar" }, mode: params.publish_dir_mode,
        pattern: "*_compound_hets.vcf.gz"

    input:
    tuple val(meta), path(vcf), path(ped)

    output:
    tuple val(meta), path("*_compound_hets.vcf.gz"), emit: compound_het_vcf
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # slivar writes the bgzipped VCF itself (the .gz name selects it), so its
    # own exit status fails the task: a crash must not read as "no compound
    # hets". Its stderr stays in the task log.
    slivar compound-hets \\
        --allow-non-trios \\
        --vcf ${vcf} \\
        --ped ${ped} \\
        --out-vcf ${meta.id}_compound_hets.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_compound_hets.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

/*
    COMPHET_PHASE — the phase of each compound-het pair, from the reads
    bin/comphet_phase.py (the rule scripts/31-slivar.sh runs) counts, for each
    pair SLIVAR found, the reads and read pairs that cover both positions:
    alt with alt or ref with ref is cis (one copy of the gene carries both),
    alt with ref is trans (both copies are hit). It writes
    <id>_compound_het_pairs.tsv with the verdict and the counts. bam and bai
    are [] for a sample without a BAM: every pair is then unknown, with the
    reason. pysam is in PYPGX_IMAGE (conf/containers.config).

    Equivalent to: scripts/31-slivar.sh, the phase of the pairs
*/
process COMPHET_PHASE {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/slivar" }, mode: params.publish_dir_mode,
        pattern: "*_compound_het_pairs.tsv"

    input:
    // Fixed names: a samplesheet path may hold any character, and these go on a command line.
    tuple val(meta), path(vcf, stageAs: 'compound_hets.vcf.gz'), path(bam, stageAs: 'sample.bam'), path(bai, stageAs: 'sample.bam.bai')

    output:
    tuple val(meta), path("${meta.id}_compound_het_pairs.tsv"), emit: pairs
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads = bam ? "--bam ${bam}" : "--skip 'no BAM for this sample in the samplesheet'"
    """
    comphet_phase.py --vcf ${vcf} --out ${meta.id}_compound_het_pairs.tsv ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'GENE\\tCHROM\\tPOS_1\\tREF_1\\tALT_1\\tPOS_2\\tREF_2\\tALT_2\\tDISTANCE\\tPHASE\\tCIS_READS\\tTRANS_READS\\tPHASE_NOTE\\n' > ${meta.id}_compound_het_pairs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
