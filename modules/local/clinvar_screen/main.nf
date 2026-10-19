/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ClinVar Pathogenic Screen — intersect sample VCF with ClinVar pathogenic/LP variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Finds known disease-causing variants the person carries by intersecting
    PASS variants against the ClinVar pathogenic subset. Writes
    <id>_clinvar_hits.vcf: the sample's matching records the person carries
    (genotype with an ALT allele) with ClinVar's ID, GENEINFO, CLNSIG and
    CLNREVSTAT copied on; <id>_clinvar_hits.tsv: the same hits, one row each;
    and prints the hits per ClinVar review-status star tier.

    Both files are read only on the contigs they share with each other and
    with the reference: bcftools norm stops at the first record whose contig
    the reference lacks (a vendor scaffold, ClinVar's NT_ contigs), and a
    record on a contig the other file lacks cannot match anyway.

    Two neighbouring records that can be one allele are also matched as their
    combined allele (bin/clinvar_join_adjacent.awk): a matched combined
    record is a hit of its own (INFO/PGP_JOINED_FROM) and its parts are
    marked PGP_PART_OF.

    A hit whose call has a genotype quality below ext.min_gq (default 30) or
    a depth below ext.min_dp (default 10) is a low-confidence call; a floor
    other than the default goes into the hits VCF as
    ##PGP_CALL_FLOOR=GQ=..,DP=.., which bin/clinvar_hits.awk and the report
    read (without the line they use the default). Set it in a config with
    withName: 'CLINVAR_SCREEN' { ext.min_gq = 20; ext.min_dp = 8 }
    (scripts/06 reads CLINVAR_MIN_GQ and CLINVAR_MIN_DP).

    CLINVAR_CDNA_CHECK, below, runs scripts/06's cDNA read check
    (bin/clip_check.py) on these hits: <id>_clinvar_cdna_check.tsv says for
    each hit whether its alt reads look like spliced cDNA. It needs the
    sample's BAM and --gencode_exons; without either it writes every hit as
    not_checked, with the reason, as scripts/06 does.

    Equivalent to: scripts/06-clinvar-screen.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINVAR_SCREEN {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/clinvar" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(clinvar)
    path(clinvar_index)
    path(reference)
    path(reference_fai)

    output:
    // The emit keeps its old name (isec_dir) so workflows/pgx.nf needs no change;
    // it now carries the annotated hits VCF, not the isec directory.
    tuple val(meta), path("${meta.id}_clinvar_hits.vcf"), emit: isec_dir
    tuple val(meta), path("${meta.id}_clinvar_hits.tsv"), emit: hits_tsv
    tuple val(meta), path("*_pass.vcf.gz"),      emit: pass_vcf
    tuple val(meta), path("*_pass.vcf.gz.tbi"),  emit: pass_vcf_index
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_gq = task.ext.min_gq != null ? task.ext.min_gq : 30
    def min_dp = task.ext.min_dp != null ? task.ext.min_dp : 10
    if (!(min_gq.toString() ==~ /[0-9]+/) || !(min_dp.toString() ==~ /[0-9]+/)) {
        error "CLINVAR_SCREEN: ext.min_gq and ext.min_dp must be whole numbers, got '${min_gq}' and '${min_dp}'"
    }
    def floor_line = (min_gq.toString() == '30' && min_dp.toString() == '10') ? '' : "##PGP_CALL_FLOOR=GQ=${min_gq},DP=${min_dp}"
    """
    set -euo pipefail

    # Step 0: the contigs both files hold records on, that the reference has too.
    # No contig in common would give zero hits that look clean, so stop and say
    # which file is named the other way.
    cut -f1 ${reference_fai} | sort -u > ref_contigs.txt
    bcftools index -s ${vcf} | awk '\$3 > 0 { print \$1 }' | sort -u > sample_contigs.txt
    bcftools index -s ${clinvar} | awk '\$3 > 0 { print \$1 }' | sort -u > clinvar_contigs.txt
    comm -12 clinvar_contigs.txt ref_contigs.txt > clinvar_keep.txt
    comm -12 sample_contigs.txt clinvar_keep.txt > sample_keep.txt
    if [ ! -s sample_keep.txt ]; then
        echo "ERROR: ${vcf} and ${clinvar} have no contig name in common." >&2
        echo "  Sample contigs:  \$(awk 'NR <= 3' sample_contigs.txt | paste -sd' ' -) ..." >&2
        echo "  ClinVar contigs: \$(awk 'NR <= 3' clinvar_contigs.txt | paste -sd' ' -) ..." >&2
        if ! grep -q '^chr' sample_contigs.txt; then
            echo "  The sample VCF is not chr-named: rename its contigs (docs/vcf-first.md)." >&2
        elif ! grep -q '^chr' clinvar_contigs.txt; then
            echo "  The ClinVar file is not chr-named: use clinvar_pathogenic_chr.vcf.gz from scripts/setup.sh." >&2
        fi
        exit 1
    fi
    # Targets with a constant end: a header without contig lengths cannot shorten them
    awk -v OFS='\\t' '{ print \$1, 1, 2147483647 }' clinvar_keep.txt > clinvar_targets.tsv
    awk -v OFS='\\t' '{ print \$1, 1, 2147483647 }' sample_keep.txt > sample_targets.tsv
    comm -23 sample_contigs.txt sample_keep.txt > sample_dropped.txt

    # Step 1: ClinVar split and left-aligned. Reuse the copy step 06 (bash) keeps beside
    # the source file when it is newer than the source; otherwise normalise here.
    SRC=\$(readlink -f ${clinvar})
    NORM_BESIDE="\${SRC%.vcf.gz}.norm.vcf.gz"
    if [ -f "\${NORM_BESIDE}" ] && [ -f "\${NORM_BESIDE}.tbi" ] && [ "\${NORM_BESIDE}" -nt "\${SRC}" ]; then
        CLINVAR_NORM="\${NORM_BESIDE}"
    else
        bcftools view -T clinvar_targets.tsv -Ou ${clinvar} \\
            | bcftools norm -m -any -c w -f ${reference} -Oz -o clinvar_norm.vcf.gz -
        bcftools index -t clinvar_norm.vcf.gz
        CLINVAR_NORM=clinvar_norm.vcf.gz
    fi

    # Step 2: PASS records, or '.,PASS' when the caller never writes PASS
    HAS_PASS=\$( (bcftools view -H -f PASS ${vcf} || true) | head -n 1 | wc -l)
    if [ "\${HAS_PASS}" -gt 0 ]; then
        FILTER=PASS
    else
        FILTER=.,PASS
        echo "NOTICE: ${vcf} has no PASS record; using records with FILTER '.' or PASS"
    fi
    if [ -s sample_dropped.txt ]; then
        awk -v OFS='\\t' '{ print \$1, 1, 2147483647 }' sample_dropped.txt > sample_dropped.tsv
        LEFT_OUT=\$(bcftools view -H -f "\${FILTER}" -R sample_dropped.tsv ${vcf} | wc -l | tr -d ' ')
        echo "NOTICE: \${LEFT_OUT} '\${FILTER}' records left out: they lie on \$(grep -c . sample_dropped.txt) contigs that ClinVar or the reference lacks (first ones: \$(awk 'NR <= 5' sample_dropped.txt | paste -sd, -))"
    fi
    bcftools view -f "\${FILTER}" -T sample_targets.tsv -Ou ${vcf} | \\
        bcftools norm -m -any -c w -f ${reference} -Oz -o ${meta.id}_pass.vcf.gz -
    bcftools index -t ${meta.id}_pass.vcf.gz
    if [ "\$(bcftools index -n ${meta.id}_pass.vcf.gz)" -eq 0 ]; then
        echo "ERROR: no records left after the '\${FILTER}' filter on the shared contigs" >&2
        exit 1
    fi

    # Step 2b: the combined alleles. A caller can write one complex allele (a
    # delins) as two records a few bases apart on the same haplotype; matched
    # one by one, the parts can hit single-submitter ClinVar records and miss
    # the better reviewed record for the whole allele. bcftools norm -m +any
    # joins only records at the same position, so clinvar_join_adjacent.awk
    # (on the task PATH) builds the allele of each pair of neighbours that can
    # be one allele, with the reference bases between the parts
    # (+fill-from-fasta), left-aligned the way ClinVar is. The candidates
    # matched below are the split calls plus these.
    JOIN_AWK="\$(command -v clinvar_join_adjacent.awk)"
    bcftools view ${meta.id}_pass.vcf.gz \\
        | awk -v mode=prepare -f "\${JOIN_AWK}" \\
        | bcftools +fill-from-fasta --no-version -Ov - -- -c INFO/PGP_JOIN_SEQ -f ${reference} \\
        | awk -v mode=finish -f "\${JOIN_AWK}" \\
        | bcftools norm --no-version -c w -d exact -f ${reference} -Oz -o joined.vcf.gz -
    bcftools index -t joined.vcf.gz
    bcftools concat --no-version -a -Oz -o candidates.vcf.gz ${meta.id}_pass.vcf.gz joined.vcf.gz
    bcftools index -t candidates.vcf.gz
    echo "Combined alleles of neighbouring records: \$(bcftools index -n joined.vcf.gz)"

    # Step 3: sample records matching a ClinVar allele, with ClinVar's ID, gene,
    # significance and review status copied on, kept only where the genotype
    # carries an ALT allele. The match is by allele, so after the split above a
    # 0/0 or ./. record, the 0/0 half of a multiallelic record and a
    # reference-only ALT '.' row would otherwise be listed as hits. The test is a
    # non-zero allele index anywhere in GT, not GT="alt", which drops a half
    # call such as ./1 that does carry the ALT allele.
    # (annotate -a needs an indexed target, so the shared records go to a file first.)
    # --no-version: no command line, with ClinVar's absolute path, in the header.
    # Last, the tag pass marks each part of a matched combined allele
    # (PGP_PART_OF), and two hits in one gene that can be one allele but
    # matched no combined record (PGP_ONE_ALLELE).
    bcftools isec --no-version -n=2 -w1 -Oz -o shared.vcf.gz candidates.vcf.gz "\${CLINVAR_NORM}"
    bcftools index -t shared.vcf.gz
    bcftools annotate --no-version -a "\${CLINVAR_NORM}" --pair-logic exact \\
        -c ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT -Ou shared.vcf.gz \\
        | bcftools view --no-version -i 'GT~"[1-9]"' -Ov -o hits_untagged.vcf
    # A floor other than the default, for bin/clinvar_hits.awk and the report.
    awk -v mode=tag -f "\${JOIN_AWK}" hits_untagged.vcf hits_untagged.vcf \\
        | awk -v l="${floor_line}" 'l != "" && /^#CHROM/ { print l } { print }' > ${meta.id}_clinvar_hits.vcf

    # The same hits as a table, one row per hit
    printf 'chrom\\tpos\\tref\\talt\\tgenotype\\tclinvar_id\\tgeneinfo\\tclnsig\\tclnrevstat\\n' > ${meta.id}_clinvar_hits.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\t%ID\\t%INFO/GENEINFO\\t%INFO/CLNSIG\\t%INFO/CLNREVSTAT\\n' \\
        ${meta.id}_clinvar_hits.vcf >> ${meta.id}_clinvar_hits.tsv

    HITS=\$(grep -c -v '^#' ${meta.id}_clinvar_hits.vcf || true)
    echo "ClinVar pathogenic hits: \${HITS}"
    # Grouped by ClinVar review stars, as scripts/06 prints them
    # (bin/clinvar_hits.awk is on the task PATH).
    if [ "\${HITS}" -gt 0 ]; then
        echo "By review status (ClinVar stars):"
        awk -f "\$(command -v clinvar_hits.awk)" ${meta.id}_clinvar_hits.vcf | cut -f1 | sort -nr | uniq -c \\
            | awk '{printf "  %d star%s: %d\\n", \$2, (\$2 == 1 ? "" : "s"), \$1}'
        LOW=\$(awk -f "\$(command -v clinvar_hits.awk)" ${meta.id}_clinvar_hits.vcf | awk -F'\\t' '\$13 == "yes"' | wc -l | tr -d ' ')
        echo "Low-confidence calls among them (GQ below ${min_gq} or DP below ${min_dp}): \${LOW}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_clinvar_hits.vcf ${meta.id}_clinvar_hits.tsv
    touch ${meta.id}_pass.vcf.gz
    touch ${meta.id}_pass.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

/*
    CLINVAR_CDNA_CHECK — step 06's cDNA read check of each ClinVar hit
    Spliced cDNA or expression-vector DNA (a lab contaminant) gives alt reads
    clipped exactly at an exon edge, with the rest of the read on the next
    exon; bin/clip_check.py counts them for each hit and flags the hit when
    most alt reads are clipped that way and few ref reads are. It reads only
    the reads at the hits. pysam is in PYPGX_IMAGE (conf/containers.config).
    bam and bai are [] for a sample without a BAM, exons [] without
    --gencode_exons: the table then says "not checked" and why.

    Equivalent to: scripts/06-clinvar-screen.sh, step 4
*/
process CLINVAR_CDNA_CHECK {
    tag "$meta.id"
    label 'process_single'

    // The table only: CLINVAR_SCREEN's versions.yml is the one in clinvar/.
    publishDir { "${params.outdir}/${meta.id}/clinvar" }, mode: params.publish_dir_mode, pattern: '*_clinvar_cdna_check.tsv'

    input:
    // Fixed names: a samplesheet path may hold any character, and these go on a command line.
    tuple val(meta), path(hits_tsv), path(bam, stageAs: 'sample.bam'), path(bai, stageAs: 'sample.bam.bai')
    path(exons, stageAs: 'exons.bed')

    output:
    tuple val(meta), path("${meta.id}_clinvar_cdna_check.tsv"), emit: cdna_check
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads = !bam ? "--skip 'no BAM for this sample in the samplesheet'" :
                !exons ? "--skip 'no GENCODE exons (--gencode_exons)'" :
                "--bam ${bam} --exons ${exons}"
    """
    clip_check.py --hits ${hits_tsv} --out ${meta.id}_clinvar_cdna_check.tsv ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_clinvar_cdna_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
