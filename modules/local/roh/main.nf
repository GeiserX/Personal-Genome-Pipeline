/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ROH — Runs of Homozygosity (consanguinity screening)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Detects long stretches of homozygous DNA that indicate shared ancestry.
    The summary keeps chr1 to chr22 only and leaves out the part of each segment in a
    centromere or heterochromatin band (masked_bands below), where few
    markers map. No population allele frequencies are used (--AF-dflt 0.4).

    Equivalent to: scripts/11-roh-analysis.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ROH {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/roh" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*_roh.txt"),         emit: roh_regions
    tuple val(meta), path("*_roh_summary.txt"), emit: roh_summary
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Auto-detect chip data: if FORMAT/PL is absent, use -G30 (genotype-only mode)
    HAS_PL=\$(bcftools view -h ${vcf} | grep -c '##FORMAT=<ID=PL' || true)

    ROH_FLAGS="--AF-dflt 0.4"
    if [ "\${HAS_PL}" -eq 0 ]; then
        ROH_FLAGS="\${ROH_FLAGS} -G30"
    fi

    # Only PASS calls (and records with no filter, as chip VCFs have), as in
    # scripts/11-roh-analysis.sh. pipefail: a failed view must fail the task.
    set -o pipefail
    bcftools view -f PASS,. -Ou ${vcf} \\
        | bcftools roh \${ROH_FLAGS} -o ${meta.id}_roh.txt -

    # Summary: autosomal segments with 5 Mb or more outside GRCh38's
    # centromere and heterochromatin bands (the acen, gvar and stalk rows of
    # UCSC's cytoBand.hg38.txt, merged; the line of scripts/11-roh-analysis.sh
    # and bin/collect_summary.py), the threshold of docs/11-roh-analysis.md;
    # length_bp and length_mb are of that part. The bash step writes the same
    # file.
    MASKED_BANDS="chr1:121700000-143200000 chr2:91800000-96000000 chr3:87800000-98600000 chr4:48200000-51800000 chr5:46100000-51400000 chr6:58500000-62600000 chr7:58100000-62100000 chr8:43200000-47200000 chr9:42200000-61500000 chr10:38000000-41600000 chr11:51000000-55800000 chr12:33200000-37800000 chr13:0-18900000 chr14:0-18200000 chr15:0-20500000 chr16:35300000-47000000 chr17:22700000-27400000 chr18:15400000-21500000 chr19:19900000-31900000 chr20:25700000-30400000 chr21:0-13000000 chr22:0-17400000 chrX:58100000-63800000 chrY:10300000-10600000 chrY:26600000-57227415"
    echo "# ROH Summary for ${meta.id}" > ${meta.id}_roh_summary.txt
    echo "# Segments >=5MB on autosomes (excludes chrX/chrY)" >> ${meta.id}_roh_summary.txt
    printf 'chrom\\tstart\\tend\\tlength_bp\\tlength_mb\\n' >> ${meta.id}_roh_summary.txt
    awk -v bands="\$MASKED_BANDS" '
        BEGIN { n = split(bands, b, " "); for (i = 1; i <= n; i++) { split(b[i], x, /[:-]/); c[i] = x[1]; s[i] = x[2]; e[i] = x[3] } }
        \$1 == "RG" && \$3 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/ {
          m = 0
          for (i = 1; i <= n; i++) if (c[i] == \$3) { o = (\$5 < e[i] ? \$5 : e[i]) - (\$4 - 1 > s[i] ? \$4 - 1 : s[i]); if (o > 0) m += o }
          if (m > \$6) m = \$6
          l = \$6 - m
          if (l >= 5000000) printf "%s\\t%s\\t%s\\t%d\\t%.1f\\n", \$3, \$4, \$5, l, l / 1e6
        }' ${meta.id}_roh.txt >> ${meta.id}_roh_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_roh.txt ${meta.id}_roh_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
