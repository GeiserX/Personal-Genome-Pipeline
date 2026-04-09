/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SURVIVOR_MERGE — Merge SV calls from multiple callers (bcftools heuristic)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Finds structural variants called by 2+ callers at overlapping positions
    using a simplified bcftools-based approach with 1kb position binning.

    EXPERIMENTAL: This uses a heuristic position-binning approach.
    For production use, consider SURVIVOR or Jasmine with proper multi-sample
    VCF merging.

    Equivalent to: scripts/22-survivor-merge.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SURVIVOR_MERGE {
    tag "$meta.id"
    label 'process_low'

    container 'staphb/bcftools:1.21'

    publishDir "${params.outdir}/${meta.id}/sv_merged", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sv_vcfs)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_sv_consensus.vcf.gz"),     emit: merged_vcf
    tuple val(meta), path("${meta.id}_sv_consensus.vcf.gz.tbi"), emit: merged_vcf_index
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Tag each VCF with a caller index, extract SV positions
    CALLER_IDX=0
    for VCF_FILE in ${sv_vcfs}; do
        CALLER_IDX=\$((CALLER_IDX + 1))
        bcftools view -f PASS,. "\$VCF_FILE" 2>/dev/null | \\
            grep -v '^#' | \\
            awk -F'\\t' -v caller=\$CALLER_IDX '{
                chrom=\$1; pos=\$2; info=\$8;
                end=pos;
                if(match(info, /END=[0-9]+/)) end=substr(info, RSTART+4, RLENGTH-4);
                svtype="UNK";
                if(match(info, /SVTYPE=[A-Z]+/)) svtype=substr(info, RSTART+7, RLENGTH-7);
                bin=int(pos/1000);
                printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t.\\tN\\t<%s>\\t.\\tPASS\\tSVTYPE=%s;END=%s\\n",
                    chrom, bin, svtype, caller, pos, chrom, pos, svtype, svtype, end;
            }' || true
    done > all_sv_tagged.tsv

    # Find bins seen by 2+ distinct callers
    awk -F'\\t' '{
        key=\$1"_"\$2"_"\$3;
        caller=\$4;
        if(!(key in seen)) {
            seen[key]=caller;
            line[key]=\$6"\\t"\$7"\\t.\\tN\\t"\$10"\\t.\\tPASS\\t"\$13;
        } else if(index(seen[key], caller) == 0) {
            seen[key]=seen[key]"|"caller;
        }
    } END {
        for(k in seen) {
            n=split(seen[k], a, "|");
            if(n >= 2) print line[k];
        }
    }' all_sv_tagged.tsv | \\
        sort -k1,1V -k2,2n > consensus_raw.txt

    # Build a valid VCF with contig headers from reference FAI
    {
        echo '##fileformat=VCFv4.2'
        echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">'
        echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">'
        if [ -f "${reference_fai}" ]; then
            awk -F'\\t' '{printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2}' "${reference_fai}"
        fi
        printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n'
        cat consensus_raw.txt
    } > ${meta.id}_sv_consensus_unsorted.vcf

    # Sort, compress, and index
    bcftools sort ${meta.id}_sv_consensus_unsorted.vcf -Oz \\
        -o ${meta.id}_sv_consensus.vcf.gz
    bcftools index -t ${meta.id}_sv_consensus.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_consensus.vcf.gz
    touch ${meta.id}_sv_consensus.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
