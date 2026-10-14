/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VCF_PRECHECK — Look at each input VCF once, before any analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Reports three facts per sample; main.nf decides what to do with them,
    so every stop names the sample and says how to fix the file. The
    sample count is the exception: more than one sample stops here.

    FILTER values. ClinVar screen, clinical filter and slivar keep only
    FILTER=PASS records. A VCF from a caller that leaves FILTER as '.'
    (unfiltered GATK HaplotypeCaller, FreeBayes) would give zero records
    everywhere and a report with zero hits. FILTER_STATUS:
      pass        — at least one PASS record; the VCF is used as given
      no_pass     — no PASS record; main.nf stops the run naming the sample
                    (also with --allow_unfiltered when no record has FILTER '.',
                    since the relaxed copy would still have no PASS record)
      unfiltered  — no PASS record, at least one FILTER '.' record, and
                    --allow_unfiltered is set: a copy with
                    FILTER '.' rewritten to PASS is emitted, which is what
                    `bcftools view -f .,PASS` would select (records with any
                    other FILTER value stay excluded)
    The fallback applies only when the file has no PASS record at all, and
    the log says so. A file with some PASS records is never relaxed.

    Contig names. The steps expect GRCh38 names with chr (chr1, chrM). With
    Ensembl names (1, MT) the mito haplogroup file comes out empty and chrX
    segments leak into the autosomal ROH summary, both with exit 0.
    CONTIG_STYLE is chr when any contig holding records starts with chr,
    other when none does, unknown when the index lists none; CONTIGS_SEEN
    holds the first five, for the message.

    gVCF. PharmCAT refuses a gVCF, and it decides by the file name too
    (.g.vcf, .genomic.vcf). GVCF:
      blocks  — a ##GVCFBlock header line, or a reference-only record (ALT
                <*>, <NON_REF> or '.') that carries INFO/END
      name    — no blocks, but the file name contains .g.vcf or .genomic.vcf
      none    — neither

    Samples. One samplesheet row is one person. A joint-called VCF holds
    one genotype column per person: slivar keeps the first column's name
    and the other steps read every column, so the results would mix people
    or cover one of them without saying so. A VCF with more than one
    sample stops here, with the sample count, the first five names and the
    bcftools command that keeps one of them.

    Build. The pipeline is GRCh38 only. A ##contig line for chr1-22, chrX or
    chrY (or 1-22, X, Y) whose length is not GRCh38's stops here, naming the
    contigs: a GRCh37 (hg19) file named the chr way would otherwise run to
    the end with every coordinate wrong. The lengths are those of the
    reference setup.sh installs (NCBI's no-ALT analysis set), the same for
    every GRCh38 build. A header that gives no length for chr1-22, X or Y
    (no ##contig lines, or lines without lengths) cannot show the build, so
    the REF bases of its first 1000 records (less any on a contig the
    reference lacks) are compared with the reference, as validate-setup.sh
    does: bcftools norm -c x drops a record whose REF differs, and the
    records it keeps give the count. On GRCh38 none or almost none differ;
    on another build most do. More than 5% different stops here, naming the
    count and the rate (a line drawn, not measured). When no record is left
    to compare, or bcftools norm fails, the task log says so and the run
    goes on.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VCF_PRECHECK {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path reference
    path reference_fai

    output:
    tuple val(meta), env('FILTER_STATUS'), env('FILTER_COUNTS'),
                     env('CONTIG_STYLE'), env('CONTIGS_SEEN'), env('GVCF'),  emit: status
    tuple val(meta), path("${meta.id}.unfiltered_as_pass.vcf.gz"),
                     path("${meta.id}.unfiltered_as_pass.vcf.gz.tbi"), emit: relaxed, optional: true
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def allow_unfiltered = params.allow_unfiltered ? 'true' : 'false'
    """
    # Contigs that hold records, from the index (no pass over the file)
    bcftools index -s ${vcf} | awk '\$3 > 0 { print \$1 }' > contigs.txt || true
    CONTIGS_SEEN=\$(awk 'NR <= 5' contigs.txt | paste -sd, -)
    if grep -q '^chr' contigs.txt; then
        CONTIG_STYLE=chr
    elif [ -s contigs.txt ]; then
        CONTIG_STYLE=other
    else
        CONTIG_STYLE=unknown
    fi

    # gVCF: by name (PharmCAT's own rule), by header, or by reference blocks
    GVCF=none
    case "${vcf.name}" in
        *.g.vcf*|*.genomic.vcf*) GVCF=name ;;
    esac
    bcftools view -h ${vcf} > header.txt

    # One sample per row: the #CHROM columns after FORMAT are the samples
    awk -F'\\t' '/^#CHROM/ { for (i = 10; i <= NF; i++) print \$i }' header.txt > samples.txt
    N_SAMPLES=\$(wc -l < samples.txt | tr -d ' ')
    if [ "\${N_SAMPLES}" -gt 1 ]; then
        NAMES=\$(awk 'NR <= 5' samples.txt | paste -sd, - | sed 's/,/, /g')
        if [ "\${N_SAMPLES}" -gt 5 ]; then NAMES="\${NAMES}, ..."; fi
        # Shell-quoted: a sample name is free text, and the command is for copying
        FIRST=\$(awk 'NR == 1' samples.txt)
        FIRST_Q=\$(printf '%q' "\${FIRST}")
        VCF_Q=\$(printf '%q' "${vcf.name}")
        {
            echo "ERROR: Sample '${meta.id}': ${vcf.name} holds \${N_SAMPLES} samples (\${NAMES})."
            echo "One samplesheet row is one sample, and the steps would mix their genotypes. Keep one"
            echo "sample's column, index the new file, and give each sample its own row, for example:"
            echo "    bcftools view -s \${FIRST_Q} -a -c 1 -Oz -o \${FIRST_Q}.vcf.gz \${VCF_Q}"
            echo "    bcftools index -t \${FIRST_Q}.vcf.gz"
            echo "(-a drops the ALT alleles that sample does not carry, -c 1 the sites where it carries none:"
            echo "the result is a variant-only VCF. Split a gvcf for that row with -s alone, so its reference"
            echo "blocks stay.)"
        } >&2
        exit 1
    fi

    # Build: every ##contig length of chr1-22, X and Y must be GRCh38's
    awk -v grch38='1:248956422 2:242193529 3:198295559 4:190214555 5:181538259 6:170805979 7:159345973 8:145138636 9:138394717 10:133797422 11:135086622 12:133275309 13:114364328 14:107043718 15:101991189 16:90338345 17:83257441 18:80373285 19:58617616 20:64444167 21:46709983 22:50818468 X:156040895 Y:57227415' '
        BEGIN { n = split(grch38, a, " "); for (i = 1; i <= n; i++) { split(a[i], kv, ":"); want[kv[1]] = kv[2] } }
        /^##contig=</ {
            id = ""; len = ""
            if (match(\$0, /[<,]ID=[^,>]+/)) id = substr(\$0, RSTART + 4, RLENGTH - 4)
            if (match(\$0, /[<,]length=[0-9]+/)) len = substr(\$0, RSTART + 8, RLENGTH - 8)
            key = id; sub(/^chr/, "", key)
            if ((key in want) && len != "") {
                n_len++
                if (len != want[key]) print id " " len " " want[key]
            }
        }
        END { print n_len + 0 > "build_lengths.txt" }' header.txt > build_mismatch.txt
    if [ -s build_mismatch.txt ]; then
        N_WRONG=\$(wc -l < build_mismatch.txt | tr -d ' ')
        {
            echo "ERROR: Sample '${meta.id}': ${vcf.name} is not on GRCh38. Its header gives \${N_WRONG} contig length(s) that"
            echo "GRCh38 does not have:"
            awk 'NR <= 5 { printf "    %s length %s (GRCh38: %s)\\n", \$1, \$2, \$3 }' build_mismatch.txt
            if grep -Eq '^(chr)?1 249250621 ' build_mismatch.txt; then
                echo "A chr1 length of 249250621 is GRCh37 (hg19)."
            fi
            echo "The pipeline needs GRCh38: every position would be read on the wrong build. Renaming the"
            echo "contigs does not fix that. Get a GRCh38 VCF from your provider, or start from the FASTQ"
            echo "so the pipeline aligns and calls the reads on GRCh38 (docs/vendor-guide.md, Genome Build)."
        } >&2
        exit 1
    fi

    # Build, when the header gives no length for chr1-22, X or Y: compare the
    # REF bases of the first 1000 records with the reference, leaving out
    # those on a contig the reference lacks (bcftools norm stops on one).
    # bcftools head stops reading by itself, so nothing dies of SIGPIPE. The
    # spot VCF declares the reference's contigs, which bcftools norm needs to
    # write a record.
    if [ "\$(cat build_lengths.txt)" -eq 0 ]; then
        bcftools head -h 0 -n 1000 ${vcf} \\
            | awk -F'\\t' 'NR == FNR { ref[\$1] = 1; next } (\$1 in ref)' ${reference_fai} - \\
            > spot_records.txt
        awk -F'\\t' 'NR == FNR { c = c "##contig=<ID=" \$1 ",length=" \$2 ">\\n"; next }
                     /^##contig=/ { next } /^#CHROM/ { printf "%s", c } { print }' ${reference_fai} header.txt > spot.vcf
        cat spot_records.txt >> spot.vcf
        N_SPOT=\$(wc -l < spot_records.txt | tr -d ' ')
        if [ "\${N_SPOT}" -eq 0 ]; then
            echo "${meta.id}: the header gives no contig lengths and no record is on a reference contig: REF bases not compared"
        elif bcftools norm -c x -f ${reference} -o spot_kept.vcf spot.vcf 2> spot_norm.log; then
            N_KEPT=\$(grep -vc '^#' spot_kept.vcf || true)
            N_DIFF=\$((N_SPOT - N_KEPT))
            PCT=\$(awk -v d="\${N_DIFF}" -v n="\${N_SPOT}" 'BEGIN { printf "%.1f", 100 * d / n }')
            echo "${meta.id}: the header gives no contig lengths; REF bases differ from the reference in \${N_DIFF} of the first \${N_SPOT} records (\${PCT}%)"
            if [ \$((N_DIFF * 20)) -gt "\${N_SPOT}" ]; then
                {
                    echo "ERROR: Sample '${meta.id}': ${vcf.name} is not on GRCh38. Its header gives no contig lengths, so its"
                    echo "REF bases were compared with the reference: they differ in \${N_DIFF} of the first \${N_SPOT} records (\${PCT}%)."
                    echo "On GRCh38 none or almost none differ; more than 5% stops the run. The file is on another build"
                    echo "(GRCh37, hg19, for example), and every position would be read on the wrong one. Renaming the"
                    echo "contigs does not fix that. Get a GRCh38 VCF from your provider, or start from the FASTQ so the"
                    echo "pipeline aligns and calls the reads on GRCh38 (docs/vendor-guide.md, Genome Build)."
                } >&2
                exit 1
            fi
        else
            echo "WARNING: ${meta.id}: could not compare the REF bases of ${vcf.name} with the reference (bcftools norm failed):" >&2
            head -n 5 spot_norm.log >&2
        fi
    fi

    if grep -q '^##GVCFBlock' header.txt; then
        GVCF=blocks
    fi
    # INFO/END can be read only when the header defines it
    QUERY='%FILTER\\n'
    if grep -q '^##INFO=<ID=END,' header.txt; then
        QUERY='%FILTER\\t%ALT\\t%INFO/END\\n'
    fi

    # One pass: FILTER counts, and reference-only records that carry END
    bcftools query -f "\${QUERY}" ${vcf} \\
        | awk -F'\\t' '{ if (\$1 == "PASS") p++; else if (\$1 == ".") d++; else o++ }
                      NF >= 3 && \$3 != "." && (\$2 == "<*>" || \$2 == "<NON_REF>" || \$2 == ".") { b++ }
                      END { printf "%d %d %d %d\\n", p, d, o, b }' \\
        > filter_counts.txt
    read -r N_PASS N_DOT N_OTHER N_BLOCKS < filter_counts.txt
    FILTER_COUNTS="PASS=\${N_PASS} .=\${N_DOT} other=\${N_OTHER}"
    echo "${meta.id}: FILTER counts \${FILTER_COUNTS}"
    if [ "\${N_BLOCKS}" -gt 0 ]; then
        GVCF=blocks
    fi
    echo "${meta.id}: contigs \${CONTIG_STYLE} (\${CONTIGS_SEEN}), gVCF \${GVCF}, reference-block records \${N_BLOCKS}"

    if [ "\${N_PASS}" -gt 0 ]; then
        FILTER_STATUS=pass
    elif [ "${allow_unfiltered}" = "true" ] && [ "\${N_DOT}" -gt 0 ]; then
        FILTER_STATUS=unfiltered
        echo "WARNING: ${meta.id} has no PASS record; --allow_unfiltered set, treating FILTER '.' as PASS" >&2
        bcftools view ${vcf} \\
            | awk -F'\\t' -v OFS='\\t' '/^#/ { print; next } \$7 == "." { \$7 = "PASS" } { print }' \\
            | bcftools view -Oz -o ${meta.id}.unfiltered_as_pass.vcf.gz
        bcftools index -t ${meta.id}.unfiltered_as_pass.vcf.gz
    else
        FILTER_STATUS=no_pass
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    FILTER_STATUS=pass
    FILTER_COUNTS="stub"
    CONTIG_STYLE=chr
    CONTIGS_SEEN=stub
    GVCF=none

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
