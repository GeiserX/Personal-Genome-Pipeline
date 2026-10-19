/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTSV — Annotate structural variants with ACMG pathogenicity classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Classifies each SV into ACMG class 1-5 and adds gene/disease annotations.

    The biocontainer holds AnnotSV's code but none of its annotation data, so
    the annotation directory (setup.sh's annotsv_annotations/, which holds
    Annotations_Human/) comes in through --annotsv_annotations. setup.sh adds
    gnomAD-SV's common SVs to its benign set and mane_select_refseq.txt, the
    MANE Select transcripts passed as -txFile; the task stops without them.

    Before AnnotSV, bin/sv_svlen_fix.py gives sequence-resolved records their
    net length change (SVLEN, and for a net deletion END = POS-1+len(REF)-
    len(ALT), one base less than AnnotSV's own default END, so its POS-1..END
    interval, which counts the anchor base, is as long as the change) and bin/cnvpytor_to_vcf.py adds the Delly (PASS) and
    CNVpytor (e-val < 0.01, q0 < 0.5) deletions and duplications no Manta
    call covers, each record with its caller in INFO SVCALLER; either input
    may be [] when that caller did not run. -bracketedAltMode DEL,DUP,INS
    keeps a breakend pair MANTA left on purpose (a lone inversion junction of
    1 Mb or more) from being read as a whole-span inversion. The class 4 and
    5 rows are written to their own table, found by column name.

    Equivalent to: scripts/05-annotsv.sh, which reads the same duphold-filtered
    calls when step 15 wrote them, and Manta's calls (with a note) otherwise.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANNOTSV {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/annotsv" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sv_vcf), path(delly_vcf), path(cnvpytor_vcf)  // delly_vcf, cnvpytor_vcf: VCF or []
    path(annotations_dir)

    output:
    tuple val(meta), path("${meta.id}_sv_annotated.tsv"), emit: annotated_tsv
    tuple val(meta), path("${meta.id}_sv_class4_5.tsv"),  emit: class4_5_tsv
    tuple val(meta), path("${meta.id}_annotsv_input.vcf"), emit: input_vcf
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def delly_arg    = delly_vcf    ? "--delly ${delly_vcf}"       : ''
    def cnvpytor_arg = cnvpytor_vcf ? "--cnvpytor ${cnvpytor_vcf}" : ''
    """
    for f in ${annotations_dir}/mane_select_refseq.txt \\
             ${annotations_dir}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed; do
        if [ ! -s "\$f" ]; then
            echo "ERROR: \$f not found: AnnotSV's data lacks the MANE Select list or gnomAD-SV's common SVs." >&2
            echo "  Run ./scripts/setup.sh --annotsv-data <genome_dir> to add them, or see docs/05-annotsv.md." >&2
            exit 1
        fi
    done

    sv_svlen_fix.py ${sv_vcf} \\
        | cnvpytor_to_vcf.py --base - --base-name manta ${delly_arg} ${cnvpytor_arg} \\
        > ${meta.id}_annotsv_input.vcf

    AnnotSV \\
        -SVinputFile ${meta.id}_annotsv_input.vcf \\
        -outputFile ${meta.id}_sv_annotated.tsv \\
        -outputDir . \\
        -genomeBuild GRCh38 \\
        -annotationsDir ${annotations_dir} \\
        -annotationMode both \\
        -txFile ${annotations_dir}/mane_select_refseq.txt \\
        -bracketedAltMode DEL,DUP,INS

    # AnnotSV writes no table when it had no SV to annotate.
    [ -s ${meta.id}_sv_annotated.tsv ] || : > ${meta.id}_sv_annotated.tsv
    awk -F'\\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[\$i] = i
                           if (!("ACMG_class" in c) || !("Annotation_mode" in c)) { print "no ACMG_class or Annotation_mode column" > "/dev/stderr"; exit 2 }
                           print; next }
                 \$c["Annotation_mode"] == "full" && (\$c["ACMG_class"] == "4" || \$c["ACMG_class"] == "5")' \\
        ${meta.id}_sv_annotated.tsv > ${meta.id}_sv_class4_5.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_annotated.tsv ${meta.id}_sv_class4_5.tsv ${meta.id}_annotsv_input.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
