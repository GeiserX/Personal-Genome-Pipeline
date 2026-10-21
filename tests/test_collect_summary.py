#!/usr/bin/env python3
"""bin/collect_summary.py and bin/render_report.py on synthetic sample folders.

Checks:
  1. a bash output folder: the summary validates against
     tests/schema/summary.schema.json, and the text and HTML reports show the
     same ClinVar count, PRS rows and QC values;
  2. a Nextflow output folder (roh/, coverage/, hla/, pharmcat/): ROH, depth and
     HLA are read (the old reports looked for ROH under vcf/ only);
  3. heteroplasmy counts only chrM calls with an allele fraction of 0.05 to 0.95;
  4. the CPSR breakdown is read from classification.tsv.gz by column name;
  5. a result older than the latest run-all.sh run whose step was skipped is
     marked stale in the summary and in both reports; one whose step was ok, or
     that is newer than the run, is not;
  6. bin/clinvar_hits.awk and collect_summary.py give every review status the
     same number of stars; 6b. and the same combined-allele note to a combined
     record, its parts and two hits that may be one allele, which both reports
     show; 6c. step 06's cDNA read check reaches each ClinVar hit by chrom,
     pos, ref and alt, and a hit the table lacks is not checked;
  7. a corrupt file (not gzip, or gzip with damaged data) makes its section
     'unreadable' instead of stopping the report;
  8. the secondary-findings tier lists the ClinVar hits and the rare
     HIGH-impact clinical records in ACMG SF v3.3 genes, and nothing else;
  9. haplocheck's contamination status and Yleaf's Y haplogroup (or
     'insufficient markers') reach both reports;
 10. the step 33 verdict with the declared sex as somalier's pedigree:
     somalier keeps the pedigree's sex when the reads cannot tell, which is
     read as unknown (not checked), not as a call that agrees; a sex somalier
     changed from the pedigree's is a mismatch; a call that agrees with
     enough chrX sites is ok; a pedigree female with chrY reads but chrX that
     cannot tell is unknown, not an aneuploidy; without a pedigree the rows
     read as before;
 11. CYP2D6 by step 36's rule: pypgx 'Indeterminate' is not a call and is
     not counted in 'Genes called'; a Cyrius genotype with
     Filter=CYP2D6_depth_unreliable is shown as not usable and is not a
     call; step 36's verdict is shown when its table is there, and ignored
     when the table is older than a caller's result;
 12. repeat expansions: Stranger's STR_STATUS is read, every locus that is
     not normal is listed (also outside the five key loci), a record
     without a status is counted apart, not flagged; a Stranger file older
     than ExpansionHunter's is not read and both reports say to rerun step
     9b; a stale card says Stale, not Complete;
 13. HLA: each allele keeps its own T1K quality; a one-allele row
     ('.', 0, -1) is not low confidence; a quality of 0 or below is, and on
     HLA-A or HLA-B the reports say the gene is not passed to PharmCAT;
 14. a section whose step failed in the latest run says Failed, not Not
     run, in the summary and both reports, also in the CYP2D6 rows;
 15. the CPIC card counts the genes called without a function phenotype;
 16. the mosaic_* rows of the step 33 table (bin/mosaic_sex_loss.py) make
     their own section, with the numbers as numbers; a table without them
     (written before the estimate) leaves the section missing, not broken; a
     mosaic_call that is not one of the four is unreadable.
 17. step 35b (Paraphase): the named genes in order, the other regions
     counted, the variant rows, SMN1 beside step 35's Parascopy call (agree,
     disagree, or no line without that table) in both reports, and a folder
     without the step lists it under Steps Not Run.
 18. step 42's traits table (bin/traits_report.py's own output): the counts
     per status leave the APOE rows out, the rows and the APOE lines reach
     both reports (the APOE lines last, after the warning, folded away in
     the HTML), a row with ancestry yes is marked, a not-covered row is
     listed as not read; a table without the reading column, or with a row
     cut short, is unreadable;
     without the file the section is missing and listed as not run.
 19. step 40's summary table: its counts and '--' reasons are read; a table
     without a whole-number 'called' row, or with a '--' reason without a
     whole-number count, is unreadable, not a zero or a reason left out; an
     older table whose step was skipped in the latest run is stale.
 21. ClinVar call quality and frequency (step 06's floor, step 13's VEP): a
     hit below GQ 30 or DP 10 is marked low confidence and left out of the
     count and the stars, in the summary and in bin/clinvar_hits.awk alike,
     and a floor of 0 marks nothing; each hit carries the gnomAD AF read from
     a bgzipped, tabix-indexed VEP file (two BGZF blocks, two contigs), a hit
     above 1% is common, one VEP lacks is unknown, and without VEP output the
     reports say why; both reports print the AF, the call mark and the
     frequency caveat on every hit;
 22. CPSR beside ClinVar (step 17): the P/LP rows by final or own class, the
     rows where CLINVAR_CLASSIFICATION and CPSR_CLASSIFICATION disagree on
     P/LP with CPSR's ACMG codes, the final P/LP rows CPSR's own rules class
     benign, the rows below the depth floor (DP_CONTROL), the cDNA read
     check of step 17 on the P/LP rows, and CPSR's bundled ClinVar release
     beside step 6's; a table with one class column shows no disagreement;
     bin/clip_check.py --cpsr checks the P/LP rows only (a fake pysam: the
     spliced reads at a P/LP row are flagged, genomic reads at the same row
     are clear, and CPSR's "13" finds the BAM's "chr13").
 23. step 32's summary beyond the diplotypes: pypgx's copy-number calls for
     its BAM genes and UGT2B17's copy number from CNVpytor over its GENCODE
     span (no copy, one copy, two copies with no call; not assessed without
     the gene lines or without CNVpytor) in a PGx copy number block; GSTT1
     not assessable on this reference; a call resting on a GQ 12 defining
     variant flagged; the CYP2D6 depth check with its MAPQ shares; the
     probable CYP2D6 call of pypgx alone in both reports, and none after a
     failed depth check, beside a usable Cyrius call, or from a summary
     written before the evidence columns.
 20. step 39's couple table (bin/couple_report.py's columns): the partner,
     the shared and X-linked genes with each side's findings from this
     sample's point of view, the one-partner count; a table naming another
     couple, or without a column, is unreadable; a table without rows names
     no partner; a single sample without the step has no "not run" line for
     it (the screen needs two people).
 24. the clinical filter's HIGH-impact list (bin/csq_rank.awk's rule, step 23):
     a HIGH consequence on a non-MANE transcript whose MANE Select one is
     intronic is not listed; without a MANE transcript the canonical one
     ranks; a heterozygote at GQ 12 with AD 3,2 and a 1/1 at GQ 11 with AD
     12,5 are weak calls, a 1/1 at GQ 50 with AD 0,30 is not, and a floor of 0
     in ##PGP_CALL_FLOOR marks nothing; an indel with no MAX_AF is listed as
     frequency unknown, a common one not at all; neither a weak call nor an
     unknown frequency counts as rare; the summary table's FREQ_TIER and
     CALL_CHECK columns and step 31's FREQ_TIER are counted; both reports
     show the counts and the marks. bin/csq_rank.awk, the program steps 23
     and 31 and their modules run, tags the same records the same way
     (PGP_IMPACT, PGP_TX, PGP_FREQ, PGP_CALL), drops the common one with
     keep=noncommon, and agrees with the summary on every record.

Run: python3 tests/test_collect_summary.py
"""
import gzip
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import types
import zlib

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests", "schema"))
import clip_check  # noqa: E402
import collect_summary  # noqa: E402
import render_report  # noqa: E402
import validate  # noqa: E402

FAILS = 0
DAY = 86400


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:400]}")
    if not ok:
        FAILS += 1


def put(path, text, gz=False, age_days=0):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        with open(path, "w") as f:
            f.write(text)
    if age_days:
        t = time.time() - age_days * DAY
        os.utime(path, (t, t))


VCF_HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
HITS = VCF_HEAD + (
    "chr1\t10\t1\tC\tT\t50\tPASS\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\tGT\t0/1\n"
    "chr1\t20\t2\tT\tG\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Likely_pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT\t1/1\n"
    "chr2\t5\t3\tA\tG\t50\tPASS\tGENEINFO=GENEC:33;CLNSIG=Pathogenic;CLNREVSTAT=no_assertion_criteria_provided\tGT\t0|1\n")
# Step 6 hits with the combined-allele tags of bin/clinvar_join_adjacent.awk:
# a 2-star combined record, its two 1-star parts, a single SNV, and two hits in
# one gene that matched no combined record.
_CV1 = "GENEINFO=SYNTH1:1;CLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided,_single_submitter"
COMBINED_HITS = VCF_HEAD + (
    f"chr22\t30000037\t1\tAGC\tA\t50\tPASS\t{_CV1};PGP_PART_OF=chr22:30000038:GCCAT:CATGGG\tGT\t0/1\n"
    "chr22\t30000038\t3\tGCCAT\tCATGGG\t.\tPASS\tPGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG;"
    "GENEINFO=SYNTH1:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\tGT\t0/1\n"
    f"chr22\t30000042\t2\tT\tTGGG\t50\tPASS\t{_CV1};PGP_PART_OF=chr22:30000038:GCCAT:CATGGG\tGT\t0/1\n"
    "chr22\t30000210\t4\tC\tT\t50\tPASS\tGENEINFO=SYNTH2:2;CLNSIG=Pathogenic;"
    "CLNREVSTAT=criteria_provided,_single_submitter\tGT\t0/1\n"
    f"chr22\t30000300\t5\tA\tG\t50\tPASS\t{_CV1.replace('SYNTH1', 'SYNTH3')};PGP_ONE_ALLELE=chr22:30000305:C:T\tGT\t0/1\n"
    f"chr22\t30000305\t6\tC\tT\t50\tPASS\t{_CV1.replace('SYNTH1', 'SYNTH3')};PGP_ONE_ALLELE=chr22:30000300:A:G\tGT\t0/1\n")
MITO = VCF_HEAD.replace("##fileformat=VCFv4.2\n", "") + "".join(
    f"chrM\t{i}\t.\tA\tG\t.\tPASS\t.\tGT:AF\t0/1:{af}\n" for i, af in enumerate(("0.02", "0.30", "0.99", "0.05", "0.949"), 1))


def bash_folder(d, s="S"):
    put(f"{d}/vcf/{s}.vcf.gz", VCF_HEAD + "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "chr1\t2\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1\nchr1\t3\t.\tC\tT\t5\tRefCall\t.\tGT\t0/0\n", gz=True)
    put(f"{d}/clinvar/{s}_clinvar_hits.vcf", HITS)
    put(f"{d}/mito/{s}_chrM_filtered.vcf.gz", MITO, gz=True)
    put(f"{d}/vcf/{s}_roh.txt", "# RG\tsample\tchrom\tstart\tend\tlength\n"
        "RG\tS\tchr1\t100\t7000100\t7000000\t3\t50\nRG\tS\tchrX\t1\t9000000\t9000000\t3\t50\n")
    put(f"{d}/mosdepth/{s}.mosdepth.summary.txt",
        "chrom\tlength\tbases\tmean\tmin\tmax\nchr1\t10\t300\t30.0\t0\t50\ntotal\t10\t300\t29.5\t0\t50\n")
    put(f"{d}/indexcov/indexcov-indexcov.ped", "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\t"
        "CNchrX\tCNchrY\nS\tS\t-9\t-9\t2\t-9\t2.0\t0.0\n")
    put(f"{d}/prs/{s}_prs_summary.tsv", "Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\n"
        "Coronary artery disease\tPGS000018\t0.12\t100\t200\nType 2 diabetes\tPGS000014\t-0.4\t50\t90\n")
    put(f"{d}/cpsr/{s}.cpsr.grch38.classification.tsv.gz",
        "SAMPLE_ID\tGENOMIC_CHANGE\tCPSR_CLASSIFICATION\tCLASSIFICATION\n"
        "S\ta\tVUS\tVUS\nS\tb\tVUS\tLikely_Pathogenic\nS\tc\tB\tVUS\n", gz=True)
    put(f"{d}/slivar/{s}_slivar_summary.tsv", "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\nchr1\t1\tA\tG\tHIGH\tX\n")
    put(f"{d}/hla_t1k/{s}_hla_genotype.tsv", "HLA-A\t2\tA*01:01:01\t30\t60\tA*02:01:01\t20\t55\t\n")
    put(f"{d}/cpic/{s}_phenotypes.tsv", "Gene\tDiplotype\tPhenotype\tStatus\n"
        "CYP2D6\t*1/*4\tIntermediate Metabolizer\tnon-normal\nCYP2C19\t*1/*1\tNormal Metabolizer\tnormal\n"
        "VKORC1\trs9923231 reference (C)/rs9923231 reference (C)\t-1639 GG\tunclassified\n")
    put(f"{d}/cpic/{s}_cpic_recommendations.txt", "Pharmacogenomic Drug Recommendations\n")
    put(f"{d}/pypgx/{s}_pypgx_summary.tsv", "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\nCYP2D6\t*4/*1\tIM\t.\tbam\n")
    put(f"{d}/cyrius/{s}_cyp2d6.tsv", "Sample\tGenotype\tFilter\nS\t*1/*4\tPASS\n")


def render(d, s="S"):
    summ = collect_summary.collect(s, d)
    return summ, render_report.text_report(summ), render_report.html_report(summ)


def bgzf_block(data):
    """One BGZF block (a gzip member with the BC extra field) holding data."""
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    comp = c.compress(data) + c.flush()
    return (struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, len(comp) + 25) + comp
            + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data)))


def write_tabixed_vcf(path, header, blocks):
    """A bgzipped VCF of one BGZF block per item of blocks (lists of record
    lines, the header in the first block) and its .tbi with the linear index
    only (no bins: bin/collect_summary.py reads the linear index)."""
    raw, lin, coff = b"", {}, 0
    for i, lines in enumerate(blocks):
        data = (header if i == 0 else "").encode()
        for line in lines:
            chrom, pos = line.split("\t")[0], int(line.split("\t")[1])
            w, ix = (pos - 1) >> 14, lin.setdefault(chrom, [])
            while len(ix) <= w:
                ix.append(None)
            if ix[w] is None:
                ix[w] = (coff << 16) | len(data)
            data += (line + "\n").encode()
        blk = bgzf_block(data)
        raw += blk
        coff += len(blk)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(raw + bgzf_block(b""))
    names = list(lin)
    nm = b"".join(n.encode() + b"\0" for n in names)
    tbi = b"TBI\x01" + struct.pack("<8i", len(names), 2, 1, 2, 0, ord("#"), 0, len(nm)) + nm
    for n in names:
        ix, last = lin[n], 0
        for j, v in enumerate(ix):  # an empty window takes the offset before it, as htslib fills it
            ix[j] = last = v if v is not None else last
        tbi += struct.pack("<i", 0) + struct.pack("<i", len(ix)) + struct.pack(f"<{len(ix)}Q", *ix)
    with gzip.open(path + ".tbi", "wb") as f:
        f.write(tbi)
    return {n: tuple(v) for n, v in lin.items()}


CSQ_HEAD = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
            'Format: Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF|gnomADg_AF">\n')
FLOOR_HITS = ("##fileformat=VCFv4.2\n##PGP_CALL_FLOOR=GQ=30,DP=10\n" + VCF_HEAD.split("\n", 1)[1] + (
    "chr1\t100\t1\tC\tT\t50\tPASS\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel"
    "\tGT:GQ:DP\t0/1:12:6\n"
    "chr1\t40000\t2\tG\tA\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel"
    "\tGT:GQ:DP\t0/1:50:40\n"
    "chr2\t20000\t3\tCA\tC\t50\tPASS\tGENEINFO=GENEC:33;CLNSIG=Likely_pathogenic;"
    "CLNREVSTAT=criteria_provided,_single_submitter\tGT:DP:GQ\t0/1:35:45\n"
    "chr3\t10\t4\tA\tG\t50\tPASS\tGENEINFO=GENED:44;CLNSIG=Pathogenic;"
    "CLNREVSTAT=criteria_provided,_single_submitter\tGT\t1/1\n"))


def clinvar_quality(work, schema):
    """21. the call floor and the gnomAD frequency on step 06's hits."""
    d = os.path.join(work, "clinvar-quality", "S")
    put(f"{d}/clinvar/S_clinvar_hits.vcf", FLOOR_HITS)
    vep_head = "##fileformat=VCFv4.2\n" + CSQ_HEAD + VCF_HEAD.split("\n", 1)[1]
    lin = write_tabixed_vcf(f"{d}/vep/S_vep.vcf.gz", vep_head, [
        ["chr1\t100\t.\tC\tT\t50\tPASS\tCSQ=T|missense_variant|MODERATE|GENEA|0.03|0.02\tGT\t0/1",
         "chr1\t40000\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|GENEB|0.0002|\tGT\t0/1"],
        ["chr2\t50\t.\tT\tG\t50\tPASS\tCSQ=G|intron_variant|MODIFIER|GENEC|0.5|0.5\tGT\t0/1",
         "chr2\t20000\t.\tCA\tC\t50\tPASS\tCSQ=-|frameshift_variant|HIGH|GENEC|0.0001|0.00005\tGT\t0/1"]])
    summ, txt, html = render(d)
    cv = summ["sections"]["clinvar"]["data"]
    by = {h["gene"]: h for h in cv["hits"]}
    check("call floor: the GQ 12, DP 6 hit is low confidence, the others are not (no GQ/DP is not low)",
          [g for g in ("GENEA", "GENEB", "GENEC", "GENED") if by[g]["low_confidence"]] == ["GENEA"]
          and (by["GENEA"]["gq"], by["GENEA"]["dp"]) == (12, 6) and (by["GENEC"]["gq"], by["GENEC"]["dp"]) == (45, 35)
          and (by["GENED"]["gq"], by["GENED"]["dp"]) == (None, None), by)
    check("call floor: the low-confidence hit is listed but left out of the count and the stars",
          cv["count"] == 3 and cv["low_confidence"] == 1 and len(cv["hits"]) == 4
          and cv["by_stars"] == {"4": 0, "3": 1, "2": 0, "1": 2, "0": 0}, cv)
    awk = subprocess.run(["awk", "-f", os.path.join(REPO, "bin", "clinvar_hits.awk"), f"{d}/clinvar/S_clinvar_hits.vcf"],
                         capture_output=True, text=True).stdout
    awk_low = {r.split("\t")[6]: tuple(r.split("\t")[10:13]) for r in awk.splitlines()}
    check("call floor: bin/clinvar_hits.awk gives the same DP, GQ and mark",
          awk_low == {"GENEA": ("6", "12", "yes"), "GENEB": ("40", "50", "no"), "GENEC": ("35", "45", "no"),
                      "GENED": (".", ".", "no")}, awk_low)
    zero = FLOOR_HITS.replace("##PGP_CALL_FLOOR=GQ=30,DP=10", "##PGP_CALL_FLOOR=GQ=0,DP=0")
    z = os.path.join(work, "clinvar-floor-zero", "S")
    put(f"{z}/clinvar/S_clinvar_hits.vcf", zero)
    cz = collect_summary.collect("S", z)["sections"]["clinvar"]["data"]
    awk_z = subprocess.run(["awk", "-f", os.path.join(REPO, "bin", "clinvar_hits.awk"), f"{z}/clinvar/S_clinvar_hits.vcf"],
                           capture_output=True, text=True).stdout
    check("call floor: a floor of GQ 0 and DP 0 marks nothing, so the mark above depends on the floor",
          cz["low_confidence"] == 0 and cz["count"] == 4 and "\tyes" not in awk_z, (cz["low_confidence"], awk_z))
    check("call floor: a hits file without the floor line uses GQ 30 and DP 10",
          collect_summary.call_floor(f"{work}/cdna/S/clinvar/S_clinvar_hits.vcf") == {"GQ": 30, "DP": 10})
    af = {g: (by[g]["af"], by[g]["common"]) for g in by}
    check("gnomAD AF: read through the tabix index from both BGZF blocks; above 1% is common; no VEP record is unknown",
          af == {"GENEA": (0.03, True), "GENEB": (0.0002, False), "GENEC": (0.0001, False), "GENED": (None, False)}
          and (by["GENEA"]["gnomad_exome_af"], by["GENEA"]["gnomad_genome_af"]) == (0.03, 0.02), af)
    got_lin = collect_summary.tabix_linear_index(f"{d}/vep/S_vep.vcf.gz")
    check("tabix: the linear index is read as written (chr2 starts in the second BGZF block)",
          got_lin == lin and got_lin["chr2"][0] >> 16 > 0 and len(got_lin["chr1"]) == 3, (got_lin, lin))
    check("gnomAD AF: the summary names the VEP file", cv["frequency"]["source"] == "vep/S_vep.vcf.gz"
          and cv["frequency"]["note"] is None, cv["frequency"])
    check("both reports: the AF and the common mark on the hit, unknown on the hit VEP lacks",
          "frequency: gnomAD AF 0.03; common: unlikely to be a high-penetrance finding" in txt
          and "frequency: gnomAD AF unknown" in txt
          and "<td>gnomAD AF 0.03; common: unlikely to be a high-penetrance finding</td>" in html
          and "<th>gnomAD frequency</th>" in html, txt)
    check("both reports: the low-confidence mark on the hit and the count beside the headline",
          "call: low-confidence call (GQ 12, DP 6): check the reads" in txt
          and "Low-confidence calls (GQ below 30 or DP below 10), listed but not counted above: 1" in txt
          and "<td>low-confidence call (GQ 12, DP 6): check the reads</td>" in html, txt)
    check("both reports: the frequency caveat covers heterozygous hits too",
          all("a heterozygous hit in a dominant gene as well as a homozygous one" in t for t in (txt, html)))
    errs = validate.validate(schema, json.loads(json.dumps(summ)))
    check("call floor and AF: the summary validates", not errs, errs)
    os.remove(f"{d}/vep/S_vep.vcf.gz.tbi")
    cn = collect_summary.collect("S", d)["sections"]["clinvar"]["data"]
    check("gnomAD AF: a VEP file without its .tbi leaves every AF unknown and says why",
          all(h["af"] is None for h in cn["hits"]) and "no .tbi" in cn["frequency"]["note"], cn["frequency"])
    shutil.rmtree(f"{d}/vep")
    summ_n, txt_n, _ = render(d)
    check("gnomAD AF: without step 13 the text report says the frequency is unknown and why",
          "gnomAD AF: unknown for every hit (no VEP output (step 13) in this folder)" in txt_n, txt_n)


CPSR_COLS = ("VAR_ID\tGENOMIC_CHANGE\tSYMBOL\tCLINVAR_CLASSIFICATION\tCPSR_CLASSIFICATION\tACMG_CODE\t"
             "ASSERTION_AUTHORITY\tCLASSIFICATION\tDP_CONTROL\n")
CPSR_ROWS = (
    "13_32340300_G_A\t13:g.32340300G>A\tBRCA2\tPathogenic\tBenign\tBA1\tClinVar\tPathogenic\t40\n"
    "17_43045700_C_T\t17:g.43045700C>T\tBRCA1\tLikely Pathogenic\tLikely Benign\tBS1\tClinVar\tLikely Pathogenic\t6\n"
    "2_47800000_A_G\t2:g.47800000A>G\tMSH2\t\tLikely Pathogenic\tPVS1,PM2\tCPSR\tLikely Pathogenic\t30\n"
    "5_112800000_T_C\t5:g.112800000T>C\tAPC\tBenign\tPathogenic\tPVS1,PM2\tClinVar\tBenign\t30\n"
    "3_100_A_G\t3:g.100A>G\tMLH1\tVUS\tVUS\tPM2\tCPSR\tVUS\t30\n")


class FakeRead:
    """The pysam.AlignedSegment attributes bin/clip_check.py reads."""

    def __init__(self, name, start, cigar, seq):
        self.query_name, self.reference_start, self.cigartuples, self.query_sequence = name, start, cigar, seq
        self.is_unmapped = self.is_secondary = self.is_duplicate = self.is_qcfail = self.is_read2 = False
        self.is_read1 = True

    def has_tag(self, tag):
        return False


def fake_pysam(pos0, exon_end):
    """A pysam stand-in with two BAMs on contig chr13: 'spliced' has 8 alt
    reads clipped at the exon end, 'genomic' 8 alt reads running on; both
    have 6 genomic ref reads."""
    def reads(kind):
        out = [FakeRead(f"ref{i}", pos0 - 50 + i, [(0, 100)], "N" * (50 - i) + "G" + "N" * (49 + i)) for i in range(6)]
        for i in range(8):
            if kind == "spliced":
                aligned = 70 + i
                s = exon_end - aligned
                out.append(FakeRead(f"alt{i}", s, [(0, aligned), (4, 100 - aligned)],
                                    "N" * (pos0 - s) + "A" + "N" * (100 - (pos0 - s) - 1)))
            else:
                out.append(FakeRead(f"alt{i}", pos0 - 50 + i, [(0, 100)], "N" * (50 - i) + "A" + "N" * (49 + i)))
        return out

    class AlignmentFile:
        references = ("chr13",)

        def __init__(self, path):
            self.kind = os.path.basename(path)

        def fetch(self, contig, start, end):
            return reads(self.kind) if contig == "chr13" else []
    return types.SimpleNamespace(AlignmentFile=AlignmentFile)


def cpsr_beside_clinvar(work, schema):
    """22. CPSR's own class beside ClinVar's, its depth, its read check and release."""
    d = os.path.join(work, "cpsr-compare", "S")
    put(f"{d}/cpsr/S.cpsr.grch38.classification.tsv.gz", CPSR_COLS + CPSR_ROWS, gz=True)
    put(f"{d}/cpsr/S.cpsr.grch38.html", "<html></html>\n")
    put(f"{d}/clinvar/S_clinvar_hits.vcf", FLOOR_HITS)
    put(f"{d}/cpsr/S_cpsr_cdna_check.tsv", "\t".join(clip_check.COLUMNS) + "\n"
        "13\t32340300\tG\tA\t12\t0\t10\t0\tclear\t\n"
        f"2\t47800000\tA\tG\t9\t8\t11\t0\tflagged\t{clip_check.FLAG_NOTE}\n")
    put(f"{d}/run_manifest.tsv", "data\tclinvar_file_date\t20261001\ndata\tcpsr_clinvar_release\t20250309\n")
    summ, txt, html = render(d)
    c = summ["sections"]["cpsr"]["data"]
    plp = [x["gene"] for x in c["plp"]]
    dis = {x["gene"]: (x["clinvar"], x["cpsr"], x["acmg_codes"]) for x in c["disagreements"]}
    check("CPSR: the P/LP rows by final or own class, not the VUS row", plp == ["BRCA2", "BRCA1", "MSH2", "APC"], plp)
    check("CPSR: both directions of a ClinVar and CPSR disagreement, with CPSR's codes",
          dis == {"BRCA2": ("Pathogenic", "Benign", "BA1"), "BRCA1": ("Likely Pathogenic", "Likely Benign", "BS1"),
                  "APC": ("Benign", "Pathogenic", "PVS1,PM2")}, dis)
    check("CPSR: two final P/LP rows that CPSR's own rules class benign or likely benign",
          c["final_plp_own_benign"] == 2 and c["classification"] == {"Pathogenic": 1, "Likely Pathogenic": 2,
                                                                     "Benign": 1, "VUS": 1}, c)
    check("CPSR: the DP 6 row is low confidence against step 6's depth floor",
          c["low_confidence"] == 1 and [x["gene"] for x in c["plp"] if x["low_confidence"]] == ["BRCA1"], c["plp"])
    cd = {x["gene"]: x["cdna"]["check"] for x in c["plp"]}
    check("CPSR: step 17's read check on each P/LP row, a row the table lacks not checked",
          cd == {"BRCA2": "clear", "BRCA1": "not_checked", "MSH2": "flagged", "APC": "not_checked"}
          and c["cdna_flagged"] == 1, cd)
    check("CPSR: both reports print the disagreement with both classes and the codes",
          all("BRCA2 13:32340300 G>A: final Pathogenic (ClinVar), ClinVar Pathogenic, CPSR's own Benign (BA1)"
              in t.replace("&#x27;", "'").replace("&gt;", ">") for t in (txt, html)), txt)
    check("CPSR: both reports print the two ClinVar releases",
          all("ClinVar release in CPSR's bundle: 20250309; step 6 used the ClinVar file of 20261001"
              in t.replace("&#x27;", "'") for t in (txt, html))
          and summ["databases"]["cpsr_clinvar_release"] == "20250309", txt)
    check("CPSR: both reports print the low-confidence row and the cDNA flag",
          "Rows where ClinVar and CPSR's own class disagree on P/LP: 3" in txt
          and "low-confidence call (DP 6): check the reads" in txt
          and "cDNA read check of the P/LP rows (alt reads clipped at exon edges): 1 flagged" in txt, txt)
    errs = validate.validate(schema, json.loads(json.dumps(summ)))
    check("CPSR comparison: the summary validates", not errs, errs)
    # Red control: one class column, nothing to compare.
    put(f"{d}/cpsr/S.cpsr.grch38.classification.tsv.gz",
        "VAR_ID\tSYMBOL\tCLASSIFICATION\n13_32340300_G_A\tBRCA2\tPathogenic\n5_112800000_T_C\tAPC\tBenign\n", gz=True)
    summ1, txt1, _ = render(d)
    c1 = summ1["sections"]["cpsr"]["data"]
    check("CPSR: a table with one class column has no disagreement and says so",
          c1["disagreements"] == [] and c1["own_class_column"] is None
          and "has no CPSR_CLASSIFICATION column" in txt1, c1)

    # bin/clip_check.py --cpsr, with a fake pysam
    pos0 = 32340299
    exon_end = pos0 + 30
    put(f"{d}/exons.bed", f"chr13\t{pos0 - 30}\t{exon_end}\nchr13\t{pos0 + 3000}\t{pos0 + 3150}\n")
    put(f"{d}/cls.tsv.gz", CPSR_COLS + CPSR_ROWS, gz=True)
    keep = sys.modules.get("pysam")
    sys.modules["pysam"] = fake_pysam(pos0, exon_end)
    try:
        got = {}
        for kind in ("spliced", "genomic"):
            out = f"{d}/{kind}.tsv"
            with open(os.devnull, "w") as null:
                so, sys.stdout = sys.stdout, null
                try:
                    rc = clip_check.main(["--cpsr", f"{d}/cls.tsv.gz", "--bam", f"{d}/{kind}", "--exons",
                                          f"{d}/exons.bed", "--out", out])
                finally:
                    sys.stdout = so
            rows = [line.rstrip("\n").split("\t") for line in open(out)][1:]
            got[kind] = (rc, [(r[0], r[1], r[8]) for r in rows])
    finally:
        if keep is None:
            sys.modules.pop("pysam", None)
        else:
            sys.modules["pysam"] = keep
    check("clip_check --cpsr: the spliced reads at a P/LP row are flagged, CPSR's contig name kept",
          got["spliced"][0] == 0 and got["spliced"][1][0] == ("13", "32340300", "flagged"), got)
    check("clip_check --cpsr: genomic reads at the same row are clear",
          got["genomic"][1][0] == ("13", "32340300", "clear"), got)
    check("clip_check --cpsr: only the four P/LP rows are checked (not the VUS row)",
          [r[0] for r in got["spliced"][1]] == ["13", "17", "2", "5"], got)


PYPGX_HEAD = ("Gene\tDiplotype\tPhenotype\tCNV_call\tSource\tStatus\tDefining_variants\tMin_GQ\tWeakest_variant"
              "\tDefining_filter\tQuality\n")
PYPGX_ROWS = (
    "CYP2D6\t*1/*2\tNormal Metabolizer\tNormal\tBAM\tcalled\t1/1 in the VCF\t50\tchr22:42126611 C>G GQ 50 AD 20,22"
    "\tPASS\tok\n"
    "CYP2A6\t*1/*1\tno pypgx phenotype rule\tNormal\tBAM\tcalled\tnone\tN/A\tN/A\tnone\tN/A\n"
    "GSTM1\t*0/*0\tno pypgx phenotype rule\tWholeDel2\tBAM\tcalled\tnone\tN/A\tN/A\tnone\tN/A\n"
    "GSTT1\tN/A\tN/A\tN/A\tBAM\tnot callable on this reference\tN/A\tN/A\tN/A\tN/A\tN/A\n"
    "CYP2C19\t*1/*2\tIntermediate Metabolizer\tN/A\tVCF\tcalled\t2/2 in the VCF\t12"
    "\tchr10:94781859 G>A GQ 12 AD 26,3\tPASS\tlow-confidence: check the reads\n")
DEPTH_OK = ("metric\tvalue\nstatus\tok\nmessage\tCYP2D6 depth is fit for a copy-number call\ngene_depth_all\t18.53\n"
            "gene_mapq1_share\t0.91\ngene_mapq20_share\t0.84\nbody_depth_all\t30.10\nbody_mapq1_share\t0.99\n"
            "body_mapq20_share\t0.97\n")
GENES_GTF = ('chr4\tHAVANA\tgene\t68516301\t68627039\t.\t-\t.\tgene_id "ENSG00000197888.11"; gene_type "protein_coding"; '
             'gene_name "UGT2B17"; level 2;\n'
             'chr4\tHAVANA\tgene\t68670652\t68694000\t.\t-\t.\tgene_id "ENSG1"; gene_name "UGT2B15"; level 2;\n')


def clinical_mane_quality(work, schema):
    """24. the clinical filter: MANE first, weak calls, unknown frequency."""
    d = os.path.join(work, "clinical-mane", "S")
    fmt = "Allele|Consequence|IMPACT|SYMBOL|Feature|CANONICAL|MANE_SELECT|MAX_AF"
    head = ("##fileformat=VCFv4.2\n"
            f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {fmt}">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    body = (
        # TTN: HIGH on a non-MANE transcript, intronic on the MANE Select one
        "chr2\t100\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|TTN|ENST2|||0.0001,"
        "T|intron_variant|MODIFIER|TTN|ENST1|YES|NM_1|0.0001\tGT:GQ:DP:AD\t0/1:50:30:15,15\n"
        # BRCA2: a heterozygote at GQ 12, AD 3,2
        "chr13\t200\t.\tA\tAT\t50\tPASS\tCSQ=AT|frameshift_variant|HIGH|BRCA2|ENST3|YES|NM_2|0.0001"
        "\tGT:GQ:DP:AD\t0/1:12:.:3,2\n"
        # LDLR: 1/1 at GQ 11 on 12 ref and 5 alt reads
        "chr19\t300\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|LDLR|ENST4|YES|NM_3|0.0001"
        "\tGT:GQ:DP:AD\t1/1:11:17:12,5\n"
        # MYBPC3: 1/1 at GQ 50 on 30 alt reads
        "chr11\t400\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|MYBPC3|ENST5|YES|NM_4|0.0001"
        "\tGT:GQ:DP:AD\t1/1:50:30:0,30\n"
        # MYH7: a multi-base deletion VEP found no frequency for
        "chr14\t500\t.\tCTTT\tC\t50\tPASS\tCSQ=-|frameshift_variant|HIGH|MYH7|ENST6|YES|NM_5|"
        "\tGT:GQ:DP:AD\t0/1:50:30:15,15\n"
        # KCNQ1: a common stop-gain (in the clinical VCF through the ClinVar tier)
        "chr11\t600\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|KCNQ1|ENST7|YES|NM_6|0.3"
        "\tGT:GQ:DP:AD\t0/1:50:30:15,15\n"
        # SCN5A: no MANE transcript; HIGH on the canonical one
        "chr3\t700\t.\tC\tT\t50\tPASS\tCSQ=T|missense_variant|MODERATE|SCN5A|ENST8|||0.0001,"
        "T|stop_gained|HIGH|SCN5A|ENST9|YES||0.0001\tGT:GQ:DP:AD\t0/1:50:30:15,15\n")
    put(f"{d}/clinical/S_clinical.vcf.gz", head + body, gz=True)
    # bin/csq_rank.awk on the same records
    awk_file = os.path.join(REPO, "bin", "csq_rank.awk")

    def ranked(keep):
        p_ = subprocess.run(["awk", "-v", f"keep={keep}", "-f", awk_file], input=head + body,
                            capture_output=True, text=True, env=dict(os.environ, LC_ALL="C"))
        tags = {}
        for line in p_.stdout.splitlines():
            if not line.startswith("#"):
                r = line.split("\t")
                tags[int(r[1])] = dict(kv.split("=", 1) for kv in r[7].split(";") if kv.startswith("PGP_"))
        return p_.returncode, tags, p_.stdout
    rc, tags, out = ranked("all")
    want = {100: ("MODIFIER", "MANE", "rare", "ok"), 200: ("HIGH", "MANE", "rare", "weak"),
            300: ("HIGH", "MANE", "rare", "weak"), 400: ("HIGH", "MANE", "rare", "ok"),
            500: ("HIGH", "MANE", "unknown", "ok"), 600: ("HIGH", "MANE", "common", "ok"),
            700: ("HIGH", "CANONICAL", "rare", "ok")}
    got = {p_: (t.get("PGP_IMPACT"), t.get("PGP_TX"), t.get("PGP_FREQ"), t.get("PGP_CALL")) for p_, t in tags.items()}
    check("csq_rank.awk: impact on the MANE transcript, else the canonical one; frequency tier; call check",
          rc == 0 and got == want, (rc, got, out[-600:]))
    check("csq_rank.awk: what made a call weak (GQ12; GQ11,AB0.29)",
          (tags.get(200, {}).get("PGP_CALL_WHY"), tags.get(300, {}).get("PGP_CALL_WHY")) == ("GQ12", "GQ11,AB0.29"), tags)
    check("csq_rank.awk: the header names the floor it used",
          "##PGP_CALL_FLOOR=GQ=30,DP=10,HET_AB=0.2-0.8,HOM_AB=0.8\n" in out, out[:900])
    rc, tags_nc, _ = ranked("noncommon")
    check("csq_rank.awk: keep=noncommon drops the common record only",
          rc == 0 and sorted(tags_nc) == [100, 200, 300, 400, 500, 700], sorted(tags_nc))
    agree = collect_summary.csq_high_impact_genes(f"{d}/clinical/S_clinical.vcf.gz")
    check("csq_rank.awk and the summary agree on which records are HIGH and which calls are weak",
          {p_ for p_, t in tags.items() if t.get("PGP_IMPACT") == "HIGH"} == {k[1] for k in agree}
          and {p_ for p_, t in tags.items() if t.get("PGP_CALL") == "weak"} == {k[1] for k, x in agree.items() if x["weak"]},
          (tags, agree))
    summ, txt, html = render(d)
    sf = summ["sections"]["clinical"]["data"]["acmg_sf"]
    by = {x["gene"]: x for x in sf["high_impact"]}
    check("clinical MANE: a HIGH consequence off the MANE Select transcript is not listed (TTN)",
          "TTN" not in by, sorted(by))
    check("clinical MANE: without a MANE transcript the canonical HIGH one is listed (SCN5A)",
          by.get("SCN5A", {}).get("transcript") == "CANONICAL" and by["SCN5A"]["consequence"] == "stop_gained",
          by.get("SCN5A"))
    check("clinical MANE: the common stop-gain is not listed as rare HIGH (KCNQ1)", "KCNQ1" not in by, sorted(by))
    check("clinical call: the GQ 12 heterozygote with AD 3,2 is a weak call (BRCA2)",
          by.get("BRCA2", {}).get("weak_call") is True and "GQ 12 < 30" in (by["BRCA2"].get("call_note") or ""),
          by.get("BRCA2"))
    check("clinical call: the GQ 11 1/1 on AD 12,5 is weak by GQ and by its alt reads (LDLR)",
          by.get("LDLR", {}).get("weak_call") is True and "GQ 11" in by["LDLR"]["call_note"]
          and "0.29 of the depth on a homozygous call" in by["LDLR"]["call_note"], by.get("LDLR"))
    check("clinical call: the GQ 50 1/1 on AD 0,30 is not weak (MYBPC3)",
          by.get("MYBPC3", {}).get("weak_call") is False, by.get("MYBPC3"))
    check("clinical frequency: the deletion with no MAX_AF is listed as frequency unknown (MYH7)",
          by.get("MYH7", {}).get("frequency_unknown") is True and by["MYH7"]["frequency"] == "unknown", by.get("MYH7"))
    check("clinical: only the sound rare calls count as rare HIGH (MYBPC3, SCN5A); 2 weak, 1 unknown frequency",
          (sf.get("high_impact_rare"), sf.get("high_impact_weak"), sf.get("high_impact_frequency_unknown")) == (2, 2, 1), sf)
    check("clinical: the text report counts 2 rare and lists the others apart, with their marks",
          "0 ClinVar P/LP, 2 rare HIGH impact (listed apart: 2 weak call(s), 1 of unknown frequency)" in txt
          and "weak call (GQ 12 < 30): check the reads" in txt and "frequency unknown, not counted as rare" in txt, txt)
    check("clinical: the HTML report marks the weak call and the unknown frequency",
          "weak call (GQ 12 &lt; 30): check the reads" in html and "frequency unknown, not counted as rare" in html)
    # A floor of 0 written by bin/csq_rank.awk marks nothing.
    put(f"{d}/clinical/S_clinical.vcf.gz", head.replace(
        "#CHROM", "##PGP_CALL_FLOOR=GQ=0,DP=0,HET_AB=0-1,HOM_AB=0\n#CHROM") + body, gz=True)
    sf0 = collect_summary.collect("S", d)["sections"]["clinical"]["data"]["acmg_sf"]
    check("clinical call: with ##PGP_CALL_FLOOR GQ=0 and no allele-balance floor no call is weak",
          sf0.get("high_impact_weak") == 0 and sf0.get("high_impact_rare") == 4, sf0)
    # The summary tables of steps 23 and 31.
    cols = ("CHROM\tPOS\tREF\tALT\tGT\tIMPACT\tGENE\tConsequence\tMAX_AF\tCADD_PHRED\tREVEL\tAM_CLASS\tTRANSCRIPT\t"
            "WORST_IMPACT\tWORST_Consequence\tFREQ_TIER\tGQ\tDP\tAD\tCALL_CHECK\tCALL_WHY\n")
    put(f"{d}/clinical/S_clinical_summary.tsv", cols
        + "chr13\t200\tA\tAT\t0/1\tHIGH\tBRCA2\tframeshift_variant\t0.0001\t.\t.\t.\tMANE\tHIGH\tframeshift_variant\trare\t12\t.\t3,2\tweak\tGQ12\n"
        + "chr14\t500\tCTTT\tC\t0/1\tHIGH\tMYH7\tframeshift_variant\t.\t.\t.\t.\tMANE\tHIGH\tframeshift_variant\tunknown\t50\t30\t15,15\tok\t.\n"
        + "chr2\t100\tC\tT\t0/1\tMODIFIER\tTTN\tintron_variant\t0.0001\t.\t.\t.\tMANE\tHIGH\tstop_gained\trare\t50\t30\t15,15\tok\t.\n")
    put(f"{d}/slivar/S_slivar_summary.tsv",
        "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\tConsequence\tExisting_variation\tGT\tMAX_AF\tgnomADe_AF\tgnomADg_AF\t"
        "FREQ_TIER\tTRANSCRIPT\tMANE_Consequence\n"
        "chr14\t500\tCTTT\tC\tHIGH\tMYH7\tframeshift_variant\t.\t0/1\t.\t.\t.\tunknown\tMANE\tframeshift_variant\n"
        "chr13\t200\tA\tAT\tHIGH\tBRCA2\tframeshift_variant\t.\t0/1\t0.0001\t0.0001\t.\trare\tMANE\tframeshift_variant\n")
    summ, txt, html = render(d)
    cl, sl = summ["sections"]["clinical"]["data"], summ["sections"]["slivar"]["data"]
    check("clinical table: one row of unknown frequency and one weak call counted; slivar: one of unknown frequency",
          (cl.get("frequency_unknown"), cl.get("weak_calls"), sl.get("frequency_unknown")) == (1, 1, 1), (cl, sl))
    check("clinical table: the text report prints both counts and the slivar one",
          "Frequency unknown (not counted as rare): 1" in txt and "Weak genotype calls (check the reads): 1" in txt
          and "Frequency unknown (prioritized, not counted as rare): 1" in txt, txt)
    check("clinical table: the HTML cards print the counts",
          "Frequency unknown (not rare)" in html and "Weak genotype calls" in html)
    check("clinical MANE section validates", not validate.validate(schema, json.loads(json.dumps(summ))),
          validate.validate(schema, json.loads(json.dumps(summ))))


def pypgx_section(work, schema):
    """23. step 32's summary beyond the diplotypes, and the probable CYP2D6 call."""
    import html as htmlmod
    root = os.path.join(work, "pgxcn")
    d = os.path.join(root, "S")
    put(f"{d}/pypgx/S_pypgx_summary.tsv", PYPGX_HEAD + PYPGX_ROWS)
    put(f"{d}/pypgx/S_cyp2d6_depth_check.tsv", DEPTH_OK)
    put(f"{root}/reference/gencode.v50.basic.genes.gtf", GENES_GTF)
    put(f"{d}/cnvpytor/S_cnvs.txt", "deletion\tchr4:68500001-68650000\t150000\t0.0300\t1.2e-30\t1e-20\t1\t1\t0.01\t0\t9\n"
        "deletion\tchr4:68660001-68700000\t40000\t0.5\t1e-5\t1\t1\t1\t0\t0\t9\n")
    summ, txt, html_ = render(d)
    html_ = htmlmod.unescape(html_)
    pd = summ["sections"]["pypgx"]["data"]
    check("pypgx: the CNV call of the BAM genes, GSTT1's N/A left out",
          pd["cnv_call"] == {"CYP2D6": "Normal", "CYP2A6": "Normal", "GSTM1": "WholeDel2"}, pd.get("cnv_call"))
    check("pypgx: GSTT1 is not callable on this reference", pd["not_callable"] == ["GSTT1"], pd.get("not_callable"))
    check("pypgx: GSTT1 is not counted as called", pd["genes_called"] == 4, pd["genes_called"])
    check("pypgx: the GQ 12 call is flagged with its weakest variant",
          pd["low_gq"] == [{"gene": "CYP2C19", "diplotype": "*1/*2", "min_gq": 12,
                            "weakest": "chr10:94781859 G>A GQ 12 AD 26,3"}], pd.get("low_gq"))
    ugt = summ["sections"]["cnvpytor"]["data"]["ugt2b17"]
    check("UGT2B17: the deletion over its GENCODE span leaves no copy",
          ugt["call"] == "no copy (deleted on both chromosomes)" and ugt["span"] == "chr4:68516301-68627039"
          and ugt["gencode_release"] == "50", ugt)
    cn_lines = ["GSTM1    WholeDel2 (pypgx, read depth)",
                "UGT2B17  no copy (deleted on both chromosomes) (CNVpytor deletion chr4:68500001-68650000, "
                "depth level 0.03, e-value 1.2e-30, covers 100% of the gene)"]
    check("text report: the PGx copy number block", "## PGx copy number" in txt and all(x in txt for x in cn_lines), txt)
    check("HTML report: the PGx Copy Number card", "PGx Copy Number" in html_ and "WholeDel2 (pypgx, read depth)" in html_
          and "no copy (deleted on both chromosomes) (CNVpytor deletion chr4:68500001-68650000" in html_)
    lines = ["Not assessable on this reference: GSTT1 (GRCh38 has it only on the ALT contig chr22_KI270879v1_alt",
             "Low-confidence: CYP2C19 *1/*2 rests on chr10:94781859 G>A GQ 12 AD 26,3: check the reads",
             "CYP2D6 depth check: ok; share of the depth at MAPQ >= 1 / >= 20: Cyrius's region with REP6 0.91 / 0.84, "
             "gene body without REP6 0.99 / 0.97",
             "Probable (one caller: pypgx): *1/*2, Normal Metabolizer; not passed to PharmCAT"]
    for x in lines:
        h = x.replace("Probable (one caller: pypgx): ", "")
        check(f"both reports: {x[:60]}", x in txt and h in html_, txt if x not in txt else "HTML lacks it")
    check("HTML report: the probable call's label", "Probable (one caller: pypgx)" in html_)
    check("CYP2D6: the probable call in the summary",
          summ["cyp2d6"]["probable"] == {"diplotype": "*1/*2", "phenotype": "Normal Metabolizer",
                                         "label": "probable (one caller: pypgx)"}, summ["cyp2d6"])
    check("pypgx section summary validates", not validate.validate(schema, json.loads(json.dumps(summ))),
          validate.validate(schema, json.loads(json.dumps(summ))))

    # red controls for the probable call
    put(f"{d}/pypgx/S_cyp2d6_depth_check.tsv", DEPTH_OK.replace("status\tok", "status\tunreliable"))
    summ, txt, _ = render(d)
    check("CYP2D6: a failed depth check gives no probable call", summ["cyp2d6"]["probable"] is None
          and "Probable (one caller" not in txt, summ["cyp2d6"])
    put(f"{d}/pypgx/S_cyp2d6_depth_check.tsv", DEPTH_OK)
    cy_head = "Sample\tGenotype\tFilter\n"
    put(f"{d}/cyrius/S_cyp2d6.tsv", cy_head + "S\t*1/*2\tPASS\n")
    check("CYP2D6: beside a usable Cyrius call there is no probable call", render(d)[0]["cyp2d6"]["probable"] is None)
    put(f"{d}/cyrius/S_cyp2d6.tsv", cy_head + "S\tNone\tNot_assigned_to_haplotypes\n")
    check("CYP2D6: a Cyrius no-call leaves pypgx the only caller: probable",
          (render(d)[0]["cyp2d6"]["probable"] or {}).get("diplotype") == "*1/*2")
    os.remove(f"{d}/cyrius/S_cyp2d6.tsv")
    put(f"{d}/pypgx/S_pypgx_summary.tsv", "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\n"
        "CYP2D6\t*1/*2\tNormal Metabolizer\tNormal\tBAM\n")
    summ, txt, _ = render(d)
    check("CYP2D6: a summary without the evidence columns gives no probable call, no flags",
          summ["cyp2d6"]["probable"] is None and summ["sections"]["pypgx"]["data"]["low_gq"] == []
          and summ["sections"]["pypgx"]["data"]["not_callable"] == [], summ["cyp2d6"])

    # UGT2B17: one copy, no call, no gene lines, no CNVpytor
    put(f"{d}/cnvpytor/S_cnvs.txt", "deletion\tchr4:68530001-68600000\t70000\t0.4800\t3e-9\t1\t1\t1\t0\t0\t9\n")
    u = render(d)[0]["sections"]["cnvpytor"]["data"]["ugt2b17"]
    check("UGT2B17: a level 0.48 deletion leaves one copy", u["call"] == "one copy", u)
    put(f"{d}/cnvpytor/S_cnvs.txt", "deletion\tchr4:1000001-1100000\t100000\t0.0\t1e-9\t1\t1\t1\t0\t0\t9\n")
    u = render(d)[0]["sections"]["cnvpytor"]["data"]["ugt2b17"]
    check("UGT2B17: no CNVpytor call over the gene reads two copies, and says why",
          u["call"] == "two copies" and u["note"] == "no CNVpytor call over the gene", u)
    os.remove(f"{root}/reference/gencode.v50.basic.genes.gtf")
    summ, txt, _ = render(d)
    check("UGT2B17: without the GENCODE gene lines it is not assessed, and the report says why",
          summ["sections"]["cnvpytor"]["data"]["ugt2b17"]["call"] == "not assessed"
          and "UGT2B17  not assessed (no GENCODE gene lines" in txt, txt)
    os.remove(f"{d}/cnvpytor/S_cnvs.txt")
    txt = render(d)[1]
    check("UGT2B17: without CNVpytor the block says it did not run",
          "UGT2B17  not assessed (CNVpytor, step 18, did not run)" in txt, txt)


def main():
    work = tempfile.mkdtemp()
    schema = json.load(open(os.path.join(REPO, "tests", "schema", "summary.schema.json")))
    try:
        # 1. bash folder: schema, and both reports agree
        d = os.path.join(work, "bash", "S")
        bash_folder(d)
        summ, txt, html = render(d)
        errs = validate.validate(schema, json.loads(json.dumps(summ)))
        check("bash folder: the summary validates against the schema", not errs, errs)
        sec = summ["sections"]
        check("CPIC: an unclassified gene (VKORC1 -1639 GG) is not counted as non-normal",
              sec["cpic"]["data"]["non_normal"] == 1, sec["cpic"]["data"])
        check("ClinVar: 3 hits, best-reviewed first", sec["clinvar"]["data"]["count"] == 3
              and [h["gene"] for h in sec["clinvar"]["data"]["hits"]] == ["GENEB", "GENEA", "GENEC"],
              sec["clinvar"]["data"])
        check("ClinVar: count by stars", sec["clinvar"]["data"]["by_stars"] == {"4": 0, "3": 1, "2": 0, "1": 1, "0": 1},
              sec["clinvar"]["data"]["by_stars"])
        # Step 06's cDNA read check, joined to each hit by chrom, pos, ref and alt;
        # a hit the table lacks is not checked, never clear.
        cd = os.path.join(work, "cdna", "S")
        put(f"{cd}/clinvar/S_clinvar_hits.vcf", HITS)
        put(f"{cd}/clinvar/S_clinvar_cdna_check.tsv",
            "chrom\tpos\tref\talt\talt_reads\talt_clipped\tref_reads\tref_clipped\tcdna_check\tcdna_note\n"
            "chr1\t10\tC\tT\t9\t8\t11\t0\tflagged\tpossible cDNA/transgene contamination: check the reads\n"
            "chr1\t20\tT\tA\t9\t0\t0\t0\tclear\t\n")
        cv = collect_summary.collect("S", cd)["sections"]["clinvar"]["data"]
        got = {h["gene"]: (h["cdna"]["check"], h["cdna"]["note"]) for h in cv["hits"]}
        check("cDNA check: one flagged hit, the others not checked (GENEB's alt differs from the table's row)",
              cv["cdna_flagged"] == 1 and got["GENEA"][0] == "flagged" and got["GENEA"][1].startswith("possible cDNA")
              and got["GENEB"] == ("not_checked", "not in the read-check table")
              and got["GENEC"] == ("not_checked", "not in the read-check table"), got)
        check("text report: the ClinVar count line", "  Pathogenic/Likely Pathogenic hits: 3\n" in txt)
        check("HTML report: the same ClinVar count in its badge, yellow (never green) with a hit",
              'class="badge badge-yellow">3</span>' in html and 'badge-green">3<' not in html)
        rows = [l for l in html.splitlines() if l.strip().startswith("<tr><td>chr")]
        check("HTML report: one row per hit, each with a Stars cell and an empty combined-allele cell",
              len(rows) == 3 and rows[0].strip().endswith("<td>3</td><td></td><td>gnomAD AF unknown</td><td></td></tr>"),
              rows)
        for r in sec["prs"]["data"]["scores"]:
            check(f"PRS row {r['pgs_id']} in both reports", r["pgs_id"] in txt and f"<td>{r['pgs_id']}</td>" in html)
        check("QC: mean depth 29.5x in both reports", "Mean depth: 29.5x" in txt and "29.5x" in html,
              sec["coverage"]["data"])
        check("QC: inferred sex female", sec["sex_check"]["data"]["inferred_sex"] == "female")
        check("variants: 3 records, 2 PASS, 2 SNP lines, 1 indel",
              sec["variants"]["data"] == {"total": 3, "pass": 2, "snps": 2, "indels": 1}, sec["variants"]["data"])
        check("CYP2D6: pypgx and Cyrius agree (allele order ignored)", summ["cyp2d6"]["agree"] is True, summ["cyp2d6"])
        # 15. CPIC unclassified
        check("CPIC: the unclassified gene is counted apart (VKORC1 -1639 GG)",
              sec["cpic"]["data"].get("unclassified") == 1, sec["cpic"]["data"])
        check("CPIC: both reports show the unclassified count",
              "Called without a function phenotype (no drug guidance): 1" in txt
              and 'Genes called without a function phenotype (no drug guidance)</span><span class="value">1<' in html,
              [l for l in txt.splitlines() if "Not called" in l])

        # 3. heteroplasmy floor
        check("heteroplasmy: only AF 0.05 to 0.95 counts (0.30, 0.05, 0.949 of 5 PASS)",
              sec["mito"]["data"]["heteroplasmic"] == 3 and sec["mito"]["data"]["pass"] == 5, sec["mito"]["data"])
        # 4. CPSR by column name
        check("CPSR: breakdown read from the CLASSIFICATION column by name",
              sec["cpsr"]["data"].get("classification") == {"VUS": 2, "Likely_Pathogenic": 1}, sec["cpsr"]["data"])

        # 2. Nextflow folder
        n = os.path.join(work, "nf", "S")
        put(f"{n}/roh/S_roh.txt", "RG\tS\tchr2\t1\t6000001\t6000000\t3\t50\n")
        put(f"{n}/roh/S_roh_summary.txt", "summary\n")
        put(f"{n}/coverage/S.mosdepth.summary.txt", "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t1\t1\t31.25\t0\t9\n")
        put(f"{n}/hla/S_hla_genotype.tsv", "HLA-B\t2\tB*07:02:01\t30\t60\tB*08:01:01\t20\t55\t\n")
        put(f"{n}/pharmcat/S.report.json", json.dumps({"pharmcatVersion": "3.2.0", "genes": {}}))
        summ_n, txt_n, _ = render(n)
        s2 = summ_n["sections"]
        check("Nextflow folder: ROH read from roh/", s2["roh"]["state"] == "ok"
              and s2["roh"]["data"]["autosomal_over_5mb"] == [{"region": "chr2:1-6000001", "mb": 6.0}], s2["roh"])
        check("Nextflow folder: ROH shown in the text report", "Autosomal ROH > 5MB: 1" in txt_n)
        check("Nextflow folder: depth from coverage/", s2["coverage"]["data"]["mean_depth"] == 31.25)
        check("Nextflow folder: HLA from hla/", s2["hla"]["data"]["loci"][0]["alleles"] == ["B*07:02:01", "B*08:01:01"])
        check("Nextflow folder: PharmCAT from pharmcat/", s2["pharmcat"]["data"].get("version") == "3.2.0")
        check("Nextflow folder: the summary validates", not validate.validate(schema, json.loads(json.dumps(summ_n))))

        # 5. stale results
        st = os.path.join(work, "stale", "S")
        bash_folder(st)
        put(f"{st}/slivar/S_slivar_summary.tsv", "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\nchr1\t1\tA\tG\tHIGH\tX\n", age_days=30)
        put(f"{st}/vcf/S_roh.txt", "RG\tS\tchr1\t1\t6000001\t6000000\t3\t50\n", age_days=30)
        put(f"{st}/logs/run_status.tsv", f"meta\tstarted_epoch\t{time.time() - DAY}\nmeta\tdeclared_sex\tmale\n"
            "step\t31\tskipped (needs VEP, step 13 skipped)\nstep\t11\tok\n")
        summ_s, txt_s, html_s = render(st)
        ss = summ_s["sections"]
        check("stale: step 31 skipped, slivar summary from an earlier run -> stale",
              ss["slivar"]["state"] == "stale" and "skipped" in (ss["slivar"]["note"] or ""), ss["slivar"])
        check("stale: marked in the text report", "[STALE: from an earlier run" in txt_s)
        check("stale: marked in the HTML report", '<div class="stale">STALE: from an earlier run' in html_s)
        check("stale: step 11 ok in that run, an old ROH file is current", ss["roh"]["state"] == "ok", ss["roh"])
        check("stale: a file written after the run started is current", ss["clinvar"]["state"] == "ok")
        check("stale: declared sex from the run status, mismatch flagged",
              summ_s["run"]["declared_sex"] == "male" and "DOES NOT MATCH" in txt_s)
        check("no run status: nothing is stale",
              all(x["state"] != "stale" for x in summ["sections"].values()))

        # 6. stars: awk and python agree
        statuses = list(collect_summary.STARS) + ["no_assertion_criteria_provided", "no_classification_provided",
                                                  "something_new", ""]
        vcf = os.path.join(work, "stars.vcf")
        put(vcf, VCF_HEAD + "".join(f"chr1\t{i}\t.\tA\tG\t50\tPASS\tCLNREVSTAT={st_}\tGT\t0/1\n"
                                    for i, st_ in enumerate(statuses, 1)))
        out = subprocess.run(["awk", "-f", os.path.join(REPO, "bin", "clinvar_hits.awk"), vcf],
                             capture_output=True, text=True).stdout.split("\n")
        awk_stars = [int(l.split("\t")[0]) for l in out if l]
        py_stars = [collect_summary.STARS.get(x, 0) for x in statuses]
        check("bin/clinvar_hits.awk and collect_summary.py give the same stars", awk_stars == py_stars,
              f"awk {awk_stars} python {py_stars}")

        # 6b. the combined-allele note (step 6's PGP_* tags): awk and python
        # word it the same, and both reports show it
        cv = os.path.join(work, "combined", "S")
        put(f"{cv}/clinvar/S_clinvar_hits.vcf", COMBINED_HITS)
        out = subprocess.run(["awk", "-f", os.path.join(REPO, "bin", "clinvar_hits.awk"),
                              f"{cv}/clinvar/S_clinvar_hits.vcf"], capture_output=True, text=True).stdout
        awk_notes = {(l.split("\t")[2]): l.split("\t")[9] for l in out.split("\n") if l}
        summ_cv, txt_cv, html_cv = render(cv)
        py_notes = {str(h["pos"]): h["combined_allele"] for h in summ_cv["sections"]["clinvar"]["data"]["hits"]}
        # All three are unphased het calls: the genotype cannot say the parts
        # share a haplotype, so the note says so (a phased pair below does not).
        unph = "; unphased, so the parts may sit on different haplotypes"
        want = {"30000037": "part of a combined allele: chr22:30000038 GCCAT>CATGGG" + unph,
                "30000038": "combined allele of chr22:30000037 AGC>A + chr22:30000042 T>TGGG" + unph,
                "30000042": "part of a combined allele: chr22:30000038 GCCAT>CATGGG" + unph,
                "30000210": ".",
                "30000300": "possibly one allele with chr22:30000305 C>T: check ClinVar for a combined record",
                "30000305": "possibly one allele with chr22:30000300 A>G: check ClinVar for a combined record"}
        check("combined allele: the note of each hit in the summary", py_notes == want, py_notes)
        check("combined allele: bin/clinvar_hits.awk gives the same notes", awk_notes == want, awk_notes)
        check("combined allele: the 2-star combined record is listed before its 1-star parts",
              [h["pos"] for h in summ_cv["sections"]["clinvar"]["data"]["hits"]][0] == 30000038,
              summ_cv["sections"]["clinvar"]["data"]["hits"])
        check("combined allele: the text report shows the notes",
              f"(part of a combined allele: chr22:30000038 GCCAT>CATGGG{unph})" in txt_cv
              and f"(combined allele of chr22:30000037 AGC>A + chr22:30000042 T>TGGG{unph})" in txt_cv, txt_cv)
        cv_rows = [l for l in html_cv.splitlines() if l.strip().startswith("<tr><td>chr22")]
        check("combined allele: the HTML table has the column and the notes, an SNV hit's cell empty",
              "<th>Stars</th><th>Combined allele</th>" in html_cv
              and any(f"<td>part of a combined allele: chr22:30000038 GCCAT&gt;CATGGG{unph}</td>" in r for r in cv_rows)
              and any("<td>30000210</td>" in r and r.strip().endswith("<td></td></tr>") for r in cv_rows), cv_rows)
        # The same hits phased (0|1, one phase set) or homozygous: both
        # parts are on one haplotype, so the notes carry no phase caveat.
        for gt_ in ("0|1", "1/1"):
            put(f"{cv}/clinvar/S_clinvar_hits.vcf", COMBINED_HITS.replace("\tGT\t0/1\n", f"\tGT\t{gt_}\n"))
            out = subprocess.run(["awk", "-f", os.path.join(REPO, "bin", "clinvar_hits.awk"),
                                  f"{cv}/clinvar/S_clinvar_hits.vcf"], capture_output=True, text=True).stdout
            awk_g = {(l.split("\t")[2]): l.split("\t")[9] for l in out.split("\n") if l}
            summ_g, _, _ = render(cv)
            py_g = {str(h["pos"]): h["combined_allele"] for h in summ_g["sections"]["clinvar"]["data"]["hits"]}
            want_g = {k: v.replace(unph, "") for k, v in want.items()}
            check(f"combined allele: GT {gt_} gets no phase caveat, in the summary and in bin/clinvar_hits.awk",
                  py_g == want_g and awk_g == want_g, (py_g, awk_g))
        put(f"{cv}/clinvar/S_clinvar_hits.vcf", COMBINED_HITS)
        check("combined allele: the summary validates", not validate.validate(schema, json.loads(json.dumps(summ_cv))),
              validate.validate(schema, json.loads(json.dumps(summ_cv))))

        # 7. a corrupt file
        c = os.path.join(work, "corrupt", "S")
        put(f"{c}/vcf/S.vcf.gz", "not gzip\n")
        put(f"{c}/clinvar/S_clinvar_hits.vcf", HITS)
        summ_c, txt_c, _ = render(c)
        check("corrupt VCF: section unreadable, the report still renders",
              summ_c["sections"]["variants"]["state"] == "unreadable" and "hits: 3" in txt_c,
              summ_c["sections"]["variants"])
        # a gzip whose header is fine but whose compressed data is damaged
        # (a truncated copy, a bad disk) raises zlib.error, not OSError
        z = os.path.join(work, "corrupt-deflate", "S")
        os.makedirs(f"{z}/vcf")
        blob = bytearray(gzip.compress((VCF_HEAD + "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n" * 2000).encode()))
        blob[20:40] = b"\xff" * 20
        with open(f"{z}/vcf/S.vcf.gz", "wb") as f:
            f.write(bytes(blob))
        put(f"{z}/clinvar/S_clinvar_hits.vcf", HITS)
        try:
            summ_z, txt_z, _ = render(z)
            state_z = summ_z["sections"]["variants"]
        except Exception as e:  # the pre-fix behaviour: the whole report dies
            state_z, txt_z = f"raised {type(e).__name__}: {e}", ""
        check("damaged gzip data: section unreadable, the report still renders",
              isinstance(state_z, dict) and state_z["state"] == "unreadable" and "hits: 3" in txt_z, state_z)

        # 8. ACMG SF tier: BRCA2 (on the list) and GENEA (not) in ClinVar; TTN
        # (on the list) and GENEX (not) with HIGH impact; MYH7 on the list
        # but MODERATE only
        a = os.path.join(work, "acmg", "S")
        put(f"{a}/clinvar/S_clinvar_hits.vcf", VCF_HEAD + (
            "chr13\t32340300\t1\tG\tA\t50\tPASS\tGENEINFO=BRCA2:675;CLNSIG=Pathogenic;"
            "CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n"
            "chr1\t10\t2\tC\tT\t50\tPASS\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;"
            "CLNREVSTAT=criteria_provided,_single_submitter\tGT\t0/1\n"))
        csq = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
               'Format: Allele|Consequence|IMPACT|SYMBOL">\n')
        put(f"{a}/clinical/S_clinical.vcf.gz", "##fileformat=VCFv4.2\n" + csq + VCF_HEAD.split("\n", 1)[1] + (
            "chr2\t178500000\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|TTN,T|intron_variant|MODIFIER|TTN-AS1\tGT\t0/1\n"
            "chr3\t100\t.\tA\tG\t50\tPASS\tCSQ=G|frameshift_variant|HIGH|GENEX\tGT\t1/1\n"
            "chr14\t23400000\t.\tG\tA\t50\tPASS\tCSQ=A|missense_variant|MODERATE|MYH7\tGT\t0/1\n"), gz=True)
        summ_a, txt_a, html_a = render(a)
        sf = summ_a["sections"]["clinical"]["data"].get("acmg_sf") or {}
        check("ACMG SF: v3.3, 84 genes", sf.get("version") == "ACMG SF v3.3" and sf.get("genes_on_list") == 84, sf)
        check("ACMG SF: the ClinVar hit in BRCA2 only", [x["gene"] for x in sf.get("clinvar_hits", [])] == ["BRCA2"], sf)
        check("ACMG SF: the HIGH-impact TTN record only (not GENEX, not MODERATE MYH7)",
              [(x["gene"], x["consequence"]) for x in sf.get("high_impact", [])] == [("TTN", "stop_gained")], sf)
        check("ACMG SF: the text report lists BRCA2 and TTN",
              "1 ClinVar P/LP, 1 rare HIGH impact" in txt_a and "    BRCA2" in txt_a and "    TTN" in txt_a, txt_a)
        check("ACMG SF: the HTML report has the card", "Secondary-Findings Genes (ACMG SF v3.3)" in html_a
              and "<td>TTN</td>" in html_a and "<td>GENEX</td>" not in html_a)

        # 9. haplocheck and Yleaf
        h = os.path.join(work, "mito", "S")
        put(f"{h}/mito/S_haplogroup.txt", '"SampleID"\t"Haplogroup"\t"Rank"\n"S"\t"H1a"\t"0.95"\n')
        put(f"{h}/mito/S_haplocheck.txt", '"Sample"\t"Contamination Status"\t"Contamination Level"\t"Distance"\n'
            '"S"\t"YES"\t"0.12"\t"5"\n')
        put(f"{h}/y_haplogroup/S_y_haplogroup.txt", "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\t"
            "QC-score\tQC-1\tQC-2\tQC-3\nS_sorted\tR-M269\tM269\t100\t42\t0.97\t1\t1\t0.97\n")
        summ_h, txt_h, html_h = render(h)
        hd = summ_h["sections"]["haplogroup"]["data"]
        check("haplocheck: status and level read", hd.get("contamination_status") == "YES"
              and hd.get("contamination_level") == "0.12", hd)
        check("haplocheck: the line in both reports", "Contamination (haplocheck): YES, two mtDNA haplogroups" in txt_h
              and "YES, two mtDNA haplogroups" in html_h, txt_h)
        check("Yleaf: the Y haplogroup in both reports", "R-M269 (42 markers, QC-score 0.97)" in txt_h
              and "Y-Chromosome Haplogroup" in html_h and "R-M269" in html_h, txt_h)
        put(f"{h}/y_haplogroup/S_y_haplogroup.txt", "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\t"
            "QC-score\tQC-1\tQC-2\tQC-3\nS_sorted\tNA\t\t100\t3\t0\t0\t0\t0\n")
        summ_n, txt_n, _ = render(h)
        check("Yleaf: Hg NA reads as insufficient markers",
              summ_n["sections"]["y_haplogroup"]["data"].get("haplogroup") == "insufficient markers"
              and "Haplogroup: insufficient markers (3 markers" in txt_n, summ_n["sections"]["y_haplogroup"])
        os.remove(f"{h}/mito/S_haplocheck.txt")
        check("no haplocheck file: the report says it was not checked",
              "Contamination (haplocheck): not checked" in render(h)[1])

        # 9b. step 41: the blood group table
        b = os.path.join(work, "blood", "S")
        put(f"{b}/bloodtype/S_bloodtype.tsv", "system\tphenotype\tgenotype\tconfidence\tnote\n"
            "ABO\tA1\tA1/O.01\thigh\t\nRhD\tD-\tRHD deleted/RHD deleted (RHD copies by read depth 0.02)\thigh\t\n"
            "Kell\tK-k+ Kp(a-b+) Js(a-b+)\tk/k Kpb/Kpb Jsb/Jsb\tmedium\tx\n")
        summ_b, txt_b, html_b = render(b)
        bd = summ_b["sections"]["bloodtype"]
        check("blood groups: read, with ABO and RhD picked out", bd["state"] == "ok" and bd["data"].get("abo") == "A1"
              and bd["data"].get("rhd") == "D-" and len(bd["data"].get("systems", [])) == 3, bd)
        check("blood groups: step 41 is a run-all step (a result older than the run can be stale)",
              bd["step"] == "41" and "41" in collect_summary.RUN_ALL_STEPS, bd["step"])
        check("blood groups: the summary still passes the schema",
              not validate.validate(schema, json.loads(json.dumps(summ_b))),
              validate.validate(schema, json.loads(json.dumps(summ_b))))
        put(f"{b}/bloodtype/S_bloodtype.tsv", "system\tantigens\nABO\tA\n")
        bad_b = render(b)[0]["sections"]["bloodtype"]
        check("blood groups: a table without the columns is unreadable, not empty",
              bad_b["state"] == "unreadable" and "phenotype" in (bad_b["note"] or ""), bad_b)
        check("blood groups: a sample without step 41 lists it as not run",
              "Blood groups (predicted)" in render(d)[0]["not_run"], render(d)[0]["not_run"])

        # 10. step 33 with the declared sex as somalier's pedigree
        head = ("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\toriginal_pedigree_sex\t"
                "gt_depth_mean\tn_hom_ref\tn_het\tn_hom_alt\tp_middling_ab\tX_depth_mean\tX_n\tX_hom_ref\t"
                "X_het\tX_hom_alt\tY_depth_mean\tY_n\n")

        def qc(sex, ped, x_n, x_het, x_hom_alt, declared, y_depth="0.0", mid="0.010"):
            path = os.path.join(work, "somalier.samples.tsv")
            put(path, head + f"S\tS\t-9\t-9\t{sex}\t-9\t{ped}\t30.0\t60\t40\t20\t{mid}\t15.0\t{x_n}\t0\t"
                f"{x_het}\t{x_hom_alt}\t{y_depth}\t5\n")
            return dict(collect_summary.sample_qc_table("S", path, declared_sex=declared))

        t = qc("1", "male", 2, 0, 2, "male")
        check("pedigree male kept on 2 chrX sites: unknown, not checked (not ok)",
              t["inferred_sex"] == "unknown" and t["sex_check"] == "not_checked", t)
        t = qc("1", "female", 39, 0, 39, "female")
        check("pedigree female, somalier set male from 39 chrX sites: mismatch",
              t["inferred_sex"] == "male" and t["sex_check"] == "mismatch", t)
        t = qc("1", "male", 39, 0, 39, "male")
        check("pedigree male, 39 homozygous chrX sites agree: ok", t["sex_check"] == "ok", t)
        t = qc("1", "male", 39, 0, 39, "male", mid="0.080")
        check("pedigree male kept on a sample somalier calls low quality (8% middling allele balance): not checked",
              t["sex_check"] == "not_checked", t)
        t = qc("1", "male", 39, 0, 39, "male", mid="nan")
        check("pedigree male kept, middling allele balance nan (no autosomal call): not checked",
              t["sex_check"] == "not_checked", t)
        t = qc("-2", "female", 4, 0, 4, "female", y_depth="14.0")
        check("pedigree female with chrY reads, chrX cannot tell: unknown, not the aneuploidy line",
              t["sex_check"] == "not_checked" and "could not tell" in t["sex_check_reason"], t)
        t = qc("-2", "female", 40, 20, 20, "female", y_depth="14.0")
        check("pedigree female, chrX heterozygous and chrY reads: the aneuploidy line",
              t["sex_check"] == "not_checked" and "chrY has reads" in t["sex_check_reason"], t)
        t = qc("1", "-9", 2, 0, 2, "male")
        check("no pedigree (-9): somalier's sex column is read as before",
              t["inferred_sex"] == "male" and t["sex_check"] == "ok", t)

        # 11. CYP2D6 with a failed depth check: pypgx Indeterminate, Cyrius's genotype filtered
        y = os.path.join(work, "cyp", "S")
        put(f"{y}/pypgx/S_pypgx_summary.tsv", "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\n"
            "CYP2D6\tIndeterminate\tIndeterminate\t.\tbam\nCYP2C19\t*1/*2\tIM\t.\tbam\nCYP2C9\tindeterminate\t.\t.\tbam\n",
            age_days=2)
        cy_head = "Sample\tGenotype\tFilter\n"
        put(f"{y}/cyrius/S_cyp2d6.tsv", cy_head + "S\t*5/*5\tCYP2D6_depth_unreliable\n", age_days=2)
        summ_y, txt_y, html_y = render(y)
        cy = summ_y["cyp2d6"]
        check("pypgx: Indeterminate (any case) is not counted as called",
              summ_y["sections"]["pypgx"]["data"]["genes_called"] == 1
              and "Genes called: 1/3" in txt_y, summ_y["sections"]["pypgx"]["data"])
        check("CYP2D6: a Cyrius genotype that failed its Filter is shown as not usable",
              cy["calls"]["Cyrius"] == "*5/*5 (not usable: Filter CYP2D6_depth_unreliable)", cy["calls"])
        check("CYP2D6: neither call is usable, so no agreement is claimed (not 'disagree')", cy["agree"] is None, cy)
        check("CYP2D6: the text report says so", "not usable: Filter CYP2D6_depth_unreliable" in txt_y
              and "pypgx and Cyrius agree: fewer than two usable calls" in txt_y, txt_y)
        check("CYP2D6: the HTML card says so", "not usable: Filter CYP2D6_depth_unreliable" in html_y
              and "fewer than two usable calls" in html_y)
        check("CYP2D6: no step 36 table, no verdict", cy.get("consensus") is None and "Step 36" not in txt_y, cy)
        put(f"{y}/pgx_consensus/S_pgx_consensus.tsv", "Gene\tResult\tOutside_call\tReason\tEvidence\n"
            "HLA-A\tnot typed\tno\tHLA typing (step 08) did not run\t-\n"
            "CYP2D6\tindeterminate\tno\tthe CYP2D6 depth check failed\tpypgx: no call (Indeterminate)\n")
        summ_y, txt_y, html_y = render(y)
        v36 = summ_y["cyp2d6"].get("consensus") or {}
        check("CYP2D6: step 36's verdict from its table",
              v36.get("result") == "indeterminate" and v36.get("passed_to_pharmcat") is False, summ_y["cyp2d6"])
        check("CYP2D6: both reports print step 36's verdict and reason",
              "Step 36: indeterminate, not passed to PharmCAT (the CYP2D6 depth check failed)" in txt_y
              and "Not passed to PharmCAT: the CYP2D6 depth check failed" in html_y, txt_y)
        check("CYP2D6 summary validates", not validate.validate(schema, json.loads(json.dumps(summ_y))),
              validate.validate(schema, json.loads(json.dumps(summ_y))))
        put(f"{y}/pgx_consensus/S_pgx_consensus.tsv", open(f"{y}/pgx_consensus/S_pgx_consensus.tsv").read(), age_days=5)
        check("CYP2D6: a step 36 table older than the callers' results is not read",
              render(y)[0]["cyp2d6"].get("consensus") is None)

        # 12. repeat expansions with Stranger: invented loci, none of the five
        # key loci, so a flagged one can only come from STR_STATUS
        x = os.path.join(work, "str", "S")
        fmt = "GT:REPCN"
        eh_rec = ("chr4\t1000\t.\tC\t<STR20>\t.\tPASS\tEND=1;REPID=LOCUS_A;VARID=LOCUS_A{st}\t" + fmt + "\t0/1:7/20\n"
                  "chr14\t2000\t.\tG\t<STR75>\t.\tPASS\tEND=1;REPID=LOCUS_B;VARID=LOCUS_B{st2}\t" + fmt + "\t0/1:20/75\n"
                  "chr9\t100\t.\tA\t<STR9>\t.\tPASS\tEND=1;REPID=NOTINCAT;VARID=NOTINCAT\t" + fmt + "\t0/1:5/9\n")
        put(f"{x}/expansion_hunter/S_eh.vcf", VCF_HEAD + eh_rec.format(st="", st2=""), age_days=1)
        put(f"{x}/expansion_hunter/S_eh_stranger.vcf",
            VCF_HEAD + eh_rec.format(st=";STR_STATUS=normal", st2=";STR_STATUS=full_mutation"))
        summ_x, txt_x, html_x = render(x)
        ed = summ_x["sections"]["expansions"]["data"]
        check("Stranger: read when it is there", summ_x["sections"]["expansions"]["source"]
              == "expansion_hunter/S_eh_stranger.vcf" and ed.get("stranger") is True, summ_x["sections"]["expansions"])
        check("Stranger: the full_mutation outside the five key loci is listed, the normal one is not",
              ed.get("flagged") == [{"locus": "LOCUS_B", "repeat_count": "20/75", "status": "full_mutation"}], ed)
        check("Stranger: a record without STR_STATUS is counted apart, not flagged", ed.get("no_status") == 1, ed)
        check("Stranger: the five key loci stay listed", [k["locus"] for k in ed["key_loci"]] == collect_summary.EH_LOCI, ed)
        check("Stranger: both reports list LOCUS_B full_mutation and the short-read caveat",
              "LOCUS_B  20/75  full_mutation" in txt_x and "<td>LOCUS_B</td><td>20/75</td><td>full_mutation</td>" in html_x
              and "can be wrong at some loci" in txt_x and "can be wrong at some loci" in html_x, txt_x)
        check("Stranger summary validates", not validate.validate(schema, json.loads(json.dumps(summ_x))))
        put(f"{x}/expansion_hunter/S_eh_stranger.vcf", open(f"{x}/expansion_hunter/S_eh_stranger.vcf").read(), age_days=3)
        summ_x2, txt_x2, html_x2 = render(x)
        check("Stranger: a file older than ExpansionHunter's is not read",
              summ_x2["sections"]["expansions"]["source"] == "expansion_hunter/S_eh.vcf"
              and not summ_x2["sections"]["expansions"]["data"].get("stranger")
              and summ_x2["sections"]["expansions"]["data"].get("stranger_outdated") is True,
              summ_x2["sections"]["expansions"])
        check("Stranger: an outdated file makes both reports say to rerun step 9b, not that it did not run",
              "rerun step 9b" in txt_x2 and "did not run" not in txt_x2
              and "rerun step 9b" in html_x2 and "not run" not in html_x2.split("<h2>Repeat Expansions</h2>", 1)[-1]
              .split("<h2>", 1)[0], txt_x2)
        check("Stranger outdated summary validates", not validate.validate(schema, json.loads(json.dumps(summ_x2))))
        put(f"{x}/logs/run_status.tsv", f"meta\tstarted_epoch\t{time.time() - 0.5 * DAY}\nstep\t09\tskipped\n")
        _, _, html_x3 = render(x)
        card_eh = html_x3.split("<h2>Repeat Expansions</h2>", 1)[-1].split("<h2>", 1)[0]
        check("Repeat card of a stale result says Stale, not Complete",
              "Stale" in card_eh and "Complete" not in card_eh, card_eh)

        # 13. HLA quality per allele
        q = os.path.join(work, "hla", "S")
        put(f"{q}/hla_t1k/S_hla_genotype.tsv",
            "HLA-A\t1\tA*01:01:01\t30\t60\t.\t0\t-1\n"
            "HLA-B\t2\tB*57:01:01\t10\t0\tB*08:01:01\t20\t40\n"
            "HLA-C\t2\tC*07:01:01\t5\t-1\tC*07:02:01\t5\t30\n")
        summ_q, txt_q, html_q = render(q)
        loci = {l["gene"]: l for l in summ_q["sections"]["hla"]["data"]["loci"]}
        check("HLA: a one-allele row keeps the allele's quality only and is not low confidence",
              loci["HLA-A"]["alleles"] == ["A*01:01:01"] and loci["HLA-A"]["quality"] == ["60"]
              and loci["HLA-A"].get("low_confidence") is False, loci["HLA-A"])
        check("HLA: HLA-B with a quality 0 allele is low confidence and withheld from PharmCAT",
              loci["HLA-B"]["quality"] == ["0", "40"] and loci["HLA-B"].get("low_confidence") is True
              and loci["HLA-B"].get("withheld_from_pharmcat") is True, loci["HLA-B"])
        check("HLA: HLA-C with a quality -1 allele is low confidence, but PharmCAT never gets HLA-C",
              loci["HLA-C"].get("low_confidence") is True and loci["HLA-C"].get("withheld_from_pharmcat") is False,
              loci["HLA-C"])
        check("HLA: the text report shows each quality and says only HLA-B is not passed to PharmCAT",
              "A*01:01:01 (quality 60)\n" in txt_q and "B*57:01:01 (quality 0) / B*08:01:01 (quality 40); low confidence"
              in txt_q and txt_q.count("low confidence (T1K quality 0 or below), not passed to PharmCAT") == 1
              and "C*07:02:01 (quality 30); low confidence (T1K quality 0 or below)\n" in txt_q, txt_q)
        check("HLA: the HTML card shows the same, and the short-read note",
              "B*57:01:01 (quality 0)" in html_q and html_q.count("not passed to PharmCAT</span>") == 1
              and "HLA typing from short-read WGS is approximate" in html_q)
        check("HLA summary validates", not validate.validate(schema, json.loads(json.dumps(summ_q))))

        # 17. step 35b: the Paraphase genes, SMN1 beside Parascopy, both reports
        q = os.path.join(work, "paralogs", "S")
        put(f"{q}/paralogs/S_paraphase_copy_number.tsv",
            "region\tgene\tcopy_number\tstatus\tdetail\ttotal_cn\tregion_depth\tgenome_depth\n"
            "smn1\tSMN1\t2\tcalled\t.\t4\t120\t31\n"
            "smn1\tSMN2\t.\tno call\tone haplotype seen\t4\t120\t31\n"
            "hba\tHBA1\t2\tcalled\tgenotype -a/aa (3 alpha-globin genes); 4.2 kb deletion (-a4.2): HBA2 lost\t3\t90\t31\n"
            "hba\tHBA2\t1\tcalled\tgenotype -a/aa (3 alpha-globin genes); 4.2 kb deletion (-a4.2): HBA2 lost\t3\t90\t31\n"
            "TNXB\tTNXB\t4\tcalled\tcopies of the whole gene family\t4\t60\t31\n"
            "SSX2\tSSX2,SSX2B\t.\tno call\t.\t.\t0\t31\n")
        put(f"{q}/paralogs/S_paraphase_variants.tsv", "region\tchrom\tpos\tref\talt\tfilter\thaplotypes\tsupport\n"
            "hba\tchr16\t172900\tC\tT\tPASS\thba_hba2hap1\thba_hba2hap1:20/21\n")
        put(f"{q}/paralogs/S_smn_copy_number.tsv", "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\t"
            "psCN_filter\tpsCN\tpsCN_qual\thomologous_regions\nchr5\t70925000\t70953000\tSMN1\tPASS\t3\t40\t"
            "PASS\t1,2\t33,30\tchr5:70049000-70077000:-\n")
        summ_p, txt_p, html_p = render(q)
        pd = summ_p["sections"]["paralogs"]
        check("paralogs: read, HBA1 then HBA2 kept in order, two other regions counted (one called)",
              pd["state"] == "ok" and [g["gene"] for g in pd["data"]["genes"]] == ["SMN1", "SMN2", "HBA1", "HBA2"]
              and pd["data"]["other_regions"] == 2 and pd["data"]["other_regions_called"] == 1, pd)
        check("paralogs: SMN1 Paraphase 2 against Parascopy 1 is a disagreement",
              pd["data"]["smn1"]["paraphase"] == 2 and pd["data"]["smn1"]["parascopy"]["cn"] == 1
              and pd["data"]["smn1"]["agree"] == "disagree", pd["data"]["smn1"])
        check("paralogs: the text report lists HBA2 1 with the deletion, and the disagreement",
              "  HBA2: 1 (called); genotype -a/aa (3 alpha-globin genes); 4.2 kb deletion (-a4.2): HBA2 lost" in txt_p
              and "SMN1 by both tools: Paraphase 2, Parascopy 1 (quality 33, filter PASS): they DISAGREE" in txt_p
              and "Variants on phased haplotypes: 1" in txt_p, txt_p)
        card_p = html_p.split("<h2>Paralog Genes (Paraphase)</h2>", 1)[-1].split("<h2>", 1)[0]
        check("paralogs: the HTML card has the HBA2 row and the SMN1 line",
              "<td>HBA2</td><td>1</td>" in card_p and "they DISAGREE" in card_p and "PacBio HiFi reads only" in card_p,
              card_p)
        check("paralogs summary validates", not validate.validate(schema, json.loads(json.dumps(summ_p))))
        put(f"{q}/paralogs/S_smn_copy_number.tsv", "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\t"
            "psCN_filter\tpsCN\tpsCN_qual\thomologous_regions\nchr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\t"
            "PASS\t2,2\t33,30\tchr5:70049000-70077000:-\n")
        check("paralogs: Parascopy 2 agrees", "Parascopy 2 (quality 33, filter PASS): the two agree" in render(q)[1])
        os.remove(f"{q}/paralogs/S_smn_copy_number.tsv")
        _, txt_p2, _ = render(q)
        check("paralogs: without step 35's table there is no SMN1 comparison line", "SMN1 by both tools" not in txt_p2,
              txt_p2)
        check("paralogs: a folder without step 35b lists it under Steps Not Run",
              "  - Paralog genes (Paraphase, PacBio HiFi)" in render(os.path.join(work, "mito", "S"))[1])

        # 14. a step that failed in the latest run
        f_ = os.path.join(work, "failed", "S")
        put(f"{f_}/clinvar/S_clinvar_hits.vcf", HITS)
        put(f"{f_}/logs/run_status.tsv", f"meta\tstarted_epoch\t{time.time() - DAY}\nstep\t10\tfailed\nstep\t11\tok\n"
                                          "step\t32\tfailed\n")
        summ_f, txt_f, html_f = render(f_)
        tel = summ_f["sections"]["telomere"]
        check("failed: no telomere file and step 10 failed -> state failed with a note",
              tel["state"] == "failed" and "step 10 failed" in (tel["note"] or ""), tel)
        check("failed: a step that never ran is still missing", summ_f["sections"]["roh"]["state"] == "missing")
        card_t = html_f.split("<h2>Telomere content (relative)</h2>", 1)[-1].split("<h2>", 1)[0]
        check("failed: the HTML card says Failed, not Not run", 'badge-red">Failed<' in card_t and "Not run" not in card_t,
              card_t)
        check("failed: the text report says FAILED and lists it under Steps Not Run",
              "[FAILED: step 10 failed" in txt_f and "  - Telomere content (TelomereHunter) (step 10 failed)" in txt_f, txt_f)
        card_c = html_f.split("<h2>CYP2D6 Across Callers</h2>", 1)[-1].split("<h2>", 1)[0]
        check("failed: a failed pypgx step shows as failed in the CYP2D6 rows of both reports, with the FAILED note",
              summ_f["cyp2d6"]["calls"]["pypgx"] == "failed (step 32)" and "failed (step 32)" in txt_f
              and "failed (step 32)" in card_c and "FAILED: step 32 failed" in card_c, card_c)
        check("failed summary validates", not validate.validate(schema, json.loads(json.dumps(summ_f))),
              validate.validate(schema, json.loads(json.dumps(summ_f))))

        # 16. mosaic loss of chrY or chrX, from the step 33 table
        mo = os.path.join(work, "mosaic", "S")
        base = "key\tvalue\nsample\tS\ninferred_sex\tmale\nsex_check\tok\n"
        put(f"{mo}/qc/S_sample_qc.tsv", base)
        sec_m = collect_summary.collect("S", mo)["sections"]["mosaic_sex_loss"]
        check("mosaic: a table without the rows leaves the section missing", sec_m["state"] == "missing", sec_m)
        rows = {"sex": "male", "sex_from": "reads", "chrom": "chrY", "depth_ratio": "0.351",
                "depth_ratio_expected": "0.5", "depth_blocks": "109", "loss_depth": "0.298",
                "loss_depth_low": "0.293", "loss_depth_high": "0.303", "x_het_sites": "", "loss_af": "",
                "loss_af_low": "", "loss_af_high": "", "floor": "0.100", "call": "detected",
                "reason": "the low end of the 95% range is above 10%: a share of the cells lost the Y"}
        put(f"{mo}/qc/S_sample_qc.tsv", base + "".join(f"mosaic_{k}\t{v}\n" for k, v in rows.items()))
        summ_m = collect_summary.collect("S", mo)
        sec_m = summ_m["sections"]["mosaic_sex_loss"]
        check("mosaic: the rows make an ok section from step 33", sec_m["state"] == "ok" and sec_m["step"] == "33"
              and sec_m["source"] == "qc/S_sample_qc.tsv", sec_m)
        check("mosaic: numbers are numbers, blanks are None, the call and reason stay text",
              sec_m["data"]["loss_depth"] == 0.298 and sec_m["data"]["loss_af"] is None
              and sec_m["data"]["x_het_sites"] is None and sec_m["data"]["call"] == "detected"
              and sec_m["data"]["chrom"] == "chrY", sec_m["data"])
        check("mosaic: the sample QC section still reads the same table",
              summ_m["sections"]["sample_qc"]["data"].get("sex_check") == "ok", summ_m["sections"]["sample_qc"])
        check("mosaic summary validates", not validate.validate(schema, json.loads(json.dumps(summ_m))),
              validate.validate(schema, json.loads(json.dumps(summ_m))))
        put(f"{mo}/qc/S_sample_qc.tsv", base + "".join(f"mosaic_{k}\t{v}\n" for k, v in
                                                       dict(rows, loss_depth="inf", loss_depth_high="-inf").items()))
        summ_i = collect_summary.collect("S", mo)
        sec_m = summ_i["sections"]["mosaic_sex_loss"]
        check("mosaic: inf is no value (JSON null), never infinity", sec_m["data"]["loss_depth"] is None
              and sec_m["data"]["loss_depth_high"] is None and "Infinity" not in json.dumps(summ_i), sec_m["data"])
        put(f"{mo}/qc/S_sample_qc.tsv", base + "mosaic_call\tmaybe\n")
        sec_m = collect_summary.collect("S", mo)["sections"]["mosaic_sex_loss"]
        check("mosaic: an unknown mosaic_call is unreadable", sec_m["state"] == "unreadable"
              and "mosaic_call" in (sec_m["note"] or ""), sec_m)
        # 18. the traits panel (step 42)
        import traits_report
        t = os.path.join(work, "traits", "S")
        rows = traits_report.read_table(os.path.join(REPO, "assets", "traits.tsv"))
        by_id = {r["rsid"]: r for r in rows}
        lct, ear = by_id["rs4988235"], by_id["rs17822931"]
        g = os.path.join(work, "traits_g.tsv")
        put(g, f"{lct['chrom']}\t{lct['pos']}\t{lct['ref']}\t{lct['alt']}\t{lct['ref']}/{lct['alt']}\n"
               f"{ear['chrom']}\t{ear['pos']}\t{ear['ref']}\t.\t{ear['ref']}/{ear['ref']}\n"
               "chr19\t44908684\tT\t.\tT/T\nchr19\t44908822\tC\t.\tC/C\n")
        os.makedirs(f"{t}/traits", exist_ok=True)
        traits_report.write_tsv(f"{t}/traits/S_traits.tsv", traits_report.report(rows, traits_report.read_genotypes(g)))
        summ_t, txt_t, html_t = render(t)
        td = summ_t["sections"]["traits"]
        check("traits: section ok, counts without the APOE rows",
              td["state"] == "ok" and (td["data"]["called"], td["data"]["not_covered"]) == (2, len(rows) - 4), td)
        check("traits: the APOE rows and the APOE line are kept apart",
              len(td["data"]["apoe"]) == 3 and not any(r["gene"] == "APOE" for r in td["data"]["rows"]), td["data"]["apoe"])
        check("traits: the text report has the read row with its genotype, reading and PMID, marked ancestry",
              f"(MCM6 rs4988235): {lct['ref']}/{lct['alt']} [ancestry]. {lct['het']} (PMID {lct['pmids']})" in txt_t, txt_t)
        check("traits: the text report lists a not-covered row as not read",
              "Not read (no genotype in the gVCF): " in txt_t and "rs671" in txt_t.split("Not read (no genotype in the gVCF): ", 1)[-1])
        apoe_at = txt_t.find(render_report.APOE_WARNING)
        check("traits: the APOE line comes after the warning in the text report, e3/e3",
              apoe_at > 0 and "(rs429358+rs7412): e3/e3." in txt_t[apoe_at:] and "e3/e3" not in txt_t[:apoe_at], txt_t)
        card_tr = html_t.split("<h2>Traits</h2>", 1)[-1].split("<h2>", 1)[0]
        check("traits: the HTML card has the row and the ancestry badge",
              f"rs4988235</span></td><td>{lct['ref']}/{lct['alt']} " in card_tr and 'badge-yellow">ancestry<' in card_tr, card_tr)
        folded = card_tr.split("<details>", 1)
        check("traits: in the HTML the APOE lines are inside the folded block only",
              len(folded) == 2 and "rs429358+rs7412" in folded[1] and "rs429358" not in folded[0], card_tr)
        check("traits: the summary validates", not validate.validate(schema, json.loads(json.dumps(summ_t))))
        put(f"{t}/traits/S_traits.tsv", "rsid\tgene\tstatus\nrs1\tX\tcalled\n")
        check("traits: a table without the reading column is unreadable",
              collect_summary.collect("S", t)["sections"]["traits"]["state"] == "unreadable")
        put(f"{t}/traits/S_traits.tsv", "rsid\tgene\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids\n"
            "rs1\tX\tT\tA/A\tcalled\n")
        check("traits: a row cut short after its status is unreadable, never a called trait without a reading",
              collect_summary.collect("S", t)["sections"]["traits"]["state"] == "unreadable")
        os.remove(f"{t}/traits/S_traits.tsv")
        summ_m, _, html_m = render(t)
        check("traits: no file -> missing, listed under Steps Not Run, the card says Not run",
              summ_m["sections"]["traits"]["state"] == "missing" and "Traits panel" in summ_m["not_run"]
              and ">Not run<" in html_m.split("<h2>Traits</h2>", 1)[-1].split("<h2>", 1)[0], summ_m["not_run"])
        # 19. step 40, the raw-data export
        e_ = os.path.join(work, "export", "S")
        put(f"{e_}/consumer_export/S_consumer_export.tsv", "key\tvalue\nsites\t10\ncalled\t7\nno_call\t3\n"
            "from_reference_blocks\t5\ntemplate_sites\t12\nfile_23andme\tS_23andme_v5.zip\n"
            "file_ancestrydna\tS_ancestrydna_layout.zip\nno_call: no-call in the gVCF\t3\n")
        ce = collect_summary.collect("S", e_)["sections"]["consumer_export"]
        check("step 40: counts and reasons read", ce["state"] == "ok" and ce["step"] == "40"
              and ce["data"].get("called") == 7 and ce["data"].get("sites") == 10
              and ce["data"].get("no_call_reasons") == {"no-call in the gVCF": 3}, ce)
        put(f"{e_}/consumer_export/S_consumer_export.tsv", "key\tvalue\nsites\t10\ncalled\t\nno_call\t3\n")
        ce = collect_summary.collect("S", e_)["sections"]["consumer_export"]
        check("step 40: a table without a whole-number 'called' row is unreadable", ce["state"] == "unreadable", ce)
        put(f"{e_}/consumer_export/S_consumer_export.tsv", "key\tvalue\nsites\t10\ncalled\t7\nno_call\t3\n"
            "no_call: filtered call\tthree\n")
        ce = collect_summary.collect("S", e_)["sections"]["consumer_export"]
        check("step 40: a '--' reason without a whole-number count is unreadable, not left out",
              ce["state"] == "unreadable" and "no_call: filtered call" in (ce["note"] or ""), ce)
        put(f"{e_}/consumer_export/S_consumer_export.tsv", "key\tvalue\nsites\t10\ncalled\t7\nno_call\t3\n", age_days=2)
        put(f"{e_}/logs/run_status.tsv", f"meta\tstarted_epoch\t{time.time() - DAY}\nstep\t40\tskipped (opt-in: add consumer_export to TOOLS)\n")
        summ_e = collect_summary.collect("S", e_)
        ce = summ_e["sections"]["consumer_export"]
        check("step 40: an older table whose step was skipped is stale", ce["state"] == "stale", ce)
        check("step 40 summary validates", not validate.validate(schema, json.loads(json.dumps(summ_e))),
              validate.validate(schema, json.loads(json.dumps(summ_e))))

        # 20. step 39, the couple table
        import couple_report
        c_ = os.path.join(work, "couple", "B")
        head = "\t".join(couple_report.TABLE_COLUMNS) + "\n"

        def crow(gene, result, risk, fa, fb):
            r = {"gene": gene, "inheritance": "AR", "lists": "ACMG", "condition": f"{gene} disease", "result": result,
                 "offspring_risk": risk, "partner_a": "A", "findings_a": fa, "partner_b": "B", "findings_b": fb,
                 "explanation": "x"}
            return "\t".join(r[k] for k in couple_report.TABLE_COLUMNS) + "\n"
        put(f"{c_}/couple/B_couple_carrier.tsv", head + crow("NAGA", "shared", "1 in 4 per pregnancy", "a-snv", "b-snv")
            + crow("FMR1", "x_linked", "A: sons 1 in 2", "a-repeat", ".") + crow("UGT1A1", "one_partner", "none shared", ".", "b-only"))
        summ_c = collect_summary.collect("B", c_)
        cs = summ_c["sections"]["couple"]
        check("step 39: partner A, one shared gene with B's own findings first, one X-linked, one partner-only",
              cs["state"] == "ok" and cs["step"] == "39" and cs["data"]["partner"] == "A"
              and [x["gene"] for x in cs["data"]["shared"]] == ["NAGA"] and cs["data"]["shared"][0]["own"] == "b-snv"
              and cs["data"]["shared"][0]["partner_findings"] == "a-snv" and cs["data"]["x_linked"][0]["gene"] == "FMR1"
              and cs["data"]["one_partner"] == 1, cs)
        txt_c, html_c = render_report.text_report(summ_c), render_report.html_report(summ_c)
        check("step 39: both reports give the shared gene with its chance and the counselling note",
              all("Shared: NAGA (NAGA disease): 1 in 4 per pregnancy" in t and "genetic counsellor" in t for t in (txt_c, html_c)),
              txt_c)
        check("step 39 summary validates", not validate.validate(schema, json.loads(json.dumps(summ_c))),
              validate.validate(schema, json.loads(json.dumps(summ_c))))
        put(f"{c_}/couple/B_couple_carrier.tsv", head + crow("NAGA", "shared", "1 in 4", "x", "y").replace("\tB\t", "\tC\t"))
        cs = collect_summary.collect("B", c_)["sections"]["couple"]
        check("step 39: a table of another couple is unreadable", cs["state"] == "unreadable", cs)
        put(f"{c_}/couple/B_couple_carrier.tsv", "gene\tresult\nNAGA\tshared\n")
        cs = collect_summary.collect("B", c_)["sections"]["couple"]
        check("step 39: a table without the columns is unreadable", cs["state"] == "unreadable", cs)
        put(f"{c_}/couple/B_couple_carrier.tsv", head)
        cs = collect_summary.collect("B", c_)["sections"]["couple"]
        check("step 39: a table without rows names no partner and counts nothing",
              cs["state"] == "ok" and cs["data"]["partner"] is None and cs["data"]["shared"] == [], cs)
        check("step 39: a single sample without the step has no not-run line for it",
              render(d)[0]["sections"]["couple"]["state"] == "missing"
              and not any("Couple" in x for x in render(d)[0]["not_run"]), render(d)[0]["not_run"])
        clinvar_quality(work, schema)
        cpsr_beside_clinvar(work, schema)
        pypgx_section(work, schema)
        clinical_mane_quality(work, schema)
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
