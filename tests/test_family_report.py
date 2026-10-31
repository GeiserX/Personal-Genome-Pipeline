#!/usr/bin/env python3
"""bin/family_report.py (step 38 and the Nextflow FAMILY workflow) on
synthetic pedigrees, somalier tables and joint VCFs.

Checks:
  1. prepare: a trio with a grandparent named but not sequenced writes
     members.tsv, somalier's pedigree (sequenced members only) and WhatsHap's
     trios; a pedigree problem stops it with exit 2 and names the problem
     (no header, an unknown column, a female father, a person who is their
     own parent, a loop, a second row, a sex the samplesheet contradicts, a
     gVCF whose sample is another, fewer than two sequenced samples);
  2. the expected relatedness of parent and child, full and half siblings,
     grandparent, first cousins and partners;
  3. check: somalier numbers that fit the trio pass; the same numbers with
     two names swapped in the pedigree exit 1 and name the pair; too few
     sites is not_checked, not a pass and not a stop; a first-degree pair
     with many opposite homozygotes is not parent and child;
  4. report: the allele from each parent where the genotypes tell, an
     all-heterozygous site left ambiguous, a planted new allele listed as a
     de novo candidate (and an opposite-homozygote site not), a son's chrX
     outside the PARs taken from the mother only, inside a PAR from both,
     the allele of an unsequenced father inferred from a duo, compound
     heterozygotes in trans and in cis, WhatsHap's phase compared with the
     genotype origin, and the same bytes on a second run.

Run: python3 tests/test_family_report.py
"""
import gzip
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "family_report.py")
PAR = os.path.join(REPO, "assets", "par_grch38.bed")
sys.path.insert(0, os.path.join(REPO, "bin"))
import family_report as fr  # noqa: E402

FAILS = []


def check(desc, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + desc + ("" if ok else "  " + detail))
    if not ok:
        FAILS.append(desc)


def run(*args):
    p = subprocess.run([sys.executable, SCRIPT] + list(args), capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


def gvcf(d, name, sample):
    path = os.path.join(d, name)
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
        f.write("chr1\t1\t.\tA\t<*>\t0\t.\tEND=10\tGT:GQ\t0/0:50\n")
    return path


def tsv(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f]


TRIO = ("sample\tsex\tfather\tmother\tpartner\n"
        "kid\tmale\tdad\tmum\t.\n"
        "dad\tmale\tgrandpa\t.\tmum\n"
        "mum\tfemale\t.\t.\tdad\n")


def prepare(d, ped_text, samples=("kid", "dad", "mum"), extra=()):
    ped = write(os.path.join(d, "ped.tsv"), ped_text)
    args = ["prepare", "--pedigree", ped, "--out-dir", os.path.join(d, "out")]
    for s in samples:
        args += ["--gvcf", "%s=%s" % (s, gvcf(d, s + ".g.vcf.gz", s))]
    return run(*(args + list(extra)))


def test_prepare():
    with tempfile.TemporaryDirectory() as d:
        rc, out = prepare(d, TRIO)
        check("prepare: a trio passes (exit %d)" % rc, rc == 0, out)
        m = {r["sample"]: r for r in tsv(os.path.join(d, "out", "members.tsv"))}
        check("members.tsv names the grandparent, not sequenced, male",
              m.get("grandpa", {}).get("sequenced") == "no" and m["grandpa"]["sex"] == "male"
              and m["grandpa"]["in_pedigree"] == "no", str(m.get("grandpa")))
        check("members.tsv: the three are sequenced",
              all(m[s]["sequenced"] == "yes" for s in ("kid", "dad", "mum")), str(m))
        ped = open(os.path.join(d, "out", "family.ped")).read().splitlines()
        check("somalier's pedigree: three rows, the kid with both parents, the unsequenced grandpa as -9",
              sorted(ped) == sorted(["family\tdad\t-9\t-9\t1\t-9", "family\tkid\tdad\tmum\t1\t-9",
                                     "family\tmum\t-9\t-9\t2\t-9"]), str(ped))
        trios = open(os.path.join(d, "out", "trios.ped")).read().splitlines()
        check("WhatsHap's trios: the kid", trios == ["family\tkid\tdad\tmum\t1\t0"], str(trios))
        check("prepare says what it will check", "kid and mum: declared parent-child" in out, out)

        # a duo: no trio for WhatsHap
        rc, out = prepare(d, TRIO, samples=("kid", "mum"))
        check("prepare: a duo passes", rc == 0, out)
        check("a duo has no trio", open(os.path.join(d, "out", "trios.ped")).read() == "", out)

    bad = [
        ("no header", "kid\tmale\tdad\tmum\t.\n", "header with a 'sample' column"),
        ("an unknown column", "sample\tsex\tdad\n" "kid\tmale\tx\n", "unknown column dad"),
        ("a female father", TRIO.replace("dad\tmale", "dad\tfemale"), "'dad' is the father of 'kid' but its sex is female"),
        ("own parent", TRIO.replace("kid\tmale\tdad", "kid\tmale\tkid"), "'kid' is named as its own father"),
        ("a loop", TRIO.replace("dad\tmale\tgrandpa", "dad\tmale\tkid"), "among its own ancestors"),
        ("a second row", TRIO + "kid\tmale\t.\t.\t.\n", "sample 'kid' has a second row"),
        ("same father and mother", TRIO.replace("kid\tmale\tdad\tmum", "kid\tmale\tdad\tdad"), "same person, 'dad', as father and mother"),
        ("a bad sex", TRIO.replace("kid\tmale", "kid\tboy"), "sex 'boy'"),
        ("partners that disagree", TRIO.replace("mum\tfemale\t.\t.\tdad", "mum\tfemale\t.\t.\tkid"), "'dad' names 'mum' as partner, but 'mum' names 'kid'"),
    ]
    for desc, text, msg in bad:
        with tempfile.TemporaryDirectory() as d:
            rc, out = prepare(d, text)
            check("prepare stops on %s (exit %d) and says so" % (desc, rc), rc == 2 and msg in out, out)
    with tempfile.TemporaryDirectory() as d:
        rc, out = prepare(d, TRIO, extra=["--declared-sex", "mum=male"])
        check("prepare stops when the samplesheet's sex contradicts the pedigree's",
              rc == 2 and "the pedigree says 'mum' is female, the samplesheet says male" in out, out)
        rc, out = prepare(d, TRIO, samples=("kid",))
        check("prepare stops with one sequenced sample", rc == 2 and "needs two or more" in out, out)
        ped = write(os.path.join(d, "ped.tsv"), TRIO)
        rc, out = run("prepare", "--pedigree", ped, "--out-dir", os.path.join(d, "o2"),
                      "--gvcf", "kid=" + gvcf(d, "a.g.vcf.gz", "kid"), "--gvcf", "dad=" + gvcf(d, "b.g.vcf.gz", "HG003"))
        check("prepare stops on a gVCF whose sample is another, with the rename command",
              rc == 2 and "holds the sample HG003, not 'dad'" in out and "bcftools reheader" in out, out)


def test_kinship():
    ped = {
        "gp": (None, None), "gm": (None, None),
        "f": ("gp", "gm"), "u": ("gp", "gm"), "m": (None, None), "w": (None, None),
        "a": ("f", "m"), "b": ("f", "m"), "h": ("f", "w"), "c": ("u", "w2"), "w2": (None, None),
    }
    members = {k: {"father": v[0], "mother": v[1], "sex": None, "partner": None} for k, v in ped.items()}
    members["m"]["partner"] = "f"
    phi = fr.kinship_table(members)
    cases = [("a", "f", 0.5, "parent-child"), ("a", "b", 0.5, "full siblings"), ("a", "h", 0.25, "half siblings"),
             ("a", "gp", 0.25, "second degree"), ("a", "c", 0.125, "third degree"), ("f", "m", 0.0, "partners"),
             ("m", "w", 0.0, "unrelated")]
    for x, y, want, label in cases:
        got = 2 * phi(x, y)
        lab = fr.declared_relationship(x, y, members, got)
        check("expected relatedness of %s and %s is %.3f, %s (got %.3f, %s)" % (x, y, want, label, got, lab),
              abs(got - want) < 1e-9 and lab == label)


PAIRS_HEAD = ("#sample_a\tsample_b\trelatedness\tibs0\tibs2\tconcordance\thets_a\thets_b\thets_ab\tshared_hets"
              "\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n")


def pairs_file(path, rows):
    with open(path, "w") as f:
        f.write(PAIRS_HEAD)
        for a, b, rel, ibs0, n in rows:
            f.write("\t".join([a, b, str(rel), str(ibs0), "100", "0.9", "1", "1", "1", "1", "1", "1", "1", str(n),
                               "0", "0", "-1"]) + "\n")
    return path


def test_check():
    data = [("dad", "kid", 0.49, 2, 5000), ("kid", "mum", 0.51, 1, 5000), ("dad", "mum", 0.01, 600, 5000)]
    with tempfile.TemporaryDirectory() as d:
        prepare(d, TRIO)
        members = os.path.join(d, "out", "members.tsv")
        pairs = pairs_file(os.path.join(d, "pairs.tsv"), data)
        out_tsv = os.path.join(d, "rel.tsv")
        rc, out = run("check", "--members", members, "--pairs", pairs, "--out", out_tsv)
        rows = tsv(out_tsv)
        check("check: the trio's numbers pass (exit %d)" % rc, rc == 0 and all(r["verdict"] == "ok" for r in rows), out)
        check("check writes the three pairs", len(rows) == 3, str(rows))

        # The same data, the pedigree with kid and dad swapped: dad is now
        # the child of kid and mum, and dad-mum (unrelated) is declared parent-child.
        swapped = TRIO.replace("kid\tmale\tdad\tmum", "@").replace("dad\tmale\tgrandpa\t.\tmum", "kid\tmale\tgrandpa\t.\tmum") \
            .replace("@", "dad\tmale\tkid\tmum").replace("mum\tfemale\t.\t.\tdad", "mum\tfemale\t.\t.\tkid")
        d2 = os.path.join(d, "swap")
        os.makedirs(d2)
        rc0, out0 = prepare(d2, swapped)
        check("the swapped pedigree is itself valid", rc0 == 0, out0)
        rc, out = run("check", "--members", os.path.join(d2, "out", "members.tsv"), "--pairs", pairs,
                      "--out", os.path.join(d2, "rel.tsv"))
        check("check: two names swapped in the pedigree stop it (exit %d)" % rc, rc == 1, out)
        check("and the message names the pair that disagrees",
              "dad and mum: declared parent-child (first degree), but somalier measures unrelated" in out, out)
        check("and says what it may be", "sample swap" in out, out)

        few = pairs_file(os.path.join(d, "few.tsv"), [(a, b, r, i, 300) for a, b, r, i, _ in data])
        rc, out = run("check", "--members", members, "--pairs", few, "--out", out_tsv)
        rows = tsv(out_tsv)
        check("check: 300 sites is not_checked and does not stop (exit %d)" % rc,
              rc == 0 and all(r["verdict"] == "not_checked" for r in rows), out)
        sib = pairs_file(os.path.join(d, "sib.tsv"), [("dad", "kid", 0.5, 400, 5000)] + data[1:])
        rc, out = run("check", "--members", members, "--pairs", sib, "--out", out_tsv)
        check("check: first degree with 8%% opposite homozygotes is not parent and child (exit %d)" % rc,
              rc == 1 and "siblings do not" in out, out)
        lost = pairs_file(os.path.join(d, "lost.tsv"), data[:2])
        rc, out = run("check", "--members", members, "--pairs", lost, "--out", out_tsv)
        check("check: a pair somalier did not compare stops it", rc == 2 and "no row for dad and mum" in out, out)


HEAD = ("##fileformat=VCFv4.2\n"
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. "
        "Format: Allele|Consequence|IMPACT|SYMBOL\">\n")


def joint(path, rows, samples=("dad", "kid", "mum"), phased=False):
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(samples))
        for chrom, pos, ref, alt, gts in rows:
            fmt = "GT:DP:AD:GQ" + (":PS" if phased else "")
            cells = []
            for g in gts:
                cell = g if ":" in g else "%s:30:15,15:50" % g
                cells.append(cell + (":1" if phased else ""))
            f.write("\t".join([chrom, str(pos), ".", ref, alt, "50", ".", ".", fmt] + cells) + "\n")
    return path


def read_inh(path):
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        return {(r[0], int(r[1])): dict(zip(head, r)) for r in (l.rstrip("\n").split("\t") for l in f)}


def test_report():
    #            dad      kid      mum
    rows = [
        ("chr1", 100, "A", "G", ["0/0", "0/1", "0/1"]),          # from mum
        ("chr1", 200, "C", "T", ["0/1", "0/1", "0/1"]),          # ambiguous
        ("chr1", 300, "G", "A", ["0/0:30:30,0:50", "0/1:30:14,16:50", "0/0:30:30,0:50"]),  # planted de novo
        ("chr1", 400, "T", "C", ["1/1", "0/0", "0/0"]),          # opposite homozygote: inconsistent, not de novo
        ("chr1", 500, "A", "C", ["0/1", "1/1", "0/1"]),          # both alt
        ("chr2", 100, "A", "G", ["0/1", "0/1", "0/0"]),          # gene G1, from dad
        ("chr2", 200, "C", "T", ["0/0", "0/1", "0/1"]),          # gene G1, from mum -> in trans
        ("chr3", 100, "A", "G", ["0/0", "0/1", "0/1"]),          # gene G2, from mum
        ("chr3", 200, "C", "T", ["0/0", "0/1", "1/1"]),          # gene G2, from mum -> in cis
        ("chrX", 74000000, "G", "A", ["0/0", "1/1", "0/1"]),     # son's X: from mum only
        ("chrX", 74000100, "G", "A", ["1/1", "0/0", "0/1"]),     # son's X: dad's allele does not count
        ("chrX", 500000, "C", "T", ["0/1", "0/1", "0/0"]),       # PAR1: autosomal, from dad
        ("chrM", 100, "A", "G", ["0/0", "0/1", "0/0"]),          # left out
    ]
    with tempfile.TemporaryDirectory() as d:
        prepare(d, TRIO)
        out = os.path.join(d, "out")
        members = os.path.join(out, "members.tsv")
        pairs = pairs_file(os.path.join(d, "pairs.tsv"),
                           [("dad", "kid", 0.49, 2, 5000), ("kid", "mum", 0.51, 1, 5000), ("dad", "mum", 0.01, 600, 5000)])
        rel = os.path.join(out, "relationships.tsv")
        run("check", "--members", members, "--pairs", pairs, "--out", rel)
        vcf = joint(os.path.join(d, "family.vcf.gz"), rows)
        cand = write(os.path.join(d, "kid_clinical.vcf"), HEAD +
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkid\n"
                     "chr2\t100\t.\tA\tG\t50\tPASS\tCSQ=G|missense_variant|MODERATE|G1\tGT\t0/1\n"
                     "chr2\t200\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|G1\tGT\t0/1\n"
                     "chr3\t100\t.\tA\tG\t50\tPASS\tCSQ=G|missense_variant|MODERATE|G2\tGT\t0/1\n"
                     "chr3\t200\t.\tC\tT\t50\tPASS\tCSQ=T|missense_variant|MODERATE|G2\tGT\t0/1\n")
        args = ["report", "--members", members, "--relationships", rel, "--vcf", vcf, "--candidates",
                "kid=" + cand, "--par-bed", PAR, "--out-dir", out]
        rc, log = run(*args)
        check("report runs (exit %d)" % rc, rc == 0, log)
        inh = read_inh(os.path.join(out, "kid_inheritance.tsv.gz"))
        def got(k):
            r = inh.get(k, {})
            return (r.get("from_father"), r.get("from_mother"), r.get("how"))
        check("chr1:100 the G came from the mother", got(("chr1", 100)) == ("A", "G", "genotype"), str(got(("chr1", 100))))
        check("chr1:200 all heterozygous: ambiguous", got(("chr1", 200)) == (".", ".", "ambiguous"), str(got(("chr1", 200))))
        check("chr1:400 opposite homozygote: inconsistent", got(("chr1", 400))[2] == "inconsistent", str(got(("chr1", 400))))
        check("chr1:500 the alt from both", got(("chr1", 500)) == ("C", "C", "genotype"), str(got(("chr1", 500))))
        check("chrX:74000000 the son's X from the mother only", got(("chrX", 74000000)) == (".", "A", "genotype"),
              str(got(("chrX", 74000000))))
        check("chrX:74000100 the father's X allele does not count for a son", got(("chrX", 74000100)) == (".", "G", "genotype"),
              str(got(("chrX", 74000100))))
        check("chrX:500000 (PAR1) from the father, as an autosome", got(("chrX", 500000)) == ("T", "C", "genotype"),
              str(got(("chrX", 500000))))
        check("chrM is left out", ("chrM", 100) not in inh, str(inh.get(("chrM", 100))))
        dn = tsv(os.path.join(out, "de_novo.tsv"))
        check("the planted new allele is the one de novo candidate, quality pass",
              [(r["chrom"], r["pos"], r["alt"], r["quality"]) for r in dn] == [("chr1", "300", "A", "pass")], str(dn))
        ch = tsv(os.path.join(out, "compound_het.tsv"))
        rel_of = {r["gene"]: (r["a_from"], r["b_from"], r["phase"]) for r in ch}
        check("G1: one from each parent, in trans", rel_of.get("G1") == ("father", "mother", "in trans"), str(ch))
        check("G2: both from the mother, in cis", rel_of.get("G2") == ("mother", "mother", "in cis"), str(ch))
        md = open(os.path.join(out, "family_report.md")).read()
        check("the report states the de novo count and the pairs", "De novo candidates (an allele neither parent has): 1, 1" in md
              and "in trans (one from each parent, both copies of the gene hit) 1, in cis (both from one parent) 1" in md, md)

        first = open(os.path.join(out, "kid_inheritance.tsv.gz"), "rb").read()
        run(*args)
        check("a second run writes the same bytes", open(os.path.join(out, "kid_inheritance.tsv.gz"), "rb").read() == first)

        # A low-quality copy of the planted site: still listed, quality low
        low = joint(os.path.join(d, "low.vcf.gz"), [("chr1", 300, "G", "A", ["0/0:5:5,0:10", "0/1:30:14,16:50", "0/0:30:30,0:50"])])
        run("report", "--members", members, "--relationships", rel, "--vcf", low, "--par-bed", PAR, "--out-dir", os.path.join(d, "low"))
        dn = tsv(os.path.join(d, "low", "de_novo.tsv"))
        check("a candidate with a parent at depth 5 is listed as low", [r["quality"] for r in dn] == ["low"], str(dn))

        # WhatsHap's phase: kid's first haplotype holds the father's allele at chr2:100
        ph = joint(os.path.join(d, "phased.vcf.gz"),
                   [("chr2", 100, "A", "G", ["0|1", "1|0", "0|0"]), ("chr2", 200, "C", "T", ["0|0", "0|1", "0|1"])], phased=True)
        rc, log = run("report", "--members", members, "--relationships", rel, "--vcf", vcf, "--phased", ph,
                      "--par-bed", PAR, "--out-dir", os.path.join(d, "ph"))
        md = open(os.path.join(d, "ph", "family_report.md")).read()
        check("WhatsHap's phase agrees with the genotype origin at both sites",
              "WhatsHap phased 2 heterozygous sites of kid by inheritance. At 2 of them the genotypes also give the "
              "origin, and the phase agrees with it at 2 (100.0%)" in md, md)

    # A duo: the father is not sequenced; his allele is inferred
    with tempfile.TemporaryDirectory() as d:
        prepare(d, TRIO, samples=("kid", "mum"))
        out = os.path.join(d, "out")
        pairs = pairs_file(os.path.join(d, "pairs.tsv"), [("kid", "mum", 0.5, 1, 5000)])
        rel = os.path.join(out, "relationships.tsv")
        rc, log = run("check", "--members", os.path.join(out, "members.tsv"), "--pairs", pairs, "--out", rel)
        check("a duo passes the check", rc == 0, log)
        vcf = joint(os.path.join(d, "duo.vcf.gz"),
                    [("chr1", 100, "A", "G", ["0/1", "0/0"]), ("chr1", 200, "A", "G", ["0/1", "0/1"]),
                     ("chr1", 300, "A", "G", ["0/1", "1/1"])], samples=("kid", "mum"))
        rc, log = run("report", "--members", os.path.join(out, "members.tsv"), "--relationships", rel, "--vcf", vcf,
                      "--par-bed", PAR, "--out-dir", out)
        inh = read_inh(os.path.join(out, "kid_inheritance.tsv.gz"))
        check("duo: the G the mother lacks came from the unsequenced father",
              (inh[("chr1", 100)]["from_father"], inh[("chr1", 100)]["father_gt"]) == ("G", "."), str(inh[("chr1", 100)]))
        check("duo: child and mother both heterozygous is ambiguous", inh[("chr1", 200)]["how"] == "ambiguous")
        md = open(os.path.join(out, "family_report.md")).read()
        check("duo: the report states the inferred share of the father, and calls no de novo",
              "The father is not sequenced. At 2 of these sites (66.7%)" in md and "needs both parents sequenced" in md, md)


FIX = os.path.join(REPO, "tests", "fixtures", "families")
TWO = os.path.join(FIX, "two_families.tsv")


def plan(*args):
    return run("plan", *args)


def plan_rows(path):
    return {(r["analysis"], r["members"]): r for r in tsv(path)}


def by_analysis(path):
    out = {}
    for r in tsv(path):
        out.setdefault(r["analysis"], []).append(r)
    return out


def test_families_file():
    with tempfile.TemporaryDirectory() as d:
        # the two-family file splits into two pedigrees
        for fam, want in (("family_a", ["child_a", "father_a", "mother_a"]), ("family_b", ["partner_x", "partner_y"])):
            rc, out = plan("--families", TWO, "--family", fam, "--write-pedigree", d)
            ped = os.path.join(d, fam, "pedigree.tsv")
            rows = tsv(ped) if os.path.exists(ped) else []
            check("plan --write-pedigree: %s gets its own pedigree (exit %d)" % (fam, rc),
                  rc == 0 and [r["sample"] for r in rows] == want, out + str(rows))
            check("%s's pedigree has the pedigree columns, no family column" % fam,
                  rows and list(rows[0]) == list(fr.PED_COLUMNS), str(rows[:1]))
        rows, errors = fr.read_pedigree(os.path.join(d, "family_a", "pedigree.tsv"))
        members, more = fr.build_members(rows, {}, "x")
        check("the written pedigree reads back: child_a's parents are father_a and mother_a",
              not errors + more and (members["child_a"]["father"], members["child_a"]["mother"]) == ("father_a", "mother_a"),
              str(errors + more))

        # a file without the column gives the same members.tsv as before the column existed
        args = ["prepare", "--pedigree", os.path.join(REPO, "assets", "pedigree.example.tsv"), "--out-dir",
                os.path.join(d, "old")]
        for s in ("child_a", "father_a", "mother_a"):
            args += ["--gvcf", "%s=%s" % (s, gvcf(d, s + ".g.vcf.gz", s))]
        rc, out = run(*args)
        got = open(os.path.join(d, "old", "members.tsv"), "rb").read() if rc == 0 else b""
        check("a pedigree without the family column: members.tsv is byte for byte the golden file",
              got == open(os.path.join(FIX, "members.golden.tsv"), "rb").read(), out + repr(got))
        # the same family from the families file: the same members.tsv
        args[2], args[4] = TWO, os.path.join(d, "new")
        rc, out = run(*(args + ["--family", "family_a"]))
        check("the same family from the families file gives the same members.tsv",
              rc == 0 and open(os.path.join(d, "new", "members.tsv"), "rb").read() == got, out)
        rc, out = run(*args)
        check("prepare on a two-family file without --family stops and lists the families",
              rc == 2 and "holds 2 families (family_a, family_b)" in out, out)
        rc, out = plan("--families", TWO)
        check("plan on a two-family file without --family or --sample stops and lists the families",
              rc == 2 and "holds 2 families (family_a, family_b)" in out, out)
        rc, out = plan("--families", TWO, "--family", "family_c")
        check("an unknown --family is named, with the file's families",
              rc == 2 and "no family 'family_c'; its families are: family_a, family_b" in out, out)
        rc, out = plan("--families", TWO, "--sample", "partner_y", "--which")
        check("--which prints the family of the sample", rc == 0 and out.strip().splitlines()[-1] == "family_b", out)
        rc, out = plan("--families", TWO, "--sample", "nobody", "--which")
        check("--which: a sample in no family prints no family and exits 0",
              rc == 0 and "family_" not in out, out)

    head = "family\tsample\tsex\tfather\tmother\tpartner\n"
    bad = [
        ("a sample in two families", head + "fa\tkid\tmale\t.\t.\t.\nfb\tkid\tmale\t.\t.\t.\nfb\tmum\tfemale\t.\t.\t.\n",
         "line 3: sample 'kid' is in family 'fb' here and in family 'fa' on line 2"),
        ("a family name breaking the sample-name rule", head + "f/a\tkid\tmale\t.\t.\t.\n",
         "line 2: family name 'f/a' may hold only letters"),
        ("a family name equal to a sample", head + "kid\tkid\tmale\t.\t.\t.\nkid\tmum\tfemale\t.\t.\t.\n",
         "the family name 'kid' is also a person's name (line 2)"),
        ("a parent in another family", head + "fa\tkid\tmale\tdad\t.\t.\nfb\tdad\tmale\t.\t.\t.\n",
         "line 2: the father of 'kid', 'dad', is in family 'fb' (line 3)"),
        ("a partner in another family", head + "fa\tdad\tmale\t.\t.\tmum\nfb\tmum\tfemale\t.\t.\tdad\n",
         "line 2: the partner of 'dad', 'mum', is in family 'fb' (line 3)"),
        ("the column filled on some rows only", head + "fa\tkid\tmale\t.\t.\t.\n.\tmum\tfemale\t.\t.\t.\n",
         "line 3: no family name, but other rows have one (line 2: 'fa')"),
    ]
    for desc, text, msg in bad:
        with tempfile.TemporaryDirectory() as d:
            rc, out = plan("--families", write(os.path.join(d, "f.tsv"), text), "--family", "fa")
            check("plan stops on %s (exit %d), naming the file and line" % (desc, rc),
                  rc == 2 and msg in out and "f.tsv" in out, out)
    with tempfile.TemporaryDirectory() as d:
        rc, out = plan("--families", TWO, "--sample", "partner_x", "--sex", "female")
        check("a sex that contradicts the file stops plan (exit %d)" % rc,
              rc == 2 and "says 'partner_x' is male, the sex given to run-all.sh says female" in out, out)
        rc, out = plan("--families", TWO, "--sample", "partner_x", "--sex", "male")
        check("the sex the file gives passes", rc == 0, out)


def test_plan():
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "plan.tsv")
        # one member with a gVCF: everything waits, naming the missing one
        rc, log = plan("--families", TWO, "--sample", "partner_x", "--have-gvcf", "partner_x", "--out", out)
        p = by_analysis(out) if rc == 0 else {}
        check("plan writes the columns family, analysis, members, decision, reason",
              rc == 0 and open(out).readline().rstrip("\n").split("\t") == list(fr.PLAN_COLUMNS), log)
        check("one member with a gVCF: the joint call waits for partner_y",
              p.get("joint_call", [{}])[0].get("decision") == "waiting"
              and p["joint_call"][0]["reason"] == "partner_y has no gVCF yet: run run-all.sh partner_y first", str(p))
        check("and the couple screen waits too", p.get("couple_carrier", [{}])[0].get("decision") == "waiting", str(p))
        check("every row is family_b's", all(r["family"] == "family_b" for rs in p.values() for r in rs), str(p))

        # the trio: phasing runs
        rc, log = plan("--families", TWO, "--family", "family_a", "--have-gvcf", "child_a", "father_a", "mother_a",
                       "--have-clinvar", "child_a", "father_a", "mother_a", "--ibd-panel", "yes", "--out", out)
        p = by_analysis(out)
        check("a trio: the joint call runs on the three", p["joint_call"][0]["decision"] == "runs"
              and p["joint_call"][0]["members"] == "child_a,father_a,mother_a", str(p))
        check("a trio: phasing by inheritance runs for child_a",
              (p["trio_phasing"][0]["decision"], p["trio_phasing"][0]["members"]) == ("runs", "child_a"), str(p))
        check("a trio: the chromosome map runs", p["inheritance_map"][0]["decision"] == "runs", str(p))
        check("with the panel the shared segments run", p["shared_segments"][0]["decision"] == "runs", str(p))
        check("the parents, partners with ClinVar hits, get the couple screen",
              (p["couple_carrier"][0]["decision"], p["couple_carrier"][0]["members"]) == ("runs", "father_a,mother_a"),
              str(p))

        # a duo: phasing skipped, the map runs; no panel: segments skipped; a partner without ClinVar: couple skipped
        rc, log = plan("--families", TWO, "--family", "family_a", "--have-gvcf", "child_a", "mother_a",
                       "--have-clinvar", "mother_a", "--ibd-panel", "no", "--out", out)
        p = by_analysis(out)
        check("a duo: the joint call runs on two and names the one left out",
              p["joint_call"][0]["decision"] == "runs" and p["joint_call"][0]["members"] == "child_a,mother_a"
              and "left out: father_a has no gVCF yet" in p["joint_call"][0]["reason"], str(p))
        check("a duo: phasing by inheritance is skipped", p["trio_phasing"][0]["decision"] == "skipped"
              and "no child has both parents sequenced" in p["trio_phasing"][0]["reason"], str(p))
        check("a duo: the chromosome map runs", p["inheritance_map"][0]["decision"] == "runs", str(p))
        check("--ibd-panel no: the shared segments are skipped", p["shared_segments"][0]["decision"] == "skipped"
              and "setup.sh --ancestry-fine-panel" in p["shared_segments"][0]["reason"], str(p))
        check("a partner without a gVCF: the couple screen waits for him",
              p["couple_carrier"][0]["decision"] == "waiting" and "father_a has no gVCF yet" in p["couple_carrier"][0]["reason"],
              str(p))
        rc, log = plan("--families", TWO, "--family", "family_b", "--have-gvcf", "partner_x", "partner_y",
                       "--have-clinvar", "partner_x", "--out", out)
        p = by_analysis(out)
        check("a partner with no ClinVar hits: the couple screen is skipped, naming her",
              p["couple_carrier"][0]["decision"] == "skipped" and "partner_y has no ClinVar hits" in p["couple_carrier"][0]["reason"],
              str(p))
        check("partners only: no child, so no phasing and no map",
              p["trio_phasing"][0]["decision"] == "skipped" and p["inheritance_map"][0]["decision"] == "skipped", str(p))

        # prepare writes the same table as plan, for the same members
        rc, log = plan("--families", TWO, "--family", "family_a", "--have-gvcf", "child_a", "father_a", "mother_a",
                       "--have-clinvar", "father_a", "--ibd-panel", "yes", "--out", out)
        args = ["prepare", "--pedigree", TWO, "--family", "family_a", "--have-clinvar", "father_a", "--ibd-panel", "yes",
                "--pedigree-out", os.path.join(d, "ped.tsv"), "--out-dir", os.path.join(d, "prep")]
        for s in ("child_a", "father_a", "mother_a"):
            args += ["--gvcf", "%s=%s" % (s, gvcf(d, s + ".g.vcf.gz", s))]
        rc2, log2 = run(*args)
        pp = os.path.join(d, "prep", "family_plan.tsv")
        check("prepare's family_plan.tsv is plan's, byte for byte",
              rc == 0 and rc2 == 0 and os.path.exists(pp) and open(pp, "rb").read() == open(out, "rb").read(),
              log + log2)
        check("prepare --pedigree-out writes the family's own pedigree",
              rc2 == 0 and [r["sample"] for r in tsv(os.path.join(d, "ped.tsv"))] == ["child_a", "father_a", "mother_a"],
              log2)


if __name__ == "__main__":
    test_prepare()
    test_kinship()
    test_check()
    test_report()
    test_families_file()
    test_plan()
    if FAILS:
        print("%d check(s) failed" % len(FAILS))
        sys.exit(1)
    print("all checks passed")
