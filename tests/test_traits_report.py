#!/usr/bin/env python3
"""bin/traits_report.py (step 42) on synthetic bcftools output, on assets/traits.tsv and on
assets/genotype_rules.tsv; and scripts/ci/check-rule-alleles.sh on the rules.

Checks:
  1. a reference block expanded at a position (ALT '.', genotype G/G) reads as
     ref/ref with the table's hom_ref text; ref/alt and alt/alt read het and
     hom_alt, phased (C|T) as unphased, the reference allele printed first;
  2. a planted alternate allele changes the reading: the same position read
     as ref/ref and then as ref/alt gives two different readings, each the
     table's own;
  3. a position with no line is not_covered (never read as the reference); an
     allele the table does not describe and a deletion starting at the
     position are other_allele, not read; a half call that reaches the script
     (the bcftools genotype filter of step 42 drops it first, so in a run its
     position is not_covered) is other_allele too;
  4. a gVCF REF that is not the table's reference base stops the run (exit 1,
     the row named), so a gVCF on another build is not read;
  5. the APOE line, now read from the APOE rule: for every pair of rs429358 and
     rs7412 genotypes, one site not covered and an allele the rule does not
     describe, the line is byte for byte the line the hand code of main
     (6056452) wrote, by its sha256; it comes right after the second APOE row,
     once; control: one changed character in a reading changes the hash;
     report() without rules reads the repository's rules table;
  6. the rules table passes the schema check, and planted problems are each
     caught: a missing and an extra combination line, a haplotype allele that
     is neither ref nor alt, diplotypes out of haplotype order, a site that
     disagrees with its traits table row, no APOE rule, a line of an unknown
     kind;
  7. the rules on synthetic gVCF records, for each rule: reference blocks only
     (homozygous reference, with the block's GQ; TAS2R38's reference haplotype
     is AAI, which the rule does not read); APOE's double heterozygote
     (e2/e4, maybe e1/e3) without phase; TAS2R38 homozygous PAV/PAV and a
     single heterozygous site, which need no phase; the triple heterozygote
     unphased (WhatsHap ran and no read pair spans the sites, or no phasing
     ran), phased in one block (PAV/AVI, read-backed), phased in two blocks
     (unphased), phased into a haplotype the rule does not describe; a site
     not covered;
  8. --phase-sites writes only the heterozygous sites of a rule that needs
     phase and has two or more of them, as a sorted VCF; WhatsHap's phased VCF
     is read back (GT with |, PS), its unphased records left out;
  9. the command line writes the traits table (one row per table row plus the
     APOE line, the columns in order) and the rules table (one line per rule),
     and exits 1 on a traits table or a rules table that fails the schema check;
 10. scripts/ci/check-rule-alleles.sh passes on the rules and fails, naming the
     site and the mistake, on planted copies: ref and alt swapped, the
     minus-strand bases, an rsID the dbSNP table lacks, a wrong position; and
     on a rules file with no site, and on a dbSNP table of another build than
     versions.env pins.

Run: python3 tests/test_traits_report.py
"""
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import traits_report as tr  # noqa: E402

TABLE = os.path.join(REPO, "assets", "traits.tsv")
RULES = os.path.join(REPO, "assets", "genotype_rules.tsv")
SCRIPT = os.path.join(REPO, "bin", "traits_report.py")
CHECK = os.path.join(REPO, "scripts", "ci", "check-rule-alleles.sh")
FAILS = 0

# The sha256 of the APOE line (its OUT_COLUMNS joined by tabs) the hand code of
# main 6056452 wrote, per rs429358 and rs7412 genotype ("." = no line).
APOE_GOLDEN = {
    ("T/T", "C/C"): "e3d7512b649f1d23f5acbcaa08da19bd6a009d5b39ac57cb6856591c50e7cfda",
    ("T/T", "C/T"): "62e910d8ec7ac45ab15e9c9ed85c93b128ac41703d1149a332082dc0d7b5f389",
    ("T/T", "T/T"): "41b4ad51c67028bb19af1b7ab7e4f414f5c0354e4c4a7c3c19309186b429f960",
    ("T/C", "C/C"): "08c5ed9e2898f2ae29cbd49d5ac97394758fa852f9d35e3c7173b2aa8d3111ad",
    ("T/C", "C/T"): "555c9cc0e8b59b10bdc5bdc2e33d93fb4ccbee1e6ea12916b17e73ffd6fdb8d8",
    ("T/C", "T/T"): "043f51b41486f217cac3a71e705e0e29c1a1ced9745e60abde2230df2fe901db",
    ("C/C", "C/C"): "ca7f36540cf1269b84649b4ad6ace19e0abc03eab7984d495de584540bb11dcc",
    ("C/C", "C/T"): "db499d1bfebedf77ed34a9d0525030d34be6336e685758a8efd53125197e3d00",
    ("C/C", "T/T"): "16e5214f90c5f2c401e5849dffd76a08dc3fc5ba2deb324d7fb02c76a92113be",
    ("T/C", "."): "34624db16a92981c0882edb3d8ce8d291162aa78b7c5208eba2360141e6ccbd8",
    ("T/G", "C/C"): "1d0aee0b8d241561ebd7477d6f4a8be274a48157a297abecd444c6c6a0a6abac",
}
# TAS2R38's three sites in the rule's order (codons 49, 262, 296).
T1, T2, T3 = ("chr7", "141973545"), ("chr7", "141972905"), ("chr7", "141972804")


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def recs(*lines):
    """bcftools' CHROM POS REF ALT GENOTYPE [GQ] lines -> read_genotypes' records."""
    d = tempfile.mkdtemp()
    try:
        p = os.path.join(d, "g.tsv")
        write(p, "".join("\t".join(map(str, l)) + "\n" for l in lines))
        return tr.read_genotypes(p)
    finally:
        shutil.rmtree(d)


def phased_vcf(*records):
    """WhatsHap-style VCF records (chrom, pos, ref, alt, GT[:PS], GT as allele indexes) -> read_phased's dict."""
    d = tempfile.mkdtemp()
    try:
        p = os.path.join(d, "p.vcf")
        body = "".join(f"{c}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t{'GT:PS' if ':' in g else 'GT'}\t{g}\n"
                       for c, pos, ref, alt, g in records)
        write(p, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + body)
        return tr.read_phased(p)
    finally:
        shutil.rmtree(d)


def apoe_line(out):
    return next(o for o in out if o["rsid"] == "rs429358+rs7412")


def line_hash(o):
    return hashlib.sha256("\t".join(str(o[c]) for c in tr.OUT_COLUMNS).encode()).hexdigest()


def run(*args):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)


def dbsnp_checks(rules):
    """10. scripts/ci/check-rule-alleles.sh on the rules and on planted copies. Runs first,
    so a rules table that breaks the other checks still shows what dbSNP says."""
    work = tempfile.mkdtemp()
    try:
        def dbsnp(path, *extra):
            r = subprocess.run(["bash", CHECK, *extra, path], capture_output=True, text=True)
            return r.returncode, r.stdout + r.stderr

        rc, log = dbsnp(RULES)
        print("  " + log.strip().replace("\n", "\n  "))
        n = sum(len(v["sites"]) for v in rules.values())
        check(f"check-rule-alleles.sh: every rule site matches dbSNP ({n} sites)",
              rc == 0 and f"{n} rule sites checked" in log, log)
        text = open(RULES).read()
        site = "site\tTAS2R38\trs1726866\tchr7\t141972905\t"

        def planted(desc, new, want):
            f = os.path.join(work, "planted.tsv")
            write(f, text.replace(site + "G\tA", site + new))
            rc, log = dbsnp(f)
            check(f"check-rule-alleles.sh control: {desc} fails, naming it", rc == 1 and want in log, (rc, log))

        planted("ref and alt swapped", "A\tG", "rs1726866: ref and alt are swapped")
        planted("the minus-strand bases", "C\tT", "rs1726866: the alleles are the minus-strand bases")
        f = os.path.join(work, "planted.tsv")
        write(f, text.replace("\trs1726866\tchr7\t141972905\t", "\trs1726867\tchr7\t141972905\t"))
        rc, log = dbsnp(f)
        check("check-rule-alleles.sh control: an rsID the dbSNP table lacks fails",
              rc == 1 and "rs1726867: not in the pinned dbSNP table" in log, (rc, log))
        write(f, text.replace("\trs1726866\tchr7\t141972905\t", "\trs1726866\tchr7\t141972906\t"))
        rc, log = dbsnp(f)
        check("check-rule-alleles.sh control: a wrong position fails",
              rc == 1 and "dbSNP at chr7:141972905" in log, (rc, log))
        write(f, "# no sites\n")
        rc, log = dbsnp(f)
        check("check-rule-alleles.sh control: a rules file with no site fails", rc == 1 and "no site line" in log,
              (rc, log))
        root = os.path.join(work, "root")
        os.makedirs(os.path.join(root, "scripts", "ci"))
        shutil.copy(os.path.join(REPO, "versions.env"), root)
        copy = os.path.join(root, "scripts", "ci", "check-rule-alleles.sh")
        write(copy, open(CHECK).read().replace("\n#: # build ", "\n#: # build x"))
        r = subprocess.run(["bash", copy, RULES], capture_output=True, text=True)
        check("check-rule-alleles.sh control: a dbSNP table of another build than versions.env pins fails",
              r.returncode == 1 and "versions.env pins" in r.stderr, (r.returncode, r.stdout, r.stderr))
    finally:
        shutil.rmtree(work)


def main():
    rules = tr.read_rules(RULES)
    dbsnp_checks(rules)
    rows = tr.read_table(TABLE)
    by_id = {r["rsid"]: r for r in rows}
    lct = by_id["rs4988235"]   # chr2:135851076 G>A
    c, p, ref, alt = lct["chrom"], lct["pos"], lct["ref"], lct["alt"]
    at = lambda *ls: recs(*ls).get((c, p), [])  # noqa: E731

    # 1. genotypes and their readings
    g = tr.call(lct, at((c, p, ref, ".", f"{ref}/{ref}")))
    check("expanded reference block: ref/ref, hom_ref reading", g == (f"{ref}/{ref}", "called", lct["hom_ref"]), g)
    g = tr.call(lct, at((c, p, ref, alt, f"{alt}/{ref}")))
    check("ref/alt: het reading, reference allele first", g == (f"{ref}/{alt}", "called", lct["het"]), g)
    g = tr.call(lct, at((c, p, ref, alt, f"{ref}|{alt}")))
    check("phased ref|alt reads as ref/alt", g == (f"{ref}/{alt}", "called", lct["het"]), g)
    g = tr.call(lct, at((c, p, ref, alt, f"{alt}/{alt}")))
    check("alt/alt: hom_alt reading", g == (f"{alt}/{alt}", "called", lct["hom_alt"]), g)

    # 2. a planted alternate allele changes the reading
    before = tr.call(lct, at((c, p, ref, ".", f"{ref}/{ref}")))[2]
    after = tr.call(lct, at((c, p, ref, alt, f"{ref}/{alt}")))[2]
    check("planted alt: the reading changes from hom_ref to het", before == lct["hom_ref"] and after == lct["het"]
          and before != after, (before, after))

    # 3. not covered and other alleles
    check("no line: not_covered", tr.call(lct, [])[1] == "not_covered")
    other = next(b for b in "ACGT" if b not in (ref, alt))
    g = tr.call(lct, at((c, p, ref, other, f"{ref}/{other}")))
    check("an allele the table does not describe: other_allele, not read", g[1] == "other_allele"
          and g[2] == tr.OTHER_ALLELE, g)
    g = tr.call(lct, at((c, p, ref, alt, f"./{alt}")))
    check("a half call: other_allele", g[1] == "other_allele", g)
    g = tr.call(lct, at((c, p, ref + "T", ref, f"{ref}T/{ref}")))
    check("a deletion starting at the position: other_allele", g[1] == "other_allele", g)
    g = tr.call(lct, at((c, p, ref + "T", ref, f"{ref}T/{ref}"), (c, p, ref, alt, f"{ref}/{alt}")))
    check("a SNV record beside an indel at the position is the one read", g == (f"{ref}/{alt}", "called", lct["het"]), g)

    # 4. a REF that is not the GRCh38 base stops
    try:
        tr.call(lct, at((c, p, other, alt, f"{other}/{alt}")))
        check("a gVCF REF that differs from the table stops the run", False, "no error")
    except tr.TableError as e:
        check("a gVCF REF that differs from the table stops the run, naming the row", "rs4988235" in str(e), e)

    # 5. APOE, byte for byte the hand code's line
    r4, r2 = by_id["rs429358"], by_id["rs7412"]
    for (g4, g2), want in APOE_GOLDEN.items():
        lines = [(r4["chrom"], r4["pos"], "T", g4.split("/")[1], g4)]
        if g2 != ".":
            lines.append((r2["chrom"], r2["pos"], "C", "T", g2))
        out = tr.report(rows, recs(*lines), rules)
        line = apoe_line(out)
        check(f"APOE line byte for byte as before: rs429358 {g4}, rs7412 {g2} -> {line['genotype']} {line['status']}",
              line_hash(line) == want, line)
    out = tr.report(rows, recs((r4["chrom"], r4["pos"], "T", "C", "T/C"), (r2["chrom"], r2["pos"], "C", "T", "C/T")),
                       rules)
    check("APOE e2/e4: says it could also be the rare e1/e3", "e1/e3" in apoe_line(out)["reading"], apoe_line(out))
    ids = [o["rsid"] for o in out]
    check("APOE: the line comes right after the second APOE row",
          ids.index("rs429358+rs7412") == max(ids.index("rs429358"), ids.index("rs7412")) + 1, ids)
    check("APOE: one line only", ids.count("rs429358+rs7412") == 1, ids)
    changed = dict(rules["APOE"], combinations=[dict(x, reading=x["reading"].replace("e3/e3.", "e3/e3!"))
                                                for x in rules["APOE"]["combinations"]])
    out = tr.report(rows, recs((r4["chrom"], r4["pos"], "T", "C", "T/C"), (r2["chrom"], r2["pos"], "C", "T", "C/C")),
                       dict(rules, APOE=changed))
    check("control: one changed character in an APOE reading changes the hash",
          line_hash(apoe_line(out)) != APOE_GOLDEN[("T/C", "C/C")], apoe_line(out))
    e4 = apoe_line(tr.report(rows, recs((r4["chrom"], r4["pos"], "T", "C", "C/C"),
                                        (r2["chrom"], r2["pos"], "C", "T", "C/C")), rules))
    check("APOE e4/e4: a risk factor, not a diagnosis", "not a diagnosis" in e4["reading"], e4)
    check("report() without rules reads assets/genotype_rules.tsv (other tests call it so)",
          line_hash(apoe_line(tr.report(rows, {}))) == line_hash(apoe_line(tr.report(rows, {}, rules))), rows[:0])

    # 6. the rules table's schema
    errs = tr.check_rules(rules, rows)
    check(f"assets/genotype_rules.tsv passes the schema check ({len(rules)} rules)", not errs, errs)
    check("the APOE and TAS2R38 rules are there", {"APOE", "TAS2R38"} <= set(rules), list(rules))

    def broken(desc, mutate, want, table=rows):
        rs = {k: dict(v, sites=[dict(s) for s in v["sites"]], haplotypes=[dict(h) for h in v["haplotypes"]],
                      combinations=[dict(x) for x in v["combinations"]]) for k, v in rules.items()}
        mutate(rs)
        e = tr.check_rules(rs, table)
        check(f"control: {desc} is caught", any(want in x for x in e), e)

    broken("a missing combination line", lambda rs: rs["TAS2R38"]["combinations"].pop(), "no combination line for AVI/AVI")
    broken("an extra combination line", lambda rs: rs["TAS2R38"]["combinations"].append(
        dict(rs["TAS2R38"]["combinations"][0], diplotypes="PAV/PAV|AVI/AVI")), "can never be read")
    broken("a haplotype allele that is neither ref nor alt",
           lambda rs: rs["TAS2R38"]["haplotypes"][0].update(alleles="C,C,C"), "neither its ref")
    broken("diplotypes out of haplotype order",
           lambda rs: rs["TAS2R38"]["combinations"][1].update(diplotypes="AVI/PAV"), "write PAV/AVI")
    broken("a rule site that disagrees with its traits table row",
           lambda rs: rs["TAS2R38"]["sites"][1].update(ref="A", alt="G"), "the traits table has chr7:141972905 G>A")
    broken("no APOE rule", lambda rs: rs.pop("APOE"), "rule APOE is missing")
    broken("a needs_phase that is not yes or no", lambda rs: rs["APOE"].update(needs_phase="maybe"), "needs_phase")
    work = tempfile.mkdtemp()
    try:
        bad = os.path.join(work, "r.tsv")
        write(bad, open(RULES).read() + "allele\tAPOE\tx\n")
        try:
            tr.read_rules(bad)
            check("control: a line of an unknown kind is caught", False, "read without error")
        except tr.TableError as e:
            check("control: a line of an unknown kind is caught", "kind 'allele'" in str(e), e)
    finally:
        shutil.rmtree(work)

    # 7. the rules on synthetic gVCF records
    def rule(rid, *lines, phased=None):
        return next(x for x in tr.rule_lines(rules, recs(*lines), phased) if x["rule"] == rid)

    a4, a2 = (r4["chrom"], r4["pos"]), (r2["chrom"], r2["pos"])
    l = rule("APOE", (*a4, "T", ".", "T/T", "37"), (*a2, "C", ".", "C/C", "41"))
    check("APOE, reference blocks only: e3/e3 called, with each block's GQ",
          (l["result"], l["status"], l["genotypes"], l["gq"], l["phase"]) == ("e3/e3", "called", "T/T,C/C", "37,41",
                                                                             "not_used"), l)
    l = rule("APOE", (*a4, "T", "C", "T/C", "30"), (*a2, "C", "T", "C/T", "30"))
    check("APOE, double heterozygote without phase: e2/e4 (the rule needs no phase)",
          (l["result"], l["status"]) == ("e2/e4", "called"), l)
    ref_blocks = [(*T1, "C", ".", "C/C", "50"), (*T2, "G", ".", "G/G", "50"), (*T3, "T", ".", "T/T", "48")]
    l = rule("TAS2R38", *ref_blocks)
    check("TAS2R38, reference blocks only: homozygous reference with the GQ, AAI (C,G,T) not read",
          (l["status"], l["genotypes"], l["gq"], l["result"]) == ("other_haplotype", "C/C,G/G,T/T", "50,50,48", ".")
          and "C,G,T" in l["reading"], l)
    pav = [(*T1, "C", "G", "G/G", "40"), (*T2, "G", ".", "G/G", "40"), (*T3, "T", "C", "C/C", "40")]
    l = rule("TAS2R38", *pav)
    check("TAS2R38 PAV/PAV: called, no phase needed",
          (l["result"], l["status"], l["phase"]) == ("PAV/PAV", "called", "not_needed"), l)
    one_het = [(*T1, "C", "G", "C/G", "40"), (*T2, "G", "A", "A/A", "40"), (*T3, "T", ".", "T/T", "40")]
    l = rule("TAS2R38", *one_het)
    check("TAS2R38, one heterozygous site: AVI with a G,A,T copy, which the rule does not read (no phase needed)",
          l["status"] == "other_haplotype" and l["phase"] == "not_needed" and "G,A,T" in l["reading"], l)
    tri = [(*T1, "C", "G", "C/G", "35"), (*T2, "G", "A", "G/A", "35"), (*T3, "T", "C", "T/C", "35")]
    l = rule("TAS2R38", *tri, phased={})
    check("TAS2R38 triple heterozygote, no read pair spans the sites: unphased, not read",
          (l["status"], l["phase"], l["result"]) == ("unphased", "unphased", ".")
          and "no read or read pair" in l["reading"], l)
    l = rule("TAS2R38", *tri)
    check("TAS2R38 triple heterozygote, no phasing ran: unphased, says why",
          l["status"] == "unphased" and "no phasing ran" in l["reading"], l)
    one_block = phased_vcf(("chr7", T3[1], "T", "C", "1|0:141972804"), ("chr7", T2[1], "G", "A", "0|1:141972804"),
                           ("chr7", T1[1], "C", "G", "1|0:141972804"))
    l = rule("TAS2R38", *tri, phased=one_block)
    check("TAS2R38 triple heterozygote phased in one block: PAV/AVI, read-backed",
          (l["result"], l["status"], l["phase"]) == ("PAV/AVI", "called", "read_backed"), l)
    two_blocks = phased_vcf(("chr7", T3[1], "T", "C", "1|0:141972804"), ("chr7", T2[1], "G", "A", "0|1:141972804"),
                            ("chr7", T1[1], "C", "G", "1|0:141973545"))
    l = rule("TAS2R38", *tri, phased=two_blocks)
    check("TAS2R38 triple heterozygote phased in two blocks: unphased", l["status"] == "unphased", l)
    odd = phased_vcf(("chr7", T3[1], "T", "C", "0|1:141972804"), ("chr7", T2[1], "G", "A", "0|1:141972804"),
                     ("chr7", T1[1], "C", "G", "1|0:141972804"))
    l = rule("TAS2R38", *tri, phased=odd)
    check("TAS2R38 phased into G,G,T and C,A,C: a haplotype the rule does not describe",
          l["status"] == "other_haplotype" and "G,G,T and C,A,C" in l["reading"], l)
    l = rule("TAS2R38", *tri[:2])
    check("TAS2R38, a site not covered: not_covered, names the three sites",
          l["status"] == "not_covered" and "every one of rs713598, rs1726866 and rs10246939" in l["reading"], l)

    # 8. the sites to phase, and WhatsHap's VCF read back
    s = tr.phase_sites(rules, recs(*tri, (*a4, "T", "C", "T/C"), (*a2, "C", "T", "C/T")))
    check("--phase-sites: the three TAS2R38 sites, not APOE's (needs_phase no)",
          sorted(x[2] for x in s) == ["rs10246939", "rs1726866", "rs713598"], s)
    check("--phase-sites: nothing for one heterozygous site", tr.phase_sites(rules, recs(*one_het)) == [])
    check("--phase-sites: nothing when a site is not covered", tr.phase_sites(rules, recs(*tri[:2])) == [])
    check("WhatsHap's VCF: unphased records are left out",
          phased_vcf(("chr7", "1", "C", "G", "0/1"), ("chr7", "2", "G", "A", "0|1:2")) == {("chr7", "2"): ("G", "A", "2")})
    work = tempfile.mkdtemp()
    try:
        g = os.path.join(work, "g.tsv")
        write(g, "".join("\t".join(x) + "\n" for x in tri))
        vcf = os.path.join(work, "in.vcf")
        res = run("--rules", RULES, "--genotypes", g, "--sample", "S", "--phase-sites", vcf)
        body = [x.split("\t") for x in open(vcf).read().splitlines() if not x.startswith("#")] if res.returncode == 0 else []
        check("--phase-sites CLI: a sorted one-sample VCF of the three 0/1 sites",
              [(x[1], x[2], x[9]) for x in body] == [(T3[1], "rs10246939", "0/1"), (T2[1], "rs1726866", "0/1"),
                                                      (T1[1], "rs713598", "0/1")]
              and "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS" in open(vcf).read(), (res.stderr, body))

        # 9. the command line
        write(g, f"{c}\t{p}\t{ref}\t{alt}\t{ref}/{alt}\t30\n" + "".join("\t".join(x) + "\n" for x in tri))
        out, rout = os.path.join(work, "S_traits.tsv"), os.path.join(work, "S_genotype_rules.tsv")
        ph = os.path.join(work, "phased.vcf")
        write(ph, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
              f"chr7\t{T3[1]}\t.\tT\tC\t.\tPASS\t.\tGT:PS\t1|0:9\nchr7\t{T2[1]}\t.\tG\tA\t.\tPASS\t.\tGT:PS\t0|1:9\n"
              f"chr7\t{T1[1]}\t.\tC\tG\t.\tPASS\t.\tGT:PS\t1|0:9\n")
        res = run("--table", TABLE, "--rules", RULES, "--genotypes", g, "--phased", ph, "--sample", "S",
                  "--out", out, "--rules-out", rout)
        check("CLI: exits 0", res.returncode == 0, res.stderr)
        lines = open(out).read().splitlines() if res.returncode == 0 else [""]
        check("CLI: the header is the output columns", lines[0].split("\t") == tr.OUT_COLUMNS, lines[0])
        check("CLI: one row per table row plus the APOE line", len(lines) - 1 == len(rows) + 1, len(lines))
        row = next((l.split("\t") for l in lines[1:] if l.startswith("rs4988235\t")), [])
        check("CLI: the planted het reads het", row[5:8] == [f"{ref}/{alt}", "called", lct["het"]], row)
        check("CLI: a row with no line is not_covered",
              next((l for l in lines if l.startswith("rs671\t")), "").split("\t")[6:7] == ["not_covered"])
        rl = open(rout).read().splitlines() if res.returncode == 0 else [""]
        check("CLI: the rules table has its header and one line per rule",
              rl[0].split("\t") == tr.RULE_OUT_COLUMNS and len(rl) - 1 == len(rules), rl[:1])
        t = next((x.split("\t") for x in rl if x.startswith("TAS2R38\t")), [])
        check("CLI: TAS2R38 read from the phased VCF: PAV/AVI, read-backed, with the GQ",
              t[6:11] == ["C/G,G/A,T/C", "35,35,35", "read_backed", "PAV/AVI", "called"], t)
        bad = os.path.join(work, "bad.tsv")
        write(bad, open(TABLE).read().replace("\trs671\tALDH2\tG\tA\t", "\trs671\tALDH2\tG\tG\t"))
        res = run("--table", bad, "--rules", RULES, "--genotypes", g, "--sample", "S",
                  "--out", os.path.join(work, "x.tsv"), "--rules-out", os.path.join(work, "y.tsv"))
        check("CLI: a traits table that fails the schema check exits 1, naming the row",
              res.returncode == 1 and "rs671" in res.stderr and not os.path.exists(os.path.join(work, "x.tsv")),
              (res.returncode, res.stderr))
        badr = os.path.join(work, "badr.tsv")
        write(badr, open(RULES).read().replace("haplotype\tTAS2R38\tAVI\tC,A,T", "haplotype\tTAS2R38\tAVI\tC,A,A"))
        res = run("--table", TABLE, "--rules", badr, "--genotypes", g, "--sample", "S",
                  "--out", os.path.join(work, "x.tsv"), "--rules-out", os.path.join(work, "y.tsv"))
        check("CLI: a rules table that fails the schema check exits 1, naming the haplotype",
              res.returncode == 1 and "AVI" in res.stderr and not os.path.exists(os.path.join(work, "y.tsv")),
              (res.returncode, res.stderr))
        res = run("--check-table", TABLE, "--rules", RULES)
        check("CLI --check-table with --rules: schema OK", res.returncode == 0 and "schema OK" in res.stdout,
              (res.returncode, res.stderr))

    finally:
        shutil.rmtree(work)

    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
