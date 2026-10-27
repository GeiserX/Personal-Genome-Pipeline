#!/usr/bin/env python3
"""bin/ancient_report.py (step 44) on synthetic files.

Checks:
  1. dates: years before 1950 become BCE and CE the way the AADR writes
     them (4250 BP is the middle of its 2600-2000 BCE; 1950 BP is 1 BCE,
     1949 BP is 1 CE), and a population's span is written oldest first, with
     one era when both ends share it;
  2. ranking: on a synthetic .evec the nearest ancient group is the one whose
     mean position is closest on the weighted PCs, a group far on PC1 but
     with a small eigenvalue on PC3 is ranked by the weighted distance, and
     each row carries its dates, country and member count;
  3. axes: a run whose PC2 has the opposite sign is turned back before the
     person is placed, and a run whose present-day coordinates differ stops
     with the reason instead of placing the person on other axes;
  4. liftover: a site inside a block maps with the block's offset, a site on
     a minus-strand chain maps to the mirrored position with its alleles
     complemented, a site in a gap or covered by two chains is dropped;
  5. genotypes: a bcftools genotype gives the copies of the .snp file's first
     allele, phased or not; a half call, an indel and an allele that is
     neither of the two are missing; person-geno appends one column, 9 where
     the person has no genotype;
  6. select: the present-day reference people come from MODERN_GROUPS, the
     ancient ones need a date, a passing assessment, enough sites and a place
     in West Eurasia, one genetic id per individual;
  7. the report command writes the table, the summary, the text and an SVG
     with the person marked;
  8. bin/render_report.py shows them: an HTML card with the table and the
     plot, and a section of the text report.

Run: python3 tests/test_ancient_report.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import ancient_report as ar  # noqa: E402

SCRIPT = os.path.join(REPO, "bin", "ancient_report.py")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


# --- 1. dates ---------------------------------------------------------------------
check("4250 BP is 2300 BCE, the AADR's midpoint of 2600-2000 BCE", ar.year_label(4250) == "2300 BCE", ar.year_label(4250))
check("8025 BP is 6075 BCE", ar.year_label(8025) == "6075 BCE", ar.year_label(8025))
check("1950 BP is 1 BCE", ar.year_label(1950) == "1 BCE", ar.year_label(1950))
check("1949 BP is 1 CE", ar.year_label(1949) == "1 CE", ar.year_label(1949))
check("450 BP is 1500 CE", ar.year_label(450) == "1500 CE", ar.year_label(450))
check("a span in one era is written once", ar.range_label(7250, 6850) == "5300-4900 BCE", ar.range_label(7250, 6850))
check("a span across the eras names both", ar.range_label(2250, 1750) == "300 BCE-200 CE", ar.range_label(2250, 1750))
check("a span of one date is that date", ar.range_label(3000, 3000) == "1050 BCE", ar.range_label(3000, 3000))

# --- 2. ranking -------------------------------------------------------------------
# Present-day people M0..M19 build the axes; three ancient groups.
EVALS = [9.0, 4.0, 0.01, 0.01]


def evec(rows, evals=EVALS):
    out = ["#eigvals: " + " ".join(f"{v:.3f}" for v in evals)]
    for i, (label, pcs) in rows.items():
        out.append(f"{i} " + " ".join(f"{v:.6f}" for v in pcs) + f" {label}")
    return "\n".join(out) + "\n"


MODERNS = {f"M{i}": ("French" if i % 2 else "Sardinian",
                     [i * 0.01 - 0.1, (i % 5) * 0.02 - 0.04, (i % 3) * 0.01 - 0.01, (i % 4) * 0.01 - 0.015])
           for i in range(20)}
PANEL = dict(MODERNS)
PANEL.update({
    "A1": ("Ancient", [0.10, 0.10, 0.0, 0.0]), "A2": ("Ancient", [0.12, 0.10, 0.0, 0.0]),   # Steppe
    "B1": ("Ancient", [-0.20, 0.00, 0.0, 0.0]),                                              # Farmer
    "C1": ("Ancient", [0.00, 0.00, 0.9, 0.0]),   # close on PC1/PC2, far on PC3 (eigenvalue 0.01)
})
ANCIENTS = {
    "A1": {"id": "A1", "group": "Steppe_EBA", "date_bp": "4800", "country": "Russia", "locality": "Samara"},
    "A2": {"id": "A2", "group": "Steppe_EBA", "date_bp": "4600", "country": "Russia", "locality": "Samara"},
    "B1": {"id": "B1", "group": "Anatolia_N", "date_bp": "8200", "country": "Turkey", "locality": "Barcin"},
    "C1": {"id": "C1", "group": "Odd_C", "date_bp": "3000", "country": "Italy", "locality": "Somewhere"},
}
person = [0.08, 0.09, 0.0, 0.0]
rows, n = ar.nearest_populations(person, PANEL, ANCIENTS, EVALS, 4)
check("three ancient groups ranked", n == 3 and len(rows) == 3, (n, rows))
order = [r["group"] for r in rows]
check("the nearest group is first (Steppe_EBA)", order[0] == "Steppe_EBA", order)
check("a group far only on a PC with a tiny eigenvalue is second, not last (Odd_C)", order[1] == "Odd_C", order)
check("the farthest group is last (Anatolia_N)", order[-1] == "Anatolia_N", order)
st = rows[0]
check("the row has the member count, dates and place", st["n"] == 2 and st["dates"] == "2850-2650 BCE"
      and st["country"] == "Russia" and st["locality"] == "Samara", st)
check("ranks are 1..n", [r["rank"] for r in rows] == [1, 2, 3], rows)
unweighted = ar.nearest_populations(person, PANEL, ANCIENTS, [1.0, 1.0, 1.0, 1.0], 4)[0]
check("control: without the eigenvalue weights Odd_C falls behind Anatolia_N",
      [r["group"] for r in unweighted] == ["Steppe_EBA", "Anatolia_N", "Odd_C"], [r["group"] for r in unweighted])

# --- 3. axes ----------------------------------------------------------------------
run = {k: (lab, [pcs[0], -pcs[1], pcs[2], pcs[3]]) for k, (lab, pcs) in MODERNS.items()}
run["S"] = ("Person", [0.08, -0.09, 0.0, 0.0])
signs = ar.axis_signs(PANEL, run, 2)
check("a flipped PC2 is found and turned back", signs == [1.0, -1.0], signs)
bad = {k: (lab, [pcs[0] * (1 if int(k[1:]) % 3 else -1)] + pcs[1:]) for k, (lab, pcs) in MODERNS.items()}
try:
    ar.axis_signs(PANEL, bad, 2)
    check("other axes stop the report", False, "no error")
except ar.Fail as e:
    check("other axes stop the report, naming the PC", "PC1" in str(e), str(e))

# --- 4. liftover --------------------------------------------------------------------
tmp = tempfile.mkdtemp()
try:
    chain = write(os.path.join(tmp, "test.chain"), "\n".join([
        # chr1 source 1000-1100 -> target 2000-2100 (+), gap of 50, then 1150-1200 -> 2150-2200
        "chain 100 chr1 100000 + 1000 1200 chr1 100000 + 2000 2200 1",
        "100 50 50", "50", "",
        # chr2 source 500-600 -> target on the minus strand of a 10,000 bp chr2, reverse start 300
        "chain 90 chr2 50000 + 500 600 chr2 10000 - 300 400 2",
        "100", "",
        # chr3: two chains over 100-200
        "chain 80 chr3 50000 + 100 200 chr3 50000 + 100 200 3", "100", "",
        "chain 70 chr3 50000 + 150 250 chr3 50000 + 5150 5250 4", "100", ""]) + "\n")
    ix = ar.index_blocks(ar.read_chain(chain))
    check("a site in the first block moves by its offset", ar.lift_point(ix, "chr1", 1051) == ("chr1", 2051, "+"),
          ar.lift_point(ix, "chr1", 1051))
    check("a site in the second block moves by its offset", ar.lift_point(ix, "chr1", 1160) == ("chr1", 2160, "+"),
          ar.lift_point(ix, "chr1", 1160))
    check("a site in the gap is dropped", ar.lift_point(ix, "chr1", 1120) is None, ar.lift_point(ix, "chr1", 1120))
    # 0-based 510 -> reverse 310 -> forward 10000 - 1 - 310 = 9689 (0-based) = 9690 1-based
    check("a minus-strand site is mirrored", ar.lift_point(ix, "chr2", 511) == ("chr2", 9690, "-"),
          ar.lift_point(ix, "chr2", 511))
    check("a site two chains cover is dropped", ar.lift_point(ix, "chr3", 171) is None, ar.lift_point(ix, "chr3", 171))
    check("a site one of them covers maps", ar.lift_point(ix, "chr3", 121) == ("chr3", 121, "+"),
          ar.lift_point(ix, "chr3", 121))
    # The lift command against a small FASTA.
    seq1 = ["N"] * 3000
    seq1[2050] = "G"   # 1-based 2051
    seq1[2159] = "T"   # 1-based 2160
    seq2 = ["N"] * 10000
    seq2[9689] = "C"   # 1-based 9690: the complement of the source G
    fa = os.path.join(tmp, "ref.fa")
    with open(fa, "w") as f, open(fa + ".fai", "w") as fai:
        off = 0
        for name, seq in (("chr1", seq1), ("chr2", seq2)):
            head = f">{name}\n"
            f.write(head)
            off += len(head)
            s = "".join(seq)
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")
            fai.write(f"{name}\t{len(s)}\t{off}\t60\t61\n")
            off += len(s) + (len(s) + 59) // 60
    snp = write(os.path.join(tmp, "src.snp"), "\n".join([
        "rs1 1 0.0 1051 G A",      # ref G at target: kept
        "rs2 1 0.0 1160 C T",      # second allele is the target base: kept
        "rs3 1 0.0 1120 A G",      # gap: dropped
        "rs4 2 0.0 511 G T",       # minus strand: alleles C/A, target C: kept
        "rs5 1 0.0 1051 A C",      # neither allele is G: dropped
        "rs6 23 0.0 1051 G A",     # X: dropped
        "rs7 1 0.0 1051 G X",      # monomorphic marker: dropped
    ]) + "\n")
    out = os.path.join(tmp, "sites.tsv")
    rc = ar.main(["lift", "--snp", snp, "--chain", chain, "--fasta", fa, "--out", out])
    got = open(out).read().splitlines() if os.path.exists(out) else []
    check("lift keeps rs1, rs2 and the complemented rs4, in GRCh38 order",
          rc == 0 and got == ["snp\tchrom\tpos\tallele1\tallele2", "rs1\tchr1\t2051\tG\tA",
                              "rs2\tchr1\t2160\tC\tT", "rs4\tchr2\t9690\tC\tA"], (rc, got))

    # --- 5. genotypes -----------------------------------------------------------------
    check("A/A is two copies of allele1 A", ar.dosage("A/A", "A", "G") == 2)
    check("A|G is one copy", ar.dosage("A|G", "A", "G") == 1)
    check("G/G is none", ar.dosage("G/G", "A", "G") == 0)
    check("a half call is missing", ar.dosage("A", "A", "G") is None)
    check("an allele outside the pair is missing", ar.dosage("A/T", "A", "G") is None)
    check("an indel is missing", ar.dosage("A/AT", "A", "G") is None)
    msnp = write(os.path.join(tmp, "moderns.snp"),
                 "rs1 1 0.0 1051 G A\nrs2 1 0.0 1160 C T\nrsX 23 0.0 5 A G\nrs4 2 0.0 511 G T\n")
    mgeno = write(os.path.join(tmp, "moderns.geno"), "012\n220\n999\n101\n")
    gts = write(os.path.join(tmp, "gts.tsv"), "chr1\t2051\tG\tA\tG/A\nchr1\t2160\tT\t.\tT/T\nchr2\t9690\tC\tG\tC/G\n")
    merged = os.path.join(tmp, "merged.geno")
    stats = os.path.join(tmp, "stats.tsv")
    rc = ar.main(["person-geno", "--sites", out, "--genotypes", gts, "--snp", msnp, "--moderns-geno", mgeno,
                  "--out", merged, "--stats", stats])
    got = open(merged).read().splitlines() if os.path.exists(merged) else []
    check("person-geno appends one column: 1, 0, 9 for X, 9 for an allele outside the pair",
          rc == 0 and got == ["0121", "2200", "9999", "1019"], (rc, got))
    st = dict(l.split("\t") for l in open(stats).read().splitlines()) if os.path.exists(stats) else {}
    check("its counts", st == {"panel_sites": "4", "sites_genotyped": "3", "other_allele": "1", "sites_used": "2"}, st)
    short = write(os.path.join(tmp, "short.geno"), "012\n220\n")
    rc = ar.main(["person-geno", "--sites", out, "--genotypes", gts, "--snp", msnp, "--moderns-geno", short,
                  "--out", merged, "--stats", stats])
    check("a genotype file shorter than the .snp file stops", rc == 1, rc)

    # --- 6. select ------------------------------------------------------------------------
    head = ["Genetic ID (long text)", "Persistent Genetic ID", "Individual ID", '"Date mean in BP in years before 1950"',
            "Date standard deviation in BP [x]", "Full Date One of two formats", "Group ID", "Locality",
            "Political Entity", "Latitude", "Longitude", "Data type",
            "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)",
            "SNPs hit on autosomal targets (Computed using easystats on HO snpset)", "ASSESSMENT",
            "ASSESSMENT WARNINGS (long)"]

    def arow(gid, ind, bp, group, lat, lon, dtype, snps, ass):
        return [gid, gid, ind, str(bp), "50", "x", group, "Loc", "Land", str(lat), str(lon), dtype, str(snps), "0", ass, ""]

    rows_ = [arow(f"F{i}.DG", f"F{i}", 0, "French", 46, 2, "Shotgun.diploid", 1100000, "PROVISIONAL_PASS")
             for i in range(60)]
    rows_ += [
        arow("F0.SG", "F0", 0, "French", 46, 2, "Shotgun.diploid", 900000, "PASS"),      # same individual, fewer sites
        arow("Y1.DG", "Y1", 0, "Yakut", 62, 129, "Shotgun.diploid", 1100000, "PASS"),    # not a reference group
        arow("H1.HO", "H1", 0, "French", 46, 2, "HumanOrigins", 600000, "PASS"),         # not shotgun diploid
        arow("A1.AG", "A1", 4800, "Russia_Steppe", 53, 50, "1240k", 500000, "PASS"),
        arow("A2.AG", "A2", 4800, "Russia_Steppe", 53, 50, "1240k", 50000, "PASS"),      # too few sites
        arow("A3.AG", "A3", 4800, "Russia_Steppe", 53, 50, "1240k", 500000, "Questionable"),
        arow("A4.AG", "A4", 4800, "China_N", 40, 116, "1240k", 500000, "PASS"),          # outside West Eurasia
        arow("A5.AG", "A5", 4800, "Unknown", "..", "..", "1240k", 500000, "PASS"),       # no place
    ]
    anno = write(os.path.join(tmp, "x.anno"), "\n".join("\t".join(r) for r in [head] + rows_) + "\n")
    ind = write(os.path.join(tmp, "x.ind"), "".join(f"{r[0]} U {r[6]}\n" for r in rows_))
    oind, plist, anc = (os.path.join(tmp, n) for n in ("panel.ind", "poplist.txt", "ancients.tsv"))
    rc = ar.main(["select", "--anno", anno, "--ind", ind, "--out-ind", oind, "--poplist", plist, "--ancients", anc])
    labels = dict(l.split("\t")[0::2] for l in open(oind).read().splitlines()) if rc == 0 else {}
    check("select: the 60 French shotgun genomes are the reference, one id per person (F0.DG)",
          rc == 0 and sum(1 for v in labels.values() if v == "French") == 60 and labels.get("F0.SG") == "Ignore", labels)
    check("select: Yakut and a Human Origins genome are ignored",
          labels.get("Y1.DG") == "Ignore" and labels.get("H1.HO") == "Ignore", labels)
    check("select: only A1 is projected", [k for k, v in labels.items() if v == "Ancient"] == ["A1.AG"], labels)
    check("select: the population list", open(plist).read() == "French\n" if rc == 0 else False)
    check("select: ancients.tsv holds A1 with its group",
          rc == 0 and open(anc).read().splitlines()[1].startswith("A1.AG\tRussia_Steppe\t4800"), open(anc).read() if rc == 0 else "")

    # --- 7. the report command ----------------------------------------------------------
    pe = write(os.path.join(tmp, "panel.evec"), evec(PANEL))
    re_ = write(os.path.join(tmp, "run.evec"), evec(run))
    at = os.path.join(tmp, "anc.tsv")
    with open(at, "w") as f:
        f.write("id\tgroup\tdate_bp\tdate_sd\tfull_date\tlocality\tcountry\tlat\tlon\tsnps\n")
        for r in ANCIENTS.values():
            f.write(f"{r['id']}\t{r['group']}\t{r['date_bp']}\t50\tx\t{r['locality']}\t{r['country']}\t50\t40\t500000\n")
    outs = {k: os.path.join(tmp, f"S_{k}") for k in ("ancient.tsv", "info.tsv", "ancient.txt", "pca.svg")}
    p = subprocess.run([sys.executable, SCRIPT, "report", "--panel-evec", pe, "--run-evec", re_, "--ancients", at,
                        "--sample", "S", "--release", "test", "--stats", stats,
                        "--out-tsv", outs["ancient.tsv"], "--out-info", outs["info.tsv"],
                        "--out-txt", outs["ancient.txt"], "--out-svg", outs["pca.svg"]], capture_output=True, text=True)
    check("report exits 0", p.returncode == 0, p.stderr)
    tsv = open(outs["ancient.tsv"]).read().splitlines() if p.returncode == 0 else []
    check("report: the table's first row is Steppe_EBA, from the flipped run turned back",
          len(tsv) == 4 and tsv[1].split("\t")[1] == "Steppe_EBA", tsv)
    info = dict(l.split("\t") for l in open(outs["info.tsv"]).read().splitlines()[1:]) if p.returncode == 0 else {}
    check("report: the person's PC2 in the panel's sign", info.get("pc2") == "0.09000", info)
    check("report: the sites used", info.get("sites_used") == "2", info)
    svg = open(outs["pca.svg"]).read() if p.returncode == 0 else ""
    check("report: the SVG is an SVG with a dark background and the person's label",
          svg.startswith("<svg") and "#15171c" in svg and ">S</text>" in svg and svg.rstrip().endswith("</svg>"), svg[:200])
    txt = open(outs["ancient.txt"]).read() if p.returncode == 0 else ""
    check("report: the text carries the caveat and says qpAdm is not run", ar.CAVEAT in txt and "qpAdm" in txt, txt)
    # --- 8. both reports read the outputs ----------------------------------------------
    sd = os.path.join(tmp, "S")
    os.makedirs(os.path.join(sd, "ancient"))
    for k, f in outs.items():
        shutil.copy(f, os.path.join(sd, "ancient", f"S_{k.replace('pca.svg', 'ancient_pca.svg').replace('info.tsv', 'ancient_info.tsv')}"))
    html_out, txt_out = os.path.join(tmp, "S.html"), os.path.join(tmp, "S.txt")
    p = subprocess.run([sys.executable, os.path.join(REPO, "bin", "render_report.py"), "--sample", "S", "--sample-dir", sd,
                        "--json", os.path.join(tmp, "summary.json"), "-o", html_out, "-o", txt_out],
                       capture_output=True, text=True)
    check("render_report exits 0", p.returncode == 0, p.stderr)
    page = open(html_out).read() if os.path.exists(html_out) else ""
    check("the HTML report has the card, its first group and the plot as an image",
          "<h2>Ancient Populations</h2>" in page and "Steppe_EBA" in page and "data:image/svg+xml;base64," in page, page[-400:])
    text = open(txt_out).read() if os.path.exists(txt_out) else ""
    check("the text report has the section with the nearest group and the caveat",
          "Ancient Populations (AADR PCA, step 44)" in text and "1. Steppe_EBA (" in text and "not descent" in text, text[-600:])

    p = subprocess.run([sys.executable, SCRIPT, "report", "--panel-evec", pe, "--run-evec", pe, "--ancients", at,
                        "--sample", "S", "--release", "test", "--out-tsv", outs["ancient.tsv"], "--out-info",
                        outs["info.tsv"], "--out-txt", outs["ancient.txt"], "--out-svg", outs["pca.svg"]],
                       capture_output=True, text=True)
    check("report: a run without the person stops", p.returncode == 1 and "no row for S" in p.stderr, p.stderr)
finally:
    shutil.rmtree(tmp)

print(f"\n{FAILS} failure(s)")
sys.exit(1 if FAILS else 0)
