#!/usr/bin/env python3
"""bin/copy_number.py and bin/cn_loci.tsv: one depth-ratio engine for every
copy-number locus.

  1. the catalog: every row loads with its columns checked; a broken row
     (a bad kind, a segment without a span, a class call without an
     interval, a reused segment name) is refused;
  2. one owner per locus: no row's targets or controls overlap any region
     (main or paralog copy) of a model the pinned Parascopy data holds
     (assets/parascopy_models_v1.7.txt), except the pairs EXEMPT lists; a
     planted AMY1C row and a planted row on C4A's copy (C4B) are caught;
  3. parity: the depths of tests/test_hba_copy_number.py and
     tests/test_lpa_kiv2.py give, through the engine, the very tables the
     code before the engine wrote (GOLD, frozen from it);
  4. synthetic bedcov tables at 0, 1, 2, 3 and 4 copies of a made-up gene
     row are each called right; a MAPQ 0 share above the row's ceiling is
     not assessed; a kind=interval row with pooled depths works the same;
  5. a control segment with no read while the other controls are covered
     raises BrokenControl, and `call` exits 3 naming BROKEN_CONTROL; with no
     read anywhere it is 'not assessed', not an error;
  6. the CLI: `bed` prints every bedcov locus, `call --json` writes one block;
  7. the wiring: steps 35, 41 and 46 and their Nextflow modules call
     copy_number.py, never hba_copy_number.py or lpa_kiv2.py, with the
     bedcov options of the catalog row. UGT2B17 has no row: pypgx calls it
     in BAM mode (assets/pgx_genes.tsv), and tests/test_pypgx_genes.sh
     checks that no catalog row names a gene pypgx calls from the BAM;
  8. the gene-dosage rows (HP, LCE3BC, CCL3L1): made-up depths at 0
     to 4 copies are each called right; the block of step 35 that counts
     them, cut out and run on tests/fixtures/cn/, writes the known count of
     each, and a row shifted by 1 kb no longer matches its fixture; the
     PARASCOPY module runs the same list; the CCL3L1 reading names no
     disease (a planted 'HIV' is caught);
  9. the interval rows (scripts/ci/build-cn-intervals.py): a small ClinGen
     BED with a 'chrx' row converts to rows on chrX that load, the same rows
     left as 'chrx' are refused, and an unknown contig stops the build; the
     committed rows are the ClinGen regions plus the named hand rows, with
     no target in a breakpoint row; made-up depths call a 22q11.2 A-D
     deletion (1 copy), a 17p12 duplication (3) and 2 for both unaltered;
     AZFc on a male at half its DAZ depth is a deletion, AZFa a male's one
     copy, every chrY row on a female not applicable and without a sex not
     assessed; `intervals` over one bedcov table writes one line per row and
     a row with a broken control does not stop the others; the interval
     block of step 35, cut out and run, writes that table, and the PARASCOPY
     module runs the same commands.

Run: python3 tests/test_copy_number.py
"""
import json
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "copy_number.py")
CATALOG = os.path.join(REPO, "bin", "cn_loci.tsv")
MODELS = os.path.join(REPO, "assets", "parascopy_models_v1.7.txt")
sys.path.insert(0, os.path.join(REPO, "bin"))
import copy_number as cn  # noqa: E402
import hba_copy_number as hba  # noqa: E402
import lpa_kiv2 as lpa  # noqa: E402

FAILS = 0
# Frozen from bin/hba_copy_number.py and bin/lpa_kiv2.py as they were before the
# engine (main at 6056452): inputs [depths..., MAPQ 0 share] and the dict assess() gave.
GOLD = json.loads(''.join((
    '{"hba":[[[30,30,30,30,0.25],{"hba2_depth":30,"hba1_depth":30,"flank_left_depth":30,"flank_right_depth":30,"map'
    'q0_share":"0.250","estimate":"4.00","alpha_globin_genes":"4","status":"called","message":"four alpha-globin ge'
    'nes (aa/aa): no deletion seen"}],[[15,30,30,30,0.25],{"hba2_depth":15,"hba1_depth":30,"flank_left_depth":30,"f'
    'lank_right_depth":30,"mapq0_share":"0.250","estimate":"3.00","alpha_globin_genes":"3","status":"called","messa'
    'ge":"three alpha-globin genes: one gene lost on one chromosome (-a/aa, as -a3.7 or -a4.2), a silent carrier; u'
    'nconfirmed: this short-read count has been checked on one real genome only (HG002, four genes), so a count und'
    'er four needs a clinical test before it means anything"}],[[15,15,30,30,0.25],{"hba2_depth":15,"hba1_depth":15'
    ',"flank_left_depth":30,"flank_right_depth":30,"mapq0_share":"0.250","estimate":"2.00","alpha_globin_genes":"2"'
    ',"status":"called","message":"two alpha-globin genes: alpha-thalassaemia trait, either both genes lost on one '
    "chromosome (--/aa, as --SEA; a child with a partner who carries the same can have Hb Bart's) or one on each (-"
    'a/-a); depth cannot tell these apart; unconfirmed: this short-read count has been checked on one real genome o'
    'nly (HG002, four genes), so a count under four needs a clinical test before it means anything"}],[[45,30,30,30'
    ',0.25],{"hba2_depth":45,"hba1_depth":30,"flank_left_depth":30,"flank_right_depth":30,"mapq0_share":"0.250","es'
    'timate":"5.00","alpha_globin_genes":"5","status":"called","message":"5 alpha-globin genes: more than four, a t'
    'riplication (as aaa/aa)"}],[[22.5,30,30,30,0.25],{"hba2_depth":22.5,"hba1_depth":30,"flank_left_depth":30,"fla'
    'nk_right_depth":30,"mapq0_share":"0.250","estimate":"3.50","alpha_globin_genes":".","status":"not assessed","m'
    'essage":"not assessed: the estimate 3.50 lies between 3 and 4 alpha-globin genes"}],[[5,5,6,6,0.25],{"hba2_dep'
    'th":5,"hba1_depth":5,"flank_left_depth":6,"flank_right_depth":6,"mapq0_share":"0.250","estimate":".","alpha_gl'
    'obin_genes":".","status":"not assessed","message":"not assessed: flank depth under 8x (6.0x and 6.0x), too low'
    ' to compare"}],[[30,30,30,15,0.25],{"hba2_depth":30,"hba1_depth":30,"flank_left_depth":30,"flank_right_depth":'
    '15,"mapq0_share":"0.250","estimate":".","alpha_globin_genes":".","status":"not assessed","message":"not assess'
    'ed: the two flanks disagree (30.0x and 15.0x); a deletion larger than --SEA, or uneven coverage, would skew th'
    'e count"}],[[15.5,25.68,30.11,31.25,0],{"hba2_depth":15.5,"hba1_depth":25.68,"flank_left_depth":30.11,"flank_r'
    'ight_depth":31.25,"mapq0_share":"0.000","estimate":"2.68","alpha_globin_genes":".","status":"not assessed","me'
    'ssage":"not assessed: reads of mapping quality 0 make up 0.0% of the depth over HBA1 and HBA2, so the aligner '
    'left out the reads that fit both genes (as novoalign does by default) and the count would come out low; realig'
    'n with minimap2 (step 02) or BWA-MEM"}],[[30,30,30,30,0.01],{"hba2_depth":30,"hba1_depth":30,"flank_left_depth'
    '":30,"flank_right_depth":30,"mapq0_share":"0.010","estimate":"4.00","alpha_globin_genes":".","status":"not ass'
    'essed","message":"not assessed: reads of mapping quality 0 make up 1.0% of the depth over HBA1 and HBA2, so th'
    'e aligner left out the reads that fit both genes (as novoalign does by default) and the count would come out l'
    'ow; realign with minimap2 (step 02) or BWA-MEM"}],[[0,0,30,30,0],{"hba2_depth":0,"hba1_depth":0,"flank_left_de'
    'pth":30,"flank_right_depth":30,"mapq0_share":"0.000","estimate":"0.00","alpha_globin_genes":"0","status":"call'
    'ed","message":"no alpha-globin gene: the pattern of Hb Bart\'s hydrops fetalis (--/--); unconfirmed: this short'
    '-read count has been checked on one real genome only (HG002, four genes), so a count under four needs a clinic'
    'al test before it means anything"}],[[7.4,7.6,30,30,0.25],{"hba2_depth":7.4,"hba1_depth":7.6,"flank_left_depth'
    '":30,"flank_right_depth":30,"mapq0_share":"0.250","estimate":"1.00","alpha_globin_genes":"1","status":"called"'
    ',"message":"one alpha-globin gene: the pattern of HbH disease (--/-a); unconfirmed: this short-read count has '
    'been checked on one real genome only (HG002, four genes), so a count under four needs a clinical test before i'
    't means anything"}],[[30.2,29.1,29.4,31.9,0.2],{"hba2_depth":30.2,"hba1_depth":29.1,"flank_left_depth":29.4,"f'
    'lank_right_depth":31.9,"mapq0_share":"0.200","estimate":"3.87","alpha_globin_genes":"4","status":"called","mes'
    'sage":"four alpha-globin genes (aa/aa): no deletion seen"}]],"lpa":[[[12,30,0,0,0.9],{"status":"estimated","ki'
    'v2_copies":"12.0","ci95_low":"12.0","ci95_high":"12.0","class":"low","kiv2_depth":"30.00","flank_depth":"30.00'
    '","mapq0_share":"0.900","message":"about 12 KIV-2 repeats in total (95% interval 12.0 to 12.0): fewer KIV-2 re'
    'peats than most people (the whole 95% interval is under 30; about 70% of the 1000 Genomes samples had 30 to 45'
    '): a low count goes with small apo(a) isoforms and a higher Lp(a) level"}],[[20,30,0.1,0.15,0.9],{"status":"es'
    'timated","kiv2_copies":"19.8","ci95_low":"17.9","ci95_high":"22.0","class":"low","kiv2_depth":"47.45","flank_d'
    'epth":"30.43","mapq0_share":"0.900","message":"about 20 KIV-2 repeats in total (95% interval 17.9 to 22.0): fe'
    'wer KIV-2 repeats than most people (the whole 95% interval is under 30; about 70% of the 1000 Genomes samples '
    'had 30 to 45): a low count goes with small apo(a) isoforms and a higher Lp(a) level"}],[[29,30,0.03,0.03,0.9],'
    '{"status":"estimated","kiv2_copies":"29.0","ci95_low":"28.3","ci95_high":"29.7","class":"low","kiv2_depth":"66'
    '.56","flank_depth":"30.09","mapq0_share":"0.900","message":"about 29 KIV-2 repeats in total (95% interval 28.3'
    ' to 29.7): fewer KIV-2 repeats than most people (the whole 95% interval is under 30; about 70% of the 1000 Gen'
    'omes samples had 30 to 45): a low count goes with small apo(a) isoforms and a higher Lp(a) level"}],[[29,30,0.'
    '25,0.25,0.9],{"status":"estimated","kiv2_copies":"28.8","ci95_low":"23.5","ci95_high":"35.1","class":"typical"'
    ',"kiv2_depth":"67.50","flank_depth":"30.72","mapq0_share":"0.900","message":"about 29 KIV-2 repeats in total ('
    '95% interval 23.5 to 35.1): within or overlapping the range most people have (30 to 45 for about 70% of the 10'
    '00 Genomes samples)"}],[[47,30,0.25,0.25,0.9],{"status":"estimated","kiv2_copies":"46.6","ci95_low":"38.3","ci'
    '95_high":"56.6","class":"typical","kiv2_depth":"106.69","flank_depth":"30.72","mapq0_share":"0.900","message":'
    '"about 47 KIV-2 repeats in total (95% interval 38.3 to 56.6): within or overlapping the range most people have'
    ' (30 to 45 for about 70% of the 1000 Genomes samples)"}],[[60,30,0.03,0.03,0.9],{"status":"estimated","kiv2_co'
    'pies":"59.9","ci95_low":"58.6","ci95_high":"61.4","class":"high","kiv2_depth":"133.11","flank_depth":"30.09","'
    'mapq0_share":"0.900","message":"about 60 KIV-2 repeats in total (95% interval 58.6 to 61.4): more KIV-2 repeat'
    's than most people (the whole 95% interval is over 45): a high count goes with large apo(a) isoforms and, on a'
    'verage, a lower Lp(a) level"}],[[20,8,0,0,0.9],{"status":"too little depth","kiv2_copies":".","ci95_low":".","'
    'ci95_high":".","class":".","kiv2_depth":"12.57","flank_depth":"8.00","mapq0_share":"0.900","message":"too litt'
    'le depth: 8.0x over the single-copy LPA exons, under the 10x the count needs"}],[[20,30,0,0,0],{"status":"not '
    'assessed","kiv2_copies":".","ci95_low":".","ci95_high":".","class":".","kiv2_depth":"47.14","flank_depth":"30.'
    '00","mapq0_share":"0.000","message":"not assessed: reads of mapping quality 0 make up 0.0% of the depth over t'
    'he KIV-2 exons, so the aligner left out the reads that fit several KIV-2 copies (as novoalign does by default)'
    ' and the count would come out low; realign with minimap2 (step 02) or BWA-MEM"}],[[-1.5,30,0,0,0.9],{"status":'
    '"not assessed","kiv2_copies":".","ci95_low":".","ci95_high":".","class":".","kiv2_depth":"1.07","flank_depth":'
    '"30.00","mapq0_share":"0.900","message":"not assessed: the KIV-2 exons have 1.1x against 30.0x over the single'
    '-copy exons, too little for any repeat count (a deletion or a mapping problem)"}],[[38,12,0.3,0.2,0.5],{"statu'
    's":"estimated","kiv2_copies":"38.0","ci95_low":"31.0","ci95_high":"46.5","class":"typical","kiv2_depth":"34.95'
    '","flank_depth":"12.23","mapq0_share":"0.500","message":"about 38 KIV-2 repeats in total (95% interval 31.0 to'
    ' 46.5): within or overlapping the range most people have (30 to 45 for about 70% of the 1000 Genomes samples)"'
    '}],[[-2,30,0,0,0.9],{"status":"not assessed","kiv2_copies":".","ci95_low":".","ci95_high":".","class":".","kiv'
    '2_depth":"0.00","flank_depth":"30.00","mapq0_share":"0.000","message":"not assessed: reads of mapping quality '
    '0 make up 0.0% of the depth over the KIV-2 exons, so the aligner left out the reads that fit several KIV-2 cop'
    'ies (as novoalign does by default) and the count would come out low; realign with minimap2 (step 02) or BWA-ME'
    'M"}]]}'
)))


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS += 1


HEADER = "\t".join(cn.COLUMNS)


def row(**kw):
    """A made-up bedcov row: one 1 kb gene between two 10 kb flanks, called
    by rounding, with a MAPQ 0 ceiling of 0.5."""
    r = {"locus": "TESTG", "kind": "gene", "chrom": "chr3", "targets": "TESTG:50000-51000",
         "controls": "fl_a:30000-40000;fl_b:60000-70000", "details": ".", "target_depth": "sum",
         "control_depth": "flanks", "scale": "2", "offset": "0", "bedcov": ".,-Q 1", "mapq0_min": ".",
         "mapq0_max": "0.5", "min_control_depth": "8", "max_control_ratio": "1.5", "call": "round:0.3",
         "interval": ".", "readings": "0=none|1=one|2=two|more={n} copies", "input": "bedcov",
         "owner": "cn-depth", "cross_check": "none", "pmids": "."}
    r.update(kw)
    return r


def write_catalog(path, rows):
    with open(path, "w") as f:
        f.write(HEADER + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in cn.COLUMNS) + "\n")


def bedcov(path, prow, depth):
    with open(path, "w") as f:
        for c, s, e, n in cn.regions(prow):
            f.write(f"{c}\t{s}\t{e}\t{n}\t{round(depth[n] * (e - s))}\n")


# The only rows allowed to overlap a Parascopy model, each for a stated reason.
# A new pair here needs the same: a cross_check in the catalog does not exempt a row.
EXEMPT = {
    # Step 41 runs by default and Parascopy is opt-in, so the RhD deletion stays a
    # depth row; RHCE's model holds RHD as a paralog copy and is its cross-check.
    ("RHD", "RHCE"),
}
DOSAGE = ["HP", "LCE3BC", "CCL3L1"]
# The counts tests/fixtures/cn/<locus>.bedcov.tsv was made for.
FIXTURE_COPIES = {"HP": 3, "LCE3BC": 1, "CCL3L1": 5}
# Words a reading of CCL3L1 must not carry: its disease links are disputed.
DENY = re.compile(r"hiv|aids|infect|disease|risk|suscept|progress|malaria|hepatitis|arthritis|lupus|kawasaki|"
                  r"psoria|diabet|crohn|vaccin", re.I)


def overlaps(loc_row, models):
    """[(model, segment)] of a parsed row's targets and controls that overlap any
    region of a model (1-based, both ends included, against 0-based segments)."""
    hits = []
    for n, s, e, _ in loc_row["_targets"] + loc_row["_controls"]:
        for name, regs in models.items():
            if any(c == loc_row["chrom"] and s < me and ms - 1 < e for c, ms, me, _ in regs):
                hits.append((name, n))
    return hits


def owner_violations(cat, models):
    return [(loc, m, seg) for loc, r in cat.items() for m, seg in overlaps(r, models) if (loc, m) not in EXEMPT]


def main():
    # 1. the catalog
    cat = cn.load_catalog(CATALOG)
    check("catalog: the six gene and dosage loci first, then the interval rows",
          list(cat)[:6] == ["HBA", "LPA_KIV2", "RHD"] + DOSAGE
          and all(r["kind"] == "interval" and r["_passes"] == ["."] for r in list(cat.values())[6:]), list(cat)[:8])
    for loc, r in cat.items():
        check(f"catalog {loc}: owner cn-depth, cross_check named, PMIDs numbers or '.'",
              r["owner"] == "cn-depth" and re.fullmatch(r"none|[a-z]+:[A-Za-z0-9_]+(,[a-z]+:[A-Za-z0-9_]+)*",
                                                        r["cross_check"])
              and re.fullmatch(r"\.|\d+(,\d+)*", r["pmids"]), r)
        if r["input"] == "bedcov":
            passes = r["_passes"]
            check(f"catalog {loc}: bedcov passes ('.' or samtools options; a second pass adds -Q 1)",
                  all(re.fullmatch(r"\.|(-[A-Za-z]( \d+)?)( -[A-Za-z]( \d+)?)*", p) for p in passes)
                  and (len(passes) == 1 or passes[1].replace(".", "").replace("-Q 1", "").strip()
                       == passes[0].replace(".", "").strip()), passes)
            segs = r["_targets"] + r["_controls"] + r["_details"]
            check(f"catalog {loc}: every segment is a span of under 1 Mb", all(0 < e - s < 1_000_000 for _, s, e, _ in segs))
            ctl = r["_controls"]
            check(f"catalog {loc}: no control overlaps a target",
                  not any(s < te and ts < e for _, s, e, _ in ctl for _, ts, te, _ in r["_targets"]))
            readings = r["_readings"]
            mode, vals = r["_call"]
            want = {"low", "typical", "high"} if mode == "class" else {str(i) for i in range(int(len(vals) + 1))} \
                if mode == "cut" else None
            check(f"catalog {loc}: a reading for every call", want is None or want <= set(readings), readings)
    with tempfile.TemporaryDirectory() as t:
        p = os.path.join(t, "c.tsv")
        for desc, kw, err in (("a bad kind", {"kind": "exon"}, "kind"),
                              ("a segment with no span", {"targets": "TESTG:51000-50000"}, "empty span"),
                              ("a class call without an interval", {"call": "class:30,45"}, "needs interval cv"),
                              ("a MAPQ 0 guard with one bedcov pass", {"bedcov": "."}, "two bedcov passes"),
                              ("a segment name used twice", {"details": "fl_a:1-2"}, "used twice"),
                              ("a NaN call threshold", {"call": "round:nan"}, "is not round"),
                              ("an infinite cut threshold", {"call": "cut:0.5,inf"}, "is not round"),
                              ("a NaN number column", {"mapq0_max": "nan"}, "not a number")):
            write_catalog(p, [row(**kw)])
            try:
                cn.load_catalog(p)
                got = "loaded"
            except cn.CatalogError as e:
                got = str(e)
            check(f"catalog: {desc} is refused", err in got, got)
        write_catalog(p, [row(), row(locus="TESTH", targets="TESTG:80000-81000", controls="fl_c:1-2;fl_d:3-4")])
        try:
            cn.load_catalog(p)
            got = "loaded"
        except cn.CatalogError as e:
            got = str(e)
        check("catalog: one segment name in two loci is refused", "used by two loci" in got, got)

        g0 = cn.parse_row(row())
        b = os.path.join(t, "dup.bedcov")
        bedcov(b, g0, {"fl_a": 30.0, "fl_b": 30.0, "TESTG": 30.0})
        with open(b) as f:
            first = f.readline()
        with open(b, "a") as f:
            f.write(first.rsplit("\t", 1)[0] + "\t0\n")
        try:
            cn.read_bedcov(b, g0)
            got = "read"
        except ValueError as e:
            got = str(e)
        check("bedcov: a region listed twice is refused, not overwritten", "appears twice" in got, got)
        with open(b, "w") as f:
            f.write(first.split("\t")[0] + "\t50000\t50000\tTESTG\t0\n")
        try:
            cn.read_bedcov(b, g0)
            got = "read"
        except ValueError as e:
            got = str(e)
        except ZeroDivisionError as e:
            got = f"ZeroDivisionError: {e}"
        check("bedcov: an empty region is refused with its name", "is empty or reversed" in got, got)

    # 2. one owner per locus
    models = {}
    with open(MODELS) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            name, reg, what = line.split()
            c, span = reg.split(":")
            s, e = span.split("-")
            models.setdefault(name, []).append((c, int(s), int(e), what))
    check("Parascopy model list: 167 models, AMY1C, C4A, CFH, CFHR4, FCGR3A, CHRNA7, APOBEC3A, GYPA, RHCE among them",
          len(models) == 167 and {"AMY1C", "C4A", "CFH", "CFHR4", "FCGR3A", "CHRNA7", "APOBEC3A", "GYPA", "RHCE"}
          <= set(models), len(models))
    check("Parascopy model list: one main region per model, and the paralog copies (over 1,000 regions in all)",
          all(sum(w == "main" for *_, w in r) == 1 and all(w in ("main", "copy") for *_, w in r)
              for r in models.values()) and sum(len(r) for r in models.values()) > 1000,
          sum(len(r) for r in models.values()))
    for loc, r in cat.items():
        for x in r["cross_check"].split(","):
            if x.startswith("parascopy:"):
                check(f"catalog {loc}: cross_check {x} is a pinned model", x.split(":", 1)[1] in models, x)
    v = owner_violations(cat, models)
    check("owner: no row counts a gene a pinned Parascopy model owns", not v, v)
    stale = [x for x in EXEMPT if x[0] not in cat or x[1] not in {m for m, _ in overlaps(cat[x[0]], models)}]
    check("owner: every exemption is still needed (its row overlaps its model)", not stale, stale)
    planted = dict(cat)
    planted["AMY1"] = cn.parse_row(row(locus="AMY1", chrom="chr1", targets="AMY1C:103745500-103760000",
                                       controls="amy_l:103000000-103010000;amy_r:104000000-104010000"))
    v = owner_violations(planted, models)
    check("owner: a planted AMY1C row is caught", v == [("AMY1", "AMY1C", "AMY1C")], v)
    planted = dict(cat)
    planted["C4B"] = cn.parse_row(row(locus="C4B", chrom="chr6", targets="C4B:32014800-32035400",
                                      controls="c4_l:31000000-31010000;c4_r:33500000-33510000"))
    v = owner_violations(planted, models)
    check("owner: a planted row on C4A's paralog copy (C4B) is caught", v == [("C4B", "C4A", "C4B")], v)
    planted = dict(cat)
    planted["GYPA"] = cn.parse_row(row(locus="GYPA", chrom="chr4", cross_check="parascopy:GYPA",
                                       targets="GYPA_T:144140000-144150000",
                                       controls="gy_l:140000000-140010000;gy_r:150000000-150010000"))
    v = owner_violations(planted, models)
    check("owner: naming the model in cross_check does not exempt a row", ("GYPA", "GYPA", "GYPA_T") in v, v)
    # DEFB4 has no row: the beta-defensin repeat (DEFB4B to DEFB103B) lies in SPAG11A's paralog copies.
    planted = dict(cat)
    planted["DEFB4"] = cn.parse_row(row(locus="DEFB4", chrom="chr8", targets="DEFB_UNIT_B:7414854-7430348",
                                        controls="defb_l:7622000-7662000;defb_r:8250000-8290000"))
    v = owner_violations(planted, models)
    check("owner: a planted DEFB4 row on the beta-defensin repeat is caught (SPAG11A's model)",
          v == [("DEFB4", "SPAG11A", "DEFB_UNIT_B")], v)
    check("owner: RHD overlaps RHCE's model only (a paralog copy), the exempt pair",
          {m for m, _ in overlaps(cat["RHD"], models)} == {"RHCE"}, overlaps(cat["RHD"], models))

    # 3. parity with the code before the engine
    for (h2, h1, left, right, share), want in GOLD["hba"]:
        m = {"flank_left": left, "HBA2": h2, "HBA1": h1, "flank_right": right}
        q = {n: v * (1 - share) if n in ("HBA2", "HBA1") else v for n, v in m.items()}
        got = hba.assess(m, q)
        check(f"parity HBA {h2}/{h1} flanks {left}/{right} MAPQ0 {share}: {want['status']} {want['alpha_globin_genes']}",
              got == want, (got, want))
    kiv = [n for *_, n in lpa.REGIONS if n.startswith("kiv2_")]
    fl = [n for *_, n in lpa.REGIONS if n.startswith("flank_")]
    for (n, flank, ks, fs, share), want in GOLD["lpa"]:
        k = flank * (n + 2) / 14
        m = {x: k * (1 + (ks if i % 2 == 0 else -ks)) for i, x in enumerate(kiv)}
        m.update({x: flank * (1 + (fs if i % 2 == 0 else -fs)) for i, x in enumerate(fl)})
        q = {x: v * (1 - share) if x.startswith("kiv2_") else v for x, v in m.items()}
        got = lpa.assess(m, q)
        check(f"parity LPA {n} copies at {flank}x: {want['status']} {want['kiv2_copies']}", got == want, (got, want))

    # 4. synthetic copies, the MAPQ 0 ceiling, a kind=interval row
    g = cn.parse_row(row())
    for copies in (0, 1, 2, 3, 4):
        m = {"fl_a": 30.0, "fl_b": 30.0, "TESTG": 15.0 * copies}
        b = cn.assess(g, m, {n: v * 0.9 for n, v in m.items()})
        check(f"synthetic: {copies} copies called {copies}", (b["status"], b["copies"]) == ("called", copies), b)
    m = {"fl_a": 30.0, "fl_b": 30.0, "TESTG": 30.0}
    b = cn.assess(g, m, {"fl_a": 30.0, "fl_b": 30.0, "TESTG": 12.0})
    check("synthetic: 60% of the target from MAPQ 0, over the 0.5 ceiling: not assessed",
          (b["status"], b["reason"]) == ("not assessed", "high_mapq0_share"), b)
    b = cn.assess(g, m, {"fl_a": 30.0, "fl_b": 30.0, "TESTG": 16.0})
    check("synthetic: 47% from MAPQ 0, under the ceiling: called 2", (b["status"], b["copies"]) == ("called", 2), b)
    m = {"fl_a": 30.0, "fl_b": 30.0, "TESTG": 37.5}
    b = cn.assess(g, m, m)
    check("synthetic: 2.5 copies is between counts", b["reason"] == "between_counts", b)
    z = cn.parse_row(row(min_control_depth="."))
    try:
        b = cn.assess(z, {"fl_a": 0.0, "fl_b": 0.0, "TESTG": 30.0}, {"fl_a": 0.0, "fl_b": 0.0, "TESTG": 30.0})
        got = b["reason"]
    except ZeroDivisionError as e:
        got = f"ZeroDivisionError: {e}"
    check("synthetic: no control depth with no minimum set is low_control_depth, not a crash",
          got == "low_control_depth", got)
    try:
        b = cn.assess(z, {"fl_a": 0.0, "fl_b": 30.0, "TESTG": 30.0}, {"fl_a": 0.0, "fl_b": 30.0, "TESTG": 30.0})
        got = b["reason"]
    except ZeroDivisionError as e:
        got = f"ZeroDivisionError: {e}"
    check("synthetic: a flank with no read, no minimum set, is controls_disagree, not a crash",
          got == "controls_disagree", got)
    iv = cn.parse_row(row(kind="interval", target_depth="pooled", control_depth="pooled", mapq0_max=".",
                          bedcov=".", targets="seg1:50000-50500;seg2:50500-51000"))
    b = cn.assess(iv, {"fl_a": 20.0, "fl_b": 40.0, "seg1": 45.0, "seg2": 45.0})
    check("interval row: pooled 45x over pooled 30x controls is 3 copies", (b["status"], b["copies"]) == ("called", 3), b)
    b = cn.assess(cat["RHD"], {n: 30.0 for *_, n in cn.regions(cat["RHD"])})
    check("a row without a MAPQ 0 guard, given no MAPQ >= 1 depths, reports no MAPQ 0 share",
          b["status"] == "called" and b["mapq0_share"] is None, b)
    try:
        cn.assess(g, m)
        got = "no error"
    except ValueError as e:
        got = str(e)
    check("a row with a MAPQ 0 guard, given no MAPQ >= 1 depths, is refused", "MAPQ >= 1 depths are needed" in got, got)

    # 5. a broken control
    bad = {"fl_a": 0.0, "fl_b": 30.0, "TESTG": 30.0}
    try:
        cn.assess(g, bad, bad)
        got = "no error"
    except cn.BrokenControl as e:
        got = str(e)
    check("broken control: fl_a with no read beside fl_b at 30x raises BROKEN_CONTROL",
          got.startswith("BROKEN_CONTROL: control segment fl_a of TESTG"), got)
    b = cn.assess(g, {"fl_a": 0.0, "fl_b": 0.0, "TESTG": 0.0}, {})
    check("no read anywhere: not assessed (too little depth), not an error",
          (b["status"], b["reason"]) == ("not assessed", "low_control_depth"), b)
    b = cn.assess(g, {"fl_a": 0.0, "fl_b": 5.0, "TESTG": 5.0}, {})
    check("one flank empty beside a thin one (5x): too little depth, not an error", b["reason"] == "low_control_depth", b)
    lrow = cat["LPA_KIV2"]
    with tempfile.TemporaryDirectory() as t:
        a, q, j = (os.path.join(t, x) for x in ("a.tsv", "q.tsv", "o.json"))
        depth = {n: 30.0 for *_, n in cn.regions(lrow)}
        depth["flank_e19"] = 0.0
        bedcov(a, lrow, depth)
        bedcov(q, lrow, depth)
        out = subprocess.run([sys.executable, SCRIPT, "call", "--locus", "LPA_KIV2", "--bedcov", a,
                              "--bedcov-mapq1", q, "--json", j], capture_output=True, text=True)
        check("broken control: `call` exits 3 with BROKEN_CONTROL and writes nothing",
              out.returncode == 3 and "BROKEN_CONTROL: control segment flank_e19 of LPA_KIV2" in out.stderr
              and not os.path.exists(j), out)

        # 6. the CLI
        out = subprocess.run([sys.executable, SCRIPT, "bed"], capture_output=True, text=True)
        names = [x.split("\t")[3] for x in out.stdout.splitlines()]
        want = sum(len(cn.regions(r)) for r in cat.values() if r["input"] == "bedcov")
        check("bed: every segment of the bedcov loci, once", out.returncode == 0 and len(names) == len(set(names)) == want,
              (out.returncode, len(names), want))
        depth = {n: 30.0 for *_, n in cn.regions(lrow)}
        for n in depth:
            if n.startswith("kiv2_"):
                depth[n] = 30.0 * 22 / 14
        bedcov(a, lrow, depth)
        bedcov(q, lrow, {n: v * (0.1 if n.startswith("kiv2_") else 1) for n, v in depth.items()})
        tab = os.path.join(t, "t.tsv")
        out = subprocess.run([sys.executable, SCRIPT, "call", "--locus", "LPA_KIV2", "--bedcov", a,
                              "--bedcov-mapq1", q, "--json", j, "--table", tab], capture_output=True, text=True)
        block = json.load(open(j)) if os.path.exists(j) else {}
        check("call --json: one block, 20 KIV-2 repeats, called",
              out.returncode == 0 and block.get("locus") == "LPA_KIV2" and block.get("status") == "called"
              and abs(block.get("estimate", 0) - 20) < 0.01, (out, block))
        tab_rows = dict(x.split("\t", 1) for x in open(tab).read().splitlines()[1:]) if os.path.exists(tab) else {}
        check("call --table: step 46's table", tab_rows.get("kiv2_copies") == "20.0" and list(tab_rows) == list(lpa.ORDER),
              tab_rows)
        out = subprocess.run([sys.executable, SCRIPT, "call", "--locus", "LPA_KIV2", "--bedcov", a, "--json", j],
                             capture_output=True, text=True)
        check("call: a row with a MAPQ 0 guard needs --bedcov-mapq1", out.returncode == 2 and "--bedcov-mapq1" in out.stderr,
              out)
        empty = os.path.join(t, "empty.tsv")
        open(empty, "w").close()
        for loc, key in (("HBA", "alpha_globin_genes"), ("LPA_KIV2", "kiv2_copies")):
            if os.path.exists(tab):
                os.remove(tab)
            out = subprocess.run([sys.executable, SCRIPT, "call", "--locus", loc, "--bedcov", empty,
                                  "--bedcov-mapq1", empty, "--table", tab], capture_output=True, text=True)
            tab_rows = dict(x.split("\t", 1) for x in open(tab).read().splitlines()[1:]) if os.path.exists(tab) else {}
            check(f"call {loc} on two empty tables (no reads on the contig): exit 0, not assessed",
                  out.returncode == 0 and tab_rows.get("status") == "not assessed" and tab_rows.get(key) == "."
                  and "no chr" in tab_rows.get("message", ""), (out, tab_rows))

    # 7. the wiring
    def read(p):
        with open(os.path.join(REPO, p)) as f:
            return f.read()
    for path, loci in (("scripts/35-paralogs.sh", ["HBA"]), ("modules/local/parascopy/main.nf", ["HBA"]),
                       ("scripts/41-bloodtype.sh", ["RHD"]), ("modules/local/bloodtype/main.nf", ["RHD"]),
                       ("scripts/46-lpa.sh", ["LPA_KIV2"]), ("modules/local/lpa/main.nf", ["LPA_KIV2"])):
        text = read(path)
        for loc in loci:
            check(f"{path}: copy_number.py bed --locus {loc}", f"copy_number.py bed --locus {loc}" in text)
            if loc != "RHD":   # step 41's report imports the engine
                check(f"{path}: copy_number.py call --locus {loc}", f"copy_number.py call --locus {loc}" in text)
            passes = ["" if p == "." else p for p in cat[loc]["_passes"]]
            if loc == "HBA":   # and the dosage loop's two passes (section 8), and the interval pass (section 9)
                passes += ["", "-Q 1", ""]
            # the commands only: options, then a path, a quoted argument or a variable
            calls = [o for o, arg in re.findall(r"samtools bedcov((?: -[A-Za-z](?: \d+)?)*) (\S+)", text)
                     if re.match(r'["$/]|\S+\.bed\b', arg)]
            check(f"{path}: its bedcov passes are the catalog's ({' and '.join(repr(p) for p in passes)})",
                  sorted(c.strip() for c in calls) == sorted(passes), calls)
        check(f"{path}: no call of the old per-gene scripts",
              not re.search(r"(hba_copy_number|lpa_kiv2)\.py (bed|check|estimate)", text))

    dosage(cat, read)
    intervals(cat, read)

    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


def dosage_depths(r, copies, q1_share=0.5):
    """Mean depths of a dosage row at a copy number: controls at 30x (two
    copies), the target copies spread evenly over the row's targets (the
    reference copies), and the MAPQ >= 1 depths with q1_share of the targets."""
    tg = {n for n, *_ in r["_targets"]}
    m = {n: 15.0 * copies / len(tg) if n in tg else 30.0 for *_, n in cn.regions(r)}
    return m, {n: v * q1_share if n in tg else v for n, v in m.items()}


# run_in and cpath for the dosage block of step 35: copy_number.py runs for
# real, samtools bedcov prints tests/fixtures/cn/<locus>[.q1].bedcov.tsv.
STUB = r"""
run_in() {
  while [ $# -gt 0 ]; do case $1 in -v|--cpus|--memory) shift 2 ;; *) break ;; esac; done
  shift
  case "$1 $2" in
    "python3 /pgp-bin/copy_number.py") shift 2; python3 "${REPO}/bin/copy_number.py" "$@" ;;
    "samtools bedcov") shift 2; q=""; if [ "$1" = -Q ]; then q=.q1; shift 2; fi
      d=$(dirname "$1"); cat "${FIX}/${d##*/}${q}.bedcov.tsv" ;;
    *) return 1 ;;
  esac
}
cpath() { printf '%s' "$1"; }
"""


def dosage(cat, read):
    step, module = read("scripts/35-paralogs.sh"), read("modules/local/parascopy/main.nf")
    m = re.search(r'^DOSAGE_LOCI="([^"]*)"$', step, re.M)
    check("step 35: DOSAGE_LOCI lists the dosage rows", m is not None and m.group(1).split() == DOSAGE,
          m and m.group(1))
    m = re.search(r'^ *dosage_loci="([^"]*)"$', module, re.M)
    check("PARASCOPY module: dosage_loci is step 35's list", m is not None and m.group(1).split() == DOSAGE,
          m and m.group(1))
    for path, text, q in (("scripts/35-paralogs.sh", step, '"'), ("modules/local/parascopy/main.nf", module, '"\\')):
        for cmd in ("bed", "call"):
            check(f"{path}: copy_number.py {cmd} --locus on each dosage row",
                  f"copy_number.py {cmd} --locus {q}$loc" in text)
    for loc in DOSAGE:
        r = cat[loc]
        check(f"dosage {loc}: two bedcov passes, all reads then MAPQ >= 1 (what the loop runs)",
              r["_passes"] == [".", "-Q 1"], r["_passes"])
        check(f"dosage {loc}: a whole-number call with a reading for 0 and 1 and for more",
              r["_call"][0] == "round" and {"0", "1", "more"} <= set(r["_readings"]), r["_readings"])
        check(f"dosage {loc}: a PMID", re.fullmatch(r"\d+(,\d+)*", r["pmids"]) is not None, r["pmids"])
        for copies in (0, 1, 2, 3, 4):
            m, q = dosage_depths(r, copies)
            b = cn.assess(r, m, q)
            check(f"dosage {loc}: {copies} copies called {copies}, with its reading",
                  (b["status"], b["copies"]) == ("called", copies) and b["reading"] == cn.reading(r, copies)
                  and "{n}" not in b["reading"], b)
        m, q = dosage_depths(r, 2.5)
        check(f"dosage {loc}: 2.5 copies is between counts", cn.assess(r, m, q)["reason"] == "between_counts")
    hp = cat["HP"]
    m, q = dosage_depths(hp, 4)
    m["HP_UNIT1"], m["HP_UNIT2"] = 60.0, 0.0
    check("dosage HP: reads of both copies on one reference copy still add up to 4",
          cn.assess(hp, m, q)["copies"] == 4)

    readings = " ".join(cat["CCL3L1"]["_readings"].values())
    check("dosage CCL3L1: no reading names a disease", not DENY.search(readings), DENY.findall(readings))
    planted = readings + " more copies go with slower HIV progression"
    check("dosage CCL3L1: a planted 'HIV' reading is caught", DENY.search(planted) is not None)

    # The block of step 35, cut out and run on the fixtures.
    blk = re.search(r'\n(DOSAGE_LOCI=.*?\nrm -rf "\$DOSAGE_WORK"\n)', step, re.S)
    check("step 35 holds the dosage block", blk is not None)
    fix = os.path.join(REPO, "tests", "fixtures", "cn")
    with tempfile.TemporaryDirectory() as t:
        env = dict(os.environ, REPO=REPO, FIX=fix, OUTDIR=t, SAMPLE="S", PGP_ROOT=REPO, PARASCOPY_IMAGE="img",
                   BAM=os.path.join(t, "S.bam"))
        out = subprocess.run(["bash", "-euo", "pipefail", "-c", STUB + (blk.group(1) if blk else "false\n")],
                             capture_output=True, text=True, env=env)
        tab = os.path.join(t, "S_cn_dosage.tsv")
        rows = {}
        if os.path.exists(tab):
            lines = open(tab).read().splitlines()
            rows = {x.split("\t")[0]: x.split("\t") for x in lines[1:]}
            check("step 35 dosage table: its header", lines[0] == "locus\tstatus\tcopies\testimate\treason\treading",
                  lines[0])
        check("step 35 dosage block: exit 0, one line per dosage row", out.returncode == 0 and list(rows) == DOSAGE,
              (out.returncode, out.stderr[-600:], list(rows)))
        for loc, n in FIXTURE_COPIES.items():
            got = rows.get(loc, [])
            check(f"step 35 on tests/fixtures/cn/{loc}: called {n}, with its reading",
                  got[1:3] == ["called", str(n)] and got[5:] == [cn.reading(cat[loc], n)], got)

    # A row shifted by 1 kb no longer matches the fixture made for its place.
    for loc in DOSAGE:
        r = dict(cat[loc])
        r["targets"] = ";".join(f"{n}:{s + 1000}-{e + 1000}" for n, s, e, _ in r["_targets"])
        shifted = cn.parse_row({k: v for k, v in r.items() if not k.startswith("_")})
        try:
            cn.read_bedcov(os.path.join(fix, f"{loc}.bedcov.tsv"), shifted)
            got = "read"
        except ValueError as e:
            got = str(e)
        check(f"dosage {loc}: shifted 1 kb, its fixture is refused", "regions differ" in got, got)


# Section 9: the interval rows.
BUILD = os.path.join(REPO, "scripts", "ci", "build-cn-intervals.py")
SMALL_BED = """track name="test"
#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb
chr22\t18159724\t18924717\tBP-A\t0\t.\t18159724\t18924717\t0,0,0
chr22\t18924718\t21111383\t22q11.2 recurrent (DGS) region (proximal, A-D) (includes TBX1)\t0\t.\t18924718\t21111383\t255,111,0
chr22\t21111384\t21562827\tBP-D\t0\t.\t21111384\t21562827\t0,0,0
chrx\t6527561\t6537770\tXp22.31 distal breakpoint\t0\t.\t6527561\t6537770\t0,0,0
chrx\t6537771\t8156913\tXp22.31 recurrent region (includes STS)\t0\t.\t6537771\t8156913\t255,111,0
"""
HAND = ["1p36_terminal", "Xq28_MECP2", "GJB6_D13S1830", "CFTR_dele2_3", "HBB_exon3_del", "HBB_619bp_del",
        "HBD_HBB_Sicilian_del", "AZFa", "AZFb", "AZFc_DAZ"]
STUB_IV = r"""
run_in() {
  while [ $# -gt 0 ]; do case $1 in -v|--cpus|--memory) shift 2 ;; *) break ;; esac; done
  shift
  case "$1 $2" in
    "python3 /pgp-bin/copy_number.py") shift 2; python3 "${REPO}/bin/copy_number.py" "$@" ;;
    "samtools bedcov") shift 2; python3 "${FAKECOV}" "$1" ;;
    *) return 1 ;;
  esac
}
cpath() { printf '%s' "$1"; }
"""
# fakecov BED: samtools bedcov over BED at 30x, the segments named in FAKE_HALF at 15x, in FAKE_ZERO at 0.
FAKECOV = r"""
import os, sys
half = set(os.environ.get("FAKE_HALF", "").split(","))
zero = set(os.environ.get("FAKE_ZERO", "").split(","))
for line in open(sys.argv[1]):
    c = line.rstrip("\n").split("\t")
    d = 0.0 if c[3] in zero else 15.0 if c[3] in half else 30.0
    print("\t".join(c + [str(int(d * (int(c[2]) - int(c[1]))))]))
"""


def depths(r, target=30.0, control=30.0):
    tg = {n for n, *_ in r["_targets"]}
    return {n: target if n in tg else control for *_, n in cn.regions(r)}


def intervals(cat, read):
    import importlib.util
    spec = importlib.util.spec_from_file_location("build_cn_intervals", BUILD)
    bld = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bld)
    rows = bld.build(SMALL_BED.splitlines(True), [], [], [], [])
    check("converter: the two region rows of the small BED, its chrx row on chrX",
          [(r["locus"], r["chrom"]) for r in rows[:2]] == [("22q11.2_DGS_proximal_A-D", "chr22"), ("Xp22.31", "chrX")],
          [(r["locus"], r["chrom"]) for r in rows])
    with tempfile.TemporaryDirectory() as t:
        p = os.path.join(t, "c.tsv")
        write_catalog(p, rows)
        try:
            got = len(cn.load_catalog(p))
        except cn.CatalogError as e:
            got = str(e)
        check("converter: its rows load", got == len(rows), got)
        write_catalog(p, [dict(r, chrom="chrx") if r["chrom"] == "chrX" else r for r in rows])
        try:
            cn.load_catalog(p)
            got = "loaded"
        except cn.CatalogError as e:
            got = str(e)
        check("converter: the same rows left as 'chrx' are refused", "chrom 'chrx'" in got, got)
    try:
        bld.build(SMALL_BED.replace("chrx", "chrQ").splitlines(True), [], [], [], [])
        got = "built"
    except bld.BuildError as e:
        got = str(e)
    check("converter: an unknown contig stops the build", "'chrQ' is not a GRCh38 chromosome" in got, got)
    check("converter: no target piece overlaps a breakpoint row",
          all(not (s < 18924717 and 18159724 < e) and not (s < 21562827 and 21111384 < e) for r in rows[:1]
              for _, s, e, _ in cn.parse_row(r)["_targets"]))

    iv = [r for r in cat.values() if r["kind"] == "interval" and r["_passes"] == ["."]]
    names = [r["locus"] for r in iv]
    check("interval rows: the ClinGen regions (over 70), then the hand rows", len(iv) > 80 and names[-10:] == HAND,
          (len(iv), names[-10:]))
    # HbVar's GRCh38 spans (979, 1035), 0-based start: one target each, never cut
    spans = {loc: [(s, e) for _, s, e, _ in cat[loc]["_targets"]] if loc in cat else None
             for loc in ("HBB_619bp_del", "HBD_HBB_Sicilian_del")}
    check("interval rows: the HBB 619 bp and Sicilian deletions at HbVar's GRCh38 spans",
          spans == {"HBB_619bp_del": [(5225388, 5226007)], "HBD_HBB_Sicilian_del": [(5219877, 5233280)]}, spans)
    check("interval rows: each cross-checked by CNVpytor, step 18",
          all(r["cross_check"] == "cnvpytor:cnvs" for r in iv))
    check("interval rows: the 15q13.3 rows keep CHRNA7 (a Parascopy model) out of their targets",
          all(not (32154858 < e and s < 32178018) for r in iv if r["locus"].startswith("15q13.3")
              for _, s, e, _ in r["_targets"])
          and any(r["locus"].startswith("15q13.3") for r in iv))
    hits = [(r["locus"], o["locus"]) for r in iv for _, cs, ce, _ in r["_controls"] for o in iv
            if o is not r and o["chrom"] == r["chrom"] and any(ts < ce and cs < te for _, ts, te, _ in o["_targets"])]
    check("interval rows: no row's control sits on another row's target (Xq28's left flank off Xq28_MECP2)",
          not hits, hits)

    a_d, p17 = cat["22q11.2_DGS_proximal_A-D"], cat["17p12_HNPP_CMT1A"]
    for r, target, want, change in ((a_d, 15.0, 1, "deletion"), (p17, 45.0, 3, "duplication"),
                                    (a_d, 30.0, 2, "no change"), (p17, 30.0, 2, "no change")):
        b = cn.assess(r, depths(r, target))
        line = cn.interval_line(r, b)
        check(f"interval {r['locus']}: targets at {target}x over 30x controls is {want} copies, {change}",
              (b["status"], b["copies"], line[8]) == ("called", want, change), (b, line))
    az = cat["AZFc_DAZ"]
    m = depths(az, 15.0 / 4, 15.0)
    check("AZFc on a male: four DAZ copies at a single copy's depth in all is 4 copies",
          cn.assess(az, {n: (15.0 if n in {x for x, *_ in az["_targets"]} else 15.0) for n in m}, None, "male")
          ["copies"] == 4)
    half = {n: (7.5 if n in {x for x, *_ in az["_targets"]} else 15.0) for n in m}
    b = cn.assess(az, half, None, "male")
    check("AZFc on a male at half its DAZ depth: 2 copies, a deletion (gr/gr or b2/b3)",
          (b["copies"], cn.interval_line(az, b)[8]) == (2, "deletion"), b)
    b = cn.assess(cat["AZFa"], depths(cat["AZFa"], 15.0, 15.0), None, "male")
    check("AZFa on a male: one copy, no change", (b["copies"], b["expected"]) == (1, 1), b)
    b = cn.assess(cat["AZFa"], depths(cat["AZFa"], 0.0, 15.0), None, "male")
    check("AZFa on a male with no read over it: 0 copies, a deletion",
          (b["copies"], cn.interval_line(cat["AZFa"], b)[8]) == (0, "deletion"), b)
    for loc in ("AZFa", "AZFb", "AZFc_DAZ"):
        b = cn.assess(cat[loc], depths(cat[loc], 15.0, 15.0), None, "female")
        check(f"{loc} on a female: not applicable", (b["status"], b["reason"]) == ("not applicable", "no_y"), b)
        b = cn.assess(cat[loc], depths(cat[loc], 15.0, 15.0), None, "unknown")
        check(f"{loc} with no sex: not assessed", (b["status"], b["reason"]) == ("not assessed", "sex_unknown"), b)
    xp = cat["Xp22.31"]
    b = cn.assess(xp, depths(xp, 15.0, 15.0), None, "male")
    check("Xp22.31 on a male at his one X's depth: 1 copy, no change", (b["copies"], b["expected"]) == (1, 1), b)
    b = cn.assess(xp, depths(xp, 0.0, 15.0), None, "male")
    check("Xp22.31 on a male with no read there: 0 copies, a deletion (the STS deletion)", b["copies"] == 0, b)

    with tempfile.TemporaryDirectory() as t:
        fake = os.path.join(t, "fakecov.py")
        with open(fake, "w") as f:
            f.write(FAKECOV)
        bed = os.path.join(t, "r.bed")
        out = subprocess.run([sys.executable, SCRIPT, "bed", "--kind", "interval"], capture_output=True, text=True)
        open(bed, "w").write(out.stdout)
        check("bed --kind interval: every segment of every interval row, once",
              out.returncode == 0 and len(out.stdout.splitlines()) == sum(len(cn.regions(r)) for r in iv), out.stderr)
        halfs = ",".join(n for n, *_ in a_d["_targets"])
        env = dict(os.environ, FAKE_HALF=halfs, FAKE_ZERO=p17["_controls"][0][0])
        cov = os.path.join(t, "cov.tsv")
        with open(cov, "w") as f:
            subprocess.run([sys.executable, fake, bed], stdout=f, env=env, check=True)
        tab = os.path.join(t, "iv.tsv")
        out = subprocess.run([sys.executable, SCRIPT, "intervals", "--bedcov", cov, "--sex", "female", "--out", tab],
                             capture_output=True, text=True)
        lines = open(tab).read().splitlines() if os.path.exists(tab) else []
        rows = {x.split("\t")[0]: x.split("\t") for x in lines[1:]}
        check("intervals: exit 0, the header, one line per interval row",
              out.returncode == 0 and lines[:1] == ["\t".join(cn.INTERVAL_COLUMNS)] and list(rows) == names,
              (out.returncode, out.stderr[-400:], lines[:1]))
        check("intervals: the 22q11.2 A-D deletion at 1 copy", rows.get(a_d["locus"], [""] * 11)[5:9]
              == ["called", "1", "2", "deletion"], rows.get(a_d["locus"]))
        check("intervals: a control with no read leaves 17p12 not assessed (broken_control), the rest called",
              rows.get(p17["locus"], [""] * 11)[10] == "broken_control"
              and rows.get("3q29", [""] * 11)[5:9] == ["called", "2", "2", "no change"], rows.get(p17["locus"]))
        check("intervals: chrY rows on a female are not applicable",
              all(rows.get(x, [""] * 11)[5] == "not applicable" for x in ("AZFa", "AZFb", "AZFc_DAZ")))
        out = subprocess.run([sys.executable, SCRIPT, "intervals", "--failed", "--out", tab], capture_output=True, text=True)
        lines = open(tab).read().splitlines()[1:]
        check("intervals --failed: every row not assessed, check_failed",
              out.returncode == 0 and len(lines) == len(iv) and all(x.split("\t")[10] == "check_failed" for x in lines))

        step, module = read("scripts/35-paralogs.sh"), read("modules/local/parascopy/main.nf")
        blk = re.search(r'\n(INTERVAL_OUT=.*?\nrm -rf "\$INTERVAL_WORK"\n)', step, re.S)
        check("step 35 holds the interval block", blk is not None)
        env = dict(os.environ, REPO=REPO, FAKECOV=fake, OUTDIR=t, SAMPLE="S", PGP_ROOT=REPO, PARASCOPY_IMAGE="img",
                   BAM=os.path.join(t, "S.bam"), SEX="male", FAKE_HALF=",".join(n for n, *_ in az["_targets"]))
        out = subprocess.run(["bash", "-euo", "pipefail", "-c", STUB_IV + (blk.group(1) if blk else "false\n")],
                             capture_output=True, text=True, env=env)
        tab = os.path.join(t, "S_cn_intervals.tsv")
        rows = {x.split("\t")[0]: x.split("\t") for x in open(tab).read().splitlines()[1:]} if os.path.exists(tab) else {}
        check("step 35 interval block: exit 0, one line per interval row", out.returncode == 0 and list(rows) == names,
              (out.returncode, out.stderr[-600:]))
        # AZFc's DAZ at 15x against controls at 30x on a male: 4 x 15 / 30 = 2 copies of 4.
        check("step 35 interval block on a male: AZFc at half its depth is a deletion, AZFa one copy",
              rows.get("AZFc_DAZ", [""] * 11)[6:9] == ["2", "4", "deletion"]
              and rows.get("AZFa", [""] * 11)[6:9] == ["1", "1", "no change"], (rows.get("AZFc_DAZ"), rows.get("AZFa")))
        good = os.path.join(t, "good_cn_intervals.tsv")
        if os.path.exists(tab):
            os.replace(tab, good)
        env["FAKECOV"] = "/nonexistent.py"
        out = subprocess.run(["bash", "-euo", "pipefail", "-c", STUB_IV + (blk.group(1) if blk else "false\n")],
                             capture_output=True, text=True, env=env)
        lines = open(tab).read().splitlines()[1:] if os.path.exists(tab) else []
        check("step 35 interval block with a failed bedcov: exit 0, every row check_failed",
              out.returncode == 0 and len(lines) == len(iv) and all(x.split("\t")[10] == "check_failed" for x in lines),
              (out.returncode, out.stderr[-400:]))
        for path, text in (("scripts/35-paralogs.sh", step), ("modules/local/parascopy/main.nf", module)):
            for cmd in ("copy_number.py bed --kind interval", "copy_number.py intervals --bedcov",
                        "copy_number.py intervals --failed"):
                check(f"{path}: {cmd}", cmd in text)
        # collect_summary: step 35's section carries the interval rows with a change
        import collect_summary as cs
        sd = os.path.join(t, "sample")
        os.makedirs(os.path.join(sd, "paralogs"))
        with open(os.path.join(sd, "paralogs", "S_smn_copy_number.tsv"), "w") as f:
            f.write("chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\t"
                    "homologous_regions\n")
        if os.path.exists(good):
            os.replace(good, os.path.join(sd, "paralogs", "S_cn_intervals.tsv"))
        _, data = cs.sec_parascopy(sd, "S")
        ci = data.get("cn_intervals") or {}
        check("collect_summary: step 35's section lists the interval rows with a change (AZFc) and counts the rest",
              [c["locus"] for c in ci.get("changes", [])] == ["AZFc_DAZ"] and ci.get("rows") == len(iv), ci)
        os.remove(os.path.join(sd, "paralogs", "S_smn_copy_number.tsv"))
        _, data = cs.sec_parascopy(sd, "S")
        ci = (data or {}).get("cn_intervals") or {}
        check("collect_summary: the interval rows stay when step 35 wrote no other table",
              [c["locus"] for c in ci.get("changes", [])] == ["AZFc_DAZ"] and ci.get("rows") == len(iv), ci)
        check("PARASCOPY module: the samplesheet sex goes to the interval rows",
              "def sex = meta.sex in ['male', 'female'] ? meta.sex : 'unknown'" in module and "--sex ${sex}" in module)


if __name__ == "__main__":
    sys.exit(main())
