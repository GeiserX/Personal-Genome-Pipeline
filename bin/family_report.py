#!/usr/bin/env python3
"""family_report.py: the family analysis of step 38 (scripts/38-family.sh) and
of the Nextflow FAMILY workflow (workflows/family.nf). Both run the same three
commands, in this order. Standard library only, so it runs in the plain python
image.

  family_report.py prepare --pedigree P.tsv --gvcf ID=PATH ... [--declared-sex ID=SEX ...] --out-dir D
      Reads and checks the pedigree (assets/pedigree.example.tsv has the
      format) and checks that each gVCF names its own sample. Writes
      D/members.tsv (every person of the pedigree, sequenced or not),
      D/family.ped (somalier relate's pedigree: the sequenced members) and
      D/trios.ped (WhatsHap's: the children whose two parents are
      sequenced; empty when there is none). Exits 2 on a problem.

  family_report.py check --members D/members.tsv --pairs family.pairs.tsv --out relationships.tsv
      For every pair of sequenced members, compares the relatedness the
      pedigree implies with the one somalier relate measured. Exits 1 when
      a pair disagrees: a sample swap, or a wrong line in the pedigree.

  family_report.py report --members D/members.tsv --relationships R.tsv --vcf FAMILY.vcf.gz
                          [--phased FAMILY.phased.vcf.gz] [--candidates ID=VCF ...]
                          [--par-bed par_grch38.bed] --out-dir D
      Per child with a sequenced parent: the allele each parent passed on at
      every site of the joint VCF where the genotypes tell
      (<child>_inheritance.tsv.gz), de novo candidates when both parents
      are sequenced (de_novo.tsv), and pairs of candidate variants in one
      gene in trans or in cis (compound_het.tsv), and the chromosome map of
      bin/family_segments.py (<child>_segments.tsv, <child>_dnapainter.csv,
      <child>_segments.svg). family_report.md says it in words, pair by
      pair and child by child.
"""
import argparse
import gzip
import itertools
import os
import re
import sys

COLUMNS = ("sample", "sex", "father", "mother", "partner")
MISSING = {"", ".", "0", "-9"}
SEX = {"male": "male", "m": "male", "1": "male",
       "female": "female", "f": "female", "2": "female",
       "unknown": None, "": None, ".": None, "0": None, "-9": None}
NAME = re.compile(r"^[A-Za-z0-9._-]+$")

# somalier relate needs this many sites both samples have a genotype at
# before its relatedness is read as a verdict.
MIN_SITES = 1000
# Relatedness (twice the kinship coefficient) bins, as KING draws them: the
# geometric middle between two degrees.
BINS = ((0.707, "same person"), (0.354, "first degree"),
        (0.177, "second degree"), (0.0884, "third degree"))
# A parent and child share an allele at every site: the share of sites where
# they are opposite homozygotes (IBS0) stays near 0 (genotype errors only).
PARENT_CHILD_MAX_IBS0 = 0.01
# De novo candidates: every one of the three must have these.
MIN_GQ = 20
MIN_DP = 10
HET_AB = (0.25, 0.75)
PRIMARY = {"chr%s" % c for c in list(range(1, 23)) + ["X", "Y"]}


def die(lines, code=2):
    for line in lines:
        print("ERROR: " + line, file=sys.stderr)
    sys.exit(code)


def open_text(path):
    with open(path, "rb") as f:
        magic = f.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


class GzOut:
    """A gzip file whose bytes depend on its content only (no name, mtime 0)."""

    def __init__(self, path):
        self.raw = open(path, "wb")
        self.gz = gzip.GzipFile(filename="", mode="wb", fileobj=self.raw, mtime=0)

    def write(self, b):
        self.gz.write(b)

    def close(self):
        self.gz.close()
        self.raw.close()


def missing(v):
    return v is None or v.strip() in MISSING


# --- pedigree ---------------------------------------------------------------

def read_pedigree(path):
    """Rows of the pedigree TSV as dicts, and the problems found. The first
    non-blank line is the header (a leading '#' is dropped); later lines
    starting with '#' are comments."""
    rows, errors, header = [], [], None
    with open(path) as f:
        for n, raw in enumerate(f, 1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if header is None:
                header = [c.strip().lower() for c in line.lstrip("#").split("\t")]
                if "sample" not in header:
                    return [], ["%s: the first line must be a header with a 'sample' column "
                                "(columns: %s, tab separated); it is '%s'" % (path, ", ".join(COLUMNS), line)]
                bad = [c for c in header if c not in COLUMNS]
                if bad:
                    errors.append("%s: unknown column%s %s; the columns are %s"
                                  % (path, "s" if len(bad) > 1 else "", ", ".join(bad), ", ".join(COLUMNS)))
                continue
            if line.startswith("#"):
                continue
            vals = [v.strip() for v in line.split("\t")]
            if len(vals) > len(header):
                errors.append("%s line %d: %d fields, the header has %d (fields are tab separated)"
                              % (path, n, len(vals), len(header)))
                continue
            vals += [""] * (len(header) - len(vals))
            rec = dict(zip(header, vals))
            rec["line"] = n
            rows.append(rec)
    if header is None:
        errors.append("%s is empty" % path)
    return rows, errors


def build_members(rows, declared, path):
    """Every person the pedigree names: its rows, and the parents and
    partners it names without a row of their own (not sequenced, no
    parents known). Returns (members, errors); a member is a dict with
    sex, father, mother, partner (None when not given) and row (bool)."""
    errors, members = [], {}

    def bad_name(v):
        return not NAME.match(v) or v.startswith("-")

    for r in rows:
        where = "%s line %d" % (path, r["line"])
        s = r.get("sample", "").strip()
        if missing(s):
            errors.append("%s: no sample name" % where)
            continue
        if bad_name(s):
            errors.append("%s: sample name '%s' may hold only letters, digits, '.', '_' and '-' "
                          "and may not start with '-'" % (where, s))
            continue
        if s in members:
            errors.append("%s: sample '%s' has a second row; each person has one" % (where, s))
            continue
        sex_raw = r.get("sex", "").strip().lower()
        if sex_raw not in SEX:
            errors.append("%s: sex '%s' of '%s' is not male, female or unknown" % (where, r.get("sex"), s))
            sex_raw = ""
        m = {"sex": SEX[sex_raw], "row": True}
        for k in ("father", "mother", "partner"):
            v = r.get(k, "").strip()
            if missing(v):
                m[k] = None
            elif bad_name(v):
                errors.append("%s: %s '%s' of '%s' is not a valid sample name" % (where, k, v, s))
                m[k] = None
            elif v == s:
                errors.append("%s: '%s' is named as its own %s" % (where, s, k))
                m[k] = None
            else:
                m[k] = v
        if m["father"] and m["father"] == m["mother"]:
            errors.append("%s: '%s' has the same person, '%s', as father and mother" % (where, s, m["father"]))
        d = declared.get(s)
        if d:
            if m["sex"] and m["sex"] != d:
                errors.append("%s: the pedigree says '%s' is %s, the samplesheet says %s"
                              % (where, s, m["sex"], d))
            m["sex"] = m["sex"] or d
        members[s] = m

    # People named but without a row: parents and partners not in the file.
    for s, m in list(members.items()):
        for k, sex in (("father", "male"), ("mother", "female"), ("partner", None)):
            p = m.get(k)
            if p and p not in members:
                members[p] = {"sex": sex, "father": None, "mother": None, "partner": None, "row": False}

    for s, m in members.items():
        if m["father"] and members[m["father"]]["sex"] == "female":
            errors.append("'%s' is the father of '%s' but its sex is female" % (m["father"], s))
        if m["mother"] and members[m["mother"]]["sex"] == "male":
            errors.append("'%s' is the mother of '%s' but its sex is male" % (m["mother"], s))
    fathers = {m["father"] for m in members.values() if m["father"]}
    mothers = {m["mother"] for m in members.values() if m["mother"]}
    for p in sorted(fathers & mothers):
        errors.append("'%s' is named as a father and as a mother" % p)

    for s, m in members.items():
        p = m["partner"]
        if p and members[p]["row"] and members[p]["partner"] not in (None, s):
            errors.append("'%s' names '%s' as partner, but '%s' names '%s'"
                          % (s, p, p, members[p]["partner"]))

    # No one is their own ancestor.
    for s in members:
        seen, todo = set(), [s]
        while todo:
            x = todo.pop()
            for k in ("father", "mother"):
                p = members.get(x, {}).get(k)
                if not p:
                    continue
                if p == s:
                    errors.append("'%s' is among its own ancestors: the parent columns form a loop" % s)
                    todo = []
                    break
                if p not in seen:
                    seen.add(p)
                    todo.append(p)
    return members, errors


def gvcf_samples(path):
    """Sample columns of a VCF's #CHROM line (reads the header only)."""
    with open_text(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\r\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    return None


def kinship_table(members):
    """phi(a, b): the kinship coefficient of two members, from the pedigree.
    Founders (no parent known) are unrelated; an unknown parent is a founder."""
    depth_memo, phi_memo = {}, {}

    def depth(x):
        if x not in depth_memo:
            m = members[x]
            ps = [p for p in (m["father"], m["mother"]) if p]
            depth_memo[x] = 1 + max(depth(p) for p in ps) if ps else 0
        return depth_memo[x]

    def phi(a, b):
        if a == b:
            m = members[a]
            if m["father"] and m["mother"]:
                return 0.5 * (1 + phi(m["father"], m["mother"]))
            return 0.5
        key = (a, b) if a < b else (b, a)
        if key not in phi_memo:
            if depth(a) < depth(b):
                a, b = b, a
            m = members[a]
            if not m["father"] and not m["mother"]:
                phi_memo[key] = 0.0
            else:
                phi_memo[key] = 0.5 * sum(phi(p, b) for p in (m["father"], m["mother"]) if p)
        return phi_memo[key]

    return phi


def declared_relationship(a, b, members, expected):
    ma, mb = members[a], members[b]
    if a in (mb["father"], mb["mother"]) or b in (ma["father"], ma["mother"]):
        return "parent-child"
    shared = {p for p in (ma["father"], ma["mother"]) if p} & {p for p in (mb["father"], mb["mother"]) if p}
    if len(shared) == 2:
        return "full siblings"
    if len(shared) == 1:
        return "half siblings"
    if ma["partner"] == b or mb["partner"] == a:
        return "partners" if expected == 0 else "partners, %s" % degree(expected)
    return degree(expected) if expected > 0 else "unrelated"


def degree(r):
    for lo, name in BINS:
        if r >= lo:
            return name
    return "unrelated"


def read_members(path):
    members = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            members[r["sample"]] = {
                "sex": r["sex"] or None, "father": r["father"] or None, "mother": r["mother"] or None,
                "partner": r["partner"] or None, "row": r["in_pedigree"] == "yes",
                "sequenced": r["sequenced"] == "yes"}
    return members


def cmd_prepare(a):
    declared = {}
    for kv in a.declared_sex or []:
        k, _, v = kv.partition("=")
        v = v.strip().lower()
        if v:
            if v not in ("male", "female"):
                die(["--declared-sex %s: the sex is male or female" % kv])
            declared[k] = v
    rows, errors = read_pedigree(a.pedigree)
    members, more = build_members(rows, declared, a.pedigree)
    errors += more
    gvcfs = {}
    for kv in a.gvcf or []:
        k, _, v = kv.partition("=")
        gvcfs[k] = v
    seq = sorted(s for s in gvcfs if s in members and members[s]["row"])
    if len(seq) < 2:
        errors.append("the pedigree names %d of the processed samples (%s); the family analysis needs two "
                      "or more, each with a row in %s and a gVCF"
                      % (len(seq), ", ".join(seq) or "none", a.pedigree))
    for s in seq:
        names = gvcf_samples(gvcfs[s])
        if names != [s]:
            errors.append("the gVCF of '%s' (%s) holds the sample%s %s, not '%s'. The joint call names each "
                          "column after the gVCF's own sample. Rename it: echo %s > name.txt; bcftools reheader "
                          "-s name.txt -o %s.g.vcf.gz %s; bcftools index -t %s.g.vcf.gz"
                          % (s, os.path.basename(gvcfs[s]), "s" if names and len(names) > 1 else "",
                             ", ".join(names) if names else "none (no #CHROM line)", s, s, s,
                             os.path.basename(gvcfs[s]), s))
    if errors:
        die(errors)
    for s in seq:
        members[s]["sequenced"] = True
    phi = kinship_table(members)
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "members.tsv"), "w") as f:
        f.write("sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced\n")
        for s in sorted(members):
            m = members[s]
            f.write("\t".join([s, m["sex"] or "", m["father"] or "", m["mother"] or "", m["partner"] or "",
                               "yes" if m["row"] else "no", "yes" if m.get("sequenced") else "no"]) + "\n")
    code = {"male": "1", "female": "2"}
    with open(os.path.join(a.out_dir, "family.ped"), "w") as f:
        for s in seq:
            m = members[s]
            f.write("\t".join(["family", s, m["father"] if m["father"] in seq else "-9",
                               m["mother"] if m["mother"] in seq else "-9",
                               code.get(m["sex"], "-9"), "-9"]) + "\n")
    trios = [s for s in seq if members[s]["father"] in seq and members[s]["mother"] in seq]
    with open(os.path.join(a.out_dir, "trios.ped"), "w") as f:
        for s in trios:
            m = members[s]
            f.write("\t".join(["family", s, m["father"], m["mother"], code.get(m["sex"], "0"), "0"]) + "\n")
    print("Family: %d people in the pedigree, %d sequenced (%s)" % (len(members), len(seq), ", ".join(seq)))
    for x, y in itertools.combinations(seq, 2):
        e = 2 * phi(x, y)
        print("  %s and %s: declared %s (expected relatedness %.3f)"
              % (x, y, declared_relationship(x, y, members, e), e))
    for s in seq:
        ps = [p for p in (members[s]["father"], members[s]["mother"]) if p in seq]
        if ps:
            print("  child %s: %s sequenced" % (s, " and ".join(ps)))
    print("Trios WhatsHap phases: %s" % (", ".join(trios) or "none (no child has both parents sequenced)"))


# --- relationship check -------------------------------------------------------

def read_pairs(path):
    pairs = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").lstrip("#").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            pairs[frozenset((r["sample_a"], r["sample_b"]))] = r
    return pairs


def verdict(declared, expected, rel, ibs0_rate, n):
    """(verdict, reason) for one pair."""
    if n < MIN_SITES:
        return "not_checked", ("somalier compared the two at %d sites, fewer than %d: too few to tell "
                               "(a sliced or targeted sample, or the wrong sites file)" % (n, MIN_SITES))
    want, got = degree(expected), degree(rel)
    if expected >= 0.177:
        if got != want:
            return "mismatch", "declared %s (%s), but somalier measures %s" % (declared, want, got)
        if declared == "parent-child" and ibs0_rate > PARENT_CHILD_MAX_IBS0:
            return "mismatch", ("first degree, but opposite homozygotes at %.1f%% of sites: a parent and "
                                "child share an allele at every site (siblings do not)" % (100 * ibs0_rate))
        return "ok", "%s, as declared" % got
    if rel >= 0.354:
        return "mismatch", "declared %s, but somalier measures %s" % (declared, got)
    if got != want and rel > expected:
        return "warn", "declared %s, somalier measures %s: closer than the pedigree says" % (declared, got)
    return "ok", "%s, as declared" % got


def cmd_check(a):
    members = read_members(a.members)
    pairs = read_pairs(a.pairs)
    phi = kinship_table(members)
    seq = sorted(s for s, m in members.items() if m["sequenced"])
    out, bad, errors = [], [], []
    for x, y in itertools.combinations(seq, 2):
        r = pairs.get(frozenset((x, y)))
        if r is None:
            errors.append("somalier relate has no row for %s and %s; its samples are: %s"
                          % (x, y, ", ".join(sorted({s for p in pairs for s in p}))))
            continue
        e = 2 * phi(x, y)
        d = declared_relationship(x, y, members, e)
        rel, ibs0, ibs2, n = float(r["relatedness"]), int(r["ibs0"]), int(r["ibs2"]), int(r["n"])
        rate = ibs0 / n if n else 0.0
        v, why = verdict(d, e, rel, rate, n)
        out.append([x, y, d, "%.3f" % e, "%.3f" % rel, str(ibs0), "%.4f" % rate, str(ibs2), str(n), v, why])
        if v == "mismatch":
            bad.append("%s and %s: %s" % (x, y, why))
        elif v != "ok":
            print("WARNING: %s and %s: %s" % (x, y, why), file=sys.stderr)
    if errors:
        die(errors)
    with open(a.out, "w") as f:
        f.write("sample_a\tsample_b\tdeclared\texpected_relatedness\trelatedness\tibs0\tibs0_rate\tibs2\tsites"
                "\tverdict\treason\n")
        for row in out:
            f.write("\t".join(row) + "\n")
    for row in out:
        print("%s and %s: %s; relatedness %s (expected %s), %s sites: %s"
              % (row[0], row[1], row[2], row[4], row[3], row[8], row[9]))
    if bad:
        die(["the data disagree with the pedigree:"] + ["  " + b for b in bad] +
            ["A sample swap (at the lab or in the samplesheet), or a wrong line in the pedigree. "
             "Check the pedigree and the samplesheet; somalier's own tables are beside relationships.tsv."],
            code=1)


# --- report -------------------------------------------------------------------

def read_par(path):
    par = {}
    if path:
        with open(path) as f:
            for line in f:
                c = line.split("\t")
                if len(c) >= 3:
                    par.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    return par


def in_par(par, chrom, pos):
    # BED: 0-based start, end exclusive; pos is 1-based
    return any(s < pos <= e for s, e in par.get(chrom, ()))


def parse_gt(s):
    """(alleles, phased) of a GT string; alleles is None when any is missing."""
    if not s or s == ".":
        return None, False
    phased = "|" in s
    parts = re.split(r"[/|]", s)
    if any(p == "." or p == "" for p in parts):
        return None, phased
    return tuple(int(p) for p in parts), phased


def fmt_field(fmt, vals, key):
    try:
        v = vals[fmt.index(key)]
    except (ValueError, IndexError):
        return None
    return None if v in (".", "") else v


def to_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def candidate_genes(path):
    """{(chrom, pos, ref, alt): genes} of a candidate VCF: VEP's CSQ SYMBOL,
    else snpEff's ANN gene, else INFO/GENE."""
    out, symbol, ann = {}, None, False
    with open_text(path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r"Format: ([^\"]+)", line)
                if m:
                    fields = m.group(1).split("|")
                    if "SYMBOL" in fields:
                        symbol = fields.index("SYMBOL")
                continue
            if line.startswith("##INFO=<ID=ANN"):
                ann = True
                continue
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[7].split(";"))
            genes = set()
            if symbol is not None and "CSQ" in info:
                for entry in info["CSQ"].split(","):
                    f2 = entry.split("|")
                    if len(f2) > symbol and f2[symbol]:
                        genes.add(f2[symbol])
            elif ann and "ANN" in info:
                for entry in info["ANN"].split(","):
                    f2 = entry.split("|")
                    if len(f2) > 3 and f2[3]:
                        genes.add(f2[3])
            elif info.get("GENE"):
                genes.update(g for g in info["GENE"].split(",") if g)
            if genes:
                for alt in c[4].split(","):
                    out[(c[0], int(c[1]), c[3], alt)] = genes
    return out


def resolve(child, f_alleles, m_alleles):
    """The (from_father, from_mother) allele assignments a diploid child
    genotype allows; a parent not sequenced (None) allows any allele."""
    a, b = child
    opts = set()
    for x, y in ((a, b), (b, a)):
        if (f_alleles is None or x in f_alleles) and (m_alleles is None or y in m_alleles):
            opts.add((x, y))
    return opts


class Child:
    def __init__(self, name, m, members, out_dir):
        self.name = name
        self.sex = m["sex"]
        self.father, self.mother = m["father"], m["mother"]
        self.f_seq = bool(self.father and members[self.father]["sequenced"])
        self.m_seq = bool(self.mother and members[self.mother]["sequenced"])
        self.path = os.path.join(out_dir, "%s_inheritance.tsv.gz" % name)
        self.fh = GzOut(self.path)
        self.fh.write(b"chrom\tpos\tref\talt\tchild_gt\tfather_gt\tmother_gt\tfrom_father\tfrom_mother\thow\n")
        self.n = {"sites": 0, "father": 0, "mother": 0, "both": 0, "ambiguous": 0, "inconsistent": 0,
                  "skipped_sex": 0}
        self.de_novo = []
        self.cands = []
        # WhatsHap: per (chrom, PS) block, sites where the first haplotype
        # holds the father's allele, and sites where it holds the mother's
        self.blocks = {}
        self.phased_hets = 0


def cmd_report(a):
    members = read_members(a.members)
    rels = []
    with open(a.relationships) as f:
        head = f.readline().rstrip("\n").split("\t")
        rels = [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f]
    par = read_par(a.par_bed)
    os.makedirs(a.out_dir, exist_ok=True)
    seq = sorted(s for s, m in members.items() if m["sequenced"])
    children = [Child(s, members[s], members, a.out_dir) for s in seq
                if any(p and members[p]["sequenced"] for p in (members[s]["father"], members[s]["mother"]))]
    cands = {}
    for kv in a.candidates or []:
        k, _, v = kv.partition("=")
        if k in members:
            cands[k] = candidate_genes(v)

    vcf = a.phased or a.vcf
    samples, idx = None, {}
    with open_text(vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                idx = {s: 9 + i for i, s in enumerate(samples)}
                lost = [s for s in seq if s not in idx]
                if lost:
                    die(["%s has no column for %s (its columns: %s)" % (vcf, ", ".join(lost), ", ".join(samples))])
                continue
            c = line.rstrip("\n").split("\t")
            chrom, pos = c[0], int(c[1])
            if chrom not in PRIMARY:
                continue
            alleles = [c[3]] + c[4].split(",")
            fmt = c[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else None
            if gi is None:
                continue
            calls = {}

            def call(s):
                if s not in calls:
                    vals = c[idx[s]].split(":")
                    calls[s] = (parse_gt(vals[gi] if gi < len(vals) else "."), vals)
                return calls[s]

            sexchr = chrom in ("chrX", "chrY") and not in_par(par, chrom, pos)
            for ch in children:
                (cg, cphased), cvals = call(ch.name)
                if cg is None:
                    continue
                fg = call(ch.father)[0][0] if ch.f_seq else None
                mg = call(ch.mother)[0][0] if ch.m_seq else None
                fset = set(fg) if fg is not None else None
                mset = set(mg) if mg is not None else None
                # a sequenced parent with no call here tells nothing
                f_known, m_known = fset is not None, mset is not None
                if not f_known and not m_known:
                    continue
                how, frm_f, frm_m = None, None, None
                if sexchr:
                    if ch.sex not in ("male", "female") or (chrom == "chrY" and ch.sex == "female"):
                        ch.n["skipped_sex"] += 1
                        continue
                    if ch.sex == "male":
                        if len(set(cg)) != 1:
                            ch.n["skipped_sex"] += 1   # a heterozygous call on a single X or Y
                            continue
                        allele = cg[0]
                        src = mset if chrom == "chrX" else fset
                        if (chrom == "chrX" and not m_known) or (chrom == "chrY" and not f_known):
                            continue
                        if allele in src:
                            how = "genotype"
                            if chrom == "chrX":
                                frm_m = allele
                            else:
                                frm_f = allele
                        else:
                            how = "inconsistent"
                if how is None:
                    if len(cg) != 2:
                        continue
                    opts = resolve(cg, fset, mset)
                    if not opts:
                        how = "inconsistent"
                    elif len(opts) == 1:
                        how = "genotype"
                        frm_f, frm_m = next(iter(opts))
                    else:
                        how = "ambiguous"
                ch.n["sites"] += 1
                if how == "genotype":
                    if frm_f is not None and frm_m is not None:
                        ch.n["both"] += 1
                    if frm_f is not None:
                        ch.n["father"] += 1
                    if frm_m is not None:
                        ch.n["mother"] += 1
                else:
                    ch.n[how] += 1
                gt = lambda s, known: (call(s)[1][gi] if gi < len(call(s)[1]) else ".") if known else "."
                ch.fh.write(("\t".join([chrom, c[1], c[3], c[4], cvals[gi], gt(ch.father, ch.f_seq),
                                        gt(ch.mother, ch.m_seq),
                                        alleles[frm_f] if frm_f is not None else ".",
                                        alleles[frm_m] if frm_m is not None else ".", how]) + "\n").encode())
                het = len(set(cg)) == 2 and 0 in cg
                # WhatsHap's phase of the child against the genotype origin
                if a.phased and cphased and het:
                    ch.phased_hets += 1
                    ps = fmt_field(fmt, cvals, "PS") or "."
                    if how == "genotype" and frm_f is not None:
                        b = ch.blocks.setdefault((chrom, ps), [0, 0])
                        b[0 if cg[0] == frm_f else 1] += 1
                # de novo: an allele the child carries that neither parent has
                if how == "inconsistent" and f_known and m_known:
                    # on a male's single X only the mother's alleles count, on his Y the father's
                    seen = mset if sexchr and chrom == "chrX" else fset if sexchr else fset | mset
                    novel = sorted(set(cg) - seen - {0})
                    if novel:
                        ch.de_novo.append(de_novo_row(ch, c, fmt, cvals, call, gi, alleles, novel, sexchr))
                # compound heterozygote candidates
                if het and ch.name in cands:
                    alt = alleles[cg[1] if cg[0] == 0 else cg[0]]
                    genes = cands[ch.name].get((chrom, pos, c[3], alt))
                    if genes:
                        ai = alleles.index(alt)
                        origin = "father" if frm_f == ai else "mother" if frm_m == ai else None
                        ch.cands.append((genes, "%s:%d:%s>%s" % (chrom, pos, c[3], alt), origin))
    if samples is None:
        die(["%s has no #CHROM line" % vcf])

    with open(os.path.join(a.out_dir, "de_novo.tsv"), "w") as f:
        f.write("child\tchrom\tpos\tref\talt\tchild_gt\tfather_gt\tmother_gt\tchild_dp\tchild_gq\tchild_ab"
                "\tfather_dp\tfather_gq\tmother_dp\tmother_gq\tquality\n")
        for ch in children:
            for row in ch.de_novo:
                f.write("\t".join(row) + "\n")
    comp = []
    for ch in children:
        by_gene = {}
        for genes, var, origin in ch.cands:
            for g in genes:
                by_gene.setdefault(g, []).append((var, origin))
        for g in sorted(by_gene):
            for (v1, o1), (v2, o2) in itertools.combinations(by_gene[g], 2):
                phase = "unknown" if not (o1 and o2) else ("in cis" if o1 == o2 else "in trans")
                comp.append([ch.name, g, v1, v2, o1 or ".", o2 or ".", phase])
    with open(os.path.join(a.out_dir, "compound_het.tsv"), "w") as f:
        f.write("child\tgene\tvariant_a\tvariant_b\ta_from\tb_from\tphase\n")
        for row in comp:
            f.write("\t".join(row) + "\n")
    for ch in children:
        ch.fh.close()
    # The chromosome map: which stretches each child got from each parent
    # (bin/family_segments.py), from the joint VCF's genotypes
    import family_segments
    seg = family_segments.run(members, a.vcf, par, a.out_dir)
    write_markdown(a, members, rels, children, comp, cands, seg)


def de_novo_row(ch, c, fmt, cvals, call, gi, alleles, novel, sexchr):
    fv, mv = call(ch.father)[1], call(ch.mother)[1]
    dp = [to_int(fmt_field(fmt, v, "DP")) for v in (cvals, fv, mv)]
    gq = [to_int(fmt_field(fmt, v, "GQ")) for v in (cvals, fv, mv)]
    ab = None
    ad = fmt_field(fmt, cvals, "AD")
    if ad:
        n = [to_int(x) or 0 for x in ad.split(",")]
        tot = sum(n)
        if tot:
            ab = n[novel[0]] / tot if novel[0] < len(n) else None
    het = len(set(parse_gt(cvals[gi])[0])) == 2
    good = all(x is not None and x >= MIN_GQ for x in gq) and all(x is not None and x >= MIN_DP for x in dp)
    if het and (ab is None or not HET_AB[0] <= ab <= HET_AB[1]):
        good = False
    s = lambda x: "." if x is None else str(x)
    return [ch.name, c[0], c[1], c[3], ",".join(alleles[i] for i in novel), cvals[gi],
            fv[gi] if gi < len(fv) else ".", mv[gi] if gi < len(mv) else ".",
            s(dp[0]), s(gq[0]), "." if ab is None else "%.2f" % ab, s(dp[1]), s(gq[1]), s(dp[2]), s(gq[2]),
            "pass" if good else "low"]


def pct(k, n):
    return "%.1f%%" % (100.0 * k / n) if n else "-"


def write_markdown(a, members, rels, children, comp, cands, seg=()):
    L = ["# Family analysis", ""]
    seq = sorted(s for s, m in members.items() if m["sequenced"])
    others = sorted(s for s, m in members.items() if not m["sequenced"])
    L.append("Sequenced: %s. In the pedigree, not sequenced: %s." % (", ".join(seq), ", ".join(others) or "none"))
    L += ["", "## Relationships", "",
          "Relatedness is twice the kinship coefficient: about the share of the genome the two have "
          "inherited from common ancestors (1 the same person, 0.5 parent and child or full siblings, "
          "0.25 second degree, 0 unrelated). somalier measures it at its sites; the pedigree gives the "
          "expected value.", "",
          "| Pair | Declared | Expected | Measured | Opposite homozygotes | Sites | Verdict |",
          "|---|---|---|---|---|---|---|"]
    for r in rels:
        L.append("| %s and %s | %s | %s | %s | %s | %s | %s: %s |"
                 % (r["sample_a"], r["sample_b"], r["declared"], r["expected_relatedness"], r["relatedness"],
                    pct(float(r["ibs0_rate"]), 1), r["sites"], r["verdict"], r["reason"]))
    for ch in children:
        n = ch.n
        L += ["", "## %s" % ch.name, ""]
        fa = "%s (%s)" % (ch.father, "sequenced" if ch.f_seq else "not sequenced") if ch.father else "not named"
        mo = "%s (%s)" % (ch.mother, "sequenced" if ch.m_seq else "not sequenced") if ch.mother else "not named"
        L.append("Father: %s. Mother: %s." % (fa, mo))
        L += ["", "Sites of the joint VCF where %s has a genotype and a sequenced parent does too: %d. "
              "`%s` lists them." % (ch.name, n["sites"], os.path.basename(ch.path)), ""]
        L.append("- the allele from the father is known at %d (%s)" % (n["father"], pct(n["father"], n["sites"])))
        L.append("- the allele from the mother is known at %d (%s)" % (n["mother"], pct(n["mother"], n["sites"])))
        L.append("- both parents heterozygous like the child, so no origin from the genotypes: %d (%s)"
                 % (n["ambiguous"], pct(n["ambiguous"], n["sites"])))
        L.append("- inconsistent with the declared parents (a genotype error, a deletion, or a new "
                 "mutation): %d (%s)" % (n["inconsistent"], pct(n["inconsistent"], n["sites"])))
        if n["skipped_sex"]:
            L.append("- left out on chrX or chrY outside the pseudoautosomal regions (sex unknown, or a "
                     "heterozygous call on a single X or Y): %d" % n["skipped_sex"])
        for who, seq_, k in (("father", ch.f_seq, "father"), ("mother", ch.m_seq, "mother")):
            if not seq_:
                L += ["", "The %s is not sequenced. At %d of these sites (%s) the genotypes give the allele "
                      "the %s passed on to %s: that is the half of the %s's genome %s carries. The other "
                      "half was not passed on and is not seen here."
                      % (who, n[k], pct(n[k], n["sites"]), who, ch.name, who, ch.name)]
        if ch.f_seq and ch.m_seq:
            good = sum(1 for r in ch.de_novo if r[-1] == "pass")
            L += ["", "De novo candidates (an allele neither parent has): %d, %d of them with GQ >= %d and "
                  "depth >= %d in all three and, for a heterozygous call, an allele balance of %.2f to %.2f "
                  "(`de_novo.tsv`). Most candidates of a genome are genotype errors; a candidate that "
                  "matters needs confirming by another method."
                  % (len(ch.de_novo), good, MIN_GQ, MIN_DP, HET_AB[0], HET_AB[1])]
        else:
            L += ["", "No de novo candidates are called: that needs both parents sequenced."]
        if a.phased:
            votes = sum(max(b) for b in ch.blocks.values())
            tot = sum(sum(b) for b in ch.blocks.values())
            L += ["", "WhatsHap phased %d heterozygous sites of %s by inheritance. At %d of them the "
                  "genotypes also give the origin, and the phase agrees with it at %d (%s)."
                  % (ch.phased_hets, ch.name, tot, votes, pct(votes, tot))]
    L += ["", "## Compound heterozygotes", ""]
    if not cands:
        L.append("No candidate variants were given (the clinical filter, step 23, did not run for a child), "
                 "so no pair is reported.")
    else:
        t = sum(1 for r in comp if r[-1] == "in trans")
        c_ = sum(1 for r in comp if r[-1] == "in cis")
        L.append("Pairs of the clinical filter's candidate variants a child carries heterozygous in one gene: "
                 "%d; in trans (one from each parent, both copies of the gene hit) %d, in cis (both from one "
                 "parent) %d, origin unknown %d (`compound_het.tsv`)." % (len(comp), t, c_, len(comp) - t - c_))
    L += list(seg)
    with open(os.path.join(a.out_dir, "family_report.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("prepare")
    q.add_argument("--pedigree", required=True)
    q.add_argument("--gvcf", action="append", metavar="ID=PATH")
    q.add_argument("--declared-sex", action="append", metavar="ID=SEX")
    q.add_argument("--out-dir", required=True)
    q = sub.add_parser("check")
    q.add_argument("--members", required=True)
    q.add_argument("--pairs", required=True)
    q.add_argument("--out", required=True)
    q = sub.add_parser("report")
    q.add_argument("--members", required=True)
    q.add_argument("--relationships", required=True)
    q.add_argument("--vcf", required=True)
    q.add_argument("--phased")
    q.add_argument("--candidates", action="append", metavar="ID=VCF")
    q.add_argument("--par-bed")
    q.add_argument("--out-dir", required=True)
    a = p.parse_args(argv)
    {"prepare": cmd_prepare, "check": cmd_check, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
