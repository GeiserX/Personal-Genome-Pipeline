#!/usr/bin/env python3
"""traits_report.py: the traits table and the genotype rules of step 42, read from the gVCF genotypes.

  traits_report.py --table assets/traits.tsv --rules assets/genotype_rules.tsv --genotypes G.tsv \\
      [--phased P.vcf] --sample S --out S_traits.tsv --rules-out S_genotype_rules.tsv
  traits_report.py --rules assets/genotype_rules.tsv --genotypes G.tsv --sample S --phase-sites IN.vcf
  traits_report.py --check-table assets/traits.tsv [--rules assets/genotype_rules.tsv]

G.tsv is what bcftools writes in step 42 (scripts/42-traits.sh and the
Nextflow TRAITS_SITES process run the same commands): one line per gVCF record
at a position of the traits table or of a rule site, after the reference
blocks are expanded to one record per base, as CHROM POS REF ALT GENOTYPE GQ,
the genotype as bases (C/T, or C|T when phased). A position with no line has
no genotype: it is "not_covered", never read as the reference. bcftools drops a
record with a missing allele (a no-call or a half call such as ./T) before
this script, so such a position has no line and reads as not_covered too; a
half call that still reaches here is other_allele.

Each traits table row becomes one output row with the genotype (the reference
allele first), a status (called, not_covered, other_allele) and the reading
for that genotype. A gVCF whose REF differs from the table's reference base at
a position is not GRCh38 (or the table is wrong): the run stops.

Each rule of the rules table (assets/genotype_rules.tsv) reads its sites
together and becomes one line of --rules-out: the two haplotypes the genotypes
make, named by the rule's haplotype lines, and the combination line of that
diplotype. A rule with needs_phase no takes the combination line of every
diplotype the unphased genotypes allow (APOE's double heterozygote is e2/e4
or e1/e3). A rule with needs_phase yes and two or more heterozygous sites
needs their phase: --phase-sites writes those sites as a VCF for WhatsHap,
and --phased reads WhatsHap's answer; sites not phased in one block leave the
rule "unphased", not read. The APOE rule's line also goes into the traits
table, right after its last site row, where the reports read it.

--check-table checks the schema of the table (and of the rules, with --rules)
and exits 1 on the first list of problems; tests/test_traits_table.py and
tests/test_traits_report.py run the same checks. Standard library only (it
runs in the plain python image and in the WhatsHap image).
"""
import argparse
import itertools
import os
import re
import sys

COLUMNS = ["chrom", "pos", "rsid", "gene", "ref", "alt", "trait", "hom_ref", "het", "hom_alt",
           "ancestry", "note", "pmids"]
OUT_COLUMNS = ["rsid", "gene", "chrom", "pos", "trait", "genotype", "status", "reading", "ancestry",
               "note", "pmids"]
RULE_OUT_COLUMNS = ["rule", "gene", "trait", "rsids", "chrom", "pos", "genotypes", "gq", "phase", "result",
                    "status", "reading", "ancestry", "note", "pmids"]
# The columns of each kind of line of the rules table, after the kind and the rule id.
RULE_KINDS = {"rule": ["gene", "trait", "needs_phase", "ancestry", "note", "pmids"],
              "site": ["rsid", "chrom", "pos", "ref", "alt"],
              "haplotype": ["name", "alleles"],
              "combination": ["diplotypes", "label", "reading"]}
CHROMS = {f"chr{n}" for n in range(1, 23)}
BASES = set("ACGT")

# APOE: rs429358 C (Arg112) and rs7412 C (Arg158) make e4; T and T make e2.
# The traits table must keep these two rows: the reports show them, and the
# APOE rule's line goes after them.
APOE_E4_SITE, APOE_E2_SITE = "rs429358", "rs7412"
APOE_ALLELES = {APOE_E4_SITE: ("T", "C"), APOE_E2_SITE: ("C", "T")}
# Rules whose line also goes into the traits table, where the reports read it
# (they fold the APOE lines away). Every other rule's line is in --rules-out.
TRAITS_TABLE_RULES = ("APOE",)
# The rules table beside the traits table in the repository, for report()'s callers
# that give none; the command line always takes --rules.
DEFAULT_RULES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "genotype_rules.tsv")
MAX_RULE_SITES = 8   # the completeness check walks 3^sites genotypes
NOT_COVERED = "No genotype at this position in the gVCF (no reads, or a no-call): not read."
OTHER_ALLELE = "An allele this table does not describe at this position: not read."


class TableError(Exception):
    pass


def read_table(path):
    """The rows of the traits table as dicts. Lines starting with # are comments;
    the first other line is the header, which must be COLUMNS."""
    rows, header = [], None
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cells = line.split("\t")
            if header is None:
                if cells != COLUMNS:
                    raise TableError(f"{path}:{n}: header is {cells}, expected {COLUMNS}")
                header = cells
                continue
            if len(cells) != len(COLUMNS):
                raise TableError(f"{path}:{n}: {len(cells)} columns, expected {len(COLUMNS)}")
            row = dict(zip(COLUMNS, cells))
            row["line"] = n
            rows.append(row)
    if header is None:
        raise TableError(f"{path}: no header line")
    return rows


def text_errors(where, text):
    errs = []
    if "snpedia" in text.lower():
        errs.append(f"{where}: mentions SNPedia; the table takes no text from it")
    if '"' in text or not text.isascii():
        errs.append(f"{where}: a double quote or a non-ASCII character (the reports read the table as plain TSV)")
    return errs


def site_errors(where, r):
    """Schema problems of one site: chrom, pos, rsid, ref and alt."""
    errs = []
    if r["chrom"] not in CHROMS:
        errs.append(f"{where}: chrom '{r['chrom']}' is not chr1 to chr22")
    if not re.fullmatch(r"[1-9][0-9]*", r["pos"]):
        errs.append(f"{where}: pos '{r['pos']}' is not a positive integer")
    if not re.fullmatch(r"rs[1-9][0-9]*", r["rsid"]):
        errs.append(f"{where}: rsid is not rs<number>")
    for k in ("ref", "alt"):
        if r[k] not in BASES:
            errs.append(f"{where}: {k} '{r[k]}' is not one of A, C, G, T")
    if r["ref"] == r["alt"]:
        errs.append(f"{where}: ref and alt are the same base")
    return errs


def check_rows(rows):
    """Every schema problem of the table rows, as messages (empty when it is valid)."""
    errs, ids, sites = [], set(), set()
    for r in rows:
        where = f"line {r['line']} ({r['rsid']})"
        errs += site_errors(where, r)
        if r["rsid"] in ids:
            errs.append(f"{where}: rsid appears twice")
        ids.add(r["rsid"])
        if (r["chrom"], r["pos"]) in sites:
            errs.append(f"{where}: position {r['chrom']}:{r['pos']} appears twice")
        sites.add((r["chrom"], r["pos"]))
        for k in ("gene", "trait", "hom_ref", "het", "hom_alt", "note"):
            if not r[k].strip():
                errs.append(f"{where}: {k} is empty")
        if r["ancestry"] not in ("yes", "no"):
            errs.append(f"{where}: ancestry '{r['ancestry']}' is not yes or no")
        if not re.fullmatch(r"[1-9][0-9]*(,[1-9][0-9]*)*", r["pmids"]):
            errs.append(f"{where}: pmids '{r['pmids']}' is not a comma-separated list of PubMed ids")
        errs += text_errors(where, " ".join(r[k] for k in COLUMNS))
    by_id = {r["rsid"]: r for r in rows}
    for rs, (ref, alt) in APOE_ALLELES.items():
        r = by_id.get(rs)
        if r is None:
            errs.append(f"APOE row {rs} is missing")
        elif (r["ref"], r["alt"]) != (ref, alt) or r["chrom"] != "chr19":
            errs.append(f"APOE row {rs} must be chr19 {ref}>{alt}, is {r['chrom']} {r['ref']}>{r['alt']}")
    return errs


# --- The rules table ---------------------------------------------------------------

def read_rules(path):
    """{rule id: {gene, trait, needs_phase, ancestry, note, pmids, line, sites, haplotypes,
    combinations}} in file order. Raises TableError on a line of an unknown kind or with
    the wrong column count, and on a site, haplotype or combination line before its rule line."""
    rules = {}
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cells = line.split("\t")
            kind = cells[0]
            if kind not in RULE_KINDS:
                raise TableError(f"{path}:{n}: kind '{kind}' is not one of {', '.join(RULE_KINDS)}")
            want = 2 + len(RULE_KINDS[kind])
            if len(cells) != want:
                raise TableError(f"{path}:{n}: a {kind} line has {len(cells)} columns, expected {want}")
            rid = cells[1]
            d = dict(zip(RULE_KINDS[kind], cells[2:]), line=n)
            if kind == "rule":
                if rid in rules:
                    raise TableError(f"{path}:{n}: rule {rid} has two rule lines")
                rules[rid] = dict(d, id=rid, sites=[], haplotypes=[], combinations=[])
                continue
            if rid not in rules:
                raise TableError(f"{path}:{n}: a {kind} line of rule '{rid}' before its rule line")
            rules[rid]["sites" if kind == "site" else kind + "s"].append(d)
    if not rules:
        raise TableError(f"{path}: no rule line")
    return rules


def diplotype_key(rule, a, b):
    """'a/b' with the two names in the order of the rule's haplotype lines."""
    order = [h["name"] for h in rule["haplotypes"]]
    return "/".join(sorted((a, b), key=order.index))


def allowed_diplotypes(rule, genotypes):
    """The diplotype keys, in haplotype order, whose two haplotypes give these unphased
    genotypes (one [allele, allele] per site)."""
    haps = rule["haplotypes"]
    out = []
    for i, j in itertools.combinations_with_replacement(range(len(haps)), 2):
        a, b = haps[i]["alleles"].split(","), haps[j]["alleles"].split(",")
        if all(sorted((a[k], b[k])) == sorted(g) for k, g in enumerate(genotypes)):
            out.append(f"{haps[i]['name']}/{haps[j]['name']}")
    return out


def check_rules(rules, table_rows=None):
    """Every schema problem of the rules, as messages (empty when they are valid). With
    table_rows, also: a rule site that is a traits table row has the row's chrom, pos, ref
    and alt, and the rules of TRAITS_TABLE_RULES exist with every site in the table."""
    errs = []
    by_rs = {r["rsid"]: r for r in table_rows or []}
    for rid, rule in rules.items():
        n0 = len(errs)
        where = f"rule {rid} (line {rule['line']})"
        if not re.fullmatch(r"[A-Za-z0-9_]+", rid):
            errs.append(f"{where}: the rule id is not letters, digits and _")
        for k in ("gene", "trait", "note"):
            if not rule[k].strip():
                errs.append(f"{where}: {k} is empty")
        for k in ("needs_phase", "ancestry"):
            if rule[k] not in ("yes", "no"):
                errs.append(f"{where}: {k} '{rule[k]}' is not yes or no")
        if not re.fullmatch(r"[1-9][0-9]*(,[1-9][0-9]*)*", rule["pmids"]):
            errs.append(f"{where}: pmids '{rule['pmids']}' is not a comma-separated list of PubMed ids")
        sites = rule["sites"]
        if not 2 <= len(sites) <= MAX_RULE_SITES:
            errs.append(f"{where}: {len(sites)} sites, expected 2 to {MAX_RULE_SITES}")
        seen = set()
        for s in sites:
            w = f"{where}, site line {s['line']} ({s['rsid']})"
            errs += site_errors(w, s)
            if s["rsid"] in seen or (s["chrom"], s["pos"]) in seen:
                errs.append(f"{w}: the site appears twice in the rule")
            seen |= {s["rsid"], (s["chrom"], s["pos"])}
            t = by_rs.get(s["rsid"])
            if t and (t["chrom"], t["pos"], t["ref"], t["alt"]) != (s["chrom"], s["pos"], s["ref"], s["alt"]):
                errs.append(f"{w}: the traits table has {t['chrom']}:{t['pos']} {t['ref']}>{t['alt']}, "
                            f"the rule {s['chrom']}:{s['pos']} {s['ref']}>{s['alt']}")
        if len({s["chrom"] for s in sites}) > 1:
            errs.append(f"{where}: the sites are on more than one chromosome")
        names, tuples = [], set()
        for h in rule["haplotypes"]:
            w = f"{where}, haplotype line {h['line']} ({h['name']})"
            if not re.fullmatch(r"[A-Za-z0-9*.]+", h["name"]):
                errs.append(f"{w}: the name is not letters, digits, * and .")
            if h["name"] in names:
                errs.append(f"{w}: the name appears twice")
            names.append(h["name"])
            al = h["alleles"].split(",")
            if len(al) != len(sites):
                errs.append(f"{w}: {len(al)} alleles for {len(sites)} sites")
                continue
            for a, s in zip(al, sites):
                if a not in (s["ref"], s["alt"]):
                    errs.append(f"{w}: allele {a} at {s['rsid']} is neither its ref {s['ref']} nor its alt {s['alt']}")
            if tuple(al) in tuples:
                errs.append(f"{w}: the same alleles as another haplotype")
            tuples.add(tuple(al))
        if len(names) < 2:
            errs.append(f"{where}: fewer than two haplotypes")
        keys = {}
        for c in rule["combinations"]:
            w = f"{where}, combination line {c['line']} ({c['diplotypes']})"
            for k in ("label", "reading"):
                if not c[k].strip():
                    errs.append(f"{w}: {k} is empty")
            errs += text_errors(w, " ".join((c["diplotypes"], c["label"], c["reading"])))
            if c["diplotypes"] in keys:
                errs.append(f"{w}: the diplotypes appear twice")
            keys[c["diplotypes"]] = c
            for d in c["diplotypes"].split("|"):
                pair = d.split("/")
                if len(pair) != 2 or any(p not in names for p in pair):
                    errs.append(f"{w}: '{d}' is not two of the rule's haplotypes joined by /")
                elif diplotype_key(rule, *pair) != d:
                    errs.append(f"{w}: '{d}' is not in the order of the haplotype lines "
                                f"(write {diplotype_key(rule, *pair)})")
        errs += text_errors(where, " ".join(str(rule[k]) for k in ("gene", "trait", "note")))
        if len(errs) > n0:
            continue
        # Completeness: the combination lines are exactly what the genotypes can give.
        if rule["needs_phase"] == "yes":
            want = {f"{a}/{b}" for a, b in itertools.combinations_with_replacement(names, 2)}
        else:
            want = set()
            choices = [[[s["ref"]] * 2, [s["ref"], s["alt"]], [s["alt"]] * 2] for s in sites]
            for g in itertools.product(*choices):
                key = "|".join(allowed_diplotypes(rule, g))
                if key:
                    want.add(key)
        for k in sorted(want - set(keys)):
            errs.append(f"{where}: no combination line for {k}")
        for k in sorted(set(keys) - want):
            errs.append(f"{where}: the combination line {k} can never be read "
                        + ("(a rule that needs phase reads one diplotype)" if "|" in k else "from these haplotypes"))
    if table_rows is not None:
        for rid in TRAITS_TABLE_RULES:
            if rid not in rules:
                errs.append(f"rule {rid} is missing: its line goes into the traits table, where the reports read it")
            else:
                for s in rules[rid]["sites"]:
                    if s["rsid"] not in by_rs:
                        errs.append(f"rule {rid}: site {s['rsid']} has no traits table row")
    return errs


# --- Genotypes ----------------------------------------------------------------------

def read_genotypes(path):
    """{(chrom, pos): [(REF, [allele, allele], GQ), ...]} from bcftools' CHROM POS REF ALT GT [GQ] lines."""
    recs = {}
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            c = line.split("\t")
            if len(c) < 5:
                raise TableError(f"{path}:{n}: {len(c)} columns, expected CHROM POS REF ALT GENOTYPE [GQ]")
            gq = c[5] if len(c) > 5 and c[5] else "."
            recs.setdefault((c[0], c[1]), []).append((c[2], re.split(r"[/|]", c[4]), gq))
    return recs


def site_genotype(site, recs):
    """(alleles sorted with the reference first, status, GQ) of one site from its gVCF records."""
    if not recs:
        return None, "not_covered", "."
    for r in recs:
        if r[0][:1] != site["ref"]:
            raise TableError(f"{site['rsid']} at {site['chrom']}:{site['pos']}: the gVCF has REF {r[0]}, the table "
                             f"{site['ref']} (the GRCh38 base). Is the gVCF on GRCh38?")
    snv = [r for r in recs if all(len(a) == 1 for a in r[1])]
    _, alleles, gq = (snv or recs)[0]
    order = {site["ref"]: 0, site["alt"]: 1}
    alleles = sorted(alleles, key=lambda a: (order.get(a, 2), a))
    if len(alleles) != 2 or any(a not in order for a in alleles):
        return alleles, "other_allele", gq
    return alleles, "called", gq


def call(row, recs):
    """(genotype, status, reading) of one traits table row from its gVCF records."""
    alleles, status, _ = site_genotype(row, recs)
    if status == "not_covered":
        return ".", "not_covered", NOT_COVERED
    gt = "/".join(alleles)
    if status == "other_allele":
        return gt, "other_allele", OTHER_ALLELE
    return gt, "called", row[("hom_ref", "het", "hom_alt")[alleles.count(row["alt"])]]


def read_phased(path):
    """{(chrom, pos): (allele of haplotype 1, allele of haplotype 2, phase block)} of the
    phased records of WhatsHap's VCF; unphased records are left out."""
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            fmt = dict(zip(c[8].split(":"), c[9].split(":")))
            gt, ps = fmt.get("GT", ""), fmt.get("PS", ".")
            if "|" not in gt or ps in (".", ""):
                continue
            bases = [c[3]] + c[4].split(",")
            try:
                a, b = (bases[int(i)] for i in gt.split("|"))
            except (ValueError, IndexError):
                continue
            out[(c[0], c[1])] = (a, b, ps)
    return out


def phase_sites(rules, recs):
    """The (chrom, pos, rsid, ref, alt) of the heterozygous sites of every rule that needs
    phase and has two or more of them, each site once, in rule order."""
    out, seen = [], set()
    for rule in rules.values():
        if rule["needs_phase"] != "yes":
            continue
        g = [(s, site_genotype(s, recs.get((s["chrom"], s["pos"]), []))) for s in rule["sites"]]
        if any(st != "called" for _, (_, st, _) in g):
            continue
        het = [s for s, (al, _, _) in g if al[0] != al[1]]
        if len(het) < 2:
            continue
        for s in het:
            k = (s["chrom"], s["pos"])
            if k not in seen:
                seen.add(k)
                out.append((s["chrom"], s["pos"], s["rsid"], s["ref"], s["alt"]))
    return out


def write_phase_vcf(path, sample, sites):
    """The sites as a one-sample VCF of 0/1 genotypes, sorted, for whatshap phase."""
    contigs = sorted({c for c, *_ in sites}, key=lambda c: int(c[3:]))
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in contigs:
            f.write(f"##contig=<ID={c}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
        for c, p, rs, ref, alt in sorted(sites, key=lambda s: (int(s[0][3:]), int(s[1]))):
            f.write(f"{c}\t{p}\t{rs}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/1\n")
    os.replace(tmp, path)


def rule_line(rule, recs, phased):
    """One output line of a rule. phased: read_phased's dict, or None when no phasing ran."""
    sites = rule["sites"]
    g = [site_genotype(s, recs.get((s["chrom"], s["pos"]), [])) for s in sites]
    rsids = [s["rsid"] for s in sites]
    line = {"rule": rule["id"], "gene": rule["gene"], "trait": rule["trait"], "rsids": ",".join(rsids),
            "chrom": sites[0]["chrom"], "pos": ",".join(s["pos"] for s in sites),
            "genotypes": ",".join("/".join(al) if al else "." for al, _, _ in g),
            "gq": ",".join(gq for _, _, gq in g), "ancestry": rule["ancestry"], "note": rule["note"],
            "pmids": rule["pmids"]}
    phase = "not_used" if rule["needs_phase"] == "no" else "not_needed"
    statuses = [st for _, st, _ in g]
    if any(st != "called" for st in statuses):
        need = (f"both {rsids[0]} and {rsids[1]}" if len(rsids) == 2
                else "every one of " + ", ".join(rsids[:-1]) + " and " + rsids[-1])
        return dict(line, phase=".", result=".",
                    status="not_covered" if "not_covered" in statuses else "other_allele",
                    reading=f"Not read: needs a genotype at {need}.")
    alleles = [al for al, _, _ in g]
    combos = {c["diplotypes"]: c for c in rule["combinations"]}
    if rule["needs_phase"] == "no":
        key = "|".join(allowed_diplotypes(rule, alleles))
    else:
        het = [i for i, al in enumerate(alleles) if al[0] != al[1]]
        h1, h2 = [al[0] for al in alleles], [al[1] for al in alleles]
        if len(het) >= 2:
            blocks = {(phased or {}).get((sites[i]["chrom"], sites[i]["pos"]), (None, None, None))[2] for i in het}
            if None in blocks or len(blocks) != 1:
                why = ("no read or read pair covered them together"
                       if phased is not None else "no phasing ran (no BAM, or phasing failed)")
                return dict(line, phase="unphased", result=".", status="unphased",
                            reading=f"Heterozygous at {len(het)} of the {len(sites)} sites, and {why}, so which "
                                    "alleles share a copy is not known: not read.")
            for i in het:
                h1[i], h2[i] = phased[(sites[i]["chrom"], sites[i]["pos"])][:2]
            phase = "read_backed"
        by_alleles = {tuple(h["alleles"].split(",")): h["name"] for h in rule["haplotypes"]}
        n1, n2 = by_alleles.get(tuple(h1)), by_alleles.get(tuple(h2))
        unknown = [",".join(h) for h, n in ((h1, n1), (h2, n2)) if n is None]
        if unknown:
            return dict(line, phase=phase, result=".", status="other_haplotype",
                        reading=f"A haplotype this rule does not describe ({' and '.join(unknown)} at "
                                f"{', '.join(rsids)}): not read.")
        key = diplotype_key(rule, n1, n2)
    if not key or key not in combos:
        return dict(line, phase=phase, result=".", status="other_haplotype",
                    reading="A combination of haplotypes this rule does not describe: not read.")
    c = combos[key]
    return dict(line, phase=phase, result=c["label"], status="called", reading=c["reading"])


def rule_lines(rules, recs, phased=None):
    """One output line per rule, in the order of the rules table."""
    return [rule_line(r, recs, phased) for r in rules.values()]


def report(rows, recs, rules=None, phased=None, lines=None):
    """The traits output rows: one per table row, with the line of each rule of
    TRAITS_TABLE_RULES right after its last site row. rules defaults to
    DEFAULT_RULES; lines, rule_lines' output, is computed when not given."""
    if rules is None:
        rules = read_rules(DEFAULT_RULES)
    if lines is None:
        lines = rule_lines(rules, recs, phased)
    after = {}
    for rid in TRAITS_TABLE_RULES:
        rule = rules.get(rid)
        if rule:
            members = {s["rsid"] for s in rule["sites"]}
            last = max((i for i, r in enumerate(rows) if r["rsid"] in members), default=None)
            if last is not None:
                after[last] = next(ln for ln in lines if ln["rule"] == rid)
    out = []
    for i, r in enumerate(rows):
        gt, status, reading = call(r, recs.get((r["chrom"], r["pos"]), []))
        out.append({"rsid": r["rsid"], "gene": r["gene"], "chrom": r["chrom"], "pos": r["pos"],
                    "trait": r["trait"], "genotype": gt, "status": status, "reading": reading,
                    "ancestry": r["ancestry"], "note": r["note"], "pmids": r["pmids"]})
        if i in after:
            ln = after[i]
            out.append({"rsid": ln["rsids"].replace(",", "+"), "gene": ln["gene"], "chrom": ln["chrom"],
                        "pos": ln["pos"], "trait": ln["trait"], "genotype": ln["result"], "status": ln["status"],
                        "reading": ln["reading"], "ancestry": ln["ancestry"], "note": ln["note"],
                        "pmids": ln["pmids"]})
    return out


def write_tsv(path, out, columns=OUT_COLUMNS):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for o in out:
            f.write("\t".join(str(o[k]) for k in columns) + "\n")
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", help="assets/traits.tsv")
    ap.add_argument("--rules", help="assets/genotype_rules.tsv")
    ap.add_argument("--genotypes", help="bcftools' CHROM POS REF ALT GENOTYPE GQ lines at the table and rule positions")
    ap.add_argument("--phased", help="WhatsHap's VCF of the --phase-sites sites (absent: no phasing ran)")
    ap.add_argument("--sample")
    ap.add_argument("--out", help="the sample's traits table (TSV)")
    ap.add_argument("--rules-out", help="the sample's genotype rules (TSV)")
    ap.add_argument("--phase-sites", metavar="VCF", help="write the sites WhatsHap must phase and exit")
    ap.add_argument("--check-table", metavar="TABLE", help="check the table's schema (and --rules') and exit")
    a = ap.parse_args(argv)
    try:
        if a.check_table:
            rows = read_table(a.check_table)
            errs = check_rows(rows)
            if a.rules:
                errs += check_rules(read_rules(a.rules), rows)
            for e in errs:
                print(f"ERROR: {e}", file=sys.stderr)
            if not errs:
                print(f"{a.check_table}{' and ' + a.rules if a.rules else ''}: schema OK")
            return 1 if errs else 0
        if a.phase_sites:
            if not (a.rules and a.genotypes and a.sample):
                ap.error("--phase-sites needs --rules, --genotypes and --sample")
            rules = read_rules(a.rules)
            errs = check_rules(rules)
            if errs:
                for e in errs:
                    print(f"ERROR: {a.rules}: {e}", file=sys.stderr)
                return 1
            sites = phase_sites(rules, read_genotypes(a.genotypes))
            write_phase_vcf(a.phase_sites, a.sample, sites)
            print(f"{a.sample}: {len(sites)} heterozygous rule sites to phase. Written: {a.phase_sites}")
            return 0
        if not (a.table and a.rules and a.genotypes and a.sample and a.out and a.rules_out):
            ap.error("give --table, --rules, --genotypes, --sample, --out and --rules-out "
                     "(or --phase-sites, or --check-table)")
        rows = read_table(a.table)
        rules = read_rules(a.rules)
        errs = [f"{a.table}: {e}" for e in check_rows(rows)] + [f"{a.rules}: {e}" for e in check_rules(rules, rows)]
        if errs:
            for e in errs:
                print(f"ERROR: {e}", file=sys.stderr)
            return 1
        phased = read_phased(a.phased) if a.phased else None
        recs = read_genotypes(a.genotypes)
        lines = rule_lines(rules, recs, phased)
        out = report(rows, recs, rules, phased, lines)
    except (TableError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    write_tsv(a.out, out)
    write_tsv(a.rules_out, lines, RULE_OUT_COLUMNS)
    counts = {}
    for o in out:
        counts[o["status"]] = counts.get(o["status"], 0) + 1
    print(f"{a.sample}: {len(rows)} table rows; " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items()))
          + f" (APOE line included). Written: {a.out}")
    print(f"{a.sample}: {len(lines)} genotype rules; "
          + ", ".join(f"{ln['rule']} {ln['result'] if ln['status'] == 'called' else ln['status']}" for ln in lines)
          + f". Written: {a.rules_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
