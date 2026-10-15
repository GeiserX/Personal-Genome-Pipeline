#!/usr/bin/env python3
"""export_consumer.py: step 40, a consumer raw-data file from the gVCF.

Genealogy and health services that match relatives or read health traits take
the text file a genotyping chip company gives its customers. This script
writes that file from the pipeline's own calls, for the sites of the 23andMe
v5 chip, so the genome can be used there. It writes a file and nothing else:
uploading it is the user's decision (docs/40-export-consumer.md).

  export_consumer.py build-sites --template T.txt.gz --dbsnp VCF|- --out SITES.tsv.gz
                                 --dropped DROPPED.tsv [--source TEXT]
      join the 23andMe v5 site list (GRCh37 positions and reference bases,
      columns chrom, pos, id, ref) with dbSNP's GRCh38 VCF by rsID. The GRCh38
      position of a site is dbSNP's, not a liftover. A site is dropped, with
      its reason in DROPPED, when it has no rsID (23andMe's own i-numbers),
      dbSNP has no SNV for it, dbSNP places it twice or on another
      chromosome, or neither strand of its GRCh37 base is one of dbSNP's
      alleles. MT sites take dbSNP's rCRS position in both builds (the list
      numbers them on hg19's chrM). setup.sh --consumer-sites runs it.

  export_consumer.py export --gvcf G.g.vcf.gz --sites SITES.tsv.gz --sample ID
                            [--sex male|female] [--par PAR.bed] --out-dir DIR
      genotype every site from the gVCF and write, in DIR:
        ID_23andme_v5.zip           ID_23andme_v5.txt, the 23andMe layout
        ID_ancestrydna_layout.zip   ID_ancestrydna_layout.txt, the same
                                    sites in the AncestryDNA layout
        ID_consumer_export.tsv      what was written (key, value)
      A site inside a reference block (GT 0/0) is homozygous reference; a
      site with no gVCF record, a no-call, a filtered call or one inside a
      deletion is '--' (AncestryDNA: 0 0). Both files are read back with the
      strict reader below before they are zipped. Exit 1 when a file fails
      it, or when the gVCF's reference bases disagree with the site table
      (the gVCF is not on GRCh38, or the table is wrong).

  export_consumer.py check --layout 23andme|ancestrydna FILE (.txt or .zip)
      the strict reader: exit 1 with every problem listed, 0 when the file
      keeps the layout.

Standard library only: it runs in PYTHON_IMAGE and under
tests/test_export_consumer.py.
"""
import argparse
import datetime
import gzip
import io
import os
import re
import sys
import zipfile

DOC_URL = "https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/docs/40-export-consumer.md"

# 23andMe's chromosome names, in the order its files list them.
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
CHROM_ORDER = {c: i for i, c in enumerate(CHROMS)}
# AncestryDNA's codes: 23 X, 24 Y, 25 the X pseudoautosomal regions, 26 MT.
ANCESTRY_CODE = dict({str(i): str(i) for i in range(1, 23)}, X="23", Y="24", MT="26")
# The X pseudoautosomal regions on GRCh37, 1-based inclusive (GRC).
PAR37 = [(60001, 2699520), (154931044, 155260560)]

# dbSNP's GRCh38 RefSeq accessions of the primary chromosomes.
REFSEQ38 = dict({f"NC_0000{i:02d}": str(i) for i in range(1, 23)}, NC_000023="X", NC_000024="Y", NC_012920="MT")

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
BASES = set("ACGT")

SITES_HEADER = ["rsid", "chrom", "pos37", "ref37", "pos38", "ref38", "strand"]

HEADER_23ANDME = "# rsid\tchromosome\tposition\tgenotype"
HEADER_ANCESTRY = "rsid\tchromosome\tposition\tallele1\tallele2"
ANCESTRY_FIRST = "#AncestryDNA raw data download"


def die(msg, code=2):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def open_text(path):
    """A text stream of PATH: '-' is stdin; gzip (and bgzip) is read by its magic bytes."""
    if path == "-":
        return io.TextIOWrapper(sys.stdin.buffer, errors="replace")
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def tsv_chrom(c):
    """chr1 / 1 / chrM / MT -> 23andMe's name, or None for any other contig."""
    c = c[3:] if c.startswith("chr") else c
    c = "MT" if c == "M" else c
    return c if c in CHROM_ORDER else None


def gvcf_contig(c):
    """23andMe's chromosome name -> the GRCh38 contig name of the pipeline."""
    return "chrM" if c == "MT" else f"chr{c}"


# --- build-sites --------------------------------------------------------------

def read_template(path):
    rows = []
    with open_text(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 4:
                die(f"{path} line {n}: want chrom, pos, id, ref; got {line[:80]!r}")
            chrom = tsv_chrom(c[0])
            if chrom is None or not c[1].isdigit():
                die(f"{path} line {n}: unknown chromosome or position: {line[:80]!r}")
            rows.append((c[2], chrom, int(c[1]), c[3].upper()))
    if not rows:
        die(f"{path} holds no site")
    return rows


def build_sites(args):
    template = read_template(args.template)
    want = {}
    for i, (rsid, *_rest) in enumerate(template):
        if rsid.startswith("rs") and rsid[2:].isdigit():
            want.setdefault(rsid, []).append(i)
    found = {}      # rsid -> set of (chrom, pos38, ref38, alts) for SNV records
    other = set()   # rsids dbSNP has only as a non-SNV (indel, MNV)
    with open_text(args.dbsnp) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            if len(c) < 8 or c[2] not in want:
                continue
            chrom = REFSEQ38.get(c[0].split(".")[0])
            if chrom is None:
                continue   # a patch or ALT scaffold: the primary assembly only
            ref, alts = c[3].upper(), [a.upper() for a in c[4].split(",")]
            if ";VC=SNV;" not in f";{c[7].rstrip()};" or len(ref) != 1 or any(len(a) != 1 for a in alts):
                other.add(c[2])
                continue
            found.setdefault(c[2], set()).add((chrom, int(c[1]), ref, tuple(alts)))
    kept, dropped = [], []
    for rsid, chrom, pos37, ref37 in template:
        def drop(why):
            dropped.append((rsid, chrom, pos37, why))
        if rsid not in want:
            drop("no rsID (a 23andMe internal id)")
            continue
        hits = found.get(rsid)
        if not hits:
            drop("not an SNV in dbSNP" if rsid in other else "not in dbSNP (merged or withdrawn)")
            continue
        # A pseudoautosomal site has a record on chrX and one on chrY: the
        # template's chromosome decides.
        same = [h for h in hits if h[0] == chrom]
        if not same:
            drop(f"on chromosome {sorted(hits)[0][0]} in GRCh38")
            continue
        if len({h[1] for h in same}) > 1:
            drop("more than one GRCh38 position in dbSNP")
            continue
        chrom38, pos38, ref38, alts = sorted(same)[0]
        if chrom == "MT":
            # GRCh37's MT and GRCh38's chrM are the same sequence (the rCRS),
            # but the list numbers its MT sites the hg19 chrM way (the Yoruba
            # sequence, one or two bases off): the rCRS position is dbSNP's.
            pos37, ref37 = pos38, ref38
        alleles = {ref38, *alts}
        if ref37 in alleles:
            strand = "+"
        elif COMP.get(ref37) in alleles:
            strand = "-"
        else:
            drop("its GRCh37 base is not a dbSNP allele on either strand")
            continue
        kept.append((rsid, chrom, pos37, ref37, pos38, ref38, strand))
    kept.sort(key=lambda r: (CHROM_ORDER[r[1]], r[2], r[0]))
    tmp = args.out + ".tmp"
    with gzip.open(tmp, "wt") as out:
        out.write(f"## consumer sites: 23andMe v5 chip, GRCh37 positions from the site list, GRCh38 positions from dbSNP\n")
        out.write(f"## source: {args.source or '.'}\n")
        out.write(f"## template_sites: {len(template)}\n")
        out.write(f"## kept: {len(kept)}\n")
        out.write(f"## dropped: {len(dropped)}\n")
        out.write("#" + "\t".join(SITES_HEADER) + "\n")
        for r in kept:
            out.write("\t".join(str(x) for x in r) + "\n")
    with open(args.dropped + ".tmp", "w") as out:
        out.write("rsid\tchrom\tpos37\treason\n")
        for r in dropped:
            out.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp, args.out)
    os.replace(args.dropped + ".tmp", args.dropped)
    reasons = {}
    for r in dropped:
        reasons[r[3]] = reasons.get(r[3], 0) + 1
    print(f"Site list: {len(template)} sites; kept {len(kept)}, dropped {len(dropped)}")
    for why, n in sorted(reasons.items(), key=lambda x: -x[1]):
        print(f"  dropped, {why}: {n}")
    if not kept:
        die("no site kept: is the dbSNP file the GRCh38 VCF?", 1)


# --- export -------------------------------------------------------------------

def read_sites(path):
    """The site table: a list of dicts, plus its '## key: value' lines."""
    meta, sites, cols = {}, [], None
    with open_text(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\r\n")
            if line.startswith("## "):
                k, _, v = line[3:].partition(": ")
                meta[k] = v
                continue
            if line.startswith("#"):
                cols = line[1:].split("\t")
                if cols != SITES_HEADER:
                    die(f"{path}: columns {cols}, want {SITES_HEADER}")
                continue
            if not line:
                continue
            if cols is None:
                die(f"{path}: no '#rsid ...' column line before the sites")
            c = line.split("\t")
            if len(c) != len(SITES_HEADER) or c[1] not in CHROM_ORDER or c[6] not in "+-":
                die(f"{path} line {n}: not a site row: {line[:80]!r}")
            sites.append({"rsid": c[0], "chrom": c[1], "pos37": int(c[2]), "ref37": c[3],
                          "pos38": int(c[4]), "ref38": c[5], "strand": c[6], "gt": None})
    if not sites:
        die(f"{path} holds no site")
    return meta, sites


def read_par(path):
    par = {}
    if not path:
        return par
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) >= 3:
                par.setdefault(c[0], []).append((int(c[1]) + 1, int(c[2])))
    return par


def in_par(par, contig, pos):
    return any(s <= pos <= e for s, e in par.get(contig, []))


class Tally(dict):
    def add(self, k, n=1):
        self[k] = self.get(k, 0) + n


def call_from_record(rec, pos, ref38, tally):
    """The alleles (bases, GRCh38 + strand) at POS from one gVCF record, or a no-call reason."""
    chrom, rpos, ref, alts, flt, fmt, smp = rec
    if fmt is None:
        return None, "no genotype"
    keys = fmt.split(":")
    vals = smp.split(":")
    gt = vals[keys.index("GT")] if "GT" in keys and keys.index("GT") < len(vals) else "."
    idx = re.split(r"[/|]", gt)
    if any(i in (".", "") for i in idx):
        return None, "no-call in the gVCF"
    alleles = [ref] + alts
    try:
        chosen = [alleles[int(i)] for i in idx]
    except (ValueError, IndexError):
        return None, "unreadable genotype"
    if rpos != pos:
        # The site lies inside a record that starts before it: a reference
        # block (each base is the reference) or a deletion.
        if all(i == "0" for i in idx) and alts and all(a in SYMBOLIC for a in alts):
            return [ref38] * len(idx), None
        return None, "inside a deletion or a longer allele"
    if ref[0] != ref38:
        tally.add("reference_mismatch")
        return None, "reference base differs from the site table"
    tally.add("reference_checked")
    if flt not in ("PASS", ".") and any(i != "0" for i in idx):
        return None, "filtered call"
    bases = []
    for a in chosen:
        b = a[0] if a and a[0] in BASES else None
        if b is None:
            return None, "symbolic allele"
        bases.append(b)
    return bases, None


SYMBOLIC = ("<*>", "<NON_REF>", ".")


def is_block(rec):
    return all(a in SYMBOLIC for a in rec[3])


def parse_record(line):
    """(record, end) of one gVCF line; record = (chrom, pos, ref, alts, filter, format, sample)."""
    c = line.rstrip("\n").split("\t")
    if len(c) < 8:
        return None
    end = int(c[1]) + len(c[3]) - 1
    m = re.search(r"(?:^|;)END=(\d+)", c[7])
    if m:
        end = max(end, int(m.group(1)))
    fmt, smp = (c[8], c[9]) if len(c) > 9 else (None, None)
    return (c[0], int(c[1]), c[3].upper(), [a.upper() for a in c[4].split(",")], c[6], fmt, smp), end


def nonref_gt(rec):
    if rec[5] is None:
        return False
    keys, vals = rec[5].split(":"), rec[6].split(":")
    if "GT" not in keys or keys.index("GT") >= len(vals):
        return False
    return any(i not in ("0", ".", "") for i in re.split(r"[/|]", vals[keys.index("GT")]))


def genotype_all(gvcf, sites, tally):
    """Walk the gVCF once and give each site its alleles (or a no-call reason).

    A site's call comes from the records that start last at or before it: a
    record starting at the site (a variant before a block), else a reference
    block that covers it. A site inside the span of a variant whose genotype
    carries a longer allele (a deletion, an MNV) is '--' whatever blocks
    cover it: DeepVariant starts its next block one base after the variant,
    inside the deleted bases. Only lines near a site are parsed.
    """
    by_contig = {}
    for s in sites:
        by_contig.setdefault(gvcf_contig(s["chrom"]), []).append(s)
    for lst in by_contig.values():
        lst.sort(key=lambda s: s["pos38"])
    done_contigs = set()
    cur, lst, i = None, [], 0
    pending = []   # the lines of the records that start at the last position seen
    spans = []     # (first, last) base inside a non-reference longer allele

    def resolve(upto):
        """Give the sites before UPTO their call from the pending records."""
        nonlocal i
        recs = [r for r in (parse_record(l) for l in pending) if r]
        start = recs[0][0][1] if recs else None
        while i < len(lst) and lst[i]["pos38"] < upto:
            s, p = lst[i], lst[i]["pos38"]
            i += 1
            spans[:] = [sp for sp in spans if sp[1] >= p]
            hit = [r for r, end in recs if start <= p <= end]
            if any(a <= p <= b for a, b in spans):
                s["why"] = "inside a deletion or a longer allele"
                continue
            if not hit:
                s["why"] = "no gVCF record"
                continue
            # At the record's start a variant wins over a block; past it, a block.
            hit.sort(key=lambda r: is_block(r) == (p == start))
            bases, why = call_from_record(hit[0], p, s["ref38"], tally)
            if bases is None:
                s["why"] = why
            else:
                s["gt"] = bases
                tally.add("from_reference_blocks" if is_block(hit[0]) or hit[0][1] != p else "from_variant_records")

    last_pos = 0
    with open_text(gvcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t1 = line.find("\t")
            t2 = line.find("\t", t1 + 1)
            contig, pos = line[:t1], int(line[t1 + 1:t2])
            if contig != cur:
                if cur is not None:
                    resolve(float("inf"))
                    done_contigs.add(cur)
                if contig in done_contigs:
                    die(f"{gvcf} is not sorted: {contig} appears twice", 1)
                cur, lst, i, pending, spans, last_pos = contig, by_contig.get(contig, []), 0, [], [], 0
            if not lst or i >= len(lst):
                continue
            if pos < last_pos:
                die(f"{gvcf} is not sorted: {contig}:{pos} after {last_pos}", 1)
            if pos != last_pos:
                # Records start in order, so a site before this one gets its
                # call from the pending records; only then are they parsed.
                if lst[i]["pos38"] < pos:
                    resolve(pos)
                pending = [line]
                last_pos = pos
            else:
                pending.append(line)
            t3 = line.find("\t", t2 + 1)
            t4 = line.find("\t", t3 + 1)
            if t4 - t3 > 2:
                r = parse_record(line)
                if r and not is_block(r[0]) and nonref_gt(r[0]):
                    spans.append((pos + 1, r[1]))
        if cur is not None:
            resolve(float("inf"))
    for s in sites:
        if s["gt"] is None and "why" not in s:
            s["why"] = "no gVCF record"


def finalise(sites, sex, par38, tally):
    """Ploidy rules of the chip files, and the GRCh37 strand: each site's letters."""
    for s in sites:
        bases = s["gt"]
        contig = gvcf_contig(s["chrom"])
        if bases is not None and s["chrom"] == "Y" and sex == "female":
            bases, s["why"] = None, "chrY site of a female sample"
        if bases is not None:
            haploid = s["chrom"] == "MT" or (sex == "male" and s["chrom"] in ("X", "Y")
                                             and not in_par(par38, contig, s["pos38"]))
            if haploid and len(bases) == 2:
                if bases[0] != bases[1]:
                    bases, s["why"] = None, "two alleles on a haploid site (chrM heteroplasmy, or male chrX/chrY)"
                else:
                    bases = bases[:1]
        if bases is not None and s["strand"] == "-":
            bases = [COMP[b] for b in bases]
        s["letters"] = "".join(sorted(bases)) if bases is not None else None
        if s["letters"] is None:
            tally.add(f"no_call: {s['why']}")


def write_23andme(path, sites, meta, today):
    with open(path, "w", newline="\n") as f:
        f.write("# This data file was generated by Personal Genome Pipeline (step 40) in the layout of a 23andMe raw data file.\n")
        f.write("# It was not made by 23andMe: the genotypes come from whole-genome sequencing calls, not from a chip.\n")
        f.write(f"# Written {today}. Sites: the 23andMe v5 chip list, {len(sites)} of {meta.get('template_sites', '?')} sites.\n")
        f.write("# We are using reference human assembly build 37 (also known as GRCh37); genotypes are on its + strand.\n")
        f.write("# A genotype of -- means the calls hold no genotype for that site.\n")
        f.write(f"# More: {DOC_URL}\n")
        f.write("#\n")
        f.write(HEADER_23ANDME + "\n")
        for s in sorted(sites, key=lambda s: (CHROM_ORDER[s["chrom"]], s["pos37"], s["rsid"])):
            f.write(f"{s['rsid']}\t{s['chrom']}\t{s['pos37']}\t{s['letters'] or '--'}\n")


def ancestry_code(s):
    if s["chrom"] == "X" and any(a <= s["pos37"] <= b for a, b in PAR37):
        return "25"
    return ANCESTRY_CODE[s["chrom"]]


def write_ancestry(path, sites, meta, today):
    rows = sorted(((int(ancestry_code(s)), s["pos37"], s["rsid"], s) for s in sites), key=lambda r: r[:3])
    with open(path, "w", newline="\n") as f:
        f.write(ANCESTRY_FIRST + "\n")
        f.write("#This file was generated by Personal Genome Pipeline (step 40) in the layout of an AncestryDNA raw data file.\n")
        f.write("#It was not made by AncestryDNA: the genotypes come from whole-genome sequencing calls, not from a chip,\n")
        f.write(f"#at the sites of the 23andMe v5 chip list ({len(sites)} of {meta.get('template_sites', '?')} sites). Written {today}.\n")
        f.write("#Positions are on human assembly build 37 (GRCh37), alleles on its + strand. 0 0 means no genotype.\n")
        f.write("#Chromosome codes: 1-22, 23 = X, 24 = Y, 25 = X pseudoautosomal, 26 = MT.\n")
        f.write(f"#More: {DOC_URL}\n")
        f.write(HEADER_ANCESTRY + "\n")
        for code, pos, rsid, s in rows:
            g = s["letters"]
            a, b = ("0", "0") if not g else (g[0], g[-1])
            f.write(f"{rsid}\t{code}\t{pos}\t{a}\t{b}\n")


def zip_one(txt, zpath):
    tmp = zpath + ".tmp"
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.write(txt, arcname=os.path.basename(txt))
    os.replace(tmp, zpath)
    os.remove(txt)


def export(args):
    meta, sites = read_sites(args.sites)
    par38 = read_par(args.par)
    tally = Tally()
    genotype_all(args.gvcf, sites, tally)
    finalise(sites, args.sex or "", par38, tally)
    checked, bad = tally.get("reference_checked", 0), tally.get("reference_mismatch", 0)
    if bad > max(10, (checked + bad) // 100):
        die(f"{bad} of {checked + bad} gVCF records at a site have another reference base than the site table: "
            "the gVCF is not on GRCh38, or the site table is wrong. Nothing was written.", 1)
    os.makedirs(args.out_dir, exist_ok=True)
    today = datetime.date.today().isoformat()
    base = os.path.join(args.out_dir, args.sample)
    t23, tanc = f"{base}_23andme_v5.txt", f"{base}_ancestrydna_layout.txt"
    write_23andme(t23, sites, meta, today)
    write_ancestry(tanc, sites, meta, today)
    problems = check_23andme(t23) + check_ancestry(tanc)
    if problems:
        for p in problems[:20]:
            print(f"  {p}", file=sys.stderr)
        die(f"the written files fail the strict reader ({len(problems)} problems); they are left unzipped in {args.out_dir}", 1)
    zip_one(t23, f"{base}_23andme_v5.zip")
    zip_one(tanc, f"{base}_ancestrydna_layout.zip")
    called = sum(1 for s in sites if s["letters"])
    rows = [("sites", len(sites)), ("called", called), ("no_call", len(sites) - called),
            ("from_reference_blocks", tally.get("from_reference_blocks", 0)),
            ("from_variant_records", tally.get("from_variant_records", 0)),
            ("reference_checked", checked), ("reference_mismatch", bad),
            ("template_sites", meta.get("template_sites", ".")), ("table_dropped", meta.get("dropped", ".")),
            ("sex", args.sex or "not given"), ("file_23andme", os.path.basename(base) + "_23andme_v5.zip"),
            ("file_ancestrydna", os.path.basename(base) + "_ancestrydna_layout.zip")]
    rows += sorted((k, v) for k, v in tally.items() if k.startswith("no_call: "))
    tmp = f"{base}_consumer_export.tsv.tmp"
    with open(tmp, "w") as f:
        f.write("key\tvalue\n")
        for k, v in rows:
            f.write(f"{k}\t{v}\n")
    os.replace(tmp, f"{base}_consumer_export.tsv")
    pct = 100.0 * called / len(sites)
    print(f"Sites: {len(sites)} (of {meta.get('template_sites', '?')} on the chip list; "
          f"{meta.get('dropped', '?')} dropped when the table was built)")
    print(f"Called: {called} ({pct:.1f}%), {tally.get('from_reference_blocks', 0)} of them from reference blocks; "
          f"no-call (--): {len(sites) - called}")
    for k, v in sorted(tally.items()):
        if k.startswith("no_call: "):
            print(f"  {k[9:]}: {v}")
    print(f"Reference bases checked at {checked} variant records, {bad} differ from the table")


# --- check (the strict reader) ------------------------------------------------

def read_layout_text(path):
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            names = z.namelist()
            if len(names) != 1 or not names[0].endswith(".txt"):
                return None, [f"{path}: want one .txt member, found {names}"]
            data = z.read(names[0])
    else:
        with open(path, "rb") as f:
            data = f.read()
    try:
        return data.decode("ascii"), []
    except UnicodeDecodeError:
        return None, [f"{path}: not plain ASCII"]


def _check_rows(path, text, comment_rule, header, ncols, row_rule):
    problems = []
    if "\r" in text:
        problems.append(f"{path}: carriage returns (CRLF line ends)")
    if not text.endswith("\n"):
        problems.append(f"{path}: the last line has no newline (a cut file)")
    lines = text.split("\n")[:-1] if text.endswith("\n") else text.split("\n")
    i = 0
    while i < len(lines) and comment_rule(lines[i]):
        i += 1
    if i >= len(lines) or lines[i] != header:
        problems.append(f"{path}: the column header {header!r} is missing after the comment lines")
        return problems, []
    rows = []
    for n, line in enumerate(lines[i + 1:], i + 2):
        c = line.split("\t")
        if len(c) != ncols:
            problems.append(f"{path} line {n}: {len(c)} columns, want {ncols}: {line[:60]!r}")
            continue
        err = row_rule(c)
        if err:
            problems.append(f"{path} line {n}: {err}: {line[:60]!r}")
        else:
            rows.append(c)
    if not rows:
        problems.append(f"{path}: no data row")
    return problems, rows


RSID = re.compile(r"^(rs|i)[1-9][0-9]*$")
POS = re.compile(r"^[1-9][0-9]*$")


def check_23andme(path, rows_out=None):
    text, problems = read_layout_text(path)
    if text is None:
        return problems

    def row(c):
        if not RSID.match(c[0]):
            return "rsid is not rs<number> or i<number>"
        if c[1] not in CHROM_ORDER:
            return "chromosome is not 1-22, X, Y or MT"
        if not POS.match(c[2]):
            return "position is not a positive whole number"
        g = c[3]
        if g == "--":
            return None
        if not (1 <= len(g) <= 2) or any(b not in BASES for b in g):
            return "genotype is not one or two of A, C, G, T, or --"
        if len(g) == 1 and c[1] not in ("X", "Y", "MT"):
            return "one letter on an autosome"
        if len(g) == 2 and c[1] == "MT":
            return "two letters on MT"
        return None
    problems, rows = _check_rows(path, text, lambda l: l.startswith("#") and l != HEADER_23ANDME,
                                 HEADER_23ANDME, 4, row)
    if not text.startswith("# This data file"):
        problems.append(f"{path}: the first line does not start with '# This data file'")
    problems += _order_and_ids(path, [(CHROM_ORDER[c[1]], int(c[2]), c[0]) for c in rows])
    if rows_out is not None:
        rows_out.extend(rows)
    return problems


def check_ancestry(path, rows_out=None):
    text, problems = read_layout_text(path)
    if text is None:
        return problems

    def row(c):
        if not RSID.match(c[0]):
            return "rsid is not rs<number> or i<number>"
        if c[1] not in {str(i) for i in range(1, 27)}:
            return "chromosome is not a code from 1 to 26"
        if not POS.match(c[2]):
            return "position is not a positive whole number"
        if (c[3], c[4]) == ("0", "0"):
            return None
        if c[3] not in BASES or c[4] not in BASES:
            return "an allele is not A, C, G or T (and not a 0 0 no-call)"
        return None
    problems, rows = _check_rows(path, text, lambda l: l.startswith("#"), HEADER_ANCESTRY, 5, row)
    if not text.startswith(ANCESTRY_FIRST + "\n"):
        problems.append(f"{path}: the first line is not {ANCESTRY_FIRST!r}")
    problems += _order_and_ids(path, [(int(c[1]), int(c[2]), c[0]) for c in rows])
    if rows_out is not None:
        rows_out.extend(rows)
    return problems


def _order_and_ids(path, keys):
    problems, seen = [], set()
    for a, b in zip(keys, keys[1:]):
        if b[:2] < a[:2]:
            problems.append(f"{path}: rows out of order at {b[2]} (after {a[2]})")
            break
    for k in keys:
        if k[2] in seen:
            problems.append(f"{path}: {k[2]} appears twice")
            break
        seen.add(k[2])
    return problems


def check(args):
    problems = check_23andme(args.file) if args.layout == "23andme" else check_ancestry(args.file)
    for p in problems:
        print(p)
    if problems:
        sys.exit(1)
    print(f"{args.file}: keeps the {args.layout} layout")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build-sites", help="join the 23andMe v5 site list with dbSNP's GRCh38 VCF")
    b.add_argument("--template", required=True)
    b.add_argument("--dbsnp", required=True, help="dbSNP's GRCh38 VCF (plain or gzip), or - for stdin")
    b.add_argument("--out", required=True)
    b.add_argument("--dropped", required=True)
    b.add_argument("--source", default="")
    e = sub.add_parser("export", help="write the 23andMe and AncestryDNA layouts from a gVCF")
    e.add_argument("--gvcf", required=True)
    e.add_argument("--sites", required=True)
    e.add_argument("--sample", required=True)
    e.add_argument("--sex", choices=["male", "female", ""], default="")
    e.add_argument("--par", help="the GRCh38 pseudoautosomal regions, BED (assets/par_grch38.bed)")
    e.add_argument("--out-dir", required=True)
    c = sub.add_parser("check", help="the strict reader of the two layouts")
    c.add_argument("--layout", choices=["23andme", "ancestrydna"], required=True)
    c.add_argument("file")
    args = ap.parse_args(argv)
    {"build-sites": build_sites, "export": export, "check": check}[args.cmd](args)


if __name__ == "__main__":
    main()
