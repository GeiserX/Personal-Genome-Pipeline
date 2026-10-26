#!/usr/bin/env python3
"""genome_lookup.py: answer "what is my genotype at rsX, at chr:pos, or in gene Y".

scripts/pgp-lookup.sh runs it in two passes (docs/lookup.md):

  genome_lookup.py resolve --query Q [--sites SITES.tsv.gz] [--table TSV ...]
                           [--genes GENES.gtf] [--annotated VCF]
      turn the query into the region(s) to read, one per line:
        kind  chrom  start  end  ref  source
      kind is 'site' (rsID or chr:pos; start = end) or 'gene' (GENCODE's span).
      ref is the GRCh38 base when a table gives it, else '.'.
      An rsID is looked up in the consumer site table (setup.sh
      --consumer-sites), then in the given tables (assets/traits.tsv,
      assets/bloodtype_markers.tsv), then in the ID column and VEP's
      Existing_variation of the annotated VCF (the sample's own variants
      only). A gene is looked up in GENCODE's gene lines. Exit 3 with the
      reason on stderr when nothing resolves.

  genome_lookup.py report --query Q --kind site|gene --chrom C --start N --end N
                          [--ref BASE] [--source TEXT] --gvcf G.vcf[.gz]
                          [--annotated A.vcf[.gz]] [--summary summary.json]
      print the answer. A site: the genotype from the gVCF (a position inside
      a reference block is ref/ref, a position no record covers is a no
      call), its depth and quality, and from the annotated VCF the gene,
      consequence, ClinVar significance and gnomAD frequency of the variant
      called there. A gene: every variant called in it, with its consequence,
      ClinVar and gnomAD, and how many of the gene's bases have a call. Both
      end with what the pipeline's reports (summary.json) say about the gene
      or rsID. The inputs may be whole files or the slices pgp-lookup.sh cuts
      with bcftools; records outside the region are skipped either way.

Standard library only: it runs in PYTHON_IMAGE and under
tests/test_genome_lookup.py.
"""
import argparse
import gzip
import json
import re
import sys

# CSQ fields worth printing, in order, with the label they get. Only the ones
# the annotated VCF's CSQ header defines and that hold a value are printed.
CSQ_SHOW = [
    ("Consequence", "consequence"),
    ("IMPACT", "impact"),
    ("HGVSc", "HGVSc"),
    ("HGVSp", "HGVSp"),
    ("Existing_variation", "known as"),
    ("ClinVar_CLNSIG", "ClinVar (pipeline's ClinVar file)"),
    ("ClinVar_CLNREVSTAT", "ClinVar review status"),
    ("ClinVar_CLNDN", "ClinVar condition"),
    ("CLIN_SIG", "clinical significance (VEP cache)"),
    ("gnomADg_AF", "gnomAD genomes AF"),
    ("gnomADe_AF", "gnomAD exomes AF"),
    ("MAX_AF", "highest population AF"),
]
# Keys of a summary.json entry that name a gene.
GENE_KEYS = ("gene", "genes", "Gene", "SYMBOL", "symbol")
SYMBOLIC = {"<*>", "<NON_REF>", "*"}
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def open_text(path):
    """A text stream of PATH: '-' is stdin; gzip (and bgzip) is read by its magic bytes."""
    if path == "-":
        return sys.stdin
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def norm_chrom(c):
    """chr7 / 7 -> 7; chrM / M / MT -> MT: contig names compared across files."""
    c = c[3:] if c.lower().startswith("chr") else c
    return "MT" if c.upper() in ("M", "MT") else c.upper()


def contig_name(c):
    """The pipeline's GRCh38 contig name of a normalised chromosome."""
    return "chrM" if c == "MT" else f"chr{c}"


# --- the query ----------------------------------------------------------------

def parse_query(q):
    """('rsid', 'rs123') | ('site', ('7', 117559590)) | ('gene', 'CFTR'), or raise ValueError."""
    q = q.strip()
    if re.fullmatch(r"(?i)rs\d+", q):
        return "rsid", q.lower()
    m = re.fullmatch(r"(?i)(?:chr)?([0-9]{1,2}|X|Y|M|MT)[:_-]([0-9][0-9,]*)", q)
    if m:
        pos = int(m.group(2).replace(",", ""))
        if pos < 1:
            raise ValueError(f"position must be 1 or more: {q}")
        return "site", (norm_chrom(m.group(1)), pos)
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", q):
        return "gene", q
    raise ValueError(f"not an rsID (rs123), a position (chr7:117559590) or a gene symbol (CFTR): {q!r}")


# --- resolve ------------------------------------------------------------------

def sites_rsid(path, rsid):
    """(chrom, pos38, ref38) rows of the consumer site table for RSID.
    The table (bin/export_consumer.py) opens with "## key: value" lines and a
    "#rsid ..." column line."""
    out = []
    with open_text(path) as f:
        head = None
        for line in f:
            if line.startswith("##") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if head is None:
                head = {k.lstrip("#"): i for i, k in enumerate(c)}
                if not {"rsid", "chrom", "pos38", "ref38"} <= set(head):
                    raise ValueError(f"{path}: not a consumer site table (columns {', '.join(c)})")
                continue
            if c[head["rsid"]].lower() == rsid:
                out.append((norm_chrom(c[head["chrom"]]), int(c[head["pos38"]]), c[head["ref38"]]))
    return out


def table_rsid(path, rsid):
    """(chrom, pos, ref) rows of a pipeline table (traits.tsv: chrom pos rsid ref;
    bloodtype_markers.tsv: a BED, chrom start end ... ref ... rsid) for RSID."""
    out = []
    with open_text(path) as f:
        head = None
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if head is None:
                head = {k: i for i, k in enumerate(c)}
                continue
            row = dict(zip(head, c))
            if (row.get("rsid") or "").lower() != rsid:
                continue
            if "pos" in row:
                pos = int(row["pos"])
            elif "end" in row:
                pos = int(row["end"])  # BED: end is the 1-based position of the last base
            else:
                continue
            ref = row.get("ref") or "."
            out.append((norm_chrom(row["chrom"]), pos, ref[:1] if ref != "." else "."))
    return out


def csq_format(header_lines):
    """The CSQ field names from the ##INFO=<ID=CSQ,...Format: a|b|c"> header line, or []."""
    for h in header_lines:
        if h.startswith("##INFO=<ID=CSQ,"):
            m = re.search(r'Format: ([^"]+)"', h)
            if m:
                return m.group(1).strip().split("|")
    return []


class Vcf:
    """A VCF read once: its header, CSQ field names, first sample, records.

    With chrom, start and end, only the records that overlap that region
    (reference blocks by their END) are kept, so a whole unindexed genome
    streams through in bounded memory; the reports filter on the same rule.
    """

    def __init__(self, path, chrom=None, start=None, end=None):
        self.path = path
        self.header = []
        self.samples = []
        self.records = []
        want = norm_chrom(chrom) if chrom is not None else None
        with open_text(path) as f:
            for line in f:
                if line.startswith("##"):
                    self.header.append(line.rstrip("\n"))
                elif line.startswith("#"):
                    self.samples = line.rstrip("\n").split("\t")[9:]
                else:
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 8:
                        continue
                    if want is not None:
                        if norm_chrom(c[0]) != want or not c[1].isdigit() or int(c[1]) > end:
                            continue
                        r = Record(c)
                        if r.end < start:
                            continue
                        self.records.append(r)
                    else:
                        self.records.append(Record(c))
        self.csq = csq_format(self.header)


class Record:
    def __init__(self, c):
        self.chrom = norm_chrom(c[0])
        self.pos = int(c[1])
        self.id = c[2]
        self.ref = c[3]
        self.alts = [] if c[4] == "." else c[4].split(",")
        self.qual = c[5]
        self.filter = c[6]
        self.info = {}
        for kv in c[7].split(";"):
            k, _, v = kv.partition("=")
            self.info[k] = v if _ else True
        self.fmt = {}
        if len(c) >= 10:
            self.fmt = dict(zip(c[8].split(":"), c[9].split(":")))
        end = self.info.get("END")
        self.end = int(end) if isinstance(end, str) and end.isdigit() else self.pos + max(len(self.ref), 1) - 1

    def real_alts(self):
        return [a for a in self.alts if a not in SYMBOLIC]

    def is_block(self):
        return not self.real_alts()

    def gt(self):
        """The allele indexes of GT (None for '.'), and whether it is phased."""
        g = self.fmt.get("GT", "")
        if not g:
            return [], False
        return [None if a == "." else int(a) for a in re.split(r"[/|]", g)], "|" in g

    def bases(self):
        alleles = [self.ref] + self.alts
        idx, phased = self.gt()
        return [None if i is None or i >= len(alleles) else alleles[i] for i in idx], phased

    def has_alt_call(self):
        idx, _ = self.gt()
        return any(i is not None and 0 < i <= len(self.alts) and self.alts[i - 1] not in SYMBOLIC for i in idx)

    def csq_entries(self, fields):
        raw = self.info.get("CSQ")
        if not fields or not isinstance(raw, str):
            return []
        return [dict(zip(fields, e.split("|"))) for e in raw.split(",")]

    def ids(self, fields):
        """Every rsID of the record: its ID column and VEP's Existing_variation."""
        out = {i for i in self.id.split(";") if i != "."}
        for e in self.csq_entries(fields):
            out.update(x for x in (e.get("Existing_variation") or "").split("&") if x)
        return out


def annotated_rsid(path, rsid):
    """(chrom, pos, ref) of the annotated VCF records whose ID or Existing_variation is RSID.
    Reads the whole file, so the cheap substring test comes first."""
    out = []
    fields = []
    pat = re.compile(r"(?<![0-9A-Za-z])" + re.escape(rsid) + r"(?![0-9])", re.I)
    with open_text(path) as f:
        head = []
        for line in f:
            if line.startswith("#"):
                head.append(line)
                if line.startswith("#CHROM"):
                    fields = csq_format(head)
                continue
            if not pat.search(line):
                continue
            r = Record(line.rstrip("\n").split("\t"))
            if rsid in {i.lower() for i in r.ids(fields)}:
                out.append((r.chrom, r.pos, r.ref[:1]))
    return out


def gene_spans(path, gene):
    """(chrom, start, end, source) of GENCODE's gene lines named GENE: the exact
    name first, else a case-insensitive match. More than one for a gene in both
    X and Y pseudoautosomal regions."""
    m = re.search(r"gencode\.v(\d+)", path)
    source = f"GENCODE {m.group(1)}" if m else "GENCODE"
    exact, loose = [], []
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            n = re.search(r'gene_name "([^"]+)"', c[8])
            if not n:
                continue
            span = (norm_chrom(c[0]), int(c[3]), int(c[4]), f"{source} ({n.group(1)})")
            if n.group(1) == gene:
                exact.append(span)
            elif n.group(1).upper() == gene.upper():
                loose.append(span)
    return exact or loose


def resolve(a):
    kind, val = parse_query(a.query)
    rows = []
    if kind == "site":
        rows.append(("site", val[0], val[1], val[1], ".", "position as given"))
    elif kind == "rsid":
        seen = set()
        hits = []
        if a.sites:
            hits += [(h, "the consumer site table (dbSNP positions)") for h in sites_rsid(a.sites, val)]
        for t in a.table or []:
            hits += [(h, t.rsplit("/", 1)[-1]) for h in table_rsid(t, val)]
        if not hits and a.annotated:
            hits += [(h, "the annotated VCF (a variant called in this sample)") for h in annotated_rsid(a.annotated, val)]
        for (c, p, ref), src in hits:
            if (c, p) not in seen:
                seen.add((c, p))
                rows.append(("site", c, p, p, ref or ".", src))
        if not rows:
            where = ["the consumer site table" if a.sites else "no consumer site table (setup.sh --consumer-sites builds one)"]
            where += ["the pipeline's tables"] if a.table else []
            where += ["the annotated VCF" if a.annotated else "no annotated VCF (step 13)"]
            print(f"ERROR: {val} was not found in {', '.join(where)}. An rsID the sample does not carry is only known "
                  "through the site table; look it up on dbSNP (https://www.ncbi.nlm.nih.gov/snp/) and ask for its "
                  "GRCh38 position instead, as chr:pos.", file=sys.stderr)
            return 3
    else:
        if not a.genes:
            print("ERROR: no GENCODE gene coordinates (reference/gencode.v*.basic.genes.gtf): run scripts/setup.sh, "
                  "which installs them, or ask for a position (chr:pos).", file=sys.stderr)
            return 3
        for c, s, e, src in gene_spans(a.genes, val):
            rows.append(("gene", c, s, e, ".", src))
        if not rows:
            print(f"ERROR: {val} is not a gene name in {a.genes}. Use the HGNC symbol (for example CFTR, not "
                  "its alias), an rsID or a position (chr:pos).", file=sys.stderr)
            return 3
    for r in rows:
        print("\t".join([r[0], contig_name(r[1]), str(r[2]), str(r[3]), r[4], r[5]]))
    return 0


# --- report -------------------------------------------------------------------

def pick_entry(entries, gene=None):
    """The CSQ entry to print: of GENE when given; then the most severe IMPACT,
    and among those VEP's PICK, MANE Select, the canonical transcript, else
    the first."""
    if gene:
        mine = [e for e in entries if (e.get("SYMBOL") or "").upper() == gene.upper()]
        entries = mine or entries
    if not entries:
        return {}
    best = min(IMPACT_RANK.get(e.get("IMPACT"), len(IMPACT_RANK)) for e in entries)
    entries = [e for e in entries if IMPACT_RANK.get(e.get("IMPACT"), len(IMPACT_RANK)) == best]
    for test in (lambda e: e.get("PICK") == "1", lambda e: bool(e.get("MANE_SELECT")),
                 lambda e: e.get("CANONICAL") == "YES"):
        for e in entries:
            if test(e):
                return e
    return entries[0]


def zygosity(bases):
    called = [b for b in bases if b is not None]
    if not bases or not called:
        return "no call"
    if len(called) < len(bases):
        return "half call"
    if len(bases) == 1:
        return "hemizygous"
    return "homozygous" if len(set(called)) == 1 else "heterozygous"


def fmt_gt(bases, phased):
    return ("|" if phased else "/").join("." if b is None else b for b in bases) or "."


def depth_text(r):
    parts = []
    for k, label in (("DP", "depth"), ("MIN_DP", "lowest depth in the block"), ("GQ", "genotype quality")):
        v = r.fmt.get(k)
        if v not in (None, "", "."):
            parts.append(f"{label} {v}")
    if r.qual not in (".", "") and not r.is_block():
        parts.append(f"QUAL {r.qual}")
    return ", ".join(parts) or "no depth or quality recorded"


def filtered(r):
    return r.filter not in ("PASS", ".", "")


def site_genotype(gvcf, chrom, pos, ref):
    """(status, lines) for one position of the gVCF."""
    cover = [r for r in gvcf.records if r.chrom == chrom and r.pos <= pos <= r.end]
    if not cover:
        return "no call", ["genotype: no call (no gVCF record covers this position: no reads there, "
                           "or the caller made no call)"]
    at = [r for r in cover if r.pos == pos and not r.is_block()]
    r = at[0] if at else ([x for x in cover if x.is_block()] or cover)[0]
    bases, phased = r.bases()
    lines = []
    if r.is_block():
        idx, _ = r.gt()
        if not idx or all(i is None for i in idx):
            status = "no call"
            lines.append(f"genotype: no call (reference block {contig_name(r.chrom)}:{r.pos}-{r.end} has no genotype)")
        else:
            base = ref if ref not in (None, "", ".") else (r.ref if r.pos == pos else None)
            gt = "/".join(base if i == 0 else "." for i in idx) if base else "/".join("ref" if i == 0 else "." for i in idx)
            status = "reference" if all(i == 0 for i in idx) else "no call"
            what = "homozygous reference" if status == "reference" else zygosity([0 if i == 0 else None for i in idx])
            lines.append(f"genotype: {gt} ({what}, inside the reference block {contig_name(r.chrom)}:{r.pos}-{r.end})")
    elif r.pos != pos:
        status = "inside a variant"
        lines.append(f"genotype: this position lies inside the variant called at {contig_name(r.chrom)}:{r.pos} "
                     f"{r.ref}>{','.join(r.real_alts())}, genotype {fmt_gt(bases, phased)} ({zygosity(bases)})")
    else:
        z = zygosity(bases)
        status = "no call" if z == "no call" else ("reference" if all(b == r.ref for b in bases if b) else "variant")
        lines.append(f"genotype: {fmt_gt(bases, phased)} ({z}; REF {r.ref}, ALT {','.join(r.real_alts()) or '.'})")
    lines.append(f"evidence: {depth_text(r)}")
    if filtered(r):
        lines.append(f"filter: {r.filter} (the caller did not pass this call; read it as uncertain)")
        status = "filtered"
    return status, lines


def annotation_lines(r, fields, gene=None, indent="  "):
    e = pick_entry(r.csq_entries(fields), gene)
    out = []
    sym = e.get("SYMBOL")
    if sym:
        out.append(f"{indent}gene: {sym}")
    for k, label in CSQ_SHOW:
        v = e.get(k) or ""
        if k.startswith("ClinVar_") or k == "CLIN_SIG":
            # VEP writes a comma inside a value as '&'; ClinVar writes spaces as '_'.
            v = re.sub(r"\s*,\s*", ", ", v.replace("&", ",").replace("_", " ")).strip()
        else:
            v = v.replace("&", ", ")
        if v:
            out.append(f"{indent}{label}: {v}")
    for k in ("CLNSIG", "CLNREVSTAT", "CLNDN"):
        v = r.info.get(k)
        if isinstance(v, str) and v and not e.get(f"ClinVar_{k}"):
            out.append(f"{indent}ClinVar {k}: {v.replace('_', ' ')}")
    if fields and not e:
        out.append(f"{indent}(no VEP annotation on this record)")
    return out


def symbols(r, fields):
    return {e.get("SYMBOL") for e in r.csq_entries(fields) if e.get("SYMBOL")}


def walk_mentions(obj, want_genes, want_rsid, path=()):
    """Every dict in OBJ that names one of WANT_GENES (by a gene key) or WANT_RSID (by rsid)."""
    if isinstance(obj, dict):
        hit = False
        for k in GENE_KEYS:
            v = obj.get(k)
            names = set()
            if isinstance(v, str):
                names = {x for x in re.split(r"[,;&/ ]+", v) if x}
            elif isinstance(v, list):
                names = {x for x in v if isinstance(x, str)}
            elif isinstance(v, dict):
                names = set(v)
            if {n.upper() for n in names} & want_genes:
                hit = True
        if want_rsid and str(obj.get("rsid", "")).lower() == want_rsid:
            hit = True
        if hit:
            yield obj
            return
        for v in obj.values():
            yield from walk_mentions(v, want_genes, want_rsid)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk_mentions(v, want_genes, want_rsid)


def flat(d, limit=300):
    parts = []
    for k, v in d.items():
        if isinstance(v, (str, int, float, bool)) and v not in ("", None):
            parts.append(f"{k}={v}")
        elif isinstance(v, list) and v and all(isinstance(x, (str, int, float)) for x in v):
            parts.append(f"{k}={','.join(map(str, v))}")
        elif isinstance(v, dict) and v and all(isinstance(x, (str, int, float)) for x in v.values()):
            parts.append(f"{k}={'; '.join(f'{a}: {b}' for a, b in v.items())}")
    s = ", ".join(parts)
    return s if len(s) <= limit else s[:limit - 3] + "..."


def summary_lines(path, genes, rsid):
    if not path:
        return ["  no summary.json for this sample: run scripts/generate-report.sh <sample> (on Nextflow output, "
                "GENOME_DIR=<outdir>) to add what the reports say"]
    with open(path) as f:
        summary = json.load(f)
    want = {g.upper() for g in genes if g}
    out = []
    for key, sec in (summary.get("sections") or {}).items():
        if sec.get("state") not in ("ok", "stale"):
            continue
        hits = list(walk_mentions(sec.get("data") or {}, want, rsid))
        if not hits:
            continue
        stale = " (stale: from an earlier run)" if sec.get("state") == "stale" else ""
        out.append(f"  {sec.get('title') or key} (step {sec.get('step')}){stale}:")
        for h in hits[:20]:
            out.append(f"    - {flat(h)}")
        if len(hits) > 20:
            out.append(f"    - ... and {len(hits) - 20} more in {sec.get('source') or 'the section'}")
    return out or ["  none of the reports mentions " + (", ".join(sorted(want | ({rsid} if rsid else set()))) or "it")]


def report_site(a, gvcf, ann):
    chrom, pos = norm_chrom(a.chrom), a.start
    print(f"== {a.query}: {contig_name(chrom)}:{pos} (GRCh38{'; from ' + a.source if a.source else ''})")
    status, lines = site_genotype(gvcf, chrom, pos, a.ref)
    for line in lines:
        print(f"  {line}")
    genes = set()
    rsid = a.query.lower() if re.fullmatch(r"(?i)rs\d+", a.query) else None
    if ann is None:
        print("  annotation: no annotated VCF (step 13), so no consequence, ClinVar or gnomAD")
    else:
        recs = [r for r in ann.records if r.chrom == chrom and r.pos <= pos <= r.end]
        if not recs:
            print("  annotation: no variant of this sample here in the annotated VCF (it holds only the sample's "
                  "variant calls); ClinVar and gnomAD are not read for a reference or no-call site")
        for r in recs:
            bases, phased = r.bases()
            print(f"  variant: {contig_name(r.chrom)}:{r.pos} {r.ref}>{','.join(r.alts)} "
                  f"{fmt_gt(bases, phased)} ({zygosity(bases)})")
            for line in annotation_lines(r, ann.csq, indent="    "):
                print(line)
            genes |= symbols(r, ann.csq)
            rsid = rsid or next((i.lower() for i in sorted(r.ids(ann.csq)) if i.lower().startswith("rs")), None)
    print("  in the reports:")
    for line in summary_lines(a.summary, genes, rsid):
        print(line)
    return status


def report_gene(a, gvcf, ann):
    chrom, start, end = norm_chrom(a.chrom), a.start, a.end
    gene = a.query
    print(f"== {gene}: {contig_name(chrom)}:{start}-{end} (GRCh38; {a.source or 'gene span'})")
    span = end - start + 1
    covered = 0
    for r in gvcf.records:
        if r.chrom != chrom or r.end < start or r.pos > end:
            continue
        idx, _ = r.gt()
        if idx and any(i is not None for i in idx):
            covered += min(r.end, end) - max(r.pos, start) + 1
    covered = min(covered, span)
    print(f"  gVCF: {covered:,} of {span:,} bases have a genotype call ({100.0 * covered / span:.1f}%); "
          "a variant can only be missing from the list below where there is no call")
    listed = 0
    if ann is not None:
        print("  variants called in the gene (annotated VCF):")
        for r in ann.records:
            if r.chrom != chrom or r.end < start or r.pos > end:
                continue
            syms = {s.upper() for s in symbols(r, ann.csq)}
            if ann.csq and gene.upper() not in syms:
                continue
            bases, phased = r.bases()
            ids = ",".join(sorted(r.ids(ann.csq))) or "no rsID"
            flag = f", filter {r.filter}" if filtered(r) else ""
            print(f"  - {contig_name(r.chrom)}:{r.pos} {r.ref}>{','.join(r.alts)} {fmt_gt(bases, phased)} "
                  f"({zygosity(bases)}; {ids}; {depth_text(r)}{flag})")
            for line in annotation_lines(r, ann.csq, gene=gene, indent="      "):
                print(line)
            listed += 1
    else:
        print("  variants called in the gene (gVCF; no annotated VCF from step 13, so no consequence, ClinVar or gnomAD):")
        for r in gvcf.records:
            if r.chrom != chrom or r.end < start or r.pos > end or r.is_block() or not r.has_alt_call():
                continue
            bases, phased = r.bases()
            flag = f", filter {r.filter}" if filtered(r) else ""
            print(f"  - {contig_name(r.chrom)}:{r.pos} {r.ref}>{','.join(r.real_alts())} {fmt_gt(bases, phased)} "
                  f"({zygosity(bases)}; {depth_text(r)}{flag})")
            listed += 1
    if not listed:
        print("  - none: the sample matches the reference wherever the gene has a call")
    print("  in the reports:")
    for line in summary_lines(a.summary, {gene}, None):
        print(line)
    return "listed"


def report(a):
    region = (a.chrom, a.start, a.end)
    gvcf = Vcf(a.gvcf, *region)
    ann = Vcf(a.annotated, *region) if a.annotated else None
    if a.kind == "site":
        report_site(a, gvcf, ann)
    else:
        report_gene(a, gvcf, ann)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("resolve")
    r.add_argument("--query", required=True)
    r.add_argument("--sites", help="the consumer site table (rsid, chrom, pos38, ref38)")
    r.add_argument("--table", action="append", help="a pipeline table with rsid and chrom/pos (repeatable)")
    r.add_argument("--genes", help="GENCODE gene lines (gencode.v*.basic.genes.gtf)")
    r.add_argument("--annotated", help="the annotated VCF (step 13 or 30)")
    p = sub.add_parser("report")
    p.add_argument("--query", required=True)
    p.add_argument("--kind", choices=("site", "gene"), required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--ref", default=".")
    p.add_argument("--source", default="")
    p.add_argument("--gvcf", required=True)
    p.add_argument("--annotated")
    p.add_argument("--summary")
    a = ap.parse_args(argv)
    try:
        return resolve(a) if a.cmd == "resolve" else report(a)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
