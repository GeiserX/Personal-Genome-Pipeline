#!/usr/bin/env python3
"""results_index.py: every gene symbol and rsID one run's results mention, for search.

  results_index.py build --summary summary.json [--sample-dir DIR] -o INDEX.json
  results_index.py query --index INDEX.json QUERY

Step 24 (scripts/24-html-report.sh) and the Nextflow HTML_REPORT process write
the index as <sample>_search_index.json beside the report, through
render_report.py --search-index. The report closes with an appendix built
from it, "Genes and rsIDs in these results"; the live site's search box
(bin/owngenome_live_search.py, assets/live/plugins/search.js) loads the file once
and filters it in the browser, so a query never leaves the page.

What it reads:
  - every section of summary.json in state ok or stale: each dict in its data
    that names a gene or an rsID, found by bin/genome_lookup.py's names_of(),
    the same rule owngenome-lookup.sh uses for "in the reports" (so the search and a
    lookup cannot disagree);
  - the gene column of step 23 (clinical filter) and step 31 (slivar), whose
    sections keep counts only: from the summary table on the bash path
    (clinical/<s>_clinical_summary.tsv GENE, slivar/<s>_slivar_summary.tsv
    SYMBOL and Existing_variation), else from the VCF the Nextflow report
    stages (INFO/PGP_GENE of the clinical VCF, the CSQ entries
    collect_summary.csq_rank picks in the prioritized VCF). Both forms give
    the same entry for the same variant when its most severe consequence is
    on the transcript csq_rank picks; the slivar TSV keeps the gene of the
    most severe one (split-vep -s worst), so the two can differ there.
A section that is missing, failed or unreadable gives no entry.

The index:
  {schema_version: 1, sample, generated_utc,
   terms: {TERM: [{section, anchor, title, step, state, source, count,
                   mentions: [{text, locus}]}]}}
TERM is a gene symbol in upper case or an rsID in lower case. One item per
section that mentions the term; count is how many variants (a mention with
a locus counts once per locus) or other mentions it has there, mentions the
first MENTIONS_KEPT of them. anchor is the section's id in the HTML
report (#anchor), or null for a key no report section shows.

A query (query(), and the same rule in search.js) is an rsID, a gene symbol
(a prefix also lists the terms it starts: BRCA gives BRCA1 and BRCA2) or a
position, which the index cannot answer: a position is a genotype question,
for scripts/owngenome-lookup.sh.

Standard library only.
"""
import argparse
import csv
import json
import os
import re
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import collect_summary  # noqa: E402
import genome_lookup  # noqa: E402

SCHEMA_VERSION = 1
# Mentions kept per term and section; count says how many there were.
MENTIONS_KEPT = 20
# Terms listed for a prefix query.
SUGGESTIONS = 20
POSITION_HINT = ("The search finds gene symbols and rsIDs in these results. A position is a question about your "
                 "genotype: run scripts/owngenome-lookup.sh <sample> {q}, which reads it from the gVCF on this computer.")


def walk(obj, blocked=frozenset()):
    """(dict, new names) for every dict of OBJ that names a gene or rsID no
    dict above it named: what genome_lookup.walk_mentions yields for each
    of those names."""
    if isinstance(obj, dict):
        genes, rsids = genome_lookup.names_of(obj)
        names = genes | rsids
        if names - blocked:
            yield obj, names - blocked
        for v in obj.values():
            yield from walk(v, blocked | names)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk(v, blocked)


def locus_of(d):
    """chrom:pos of a mention that carries them, else None."""
    c, p = d.get("chrom"), d.get("pos")
    if isinstance(c, str) and c and (isinstance(p, int) or (isinstance(p, str) and p.isdigit())):
        return f"{c}:{p}"
    m = re.match(r"^((?:chr)?[0-9XYMT]{1,2}):([0-9]+)\b", str(d.get("variant") or ""))
    return f"{m.group(1)}:{m.group(2)}" if m else None


def rs_ids(value):
    """The rsIDs of an Existing_variation value ('rs123&COSV1' -> ['rs123'])."""
    return sorted({x.lower() for x in re.split(r"[&,]", value or "") if genome_lookup.RSID_RE.match(x)})


def variant_text(chrom, pos, ref, alt):
    return f"{chrom}:{pos} {ref}>{alt}"


def clinical_names(path):
    """[(text, locus, genes, rsids)] of step 23's variants: the GENE column of
    the summary table, else INFO/PGP_GENE of the clinical VCF (the gene the
    table's GENE column is written from), else the genes of the CSQ entries
    csq_rank picks."""
    out = []
    if path.endswith(".tsv"):
        for r in collect_summary.read_tsv(path):
            g = r.get("GENE") or ""
            out.append((variant_text(r.get("CHROM"), r.get("POS"), r.get("REF"), r.get("ALT")),
                        f"{r.get('CHROM')}:{r.get('POS')}", {g.upper()} - {"", "."}, set()))
        return out
    return vcf_names(path, gene_tag="PGP_GENE", with_rsids=False)


def slivar_names(path):
    """[(text, locus, genes, rsids)] of step 31's prioritized variants: SYMBOL
    and Existing_variation of the summary table, else the CSQ of the VCF."""
    out = []
    if path.endswith(".tsv"):
        for r in collect_summary.read_tsv(path):
            g = r.get("SYMBOL") or ""
            out.append((variant_text(r.get("CHROM"), r.get("POS"), r.get("REF"), r.get("ALT")),
                        f"{r.get('CHROM')}:{r.get('POS')}", {g.upper()} - {"", "."},
                        set(rs_ids(r.get("Existing_variation")))))
        return out
    return vcf_names(path, gene_tag=None, with_rsids=True)


def vcf_names(path, gene_tag, with_rsids):
    fmt, out = None, []
    with collect_summary.open_text(path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r"Format: ([^\"]+)", line)
                fmt = m.group(1).split("|") if m else None
                continue
            if line.startswith("#"):
                continue
            r = line.rstrip("\n").split("\t")
            if len(r) < 8:
                continue
            info = collect_summary.info_map(r[7])
            csq = info.get("CSQ", "")
            genes = set()
            if gene_tag and info.get(gene_tag) not in (None, "", "."):
                genes = {info[gene_tag].upper()}
            elif fmt:
                _, entries, _ = collect_summary.csq_rank(fmt, csq)
                genes = {e["gene"].upper() for e in entries if e["gene"]}
            rsids = set()
            if with_rsids and fmt and "Existing_variation" in fmt:
                i = fmt.index("Existing_variation")
                for t in csq.split(","):
                    v = t.split("|")
                    if len(v) == len(fmt):
                        rsids.update(rs_ids(v[i]))
            out.append((variant_text(r[0], r[1], r[3], r[4]), f"{r[0]}:{r[1]}", genes, rsids))
    return out


# The sections whose genes are in their step's table, not in summary.json.
TABLES = {"clinical": clinical_names, "slivar": slivar_names}


def build(summary, sample_dir=None, anchors=None):
    """The index of SUMMARY (a summary.json dict). SAMPLE_DIR, when given, is
    where the step 23 and 31 tables are read from (the paths summary.json
    records); ANCHORS maps a section key to its report anchor."""
    anchors = anchors or {}
    terms = {}

    def add(term, key, sec, text, locus):
        items = terms.setdefault(term, [])
        item = next((x for x in items if x["section"] == key), None)
        if item is None:
            item = {"section": key, "anchor": anchors.get(key), "title": sec.get("title") or key,
                    "step": sec.get("step"), "state": sec.get("state"), "source": sec.get("source"),
                    "count": 0, "mentions": []}
            items.append(item)
        # One mention per variant (locus) or per text: the ClinVar list, the
        # secondary findings and a step's table can name the same variant.
        seen = item.setdefault("_seen", set())
        mark = ("locus", locus) if locus else ("text", text)
        if mark in seen:
            return
        seen.add(mark)
        item["count"] += 1
        if len(item["mentions"]) < MENTIONS_KEPT:
            item["mentions"].append({"text": text, "locus": locus})

    for key, sec in (summary.get("sections") or {}).items():
        if sec.get("state") not in ("ok", "stale"):
            continue
        for d, names in walk(sec.get("data") or {}):
            text, locus = genome_lookup.flat(d, 160), locus_of(d)
            for n in sorted(names):
                add(n, key, sec, text, locus)
        reader, src = TABLES.get(key), sec.get("source")
        if reader and src and sample_dir:
            path = src if os.path.isabs(src) else os.path.join(sample_dir, src)
            if os.path.isfile(path):
                # A table that cannot be read (truncated gzip, a bad row) adds
                # nothing, and the report is still written.
                try:
                    records = reader(path)
                except (OSError, EOFError, csv.Error, zlib.error) as e:
                    print(f"[WARN] {key}: cannot read {src} for the search index: {e}", file=sys.stderr)
                    records = []
                for text, locus, genes, rsids in records:
                    for n in sorted(genes | rsids):
                        add(n, key, sec, text, locus)
    for items in terms.values():
        for it in items:
            del it["_seen"]
    return {"schema_version": SCHEMA_VERSION, "sample": summary.get("sample"),
            "generated_utc": summary.get("generated_utc"),
            "terms": {k: terms[k] for k in sorted(terms)}}


def query(index, q):
    """{kind, term, items, suggestions, hint} for one query; kind is rsid,
    gene, site or error. A gene query that is no term, or a prefix of others,
    lists up to SUGGESTIONS terms that start with it."""
    try:
        kind, val = genome_lookup.parse_query(q)
    except ValueError as e:
        return {"kind": "error", "term": None, "items": [], "suggestions": [], "hint": str(e)}
    terms = index.get("terms") or {}
    if kind == "site":
        return {"kind": "site", "term": None, "items": [], "suggestions": [],
                "hint": POSITION_HINT.format(q=q.strip())}
    term = val.lower() if kind == "rsid" else val.upper()
    sugg = [] if kind == "rsid" else [t for t in terms if t.startswith(term) and t != term][:SUGGESTIONS]
    return {"kind": kind, "term": term, "items": terms.get(term, []), "suggestions": sugg, "hint": None}


def query_lines(res):
    if res["kind"] in ("error", "site"):
        return [res["hint"]]
    out = []
    if res["items"]:
        out.append(f"{res['term']} is in {len(res['items'])} section(s):")
        for it in res["items"]:
            stale = " (stale: from an earlier run)" if it["state"] == "stale" else ""
            out.append(f"  {it['title']} (step {it['step']}){stale}: {it['count']} mention(s)")
            out += [f"    - {m['text']}" for m in it["mentions"]]
            if it["count"] > len(it["mentions"]):
                out.append(f"    - ... and {it['count'] - len(it['mentions'])} more in {it['source']}")
    else:
        out.append(f"none of these results mentions {res['term']}")
    if res["suggestions"]:
        out.append("terms that start with it: " + ", ".join(res["suggestions"]))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--summary", required=True)
    b.add_argument("--sample-dir", help="the sample folder, for the step 23 and 31 tables")
    b.add_argument("-o", "--out", required=True)
    q = sub.add_parser("query")
    q.add_argument("--index", required=True)
    q.add_argument("query")
    a = ap.parse_args(argv)
    if a.cmd == "build":
        import render_report  # the report's section anchors; render_report imports this module
        with open(a.summary) as f:
            s = json.load(f)
        idx = build(s, a.sample_dir, render_report.SECTION_OF)
        with open(a.out, "w") as f:
            json.dump(idx, f, indent=1)
            f.write("\n")
        print(f"{len(idx['terms'])} terms")
        return 0
    with open(a.index) as f:
        idx = json.load(f)
    res = query(idx, a.query)
    for line in query_lines(res):
        print(line)
    return 0 if res["kind"] in ("rsid", "gene") else 2


if __name__ == "__main__":
    sys.exit(main())
