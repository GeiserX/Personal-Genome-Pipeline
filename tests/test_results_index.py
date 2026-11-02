#!/usr/bin/env python3
"""bin/results_index.py, its report appendix and the live search on a synthetic sample.

The sample (tests/fixtures/results_index/SYNTH, invented) has ClinVar hits in
CFTR and BRCA2, a CPIC table (CYP2C19, CYP2D6), and the tables of steps 23 and
31 in both forms: the summary TSVs the bash steps write and the VCFs the
Nextflow report stages.

Checks:
  1. every term the index finds in summary.json, genome_lookup.summary_lines
     finds too, through the one names_of() rule: with a GENE_KEYS entry
     removed both lose CYP2C19 (so the check can fail);
  2. the step 23 and 31 genes come from the TSV on the bash path and from the
     VCF on the Nextflow path, and both give the same entries; with the VCF
     reader removed the Nextflow form has none (red);
  3. Existing_variation 'rs123&COSV1' gives rs123 only, and query('RS123')
     finds it;
  4. a stale section's entries say stale; a failed or missing one gives none;
  5. query('BRCA') lists BRCA1 and BRCA2; a position gets the pgp-lookup.sh
     hint; nothing else is a query;
  6. the HTML report closes with one line per term, each link resolves to a
     section of the page, a term that is missing from the index is missing
     from the appendix (red), and the page has no script;
  7. the command line builds and queries the index (render_report.py
     --search-index writes the same index);
  8. the live search: pgp_live_search.py registers the file glob and its
     state hook; search.js's match() under node answers each fixture query
     as query() does, and has no request in it besides the one index fetch.

Run: python3 tests/test_results_index.py
"""
import copy
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402
import genome_lookup  # noqa: E402
import pgp_live_search  # noqa: E402
import render_report  # noqa: E402
import results_index  # noqa: E402

FIXTURE = os.path.join(REPO, "tests", "fixtures", "results_index", "SYNTH")
SEARCH_JS = os.path.join(REPO, "assets", "live", "plugins", "search.js")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:900]}")
    if not ok:
        FAILS += 1


def sample(work, form):
    """A copy of the fixture with the step 23 and 31 tables in one FORM:
    'tsv' (bash: summary TSVs, VCFs too) or 'vcf' (Nextflow: VCFs only)."""
    d = os.path.join(work, form, "SYNTH")
    shutil.copytree(FIXTURE, d)
    for v in ("clinical/SYNTH_clinical.vcf", "slivar/SYNTH_prioritized.vcf"):
        with open(os.path.join(d, v), "rb") as f, gzip.open(os.path.join(d, v + ".gz"), "wb") as g:
            g.write(f.read())
        os.remove(os.path.join(d, v))
    if form == "vcf":
        os.remove(os.path.join(d, "clinical/SYNTH_clinical_summary.tsv"))
        os.remove(os.path.join(d, "slivar/SYNTH_slivar_summary.tsv"))
    return d


def lookup_says(summary_path, term):
    """genome_lookup.summary_lines for one gene or rsID, as one string."""
    rs = term if term.startswith("rs") else None
    return "\n".join(genome_lookup.summary_lines(summary_path, set() if rs else {term}, rs))


def table_entries(idx, key):
    return {t: [m["text"] for it in items if it["section"] == key for m in it["mentions"]]
            for t, items in idx["terms"].items() if any(it["section"] == key for it in items)}


def main():
    work = tempfile.mkdtemp()
    try:
        tsv, vcf = sample(work, "tsv"), sample(work, "vcf")
        s_tsv, s_vcf = collect_summary.collect("SYNTH", tsv), collect_summary.collect("SYNTH", vcf)
        check("fixture: the clinical and slivar sections read the TSV on the bash form and the VCF on the other",
              (s_tsv["sections"]["clinical"]["source"], s_vcf["sections"]["clinical"]["source"],
               s_tsv["sections"]["slivar"]["source"], s_vcf["sections"]["slivar"]["source"])
              == ("clinical/SYNTH_clinical_summary.tsv", "clinical/SYNTH_clinical.vcf.gz",
                  "slivar/SYNTH_slivar_summary.tsv", "slivar/SYNTH_prioritized.vcf.gz"),
              {k: s_tsv["sections"][k]["source"] for k in ("clinical", "slivar")})
        idx = results_index.build(s_tsv, tsv, render_report.SECTION_OF)
        T = idx["terms"]
        check("index: the genes and rsIDs of the fixture, genes upper-case, rsIDs lower-case",
              sorted(T) == ["BRCA1", "BRCA2", "CFTR", "CYP2C19", "CYP2D6", "rs1000002", "rs123"], sorted(T))
        check("index: shape (schema_version 1, the sample, one item per section with its report anchor)",
              idx["schema_version"] == 1 and idx["sample"] == "SYNTH"
              and {(it["section"], it["anchor"]) for it in T["BRCA2"]}
              == {("clinvar", "clinvar"), ("carrier", "carrier-list"), ("clinical", "clinical-filter"),
                  ("slivar", "slivar")}, T["BRCA2"])
        check("index: one mention per variant in a section (the secondary findings and the table name chr13:2000)",
              [(it["section"], it["count"]) for it in T["BRCA2"]]
              == [("clinvar", 1), ("carrier", 1), ("clinical", 1), ("slivar", 1)], T["BRCA2"])

        # 1. the same mentions as pgp-lookup.sh
        spath = os.path.join(work, "summary.json")
        with open(spath, "w") as f:
            json.dump(s_tsv, f)
        from_summary = sorted({t for t, items in T.items() if any(it["section"] not in results_index.TABLES
                                                                  for it in items)})
        missed = [t for t in from_summary if "none of the reports mentions" in lookup_says(spath, t)]
        check("lookup agrees: every term the index finds in summary.json, summary_lines finds too",
              from_summary == ["BRCA2", "CFTR", "CYP2C19", "CYP2D6"] and not missed, (from_summary, missed))
        for t in from_summary:
            secs = {s_tsv["sections"][it["section"]]["title"] for it in T[t] if it["section"] not in results_index.TABLES}
            said = lookup_says(spath, t)
            check(f"lookup agrees on {t}: the same sections", all(x in said for x in secs), (secs, said))
        # A ClinVar hit's rsids (collect_summary puts them there from VEP): found by both.
        withrs = copy.deepcopy(s_tsv)
        withrs["sections"]["clinvar"]["data"]["hits"][0]["rsids"] = ["rs7000001"]
        rpath = os.path.join(work, "summary-rs.json")
        with open(rpath, "w") as f:
            json.dump(withrs, f)
        ri = results_index.build(withrs, tsv, render_report.SECTION_OF)["terms"]
        check("a ClinVar hit's rsids: the index and the lookup both find it in the ClinVar section",
              [it["section"] for it in ri.get("rs7000001", [])] == ["clinvar"]
              and "ClinVar screen (step 06)" in lookup_says(rpath, "rs7000001"),
              (ri.get("rs7000001"), lookup_says(rpath, "rs7000001")))
        check("red control: without the rsids the lookup finds nothing for it",
              "none of the reports mentions rs7000001" in lookup_says(spath, "rs7000001"))
        keep = genome_lookup.GENE_KEYS
        genome_lookup.GENE_KEYS = tuple(k for k in keep if k != "gene")
        try:
            gone = results_index.build(s_tsv, tsv, render_report.SECTION_OF)["terms"]
            said = lookup_says(spath, "CYP2C19")
        finally:
            genome_lookup.GENE_KEYS = keep
        check("red control: without the 'gene' key both the index and the lookup lose CYP2C19",
              "CYP2C19" not in gone and "none of the reports mentions CYP2C19" in said, (sorted(gone), said))

        # 2. bash and Nextflow forms
        nf = results_index.build(s_vcf, vcf, render_report.SECTION_OF)
        for key in results_index.TABLES:
            a, b = table_entries(idx, key), table_entries(nf, key)
            check(f"{key}: the TSV form and the VCF form give the same entries", a == b and a, (a, b))
        check("clinical: BRCA2 and BRCA1 from the table; slivar: CFTR, BRCA2 and both rsIDs",
              sorted(table_entries(idx, "clinical")) == ["BRCA1", "BRCA2"]
              and sorted(table_entries(idx, "slivar")) == ["BRCA2", "CFTR", "rs1000002", "rs123"],
              (table_entries(idx, "clinical"), table_entries(idx, "slivar")))
        keep = results_index.vcf_names
        results_index.vcf_names = lambda path, gene_tag, with_rsids: []
        try:
            red = results_index.build(s_vcf, vcf, render_report.SECTION_OF)
        finally:
            results_index.vcf_names = keep
        check("red control: without the VCF reader the Nextflow form has no table entries, so the check above can fail",
              not table_entries(red, "slivar") and table_entries(nf, "slivar"), table_entries(red, "slivar"))
        bash_like = copy.deepcopy(nf)
        nfc = copy.deepcopy(idx)
        for x in (bash_like, nfc):
            x["generated_utc"] = None
            for items in x["terms"].values():
                for it in items:
                    it["source"] = None
        check("bash and Nextflow: the whole index is the same, the file names aside", bash_like["terms"] == nfc["terms"],
              (bash_like["terms"], nfc["terms"]))

        # 3. Existing_variation
        check("Existing_variation: 'rs123&COSV1' gives rs123 only", results_index.rs_ids("rs123&COSV1") == ["rs123"]
              and "COSV1" not in T and "cosv1" not in T)
        r = results_index.query(idx, "RS123")
        check("query('RS123') finds rs123 in variant prioritization", r["kind"] == "rsid" and r["term"] == "rs123"
              and [it["section"] for it in r["items"]] == ["slivar"]
              and r["items"][0]["mentions"][0]["locus"] == "chr7:3000", r)

        # 4. states
        st = copy.deepcopy(s_tsv)
        st["sections"]["cpic"]["state"] = "stale"
        st["sections"]["mito"].update(state="failed", data={"variants": [{"gene": "MT-ND4", "variant": "chrM:11778 G>A"}]})
        st["sections"]["lpa"].update(state="missing", data={"variants": [{"rsid": "rs10455872", "gene": "LPA"}]})
        st["sections"]["haplogroup"].update(state="unreadable", data={"gene": "MT-CYB"})
        si = results_index.build(st, tsv, render_report.SECTION_OF)["terms"]
        check("states: a stale section's entries say stale", [it["state"] for it in si["CYP2D6"]] == ["stale"], si["CYP2D6"])
        check("states: a failed, missing or unreadable section gives no entry",
              not {"MT-ND4", "LPA", "rs10455872", "MT-CYB"} & set(si), sorted(si))
        ok = copy.deepcopy(st)
        ok["sections"]["mito"]["state"] = "ok"
        check("red control: the same data in an ok section is indexed",
              "MT-ND4" in results_index.build(ok, tsv, render_report.SECTION_OF)["terms"])

        # 5. queries
        r = results_index.query(idx, "brca")
        check("query('BRCA') lists BRCA1 and BRCA2", r["kind"] == "gene" and not r["items"]
              and r["suggestions"] == ["BRCA1", "BRCA2"], r)
        r = results_index.query(idx, "brca2")
        check("query('brca2') finds BRCA2 in four sections", [it["section"] for it in r["items"]]
              == ["clinvar", "carrier", "clinical", "slivar"] and r["suggestions"] == [], r)
        r = results_index.query(idx, "chr7:1,000")
        check("a position gets the pgp-lookup.sh hint and no entries", r["kind"] == "site" and not r["items"]
              and "pgp-lookup.sh" in r["hint"], r)
        r = results_index.query(idx, "TP53")
        check("a gene the results do not name: no entries", r["kind"] == "gene" and not r["items"], r)
        check("not a query: an error", results_index.query(idx, "two words")["kind"] == "error")

        # 6. the appendix
        page = render_report.html_report(s_tsv, idx)
        ids = set(re.findall(r'id="([^"]+)"', page))
        app = page[page.index('id="ch-terms"'):page.index('<div class="disclaimer">')]
        rows = re.findall(r"<li><strong>([^<]+)</strong>: (.*)</li>", app)
        hrefs = re.findall(r'href="#([^"]+)"', app)
        check("appendix: one line per term, in the index's order", [t for t, _ in rows] == list(T), rows)
        check("appendix: every link resolves to an id of the page", hrefs and all(h in ids for h in hrefs),
              [h for h in hrefs if h not in ids])
        check("appendix: BRCA2 links its four sections", re.search(
            r'<strong>BRCA2</strong>: <a href="#clinvar">[^<]+</a>, <a href="#carrier-list">[^<]+</a>, '
            r'<a href="#clinical-filter">[^<]+</a>, <a href="#slivar">[^<]+</a></li>', app) is not None, app[:2000])
        check("appendix: in the contents", '<a href="#ch-terms">Genes and rsIDs in these results</a>' in page)
        check("the report has no script", "<script" not in page.lower())
        less = copy.deepcopy(idx)
        less["terms"].pop("CYP2D6")
        check("red control: a term left out of the index is left out of the appendix",
              "<strong>CYP2D6</strong>" in page and "<strong>CYP2D6</strong>" not in render_report.html_report(s_tsv, less))
        bare = render_report.html_report(s_tsv)
        check("without an index the appendix is built from the summary alone (no step 23 or 31 genes)",
              "<strong>CFTR</strong>" in bare and "<strong>BRCA1</strong>" not in bare, "")

        # 7. command lines
        out = os.path.join(work, "idx.json")
        p = subprocess.run([sys.executable, os.path.join(REPO, "bin", "results_index.py"), "build", "--summary", spath,
                            "--sample-dir", tsv, "-o", out], capture_output=True, text=True)
        with open(out) as f:
            cli = json.load(f)
        check("CLI build: the same terms as build()", p.returncode == 0 and cli["terms"] == T, p.stderr)
        q = subprocess.run([sys.executable, os.path.join(REPO, "bin", "results_index.py"), "query", "--index", out,
                            "rs123"], capture_output=True, text=True)
        check("CLI query rs123: exit 0, names the section", q.returncode == 0 and "Variant prioritization" in q.stdout,
              q.stdout + q.stderr)
        q = subprocess.run([sys.executable, os.path.join(REPO, "bin", "results_index.py"), "query", "--index", out,
                            "7:117559590"], capture_output=True, text=True)
        check("CLI query of a position: exit 2 with the hint", q.returncode == 2 and "pgp-lookup.sh" in q.stdout, q.stdout)
        rr = os.path.join(work, "rr.json")
        with redirect_stdout(io.StringIO()):
            render_report.main(["--sample", "SYNTH", "--sample-dir", tsv, "--json", os.path.join(work, "s2.json"),
                                "-o", os.path.join(work, "r.html"), "--search-index", rr])
        with open(rr) as f:
            check("render_report.py --search-index writes the same terms", json.load(f)["terms"] == T)

        # 8. the live search
        run = os.path.join(work, "SYNTH")
        os.makedirs(run)
        check("live: no index file, no state", pgp_live_search.search_state(run, {}, []) == {})
        shutil.copy(rr, os.path.join(run, "SYNTH_search_index.json"))
        st8 = pgp_live_search.search_state(run, {}, [])
        check("live: the state names the file and its mtime", st8.get("file") == "SYNTH_search_index.json"
              and isinstance(st8.get("mtime"), int), st8)

        class App:
            globs, hooks = [], {}

            def run_file(self, g):
                self.globs.append(g)

            def state_hook(self, name, fn):
                self.hooks[name] = fn
        app_ = App()
        pgp_live_search.register(app_)
        check("live: the plugin serves *_search_index.json and hooks plugins.search",
              app_.globs == ["*_search_index.json"] and app_.hooks.get("search") is pgp_live_search.search_state)
        js = open(SEARCH_JS).read()
        check("search.js: one fetch (the index), no XMLHttpRequest, no sendBeacon, no WebSocket, no eval",
              js.count("fetch(") == 1 and not re.search(r"XMLHttpRequest|sendBeacon|WebSocket|eval\(|innerHTML", js))
        queries = ["RS123", "rs1000002", "brca", "BRCA2", "cyp2", "TP53", "chr7:1,000", "13-2000", "chr1:0", "two words", ""]
        node = shutil.which("node")
        check("node is on the PATH (ubuntu-latest has it)", node is not None)
        if node:
            prog = ("const m=require(process.argv[1]).match;const idx=JSON.parse(require('fs').readFileSync(process.argv[2]));"
                    "const qs=JSON.parse(process.argv[3]);console.log(JSON.stringify(qs.map(q=>{const r=m(idx,q);"
                    "return [r.kind,r.term,r.items.map(i=>i.section),r.suggestions]})))")
            p = subprocess.run([node, "-e", prog, SEARCH_JS, rr, json.dumps(queries)], capture_output=True, text=True)
            got = json.loads(p.stdout) if p.returncode == 0 else p.stderr
            want = []
            for q in queries:
                r = results_index.query(idx, q)
                want.append([r["kind"], r["term"], [i["section"] for i in r["items"]], r["suggestions"]])
            check("search.js match() answers each fixture query as query() does", got == want, (got, want))
            check("search.js: rs123 lists variant prioritization", got[0] == ["rsid", "rs123", ["slivar"], []], got)
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
