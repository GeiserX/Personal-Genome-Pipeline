#!/usr/bin/env python3
"""render_report.py: the text and HTML reports, both from one summary.

  render_report.py --sample S --sample-dir DIR [--json OUT.json] [-o REPORT.txt] [-o REPORT.html]
  render_report.py --summary IN.json [-o REPORT.txt] [-o REPORT.html]

-o writes the text report for a .txt name and the HTML report for .html.

With --sample-dir it first collects the summary (bin/collect_summary.py) and
writes it to --json (default DIR/summary.json). Every number in both reports
comes from that summary, so the two cannot disagree. A section whose result is
from an earlier run (state "stale") is shown with that date and a STALE mark.
Standard library only.
"""
import argparse
import html
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import collect_summary  # noqa: E402

E = html.escape


def stale_note(sec):
    if sec["state"] == "stale":
        return f"STALE: {sec['note']}"
    if sec["state"] == "unreadable":
        return f"could not be read: {sec['note']}"
    return ""


# --- text ------------------------------------------------------------------------

def text_report(s):
    S = s["sections"]
    out = []
    w = out.append
    w("=" * 80)
    w("  GENOMICS ANALYSIS REPORT")
    w(f"  Sample: {s['sample']}")
    w(f"  Generated: {s['generated_utc']}")
    if s["run"]["started_utc"]:
        w(f"  Latest run-all.sh run: {s['run']['started_utc']}")
    w("=" * 80)
    w("")

    def head(key, title=None):
        sec = S[key]
        if sec["state"] == "missing":
            return False
        w(f"## {title or sec['title']}")
        w("---")
        note = stale_note(sec)
        if note:
            w(f"  [{note}]")
        return sec["state"] in ("ok", "stale")

    # QC
    cov, sexc = S["coverage"], S["sex_check"]
    if cov["state"] != "missing" or sexc["state"] != "missing":
        w("## Quality control")
        w("---")
        for sec in (cov, sexc):
            if stale_note(sec):
                w(f"  [{sec['title']}: {stale_note(sec)}]")
        md = cov["data"].get("mean_depth")
        w(f"  Mean depth: {md:.1f}x" if isinstance(md, (int, float)) else "  Mean depth: not available")
        inferred = sexc["data"].get("inferred_sex")
        declared = s["run"]["declared_sex"]
        w(f"  Inferred sex: {inferred or 'not available'}")
        w(f"  Declared sex: {declared or 'not declared'}"
          + ("  ** DOES NOT MATCH the inferred sex **" if inferred and declared and inferred != declared else ""))
        w("")

    if head("variants", "Variant Calling (DeepVariant)"):
        d = S["variants"]["data"]
        w(f"  Total variants: {d['total']}")
        w(f"  PASS variants:  {d['pass']}")
        w(f"  SNPs: {d['snps']}  Indels: {d['indels']}")
        w("")

    if head("clinvar", "ClinVar Pathogenic Screen"):
        d = S["clinvar"]["data"]
        w(f"  Pathogenic/Likely Pathogenic hits: {d['count']}")
        w(f"  ClinVar file date: {s['databases']['clinvar_release'] or 'unknown'}")
        if d["count"]:
            by = d["by_stars"]
            w("  By review status: " + ", ".join(f"{k} star{'s' if k != '1' else ''} {by[k]}" for k in ("4", "3", "2", "1", "0")))
            w("")
            w("  Hits (best-reviewed first, up to 50): gene, position, genotype, significance [review status] stars")
            for h in d["hits"][:50]:
                w(f"    {h['gene']:<10} {h['chrom']}:{h['pos']} {h['ref']}>{h['alt']}  {h['genotype']}  "
                  f"{h['significance']} [{h['review_status']}] {h['stars']}*")
        w("")

    if head("pharmcat", "Pharmacogenomics (PharmCAT)"):
        d = S["pharmcat"]["data"]
        w(f"  Report: {S['pharmcat']['source']} (PharmCAT {d.get('version') or 'version unknown'})")
        if d.get("html_report"):
            w("  HTML report beside it (open in a browser for full details)")
        w("")

    if head("cpic", "CPIC Drug-Gene Recommendations"):
        d = S["cpic"]["data"]
        if d["parse_failed"]:
            w("  The PharmCAT report could not be parsed: see the CPIC report.")
        w(f"  Genes with a non-normal phenotype: {d['non_normal']}  More than one possible result: "
          f"{d.get('ambiguous', 0)}  Not called: {d['not_called']}")
        for g in d["genes"]:
            if g["status"] in ("non-normal", "ambiguous"):
                w(f"    {g['gene']:<10} {g['diplotype']:<28} {g['phenotype']}")
        for warn in d["warnings"]:
            w(f"  WARNING: {warn}")
        w(f"  Drugs per gene: {os.path.join(os.path.dirname(S['cpic']['source']), s['sample'] + '_cpic_recommendations.txt')}")
        w("")

    if head("pypgx", "pypgx Pharmacogenomics"):
        d = S["pypgx"]["data"]
        w(f"  Genes called: {d['genes_called']}/{d['genes_total']}")
        w(f"  CYP2D6 diplotype: {d['cyp2d6']}")
        if "comparison" in d:
            w(f"  PharmCAT vs pypgx conflicts: {d['comparison']['conflicts']}")
            if d["comparison"]["one_tool_only"]:
                w(f"  Scope differences (only one tool called): {d['comparison']['one_tool_only']}")
        w("")

    calls = s["cyp2d6"]["calls"]
    if any(calls.values()):
        w("## CYP2D6 across callers")
        w("---")
        for k, v in calls.items():
            w(f"  {k:<9} {v or 'not run'}")
        agree = s["cyp2d6"]["agree"]
        w("  Agreement: " + ("yes" if agree else "NO, review before acting on CYP2D6" if agree is False
                              else "fewer than two callers have a call"))
        w("")

    if head("hla", "HLA Typing (T1K)"):
        for l in S["hla"]["data"]["loci"]:
            w(f"    {l['gene']:<8} {' / '.join(l['alleles']) or 'no call'}")
        w(f"  HLA database: {s['databases']['hla_database'] or 'unknown'}")
        w("")

    if head("prs", "Polygenic Risk Scores"):
        for r in S["prs"]["data"]["scores"]:
            w(f"  {r['condition']:<35} {r['score']} ({r['matched']}/{r['total']} variants matched) {r['pgs_id']}")
        w("  NOTE: Raw PRS scores are NOT interpretable without an ancestry-matched")
        w("  reference panel; hom-ref sites are absent from the VCF. See docs/25-prs.md.")
        w("")

    if head("manta", "Structural Variants (Manta)"):
        d = S["manta"]["data"]
        w(f"  Total SVs: {d['total']}  PASS: {d['pass']}")
        w("")
    if head("delly", "Structural Variants (Delly)"):
        d = S["delly"]["data"]
        w(f"  Total SVs: {d['total']}  PASS: {d['pass']}")
        w("")
    if head("cnvpytor", "Copy Number Variants (CNVpytor)"):
        d = S["cnvpytor"]["data"]
        w(f"  Total CNVs: {d['total']} ({d['deletions']} deletions, {d['duplications']} duplications)")
        w(f"  Significant (e-val < 0.01): {d['significant']}")
        w("")
    if head("sv_consensus", "SV Consensus Merge"):
        w(f"  Consensus SVs (2+ callers): {S['sv_consensus']['data']['consensus']}")
        w("")

    if head("expansions", "Repeat Expansions (ExpansionHunter)"):
        d = S["expansions"]["data"]
        w(f"  Loci in output: {d['records']}")
        for l in d["key_loci"]:
            w(f"    {l['locus']:<8} {l['repeat_count']}")
        w("  (See docs/interpreting-results.md for disease thresholds)")
        w("")

    if head("telomere", "Telomere Length (TelomereHunter)"):
        w(f"  Telomere content: {S['telomere']['data']['tel_content']}")
        w("")

    if head("roh", "Runs of Homozygosity"):
        d = S["roh"]["data"]
        w(f"  Segments: {d['segments']}  Total: {d['total_mb']} MB  Largest: {d['largest_mb']} MB")
        w(f"  Autosomal ROH > 5MB: {len(d['autosomal_over_5mb'])}")
        for r in d["autosomal_over_5mb"]:
            w(f"    {r['region']}  {r['mb']}MB")
        w("")

    if head("haplogroup", "Mitochondrial Haplogroup"):
        w(f"  Haplogroup: {S['haplogroup']['data']['haplogroup']}")
        w("")

    if head("mito", "Mitochondrial Variants (Mutect2)"):
        d = S["mito"]["data"]
        w(f"  PASS variants: {d['pass']}")
        w(f"  Heteroplasmic (AF 0.05 to 0.95): {d['heteroplasmic']}")
        w("")

    if head("cpsr", "Cancer Predisposition Screening (CPSR)"):
        d = S["cpsr"]["data"]
        if d.get("html_report"):
            w("  HTML report in cpsr/ (open in a browser)")
        if d.get("classification"):
            w(f"  Classification breakdown ({d['classification_column']}):")
            for k, v in d["classification"].items():
                w(f"    {k}: {v}")
        elif "classification_column" in d:
            w("  Classification table has no CLASSIFICATION column.")
        w("")

    if head("clinical", "Clinical Variant Filter"):
        d = S["clinical"]["data"]
        w(f"  Clinical variants: {d['variants']} in {d['genes']} genes")
        w("  By impact: " + ", ".join(f"{k} {v}" for k, v in sorted(d["by_impact"].items())))
        w("")

    if head("slivar", "Variant Prioritization (slivar)"):
        d = S["slivar"]["data"]
        w(f"  Prioritized variants: {d['prioritized']}")
        if "compound_het_variants" in d:
            w(f"  Compound het candidates: {d['compound_het_variants']} variants across {d['compound_het_genes']} genes")
            w("  (Candidates only: unphased single-sample data. See docs/31-slivar.md)")
        w("")

    w("## Steps Not Run")
    w("---")
    if s["not_run"]:
        for t in s["not_run"]:
            w(f"  - {t}")
    else:
        w("  All major steps completed.")
    w("")
    w("## Not Assessed by This Pipeline")
    w("---")
    for t in s["not_assessed"]:
        w(f"  - {t}")
    w("")
    m = s["manifest"]
    if m:
        commit = next((r[2] for r in m if r[0] == "run" and r[1] == "git_commit"), "unknown")
        w(f"Run manifest: run_manifest.tsv (pipeline commit {commit}, "
          f"{sum(1 for r in m if r[0] == 'image')} images)")
    w("=" * 80)
    w("  DISCLAIMER: This is NOT a clinical report. Discuss findings with a")
    w("  qualified healthcare professional before making medical decisions.")
    w("=" * 80)
    return "\n".join(out) + "\n"


# --- html ------------------------------------------------------------------------

CSS = """
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         background: #f5f5f5; color: #333; line-height: 1.6; }
  .container { max-width: 1100px; margin: 0 auto; padding: 20px; }
  .header { background: linear-gradient(135deg, #1a5276 0%, #2e86c1 100%);
            color: white; padding: 30px; border-radius: 12px; margin-bottom: 24px; }
  .header h1 { font-size: 28px; margin-bottom: 8px; }
  .header .meta { opacity: 0.85; font-size: 14px; }
  .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 20px; margin-bottom: 24px; }
  .card { background: white; border-radius: 10px; padding: 24px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
  .card h2 { font-size: 18px; color: #1a5276; margin-bottom: 16px;
             padding-bottom: 8px; border-bottom: 2px solid #eee; }
  .stat { display: flex; justify-content: space-between; padding: 8px 0;
          border-bottom: 1px solid #f0f0f0; gap: 12px; }
  .stat:last-child { border-bottom: none; }
  .stat .label { color: #666; }
  .stat .value { font-weight: 600; text-align: right; }
  .badge { display: inline-block; padding: 2px 10px; border-radius: 12px;
           font-size: 13px; font-weight: 600; }
  .badge-green { background: #d5f5e3; color: #196f3d; }
  .badge-yellow { background: #fef9e7; color: #7d6608; }
  .badge-red { background: #fadbd8; color: #922b21; }
  .badge-gray { background: #eee; color: #666; }
  .stale { background: #fdebd0; color: #7e5109; border-radius: 6px; padding: 6px 10px;
           font-size: 13px; margin-bottom: 10px; }
  table { width: 100%; border-collapse: collapse; font-size: 14px; }
  th, td { padding: 8px 12px; text-align: left; border-bottom: 1px solid #eee; }
  th { background: #f8f9fa; font-weight: 600; color: #555; }
  .full-width { grid-column: 1 / -1; }
  ul.plain { margin-left: 18px; }
  .disclaimer { background: #fff3cd; border: 1px solid #ffc107; border-radius: 8px;
                padding: 16px; margin-top: 24px; font-size: 14px; }
  .footer { text-align: center; color: #999; font-size: 13px; margin-top: 24px; padding: 16px; }
  .footer table { font-size: 12px; color: #666; margin-top: 8px; }
  @media (max-width: 700px) { .grid { grid-template-columns: 1fr; } }
"""


def stat(label, value):
    return f'    <div class="stat"><span class="label">{E(str(label))}</span><span class="value">{value}</span></div>'


def badge(text, colour):
    return f'<span class="badge badge-{colour}">{E(str(text))}</span>'


def card(sec, title, body, full=False):
    """A card for one section; a missing one shows 'Not run', a stale one its note."""
    lines = [f'  <div class="card{" full-width" if full else ""}">', f"    <h2>{E(title)}</h2>"]
    if sec is not None and sec["state"] in ("stale", "unreadable"):
        lines.append(f'    <div class="stale">{E(stale_note(sec))}</div>')
    if sec is not None and sec["state"] in ("missing", "unreadable"):
        lines.append(stat("Status", badge("Not run" if sec["state"] == "missing" else "Unreadable", "gray")))
    else:
        lines += body
    lines.append("  </div>")
    return lines


def html_report(s):
    S = s["sections"]
    out = ["<!DOCTYPE html>", '<html lang="en">', "<head>", '<meta charset="UTF-8">',
           '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
           "<title>Personal Genome Pipeline Report</title>", "<style>" + CSS + "</style>", "</head>",
           "<body>", '<div class="container">', '<div class="header">', "  <h1>Genomic Analysis Report</h1>",
           '  <div class="meta">',
           f"    Sample: <strong>{E(s['sample'])}</strong> &nbsp;|&nbsp; Generated: {E(s['generated_utc'])}"
           + (f" &nbsp;|&nbsp; Latest run: {E(s['run']['started_utc'])}" if s["run"]["started_utc"] else ""),
           "  </div>", "</div>", '<div class="grid">']
    a = out.extend

    # QC
    cov, sexc = S["coverage"], S["sex_check"]
    md = cov["data"].get("mean_depth")
    inferred = sexc["data"].get("inferred_sex")
    declared = s["run"]["declared_sex"]
    sex_badge = E(inferred or "not available")
    if inferred and declared:
        sex_badge = badge(inferred, "green" if inferred == declared else "red")
    qc = [stat("Mean depth", f"{md:.1f}x" if isinstance(md, (int, float)) else "not available"),
          stat("Inferred sex (indexcov)", sex_badge),
          stat("Declared sex", E(declared or "not declared"))]
    qc_note = [f'    <div class="stale">{E(sec["title"] + ": " + stale_note(sec))}</div>'
               for sec in (cov, sexc) if stale_note(sec)]
    a(['  <div class="card">', "    <h2>Quality Control</h2>"] + qc_note + qc + ["  </div>"])

    v = S["variants"]["data"]
    a(card(S["variants"], "Variant Calling", [stat("Total variants", v.get("total")), stat("PASS variants", v.get("pass")),
                                              stat("SNPs", v.get("snps")), stat("Indels", v.get("indels"))]))

    c = S["clinvar"]["data"]
    n = c.get("count", 0)
    by = c.get("by_stars", {})
    cv_body = [stat("ClinVar matches", badge(n, "yellow" if n > 5 else "green")),
               stat("By stars (4/3/2/1/0)", E(" / ".join(str(by.get(k, 0)) for k in ("4", "3", "2", "1", "0")))),
               stat("ClinVar file date", E(s["databases"]["clinvar_release"] or "unknown"))]
    a(card(S["clinvar"], "ClinVar Screening", cv_body))

    pg = S["pharmcat"]["data"]
    cp = S["cpic"]["data"]
    pgx_body = [stat("PharmCAT report", badge("Complete", "green") if S["pharmcat"]["state"] in ("ok", "stale")
                     else badge("Not run", "gray"))]
    if pg.get("version"):
        pgx_body.append(stat("PharmCAT version", E(pg["version"])))
    if S["cpic"]["state"] in ("ok", "stale"):
        pgx_body += [stat("Genes with a non-normal phenotype", cp["non_normal"]),
                     stat("Genes with more than one possible result", cp.get("ambiguous", 0)),
                     stat("Genes not called", cp["not_called"])]
        if S["cpic"]["state"] == "stale":
            pgx_body.insert(0, f'    <div class="stale">{E(stale_note(S["cpic"]))}</div>')
    if S["pypgx"]["state"] in ("ok", "stale") and "comparison" in S["pypgx"]["data"]:
        pgx_body.append(stat("PharmCAT vs pypgx conflicts", S["pypgx"]["data"]["comparison"]["conflicts"]))
    a(card(S["pharmcat"], "Pharmacogenomics", pgx_body))

    calls = s["cyp2d6"]["calls"]
    agree = s["cyp2d6"]["agree"]
    cy_body = [stat(k, E(v or "not run")) for k, v in calls.items()]
    cy_body.append(stat("Agreement", badge("yes", "green") if agree else badge("no", "red") if agree is False
                        else badge("fewer than two calls", "gray")))
    a(["  <div class=\"card\">", "    <h2>CYP2D6 Across Callers</h2>"]
      + [f'    <div class="stale">{E(S[k]["title"] + ": " + stale_note(S[k]))}</div>'
         for k in ("cpic", "pypgx", "cyrius") if S[k]["state"] == "stale"]
      + cy_body + ["  </div>"])

    hla = S["hla"]["data"]
    a(card(S["hla"], "HLA Typing", [stat(l["gene"], E(" / ".join(l["alleles"]) or "no call")) for l in hla.get("loci", [])]
           + [stat("HLA database", E(s["databases"]["hla_database"] or "unknown"))]))

    m, dl, cn, sv = (S[k]["data"] for k in ("manta", "delly", "cnvpytor", "sv_consensus"))
    sv_body = [stat("Manta SVs (total / PASS)", f"{m.get('total', 'N/A')} / {m.get('pass', 'N/A')}"),
               stat("Delly SVs (total / PASS)", f"{dl.get('total', 'N/A')} / {dl.get('pass', 'N/A')}"),
               stat("CNVpytor CNVs", cn.get("total", "N/A")),
               stat("Consensus SVs (2+ callers)", sv.get("consensus", "N/A"))]
    notes = [f'    <div class="stale">{E(S[k]["title"] + ": " + stale_note(S[k]))}</div>'
             for k in ("manta", "delly", "cnvpytor", "sv_consensus") if stale_note(S[k])]
    a(['  <div class="card">', "    <h2>Structural Variants</h2>"] + notes + sv_body + ["  </div>"])

    cs = S["cpsr"]["data"]
    cpsr_body = [stat("CPSR report", badge("Complete", "green") if cs.get("html_report") else badge("No HTML", "gray"))]
    for k, val in (cs.get("classification") or {}).items():
        cpsr_body.append(stat(k, val))
    a(card(S["cpsr"], "Cancer Predisposition", cpsr_body))

    eh = S["expansions"]["data"]
    eh_body = [stat("ExpansionHunter", badge("Complete", "green"))]
    if eh.get("key_loci"):
        eh_body.append("    <table><tr><th>Locus</th><th>Repeat Count</th></tr>"
                       + "".join(f"<tr><td>{E(l['locus'])}</td><td>{E(str(l['repeat_count']))}</td></tr>"
                                 for l in eh["key_loci"]) + "</table>")
    a(card(S["expansions"], "Repeat Expansions", eh_body))

    r, hg, tl = S["roh"]["data"], S["haplogroup"]["data"], S["telomere"]["data"]
    anc = [stat("Mitochondrial haplogroup", E(hg.get("haplogroup", "N/A"))),
           stat("ROH total", f"{r.get('total_mb', 'N/A')} MB"),
           stat("ROH largest segment", f"{r.get('largest_mb', 'N/A')} MB"),
           stat("Autosomal ROH > 5 MB", len(r["autosomal_over_5mb"]) if "autosomal_over_5mb" in r else "N/A"),
           stat("Telomere content", E(str(tl.get("tel_content", "N/A"))))]
    notes = [f'    <div class="stale">{E(S[k]["title"] + ": " + stale_note(S[k]))}</div>'
             for k in ("haplogroup", "roh", "telomere") if stale_note(S[k])]
    a(['  <div class="card">', "    <h2>Ancestry &amp; Identity</h2>"] + notes + anc + ["  </div>"])

    mi = S["mito"]["data"]
    a(card(S["mito"], "Mitochondrial Analysis", [stat("chrM variants (PASS)", mi.get("pass")),
                                                 stat("Heteroplasmic (AF 0.05 to 0.95)", mi.get("heteroplasmic"))]))

    cl = S["clinical"]["data"]
    a(card(S["clinical"], "Clinical Variant Filter",
           [stat("Total interesting variants", cl.get("variants")), stat("Genes", cl.get("genes"))]
           + [stat(f"{k} impact", v) for k, v in sorted((cl.get("by_impact") or {}).items())]))

    sl = S["slivar"]["data"]
    sl_body = [stat("Prioritized variants", sl.get("prioritized"))]
    if "compound_het_variants" in sl:
        sl_body.append(stat("Compound het candidates", f"{sl['compound_het_variants']} variants, {sl['compound_het_genes']} genes"))
    a(card(S["slivar"], "Variant Prioritization (slivar)", sl_body))

    if S["prs"]["state"] in ("ok", "stale"):
        rows = "".join(f"<tr><td>{E(x['condition'])}</td><td>{E(x['pgs_id'])}</td><td>{E(x['score'])}</td>"
                       f"<td>{E(x['matched'])}/{E(x['total'])}</td></tr>\n" for x in S["prs"]["data"]["scores"])
        a(card(S["prs"], "Polygenic Risk Scores (raw, not percentiles)",
               ["    <table>", "      <tr><th>Condition</th><th>PGS</th><th>Score</th><th>Variants matched</th></tr>",
                rows.rstrip("\n"), "    </table>"], full=True))

    if S["cpic"]["state"] in ("ok", "stale") and (cp["non_normal"] or cp.get("ambiguous")):
        rows = "".join(f"<tr><td>{E(g['gene'])}</td><td>{E(g['diplotype'])}</td><td>{E(g['phenotype'])}</td></tr>\n"
                       for g in cp["genes"] if g["status"] in ("non-normal", "ambiguous"))
        body = ["    <table>", "      <tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th></tr>", rows.rstrip("\n"),
                "    </table>"]
        body += [f'    <div class="stale">WARNING: {E(x)}</div>' for x in cp["warnings"]]
        body.append(f"    <p>Drugs for each gene: cpic/{E(s['sample'])}_cpic_recommendations.txt</p>")
        a(card(S["cpic"], "Genes With a Non-Normal or Unresolved Phenotype (CPIC)", body, full=True))

    if S["clinvar"]["state"] in ("ok", "stale") and c.get("hits"):
        rows = [f"<tr><td>{E(h['chrom'])}</td><td>{h['pos']}</td><td>{E(h['ref'])}</td><td>{E(h['alt'])}</td>"
                f"<td>{E(h['genotype'])}</td><td>{E(h['gene'])}</td><td>{E(h['significance'])}</td>"
                f"<td>{E(h['review_status'])}</td><td>{h['stars']}</td></tr>" for h in c["hits"][:50]]
        a(card(S["clinvar"], f"ClinVar Hits (best-reviewed first, {min(n, 50)} of {n})",
               ["    <table>",
                "      <tr><th>Chr</th><th>Position</th><th>Ref</th><th>Alt</th><th>Genotype</th><th>Gene</th>"
                "<th>Significance</th><th>Review status</th><th>Stars</th></tr>"]
               + ["      " + x for x in rows] + ["    </table>"], full=True))

    lists = ['  <div class="card full-width">', "    <h2>Steps Not Run</h2>"]
    lists += (['    <ul class="plain">'] + [f"      <li>{E(t)}</li>" for t in s["not_run"]] + ["    </ul>"]
              if s["not_run"] else ["    <p>All major steps completed.</p>"])
    lists += ["    <h2 style=\"margin-top:16px\">Not Assessed by This Pipeline</h2>", '    <ul class="plain">']
    lists += [f"      <li>{E(t)}</li>" for t in s["not_assessed"]] + ["    </ul>", "  </div>"]
    a(lists)
    a(["</div>", "", '<div class="disclaimer">',
       "  <strong>Disclaimer:</strong> This report is for educational and research purposes only.",
       "  It is not a clinical diagnosis. Always discuss genomic findings with a qualified healthcare",
       "  professional before making any medical decisions. Variants of Uncertain Significance (VUS)",
       "  are not clinically actionable.", "</div>", "", '<div class="footer">',
       '  Generated by <a href="https://github.com/GeiserX/Personal-Genome-Pipeline">Personal-Genome-Pipeline</a>',
       "  &mdash; 100% local analysis, no data uploaded"])
    m = s["manifest"]
    if m:
        a(["  <details><summary>Run manifest (run_manifest.tsv)</summary>", "  <table>"]
          + [f"    <tr><td>{E(r[0])}</td><td>{E(r[1])}</td><td>{E(r[2])}</td></tr>" for r in m]
          + ["  </table></details>"])
    else:
        a(["  <p>No run_manifest.tsv in the sample folder.</p>"])
    a(["</div>", "", "</div>", "</body>", "</html>"])
    return "\n".join(out) + "\n"


def write(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write(text)
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--summary", help="read this summary JSON instead of collecting one")
    ap.add_argument("--sample")
    ap.add_argument("--sample-dir")
    ap.add_argument("--json", help="where to write the collected summary (default SAMPLE_DIR/summary.json)")
    ap.add_argument("-o", "--out", action="append", default=[],
                    help="report to write: .txt for text, .html for HTML (repeatable)")
    a = ap.parse_args(argv)
    for o in a.out:
        if not o.endswith((".txt", ".html")):
            ap.error(f"-o {o}: the name must end in .txt or .html")
    if a.summary:
        with open(a.summary) as f:
            s = json.load(f)
    elif a.sample and a.sample_dir:
        s = collect_summary.collect(a.sample, a.sample_dir)
        write(a.json or os.path.join(a.sample_dir, "summary.json"), json.dumps(s, indent=1) + "\n")
    else:
        ap.error("give --summary, or --sample and --sample-dir")
    for o in a.out:
        write(o, text_report(s) if o.endswith(".txt") else html_report(s))
    states = {}
    for sec in s["sections"].values():
        states[sec["state"]] = states.get(sec["state"], 0) + 1
    print("Sections: " + ", ".join(f"{k} {v}" for k, v in sorted(states.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
