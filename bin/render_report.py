#!/usr/bin/env python3
"""render_report.py: the text and HTML reports, both from one summary.

  render_report.py --sample S --sample-dir DIR [--json OUT.json] [-o REPORT.txt] [-o REPORT.html]
  render_report.py --summary IN.json [-o REPORT.txt] [-o REPORT.html]

-o writes the text report for a .txt name and the HTML report for .html.

With --sample-dir it first collects the summary (bin/collect_summary.py) and
writes it to --json (default DIR/summary.json). Every number in both reports
comes from that summary, so the two cannot disagree. A section whose result is
from an earlier run (state "stale") is shown with that date and a STALE mark,
and one whose step failed in the latest run (state "failed") says Failed.
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
    if sec["state"] == "failed":
        return f"FAILED: {sec['note']}"
    return ""


# Sentences both reports print, kept to what the pipeline can back.
CLINVAR_FREQ_NOTE = ("Population frequency is not checked here. Look up a homozygous P/LP hit in gnomAD "
                     "before reading much into it.")
HLA_NOTE = ("HLA typing from short-read WGS is approximate. If either HLA-A or HLA-B allele has a T1K quality "
            "of 0 or less, that gene is not passed to PharmCAT (step 36).")
EH_NOTE = ("Repeat-expansion calls from short reads can be wrong at some loci, and a long expansion "
           "cannot be sized. Disease thresholds: docs/interpreting-results.md.")
TELOMERE_UNIT = "intratelomeric reads per million reads with 48-52% GC"
TELOMERE_NOTE = (f"Relative telomere content ({TELOMERE_UNIT}), not a telomere length: compare it only with "
                 "samples of a similar age sequenced on the same platform (docs/10-telomere-analysis.md).")
PRS_ANCESTRY_NOTE = ("Most PRS scores were developed in European-ancestry populations and predict less well for "
                     "other ancestries, even with an ancestry-adjusted percentile.")


def hla_text(locus):
    """One HLA locus: each allele with its T1K quality, and the low-confidence mark."""
    q = locus.get("quality") or []
    parts = [a + (f" (quality {q[i]})" if i < len(q) and q[i] != "" else "") for i, a in enumerate(locus["alleles"])]
    txt = " / ".join(parts) or "no call"
    if locus.get("low_confidence"):
        txt += "; low confidence (T1K quality 0 or below)"
        if locus.get("withheld_from_pharmcat"):
            txt += ", not passed to PharmCAT"
    return txt


def cyp2d6_verdict(c):
    """Step 36's CYP2D6 verdict as one line, or '' without its table."""
    v = c.get("consensus")
    if not v:
        return ""
    where = "passed to PharmCAT" if v.get("passed_to_pharmcat") else "not passed to PharmCAT"
    return f"{v.get('result') or '?'}, {where} ({v.get('reason') or 'no reason given'})"


# The PharmCAT column is not a third caller: PharmCAT calls no CYP2D6 from a
# VCF, it shows the call step 36 passed to it.
CYP2D6_LABELS = {"PharmCAT": "PharmCAT (step 36's call)"}


# A section the summary lacks (one written before the section existed).
MISSING = {"title": "", "state": "missing", "note": None, "data": {}}


def freemix_text(q):
    """FREEMIX with its verdict, as both reports print it."""
    if q.get("contamination") == "not_run" or not q.get("freemix"):
        return "not run"
    warn = q.get("freemix_warn_above", "")
    few = ("; fewer than 1,000 panel markers had reads, so VerifyBamID2 ran without its marker check"
           if q.get("verifybamid2_marker_check") == "skipped" else "")
    if q.get("contamination") == "warn":
        return f"{q['freemix']} (above {warn}: possible contamination, see docs/33-sample-qc.md{few})"
    return f"{q['freemix']} (warning above {warn}{few})"


def prs_chrx_note(d):
    """A sentence when chrX rows of a score were left out (no sex given), else ''."""
    left = [f"{r['pgs_id']} ({r['chrx'].split()[0]})" for r in d.get("scores", [])
            if (r.get("chrx") or "").endswith("left out")]
    if not left:
        return ""
    return (f" Rows on chrX were left out, because the sample's sex was not given: {', '.join(left)}. "
            "Give the sex to score them (docs/25-prs.md).")


def prs_note(d):
    """The line under the PRS scores, in both reports: what the numbers can be compared with."""
    x = prs_chrx_note(d)
    if d.get("adjusted"):
        anc = d.get("ancestry") or {}
        group = anc.get("population") or next((r.get("group") for r in d["scores"] if r.get("group")), "")
        panel = f" of the {anc['panel']} reference panel" if anc.get("panel") else " of the reference panel"
        low = " The ancestry match is low-confidence: read the percentile with care." if anc.get("low_confidence") == "True" else ""
        return (f"Percentile: where the score falls among the {group or 'most similar'} samples{panel}, "
                f"the group whose genetic ancestry is most similar to this sample's (pgsc_calc).{low} "
                f"A percentile is not a risk. {PRS_ANCESTRY_NOTE} See docs/25-prs.md." + x)
    return ("Raw score only: no ancestry reference panel was installed, so there is no percentile and the sum "
            f"cannot be compared with anyone (scripts/setup.sh --ancestry-panel adds it). {PRS_ANCESTRY_NOTE} "
            "See docs/25-prs.md." + x)


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
        sec = S.get(key, MISSING)
        if sec["state"] == "missing":
            return False
        w(f"## {title or sec['title']}")
        w("---")
        note = stale_note(sec)
        if note:
            w(f"  [{note}]")
        if sec["state"] == "failed":
            w("")
        return sec["state"] in ("ok", "stale")

    # QC
    cov, sexc, sqc = S["coverage"], S["sex_check"], S.get("sample_qc", MISSING)
    if cov["state"] != "missing" or sexc["state"] != "missing" or sqc["state"] != "missing":
        w("## Quality control")
        w("---")
        for sec in (cov, sexc, sqc):
            if stale_note(sec):
                w(f"  [{sec['title']}: {stale_note(sec)}]")
        md = cov["data"].get("mean_depth")
        w(f"  Mean depth: {md:.1f}x" if isinstance(md, (int, float)) else "  Mean depth: not available")
        inferred = sexc["data"].get("inferred_sex")
        declared = s["run"]["declared_sex"]
        w(f"  Inferred sex: {inferred or 'not available'}")
        w(f"  Declared sex: {declared or 'not declared'}"
          + ("  ** DOES NOT MATCH the inferred sex **" if inferred and declared and inferred != declared else ""))
        if sqc["state"] in ("ok", "stale"):
            q = sqc["data"]
            w(f"  Sex from the reads (somalier): {q.get('inferred_sex')}"
              + ("  ** DOES NOT MATCH the declared sex **" if q.get("sex_check") == "mismatch" else ""))
            if q.get("sex_check") == "not_checked":
                w(f"    not checked: {q.get('sex_check_reason')}")
            w(f"  Contamination (VerifyBamID2 FREEMIX): {freemix_text(q)}")
            if q.get("same_person_as"):
                w(f"  ** The same person as: {q['same_person_as']} (a duplicate or a sample swap) **")
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
        w(f"  {CLINVAR_FREQ_NOTE}")
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
          f"{d.get('ambiguous', 0)}  Not called: {d['not_called']}  "
          f"Called without a function phenotype (no drug guidance): {d.get('unclassified', 0)}")
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
    verdict = cyp2d6_verdict(s["cyp2d6"])
    if any(calls.values()) or verdict:
        w("## CYP2D6 across callers")
        w("---")
        for k, v in calls.items():
            w(f"  {CYP2D6_LABELS.get(k, k):<26} {v or 'not run'}")
        agree = s["cyp2d6"]["agree"]
        w("  pypgx and Cyrius agree: " + ("yes" if agree else "NO, review before acting on CYP2D6" if agree is False
                                          else "fewer than two usable calls"))
        if verdict:
            w(f"  Step 36: {verdict}")
        w("")

    if head("hla", "HLA Typing (T1K)"):
        for l in S["hla"]["data"]["loci"]:
            w(f"    {l['gene']:<8} {hla_text(l)}")
        w(f"  HLA database: {s['databases']['hla_database'] or 'unknown'}")
        w(f"  {HLA_NOTE}")
        w("")

    if head("prs", "Polygenic Risk Scores"):
        d = S["prs"]["data"]
        for r in d["scores"]:
            pct = (f"  percentile {r['percentile']} ({r.get('group') or 'group unknown'})"
                   if r.get("percentile") else "")
            w(f"  {r['condition']:<35} {r['score']} ({r['matched']}/{r['total']} variants matched) {r['pgs_id']}{pct}")
        w("  " + prs_note(d))
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
        if not d.get("stranger"):
            w("  Stranger (step 9b) did not run: no locus is marked as normal or expanded.")
        elif d.get("flagged"):
            w(f"  Outside the normal range (Stranger): {len(d['flagged'])}")
            for l in d["flagged"]:
                w(f"    {l['locus']:<8} {l['repeat_count']}  {l['status']}")
        else:
            w("  Outside the normal range (Stranger): none")
        if d.get("stranger") and d.get("no_status"):
            w(f"  Loci without a Stranger status (not in its catalog): {d['no_status']}")
        w(f"  {EH_NOTE}")
        w("")

    if head("telomere", "Telomere content (relative, TelomereHunter)"):
        w(f"  Telomere content: {S['telomere']['data']['tel_content']} {TELOMERE_UNIT}")
        w(f"  {TELOMERE_NOTE}")
        w("")

    if head("roh", "Runs of Homozygosity"):
        d = S["roh"]["data"]
        w(f"  Segments: {d['segments']}  Total: {d['total_mb']} MB  Largest: {d['largest_mb']} MB")
        w(f"  Autosomal ROH > 5MB: {len(d['autosomal_over_5mb'])}")
        for r in d["autosomal_over_5mb"]:
            w(f"    {r['region']}  {r['mb']}MB")
        w("")

    if head("haplogroup", "Mitochondrial Haplogroup"):
        d = S["haplogroup"]["data"]
        w(f"  Haplogroup: {d['haplogroup']}")
        w(f"  Contamination (haplocheck): {contamination_text(d)}")
        w("")

    if head("y_haplogroup", "Y-Chromosome Haplogroup (Yleaf)"):
        d = S["y_haplogroup"]["data"]
        w(f"  Haplogroup: {d['haplogroup']} ({d['valid_markers']} markers, QC-score {d['qc_score']})")
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
        w(f"  Clinical variants: {d['variants']}" + (f" in {d['genes']} genes" if "genes" in d else ""))
        if d.get("by_impact"):
            w("  By impact: " + ", ".join(f"{k} {v}" for k, v in sorted(d["by_impact"].items())))
        sf = d.get("acmg_sf")
        if sf:
            w(f"  Secondary-findings genes ({sf['version']}, {sf['genes_on_list']} genes): "
              f"{len(sf['clinvar_hits'])} ClinVar P/LP, {len(sf['high_impact'])} rare HIGH impact")
            for x in sf["clinvar_hits"]:
                w(f"    {x['gene']:<8} {x['variant']}  {x['genotype']}  ClinVar {x['significance']}")
            for x in sf["high_impact"]:
                w(f"    {x['gene']:<8} {x['variant']}  {x['genotype']}  {x['consequence']}")
            w(f"    {ACMG_NOTE}")
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


ACMG_NOTE = ("A list to review with a clinician: the ACMG reporting rules per gene (for example HFE "
             "homozygous C282Y only, BTD and CYP27A1 two variants) are not applied here.")


def contamination_text(d):
    if "contamination_status" not in d:
        return "not checked (needs step 20's Mutect2 chrM calls)"
    s = d["contamination_status"]
    word = {"YES": "YES, two mtDNA haplogroups in the reads (another person's DNA?)", "NO": "no",
            "ND": "not determined"}.get(s, s)
    return f"{word} (level {d.get('contamination_level', '.')})"


def stat(label, value):
    return f'    <div class="stat"><span class="label">{E(str(label))}</span><span class="value">{value}</span></div>'


def small_note(text):
    """A sentence under a card's numbers."""
    return f'    <p style="font-size:13px;color:#555;margin-top:8px">{E(text)}</p>'


def badge(text, colour):
    return f'<span class="badge badge-{colour}">{E(str(text))}</span>'


def done_badge(sec):
    """The Status badge of a section that ran: Complete, or Stale for a result from an earlier run."""
    return badge("Stale", "yellow") if sec["state"] == "stale" else badge("Complete", "green")


def card(sec, title, body, full=False):
    """A card for one section; a missing one shows 'Not run', a failed one
    'Failed', a stale one its note."""
    lines = [f'  <div class="card{" full-width" if full else ""}">', f"    <h2>{E(title)}</h2>"]
    if sec is not None and sec["state"] in ("stale", "unreadable", "failed"):
        lines.append(f'    <div class="stale">{E(stale_note(sec))}</div>')
    if sec is not None and sec["state"] in ("missing", "unreadable", "failed"):
        lines.append(stat("Status", {"missing": badge("Not run", "gray"), "failed": badge("Failed", "red")}
                          .get(sec["state"], badge("Unreadable", "gray"))))
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
    sqc = S.get("sample_qc", MISSING)
    if sqc["state"] in ("ok", "stale"):
        q = sqc["data"]
        check = q.get("sex_check")
        qc.append(stat("Sex from the reads (somalier)",
                       badge(q.get("inferred_sex", "unknown"),
                             "green" if check == "ok" else "red" if check == "mismatch" else "gray")))
        if check == "not_checked":
            qc.append(stat("Not checked", f'<span style="font-weight:normal;font-size:13px">'
                                          f'{E(q.get("sex_check_reason", ""))}</span>'))
        colour = {"ok": "green", "warn": "yellow"}.get(q.get("contamination"), "gray")
        qc.append(stat("Contamination (FREEMIX)", badge(freemix_text(q), colour)))
        if q.get("same_person_as"):
            qc.append(stat("The same person as", badge(q["same_person_as"], "red")))
    else:
        qc.append(stat("Identity and contamination (step 33)", badge("Not run", "gray")))
    qc_note = [f'    <div class="stale">{E(sec["title"] + ": " + stale_note(sec))}</div>'
               for sec in (cov, sexc, sqc) if stale_note(sec)]
    a(['  <div class="card">', "    <h2>Quality Control</h2>"] + qc_note + qc + ["  </div>"])

    v = S["variants"]["data"]
    a(card(S["variants"], "Variant Calling", [stat("Total variants", v.get("total")), stat("PASS variants", v.get("pass")),
                                              stat("SNPs", v.get("snps")), stat("Indels", v.get("indels"))]))

    c = S["clinvar"]["data"]
    n = c.get("count", 0)
    by = c.get("by_stars", {})
    # Never green: no hit is not an all-clear, and a hit is something to look
    # at whatever its count, its genotype or its stars.
    cv_body = [stat("ClinVar matches", badge(n, "yellow" if n else "gray")),
               stat("By stars (4/3/2/1/0)", E(" / ".join(str(by.get(k, 0)) for k in ("4", "3", "2", "1", "0")))),
               stat("ClinVar file date", E(s["databases"]["clinvar_release"] or "unknown")),
               small_note(CLINVAR_FREQ_NOTE)]
    a(card(S["clinvar"], "ClinVar Screening", cv_body))

    pg = S["pharmcat"]["data"]
    cp = S["cpic"]["data"]
    pgx_body = [stat("PharmCAT report", badge("Complete", "green") if S["pharmcat"]["state"] in ("ok", "stale")
                     else badge("Not run", "gray"))]
    if pg.get("version"):
        pgx_body.append(stat("PharmCAT version", E(pg["version"])))
    if S["pypgx"]["state"] in ("ok", "stale") and "comparison" in S["pypgx"]["data"]:
        pgx_body.append(stat("PharmCAT vs pypgx conflicts", S["pypgx"]["data"]["comparison"]["conflicts"]))
    a(card(S["pharmcat"], "Pharmacogenomics", pgx_body))

    cpic_body = [stat("Status", done_badge(S["cpic"]))]
    if not cp.get("parse_failed"):
        cpic_body += [stat("Genes with a non-normal phenotype", cp.get("non_normal", 0)),
                      stat("Genes with more than one possible result", cp.get("ambiguous", 0)),
                      stat("Genes not called", cp.get("not_called", 0)),
                      stat("Genes called without a function phenotype (no drug guidance)", cp.get("unclassified", 0))]
    else:
        cpic_body.append(stat("PharmCAT report", badge("could not be parsed", "red")))
    cpic_body.append(stat("Tip", '<span style="font-weight:normal;font-size:13px">Recommendations per gene: '
                                 f'cpic/{E(s["sample"])}_cpic_recommendations.txt</span>'))
    a(card(S["cpic"], "CPIC Drug Recommendations", cpic_body))

    calls = s["cyp2d6"]["calls"]
    agree = s["cyp2d6"]["agree"]
    cy_body = [stat(CYP2D6_LABELS.get(k, k), E(v or "not run")) for k, v in calls.items()]
    cy_body.append(stat("pypgx and Cyrius agree", badge("yes", "green") if agree else badge("no", "red")
                        if agree is False else badge("fewer than two usable calls", "gray")))
    cons = s["cyp2d6"].get("consensus")
    if cons:
        cy_body.append(stat("Step 36", badge(cons.get("result") or "?", "green" if cons.get("passed_to_pharmcat")
                                             else "yellow")))
        cy_body.append(small_note(("Passed to PharmCAT: " if cons.get("passed_to_pharmcat") else "Not passed to PharmCAT: ")
                                  + (cons.get("reason") or "no reason given")))
    a(["  <div class=\"card\">", "    <h2>CYP2D6 Across Callers</h2>"]
      + [f'    <div class="stale">{E(S[k]["title"] + ": " + stale_note(S[k]))}</div>'
         for k in ("cpic", "pypgx", "cyrius") if S[k]["state"] == "stale"]
      + cy_body + ["  </div>"])

    hla = S["hla"]["data"]
    a(card(S["hla"], "HLA Typing", [stat(l["gene"], E(hla_text(l))) for l in hla.get("loci", [])]
           + [stat("HLA database", E(s["databases"]["hla_database"] or "unknown")), small_note(HLA_NOTE)]))

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
    eh_body = [stat("Status", done_badge(S["expansions"]))]
    if not eh.get("stranger"):
        eh_body.append(stat("Outside the normal range (Stranger)", badge("not run", "gray")))
    else:
        fl = eh.get("flagged") or []
        eh_body.append(stat("Outside the normal range (Stranger)", badge(len(fl), "yellow" if fl else "gray")))
    rows = [(l["locus"], l["repeat_count"], l["status"]) for l in eh.get("flagged") or []]
    rows += [(l["locus"], l["repeat_count"], "") for l in eh.get("key_loci") or []
             if l["locus"] not in {r[0] for r in rows}]
    if rows:
        # The Stranger column only when it flagged a locus: an empty cell is a locus it did not flag.
        st_col = bool(eh.get("flagged"))
        eh_body.append("    <table><tr><th>Locus</th><th>Repeat Count</th>" + ("<th>Stranger</th>" if st_col else "")
                       + "</tr>" + "".join(f"<tr><td>{E(x)}</td><td>{E(str(y))}</td>"
                                         + (f"<td>{E(z)}</td>" if st_col else "") + "</tr>" for x, y, z in rows)
                       + "</table>")
    eh_body.append(small_note(EH_NOTE))
    a(card(S["expansions"], "Repeat Expansions", eh_body))

    r, hg, tl = S["roh"]["data"], S["haplogroup"]["data"], S["telomere"]["data"]
    a(card(S["roh"], "Runs of Homozygosity",
           [stat("Status", done_badge(S["roh"])),
            stat("ROH total", f"{r.get('total_mb', 0):.1f} MB"),
            stat("ROH largest segment", f"{r.get('largest_mb', 0):.1f} MB"),
            stat("Segments", r.get("segments", 0)),
            stat("Autosomal ROH > 5 MB", len(r.get("autosomal_over_5mb", [])))]))
    a(card(S["haplogroup"], "Mitochondrial Haplogroup",
           [stat("Status", done_badge(S["haplogroup"])), stat("Haplogroup", E(hg.get("haplogroup", "."))),
            stat("Contamination (haplocheck)", E(contamination_text(hg)))]))
    if S.get("y_haplogroup", MISSING)["state"] != "missing":
        yh = S["y_haplogroup"]["data"]
        a(card(S["y_haplogroup"], "Y-Chromosome Haplogroup",
               [stat("Haplogroup", E(yh.get("haplogroup", "."))), stat("Markers", E(str(yh.get("valid_markers", ".")))),
                stat("QC-score", E(str(yh.get("qc_score", "."))))]))
    a(card(S["telomere"], "Telomere content (relative)",
           [stat("Telomere content", E(str(tl.get("tel_content", ".")))), small_note(TELOMERE_NOTE)]))

    mi = S["mito"]["data"]
    a(card(S["mito"], "Mitochondrial Analysis", [stat("chrM variants (PASS)", mi.get("pass")),
                                                 stat("Heteroplasmic (AF 0.05 to 0.95)", mi.get("heteroplasmic"))]))

    cl = S["clinical"]["data"]
    a(card(S["clinical"], "Clinical Variant Filter",
           [stat("Total interesting variants", cl.get("variants"))]
           + ([stat("Genes", cl["genes"])] if "genes" in cl else [])
           + [stat(f"{k} impact", v) for k, v in sorted((cl.get("by_impact") or {}).items())]))
    sf = cl.get("acmg_sf")
    if sf and (sf["clinvar_hits"] or sf["high_impact"]):
        rows = "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(x['variant'])}</td><td>{E(x['genotype'])}</td>"
                       f"<td>{E('ClinVar ' + x['significance'])}</td></tr>\n" for x in sf["clinvar_hits"])
        rows += "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(x['variant'])}</td><td>{E(x['genotype'])}</td>"
                        f"<td>{E(x['consequence'])}</td></tr>\n" for x in sf["high_impact"])
        a(card(S["clinical"], f"Secondary-Findings Genes ({sf['version']})",
               ["    <table>", "      <tr><th>Gene</th><th>Variant</th><th>Genotype</th><th>Evidence</th></tr>",
                rows.rstrip("\n"), "    </table>", f"    <p>{E(ACMG_NOTE)}</p>"], full=True))

    sl = S["slivar"]["data"]
    sl_body = [stat("Prioritized variants", sl.get("prioritized"))]
    if "compound_het_variants" in sl:
        sl_body.append(stat("Compound het candidates", f"{sl['compound_het_variants']} variants, {sl['compound_het_genes']} genes"))
    a(card(S["slivar"], "Variant Prioritization (slivar)", sl_body))

    if S["prs"]["state"] in ("ok", "stale"):
        d = S["prs"]["data"]
        rows = "".join(f"<tr><td>{E(x['condition'])}</td><td>{E(x['pgs_id'])}</td><td>{E(x['score'])}</td>"
                       f"<td>{E(x['matched'])}/{E(x['total'])}</td>"
                       f"<td>{E(x['percentile'] + ' (' + (x.get('group') or 'group unknown') + ')') if x.get('percentile') else 'raw score only'}</td></tr>\n"
                       for x in d["scores"])
        a(card(S["prs"], "Polygenic Risk Scores" + ("" if d.get("adjusted") else " (raw, not percentiles)"),
               ["    <table>", "      <tr><th>Condition</th><th>PGS</th><th>Score</th><th>Variants matched</th><th>Percentile</th></tr>",
                rows.rstrip("\n"), "    </table>", f"    <p>{E(prs_note(d))}</p>"], full=True))

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
    ap.add_argument("--declared-sex", help="the samplesheet's sex, when no run status or manifest records it")
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
        s = collect_summary.collect(a.sample, a.sample_dir, a.declared_sex)
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
