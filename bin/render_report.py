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
import re
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
CLINVAR_FREQ_NOTE = ("Check every P/LP hit's population frequency, a heterozygous hit in a dominant gene as well as "
                     "a homozygous one: a hit above 1% in gnomAD is unlikely to be a high-penetrance finding, and a "
                     "hit whose frequency is unknown has not been checked.")
HLA_NOTE = ("HLA typing from short-read WGS is approximate. If either HLA-A or HLA-B allele has a T1K quality "
            "of 0 or less, that gene is not passed to PharmCAT (step 36).")
EH_NOTE = ("Repeat-expansion calls from short reads can be wrong at some loci, and a long expansion "
           "cannot be sized. Disease thresholds: docs/interpreting-results.md.")
TELOMERE_UNIT = "intratelomeric reads per million reads with 48-52% GC"
TELOMERE_NOTE = (f"Relative telomere content ({TELOMERE_UNIT}), not a telomere length: compare it only with "
                 "samples of a similar age sequenced on the same platform (docs/10-telomere-analysis.md).")
ANCESTRY_FINE_NOTE = ("Each name is a group of reference populations we put together and named, not a nationality. "
                      "Continental groups are well apart; neighbouring regions are not: splits within Iberia are weak "
                      "for every method, and Greek and Southern Italian ancestry overlap. Read a regional fraction "
                      "with its bootstrap range. See docs/26b-ancestry-fine.md.")
ANCESTRY_FINE_LEVELS = (("continent", "Continental groups"), ("region", "Regional groups"))
TRAITS_NOTE = ("Each reading is what the cited studies found on average in the people they studied, not a "
               "diagnosis. A row marked ancestry reads differently in other populations; its note says how.")
APOE_WARNING = ("APOE and Alzheimer's disease: a risk factor, not a diagnosis. Many people choose not to learn "
                "their APOE result; docs/42-traits.md explains it.")
EXOMISER_NOTE = ("A ranking to bring to a clinical geneticist, not a diagnosis: a high score means a rare, likely "
                 "damaging variant in a gene linked to symptoms like these, not that it causes them, and an empty or "
                 "low list does not rule a genetic cause out. Coding and splice variants only (docs/48-exomiser.md).")
COUPLE_NOTE = ("A screen, not a diagnosis: the couple report (couple_<a>_<b>.md in the family folder) lists the "
               "genes, variant types and repeats it does not see. See a genetic counsellor before acting on a "
               "result, or on its absence (docs/39-couple.md).")
PRS_ANCESTRY_NOTE = ("Most PRS scores were developed in European-ancestry populations and predict less well for "
                     "other ancestries, even with an ancestry-adjusted percentile.")


def combined_cell(h):
    """A ClinVar hit's combined-allele note (step 6), "" when it has none."""
    note = h.get("combined_allele", ".")
    return "" if note in (".", "") else note


def cdna_text(h):
    """Step 06's cDNA read check of one ClinVar hit, as both reports print it."""
    c = h["cdna"]
    if c["check"] == "flagged":
        return (f"{c['note']} ({c['alt_clipped']} of {c['alt_reads']} alt reads and {c['ref_clipped']} of "
                f"{c['ref_reads']} ref reads clipped at an exon edge)")
    if c["check"] == "clear":
        return "clear"
    return "not checked" + (f": {c['note']}" if c["note"] else "")


def af_text(h):
    """A ClinVar hit's gnomAD allele frequency (from step 13's VEP output), as
    both reports print it; "" for a summary written before the field."""
    if "af" not in h:
        return ""
    if h["af"] is None:
        return "gnomAD AF unknown"
    txt = f"gnomAD AF {h['af']:.3g}"
    return txt + ("; common: unlikely to be a high-penetrance finding" if h.get("common") else "")


def call_text(h, floor=None):
    """The low-confidence mark of a ClinVar or CPSR row, "" for a confident one."""
    if not h.get("low_confidence"):
        return ""
    vals = ", ".join(f"{k} {h[k.lower()]}" for k in ("GQ", "DP") if h.get(k.lower()) is not None)
    return f"low-confidence call ({vals}): check the reads"


def high_rare(sf):
    """The secondary-findings HIGH-impact records counted as rare: neither a
    weak call nor of unknown frequency (a summary written before those marks
    counts every record)."""
    return sf.get("high_impact_rare", len(sf["high_impact"]))


def high_apart(sf):
    """The HIGH-impact records listed but not counted as rare, "" without any."""
    parts = []
    if sf.get("high_impact_weak"):
        parts.append(f"{sf['high_impact_weak']} weak call(s)")
    if sf.get("high_impact_frequency_unknown"):
        parts.append(f"{sf['high_impact_frequency_unknown']} of unknown frequency")
    return f" (listed apart: {', '.join(parts)})" if parts else ""


def high_marks(x, lead):
    """The marks of one HIGH-impact record (unknown frequency, weak call), each
    after lead ("  (" closes with ")"), "" for a rare record with a sound call."""
    marks = []
    if x.get("frequency_unknown"):
        marks.append("frequency unknown, not counted as rare")
    if x.get("weak_call"):
        marks.append(f"weak call ({x.get('call_note')}): check the reads")
    if not marks:
        return ""
    return lead + "; ".join(marks) + (")" if lead.endswith("(") else "")


def clinvar_floor_line(d):
    """How many ClinVar hits fall below step 6's call floor, or None for a
    summary written before the floor."""
    if "low_confidence" not in d:
        return None
    f = d.get("call_floor") or {}
    return (f"Low-confidence calls (GQ below {f.get('GQ')} or DP below {f.get('DP')}), listed but not counted "
            f"above: {d['low_confidence']}")


def clinvar_freq_line(d):
    """Where the hits' gnomAD frequencies came from, or why there are none."""
    fr = d.get("frequency")
    if not fr:
        return None
    if fr.get("note"):
        return f"gnomAD AF: unknown for every hit ({fr['note']})"
    common = sum(1 for h in d.get("hits", []) if h.get("common"))
    return f"gnomAD AF from {fr['source']}: {common} hit(s) above {fr['common_af']:.0%}"


def cpsr_row_text(x):
    """One CPSR P/LP or disagreement row: variant, gene, the classes and codes."""
    parts = [f"final {x['final']}" + (f" ({x['authority']})" if x.get("authority") else "")]
    if x.get("clinvar"):
        parts.append(f"ClinVar {x['clinvar']}")
    if x.get("cpsr"):
        parts.append(f"CPSR's own {x['cpsr']}" + (f" ({x['acmg_codes']})" if x.get("acmg_codes") else ""))
    return f"{x['gene']} {x['variant']}: " + ", ".join(parts)


def cpsr_lines(d, dbs):
    """The CPSR lines both reports share beyond the breakdown: the two ClinVar
    releases, CPSR's own benign calls among the final P/LP rows, the rows where
    CPSR and ClinVar disagree, the low-confidence rows and the cDNA flag."""
    out = []
    if "plp" not in d:
        return out
    out.append(f"ClinVar release in CPSR's bundle: {dbs.get('cpsr_clinvar_release') or 'unknown'}; "
               f"step 6 used the ClinVar file of {dbs.get('clinvar_release') or 'unknown'}")
    if d.get("own_class_column"):
        out.append(f"Final P/LP rows that CPSR's own rules class benign or likely benign: {d['final_plp_own_benign']}")
        out.append(f"Rows where ClinVar and CPSR's own class disagree on P/LP: {len(d['disagreements'])}")
        out += ["  " + cpsr_row_text(x) for x in d["disagreements"]]
    else:
        out.append("The classification table has no CPSR_CLASSIFICATION column, so no disagreement can be shown.")
    out.append(f"P/LP rows below the depth floor (DP below {d['call_floor_dp']}): {d['low_confidence']}")
    out += ["  " + cpsr_row_text(x) + f"; {call_text(x)}" for x in d["plp"] if x.get("low_confidence")]
    if "cdna_flagged" in d:
        out.append(f"cDNA read check of the P/LP rows (alt reads clipped at exon edges): {d['cdna_flagged']} flagged")
        out += ["  " + cpsr_row_text(x) + f"; reads: {cdna_text(x)}" for x in d["plp"] if x["cdna"]["check"] == "flagged"]
    return out


def hla_text(locus):
    """One HLA locus: each allele with its T1K quality, and the low-confidence mark.
    A quality-0 allele with more than collect_summary.HLA_MAX_CANDIDATES candidates
    reads 'no confident call'; the HFE row carries its label."""
    q = locus.get("quality") or []
    n = locus.get("candidates") or []
    parts = []
    for i, a in enumerate(locus["alleles"]):
        qi = q[i] if i < len(q) else ""
        if (i < len(n) and n[i] > collect_summary.HLA_MAX_CANDIDATES and qi != ""
                and collect_summary.hla_quality(qi) <= 0):
            parts.append(f"no confident call ({n[i]} candidates, quality {qi})")
        else:
            parts.append(a + (f" (quality {qi})" if qi != "" else ""))
    txt = " / ".join(parts) or "no call"
    if locus.get("low_confidence"):
        txt += "; low confidence (T1K quality 0 or below)"
        if locus.get("withheld_from_pharmcat"):
            txt += ", not passed to PharmCAT"
    if locus.get("label"):
        txt += f" [{locus['label']}]"
    return txt


HLA_DERIVED_NOTE = ("Information about alleles, not a diagnosis: each line is an association reported in the cited "
                    "study. Carrying an allele does not mean the condition or reaction will happen, and not carrying "
                    "it does not rule it out. The genotypes are unphased.")


KIR_NOTE = ("KIR genotypes as T1K typed them. No step interprets them: they matter for transplant matching and "
            "some immune conditions, not for drug guidance.")
OT_AGREE = {"yes": "agree", "no": "disagree"}


def hla_optitype_lines(ot):
    """Step 08's OptiType check, both reports: (gene, line) per class I gene,
    or one ('', line) with the reason the callers were not compared."""
    if not ot["genes"]:
        return [("", f"not compared: {ot['reason']}")]
    out = []
    for g in ot["genes"]:
        line = f"{OT_AGREE.get(g['agree'], g['agree'] or 'not compared')}: T1K {g['t1k']}, OptiType {g['optitype']}"
        if g["note"] not in ("", "-"):
            line += f" ({g['note']})"
        out.append((g["gene"], line))
    return out


def kir_lines(kir):
    """Step 08's KIR genes, both reports: (gene, alleles) per typed gene, then
    ('', line) for the genes with no allele, or the reason nothing was typed."""
    if kir["reason"]:
        return [("", kir["reason"])]
    out = [(g["gene"], " / ".join(a + (f" (quality {q})" if q != "" else "")
                                  for a, q in zip(g["alleles"], g["quality"])))
           for g in kir["typed"]]
    if kir["not_called"]:
        out.append(("", "No allele called (absent, or too few reads): " + ", ".join(kir["not_called"])))
    return out


def hla_derived_lines(d):
    """The derived HLA lines both reports print: the coeliac DQ line, the
    DQB1*02 dose and the allele table (bin/hla_alleles.tsv), present ones first."""
    if not d:
        return []
    c = d.get("coeliac") or {}
    if c.get("status") == "derived":
        dq25 = c["dq2_5"] + (f", {c['dq2_5_phase']}" if c.get("dq2_5_phase") else "")
        out = [f"Coeliac DQ: DQ2.5 {dq25}; DQ2.2 {c['dq2_2']}; DQ8 {c['dq8']}"]
    else:
        out = [f"Coeliac DQ: not determined ({c.get('reason') or 'no DQA1 and DQB1 call'})"]
    dose = c.get("dqb1_02_dose")
    out.append(f"DQB1*02 copies: {dose if dose is not None else 'not determined'}")
    order = {"present": 0, "uncertain": 1, "absent": 2, "not_determined": 3}
    for a in sorted(d.get("alleles") or [], key=lambda x: order.get(x["state"], 9)):
        state = a["state"].replace("_", " ") + (f" ({a['reason']})" if a.get("reason") else "")
        out.append(f"{a['allele']} {state}: {a['association']} ({a['citation']})")
    return out


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

# Why step 32 could not call a gene on this reference.
NOT_CALLABLE_WHY = {"GSTT1": "GRCh38 has it only on the ALT contig chr22_KI270879v1_alt, "
                             "which a reference without ALT contigs lacks"}


def probable_text(p):
    """The probable CYP2D6 call (pgx_outside_calls.probable_cyp2d6) as one line."""
    return (f"{p.get('diplotype')}" + (f", {p['phenotype']}" if p.get("phenotype") else "")
            + "; not passed to PharmCAT, so no drug guidance rests on it")


def pypgx_lines(d):
    """Step 32's lines below its counts, the same in both reports: the genes
    this reference cannot hold, the calls resting on a low-GQ variant, and the
    CYP2D6 depth check with its MAPQ shares."""
    lines = []
    for g in d.get("not_callable", []):
        lines.append(f"Not assessable on this reference: {g} ({NOT_CALLABLE_WHY.get(g, 'see docs/32-pypgx.md')})")
    for x in d.get("low_gq", []):
        lines.append(f"Low-confidence: {x['gene']} {x['diplotype']} rests on {x['weakest'] or 'a low-GQ variant'}: "
                     "check the reads")
    dc = d.get("depth_check")
    if dc:
        nm = "not measured"
        lines.append(f"CYP2D6 depth check: {dc.get('status') or 'unreadable'}; share of the depth at MAPQ >= 1 / >= 20: "
                     f"Cyrius's region with REP6 {dc.get('gene_mapq1_share', nm)} / {dc.get('gene_mapq20_share', nm)}, "
                     f"gene body without REP6 {dc.get('body_mapq1_share', nm)} / {dc.get('body_mapq20_share', nm)}")
    return lines


def pgx_copy_number(S):
    """(gene, copy number) rows, the same in both reports: pypgx's call from
    read depth for its BAM-based genes, and UGT2B17 from CNVpytor's calls
    over its GENCODE span."""
    pg, cn = S.get("pypgx", MISSING), S.get("cnvpytor", MISSING)
    if pg["state"] == "missing" and cn["state"] == "missing":
        return []
    rows = [(g, f"{v} (pypgx, read depth)") for g, v in (pg["data"].get("cnv_call") or {}).items()]
    u = cn["data"].get("ugt2b17") if cn["state"] in ("ok", "stale") else None
    if u:
        rows.append(("UGT2B17", f"{u['call']} ({u['note']})"))
    else:
        rows.append(("UGT2B17", f"not assessed (CNVpytor, step 18, {'failed' if cn['state'] == 'failed' else 'did not run' if cn['state'] == 'missing' else 'gave no UGT2B17 result'})"))
    return rows


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


MOSAIC_NOTE = ("Loss of chrY in a share of the blood cells of men, and of one chrX in women, becomes more common "
               "with age. Population studies link loss of Y with smoking and with higher rates of heart disease "
               "and some cancers; loss of X is less studied. That is a population-level association, not a "
               "diagnosis for one person. Below about 10% of the cells nothing can be told at 30x "
               "(docs/33-sample-qc.md).")
MOSAIC_CALL = {"detected": ("detected", "yellow"), "not_detected": ("not detected above 10%", "gray"),
               "gain": ("an extra copy, not a loss", "yellow"), "not_estimated": ("not estimated", "gray")}


def pct(x):
    return f"{x * 100:.1f}%" if isinstance(x, (int, float)) else "?"


def mosaic_lines(q):
    """The estimate from depth and from chrX allele fractions, one line each, as both reports print them."""
    what = "the Y" if q.get("chrom") == "chrY" else "one X"
    out = []
    if isinstance(q.get("loss_depth"), (int, float)):
        ratio = q.get("depth_ratio")
        ratio = f"{ratio:.3f}" if isinstance(ratio, (int, float)) else "?"
        out.append(f"From depth: estimated share of the cells that lost {what} {pct(q['loss_depth'])} (95% range "
                   f"{pct(q.get('loss_depth_low'))} to {pct(q.get('loss_depth_high'))}); {q.get('chrom')} to "
                   f"autosome depth ratio {ratio}, {q.get('depth_ratio_expected')} expected")
    if isinstance(q.get("loss_af"), (int, float)):
        out.append(f"From chrX allele fractions at {q.get('x_het_sites')} heterozygous sites: {pct(q['loss_af'])} "
                   f"(95% range {pct(q.get('loss_af_low'))} to {pct(q.get('loss_af_high'))})")
    return out


def read_count_text(d):
    """Step 02's two counts as one line, in both reports."""
    txt = f"{d['input_reads']} reads in, {d['primary_records']} primary records in the BAM"
    if d.get("read_count_check") == "mismatch":
        m = d["missing_reads"]
        txt += (f": {m} reads missing" if m > 0 else ": more records than reads in") + \
            " (the BAM was kept by the warn setting, see docs/troubleshooting.md)"
    return txt


def prs_chrx_note(d):
    """A sentence when chrX rows of a score were left out (no sex given), else ''."""
    left = [f"{r['pgs_id']} ({r['chrx'].split()[0]})" for r in d.get("scores", [])
            if (r.get("chrx") or "").endswith("left out")]
    if not left:
        return ""
    return (f" Rows on chrX were left out, because the sample's sex was not given: {', '.join(left)}. "
            "Give the sex to score them (docs/25-prs.md).")


PRS_NOT_APPLICABLE = "not applicable to the declared sex"


def prs_percentile(r):
    """The percentile cell of a PRS row in both reports: the percentile with
    its group, 'not applicable to the declared sex' for a score of the other
    sex (it has no percentile), or '' when there is none."""
    if r.get("sex_check") == PRS_NOT_APPLICABLE:
        return f"{PRS_NOT_APPLICABLE} ({r.get('applies_to') or 'sex'}-specific score)"
    if r.get("percentile"):
        return f"{r['percentile']} ({r.get('group') or 'group unknown'})"
    return ""


def prs_sex_note(d):
    """A sentence naming the sex-specific scores when no sex was declared, else ''."""
    named = [f"{r['pgs_id']} ({r['applies_to']})" for r in d.get("scores", [])
             if r.get("sex_check") == "no sex declared" and r.get("applies_to")]
    if not named:
        return ""
    return (f" No sex was declared, so these sex-specific scores are shown for any sex: {', '.join(named)}. "
            "Give the sex to mark a score of the other sex not applicable (docs/25-prs.md).")


def prs_note(d):
    """The line under the PRS scores, in both reports: what the numbers can be compared with."""
    x = prs_chrx_note(d) + prs_sex_note(d)
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

def couple_lines(d):
    """Step 39's lines, the same in both reports: the counts, then each shared
    and X-linked gene with the chance for a child and both partners' findings."""
    partner = d.get("partner") or "the partner"
    lines = [f"Partner: {partner}. Genes where both carry a change: {len(d.get('shared', []))}; X-linked "
             f"findings: {len(d.get('x_linked', []))}; genes where only one partner carries a change: "
             f"{d.get('one_partner', 0)}"]
    for label, key in (("Shared", "shared"), ("X-linked", "x_linked")):
        for x in d.get(key, []):
            lines.append(f"{label}: {x['gene']} ({x['condition']}): {x['risk']}. This sample: {x['own']}. "
                         f"{partner}: {x['partner_findings']}")
    return lines + [COUPLE_NOTE]


def share(x):
    """A fraction as a whole percentage, for the fine-grained ancestry section."""
    return "." if x is None else f"{100 * x:.0f}%"


def ancestry_fine_rows(d, level):
    """The groups of one level worth a line: every continental group of 1% or
    more, a regional one of 1% or more or whose bootstrap range reaches 5%."""
    return [r for r in d.get(level, []) if r["fraction"] is not None
            and (r["fraction"] >= 0.01 or (level == "region" and (r.get("high") or 0) >= 0.05))]


def ancestry_fine_spread(r):
    return f" ({share(r['low'])} to {share(r['high'])})" if r.get("low") is not None else ""


def ancestry_fine_left_out(d):
    """The groups the model leaves out, as one sentence, or ''."""
    if not d.get("left_out"):
        return ""
    # One TSV row per level ("Oceania (continental group)"); merge a group the
    # two levels share into one name, as the step's own summary does.
    kinds = {}
    for r in d["left_out"]:
        m = re.fullmatch(r"(.*) \((.*) group\)", r["group"])
        name, kind = (m.group(1), m.group(2)) if m else (r["group"], "")
        kinds.setdefault(name, [])
        if kind and kind not in kinds[name]:
            kinds[name].append(kind)
    names = [f"{n} ({' and '.join(k)} group)" if k else n for n, k in kinds.items()]
    return (f"Left out of the model: {', '.join(names)}. The panel has fewer than "
            f"{d.get('min_references', 10)} unrelated reference people for it, so its ancestry is spread over the "
            "groups shown.")


def ancestry_fine_painted(d):
    """The painted shares, or why there are none, and the chromosomes RFMix failed on."""
    failed = ", ".join(r["code"] for r in d.get("not_painted", []))
    if d.get("painting"):
        text = ", ".join(f"{r['group']} {share(r['fraction'])}" for r in d["painting"])
    else:
        text = "none" + ("" if failed else " (no chromosome long enough)")
    return text + (f"; not painted: RFMix failed on {failed}" if failed else "")


def exomiser_terms(d):
    return ", ".join(f"{t['label']} ({t['id']})" if t.get("label") else t["id"] for t in d.get("terms", [])) or "."


def exomiser_text(d, sample):
    """The Exomiser section of the text report: the terms, then each rank."""
    out = [f"  Symptoms searched: {exomiser_terms(d)}"]
    if not d["genes"]:
        out.append("  No gene passed Exomiser's filters with a variant that fits a mode of inheritance.")
    for g in d["genes"]:
        out.append(f"  {g['rank']}. {g['gene']} ({g['moi']}) score {g['exomiser_score']} (symptom match {g['phenotype_score']})")
        for v in g["variants"]:
            out.append(f"       {v['variant']} {v['genotype']} {v['consequence']} {v['hgvs']} ACMG {v['acmg']}"
                       + (f"; {v['disease']}" if v["disease"] != "." else ""))
    out.append(f"  {EXOMISER_NOTE} In words: exomiser/{sample}_exomiser.txt")
    return out


def exomiser_html(d, sample):
    """The Exomiser card's body: one table row per rank and variant."""
    rows = []
    for g in d["genes"]:
        for i, v in enumerate(g["variants"]):
            lead = (f"<td>{E(g['rank'])}</td><td>{E(g['gene'])}</td><td>{E(g['moi'])}</td><td>{E(g['exomiser_score'])}"
                    f"</td><td>{E(g['phenotype_score'])}</td>") if i == 0 else "<td></td>" * 5
            rows.append(f"      <tr>{lead}<td>{E(v['variant'])} {E(v['genotype'])}<br><span style=\"font-size:12px;"
                        f"color:#666\">{E(v['consequence'])} {E(v['hgvs'])}</span></td><td>{E(v['acmg'])}</td>"
                        f"<td>{E(v['disease'])}</td></tr>")
    body = [stat("Symptoms searched", E(exomiser_terms(d)))]
    if rows:
        body += ["    <table>", "      <tr><th>Rank</th><th>Gene</th><th>Inheritance</th><th>Score</th>"
                 "<th>Symptom match</th><th>Variant</th><th>ACMG (Exomiser)</th><th>Linked disease</th></tr>"]
        body += rows + ["    </table>"]
    else:
        body.append(stat("Ranked genes", "none passed Exomiser's filters"))
    body += [small_note(EXOMISER_NOTE), small_note(f"Exomiser's own report: exomiser/{sample}_exomiser.html")]
    return body


def traits_text(d, sample):
    """The traits section of the text report: the rows read, the ones not read,
    and the APOE lines last, after a warning."""
    out = [f"  Read: {d['called']}  No genotype at the position: {d['not_covered']}  Other allele: {d['other_allele']}"]
    for r in d["rows"]:
        if r["status"] == "not_covered":
            continue
        mark = " [ancestry]" if r["ancestry"] == "yes" else ""
        out.append(f"  {r['trait']} ({r['gene']} {r['rsid']}): {r['genotype']}{mark}. {r['reading']} (PMID {r['pmids']})")
    missing = [r["rsid"] for r in d["rows"] if r["status"] == "not_covered"]
    if missing:
        out.append(f"  Not read (no genotype in the gVCF): {', '.join(missing)}")
    out.append(f"  {TRAITS_NOTE} Notes per row: traits/{sample}_traits.tsv")
    if d.get("apoe"):
        out.append(f"  {APOE_WARNING}")
        for r in d["apoe"]:
            out.append(f"    {r['trait']} ({r['rsid']}): {r['genotype']}. {r['reading']}")
    return out


def traits_html(d, sample):
    """The traits card's body: one table row per table row read or not read,
    and the APOE lines folded away behind a warning."""
    def tr(r):
        anc = " " + badge("ancestry", "yellow") if r["ancestry"] == "yes" else ""
        return (f"<tr><td>{E(r['trait'])}</td><td>{E(r['gene'])}<br><span style=\"font-size:12px;color:#666\">"
                f"{E(r['rsid'])}</span></td><td>{E(r['genotype'])}{anc}</td><td>{E(r['reading'])}"
                f"<br><span style=\"font-size:12px;color:#666\">{E(r['note'])}</span></td><td>PMID {E(r['pmids'])}</td></tr>")
    head = "      <tr><th>Trait</th><th>Gene</th><th>Genotype</th><th>Reading</th><th>Source</th></tr>"
    body = [stat("Read / no genotype / other allele", E(f"{d['called']} / {d['not_covered']} / {d['other_allele']}")),
            "    <table>", head] + ["      " + tr(r) for r in d["rows"]] + ["    </table>", small_note(TRAITS_NOTE)]
    if d.get("apoe"):
        body += [f"    <details><summary>{E(APOE_WARNING)} Open to show it.</summary>", "    <table>", head]
        body += ["      " + tr(r) for r in d["apoe"]] + ["    </table></details>"]
    body.append(small_note(f"Every row with its note: traits/{sample}_traits.tsv"))
    return body


COMPHET_NOTE = ("A pair can hit both copies of the gene only in trans. cis: one copy carries both variants. "
                "unknown: fewer than 2 reads or read pairs cover both positions, the reads disagree, or there is no BAM; "
                "the pairs table gives the reason. A parent's genome or long reads can phase it. "
                "Every pair: slivar/<sample>_compound_het_pairs.tsv")


def comphet_phase_text(d):
    """Step 31's phase of the compound-het pairs, as text lines."""
    ph = d["compound_het_phase"]
    lines = [f"  Phase from the reads: {ph.get('trans', 0)} trans, {ph.get('cis', 0)} cis, "
             f"{ph.get('unknown', 0)} unknown of {d.get('compound_het_pairs', 0)} pairs"]
    for x in d.get("compound_het_trans") or []:
        lines.append(f"    {x['gene']:<8} {x['variants']}  in trans on {x['reads']} reads or read pairs")
    if ph.get("trans", 0) > len(d.get("compound_het_trans") or []):
        lines.append(f"    ... {ph['trans'] - len(d['compound_het_trans'])} more in the pairs table")
    return lines + [f"  ({COMPHET_NOTE})"]


def comphet_phase_html(d):
    """The same, for the slivar card."""
    ph = d["compound_het_phase"]
    out = [stat("Pairs in trans / cis / unknown", f"{ph.get('trans', 0)} / {ph.get('cis', 0)} / {ph.get('unknown', 0)}")]
    trans = d.get("compound_het_trans") or []
    if trans:
        out += ["    <table>", "      <tr><th>Gene</th><th>Pair in trans</th><th>Reads</th></tr>"]
        out += [f"      <tr><td>{E(x['gene'])}</td><td>{E(x['variants'])}</td><td>{E(x['reads'])}</td></tr>"
                for x in trans]
        out.append("    </table>")
    return out + [small_note(COMPHET_NOTE)]


ANNOTSV_NOTE = ("A class 4 or 5 SV made by one caller, a deletion on the exon edges of a transcript (the reads of "
               "a processed transcript or cDNA aligned across an intron) or one at 1% or more in AnnotSV's benign "
               "set is unlikely to be a rare pathogenic SV: check the reads before acting on it (docs/05-annotsv.md).")


def annotsv_support_text(x):
    if x["support"] is None:
        return "not known (no step 22 consensus)"
    if x["support"] == 1:
        return f"single caller ({x['caller']})"
    return f"{x['support']} callers ({', '.join(x['support_callers'])})"


def annotsv_genes_text(x):
    g = x["genes"]
    return ", ".join(g[:3]) + (f" and {len(g) - 3} more" if len(g) > 3 else "") if g else "no gene"


def annotsv_flags(x):
    """The exon-edge and population flags of one class 4-5 SV, as a list."""
    out = []
    if x["exon_edge"] == "yes":
        out.append("on exon edges: " + x["exon_edge_detail"])
    elif x["exon_edge"] == "not checked":
        out.append("exon edges not checked")
    if x["population_af"] is not None:
        out.append(f"benign-set AF {x['population_af']:g}" + (f" ({x['population_source']})" if x["population_source"] else "")
                   + (", common" if x["common"] else ""))
    return out


def annotsv_lines(d):
    """Step 05's class 4-5 lines, the same in both reports."""
    head = f"Class 4-5 SVs: {d['count']}"
    if d["count"]:
        head += (f" (single caller: {d['single_caller'] if d['consensus'] else 'not known'}; deletions on exon edges: "
                 f"{d['exon_edge'] if d['exon_edge_checked'] else 'not checked'}; common: {d['common']})")
    lines = [head]
    if d["count"] and not d["consensus"]:
        lines.append(d.get("consensus_note") or "Caller support: step 22's consensus did not run, so how many callers "
                     "saw each SV is not known.")
    return lines


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
    rdc = S.get("read_count", MISSING)
    if any(sec["state"] != "missing" for sec in (cov, sexc, sqc, rdc)):
        w("## Quality control")
        w("---")
        for sec in (cov, sexc, sqc, rdc):
            if stale_note(sec):
                w(f"  [{sec['title']}: {stale_note(sec)}]")
        md = cov["data"].get("mean_depth")
        w(f"  Mean depth: {md:.1f}x" if isinstance(md, (int, float)) else "  Mean depth: not available")
        if rdc["state"] in ("ok", "stale"):
            flag = "  ** READS LOST IN ALIGNMENT **" if rdc["data"].get("read_count_check") == "mismatch" else ""
            w(f"  Alignment: {read_count_text(rdc['data'])}{flag}")
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

    if head("mosaic_sex_loss"):
        q = S["mosaic_sex_loss"]["data"]
        w(f"  {q.get('chrom') or 'chrY or chrX'}: {MOSAIC_CALL.get(q.get('call'), (q.get('call'), ''))[0]}"
          f" ({q.get('reason')})")
        for line in mosaic_lines(q):
            w(f"  {line}")
        w(f"  {MOSAIC_NOTE}")
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
        low = clinvar_floor_line(d)
        if low:
            w(f"  {low}")
        w(f"  ClinVar file date: {s['databases']['clinvar_release'] or 'unknown'}")
        if d["hits"]:
            by = d["by_stars"]
            w("  By review status: " + ", ".join(f"{k} star{'s' if k != '1' else ''} {by[k]}" for k in ("4", "3", "2", "1", "0")))
            if "cdna_flagged" in d:
                w(f"  cDNA read check (alt reads clipped at exon edges): {d['cdna_flagged']} flagged")
            freq = clinvar_freq_line(d)
            if freq:
                w(f"  {freq}")
            w("")
            w("  Hits (best-reviewed first, up to 50): gene, position, genotype, significance [review status] stars")
            for h in d["hits"][:50]:
                note = combined_cell(h)
                w(f"    {h['gene']:<10} {h['chrom']}:{h['pos']} {h['ref']}>{h['alt']}  {h['genotype']}  "
                  f"{h['significance']} [{h['review_status']}] {h['stars']}*" + (f"  ({note})" if note else ""))
                if af_text(h):
                    w(f"               frequency: {af_text(h)}")
                if call_text(h):
                    w(f"               call: {call_text(h)}")
                if "cdna" in h:
                    w(f"               reads: {cdna_text(h)}")
        w(f"  {CLINVAR_FREQ_NOTE}")
        w("")

    if head("couple", "Couple Carrier Screen (step 39)"):
        for line in couple_lines(S["couple"]["data"]):
            w(f"  {line}")
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
        for line in pypgx_lines(d):
            w(f"  {line}")
        w("")

    cn_rows = pgx_copy_number(S)
    if cn_rows:
        w("## PGx copy number")
        w("---")
        for g, v in cn_rows:
            w(f"  {g:<8} {v}")
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
        if s["cyp2d6"].get("probable"):
            w(f"  {s['cyp2d6']['probable']['label'].capitalize()}: {probable_text(s['cyp2d6']['probable'])}")
        w("")

    if head("paralogs", "Paralog Genes (Paraphase, PacBio HiFi)"):
        for line in paralog_lines(S["paralogs"]["data"]):
            w(f"  {line}")
        w(f"  {PARALOG_NOTE}")
        w("")

    if head("parascopy", "SMN1/SMN2 Copy Number (Parascopy, short reads)"):
        for line in parascopy_lines(S["parascopy"]["data"]):
            w(f"  {line}")
        w(f"  {PARASCOPY_NOTE}")
        w("")
    unchecked = paralogs_unchecked(S)
    if unchecked:
        w("## Paralog Copy Number")
        w("---")
        w(f"  {unchecked}")
        w("")

    if head("hla", "HLA Typing (T1K)"):
        for l in S["hla"]["data"]["loci"]:
            w(f"    {l['gene']:<8} {hla_text(l)}")
        for line in hla_derived_lines(S["hla"]["data"].get("derived")):
            w(f"  {line}")
        if S["hla"]["data"].get("derived"):
            w(f"  {HLA_DERIVED_NOTE}")
        w(f"  HLA database: {s['databases']['hla_database'] or 'unknown'}")
        w(f"  {HLA_NOTE}")
        ot = S["hla"]["data"].get("optitype")
        if ot is not None:
            w("  OptiType check (HLA_QC=true, two fields):")
            for g, line in hla_optitype_lines(ot):
                w(f"    {g:<8} {line}" if g else f"    {line}")
        kir = S["hla"]["data"].get("kir")
        if kir is not None:
            w("  KIR genes (T1K, KIR=true):")
            for g, line in kir_lines(kir):
                w(f"    {g:<8} {line}" if g else f"    {line}")
            w(f"  {KIR_NOTE}")
        w("")

    if head("prs", "Polygenic Risk Scores"):
        d = S["prs"]["data"]
        for r in d["scores"]:
            pct = prs_percentile(r)
            pct = f"  {pct}" if r.get("sex_check") == PRS_NOT_APPLICABLE else (f"  percentile {pct}" if pct else "")
            w(f"  {r['condition']:<35} {r['score']} ({r['matched']}/{r['total']} variants matched) {r['pgs_id']}{pct}")
        w("  " + prs_note(d))
        w("")

    if head("ancestry_fine", "Fine-Grained Ancestry (ADMIXTURE, RFMix)"):
        d = S["ancestry_fine"]["data"]
        for level, title in ANCESTRY_FINE_LEVELS:
            rows = ancestry_fine_rows(d, level)
            if rows:
                w(f"  {title} (fraction, bootstrap range):")
                for r in rows:
                    w(f"    {r['group']:<50} {share(r['fraction']):>5}{ancestry_fine_spread(r)}")
        if ancestry_fine_left_out(d):
            w(f"  {ancestry_fine_left_out(d)}")
        w(("  Chromosome painting, share of the painted genome: " if d.get("painting") else "  Chromosome painting: ")
          + ancestry_fine_painted(d))
        w(f"  {ANCESTRY_FINE_NOTE}")
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
        filt = cnvpytor_filtered_text(d)
        if filt:
            w(f"  Significant after filters: {filt}")
            w(f"  {CNVPYTOR_NOTE}")
        w("")
    if head("sv_consensus", "SV Consensus Merge"):
        w(f"  Consensus SVs (2+ callers): {S['sv_consensus']['data']['consensus']}")
        w("")
    if head("annotsv", "Structural Variants, Class 4 and 5 (AnnotSV)"):
        d = S["annotsv"]["data"]
        for line in annotsv_lines(d):
            w(f"  {line}")
        for x in d["rows"]:
            w(f"    {x['type']:<4} {x['sv']}  class {x['acmg_class']}  {annotsv_genes_text(x)}  "
              f"support: {annotsv_support_text(x)}")
            for f in annotsv_flags(x):
                w(f"         {f}")
        if d["count"]:
            w(f"  {ANNOTSV_NOTE}")
        w("")

    if head("expansions", "Repeat Expansions (ExpansionHunter)"):
        d = S["expansions"]["data"]
        w(f"  Loci in output: {d['records']}")
        if d.get("stranger_outdated"):
            w("  Stranger's result is older than ExpansionHunter's: rerun step 9b. No locus is marked as normal or expanded.")
        elif not d.get("stranger"):
            w("  Stranger (step 9b) did not run: no locus is marked as normal or expanded.")
        elif d.get("flagged"):
            w(f"  Outside the normal range (Stranger): {len(d['flagged'])}")
            by = {x["locus"]: x for x in d.get("loci") or []}
            for l in d["flagged"]:
                w(f"    {l['locus']:<8} {str_count(by.get(l['locus'], l))}  {l['status']}")
                if l.get("rfc1_motif"):
                    w(f"      {rfc1_text(l['rfc1_motif'])}")
        else:
            w("  Outside the normal range (Stranger): none")
        if d.get("stranger") and d.get("no_status"):
            w(f"  Loci without a Stranger status (not in its catalog): {d['no_status']}")
        if d.get("unsized"):
            w(f"  Not sized from short reads, so not counted as normal: {len(d['unsized'])} ({', '.join(d['unsized'])})")
        if d.get("intermediate"):
            w(f"  In an intermediate range (bin/str_intermediate_ranges.tsv): {', '.join(d['intermediate'])}")
        rows = str_rows(d)
        if rows:
            w("  Every locus (key loci first): repeat count, Stranger status, notes")
            for locus, count, status, notes in rows:
                w(f"    {locus:<8} {count:<10} {status or '-'}" + (f"  {notes}" if notes else ""))
        if any(l.get("rfc1_motif") for l in d.get("flagged") or []):
            w(f"  {RFC1_NOTE}")
        w(f"  {EH_NOTE}")
        w("")

    if head("telomere", "Telomere content (relative, TelomereHunter)"):
        w(f"  Telomere content: {S['telomere']['data']['tel_content']} {TELOMERE_UNIT}")
        for line in telomere_lines(S["telomere"]["data"]):
            w(f"  {line}")
        w(f"  {TELOMERE_NOTE}")
        w("")

    if head("roh", "Runs of Homozygosity"):
        d = S["roh"]["data"]
        w(f"  Segments: {d['segments']}  Total: {d['total_mb']} MB  Largest: {d['largest_mb']} MB")
        w(f"  Autosomal ROH > 5MB: {len(d['autosomal_over_5mb'])}")
        for r in d["autosomal_over_5mb"]:
            w(f"    {roh_region(r)}")
        for line in roh_lines(d):
            w(f"  {line}")
        w("")

    if head("haplogroup", "Mitochondrial Haplogroup"):
        d = S["haplogroup"]["data"]
        w(f"  Haplogroup: {d['haplogroup']}")
        w(f"  Contamination (haplocheck): {contamination_text(d)}")
        for line in explore_lines(d):
            w(f"  {line}")
        w("")

    if head("y_haplogroup", "Y-Chromosome Haplogroup (Yleaf)"):
        d = S["y_haplogroup"]["data"]
        w(f"  Haplogroup: {d['haplogroup']} ({d['valid_markers']} markers, QC-score {d['qc_score']})")
        for line in explore_lines(d):
            w(f"  {line}")
        w("")

    if head("y_str", "Y-STR Markers (HipSTR)"):
        d = S["y_str"]["data"]
        w(f"  {ystr_counts(d)}")
        for x in d["markers"]:
            if x["count"]:
                w(f"  {x['panel']:<6} {x['marker']:<13} {x['count']:<6} ({x['confidence']})")
        w(f"  {YSTR_NOTE}")
        w("")

    if head("consumer_export", "Raw-Data Export (23andMe v5 layout)"):
        for line in consumer_export_lines(S["consumer_export"]["data"]):
            w(f"  {line}")
        w("")

    if head("bloodtype", "Blood Groups (predicted from the genome)"):
        for x in S["bloodtype"]["data"]["systems"]:
            w(f"  {x['system']:<9} {x['phenotype']:<34} {x['genotype']:<28} ({x['confidence']})")
        w(f"  {BLOODTYPE_NOTE}")
        w("")

    if head("lpa", "Lp(a): LPA KIV-2 Repeat Count and Risk Variants"):
        for line in lpa_lines(S["lpa"]["data"]):
            w(f"  {line}")
        w("")

    if head("chip", "Clonal Haematopoiesis (CHIP, step 49)"):
        out.extend(f"  {line}" for line in chip_lines(S["chip"]["data"], s["sample"]))
        w("")

    if head("mito", "Mitochondrial Variants (Mutect2)"):
        d = S["mito"]["data"]
        w(f"  PASS variants: {d['pass']}")
        w(f"  Heteroplasmic (AF 0.05 to 0.95, outside the D-loop poly-C tracts): {d['heteroplasmic']}")
        for line in mito_lines(d):
            w(f"  {line}")
        if d.get("polyc_heteroplasmic"):
            w(f"  {MITO_POLYC_NOTE}")
        w(f"  {MITO_NOTE}")
        w("")

    if head("traits", "Traits Panel (step 42)"):
        out.extend(traits_text(S["traits"]["data"], s["sample"]))
        w("")

    if head("viral", "Viral DNA in the Unmapped Reads (step 47)"):
        for line in viral_lines(S["viral"]["data"]):
            w(f"  {line}")
        w("")

    if head("exomiser", "Symptom-Driven Ranking (Exomiser, step 48)"):
        out.extend(exomiser_text(S["exomiser"]["data"], s["sample"]))
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
        for line in cpsr_lines(d, s["databases"]):
            w(f"  {line}")
        w("")

    if head("clinical", "Clinical Variant Filter"):
        d = S["clinical"]["data"]
        w(f"  Clinical variants: {d['variants']}" + (f" in {d['genes']} genes" if "genes" in d else ""))
        if d.get("by_impact"):
            w("  By impact: " + ", ".join(f"{k} {v}" for k, v in sorted(d["by_impact"].items())))
        if "frequency_unknown" in d:
            w(f"  Frequency unknown (not counted as rare): {d['frequency_unknown']}")
        if "weak_calls" in d:
            w(f"  Weak genotype calls (check the reads): {d['weak_calls']}")
        sf = d.get("acmg_sf")
        if sf:
            w(f"  Secondary-findings genes ({sf['version']}, {sf['genes_on_list']} genes): "
              f"{len(sf['clinvar_hits'])} ClinVar P/LP, {high_rare(sf)} rare HIGH impact" + high_apart(sf))
            for x in sf["clinvar_hits"]:
                w(f"    {x['gene']:<8} {x['variant']}  {x['genotype']}  ClinVar {x['significance']}"
                  + ("  (low-confidence call: check the reads)" if x.get("low_confidence") else ""))
            for x in sf["high_impact"]:
                w(f"    {x['gene']:<8} {x['variant']}  {x['genotype']}  {x['consequence']}" + high_marks(x, "  ("))
            w(f"    {ACMG_NOTE}")
        w("")

    if head("slivar", "Variant Prioritization (slivar)"):
        d = S["slivar"]["data"]
        w(f"  Prioritized variants: {d['prioritized']}")
        if "frequency_unknown" in d:
            w(f"  Frequency unknown (prioritized, not counted as rare): {d['frequency_unknown']}")
        if "compound_het_variants" in d:
            w(f"  Compound het candidates: {d['compound_het_variants']} variants across {d['compound_het_genes']} genes")
        if "compound_het_phase" in d:
            for line in comphet_phase_text(d):
                w(line)
        elif "compound_het_variants" in d:
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


STR_SIZING_TEXT = {
    "flanking_only": "size not determined (flanking reads only)",
    "in_repeat": "longer than a read (in-repeat reads only): shown as at least the interval's lower end",
    "low_confidence": "low confidence (fewer than 2 spanning reads, or the interval reaches the disease range)",
    "no_call": "no call",
}
RFC1_NOTE = ("RFC1: Stranger's catalog sizes the degenerate AARRG unit and flags a low count, so the motif "
             "decides. The motif count reads the reads at the locus; it cannot size a long expansion.")


def str_count(l):
    """The repeat count of a locus: an in-repeat allele as >= the lower end of
    its interval, or 'not sized' without one."""
    counts = l["repeat_count"].split("/")
    sizing = (l.get("allele_sizing") or "").split("/")
    cis = (l.get("ci") or "").split("/")
    out = []
    for i, c in enumerate(counts):
        z = sizing[i] if i < len(sizing) else ""
        ci = cis[i] if len(cis) == len(counts) else (cis[0] if len(cis) == 1 else "")
        lo = ci.split("-")[0] if ci else ""
        if z == "in_repeat":
            c = f">={lo}" if lo.isdigit() else "not sized"
        out.append(c)
    return "/".join(out)


def str_status(l):
    """Stranger's status, never 'normal' for a locus that was not sized."""
    st = l.get("status") or ""
    if st == "normal" and l.get("sizing") not in ("sized", "unknown"):
        return "not sized"
    return st


def str_notes(l):
    """What the report says beside a locus: sizing, the catalog's normal maximum, an intermediate range."""
    out = []
    if l.get("sizing") in STR_SIZING_TEXT:
        out.append(STR_SIZING_TEXT[l["sizing"]] + (f", interval {l['ci']}" if l.get("ci") not in (None, "", ".") else ""))
    if l.get("normal_max"):
        out.append(f"normal max {l['normal_max']}")
    if l.get("zone"):
        out.append(f"{l['zone']}: {l['zone_label']}; {l['zone_citation']}")
    return "; ".join(out)


def str_rows(d):
    """(locus, count, status, notes) for every locus, key loci first; a key locus missing from the output says so."""
    rows = [(l["locus"], str_count(l), str_status(l), str_notes(l)) for l in d.get("loci") or []]
    have = {r[0] for r in rows}
    missing = [(k["locus"], "not in output", "", "") for k in d.get("key_loci") or [] if k["locus"] not in have]
    return missing + rows


def rfc1_text(m):
    return f"RFC1 motif count: {m.get('verdict', '.')}: {m.get('note', '')}"


PARALOG_NOTE = ("Paraphase reads PacBio HiFi reads only. A 'no call' means it saw one haplotype and the depth "
                "did not settle one copy or two. A research estimate, not a diagnostic test: confirm a result "
                "(for example an HBA deletion or one SMN1 copy) with a clinical laboratory.")


def smn1_text(m):
    """SMN1 from Paraphase beside Parascopy (step 35), with the verdict."""
    pc = m.get("parascopy")
    if pc is None:
        return None
    pp = m.get("paraphase")
    verdict = {"agree": "the two agree", "disagree": "they DISAGREE: check both before using either"}.get(
        m.get("agree"), "no comparison: one of them has no call")
    return f"Paraphase {pp if pp is not None else 'no call'}, Parascopy {pc['text']}: {verdict}"


def paralog_lines(d):
    """Lines of the text report's paralog section."""
    out = [f"{g['gene']}: {g['copy_number']} ({g['status']})" + ("" if g["detail"] in ("", ".") else f"; {g['detail']}")
           for g in d.get("genes", [])]
    smn = smn1_text(d.get("smn1") or {})
    if smn:
        out.append(f"SMN1 by both tools: {smn}")
    if d.get("other_regions"):
        out.append(f"Other Paraphase regions: {d['other_regions']} ({d['other_regions_called']} with a call)")
    out.append(f"Variants on phased haplotypes: {sum(d.get('variants', {}).values())}")
    return out


PARASCOPY_NOTE = ("Parascopy reads short reads at SMN1 and SMN2 only; CYP21A2 has no Parascopy model and is not "
                  "assessed. A quality under 20 or a filter other than PASS leaves a copy unresolved. HBA1 and "
                  "HBA2 have no Parascopy model either: step 35 counts the alpha-globin genes together from read "
                  "depth against flanking sequence, which cannot say which gene or which chromosome lost a copy. "
                  "A research estimate, not a diagnostic test: confirm a result with a clinical laboratory (MLPA "
                  "or qPCR; gap-PCR for alpha-thalassaemia).")


def parascopy_lines(d):
    """Step 35's lines, the same in both reports."""
    smn1, smn2 = d.get("smn1") or {}, d.get("smn2") or {}
    return [f"SMN1 copies: {smn1.get('text', '.')}",
            f"SMN2 copies: {smn2.get('text', '.')}",
            f"SMN1 and SMN2 together (agCN) at exon 7: {d.get('agcn', '.')} (quality {d.get('agcn_qual', '.')}, "
            f"filter {d.get('agcn_filter', '.')})",
            f"Reading: {d.get('reading', '.')}"] + ([hba_line(d["hba"])] if d.get("hba") else [])


def hba_line(h):
    """Step 35's HBA1/HBA2 count from depth, one line."""
    if h.get("status") == "called":
        value = f"{h.get('alpha_globin_genes', '.')} ({h.get('message', '.')}; estimate {h.get('estimate', '.')})"
    else:
        value = h.get("message", "not assessed")
    return f"HBA1 and HBA2 together (alpha-globin genes, from depth): {value}"


def paralogs_unchecked(S):
    """'Paralog copy number: not checked (...)' when neither step 35 nor step
    35b left a result, so a report without them never reads as normal."""
    states = {k: S.get(k, MISSING)["state"] for k in ("parascopy", "paralogs")}
    if any(v in ("ok", "stale") for v in states.values()):
        return None
    if all(v == "missing" for v in states.values()):
        why = "neither step 35 nor step 35b ran"
    else:
        words = {"missing": "not run", "failed": "failed", "unreadable": "unreadable"}
        why = ", ".join(f"step {st} {words.get(states[k], states[k])}" for k, st in (("parascopy", "35"), ("paralogs", "35b")))
    return f"Paralog copy number: not checked ({why})"


ACMG_NOTE = ("A list to review with a clinician: the ACMG reporting rules per gene (for example HFE "
             "homozygous C282Y only, BTD and CYP27A1 two variants) are not applied here.")


def contamination_text(d):
    if "contamination_status" not in d:
        return "not checked (needs step 20's Mutect2 chrM calls)"
    s = d["contamination_status"]
    word = {"YES": "YES, two mtDNA haplogroups in the reads (another person's DNA?)", "NO": "no",
            "ND": "not determined"}.get(s, s)
    return f"{word} (level {d.get('contamination_level', '.')})"


YSTR_NOTE = ("Repeat counts read from the genome with HipSTR, in FamilyTreeDNA's Y-111 order. A count can "
             "differ from FamilyTreeDNA's by a fixed number of repeats at some markers, as labs count some repeat "
             "blocks differently: compare with a known result before trusting a difference. Multi-copy markers and "
             "repeats longer than the reads are not called. Every marker, with its reads and notes, is in "
             "y_str/<sample>_y_str.tsv; y_str/<sample>_y_str.csv is the same list for a spreadsheet.")


def ystr_counts(d):
    return (f"Markers called: {d['high']} high and {d['low']} low confidence of {len(d['markers'])}; "
            f"{d['no_call']} no call, {d['not_callable']} not callable from short reads")


def ystr_body(d):
    """The Y-STR card: the counts, then every called marker in FamilyTreeDNA's order."""
    rows = "".join(f"<tr><td>{E(x['panel'])}</td><td>{E(x['marker'])}</td><td>{E(x['count'])}</td>"
                   f"<td>{E(x['confidence'])}</td></tr>" for x in d.get("markers", []) if x["count"])
    return [f"    <p>{E(ystr_counts(d))}</p>",
            "    <table><tr><th>Panel</th><th>Marker</th><th>Repeats</th><th>Confidence</th></tr>" + rows + "</table>",
            small_note(YSTR_NOTE)]


BLOODTYPE_NOTE = ("Predicted from the DNA at the variants that define each antigen, not typed from blood. "
                  "RhD comes from the read depth over RHD; weak and partial D are not assessed. "
                  "The evidence for each site is in bloodtype/<sample>_bloodtype_evidence.tsv.")


CHIP_NOTE = ("Clonal haematopoiesis is common with age and most people with it never get a blood cancer: "
             "a call here is not a diagnosis. A call is a change in a CHIP gene's published classes, found by "
             "Mutect2 without a matched normal; one at a VAF of 35% or more is left out as possibly inherited. "
             "docs/49-chip.md says how to read it.")


def chip_lines(d, sample):
    """The CHIP section: the calls, the detection limit at this depth, the genes too shallow to read."""
    out = []
    for c in d.get("calls", []):
        out.append(f"{c['gene']} {c['protein_change']}: VAF {c['vaf']} ({c['alt_reads']} of {c['depth']} reads), {c['note']}")
    if not d.get("calls"):
        out.append("No CHIP call.")
    if d.get("median_limit") is not None:
        out.append(f"Detection limit: a clone of VAF about {d['median_limit'] * 100:.0f}% or more is found 8 times "
                   f"in 10 at this depth (median {d['median_depth']:.0f}x over the genes' coding exons). CHIP starts "
                   "at a VAF of 2%: smaller clones are likely missed, so no call does not rule one out.")
    else:
        out.append("No CHIP gene has the depth to find a clone: no result either way.")
    if d.get("not_covered"):
        out.append(f"Too shallow or not on the reference: {', '.join(d['not_covered'])}.")
    out.append(f"{d.get('other', 0)} other protein changes listed, not counted: chip/{sample}_chip.tsv.")
    out.append(CHIP_NOTE)
    return out


def bloodtype_body(d):
    """The blood group card: one row per system, then what the prediction is."""
    rows = "".join(f"<tr><td>{E(x['system'])}</td><td>{E(x['phenotype'])}</td><td>{E(x['genotype'])}</td>"
                   f"<td>{E(x['confidence'])}</td></tr>" for x in d.get("systems", []))
    return [stat("ABO", E(d.get("abo", "."))), stat("RhD", E(d.get("rhd", "."))),
            "    <table><tr><th>System</th><th>Predicted phenotype</th><th>Genotype</th><th>Confidence</th></tr>"
            + rows + "</table>", small_note(BLOODTYPE_NOTE)]


def lpa_lines(d):
    """Step 46's count with its interval, step 42's two LPA risk variants,
    the flag and the advice: the lines of the text section and the card."""
    if d.get("status") == "estimated":
        out = [f"KIV-2 repeats (both chromosomes): {d['kiv2_copies']} (95% interval {d['ci95_low']} to "
               f"{d['ci95_high']}), {d['class']}"]
    else:
        out = [f"KIV-2 repeats: {d.get('status', '.')}"]
    out.append(f"  {d.get('message', '.')}")
    for v in d.get("variants", []):
        n = v["risk_alleles"]
        out.append(f"{v['rsid']}: {v['genotype']} ({'not read' if n is None else str(n) + ' risk allele' + ('s' if n != 1 else '')})")
    if not d.get("variants"):
        out.append("rs10455872 and rs3798220: not read (step 42, the traits panel, did not run)")
    out.append(d.get("flag_text", "."))
    out.append(d.get("advice", "."))
    return out


VIRAL_NOTE = ("From the reads that do not map to the human genome, realigned to a fixed set of human viral "
              "genomes. Not a test for an infection: most adults carry herpesviruses that stay in their cells, and "
              "reagents carry viral DNA too. HHV-6 integration follows Liu et al. 2020 (depth over the viral genome "
              "at 0.30 of the human depth or more). Every genome: viral/<sample>_viral.tsv.")


def viral_lines(d):
    """Step 47's lines, the same in both reports: the HHV-6 call, then each genome with reads."""
    lines = [f"HHV-6: {d.get('hhv6', '.')}"]
    for x in d.get("found", []):
        copies = "" if x["copies_per_cell"] == "." else f", {x['copies_per_cell']} copies per cell"
        lines.append(f"{x['virus']}: {x['call']} ({x['reads']} reads, {float(x['breadth']):.0%} of the genome{copies})")
    if not d.get("found"):
        lines.append(f"No read on any of the {d.get('genomes', 0)} viral genomes.")
    return lines + [VIRAL_NOTE]


def viral_body(d):
    """The viral DNA card: the HHV-6 call, a row per genome with reads, the note."""
    rows = "".join(f"<tr><td>{E(x['virus'])}</td><td>{E(x['call'])}</td><td>{E(x['reads'])}</td>"
                   f"<td>{E(x['breadth'])}</td><td>{E(x['relative_depth'])}</td><td>{E(x['copies_per_cell'])}</td></tr>"
                   for x in d.get("found", []))
    table = ("    <table><tr><th>Virus</th><th>Call</th><th>Reads</th><th>Breadth</th><th>Depth / human</th>"
             "<th>Copies per cell</th></tr>" + rows + "</table>") if rows else small_note("No read on any viral genome.")
    return [stat("HHV-6", E(d.get("hhv6", "."))), table, small_note(VIRAL_NOTE)]


ROH_NOTE = ("Autosomes only: chrX and chrY are left out, and so is the part of a segment in a centromere or "
            "heterochromatin band (GRCh38 cytoband acen, gvar, stalk), where few markers map. No population allele "
            "frequencies were used (bcftools roh --AF-dflt 0.4): read the total and the short segments as rough; "
            "the long segments are the signal.")
CNVPYTOR_NOTE = ("Filtered: e-value below 0.01, under half of the reads with mapping quality 0 (q0) and under half "
                 "N bases (pN), no overlap with a centromere or heterochromatin band, and no chrY call when the "
                 "declared sex is female.")
MITO_POLYC_NOTE = ("Length heteroplasmy in the D-loop poly-C tracts is found in most people and is usually "
                   "technical; it is listed here and not counted as heteroplasmy.")
MITO_NOTE = ("ClinVar and gnomAD v3.1 (chrM sites) annotate the calls in step 20. MITOMAP (mitomap.org) has the "
             "curated disease status of mtDNA variants; it is a manual lookup, as it refuses scripted requests.")


def roh_lines(d):
    """Step 11's lines under the headline, the same in both reports."""
    lines = []
    if d.get("masked_mb"):
        lines.append(f"Left out in centromere and heterochromatin bands: {d['masked_mb']} MB")
    if d.get("sex_chromosome_segments"):
        lines.append(f"Segments on chrX/chrY, not counted: {d['sex_chromosome_segments']}")
    return lines + [ROH_NOTE]


def roh_region(r):
    return f"{r['region']}  {r['mb']}MB" + (f" ({r['masked_mb']} MB in a band left out)" if r.get("masked_mb") else "")


def telomere_lines(d):
    """The read count under the telomere value, and the warning when it is too thin."""
    if "intratelomeric_reads" not in d:
        return [f"Intratelomeric reads: not checked ({d['read_check']})"] if d.get("read_check") else []
    n = f"{d['intratelomeric_reads']:,}"
    if "read_floor" not in d:
        return [f"Intratelomeric reads: {n}; {d.get('read_check', 'not checked')}"]
    lines = [f"Intratelomeric reads: {n} (floor {d['read_floor']:,} at {d['depth']}x)"]
    if d.get("low_reads"):
        lines.append("WARNING: too few telomeric reads for a stable estimate; do not compare this value.")
    return lines


def cnvpytor_filtered_text(d):
    if "significant_filtered" not in d:
        return None
    dr = d.get("dropped") or {}
    parts = [f"{dr[k]} {label}" for k, label in (("q0", "q0 0.5 or more"), ("gap", "over N gaps (pN 0.5 or more)"),
                                               ("band", "over a centromere or heterochromatin band"),
                                               ("chry_female", "on chrY of a declared female")) if dr.get(k)]
    out = f"{d['significant_filtered']}" + (f" (left out: {', '.join(parts)})" if parts else "")
    if not d.get("declared_sex"):
        out += "; no declared sex, so chrY calls are kept"
    return out


def mito_variant_text(v):
    def val(x):
        return "." if x in (None, "", ".") else x
    cv = val(v.get("clinvar_significance"))
    cv = "not in ClinVar" if cv == "." else f"ClinVar {cv} ({val(v.get('clinvar_review'))}; {val(v.get('clinvar_condition'))})"
    gn = ("not in gnomAD" if val(v.get("gnomad_af_hom")) == "." else
          f"gnomAD AF hom {v['gnomad_af_hom']}, het {val(v.get('gnomad_af_het'))}")
    return f"chrM:{v['pos']} {v['ref']}>{v['alt']}  AF {val(v.get('af'))}  {cv}  {gn}" + ("  poly-C tract" if v.get("polyc") else "")


def mito_lines(d):
    """Step 20's lines under the counts, the same in both reports."""
    lines = []
    pc = d.get("polyc_heteroplasmic") or []
    if "polyc_heteroplasmic" in d:
        lines.append(f"Length heteroplasmy in the D-loop poly-C tracts ({', '.join(d.get('polyc_tracts') or [])}; "
                     f"common, usually technical): {len(pc)}")
        lines += [f"  {x['variant']}  AF {x['af']}" for x in pc]
    a = d.get("annotation") or {}
    if a.get("source"):
        lines.append(f"Annotation ({a['source']}): ClinVar {a.get('clinvar') or 'unknown'}; gnomAD {a.get('gnomad') or 'unknown'}")
        lines.append(f"ClinVar pathogenic or likely pathogenic: {a.get('clinvar_plp', 0)}")
        if a.get("variants"):
            lines.append("Heteroplasmic, or in ClinVar with a class other than benign:")
            lines += [f"  {mito_variant_text(v)}" for v in a["variants"]]
    elif a.get("note"):
        lines.append(f"ClinVar and gnomAD: not annotated ({a['note']})")
    return lines


# scripts/pgp-view.sh serves the genome viewer here (PGP_VIEW_PORT's default);
# the part after '#' is the locus it opens at.
VIEWER_URL = "http://127.0.0.1:8765/"


def viewer_links(sites):
    """One line under a findings table: each (chrom, pos, label) as a link that
    opens the reads at chrom:pos in the local genome viewer. Empty for none."""
    links = [f'<a href="{VIEWER_URL}#{E(str(c))}:{E(str(p))}" target="_blank" rel="noopener noreferrer">{E(str(t))}</a>'
             for c, p, t in sites if c and p]
    if not links:
        return []
    return ['    <p class="viewer">See the reads in the genome viewer on this computer (start it first: '
            "scripts/pgp-view.sh &lt;sample&gt;, docs/viewer.md): " + " &middot; ".join(links) + "</p>"]


def variant_site(variant):
    """(chrom, pos) of a "chrom:pos ref>alt" text, (None, None) for any other."""
    m = re.match(r"^([A-Za-z0-9._-]+):([0-9]+) ", variant or "")
    return (m.group(1), m.group(2)) if m else (None, None)


EXPLORE_NOTE = ("Nothing is uploaded by this pipeline. The docs link lists the files each service reads, "
                "which steps 12 and 37 write when they run, and what an upload shares.")


def explore_lines(d):
    """Where to explore a haplogroup further, the same in both reports: the
    links built from the branch name, then the note."""
    lines = [f"{x['name']}: {x['url']}" for x in d.get("explore", [])]
    return lines + [EXPLORE_NOTE] if lines else []


def explore_html(d):
    rows = [stat(x["name"], f'<a href="{E(x["url"])}">{E(x["url"])}</a>') for x in d.get("explore", [])]
    return rows + [small_note(EXPLORE_NOTE)] if rows else []


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


CONSUMER_NOTE = ("These files stay on this machine. Uploading one to a service shares your genotypes with it, "
                 "and with the relatives it matches you to: docs/40-export-consumer.md says what to weigh first.")


def consumer_export_lines(d):
    """Step 40's lines, the same in both reports."""
    sites, called = d.get("sites", 0), d.get("called", 0)
    pct = f" ({100.0 * called / sites:.1f}%)" if sites else ""
    lines = [f"Sites with a genotype: {called} of {sites}{pct}; '--' at {d.get('no_call', 0)}",
             f"Homozygous reference from the gVCF's reference blocks: {d.get('from_reference_blocks', 0)}",
             f"Files: consumer_export/{d.get('file_23andme', '.')} and consumer_export/{d.get('file_ancestrydna', '.')}"]
    reasons = d.get("no_call_reasons") or {}
    if reasons:
        lines.append("'--' because: " + "; ".join(f"{k} {v}" for k, v in sorted(reasons.items(), key=lambda x: -x[1])))
    return lines + [CONSUMER_NOTE]


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
    rdc = S.get("read_count", MISSING)
    if rdc["state"] in ("ok", "stale"):
        qc.insert(1, stat("Reads kept by alignment",
                          badge(read_count_text(rdc["data"]),
                                "green" if rdc["data"].get("read_count_check") == "ok" else "red")))
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
               for sec in (cov, sexc, sqc, rdc) if stale_note(sec)]
    a(['  <div class="card">', "    <h2>Quality Control</h2>"] + qc_note + qc + ["  </div>"])

    mo = S.get("mosaic_sex_loss", MISSING)
    q = mo["data"]
    label, colour = MOSAIC_CALL.get(q.get("call"), (q.get("call") or "?", "gray"))
    mo_body = [stat(q.get("chrom") or "chrY or chrX", badge(label, colour)),
               stat("Why", f'<span style="font-weight:normal;font-size:13px">{E(q.get("reason") or "")}</span>')]
    mo_body += [small_note(line) for line in mosaic_lines(q)] + [small_note(MOSAIC_NOTE)]
    a(card(mo, "Mosaic Loss of Chromosome Y or X", mo_body))

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
    if "cdna_flagged" in c:
        cv_body.insert(2, stat("Flagged by the cDNA read check",
                               badge(c["cdna_flagged"], "red" if c["cdna_flagged"] else "gray")))
    if "low_confidence" in c:
        f = c.get("call_floor") or {}
        cv_body.insert(1, stat(f"Low-confidence calls (GQ below {f.get('GQ')} or DP below {f.get('DP')}), not counted",
                               badge(c["low_confidence"], "yellow" if c["low_confidence"] else "gray")))
    if clinvar_freq_line(c):
        cv_body.insert(-1, small_note(clinvar_freq_line(c)))
    a(card(S["clinvar"], "ClinVar Screening", cv_body))
    if S.get("couple", MISSING)["state"] != "missing":
        a(card(S["couple"], "Couple Carrier Screen",
               [f"    <p>{E(line)}</p>" for line in couple_lines(S["couple"]["data"])], full=True))

    pg = S["pharmcat"]["data"]
    cp = S["cpic"]["data"]
    pgx_body = [stat("PharmCAT report", badge("Complete", "green") if S["pharmcat"]["state"] in ("ok", "stale")
                     else badge("Not run", "gray"))]
    if pg.get("version"):
        pgx_body.append(stat("PharmCAT version", E(pg["version"])))
    if S["pypgx"]["state"] in ("ok", "stale") and "comparison" in S["pypgx"]["data"]:
        pgx_body.append(stat("PharmCAT vs pypgx conflicts", S["pypgx"]["data"]["comparison"]["conflicts"]))
    a(card(S["pharmcat"], "Pharmacogenomics", pgx_body))
    if S["pypgx"]["state"] != "missing":
        pd = S["pypgx"]["data"]
        a(card(S["pypgx"], "pypgx Pharmacogenomics",
               [stat("Genes called", f"{pd.get('genes_called')}/{pd.get('genes_total')}"),
                stat("CYP2D6 diplotype", E(pd.get("cyp2d6") or "not called"))]
               + [f"    <p>{E(line)}</p>" for line in pypgx_lines(pd)]))
    cn_rows = pgx_copy_number(S)
    if cn_rows:
        a(card(None, "PGx Copy Number", [stat(g, E(v)) for g, v in cn_rows]))

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
    prob = s["cyp2d6"].get("probable")
    if prob:
        cy_body.append(stat(prob["label"].capitalize(), E(probable_text(prob))))
    a(["  <div class=\"card\">", "    <h2>CYP2D6 Across Callers</h2>"]
      + [f'    <div class="stale">{E(S[k]["title"] + ": " + stale_note(S[k]))}</div>'
         for k in ("cpic", "pypgx", "cyrius") if S[k]["state"] in ("stale", "failed")]
      + cy_body + ["  </div>"])

    hla = S["hla"]["data"]
    hla_body = [stat(l["gene"], E(hla_text(l))) for l in hla.get("loci", [])]
    if hla.get("derived"):
        hla_body += [small_note(x) for x in hla_derived_lines(hla["derived"])] + [small_note(HLA_DERIVED_NOTE)]
    hla_body += [stat("HLA database", E(s["databases"]["hla_database"] or "unknown")), small_note(HLA_NOTE)]
    if hla.get("optitype") is not None:
        hla_body += [stat(f"{g}: T1K vs OptiType", E(line)) if g else stat("OptiType check", E(line))
                     for g, line in hla_optitype_lines(hla["optitype"])]
    if hla.get("kir") is not None:
        hla_body += [stat(g, E(line)) if g else small_note(line) for g, line in kir_lines(hla["kir"])]
        hla_body.append(small_note(KIR_NOTE))
    a(card(S["hla"], "HLA Typing", hla_body))

    m, dl, cn, sv = (S[k]["data"] for k in ("manta", "delly", "cnvpytor", "sv_consensus"))
    sv_body = [stat("Manta SVs (total / PASS)", f"{m.get('total', 'N/A')} / {m.get('pass', 'N/A')}"),
               stat("Delly SVs (total / PASS)", f"{dl.get('total', 'N/A')} / {dl.get('pass', 'N/A')}"),
               stat("CNVpytor CNVs", cn.get("total", "N/A")),
               stat("CNVpytor significant (e-val < 0.01)", cn.get("significant", "N/A"))]
    if cnvpytor_filtered_text(cn):
        sv_body += [stat("CNVpytor significant after filters", E(cnvpytor_filtered_text(cn))), small_note(CNVPYTOR_NOTE)]
    sv_body.append(stat("Consensus SVs (2+ callers)", sv.get("consensus", "N/A")))
    notes = [f'    <div class="stale">{E(S[k]["title"] + ": " + stale_note(S[k]))}</div>'
             for k in ("manta", "delly", "cnvpytor", "sv_consensus") if stale_note(S[k])]
    a(['  <div class="card">', "    <h2>Structural Variants</h2>"] + notes + sv_body + ["  </div>"])
    if S.get("annotsv", MISSING)["state"] != "missing":
        av = S["annotsv"]["data"]
        body = [stat("Status", done_badge(S["annotsv"]))] + ([small_note(x) for x in annotsv_lines(av)] if av else [])
        if av.get("rows"):
            body.append("    <table><tr><th>SV</th><th>Class</th><th>Genes</th><th>Callers</th><th>Checks</th></tr>"
                        + "".join(f"<tr><td>{E(x['type'] + ' ' + x['sv'])}</td><td>{E(x['acmg_class'])}</td>"
                                  f"<td>{E(annotsv_genes_text(x))}</td><td>{E(annotsv_support_text(x))}</td>"
                                  f"<td>{E('; '.join(annotsv_flags(x)) or 'none')}</td></tr>" for x in av["rows"])
                        + "</table>")
            body.append(small_note(ANNOTSV_NOTE))
        a(card(S["annotsv"], "Structural Variants, Class 4 and 5 (AnnotSV)", body, full=True))

    cs = S["cpsr"]["data"]
    cpsr_body = [stat("CPSR report", badge("Complete", "green") if cs.get("html_report") else badge("No HTML", "gray"))]
    for k, val in (cs.get("classification") or {}).items():
        cpsr_body.append(stat(k, val))
    cpsr_body += [small_note(line.strip()) for line in cpsr_lines(cs, s["databases"])]
    a(card(S["cpsr"], "Cancer Predisposition", cpsr_body))

    eh = S["expansions"]["data"]
    eh_body = [stat("Status", done_badge(S["expansions"]))]
    if eh.get("stranger_outdated"):
        eh_body.append(stat("Outside the normal range (Stranger)",
                            badge("older than ExpansionHunter's result: rerun step 9b", "yellow")))
    elif not eh.get("stranger"):
        eh_body.append(stat("Outside the normal range (Stranger)", badge("not run", "gray")))
    else:
        fl = eh.get("flagged") or []
        eh_body.append(stat("Outside the normal range (Stranger)", badge(len(fl), "yellow" if fl else "gray")))
    if eh.get("unsized"):
        eh_body.append(stat("Not sized from short reads (not counted as normal)", badge(len(eh["unsized"]), "yellow")))
    if eh.get("intermediate"):
        eh_body.append(stat("In an intermediate range", badge(", ".join(eh["intermediate"]), "yellow")))
    for l in eh.get("flagged") or []:
        if l.get("rfc1_motif"):
            eh_body.append(small_note(rfc1_text(l["rfc1_motif"])))
            eh_body.append(small_note(RFC1_NOTE))
    # Every locus, the five key loci first.
    rows = str_rows(eh)
    if rows:
        eh_body.append("    <table><tr>" + "".join(f"<th>{c}</th>" for c in ("Locus", "Repeat Count", "Stranger", "Notes"))
                       + "</tr>"
                       + "".join("<tr>" + "".join(f"<td>{E(str(x))}</td>" for x in r) + "</tr>" for r in rows)
                       + "</table>")
    eh_body.append(small_note(EH_NOTE))
    a(card(S["expansions"], "Repeat Expansions", eh_body))

    r, hg, tl = S["roh"]["data"], S["haplogroup"]["data"], S["telomere"]["data"]
    a(card(S["roh"], "Runs of Homozygosity",
           [stat("Status", done_badge(S["roh"])),
            stat("ROH total", f"{r.get('total_mb', 0):.1f} MB"),
            stat("ROH largest segment", f"{r.get('largest_mb', 0):.1f} MB"),
            stat("Segments", r.get("segments", 0)),
            stat("Autosomal ROH > 5 MB", len(r.get("autosomal_over_5mb", [])))]
           + [small_note(roh_region(x)) for x in r.get("autosomal_over_5mb", [])]
           + [small_note(x) for x in roh_lines(r)]))
    a(card(S["haplogroup"], "Mitochondrial Haplogroup",
           [stat("Status", done_badge(S["haplogroup"])), stat("Haplogroup", E(hg.get("haplogroup", "."))),
            stat("Contamination (haplocheck)", E(contamination_text(hg)))] + explore_html(hg)))
    if S.get("paralogs", MISSING)["state"] != "missing":
        pg = S["paralogs"]["data"]
        rows = "".join(f"<tr><td>{E(g['gene'])}</td><td>{E(g['copy_number'])}</td><td>{E(g['status'])}</td>"
                       f"<td>{E(g['detail'])}</td></tr>" for g in pg.get("genes", []))
        body = [stat("Status", done_badge(S["paralogs"]))]
        if rows:
            body.append("    <table><tr><th>Gene</th><th>Copy number</th><th>Call</th><th>Detail</th></tr>"
                        + rows + "</table>")
        smn = smn1_text(pg.get("smn1") or {})
        if smn:
            body.append(stat("SMN1 by both tools", E(smn)))
        if pg.get("other_regions"):
            body.append(stat("Other Paraphase regions", f"{pg['other_regions']} ({pg['other_regions_called']} with a call)"))
        body.append(stat("Variants on phased haplotypes", sum(pg.get("variants", {}).values())))
        body.append(small_note(PARALOG_NOTE))
        a(card(S["paralogs"], "Paralog Genes (Paraphase)", body, full=True))
    if S.get("parascopy", MISSING)["state"] != "missing":
        a(card(S["parascopy"], "SMN1/SMN2 Copy Number (Parascopy)",
               [stat("Status", done_badge(S["parascopy"]))]
               + [stat(*(E(x) for x in line.split(": ", 1))) for line in parascopy_lines(S["parascopy"]["data"])]
               + [small_note(PARASCOPY_NOTE)]))
    unchecked = paralogs_unchecked(S)
    if unchecked:
        a(card(None, "Paralog Copy Number", [stat("Status", badge("Not checked", "gray")), small_note(unchecked)]))
    if S.get("y_haplogroup", MISSING)["state"] != "missing":
        yh = S["y_haplogroup"]["data"]
        a(card(S["y_haplogroup"], "Y-Chromosome Haplogroup",
               [stat("Haplogroup", E(yh.get("haplogroup", "."))), stat("Markers", E(str(yh.get("valid_markers", ".")))),
                stat("QC-score", E(str(yh.get("qc_score", "."))))] + explore_html(yh)))
    if S.get("y_str", MISSING)["state"] != "missing":
        ys = S["y_str"]
        a(card(ys, "Y-STR Markers (HipSTR)", ystr_body(ys["data"]) if ys["state"] in ("ok", "stale") else [],
               full=True))
    if S.get("consumer_export", MISSING)["state"] != "missing":
        a(card(S["consumer_export"], "Raw-Data Export (23andMe v5 layout)",
               [f"    <p>{E(line)}</p>" for line in consumer_export_lines(S["consumer_export"]["data"])]))
    if S.get("bloodtype", MISSING)["state"] != "missing":
        a(card(S["bloodtype"], "Blood Groups (predicted)", bloodtype_body(S["bloodtype"]["data"]), full=True))
    if S.get("lpa", MISSING)["state"] != "missing":
        lp = S["lpa"]["data"]
        a(card(S["lpa"], "Lp(a): LPA KIV-2 Repeat Count",
               ([stat("Genetic sign of high Lp(a)", badge("Yes", "yellow"))] if lp.get("flag") == "high" else [])
               + [f"    <p>{E(line)}</p>" for line in lpa_lines(lp)]))
    if S.get("viral", MISSING)["state"] != "missing":
        a(card(S["viral"], "Viral DNA in the Unmapped Reads", viral_body(S["viral"]["data"]), full=True))
    if S.get("chip", MISSING)["state"] != "missing":
        ch = S["chip"]["data"]   # the card shows its body only for a section that was read
        a(card(S["chip"], "Clonal Haematopoiesis (CHIP)",
               [stat("CHIP calls", E(str(len(ch.get("calls", [])))))]
               + [f"    <p>{E(line)}</p>" for line in chip_lines(ch, s["sample"])[:-1]]
               + [small_note(CHIP_NOTE)], full=True))
    a(card(S["telomere"], "Telomere content (relative)",
           [stat("Telomere content", E(str(tl.get("tel_content", "."))))]
           + [stat("Warning", badge(x[len("WARNING: "):], "yellow")) if x.startswith("WARNING: ") else small_note(x)
              for x in telomere_lines(tl)]
           + [small_note(TELOMERE_NOTE)]))

    mi = S["mito"]["data"]
    ma = mi.get("annotation") or {}
    mito_body = [stat("chrM variants (PASS)", mi.get("pass")),
                 stat("Heteroplasmic (AF 0.05 to 0.95, outside the D-loop poly-C tracts)", mi.get("heteroplasmic")),
                 stat("Length heteroplasmy in the D-loop poly-C tracts", len(mi.get("polyc_heteroplasmic") or []))]
    mito_body += [small_note(f"{x['variant']}  AF {x['af']}") for x in mi.get("polyc_heteroplasmic") or []]
    if mi.get("polyc_heteroplasmic"):
        mito_body.append(small_note(MITO_POLYC_NOTE))
    if ma.get("source"):
        plp = ma.get("clinvar_plp", 0)
        mito_body.append(stat("ClinVar pathogenic or likely pathogenic", badge(plp, "yellow" if plp else "gray")))
        if ma.get("variants"):
            mito_body.append("    <table><tr>" + "".join(f"<th>{c}</th>" for c in (
                "Variant", "AF", "ClinVar", "Review", "Condition", "gnomAD AF hom", "gnomAD AF het", "Poly-C tract"))
                + "</tr>" + "".join(
                    "<tr>" + "".join(f"<td>{E(str(x))}</td>" for x in (
                        f"chrM:{v['pos']} {v['ref']}>{v['alt']}", v.get("af"), v.get("clinvar_significance"),
                        v.get("clinvar_review"), v.get("clinvar_condition"), v.get("gnomad_af_hom"),
                        v.get("gnomad_af_het"), "yes" if v.get("polyc") else "")) + "</tr>"
                    for v in ma["variants"]) + "</table>")
            mito_body += viewer_links(("chrM", v["pos"], f"chrM:{v['pos']}") for v in ma["variants"])
        mito_body.append(small_note(f"ClinVar {ma.get('clinvar') or 'unknown'}; gnomAD {ma.get('gnomad') or 'unknown'}"))
    elif ma.get("note"):
        mito_body.append(small_note(f"ClinVar and gnomAD: not annotated ({ma['note']})"))
    mito_body.append(small_note(MITO_NOTE))
    a(card(S["mito"], "Mitochondrial Analysis", mito_body, full=bool(ma.get("variants"))))

    cl = S["clinical"]["data"]
    a(card(S["clinical"], "Clinical Variant Filter",
           [stat("Total interesting variants", cl.get("variants"))]
           + ([stat("Genes", cl["genes"])] if "genes" in cl else [])
           + [stat(f"{k} impact", v) for k, v in sorted((cl.get("by_impact") or {}).items())]
           + ([stat("Frequency unknown (not rare)", cl["frequency_unknown"])] if "frequency_unknown" in cl else [])
           + ([stat("Weak genotype calls", cl["weak_calls"])] if "weak_calls" in cl else [])))
    sf = cl.get("acmg_sf")
    if sf and (sf["clinvar_hits"] or sf["high_impact"]):
        rows = "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(x['variant'])}</td><td>{E(x['genotype'])}</td>"
                       f"<td>{E('ClinVar ' + x['significance'] + ('; low-confidence call: check the reads' if x.get('low_confidence') else ''))}</td></tr>\n"
                       for x in sf["clinvar_hits"])
        rows += "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(x['variant'])}</td><td>{E(x['genotype'])}</td>"
                        f"<td>{E(x['consequence'] + high_marks(x, '; '))}</td></tr>\n" for x in sf["high_impact"])
        a(card(S["clinical"], f"Secondary-Findings Genes ({sf['version']})",
               ["    <table>", "      <tr><th>Gene</th><th>Variant</th><th>Genotype</th><th>Evidence</th></tr>",
                rows.rstrip("\n"), "    </table>"]
               + viewer_links((*variant_site(x["variant"]), f"{x['gene']} {x['variant'].split(' ')[0]}")
                              for x in sf["clinvar_hits"] + sf["high_impact"])
               + [f"    <p>{E(ACMG_NOTE)}</p>"], full=True))

    sl = S["slivar"]["data"]
    sl_body = [stat("Prioritized variants", sl.get("prioritized"))]
    if "frequency_unknown" in sl:
        sl_body.append(stat("Frequency unknown (not rare)", sl["frequency_unknown"]))
    if "compound_het_variants" in sl:
        sl_body.append(stat("Compound het candidates", f"{sl['compound_het_variants']} variants, {sl['compound_het_genes']} genes"))
    if "compound_het_phase" in sl:
        sl_body += comphet_phase_html(sl)
    a(card(S["slivar"], "Variant Prioritization (slivar)", sl_body))

    if S["prs"]["state"] in ("ok", "stale"):
        d = S["prs"]["data"]
        rows = "".join(f"<tr><td>{E(x['condition'])}</td><td>{E(x['pgs_id'])}</td><td>{E(x['score'])}</td>"
                       f"<td>{E(x['matched'])}/{E(x['total'])}</td>"
                       f"<td>{E(prs_percentile(x) or 'raw score only')}</td></tr>\n"
                       for x in d["scores"])
        a(card(S["prs"], "Polygenic Risk Scores" + ("" if d.get("adjusted") else " (raw, not percentiles)"),
               ["    <table>", "      <tr><th>Condition</th><th>PGS</th><th>Score</th><th>Variants matched</th><th>Percentile</th></tr>",
                rows.rstrip("\n"), "    </table>", f"    <p>{E(prs_note(d))}</p>"], full=True))

    if S.get("ancestry_fine", MISSING)["state"] != "missing":
        d = S["ancestry_fine"]["data"]
        body = []
        for level, title in ANCESTRY_FINE_LEVELS:
            rows = "".join(f"<tr><td>{E(r['group'])}</td><td>{E(share(r['fraction']))}</td>"
                           f"<td>{E(ancestry_fine_spread(r).strip(' ()') or '.')}</td></tr>\n"
                           for r in ancestry_fine_rows(d, level))
            if rows:
                body += [f"    <h3>{E(title)}</h3>", "    <table>",
                         "      <tr><th>Group</th><th>Fraction</th><th>Bootstrap range</th></tr>",
                         rows.rstrip("\n"), "    </table>"]
        if ancestry_fine_left_out(d):
            body.append(f"    <p>{E(ancestry_fine_left_out(d))}</p>")
        body.append(stat("Painted genome", E(ancestry_fine_painted(d))))
        if d.get("painting_png"):
            body.append(f'    <img alt="Chromosome painting: chromosomes 1 to 22, two haplotypes each, coloured by '
                        f'continental group" style="max-width:100%" src="data:image/png;base64,{d["painting_png"]}">')
            body.append("    <p>" + " ".join(f'<span style="color:{E(x["colour"])}">&#9632;</span> {E(x["group"])}'
                                             for x in d.get("legend", [])) + " <span style=\"color:#e1e1e1\">&#9632;</span> not painted</p>")
        body.append(f"    <p>{E(ANCESTRY_FINE_NOTE)}</p>")
        a(card(S["ancestry_fine"], "Fine-Grained Ancestry", body, full=True))

    tr_sec = S.get("traits", MISSING)
    a(card(tr_sec, "Traits", traits_html(tr_sec["data"], s["sample"]) if tr_sec["state"] in ("ok", "stale") else [],
           full=True))

    ex_sec = S.get("exomiser", MISSING)
    if ex_sec["state"] != "missing":
        a(card(ex_sec, "Symptom-Driven Ranking (Exomiser)",
               exomiser_html(ex_sec["data"], s["sample"]) if ex_sec["state"] in ("ok", "stale") else [], full=True))

    if S["cpic"]["state"] in ("ok", "stale") and (cp["non_normal"] or cp.get("ambiguous")):
        rows = "".join(f"<tr><td>{E(g['gene'])}</td><td>{E(g['diplotype'])}</td><td>{E(g['phenotype'])}</td></tr>\n"
                       for g in cp["genes"] if g["status"] in ("non-normal", "ambiguous"))
        body = ["    <table>", "      <tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th></tr>", rows.rstrip("\n"),
                "    </table>"]
        body += [f'    <div class="stale">WARNING: {E(x)}</div>' for x in cp["warnings"]]
        body.append(f"    <p>Drugs for each gene: cpic/{E(s['sample'])}_cpic_recommendations.txt</p>")
        a(card(S["cpic"], "Genes With a Non-Normal or Unresolved Phenotype (CPIC)", body, full=True))

    if S["clinvar"]["state"] in ("ok", "stale") and c.get("hits"):
        # The combined-allele cell is empty, not '.', for a hit that is
        # neither a combined allele nor part of one: no note is no missing data.
        # The cDNA read check column only when step 06 wrote its table.
        reads = "cdna_flagged" in c
        # The frequency and call columns only for a summary that carries them.
        freq = "frequency" in c
        call = "low_confidence" in c
        nh = len(c["hits"])
        rows = [f"<tr><td>{E(h['chrom'])}</td><td>{h['pos']}</td><td>{E(h['ref'])}</td><td>{E(h['alt'])}</td>"
                f"<td>{E(h['genotype'])}</td><td>{E(h['gene'])}</td><td>{E(h['significance'])}</td>"
                f"<td>{E(h['review_status'])}</td><td>{h['stars']}</td>"
                f"<td>{E(combined_cell(h))}</td>"
                + (f"<td>{E(af_text(h))}</td>" if freq else "")
                + (f"<td>{E(call_text(h))}</td>" if call else "")
                + (f"<td>{E(cdna_text(h))}</td>" if reads else "") + "</tr>" for h in c["hits"][:50]]
        a(card(S["clinvar"], f"ClinVar Hits (best-reviewed first, {min(nh, 50)} of {nh})",
               ["    <table>",
                "      <tr><th>Chr</th><th>Position</th><th>Ref</th><th>Alt</th><th>Genotype</th><th>Gene</th>"
                "<th>Significance</th><th>Review status</th><th>Stars</th><th>Combined allele</th>"
                + ("<th>gnomAD frequency</th>" if freq else "")
                + ("<th>Call quality</th>" if call else "")
                + ("<th>cDNA read check</th>" if reads else "") + "</tr>"]
               + ["      " + x for x in rows] + ["    </table>"]
               + viewer_links((h["chrom"], h["pos"], f"{h['gene']} {h['chrom']}:{h['pos']}") for h in c["hits"][:50]),
               full=True))

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
