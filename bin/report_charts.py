#!/usr/bin/env python3
"""report_charts.py: the report's five charts, as data and as inline SVG.

  build(summary, notes, helpers) -> [spec, ...]
  svg(spec) -> str

One chart model for both places a chart is shown. bin/render_report.py puts
each SVG in its report section, under that section's caveat, and writes every
spec with its SVG to <sample>_charts.json (--charts-json), which the live
page's charts plugin reads (bin/pgp_live_charts.py, assets/live/plugins/charts.js).

The five charts, each drawn only when its section ran:
  clinvar   ClinVar hits per review star, 4 to 0, and a grey bar for the
            low-confidence calls, which the star counts leave out. Yellow and
            grey only, never green; with no hit a grey frame says that no hit
            is not an all-clear.
  ancestry  fine-grained ancestry: each group render_report lists, its share
            as a dot and its bootstrap range as a line.
  prs       PRS percentiles on a 0 to 100 strip, only when pgsc_calc had an
            ancestry panel. Without one a raw score has no scale, so the chart
            is the report's own note; a score of the other sex, or one with no
            percentile, says "no percentile" and is never placed.
  genome    the chromosomes 1 to 22, X and Y of GRCh38 with the runs of
            homozygosity of 5 MB or more and a tick for each ClinVar hit.
            Hits on chrM and on other contigs are counted in not_drawn.
  mito      the heteroplasmic chrM sites by position and allele fraction, the
            D-loop poly-C length variants marked apart.

A spec is {id, section, kind, title, caveat_key, caveat, axis, marks,
not_drawn}: section is the collect_summary.SECTIONS key whose data it shows;
each mark is {id, label, value, low, high, chrom, pos, end, tone, anchor,
locus, detail}, with None for what does not apply; a ClinVar tick on the
genome chart also carries the ClinVar caveat, the chart's own being ROH's. anchor is the id of the
report element the mark links to (the ClinVar hits table, or the section); locus
is "chrom:pos" or "chrom:start-end" for the genome viewer. The caveat text
is passed in by the caller (notes), never written here, so the chart and the
report say the same thing.

The SVG has no script, no style attribute and no outside reference: colours
are presentation attributes, links are <a href="#anchor"> inside the page and
each mark carries a <title> and data-mark="<id>". Text uses currentColor, so
it follows the page's colour on screen and in print. Standard library only.
"""
import html
import math
import re

SCHEMA_VERSION = 1

# GRCh38 primary assembly lengths (UCSC hg38 chrom.sizes).
GRCH38 = (("chr1", 248956422), ("chr2", 242193529), ("chr3", 198295559), ("chr4", 190214555),
          ("chr5", 181538259), ("chr6", 170805979), ("chr7", 159345973), ("chr8", 145138636),
          ("chr9", 138394717), ("chr10", 133797422), ("chr11", 135086622), ("chr12", 133275309),
          ("chr13", 114364328), ("chr14", 107043718), ("chr15", 101991189), ("chr16", 90338345),
          ("chr17", 83257441), ("chr18", 80373285), ("chr19", 58617616), ("chr20", 64444167),
          ("chr21", 46709983), ("chr22", 50818468), ("chrX", 156040895), ("chrY", 57227415))
CHROM_LEN = dict(GRCH38)
MITO_LEN = 16569
ROH_MIN_MB = 5.0       # collect_summary.ROH_MIN_BP: the runs F_ROH counts
HITS_SHOWN = 50        # render_report shows the first 50 ClinVar hits in its table
HITS_ANCHOR = "clinvar-hits"   # the id render_report puts on that table

# Colours: yellow for something to look at, grey for what is not counted,
# blue for plain information. No green anywhere: a chart is not an all-clear.
TONE = {"attention": "#e0b94a", "muted": "#8a909b", "info": "#5b9ee6"}
FRAME = "#8a909b"
LINE = "#4a505b"
W = 640

E = html.escape


def _ok(S, key):
    return S.get(key, {}).get("state") in ("ok", "stale")


def mark(mid, label, value=None, **k):
    m = {"id": mid, "label": label, "value": value, "low": None, "high": None, "chrom": None, "pos": None,
         "end": None, "tone": "info", "anchor": None, "locus": None, "detail": ""}
    m.update(k)
    return m


def spec(cid, section, kind, title, notes, caveat_key, marks, axis=None, not_drawn=None, text=None):
    out = {"id": cid, "section": section, "kind": kind, "title": title, "caveat_key": caveat_key,
           "caveat": notes[caveat_key], "axis": axis or {}, "marks": marks, "not_drawn": not_drawn or []}
    if text is not None:
        out["text"] = text
    return out


def hit_anchor(i):
    """The report's hits table for a hit it lists, the ClinVar section for one past it."""
    return HITS_ANCHOR if i < HITS_SHOWN else "clinvar"


def hit_label(h):
    return f"{h.get('gene') or '.'} {h['chrom']}:{h['pos']} {h.get('ref', '')}>{h.get('alt', '')}"


# --- the five charts ------------------------------------------------------------------------

def clinvar_chart(S, notes):
    if not _ok(S, "clinvar"):
        return None
    c = S["clinvar"]["data"]
    by = c.get("by_stars") or {}
    hits = c.get("hits") or []
    marks = []
    for k in ("4", "3", "2", "1", "0"):
        n = int(by.get(k, 0))
        idx = [i for i, h in enumerate(hits) if str(h.get("stars")) == k and not h.get("low_confidence")]
        names = ", ".join(hit_label(hits[i]) for i in idx[:10]) + (" and more" if len(idx) > 10 else "")
        marks.append(mark(f"stars-{k}", f"{k} star{'' if k == '1' else 's'}", n,
                          tone="attention" if n else "muted",
                          anchor=hit_anchor(idx[0]) if idx else "clinvar",
                          detail=f"{n} hit{'' if n == 1 else 's'} with {k} review star{'' if k == '1' else 's'}"
                          + (f": {names}" if names else "")))
    if "low_confidence" in c:
        n = int(c["low_confidence"])
        idx = [i for i, h in enumerate(hits) if h.get("low_confidence")]
        f = c.get("call_floor") or {}
        marks.append(mark("low-confidence", "Low-confidence, not counted", n, tone="muted",
                          anchor=hit_anchor(idx[0]) if idx else "clinvar",
                          detail=f"{n} call{'' if n == 1 else 's'} with GQ below {f.get('GQ')} or DP below "
                          f"{f.get('DP')}: listed in the hits, not counted in the stars"))
    return spec("clinvar-stars", "clinvar", "bars", "ClinVar hits by review stars", notes, "CLINVAR_FREQ_NOTE",
                marks, axis={"unit": "hits"})


def ancestry_chart(S, notes, rows_of, levels):
    if not _ok(S, "ancestry_fine"):
        return None
    d = S["ancestry_fine"]["data"]
    marks = []
    for level, title in levels:
        for r in rows_of(d, level):
            v = 100 * r["fraction"]
            low = 100 * r["low"] if r.get("low") is not None else None
            high = 100 * r["high"] if r.get("high") is not None else None
            marks.append(mark(f"{level}-{r.get('code') or len(marks)}", r["group"], v, low=low, high=high,
                              anchor="ancestry-fine", group=title,
                              detail=f"{title}: {v:.0f}%"
                              + (f" (bootstrap range {low:.0f}% to {high:.0f}%)" if low is not None else "")))
    if not marks:
        return None
    return spec("ancestry-fine", "ancestry_fine", "ranges", "Fine-grained ancestry, with bootstrap ranges", notes,
                "ANCESTRY_FINE_NOTE", marks, axis={"min": 0, "max": 100, "unit": "%"})


def prs_chart(S, notes, percentile_text, note_text, not_applicable):
    if not _ok(S, "prs"):
        return None
    d = S["prs"]["data"]
    if not d.get("adjusted"):
        return spec("prs", "prs", "text", "Polygenic risk scores", notes, "PRS_ANCESTRY_NOTE", [],
                    text=note_text(d))
    marks = []
    for i, r in enumerate(d.get("scores") or []):
        value = None
        if r.get("sex_check") != not_applicable and r.get("percentile"):
            try:
                value = float(r["percentile"])
            except ValueError:
                value = None
            if value is not None and not 0 <= value <= 100:
                value = None
        marks.append(mark(f"pgs-{r.get('pgs_id') or i}", f"{r.get('condition') or '.'} ({r.get('pgs_id') or '.'})",
                          value, tone="info" if value is not None else "muted", anchor="prs",
                          detail=(percentile_text(r) or "no percentile") if value is not None
                          else "no percentile" + (f": {percentile_text(r)}" if percentile_text(r) else "")))
    return spec("prs", "prs", "strip", "Polygenic risk score percentiles", notes, "PRS_ANCESTRY_NOTE", marks,
                axis={"min": 0, "max": 100, "unit": "percentile"})


REGION_RE = re.compile(r"^(chr[0-9XYM]+):([0-9]+)-([0-9]+)$")


def genome_chart(S, notes):
    roh, cv = _ok(S, "roh"), _ok(S, "clinvar")
    if not (roh or cv):
        return None
    marks, not_drawn = [], {}
    if roh:
        short = 0
        for i, r in enumerate(S["roh"]["data"].get("autosomal_over_5mb") or []):
            m = REGION_RE.match(str(r.get("region", "")))
            if not m:
                raise ValueError(f"run of homozygosity {r.get('region')!r}: not a chrom:start-end region")
            chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            if end < start:
                raise ValueError(f"run of homozygosity {r['region']!r}: the end is before the start")
            if float(r.get("mb", 0)) < ROH_MIN_MB:
                short += 1
                continue
            if chrom not in CHROM_LEN:
                not_drawn[f"runs of homozygosity on {chrom}"] = not_drawn.get(f"runs of homozygosity on {chrom}", 0) + 1
                continue
            marks.append(mark(f"roh-{i}", f"Run of homozygosity {r['region']}", float(r["mb"]), chrom=chrom,
                              pos=start, end=end, tone="info", anchor="roh", locus=f"{chrom}:{start}-{end}",
                              kind="roh", detail=f"{r['mb']} MB run of homozygosity at {r['region']}"))
        if short:
            not_drawn["runs of homozygosity under 5 MB"] = short
    if cv:
        for i, h in enumerate(S["clinvar"]["data"].get("hits") or []):
            chrom = h.get("chrom")
            if chrom == "chrM":
                not_drawn["ClinVar hits on chrM"] = not_drawn.get("ClinVar hits on chrM", 0) + 1
                continue
            if chrom not in CHROM_LEN:
                not_drawn["ClinVar hits on other contigs (alt, unplaced)"] = \
                    not_drawn.get("ClinVar hits on other contigs (alt, unplaced)", 0) + 1
                continue
            low = bool(h.get("low_confidence"))
            marks.append(mark(f"hit-{i}", hit_label(h), h.get("stars"), chrom=chrom, pos=int(h["pos"]),
                              tone="muted" if low else "attention", anchor=hit_anchor(i),
                              locus=f"{chrom}:{h['pos']}", kind="hit", caveat=notes["CLINVAR_FREQ_NOTE"],
                              detail=f"{h.get('significance') or '.'}, {h.get('stars')} star"
                              f"{'' if h.get('stars') == 1 else 's'}, {h.get('genotype') or '.'}"
                              + (", low-confidence call" if low else "")))
    return spec("genome", "roh", "ideogram", "Runs of homozygosity and ClinVar hits on the chromosomes", notes,
                "ROH_NOTE", marks, axis={"assembly": "GRCh38"},
                not_drawn=[{"what": k, "count": v} for k, v in not_drawn.items()])


POLYC_RE = re.compile(r"^chrM:([0-9]+) ")


def mito_chart(S, notes):
    if not _ok(S, "mito"):
        return None
    d = S["mito"]["data"]
    marks, not_drawn = [], []
    sites = d.get("heteroplasmic_sites")
    for x in sites or []:
        marks.append(mark(f"m{x['pos']}-{x['alt']}", f"m.{x['pos']}{x['ref']}>{x['alt']}", float(x["af"]),
                          chrom="chrM", pos=int(x["pos"]), tone="attention", anchor="mito",
                          locus=f"chrM:{x['pos']}", kind="site", detail=f"allele fraction {x['af']}"))
    if sites is None and d.get("heteroplasmic"):
        not_drawn.append({"what": "heteroplasmic sites this summary does not list (made before it listed them)",
                          "count": int(d["heteroplasmic"])})
    for x in d.get("polyc_heteroplasmic") or []:
        m = POLYC_RE.match(x.get("variant", ""))
        if not m:
            not_drawn.append({"what": f"poly-C call {x.get('variant')!r} with no position", "count": 1})
            continue
        pos = int(m.group(1))
        marks.append(mark(f"polyc-{pos}", f"{x['variant']} (D-loop poly-C tract)", float(x["af"]), chrom="chrM",
                          pos=pos, tone="muted", anchor="mito", locus=f"chrM:{pos}", kind="polyc",
                          detail=f"length heteroplasmy in a D-loop poly-C tract, allele fraction {x['af']}; "
                          "not counted as heteroplasmic"))
    return spec("mito", "mito", "scatter", "Heteroplasmic chrM sites", notes, "MITO_NOTE", marks,
                axis={"min": 0, "max": 1, "unit": "allele fraction", "floor": d.get("heteroplasmy_floor", 0.05)},
                not_drawn=not_drawn)


def build(summary, notes, helpers):
    """The specs of the five charts, in report order, for the sections that ran.
    notes: {caveat key: text}; helpers: render_report's ancestry_fine_rows,
    ANCESTRY_FINE_LEVELS, prs_percentile, prs_note and PRS_NOT_APPLICABLE."""
    S = summary["sections"]
    out = [clinvar_chart(S, notes), genome_chart(S, notes), mito_chart(S, notes),
           ancestry_chart(S, notes, helpers["ancestry_fine_rows"], helpers["ANCESTRY_FINE_LEVELS"]),
           prs_chart(S, notes, helpers["prs_percentile"], helpers["prs_note"], helpers["PRS_NOT_APPLICABLE"])]
    for x in out:
        if x is not None:
            x["svg"] = svg(x)
    return [x for x in out if x is not None]


# --- SVG ------------------------------------------------------------------------------------

def _f(x):
    if not math.isfinite(x):
        raise ValueError(f"a chart coordinate is {x}, not a number")
    return f"{x:.1f}".rstrip("0").rstrip(".")


def _text(x, y, s, size=12, anchor="start", fill="currentColor", extra=""):
    return (f'<text x="{_f(x)}" y="{_f(y)}" font-size="{size}" text-anchor="{anchor}" fill="{fill}"'
            f'{extra}>{E(s)}</text>')


def _a(m, inner):
    """A mark: a link to its anchor in the report, with its title and id."""
    title = f"{m['label']}: {m['detail']}" if m.get("detail") else m["label"]
    return (f'<a href="#{E(m["anchor"] or "")}" data-mark="{E(m["id"])}"><title>{E(title)}</title>{inner}</a>')


def _open(spec_, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" class="chart" data-chart="{E(spec_["id"])}" '
            f'viewBox="0 0 {W} {_f(h)}" width="100%" role="img" aria-label="{E(spec_["title"])}">'
            f'<title>{E(spec_["title"])}</title>')


def _wrap(s, n):
    words, lines, cur = s.split(), [], ""
    for w in words:
        if cur and len(cur) + 1 + len(w) > n:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    return lines + ([cur] if cur else [])


def _footer(spec_, y):
    parts = [f'{x["count"]} {x["what"]}' for x in spec_.get("not_drawn") or []]
    if not parts:
        return "", y
    lines = _wrap("Not drawn: " + "; ".join(parts) + ".", 100)
    return "".join(_text(8, y + 14 * i, s, 11, fill="currentColor", extra=' opacity="0.75"')
                   for i, s in enumerate(lines)), y + 14 * len(lines)


def svg_bars(sp):
    marks = sp["marks"]
    total = sum(m["value"] or 0 for m in marks)
    if total == 0:
        h = 96
        return (_open(sp, h) + f'<rect x="1" y="1" width="{W - 2}" height="{h - 2}" rx="6" fill="none" '
                f'stroke="{FRAME}" stroke-dasharray="6 4"/>'
                + _text(W / 2, 40, "No ClinVar hit in this sample.", 14, "middle")
                + _text(W / 2, 64, "No hit is not an all-clear: the screen reads only the variants ClinVar lists.",
                        12, "middle") + "</svg>")
    top = max(m["value"] or 0 for m in marks)
    left, right, row = 190, 60, 26
    h = 16 + row * len(marks) + 8
    body = []
    for i, m in enumerate(marks):
        y = 12 + i * row
        wbar = (W - left - right) * (m["value"] or 0) / top
        fill = TONE[m["tone"]]
        bar = (f'<rect x="{left}" y="{_f(y)}" width="{_f(max(wbar, 0))}" height="18" rx="3" fill="{fill}"/>'
               if m["value"] else
               f'<rect x="{left}" y="{_f(y)}" width="2" height="18" fill="{FRAME}"/>')
        body.append(_a(m, _text(left - 8, y + 13, m["label"], 12, "end") + bar
                       + _text(left + max(wbar, 2) + 6, y + 13, str(m["value"]), 12)))
    foot, h2 = _footer(sp, h + 4)
    return _open(sp, h2 + 4) + "".join(body) + foot + "</svg>"


def svg_ranges(sp):
    marks = sp["marks"]
    left, right, row = 260, 50, 22
    x0, x1 = left, W - right
    sx = lambda v: x0 + (x1 - x0) * max(0.0, min(100.0, v)) / 100   # noqa: E731
    y, body, group = 22, [], None
    axis = [f'<line x1="{x0}" y1="14" x2="{x1}" y2="14" stroke="{LINE}"/>']
    for t in (0, 25, 50, 75, 100):
        axis.append(_text(sx(t), 10, f"{t}%", 10, "middle"))
    for m in marks:
        if m.get("group") != group:
            group = m.get("group")
            y += 6
            body.append(_text(8, y + 10, group, 12, extra=' font-weight="600"'))
            y += row
        cy = y + 6
        inner = _text(left - 8, cy + 4, m["label"] if len(m["label"]) <= 40 else m["label"][:38] + "...", 12, "end")
        if m["low"] is not None and m["high"] is not None:
            inner += (f'<line class="whisker" x1="{_f(sx(m["low"]))}" y1="{_f(cy)}" x2="{_f(sx(m["high"]))}" '
                      f'y2="{_f(cy)}" stroke="{TONE["info"]}" stroke-width="2"/>')
        inner += f'<circle cx="{_f(sx(m["value"]))}" cy="{_f(cy)}" r="4.5" fill="{TONE["info"]}"/>'
        inner += _text(x1 + 6, cy + 4, f"{m['value']:.0f}%", 11)
        body.append(_a(m, inner))
        y += row
    return _open(sp, y + 8) + "".join(axis) + "".join(body) + "</svg>"


def svg_strip(sp):
    marks = sp["marks"]
    left, right, row = 260, 120, 24
    x0, x1 = left, W - right
    sx = lambda v: x0 + (x1 - x0) * v / 100   # noqa: E731
    body = [f'<line x1="{x0}" y1="14" x2="{x1}" y2="14" stroke="{LINE}"/>']
    for t in (0, 25, 50, 75, 100):
        body.append(_text(sx(t), 10, str(t), 10, "middle"))
    y = 26
    for m in marks:
        cy = y + 6
        inner = _text(left - 8, cy + 4, m["label"] if len(m["label"]) <= 40 else m["label"][:38] + "...", 12, "end")
        inner += f'<line x1="{x0}" y1="{_f(cy)}" x2="{x1}" y2="{_f(cy)}" stroke="{LINE}"/>'
        if m["value"] is None:
            inner += _text(x0 + 4, cy + 4, "no percentile", 11, extra=' opacity="0.8"')
        else:
            inner += (f'<circle cx="{_f(sx(m["value"]))}" cy="{_f(cy)}" r="5" fill="{TONE[m["tone"]]}"/>'
                      + _text(x1 + 8, cy + 4, f"{m['value']:g}", 11))
        body.append(_a(m, inner))
        y += row
    return _open(sp, y + 8) + "".join(body) + "</svg>"


def svg_text(sp):
    lines = _wrap(sp.get("text", ""), 95)
    h = 20 + 16 * len(lines)
    return (_open(sp, h) + f'<rect x="1" y="1" width="{W - 2}" height="{h - 2}" rx="6" fill="none" stroke="{FRAME}"/>'
            + "".join(_text(12, 20 + 16 * i, s, 12) for i, s in enumerate(lines)) + "</svg>")


def svg_ideogram(sp):
    left, right, row = 46, 16, 18
    x0, x1 = left, W - right
    top = max(n for _, n in GRCH38)
    sx = lambda bp: x0 + (x1 - x0) * bp / top   # noqa: E731
    body, rows = [], {c: 10 + i * row for i, (c, _) in enumerate(GRCH38)}
    for c, n in GRCH38:
        y = rows[c]
        body.append(_text(left - 6, y + 9, c[3:], 10, "end"))
        body.append(f'<rect x="{x0}" y="{y + 2}" width="{_f(sx(n) - x0)}" height="8" rx="4" fill="none" '
                    f'stroke="{LINE}"/>')
    for m in sp["marks"]:
        y = rows[m["chrom"]]
        if m.get("kind") == "roh":
            a, b = sx(m["pos"]), sx(m["end"])
            inner = (f'<rect class="roh" x="{_f(a)}" y="{y + 2}" width="{_f(max(b - a, 1.5))}" height="8" '
                     f'fill="{TONE["info"]}"/>')
        else:
            x = sx(m["pos"])
            inner = (f'<rect class="tick" x="{_f(x - 1.5)}" y="{y - 1}" width="3" height="14" '
                     f'fill="{TONE[m["tone"]]}"/>')
        body.append(_a(m, inner))
    y = 10 + len(GRCH38) * row + 6
    legend = (f'<rect x="{x0}" y="{y}" width="14" height="8" fill="{TONE["info"]}"/>'
              + _text(x0 + 20, y + 8, "run of homozygosity, 5 MB or more", 11)
              + f'<rect x="{x0 + 250}" y="{y - 2}" width="3" height="12" fill="{TONE["attention"]}"/>'
              + _text(x0 + 258, y + 8, "ClinVar hit", 11)
              + f'<rect x="{x0 + 340}" y="{y - 2}" width="3" height="12" fill="{TONE["muted"]}"/>'
              + _text(x0 + 348, y + 8, "low-confidence ClinVar call", 11))
    foot, h = _footer(sp, y + 26)
    return _open(sp, h + 4) + "".join(body) + legend + foot + "</svg>"


def svg_scatter(sp):
    left, right, tp, bottom = 44, 16, 12, 30
    h = 200
    x0, x1, y0, y1 = left, W - right, tp, h - bottom
    sx = lambda p: x0 + (x1 - x0) * (p - 1) / (MITO_LEN - 1)   # noqa: E731
    sy = lambda af: y1 - (y1 - y0) * max(0.0, min(1.0, af))   # noqa: E731
    body = [f'<rect x="{x0}" y="{y0}" width="{x1 - x0}" height="{y1 - y0}" fill="none" stroke="{LINE}"/>']
    for t in (0, 0.5, 1):
        body.append(_text(x0 - 6, sy(t) + 4, f"{t:g}", 10, "end"))
    for t, anchor in ((1, "start"), (4000, "middle"), (8000, "middle"), (12000, "middle"), (MITO_LEN, "end")):
        body.append(_text(sx(t), y1 + 14, f"m.{t}", 10, anchor))
    body.append(_text((x0 + x1) / 2, h - 2, "position on chrM; height: allele fraction", 10, "middle",
                      extra=' opacity="0.8"'))
    for lo, hi in ((303, 315), (16184, 16193)):
        body.append(f'<rect x="{_f(sx(lo) - 1)}" y="{y0}" width="{_f(max(sx(hi) - sx(lo), 2))}" '
                    f'height="{y1 - y0}" fill="{TONE["muted"]}" opacity="0.25"/>')
    if not sp["marks"]:
        body.append(_text((x0 + x1) / 2, (y0 + y1) / 2, "No heteroplasmic site between allele fraction 0.05 and 0.95",
                          12, "middle"))
    for m in sp["marks"]:
        x, y = sx(m["pos"]), sy(m["value"])
        if m.get("kind") == "polyc":
            inner = (f'<rect class="polyc" x="{_f(x - 4)}" y="{_f(y - 4)}" width="8" height="8" fill="none" '
                     f'stroke="{TONE["muted"]}" stroke-width="2"/>')
        else:
            inner = f'<circle class="site" cx="{_f(x)}" cy="{_f(y)}" r="5" fill="{TONE["attention"]}"/>'
        body.append(_a(m, inner))
    foot, h2 = _footer(sp, h + 10)
    return _open(sp, h2 + 4) + "".join(body) + foot + "</svg>"


DRAW = {"bars": svg_bars, "ranges": svg_ranges, "strip": svg_strip, "text": svg_text, "ideogram": svg_ideogram,
        "scatter": svg_scatter}


def svg(sp):
    return DRAW[sp["kind"]](sp)
