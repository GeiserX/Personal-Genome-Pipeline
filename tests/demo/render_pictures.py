#!/usr/bin/env python3
"""render_pictures.py: the three docs pictures from DEMO-001's reports.

  render_pictures.py --genome-dir DIR --out OUTDIR [--font FILE]
                     [--html-stop-before TITLES] [--cpic-head N]

reads DIR/DEMO-001 after steps 24, 27 and 28 ran on it and writes
  demo-html-report.png     the step 24 report in a 1180-pixel-wide browser,
                           from the top down to the first card whose title
                           starts with one of TITLES (comma-separated); by
                           default the first table that lists genes or
                           variants, so no gene of the invented sample shows
  demo-cpic-report.png     the step 27 text report as `cat` prints it in a
                           900-pixel-wide dark terminal, drawn from the file
                           with Pillow and a monospace font (the same file and
                           font give the same picture); --cpic-head N draws
                           `head -n N` instead
  demo-multiqc-report.png  the top 1600 x 1118 pixels of the step 28 report

Every PNG is written again from its pixels alone, so it holds no text, time
or EXIF chunk. Needs playwright (with its Chromium) and pillow;
tests/demo/retake-screenshots.sh installs both in a scratch venv.
"""
import argparse
import os
import sys

from PIL import Image, ImageDraw, ImageFont

SAMPLE = "DEMO-001"
# Cards of the step 24 report that list genes or variants, in the order the
# report has them. The picture stops above the first one present.
STOP_BEFORE = ["Secondary-Findings Genes", "Genes With a Non-Normal", "ClinVar Hits", "Steps Not Run"]
FONTS = ["/System/Library/Fonts/Menlo.ttc", "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"]


def clean_png(path):
    """Write PATH again from its pixels only (IHDR, IDAT, IEND)."""
    with Image.open(path) as im:
        im.load()
        rgb = im.convert("RGB")
    Image.frombytes("RGB", rgb.size, rgb.tobytes()).save(path, format="PNG", optimize=True)


def html_pictures(sample_dir, out, stop_titles):
    from playwright.sync_api import sync_playwright

    report = os.path.join(sample_dir, f"{SAMPLE}_report.html")
    multiqc = os.path.join(sample_dir, "multiqc", "multiqc_report.html")
    for p in (report, multiqc):
        if not os.path.isfile(p):
            raise SystemExit(f"ERROR: {p} not found: run steps 24 and 28 first")
    with sync_playwright() as pw:
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={"width": 1180, "height": 900}, device_scale_factor=1)
        page.goto("file://" + os.path.abspath(report))
        page.wait_for_load_state("load")
        cut = page.evaluate(
            """(titles) => {
                 const cards = [...document.querySelectorAll('.card')];
                 const stop = cards.find(c => { const h = c.querySelector('h2');
                   return h && titles.some(t => h.textContent.trim().startsWith(t)); });
                 if (!stop) return null;
                 return Math.floor(stop.getBoundingClientRect().top + window.scrollY);
               }""", stop_titles)
        if cut is None:
            raise SystemExit(f"ERROR: no card titled {stop_titles} in {report}; the crop has nothing to stop above")
        # The grid gap is 20 pixels: leave 12 below the last card shown.
        path = os.path.join(out, "demo-html-report.png")
        page.screenshot(path=path, full_page=True, clip={"x": 0, "y": 0, "width": 1180, "height": cut - 8})
        clean_png(path)
        print(f"  {path}: stops above the card at y={cut}")

        page = browser.new_page(viewport={"width": 1600, "height": 1118}, device_scale_factor=1)
        page.goto("file://" + os.path.abspath(multiqc))
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)   # the plots draw after load
        path = os.path.join(out, "demo-multiqc-report.png")
        page.screenshot(path=path)
        clean_png(path)
        print(f"  {path}")
        browser.close()


def cpic_picture(sample_dir, out, font_path, head):
    rel = f"{SAMPLE}/cpic/{SAMPLE}_cpic_recommendations.txt"
    with open(os.path.join(os.path.dirname(sample_dir), rel)) as f:
        lines = f.read().splitlines()
    command = f"cat {rel}"
    if head:
        lines, command = lines[:head], f"head -n {head} {rel}"
    font = ImageFont.truetype(font_path, 14)
    width, pad, bar, step = 900, 22, 30, 20
    advance = font.getlength("M")
    cols = int((width - 2 * pad) // advance)
    rows = [("$ ", command)]
    for line in lines:
        # A terminal wraps a long line at its last column.
        for i in range(0, max(len(line), 1), cols):
            rows.append(("", line[i:i + cols]))
    height = bar + pad + len(rows) * step + pad
    im = Image.new("RGB", (width, height), (30, 30, 40))
    d = ImageDraw.Draw(im)
    d.rectangle([0, 0, width, bar], fill=(44, 44, 56))
    for i, colour in enumerate([(255, 95, 87), (254, 188, 46), (40, 200, 64)]):
        x = 17 + 20 * i
        d.ellipse([x - 6, bar // 2 - 6, x + 6, bar // 2 + 6], fill=colour)
    y = bar + pad
    for prompt, text in rows:
        x = pad
        if prompt:
            d.text((x, y), prompt, font=font, fill=(80, 250, 123))
            x += advance * len(prompt)
        d.text((x, y), text, font=font, fill=(222, 222, 228))
        y += step
    path = os.path.join(out, "demo-cpic-report.png")
    im.save(path, format="PNG", optimize=True)
    clean_png(path)
    print(f"  {path}: {len(rows)} terminal lines of {cols} columns")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--font", help="monospace TrueType font for the CPIC picture (default: Menlo, else DejaVu Sans Mono)")
    ap.add_argument("--html-stop-before", default=",".join(STOP_BEFORE),
                    help="comma-separated card titles (prefixes); the HTML picture stops above the first one")
    ap.add_argument("--cpic-head", type=int, default=0, help="draw only the first N lines (head -n N)")
    ap.add_argument("--only", choices=["html", "cpic"], help="draw only these pictures (html: the step 24 and MultiQC ones)")
    a = ap.parse_args(argv)
    font = a.font or next((f for f in FONTS if os.path.isfile(f)), None)
    if not font:
        ap.error("no monospace font found: give --font")
    sample_dir = os.path.join(a.genome_dir, SAMPLE)
    os.makedirs(a.out, exist_ok=True)
    if a.only != "html":
        cpic_picture(sample_dir, a.out, font, a.cpic_head)
    if a.only != "cpic":
        html_pictures(sample_dir, a.out, [t.strip() for t in a.html_stop_before.split(",") if t.strip()])
    return 0


if __name__ == "__main__":
    sys.exit(main())
