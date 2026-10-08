#!/usr/bin/env python3
"""check_pictures.py: check the docs pictures before they are committed.

  check_pictures.py [--ocr] [--terms-file FILE] [--expect DEMO-001] PNG...
  check_pictures.py --self-test [--ocr] [--terms-file FILE]

For every PNG:
  chunks  only IHDR, PLTE, IDAT and IEND: no text, time, EXIF or colour
          profile chunk that could carry a name, a path or a date
  ocr     (--ocr, macOS only: the Vision framework through pyobjc) the text
          read from the picture contains the --expect text and none of the
          private terms: home and volume paths, private IPv4 addresses,
          tailnet names, plus one case-insensitive regular expression per
          line of --terms-file. Keep that file outside the repository: it
          lists your own names and host names.

--self-test plants one bad picture for each check (a tEXt chunk, a home path,
a private address, a line of the terms file, a picture without the expected
text) and fails unless every one of them is caught. Run it before trusting a
pass. Needs pillow; --ocr needs pyobjc-framework-Vision.
"""
import argparse
import os
import re
import struct
import sys
import tempfile
import zlib

ALLOWED = {b"IHDR", b"PLTE", b"IDAT", b"IEND"}
OCTET = r"(?:25[0-5]|2[0-4][0-9]|1[0-9][0-9]|[1-9]?[0-9])"
BUILT_IN = [
    ("home or volume path", re.compile(r"/(?:Users|home|Volumes|mnt/user)/", re.I)),
    ("private address", re.compile(r"(?<![0-9.])(?:10\." + OCTET + r"|172\.(?:1[6-9]|2[0-9]|3[01])|192\.168|"
                                   r"100\.(?:6[4-9]|[7-9][0-9]|1[01][0-9]|12[0-7]))\." + OCTET + r"\." + OCTET)),
    ("tailnet name", re.compile(r"\.ts\.net\b", re.I)),
]


def chunks(path):
    with open(path, "rb") as f:
        data = f.read()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")
    out, i = [], 8
    while i < len(data):
        n, kind = struct.unpack(">I4s", data[i:i + 8])
        out.append(kind)
        i += 12 + n
    return out


def check_chunks(path):
    bad = sorted({k.decode("latin-1") for k in chunks(path) if k not in ALLOWED})
    return [f"chunk {k}" for k in bad]


def ocr_text(path):
    """The text Vision reads, tile by tile: a tall picture is cut into
    1000-pixel bands (with overlap), as Vision scales a big image down."""
    import Vision
    from Foundation import NSURL
    from PIL import Image

    with Image.open(path) as im:
        im = im.convert("RGB")
    texts = []
    with tempfile.TemporaryDirectory() as tmp:
        top, n = 0, 0
        while top < im.height:
            tile = os.path.join(tmp, f"t{n}.png")
            im.crop((0, top, im.width, min(im.height, top + 1000))).save(tile)
            req = Vision.VNRecognizeTextRequest.alloc().init()
            req.setRecognitionLevel_(0)          # accurate
            req.setUsesLanguageCorrection_(False)
            handler = Vision.VNImageRequestHandler.alloc().initWithURL_options_(NSURL.fileURLWithPath_(tile), None)
            ok, err = handler.performRequests_error_([req], None)
            if not ok:
                raise RuntimeError(f"Vision failed on {path}: {err}")
            texts += [r.topCandidates_(1)[0].string() for r in req.results() or []]
            top += 900
            n += 1
    return "\n".join(texts)


def load_terms(path):
    """The built-in patterns plus every active line of PATH. A PATH with no
    active line is refused: a check that reads it would pass on nothing."""
    terms = list(BUILT_IN)
    if path:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    terms.append(("private term", re.compile(line, re.I)))
        if len(terms) == len(BUILT_IN):
            raise SystemExit(f"ERROR: the terms file {path} has no active pattern (only blank or # lines)")
    return terms


def check_ocr(path, expect, terms):
    text = ocr_text(path)
    flat = " ".join(text.split())
    problems = []
    if expect and expect not in flat and expect.replace("-", "") not in flat.replace("-", "").replace("—", ""):
        problems.append(f"OCR did not find {expect!r}")
    for i, (kind, rx) in enumerate(terms, 1):
        if rx.search(text) or rx.search(flat):
            # Name the kind, not the match: the match may be the private term itself.
            problems.append(f"OCR found a {kind} (pattern {i} of the list)")
    return problems, len(text)


def check(paths, ocr, expect, terms):
    failed = 0
    for p in paths:
        problems = check_chunks(p)
        note = ""
        if ocr:
            more, n = check_ocr(p, expect, terms)
            problems += more
            note = f", {n} characters read"
        if problems:
            failed += 1
            print(f"FAIL {p}: " + "; ".join(problems))
        else:
            print(f"PASS {p}: image chunks only" + (f"; {expect} found, no private term{note}" if ocr else ""))
    return failed


def plant(tmp, name, text, chunk=None):
    from PIL import Image, ImageDraw, ImageFont
    fonts = ["/System/Library/Fonts/Menlo.ttc", "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"]
    font = ImageFont.truetype(next(f for f in fonts if os.path.isfile(f)), 22)
    im = Image.new("RGB", (900, 160), (255, 255, 255))
    ImageDraw.Draw(im).multiline_text((20, 20), text, font=font, fill=(0, 0, 0), spacing=12)
    path = os.path.join(tmp, name)
    im.save(path)
    if chunk:
        with open(path, "rb") as f:
            data = f.read()
        body = b"Comment\x00" + chunk.encode()
        extra = struct.pack(">I", len(body)) + b"tEXt" + body + struct.pack(">I", zlib.crc32(b"tEXt" + body))
        with open(path, "wb") as f:
            f.write(data[:33] + extra + data[33:])   # after IHDR
    return path


def self_test(ocr, terms_file):
    terms = load_terms(terms_file)
    with tempfile.TemporaryDirectory() as tmp:
        good = plant(tmp, "good.png", "Sample: DEMO-001\nGenomic Analysis Report")
        cases = [("a clean picture passes", good, True)]
        cases.append(("a tEXt chunk is caught", plant(tmp, "text.png", "Sample: DEMO-001", chunk="written on a laptop"), False))
        if ocr:
            # Built here, so this file holds no path or address the personal-data check would flag.
            home = "/" + "/".join(["Users", "someone", "genome"])
            addr = ".".join(["192", "168", "20", "31"])
            cases += [
                ("a home path is caught", plant(tmp, "home.png", f"Sample: DEMO-001\nOutput: {home}"), False),
                ("a private address is caught", plant(tmp, "ip.png", f"Sample: DEMO-001\nhost {addr}"), False),
                ("a picture without DEMO-001 is caught", plant(tmp, "noid.png", "Sample: OTHER-002"), False),
            ]
            first = next((t[1].pattern for t in terms if t[0] == "private term"), None)
            if terms_file and first:
                word = re.sub(r"[^A-Za-z0-9 ]", "", first) or "x"
                cases.append(("the first line of the terms file is caught",
                              plant(tmp, "term.png", f"Sample: DEMO-001\nRun by {word}"), False))
        bad = 0
        for desc, path, should_pass in cases:
            ok = check([path], ocr, "DEMO-001", terms) == 0
            right = ok == should_pass
            bad += not right
            print(f"[{'PASS' if right else 'FAIL'}] self-test: {desc}")
    if bad:
        print(f"self-test: {bad} planted case(s) were not caught as expected", file=sys.stderr)
        return 1
    print("self-test: every planted case behaved as expected")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pngs", nargs="*")
    ap.add_argument("--ocr", action="store_true", help="read the text with macOS Vision and check it")
    ap.add_argument("--terms-file", help="one private term (a regular expression) per line; keep it outside the repository")
    ap.add_argument("--expect", default="DEMO-001", help="text every picture must contain (default DEMO-001)")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test(a.ocr, a.terms_file)
    if not a.pngs:
        ap.error("give the PNGs to check")
    return 1 if check(a.pngs, a.ocr, a.expect, load_terms(a.terms_file)) else 0


if __name__ == "__main__":
    sys.exit(main())
