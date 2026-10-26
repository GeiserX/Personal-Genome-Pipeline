#!/usr/bin/env python3
"""Serve the local genome viewer (igv.js) for one sample on 127.0.0.1.

scripts/pgp-view.sh finds the sample's files and starts this server. It
serves exactly these, and nothing else from the disk:

  /               index.html of this folder
  /viewer.js      the page's script, from this folder
  /igv.min.js     the pinned igv.js release pgp-view.sh checked by sha256
  /config.json    the reference, the tracks and the start locus (build_config)
  /data/<n>/<name>  (name percent-encoded in the URL) the n-th file of the config: the reference FASTA and its
                  .fai, the BAM or CRAM and its index, the VCF and the SV VCF
                  and their .tbi

It binds 127.0.0.1 only (BIND is not an option), answers byte-range requests
(igv.js reads a BAM, CRAM, FASTA or indexed VCF in pieces), refuses a Host
header that is not 127.0.0.1, localhost or [::1] (so a web page cannot reach
it through a DNS name pointed at 127.0.0.1), and sends a Content-Security-
Policy that lets the page talk to this server only: igv.js's own calls to
outside services (gene-name search, BLAT, genome lists) are blocked by the
browser, so nothing about the sample leaves the machine.
"""

import argparse
import http.server
import json
import os
import re
import sys
import urllib.parse

BIND = "127.0.0.1"
HERE = os.path.dirname(os.path.abspath(__file__))
ALLOWED_HOSTS = {"127.0.0.1", "localhost", "[::1]"}
CSP = ("default-src 'none'; script-src 'self'; connect-src 'self'; img-src 'self' data: blob:; "
       "style-src 'self' 'unsafe-inline'; font-src 'self' data:; worker-src 'self' blob:; "
       "base-uri 'none'; form-action 'none'; frame-ancestors 'none'")
LOCUS_RE = re.compile(r"^([A-Za-z0-9._-]+)(?::([0-9][0-9,]*)(?:-([0-9][0-9,]*))?)?$")


def read_fai(fai):
    """Contig names and lengths of a .fai, in order."""
    out = []
    with open(fai) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0]:
                out.append((p[0], int(p[1])))
    if not out:
        raise ValueError(f"{fai} names no contig")
    return out


def check_locus(locus, contigs):
    """The locus as igv.js reads it (chr, chr:pos or chr:start-end), checked
    against the reference; the first contig when none is given."""
    if not locus:
        return contigs[0][0]
    m = LOCUS_RE.match(locus)
    if not m:
        raise ValueError(f"'{locus}' is not a locus: give chr1, chr1:12345 or chr1:12000-13000")
    names = dict(contigs)
    if m.group(1) not in names:
        raise ValueError(f"'{m.group(1)}' is not a contig of the reference")
    for x in (m.group(2), m.group(3)):
        if x and not 1 <= int(x.replace(",", "")) <= names[m.group(1)]:
            raise ValueError(f"{x} is outside {m.group(1)} (1 to {names[m.group(1)]})")
    return locus


def index_of(names):
    """The first of these index files that exists."""
    return next((n for n in names if os.path.isfile(n)), None)


def build_config(sample, fasta, reads=None, vcf=None, sv_vcf=None, locus=None):
    """(config, files): the igv.js browser options and the list of files the
    server may send, /data/<index>/<basename> each."""
    files = []

    def url(path):
        files.append(os.path.abspath(path))
        return f"data/{len(files) - 1}/{urllib.parse.quote(os.path.basename(path))}"

    fai = fasta + ".fai"
    if not os.path.isfile(fasta) or not os.path.isfile(fai):
        raise ValueError(f"the reference needs its FASTA and .fai: {fasta}")
    contigs = read_fai(fai)
    cfg = {"sample": sample,
           "locus": check_locus(locus, contigs),
           "reference": {"id": "reference", "name": os.path.basename(fasta),
                         "fastaURL": url(fasta), "indexURL": url(fai)},
           "tracks": []}
    if reads:
        cram = reads.endswith(".cram")
        stem = os.path.splitext(reads)[0]
        idx = index_of([reads + ".crai", stem + ".crai"] if cram
                       else [reads + ".bai", stem + ".bai", reads + ".csi"])
        if not idx:
            raise ValueError(f"no index beside {reads}: run samtools index on it")
        cfg["tracks"].append({"name": f"{sample} reads", "type": "alignment", "format": "cram" if cram else "bam",
                              "url": url(reads), "indexURL": url(idx)})
    for path, name in ((vcf, f"{sample} small variants"), (sv_vcf, f"{sample} structural variants")):
        if not path:
            continue
        t = {"name": name, "type": "variant", "format": "vcf", "url": url(path)}
        idx = index_of([path + ".tbi", path + ".csi"])
        if idx:
            t["indexURL"] = url(idx)
        else:
            t["indexed"] = False
        cfg["tracks"].append(t)
    return cfg, files


def parse_range(header, size):
    """(start, end) inclusive for a single 'bytes=' range, None when there is
    none or it is not one we serve whole-file for, 'bad' when unsatisfiable."""
    m = re.match(r"^bytes=(\d*)-(\d*)$", (header or "").strip())
    if not m or (not m.group(1) and not m.group(2)):
        return None
    if m.group(1):
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else size - 1
    else:
        start, end = max(0, size - int(m.group(2))), size - 1
    if start >= size or end < start:
        return "bad"
    return start, min(end, size - 1)


class Handler(http.server.BaseHTTPRequestHandler):
    server_version = "pgp-view"
    sys_version = ""

    def log_message(self, fmt, *args):
        if not self.server.quiet:
            sys.stderr.write("%s\n" % (fmt % args))

    def _head(self, code, ctype, length, extra=()):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(length))
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Security-Policy", CSP)
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        for k, v in extra:
            self.send_header(k, v)
        self.end_headers()

    def _fail(self, code, text):
        body = (text + "\n").encode()
        self._head(code, "text/plain; charset=utf-8", len(body))
        if self.command != "HEAD":
            self.wfile.write(body)

    def _target(self):
        host = re.sub(r":[0-9]+$", "", self.headers.get("Host") or "")
        if host not in ALLOWED_HOSTS:
            return None, "the Host header must be 127.0.0.1 or localhost"
        # The config percent-encodes each file name (a space, '#' or '?' in
        # it), so the browser asks for the encoded form; the routes hold the
        # plain name.
        path = urllib.parse.unquote(self.path.split("?", 1)[0].split("#", 1)[0])
        return self.server.routes.get(path), None

    def _serve(self):
        target, err = self._target()
        if err:
            return self._fail(403, err)
        if target is None:
            return self._fail(404, "not found")
        ctype, src = target
        if isinstance(src, bytes):
            data = src
            self._head(200, ctype, len(data))
            if self.command != "HEAD":
                self.wfile.write(data)
            return
        size = os.path.getsize(src)
        rng = parse_range(self.headers.get("Range"), size)
        if rng == "bad":
            return self._head(416, "text/plain; charset=utf-8", 0, [("Content-Range", f"bytes */{size}")])
        start, end = rng if rng else (0, size - 1)
        length = end - start + 1 if size else 0
        extra = [("Content-Range", f"bytes {start}-{end}/{size}")] if rng else []
        self._head(206 if rng else 200, ctype, length, extra)
        if self.command == "HEAD" or not length:
            return
        with open(src, "rb") as f:
            f.seek(start)
            left = length
            while left > 0:
                chunk = f.read(min(1 << 20, left))
                if not chunk:
                    break
                self.wfile.write(chunk)
                left -= len(chunk)

    def do_GET(self):
        try:
            self._serve()
        except (BrokenPipeError, ConnectionResetError):
            pass

    do_HEAD = do_GET


def make_server(cfg, files, igv_js, port, quiet=False):
    """The server on 127.0.0.1:port (0 picks a free port), not started."""
    routes = {"/": ("text/html; charset=utf-8", os.path.join(HERE, "index.html")),
              "/index.html": ("text/html; charset=utf-8", os.path.join(HERE, "index.html")),
              "/viewer.js": ("text/javascript; charset=utf-8", os.path.join(HERE, "viewer.js")),
              "/igv.min.js": ("text/javascript; charset=utf-8", igv_js),
              "/config.json": ("application/json", json.dumps(cfg, indent=1).encode())}
    for i, f in enumerate(files):
        routes[f"/data/{i}/{os.path.basename(f)}"] = ("application/octet-stream", f)
    srv = http.server.ThreadingHTTPServer((BIND, port), Handler)
    srv.daemon_threads = True
    srv.routes, srv.quiet = routes, quiet
    return srv


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fasta", required=True, help="reference FASTA, with its .fai beside it")
    ap.add_argument("--reads", help="BAM or CRAM, with its index beside it")
    ap.add_argument("--vcf", help="small-variant VCF")
    ap.add_argument("--sv-vcf", help="structural-variant VCF")
    ap.add_argument("--locus", help="chr, chr:pos or chr:start-end to open at")
    ap.add_argument("--igv-js", required=True, help="the checked igv.min.js")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args(argv)
    try:
        cfg, files = build_config(a.sample, a.fasta, a.reads, a.vcf, a.sv_vcf, a.locus)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    try:
        srv = make_server(cfg, files, a.igv_js, a.port)
    except OSError as e:
        print(f"ERROR: cannot listen on {BIND}:{a.port}: {e}. Set PGP_VIEW_PORT to a free port.", file=sys.stderr)
        return 1
    print(f"Viewer for {a.sample}: http://{BIND}:{srv.server_address[1]}/#{cfg['locus']}", flush=True)
    print("Press Ctrl-C to stop it.", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
