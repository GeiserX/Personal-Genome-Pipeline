"""owngenome_live_viewer.py: the genome viewer on the live site, a plugin of bin/owngenome_live.py.

scripts/owngenome-view.sh serves one sample's reads in igv.js on 127.0.0.1,
with no login. This plugin serves the same page for every run of the site,
behind the site's password, so a position in the report opens the reads from
any device that can reach the site:

  GET /runs/<run>/viewer/              the page (assets/viewer/index.html)
  GET /runs/<run>/viewer/viewer.js     its script (assets/viewer/viewer.js)
  GET /runs/<run>/viewer/igv.min.js    the pinned igv.js (IGV_JS_VERSION in versions.env)
  GET /runs/<run>/viewer/config.json   the reference and the run's tracks
  GET /runs/<run>/viewer/data/<n>/<name>  the n-th file of that config, by byte range

The files are the ones owngenome-view.sh picks by default
(assets/viewer/serve.py sample_files): the run folder's BAM or CRAM, its VCF
and its SV VCF, each with its index, and the reference. A track whose file
or index lies outside GENOME_DIR (a symlink out of it) is left out. The page and the files carry the
viewer's Content-Security-Policy (serve.CSP): the page talks to this site
only. The report the site serves (report.html) gets its viewer links pointed
here. Opening the page writes one access-log row; the byte ranges igv.js reads
write none.

Flags it adds to owngenome_live.py (also as KEY=VALUE lines in --config;
scripts/owngenome-live.sh passes both):

  --viewer-fasta FILE    the reference FASTA, with its .fai beside it (REF_FASTA)
  --viewer-igv-js FILE   the igv.min.js scripts/owngenome-live.sh checked by sha256

Without both the viewer answers 404 with the reason, and the report keeps its
links to the local viewer. Standard library only.
"""
import importlib.util
import json
import os
import re
import sys
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
VIEWER_DIR = os.path.join(os.path.dirname(HERE), "assets", "viewer")
sys.path.insert(0, HERE)
import render_report  # noqa: E402

_spec = importlib.util.spec_from_file_location("owngenome_viewer_serve", os.path.join(VIEWER_DIR, "serve.py"))
serve = importlib.util.module_from_spec(_spec)
_keep, sys.dont_write_bytecode = sys.dont_write_bytecode, True  # no __pycache__ in assets/viewer
try:
    _spec.loader.exec_module(serve)
finally:
    sys.dont_write_bytecode = _keep

LOCAL_TEXT = "served from this computer only."
LIVE_TEXT = "served by the run site, behind its login."
LIVE_NOTE = ('    <p class="viewer">The Viewer column opens the reads at each row in this site\'s genome viewer '
             "(docs/live-site.md).</p>")
DATA_RE = re.compile(r"data/([0-9]+)/(.+)")
STATIC = {"viewer.js": "text/javascript; charset=utf-8"}

APP = None  # the App register() was given; its config holds the flags once they are read


def flag(name):
    return getattr(APP.config, name, None) if APP and APP.config else None


def ready():
    """Why the viewer cannot run, or None when it can."""
    fasta, igv = flag("viewer_fasta"), flag("viewer_igv_js")
    if not fasta or not igv:
        return "the genome viewer is not set up: start the site with --viewer-fasta and --viewer-igv-js " \
               "(scripts/owngenome-live.sh passes both)"
    if not os.path.isfile(fasta) or not os.path.isfile(fasta + ".fai"):
        return "the genome viewer's reference or its .fai is missing"
    if not os.path.isfile(igv):
        return "igv.min.js is missing"
    return None


def run_config(run, run_dir, root):
    """(config, files) of one run: serve.build_config on the run folder's files that lie inside
    GENOME_DIR. ValueError when the reference or an index is wrong."""
    top = os.path.realpath(run_dir)

    def inside(p):
        real = os.path.realpath(p)
        return all(real.startswith(t.rstrip(os.sep) + os.sep) for t in (top, root))
    reads, vcf, sv_vcf = (p if p and inside(p) else None for p in serve.sample_files(run_dir, run))
    cfg, files = serve.build_config(run, flag("viewer_fasta"), reads, vcf, sv_vcf)
    # An index beside a track can be a symlink out of GENOME_DIR too: then the track goes, with
    # its file and its index. Files 0 and 1 are the reference and its .fai (--viewer-fasta, any folder).
    files = [f if n < 2 or inside(f) else None for n, f in enumerate(files)]

    def nums(t):
        return [int(t[k].split("/")[1]) for k in ("url", "indexURL") if k in t]
    for t in cfg["tracks"]:
        if not all(files[n] for n in nums(t)):
            for n in nums(t):
                files[n] = None
    cfg["tracks"] = [t for t in cfg["tracks"] if all(files[n] for n in nums(t))]
    return cfg, files


def head_writer(req):
    def head(code, ctype, length, extra=()):
        req.send_response(code)
        req.send_header("Content-Type", ctype)
        req.send_header("Content-Length", str(length))
        req.send_header("Cross-Origin-Resource-Policy", "same-origin")
        for k, v in extra:
            req.send_header(k, v)
        req.headers_common(serve.CSP)
        req.end_headers()
    return head


def send(req, code, ctype, body):
    body = body.encode() if isinstance(body, str) else body
    head_writer(req)(code, ctype, len(body))
    if req.command != "HEAD":
        req.wfile.write(body)


def viewer(req, m):
    run, rest = urllib.parse.unquote(m.group(1)), m.group(2)
    page = rest in ("", "index.html")
    # The page asks for a login like every page of the site; the files answer 401.
    if req.live.password_hash and req.sid is None:
        req.log("viewer", status=303 if page else 401)
        if page:
            return req.send_body(303, "text/plain", "", [("Location", "/login?next=" + urllib.parse.quote(
                urllib.parse.urlsplit(req.path).path))])
        return req.deny(401)
    d = req.run_dir(run)
    if not d:
        req.log("viewer", status=404)
        return req.deny(404)
    why = ready()
    if why:
        req.log("viewer", run, status=404)
        return req.deny(404, why)
    try:
        if page:
            with open(os.path.join(VIEWER_DIR, "index.html")) as f:
                html = f.read().replace(LOCAL_TEXT, LIVE_TEXT)
            req.log("viewer", run, status=200)
            return send(req, 200, "text/html; charset=utf-8", html)
        if rest in STATIC:
            with open(os.path.join(VIEWER_DIR, rest), "rb") as f:
                return send(req, 200, STATIC[rest], f.read())
        if rest == "igv.min.js":
            return serve.send_file(req, head_writer(req), "text/javascript; charset=utf-8", flag("viewer_igv_js"))
        cfg, files = run_config(run, d, req.live.root)
        if rest == "config.json":
            return send(req, 200, "application/json", json.dumps(cfg, indent=1))
        dm = DATA_RE.fullmatch(rest)
        # The config's URLs quote the name (serve.build_config); the path is still quoted here.
        n = int(dm.group(1)) if dm else len(files)
        if n < len(files) and files[n] and os.path.basename(files[n]) == urllib.parse.unquote(dm.group(2)):
            return serve.send_file(req, head_writer(req), "application/octet-stream", files[n])
    except ValueError as e:
        req.log("viewer", run, status=409)
        return req.deny(409, f"the genome viewer cannot open this run: {e}")
    except (BrokenPipeError, ConnectionResetError):
        return None
    return req.deny(404)


def to_slash(req, m):
    run = urllib.parse.unquote(m.group(1))
    if not req.run_dir(run):
        return req.deny(404)
    return req.send_body(301, "text/plain", "", [("Location", f"/runs/{urllib.parse.quote(run)}/viewer/")])


def point_links(run, html):
    """The served report's viewer links and notes, pointed at this site's viewer of the run."""
    if ready():
        return html
    here = f"/runs/{urllib.parse.quote(run)}/viewer/"
    return (html.replace(f'href="{render_report.VIEWER_URL}#', f'href="{here}#')
            .replace(render_report.VIEWER_NOTE, LIVE_NOTE))


def register(app):
    global APP
    APP = app
    app.add_argument("--viewer-fasta", help="reference FASTA (with its .fai) for the genome viewer")
    app.add_argument("--viewer-igv-js", help="the checked igv.min.js for the genome viewer")
    app.route("GET", r"/runs/([^/]+)/viewer", to_slash)
    # auth=False: viewer() asks for the login itself, so the page sends the browser to it.
    app.route("GET", r"/runs/([^/]+)/viewer/(.*)", viewer, auth=False)
    app.report_filter(point_links)
