#!/usr/bin/env python3
"""pgp_live.py: watch a run from a browser while it goes, behind a password.

  pgp_live.py hash-password            read a password (twice from a terminal,
                                       or one line on stdin), print its hash
  pgp_live.py --root GENOME_DIR [--state-dir DIR] [--bind ADDR] [--port N]
              [--ingest-bind ADDR] [--ingest-port N] [--password-file FILE]
              [--allowed-host NAME ...] [--config FILE] ...

Two listeners, Python standard library only:

  ingest   (--ingest-bind 127.0.0.1, --ingest-port 8788) takes POST
           /ingest/<run> only: a Nextflow weblog message, or a pgp_step /
           pgp_run_finished message from scripts/run-all.sh. Each message is
           cut to the fields a progress view needs and appended to
           <state-dir>/runs/<run>/events.jsonl. It has no password: Nextflow's
           weblog cannot send one. Keep it on loopback: a local user can post
           fake progress there, but cannot read anything from it.
  browser  (--bind 127.0.0.1, --port 8787) serves the pages, after a login
           when a password hash is set:
             GET/POST /login, POST /logout
             GET /                          the runs under --root
             GET /assets/<file>, /assets/plugins/<file>
             GET /runs/<run>/               the live page (assets/live/index.html)
             GET /runs/<run>/state.json     bin/run_state.py's state, plus the
                                            headline of each finished section
                                            and what the plugins add
             GET /runs/<run>/events         the stored messages as Server-Sent
                                            Events (Last-Event-ID resumes)
             GET /runs/<run>/report.html    the run's <run>_report.html
             GET /runs/<run>/<relpath>      a file that report links to by a
                                            relative href, or one a plugin
                                            registered; nothing else from disk

The password is never stored: the server reads a PBKDF2-SHA256 hash from
--password-file or PGP_LIVE_PASSWORD_HASH, and refuses to listen on an address
other than loopback without one (exit 2). The access log
(<state-dir>/access.log, TSV) records who opened what, never a query string,
a request body, a result or the password.

Plugins: every bin/pgp_live_*.py (sorted) is imported at start and its
register(app) called; see App below. docs/live-site.md has the deploy guide.
"""
import argparse
import fnmatch
import getpass
import hashlib
import hmac
import html
import importlib.util
import ipaddress
import json
import os
import re
import secrets
import socket
import ssl
import sys
import threading
import time
import urllib.parse
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
ASSETS = os.path.join(REPO, "assets", "live")
sys.path.insert(0, HERE)
import run_state  # noqa: E402

PBKDF2_ITERATIONS = 600000
COOKIE = "pgp_live"
HEARTBEAT_SECONDS = 15
FAIL_WINDOW_SECONDS = 15 * 60
LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1"}
CSP = ("default-src 'none'; script-src 'self'; connect-src 'self'; style-src 'self'; img-src 'self' data:; "
       "form-action 'self'; base-uri 'none'; frame-ancestors 'none'")
# A run's own HTML (the report, CPSR, PharmCAT, MultiQC) carries inline styles and, for third-party
# pages, scripts and outside links. It is sent sandboxed: no script runs, nothing is fetched from
# outside, and it keeps its styles. allow-same-origin lets its links reach the run's other files
# with the session cookie; without allow-scripts that grants the page nothing else.
RUN_HTML_CSP = ("sandbox allow-same-origin allow-popups allow-popups-to-escape-sandbox; default-src 'none'; "
                "style-src 'self' 'unsafe-inline'; img-src 'self' data:; form-action 'none'; base-uri 'none'; "
                "frame-ancestors 'none'")
RUN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")  # validate_sample in scripts/lib/common.sh
ASSET_RE = re.compile(r"^(?:plugins/)?[A-Za-z0-9][A-Za-z0-9._-]*$")
TYPES = {".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8", ".css": "text/css; charset=utf-8",
         ".json": "application/json", ".txt": "text/plain; charset=utf-8", ".tsv": "text/plain; charset=utf-8",
         ".md": "text/plain; charset=utf-8", ".svg": "image/svg+xml", ".png": "image/png", ".pdf": "application/pdf"}
# Weblog fields kept. metadata.parameters and the work directory are paths: they stay in
# events.jsonl for the plugins that need them and never reach a page.
TRACE_KEYS = ("task_id", "hash", "name", "process", "tag", "status", "exit", "attempt", "submit", "start",
              "complete", "realtime", "cpus", "memory", "peak_rss", "peak_vmem", "rchar", "wchar")
PARAM_KEYS = ("outdir", "max_memory", "max_cpus")
STEP_STATES = ("started", "ok", "failed")


def valid_run(name):
    return bool(name) and name not in (".", "..") and bool(RUN_RE.match(name))


def utc_now():
    return datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --- the password ---------------------------------------------------------------------------

def hash_password(password, salt=None, iterations=PBKDF2_ITERATIONS):
    salt = salt or secrets.token_hex(16)
    dk = hashlib.pbkdf2_hmac("sha256", password.encode(), salt.encode(), iterations)
    return f"pbkdf2_sha256${iterations}${salt}${dk.hex()}"


def parse_hash(encoded):
    """(iterations, salt, digest) of a hash-password line; ValueError when it is not one."""
    parts = (encoded or "").strip().split("$")
    if len(parts) != 4 or parts[0] != "pbkdf2_sha256" or not parts[1].isdigit() or not parts[2] \
            or not re.fullmatch(r"[0-9a-f]{64}", parts[3]):
        raise ValueError("not a pgp_live.py hash-password line (pbkdf2_sha256$<iterations>$<salt>$<hash>)")
    return int(parts[1]), parts[2], parts[3]


def verify_password(password, encoded):
    iterations, salt, digest = parse_hash(encoded)
    dk = hashlib.pbkdf2_hmac("sha256", password.encode(), salt.encode(), iterations)
    return hmac.compare_digest(dk.hex(), digest)


def hash_password_main(argv):
    ap = argparse.ArgumentParser(prog="pgp_live.py hash-password",
                                 description="Print the hash of a password for --password-file.")
    ap.parse_args(argv)
    if sys.stdin.isatty():
        pw = getpass.getpass("Password: ")
        if pw != getpass.getpass("Again: "):
            print("ERROR: the two passwords differ", file=sys.stderr)
            return 2
    else:
        pw = sys.stdin.readline().rstrip("\r\n")
    if len(pw) < 8:
        print("ERROR: use a password of 8 characters or more", file=sys.stderr)
        return 2
    print(hash_password(pw))
    return 0


# --- what an ingested message keeps ---------------------------------------------------------

def scalar(v, n=200):
    if isinstance(v, bool) or v is None:
        return v
    if isinstance(v, (int, float)):
        return v
    return str(v)[:n]


def cut(ev):
    """The fields of one message the site keeps; ValueError when it is not a message."""
    if not isinstance(ev, dict):
        raise ValueError("not a JSON object")
    kind = ev.get("event")
    if not isinstance(kind, str) or not re.fullmatch(r"[a-z_]{1,40}", kind):
        raise ValueError("no 'event' name")
    out = {"event": kind}
    for k in ("runName", "runId", "utcTime"):
        if isinstance(ev.get(k), str):
            out[k] = ev[k][:100]
    tr = ev.get("trace")
    if isinstance(tr, dict):
        out["trace"] = {k: scalar(tr[k]) for k in TRACE_KEYS if k in tr and not isinstance(tr[k], (dict, list))}
    md = ev.get("metadata")
    if isinstance(md, dict):
        keep = {}
        params = md.get("parameters")
        if isinstance(params, dict):
            p = {k: scalar(params[k], 1000) for k in PARAM_KEYS if k in params and not isinstance(params[k], (dict, list))}
            if p:
                keep["parameters"] = p
        wf = md.get("workflow")
        if isinstance(wf, dict):
            if isinstance(wf.get("workDir"), str):
                keep["workflow"] = {"workDir": wf["workDir"][:1000]}
            if isinstance(wf.get("errorMessage"), str):
                out["errorMessage"] = wf["errorMessage"][:500]
        if keep:
            out["metadata"] = keep
    if isinstance(ev.get("errorMessage"), str):
        out["errorMessage"] = ev["errorMessage"][:500]
    if kind == "pgp_step":
        step, state = ev.get("step"), ev.get("state")
        if not isinstance(step, str) or not re.fullmatch(r"[0-9]{2}[a-z]?", step) or state not in STEP_STATES:
            raise ValueError("pgp_step needs 'step' (like 14 or 04b) and 'state' (started, ok or failed)")
        out.update(step=step, state=state)
    elif kind == "pgp_run_finished":
        if not isinstance(ev.get("exit"), int) or isinstance(ev.get("exit"), bool):
            raise ValueError("pgp_run_finished needs an integer 'exit'")
        out["exit"] = ev["exit"]
    return out


def step_of_process(steps):
    return {p: s["step"] for s in steps for p in s["processes"]}


def messages(events, steps):
    """What the browser is sent for each stored message, in order: (id, payload). The id is the
    message's line number in events.jsonl. No path, parameter or tag leaves here."""
    owner = step_of_process(steps)
    out = []
    for i, ev in enumerate(events, 1):
        p = {"event": ev.get("event"), "utc": ev.get("utcTime")}
        tr = ev.get("trace")
        if isinstance(tr, dict):
            name = tr.get("name") or tr.get("process") or ""
            proc, _ = run_state.task_name(name)
            p.update(process=proc, step=owner.get(proc), status=tr.get("status"), attempt=tr.get("attempt"))
        if ev.get("event") == "pgp_step":
            p.update(step=ev.get("step"), state=ev.get("state"))
        if ev.get("event") == "pgp_run_finished":
            p["exit"] = ev.get("exit")
        out.append((i, {k: v for k, v in p.items() if v is not None}))
    return out


# --- headlines --------------------------------------------------------------------------------

def _rr():
    import render_report
    return render_report


def _with_note(lines, note):
    """The first line, and the section's caveat when the helper ends with it (never re-worded)."""
    if not lines:
        return None
    return lines[0] + (f" {note}" if note and lines[-1] == note else "")


def _h_clinvar(d, summary):
    rr = _rr()
    n = d.get("count") or 0
    return f"Pathogenic/Likely Pathogenic hits: {n}. {rr.CLINVAR_FREQ_NOTE}", n > 0


def _h_carrier(d, summary):
    lines = _rr().carrier_lines(d)
    return lines[0], bool(d.get("carriers"))


def _h_couple(d, summary):
    rr = _rr()
    return _with_note(rr.couple_lines(d), rr.COUPLE_NOTE), bool(d.get("shared") or d.get("x_linked"))


def _h_roh(d, summary):
    rr = _rr()
    lines = rr.roh_lines(d)
    first = next((x for x in lines if x.startswith("**")), lines[0])
    return f"{first} {rr.ROH_NOTE}" if first != rr.ROH_NOTE else first, bool(d.get("upd_flag"))


def _h_pypgx(d, summary):
    lines = _rr().pypgx_lines(d)
    return (lines[0] if lines else None), bool(d.get("low_gq") or d.get("not_callable"))


def _h_annotsv(d, summary):
    return _rr().annotsv_lines(d)[0], bool(d.get("count"))


def _h_viral(d, summary):
    rr = _rr()
    return _with_note(rr.viral_lines(d), rr.VIRAL_NOTE), bool(d.get("found"))


def _h_archaic(d, summary):
    return _rr().archaic_lines(d)[0], False


def _h_lpa(d, summary):
    return _rr().lpa_lines(d)[0], False


HEADLINES = {"clinvar": _h_clinvar, "carrier": _h_carrier, "couple": _h_couple, "roh": _h_roh,
             "pypgx": _h_pypgx, "annotsv": _h_annotsv, "viral": _h_viral, "archaic": _h_archaic, "lpa": _h_lpa}


def headline(key, sec, summary=None):
    """(text, tone) of one summary section, tone 'info' or 'attention'; None while it has no file.
    The text comes from bin/render_report.py's helpers, so it reads as the report does; no tone says
    a section is fine (a ClinVar screen with no hit is 'info', never an all-clear)."""
    state = (sec or {}).get("state")
    if state in (None, "missing"):
        return None
    title = sec.get("title") or key
    if state in ("failed", "unreadable"):
        return f"{title}: {sec.get('note') or state}", "attention"
    if state == "stale":
        return f"{title}: {_rr().stale_note(sec)}", "attention"
    text, attention = None, False
    fn = HEADLINES.get(key)
    if fn:
        try:
            text, attention = fn(sec.get("data") or {}, summary or {})
        except (KeyError, TypeError, ValueError, AttributeError, IndexError):
            text, attention = None, False
    return (text or f"{title}: finished, see the report"), ("attention" if attention else "info")


# --- the plugin interface ---------------------------------------------------------------------

class App:
    """What a plugin (bin/pgp_live_<name>.py, function register(app)) can add:
      app.add_argument(...)            a flag of its own (before the flags are read)
      app.route(method, pattern, fn, auth=True)
                                       fn(req, match) for a path matching the regex pattern
                                       (whole path); req has send_json(obj, code=200),
                                       send_body(code, ctype, body), run_dir(name), client
      app.run_file(glob)               more files of a run folder the browser may fetch
      app.report_filter(fn)            fn(run, html) -> html, applied to report.html
      app.state_hook(name, fn)         fn(run_dir, state, events) -> dict, in state.json
                                       under plugins[name]
      app.run_kind(fn)                 fn(app) -> [{"name", "dir", "kind"}]: more runs
      app.log(action, run, section, status)  one access-log row for the current request
      app.config                       the parsed flags
    A plugin's page script goes in assets/live/plugins/<name>.js and calls
    window.PGP.register({name, onState, onEvent, panel})."""

    def __init__(self):
        self.parser = None
        self.config = None
        self.routes, self.run_globs, self.report_filters, self.state_hooks, self.run_kinds = [], [], [], [], []
        self.local = threading.local()
        self.access = None

    def add_argument(self, *a, **k):
        self.parser.add_argument(*a, **k)

    def route(self, method, pattern, fn, auth=True):
        self.routes.append((method.upper(), re.compile(pattern), fn, auth))

    def run_file(self, pattern):
        self.run_globs.append(pattern)

    def report_filter(self, fn):
        self.report_filters.append(fn)

    def state_hook(self, name, fn):
        self.state_hooks.append((name, fn))

    def run_kind(self, fn):
        self.run_kinds.append(fn)

    def log(self, action, run="-", section="-", status="-"):
        req = getattr(self.local, "req", None)
        if self.access:
            self.access.write(req.client if req else "-", req.session8 if req else "-", action, run, section, status)


def load_plugins(app, plugin_dir):
    names = sorted(f for f in os.listdir(plugin_dir) if re.fullmatch(r"pgp_live_[A-Za-z0-9_]+\.py", f)) \
        if os.path.isdir(plugin_dir) else []
    for f in names:
        spec = importlib.util.spec_from_file_location(f[:-3], os.path.join(plugin_dir, f))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        if not callable(getattr(mod, "register", None)):
            raise SystemExit(f"ERROR: plugin {f} has no register(app)")
        mod.register(app)
    return names


# --- the access log ---------------------------------------------------------------------------

class AccessLog:
    HEADER = "utc\tclient\tsession8\taction\trun\tsection\tstatus\n"

    def __init__(self, path, max_mb):
        self.path, self.max_bytes, self.lock = path, int(max_mb * 1024 * 1024), threading.Lock()
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    def write(self, *fields):
        row = "\t".join([utc_now()] + [re.sub(r"[\t\r\n]", " ", str(x if x not in (None, "") else "-"))[:200]
                                         for x in fields]) + "\n"
        with self.lock:
            try:
                if self.max_bytes and os.path.getsize(self.path) > self.max_bytes:
                    os.replace(self.path, self.path + ".1")
            except OSError:
                pass
            new = not os.path.exists(self.path)
            with open(self.path, "a") as f:
                if new:
                    f.write(self.HEADER)
                f.write(row)


# --- the stored messages ----------------------------------------------------------------------

class Store:
    """<state-dir>/runs/<run>/events.jsonl, kept in memory as well, with a condition that wakes
    every event stream when a message is added."""

    def __init__(self, state_dir):
        self.dir = os.path.join(state_dir, "runs")
        self.cond = threading.Condition()
        self.runs = {}

    def path(self, run):
        return os.path.join(self.dir, run, "events.jsonl")

    def _load(self, run):
        if run not in self.runs:
            self.runs[run] = run_state.read_events(self.path(run))
        return self.runs[run]

    def known(self):
        try:
            return sorted(r for r in os.listdir(self.dir) if valid_run(r) and os.path.isfile(self.path(r)))
        except OSError:
            return []

    def events(self, run):
        with self.cond:
            return list(self._load(run))

    def append(self, run, ev):
        with self.cond:
            lst = self._load(run)
            os.makedirs(os.path.dirname(self.path(run)), exist_ok=True)
            with open(self.path(run), "a") as f:
                f.write(json.dumps(ev, separators=(",", ":")) + "\n")
            lst.append(ev)
            self.cond.notify_all()
            return len(lst)

    def wait(self, run, after, timeout):
        with self.cond:
            lst = self._load(run)
            self.cond.wait_for(lambda: len(lst) > after, timeout)
            return list(lst[after:])


# --- the server -------------------------------------------------------------------------------

def is_loopback(addr):
    if addr in ("localhost",):
        return True
    try:
        return ipaddress.ip_address(addr.strip("[]")).is_loopback
    except ValueError:
        try:
            return all(ipaddress.ip_address(i[4][0]).is_loopback for i in socket.getaddrinfo(addr, None))
        except (OSError, ValueError):
            return False


class Live:
    """The state both listeners share."""

    def __init__(self, app, cfg, password_hash):
        self.app, self.cfg, self.password_hash = app, cfg, password_hash
        self.root = os.path.realpath(cfg.root)
        self.store = Store(cfg.state_dir)
        self.steps = run_state.read_steps(cfg.steps)
        self.sessions, self.failures, self.lock = {}, {}, threading.Lock()
        self.streams = threading.BoundedSemaphore(cfg.max_streams)
        self.headline_cache = {}
        self.allowed = set(LOOPBACK_HOSTS) | {h.lower() for h in cfg.allowed_host or []}
        if cfg.bind not in ("0.0.0.0", "::"):
            self.allowed.add(cfg.bind.strip("[]").lower())
        self.trust = [ipaddress.ip_network(c, strict=False) for c in cfg.trust_proxy or []]
        self.stopping = False

    # sessions
    def new_session(self):
        token = secrets.token_urlsafe(32)
        now = time.time()
        with self.lock:
            self.sessions[hashlib.sha256(token.encode()).hexdigest()] = {"created": now, "seen": now}
        return token

    def session(self, token):
        """The session's id (sha256 of the cookie) when it is live, else None."""
        if not token:
            return None
        sid = hashlib.sha256(token.encode()).hexdigest()
        now = time.time()
        with self.lock:
            s = self.sessions.get(sid)
            if not s:
                return None
            if now - s["created"] > self.cfg.session_hours * 3600 or now - s["seen"] > self.cfg.session_idle_minutes * 60:
                del self.sessions[sid]
                return None
            s["seen"] = now
        return sid

    def end_session(self, sid):
        with self.lock:
            self.sessions.pop(sid, None)

    def blocked(self, client):
        now = time.time()
        with self.lock:
            f = [t for t in self.failures.get(client, []) if now - t < FAIL_WINDOW_SECONDS]
            self.failures[client] = f
            return len(f) >= self.cfg.login_max_failures

    def failed(self, client):
        with self.lock:
            self.failures.setdefault(client, []).append(time.time())

    # runs
    def runs(self):
        out = {}
        try:
            names = sorted(os.listdir(self.root))
        except OSError:
            names = []
        for n in names:
            d = os.path.join(self.root, n)
            if n.startswith(".") or not valid_run(n) or not self.inside(d, self.root) or not os.path.isdir(d):
                continue
            if any(os.path.exists(os.path.join(d, x)) for x in ("logs/run_status.tsv", "nextflow", "run_manifest.tsv")):
                out[n] = {"name": n, "dir": d, "kind": "sample"}
        for n in self.store.known():
            out.setdefault(n, {"name": n, "dir": None, "kind": "sample"})
        for fn in self.app.run_kinds:
            for r in fn(self.app) or []:
                if valid_run(r.get("name", "")) and r["name"] not in out:
                    out[r["name"]] = {"name": r["name"], "dir": r.get("dir"), "kind": r.get("kind", "other")}
        return out

    @staticmethod
    def inside(path, top):
        real = os.path.realpath(path)
        return real == top or real.startswith(top.rstrip(os.sep) + os.sep)

    def run(self, name):
        return self.runs().get(name) if valid_run(name) else None

    def state(self, run):
        d = run["dir"] or os.path.join(self.root, run["name"])
        events_path = self.store.path(run["name"])
        st = run_state.collect(run["name"], d, events_path=events_path, steps_path=self.cfg.steps)
        events = self.store.events(run["name"])
        if any(e.get("event") == "pgp_run_finished" for e in events) and st["run"]["state"] == "running":
            st["run"]["state"] = "finished"
        st["events"] = len(events)  # the page paints a message ahead of state.json only past this id
        st["headlines"] = self.headlines(run["name"], d, st)
        st["plugins"] = {}
        for name, fn in self.app.state_hooks:
            st["plugins"][name] = fn(d, st, events)
        return st

    def headlines(self, sample, d, st):
        """{section: {text, tone, step}} for the sections of the steps that are done. The summary
        is read again only when the set of done steps changes."""
        done = tuple(s["step"] for s in st["steps"] if s["state"] == "done")
        key = (sample, done)
        hit = self.headline_cache.get(sample)
        if hit and hit[0] == key:
            return hit[1]
        out = {}
        if done and os.path.isdir(d):
            import collect_summary
            try:
                summary = collect_summary.collect(sample, d)
            except Exception as e:  # a reader bug must not take the page down
                print(f"pgp_live: summary of {sample} not read: {type(e).__name__}", file=sys.stderr)
                summary = None
            for s in st["steps"]:
                if s["state"] != "done" or not summary:
                    continue
                for k in s["sections"]:
                    h = headline(k, summary["sections"].get(k), summary)
                    if h:
                        out[k] = {"text": h[0], "tone": h[1], "step": s["step"]}
        self.headline_cache[sample] = (key, out)
        return out

    def run_files(self, run):
        """Relative paths of the run folder the browser may fetch: the report's relative links and
        the plugins' globs (checked when asked)."""
        d = run["dir"]
        rep = os.path.join(d, f"{run['name']}_report.html") if d else None
        out = set()
        if rep and os.path.isfile(rep):
            with open(rep, errors="replace") as f:
                for h in re.findall(r'href="([^"#?]+)', f.read()):
                    h = html.unescape(h)
                    if re.match(r"^[a-z][a-z0-9+.-]*:", h, re.I) or h.startswith("/"):
                        continue
                    rel = os.path.normpath(urllib.parse.unquote(h))
                    if not rel.startswith("..") and rel != ".":
                        out.add(rel)
        return out

    def run_file_allowed(self, run, rel):
        return rel in self.run_files(run) or any(fnmatch.fnmatch(rel, g) for g in self.app.run_globs)

    def index_page(self):
        assets = self.cfg.assets
        with open(os.path.join(assets, "index.html")) as f:
            return f.read()


class BrowserHandler(BaseHTTPRequestHandler):
    server_version = "pgp-live"
    sys_version = ""
    live = None  # set per server

    # helpers
    def log_message(self, fmt, *args):
        if not self.live.cfg.quiet:
            sys.stderr.write("pgp_live: %s %s\n" % (self.command, urllib.parse.urlsplit(self.path).path))

    @property
    def client(self):
        peer = self.client_address[0]
        if self.live.trust:
            try:
                if any(ipaddress.ip_address(peer) in n for n in self.live.trust):
                    xff = [x.strip() for x in (self.headers.get("X-Forwarded-For") or "").split(",") if x.strip()]
                    if xff:
                        return xff[-1]
            except ValueError:
                pass
        return peer

    def headers_common(self, csp=CSP):
        self.send_header("Content-Security-Policy", csp)
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Frame-Options", "DENY")

    def send_body(self, code, ctype, body, headers=(), csp=CSP):
        if isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        for k, v in headers:
            self.send_header(k, v)
        self.headers_common(csp)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def send_json(self, obj, code=200):
        self.send_body(code, "application/json", json.dumps(obj))

    def deny(self, code, text=None):
        self.send_body(code, "text/plain; charset=utf-8", (text or {401: "login required", 403: "forbidden",
                                                                   404: "not found", 429: "too many failed logins"}
                                                           .get(code, "error")) + "\n")

    def run_dir(self, name):
        r = self.live.run(name)
        return r and r["dir"]

    def cookie_token(self):
        for part in (self.headers.get("Cookie") or "").split(";"):
            k, _, v = part.strip().partition("=")
            if k == COOKIE:
                return v
        return None

    def log(self, action, run="-", section="-", status="-"):
        self.live.app.log(action, run, section, status)

    def host_ok(self):
        h = (self.headers.get("Host") or "").strip().lower()
        m = re.match(r"^(\[[^\]]*\]|[^:]+)(?::\d+)?$", h)
        return bool(m) and m.group(1).strip("[]") in {x.strip("[]") for x in self.live.allowed}

    def origin_ok(self):
        """A POST from another site is refused. Browsers say where a request comes from in
        Sec-Fetch-Site; with Referrer-Policy no-referrer they send 'Origin: null' for our own forms."""
        site = self.headers.get("Sec-Fetch-Site")
        if site:
            return site in ("same-origin", "none")
        o = self.headers.get("Origin")
        if not o or o == "null":
            return True
        return urllib.parse.urlsplit(o).netloc.lower() == (self.headers.get("Host") or "").lower()

    # dispatch
    def do_HEAD(self):
        self.handle_any("GET")

    def do_GET(self):
        self.handle_any("GET")

    def do_POST(self):
        self.handle_any("POST")

    def handle_any(self, method):
        self.live.app.local.req = self
        self.sid = None
        self.session8 = "-"
        if not self.host_ok():
            self.log("host-refused", status=403)
            return self.deny(403)
        if method == "POST" and not self.origin_ok():
            self.log("origin-refused", status=403)
            return self.deny(403)
        path = urllib.parse.urlsplit(self.path).path
        if self.live.password_hash:
            self.sid = self.live.session(self.cookie_token())
            self.session8 = self.sid[:8] if self.sid else "-"
        authed = self.sid is not None or not self.live.password_hash
        if path == "/login":
            return self.login(method)
        if path == "/logout" and method == "POST":
            if self.sid:
                self.live.end_session(self.sid)
            self.log("logout", status=303)
            return self.send_body(303, "text/plain", "", [("Location", "/login"), ("Set-Cookie", self.cookie("", 0))])
        if method == "GET" and path.startswith("/assets/"):
            return self.asset(path[len("/assets/"):])
        for m_, rx, fn, auth in self.live.app.routes:
            mt = rx.fullmatch(path)
            if mt and m_ == method:
                if auth and not authed:
                    self.log("plugin", status=401)
                    return self.deny(401)
                return fn(self, mt)
        if method != "GET":
            return self.deny(404)
        if not authed:
            page = path == "/" or path.endswith("/")
            self.log("no-session", status=303 if page else 401)
            if page:
                return self.send_body(303, "text/plain", "", [("Location", "/login?next=" + urllib.parse.quote(path))])
            return self.deny(401)
        if path == "/":
            return self.index()
        m = re.fullmatch(r"/runs/([^/]+)(/.*)?", path)
        if not m:
            return self.deny(404)
        run = self.live.run(urllib.parse.unquote(m.group(1)))
        if not run:
            self.log("run", status=404)
            return self.deny(404)
        rest = m.group(2)
        if rest is None:
            return self.send_body(301, "text/plain", "", [("Location", f"/runs/{run['name']}/")])
        return self.run_path(run, rest[1:])

    def cookie(self, value, max_age):
        secure = "; Secure" if (self.live.cfg.cookie_secure or self.live.cfg.tls_cert) else ""
        return f"{COOKIE}={value}; Path=/; HttpOnly; SameSite=Strict; Max-Age={max_age}{secure}"

    def login(self, method):
        if method == "GET":
            return self.login_page(200, "")
        client = self.client
        if not self.live.password_hash:
            return self.send_body(303, "text/plain", "", [("Location", "/")])
        if self.live.blocked(client):
            self.log("login-blocked", status=429)
            return self.login_page(429, "Too many failed logins. Wait 15 minutes.", [("Retry-After", "900")])
        n = int(self.headers.get("Content-Length") or 0)
        if n > 4096:
            return self.deny(413, "too large")
        form = urllib.parse.parse_qs(self.rfile.read(n).decode(errors="replace"))
        pw = (form.get("password") or [""])[0]
        nxt = (form.get("next") or ["/"])[0]
        if not re.fullmatch(r"/[A-Za-z0-9._~/%-]*", nxt) or nxt.startswith("//"):
            nxt = "/"
        try:
            ok = verify_password(pw, self.live.password_hash)
        except ValueError:
            ok = False
        if not ok:
            self.live.failed(client)
            self.log("login-fail", status=401)
            return self.login_page(401, "Wrong password.")
        with self.live.lock:
            self.live.failures.pop(client, None)
        token = self.live.new_session()
        self.session8 = hashlib.sha256(token.encode()).hexdigest()[:8]
        self.log("login", status=303)
        return self.send_body(303, "text/plain", "", [("Location", nxt),
                                                       ("Set-Cookie", self.cookie(token, int(self.live.cfg.session_hours * 3600)))])

    def login_page(self, code, message, headers=()):
        with open(os.path.join(self.live.cfg.assets, "login.html")) as f:
            page = f.read()
        q = urllib.parse.parse_qs(urllib.parse.urlsplit(self.path).query)
        nxt = (q.get("next") or ["/"])[0]
        page = page.replace("{{message}}", html.escape(message)).replace("{{next}}", html.escape(nxt, quote=True))
        self.send_body(code, TYPES[".html"], page, headers)

    def asset(self, name):
        if not ASSET_RE.match(name):
            return self.deny(404)
        p = os.path.join(self.live.cfg.assets, name)
        top = os.path.realpath(self.live.cfg.assets)
        if not os.path.isfile(p) or not Live.inside(p, top) or name in ("index.html", "login.html"):
            return self.deny(404)
        with open(p, "rb") as f:
            self.send_body(200, TYPES.get(os.path.splitext(p)[1], "application/octet-stream"), f.read())

    def index(self):
        runs = self.live.runs()
        items = "\n".join(f'<li><a href="/runs/{html.escape(n)}/">{html.escape(n)}</a> '
                          f'<span class="kind">{html.escape(r["kind"])}</span></li>' for n, r in sorted(runs.items()))
        self.log("index", status=200)
        page = ('<!doctype html>\n<html lang="en" data-theme="dark"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width, initial-scale=1"><title>Runs</title>'
                '<link rel="stylesheet" href="/assets/live.css"></head><body><header><h1>Runs</h1>'
                '<form method="post" action="/logout"><button type="submit">Log out</button></form></header>'
                f'<main><ul class="runs">{items or "<li>No run under this folder yet.</li>"}</ul></main></body></html>\n')
        self.send_body(200, TYPES[".html"], page)

    def run_page(self, run):
        plugins_dir = os.path.join(self.live.cfg.assets, "plugins")
        js = sorted(f for f in os.listdir(plugins_dir) if ASSET_RE.match(f) and f.endswith(".js")) \
            if os.path.isdir(plugins_dir) else []
        tags = "\n".join(f'<script src="/assets/plugins/{html.escape(f)}" defer></script>' for f in js)
        page = (self.live.index_page().replace("{{assets}}", "/assets").replace("{{run}}", html.escape(run["name"]))
                .replace("{{mode}}", "live").replace("{{plugins}}", tags))
        self.log("run", run["name"], status=200)
        self.send_body(200, TYPES[".html"], page)

    def run_path(self, run, rest):
        name = run["name"]
        if rest == "":
            return self.run_page(run)
        if rest == "state.json":
            self.log("state", name, status=200)
            return self.send_json(self.live.state(run))
        if rest == "events":
            return self.events(run)
        d = run["dir"]
        if rest == "report.html":
            p = os.path.join(d, f"{name}_report.html") if d else None
            if not p or not os.path.isfile(p) or not Live.inside(p, os.path.realpath(d)):
                self.log("report", name, status=404)
                return self.deny(404)
            with open(p, errors="replace") as f:
                page = f.read()
            for fn in self.live.app.report_filters:
                page = fn(name, page)
            self.log("report", name, status=200)
            return self.send_body(200, TYPES[".html"], page, csp=RUN_HTML_CSP)
        rel = os.path.normpath(urllib.parse.unquote(rest))
        if not d or rel.startswith("..") or os.path.isabs(rel) or not self.live.run_file_allowed(run, rel):
            self.log("file", name, rel, 404)
            return self.deny(404)
        p = os.path.join(d, rel)
        if not os.path.isfile(p) or not Live.inside(p, os.path.realpath(d)) or not Live.inside(p, self.live.root):
            self.log("file", name, rel, 404)
            return self.deny(404)
        ext = os.path.splitext(p)[1].lower()
        with open(p, "rb") as f:
            body = f.read()
        self.log("file", name, rel, 200)
        self.send_body(200, TYPES.get(ext, "application/octet-stream"), body,
                       csp=RUN_HTML_CSP if ext in (".html", ".htm", ".svg") else CSP)

    def events(self, run):
        name = run["name"]
        if not self.live.streams.acquire(blocking=False):
            self.log("events", name, status=503)
            return self.deny(503, "too many open event streams")
        try:
            try:
                after = int(self.headers.get("Last-Event-ID") or 0)
            except ValueError:
                after = 0
            self.log("events", name, status=200)
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("X-Accel-Buffering", "no")
            self.headers_common()
            self.end_headers()
            self.wfile.write(b"retry: 3000\n\n")
            self.wfile.flush()
            sent = after
            while not self.live.stopping:
                new = self.live.store.wait(name, sent, HEARTBEAT_SECONDS)
                if not new:
                    self.wfile.write(b": heartbeat\n\n")
                else:
                    all_events = self.live.store.events(name)
                    for i, payload in messages(all_events, self.live.steps)[sent:sent + len(new)]:
                        self.wfile.write(f"id: {i}\ndata: {json.dumps(payload, separators=(',', ':'))}\n\n".encode())
                    sent += len(new)
                self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError, ssl.SSLError, OSError):
            pass
        finally:
            self.live.streams.release()
        self.close_connection = True


class IngestHandler(BaseHTTPRequestHandler):
    server_version = "pgp-live-ingest"
    sys_version = ""
    live = None

    def log_message(self, fmt, *args):
        pass

    def reply(self, code, text):
        body = (text + "\n").encode()
        self.send_response(code)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        self.reply(404, "not found")

    def do_POST(self):
        path = re.sub(r"^/+", "/", urllib.parse.urlsplit(self.path).path)  # a PGP_LIVE_URL ending in '/'
        if not path.startswith("/ingest/"):
            return self.reply(404, "not found")
        if self.headers.get("Origin"):
            return self.reply(403, "a browser cannot post here")  # Nextflow and curl send no Origin
        run = urllib.parse.unquote(path[len("/ingest/"):])
        if not valid_run(run):
            return self.reply(400, "bad run name")
        try:
            n = int(self.headers.get("Content-Length"))
        except (TypeError, ValueError):
            return self.reply(411, "Content-Length required")
        if n > self.live.cfg.max_body_kb * 1024:
            self.close_connection = True
            return self.reply(413, "too large")
        raw = self.rfile.read(n)
        try:
            ev = cut(json.loads(raw.decode("utf-8")))
        except (ValueError, UnicodeDecodeError) as e:
            return self.reply(400, f"not a message: {e}"[:200])
        self.live.store.append(run, ev)
        self.reply(200, "ok")


def make_server(handler, live, bind, port, tls=None):
    cls = type(handler.__name__, (handler,), {"live": live})
    srv_cls = type("Server", (ThreadingHTTPServer,), {"daemon_threads": True, "allow_reuse_address": True,
                                                      "address_family": socket.AF_INET6 if ":" in bind else socket.AF_INET})
    srv = srv_cls((bind.strip("[]"), port), cls)
    if tls:
        ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
        ctx.load_cert_chain(*tls)
        srv.socket = ctx.wrap_socket(srv.socket, server_side=True)
    return srv


# --- the flags ----------------------------------------------------------------------------------

def parser():
    ap = argparse.ArgumentParser(prog="pgp_live.py", description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", required=True, help="GENOME_DIR: the folder holding one folder per run")
    ap.add_argument("--state-dir", help="where messages and the access log are written (default ROOT/.pgp-live)")
    ap.add_argument("--bind", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8787)
    ap.add_argument("--ingest-bind", default="127.0.0.1")
    ap.add_argument("--ingest-port", type=int, default=8788)
    ap.add_argument("--password-file", help="a file holding one hash-password line (else PGP_LIVE_PASSWORD_HASH)")
    ap.add_argument("--allowed-host", action="append", default=[],
                    help="a Host name the browser may use besides loopback and the bind address (repeatable)")
    ap.add_argument("--trust-proxy", action="append", default=[],
                    help="a reverse proxy's address or CIDR: only its X-Forwarded-For is read (repeatable)")
    ap.add_argument("--tls-cert")
    ap.add_argument("--tls-key")
    ap.add_argument("--cookie-secure", action="store_true", help="mark the cookie Secure (TLS ends at a proxy)")
    ap.add_argument("--session-hours", type=float, default=12, help="a login lasts at most this long")
    ap.add_argument("--session-idle-minutes", type=float, default=120, help="a login unused this long ends")
    ap.add_argument("--max-streams", type=int, default=32, help="open event streams at once")
    ap.add_argument("--max-body-kb", type=int, default=1024, help="largest ingested message")
    ap.add_argument("--login-max-failures", type=int, default=5, help="failed logins per client in 15 minutes")
    ap.add_argument("--access-log", help="default STATE_DIR/access.log")
    ap.add_argument("--access-log-max-mb", type=float, default=10, help="then it moves to access.log.1")
    ap.add_argument("--config", help="KEY=VALUE lines setting any of these flags (flag names with _), "
                                     "read before the flags")
    ap.add_argument("--plugin-dir", default=HERE, help=argparse.SUPPRESS)
    ap.add_argument("--assets", default=ASSETS, help=argparse.SUPPRESS)
    ap.add_argument("--steps", default=run_state.STEPS_TSV, help=argparse.SUPPRESS)
    ap.add_argument("--quiet", action="store_true", help="no request line on stderr")
    return ap


def read_config(path, ap):
    known = {a.dest: a for a in ap._actions}
    out = {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            k, sep, v = line.partition("=")
            k, v = k.strip().replace("-", "_"), v.strip()
            if not sep or k not in known or k in ("config", "help"):
                raise SystemExit(f"ERROR: {path}:{n}: not a KEY=VALUE line of a known flag")
            a = known[k]
            if isinstance(a.default, list):
                out.setdefault(k, []).append(v)
            elif a.const is True:
                out[k] = v.lower() in ("1", "true", "yes")
            else:
                out[k] = a.type(v) if a.type else v
    return out


def configure(argv, app):
    ap = parser()
    pre, _ = ap.parse_known_args(argv)
    app.parser = ap
    load_plugins(app, pre.plugin_dir)
    pre, _ = ap.parse_known_args(argv)
    if pre.config:
        ap.set_defaults(**read_config(pre.config, ap))
    cfg = ap.parse_args(argv)
    cfg.state_dir = cfg.state_dir or os.path.join(cfg.root, ".pgp-live")
    cfg.access_log = cfg.access_log or os.path.join(cfg.state_dir, "access.log")
    app.config = cfg
    return cfg


def read_password_hash(cfg):
    if cfg.password_file:
        with open(cfg.password_file) as f:
            h = f.read().strip()
    else:
        h = os.environ.get("PGP_LIVE_PASSWORD_HASH", "").strip()
    if h:
        parse_hash(h)
    return h or None


def start(argv, app=None):
    """Configure and bind both listeners; returns (live, browser server, ingest server)."""
    app = app or App()
    cfg = configure(argv, app)
    try:
        pw = read_password_hash(cfg)
    except (OSError, ValueError) as e:
        raise SystemExit(f"ERROR: password hash: {e}") from None
    if not pw and not is_loopback(cfg.bind):
        print(f"ERROR: --bind {cfg.bind} is not loopback and no password hash is set: run "
              "'pgp_live.py hash-password' and pass --password-file (or PGP_LIVE_PASSWORD_HASH)", file=sys.stderr)
        raise SystemExit(2)
    if bool(cfg.tls_cert) != bool(cfg.tls_key):
        raise SystemExit("ERROR: --tls-cert and --tls-key go together")
    if not os.path.isdir(cfg.root):
        raise SystemExit(f"ERROR: --root {cfg.root} is not a folder")
    os.makedirs(cfg.state_dir, exist_ok=True)
    live = Live(app, cfg, pw)
    app.access = AccessLog(cfg.access_log, cfg.access_log_max_mb)
    if not pw:
        print("pgp_live: no password hash set: anyone who can reach the loopback port sees the runs", file=sys.stderr)
    if not is_loopback(cfg.ingest_bind):
        print(f"pgp_live: WARNING: the ingest port listens on {cfg.ingest_bind}: anyone who reaches it can post "
              "fake progress (it serves nothing)", file=sys.stderr)
    tls = (cfg.tls_cert, cfg.tls_key) if cfg.tls_cert else None
    web = make_server(BrowserHandler, live, cfg.bind, cfg.port, tls)
    ingest = make_server(IngestHandler, live, cfg.ingest_bind, cfg.ingest_port)
    return live, web, ingest


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["hash-password"]:
        return hash_password_main(argv[1:])
    if argv[:1] == ["serve"]:
        argv = argv[1:]
    live, web, ingest = start(argv)
    threading.Thread(target=ingest.serve_forever, daemon=True).start()
    scheme = "https" if live.cfg.tls_cert else "http"
    print(f"pgp_live: browser {scheme}://{live.cfg.bind}:{web.server_address[1]}/ ; "
          f"ingest http://{live.cfg.ingest_bind}:{ingest.server_address[1]}/ingest/<run>", file=sys.stderr)
    try:
        web.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        live.stopping = True
        web.shutdown()
        ingest.shutdown()
    return 0


if __name__ == "__main__":
    sys.exit(main())
