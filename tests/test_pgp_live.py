#!/usr/bin/env python3
"""bin/pgp_live.py and assets/live/ on the synthetic DEMO-001 run (tests/fixtures/live).

Checks:
  1. hash-password: the hash verifies the password and no other; the server
     reads only a hash; a non-loopback --bind without one exits 2;
  2. login: no cookie gives a redirect to /login (pages) or 401 (data); a wrong
     password is 401 with a login-fail row; the right one sets an HttpOnly,
     SameSite=Strict cookie; an expired session is 401; after 5 failures the
     6th try is 429; the plain password is in no log, no stderr and no stored
     message;
  3. Host: a foreign Host header is 403, an --allowed-host name is served;
  4. ingest: POST /ingest on the browser port is 404; a 2 MiB body is 413,
     a body that is not JSON 400, a run named '../x' 400; metadata.parameters
     (an HPO list, a VCF path) and a task's script are not stored;
  5. the event stream sends a posted message within 2 s, Last-Event-ID
     resumes after the id given, and what it sends equals messages();
  6. state.json's step states are run_state's for the same files; ClinVar's
     headline is never 'ok' and carries the report's caveat;
  7. every browser response carries the Content-Security-Policy, and
     assets/live/ holds no outside URL, no inline script and no on* attribute;
  8. a plugin (tests/fixtures/live/plugins) adds a flag, a route, a run file
     and a state hook, all behind the session; the page gets one script tag
     per file in assets/live/plugins;
  9. run files: a file the report links to is served, an unlinked one and a
     symlink out of the root are 404;
 10. in a headless Chrome (skipped without Chrome, except in CI): the login
     form, the step map, a box that repaints after a posted message, the
     replay mode playing tests/fixtures/live/frames.json, and no policy
     violation.
Checks 2, 4, 7 and 8 are also run against a planted mistake and must fail
there: a check that cannot fail is no check.

Run: python3 tests/test_pgp_live.py
"""
import contextlib
import http.client
import io
import json
import os
import re
import select
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import urllib.parse
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import pgp_live  # noqa: E402
import run_state  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "live")
ASSETS = os.path.join(REPO, "assets", "live")
RUN = "DEMO-001"
PASSWORD = "synthetic test password 7"
FAILS = []
RESPONSES = []  # (path, headers) of every browser response


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS.append(desc)


def control(desc, still_passes):
    """A planted mistake: the check must fail on it."""
    check(f"control: {desc}", not still_passes, "the check passed on the planted mistake")


def request(port, method, path, body=None, headers=None, record=True):
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    h = {"Host": f"127.0.0.1:{port}"}
    h.update(headers or {})
    c.request(method, path, body=body, headers=h)
    r = c.getresponse()
    data = r.read()
    hdr = {k.lower(): v for k, v in r.getheaders()}
    c.close()
    if record:
        RESPONSES.append((path, hdr))
    return r.status, hdr, data


def make_root(work):
    root = os.path.join(work, "genome")
    d = os.path.join(root, RUN)
    shutil.copytree(os.path.join(FIX, RUN), d)
    os.makedirs(os.path.join(d, "clinvar"))
    os.makedirs(os.path.join(d, "extra"))
    with open(os.path.join(d, "clinvar", "hits.tsv"), "w") as f:
        f.write("gene\tvariant\nGENEA\tchr1:100 A>G\n")
    with open(os.path.join(d, "clinvar", "unlinked.tsv"), "w") as f:
        f.write("not linked\n")
    with open(os.path.join(d, "extra", "note.txt"), "w") as f:
        f.write("plugin file\n")
    outside = os.path.join(work, "outside.txt")
    with open(outside, "w") as f:
        f.write("outside the root\n")
    os.symlink(outside, os.path.join(d, "clinvar", "escape.tsv"))
    with open(os.path.join(d, f"{RUN}_report.html"), "w") as f:
        f.write('<!doctype html><html><head><style>body{color:#eee}</style></head><body>'
                '<a href="clinvar/hits.tsv">hits</a> <a href="clinvar/escape.tsv">escape</a> '
                '<a href="https://example.org/x">outside</a></body></html>\n')
    return root


def serve(argv, app=None):
    live, web, ingest = pgp_live.start(argv, app)
    for s in (web, ingest):
        threading.Thread(target=s.serve_forever, daemon=True).start()
    return live, web, ingest


def stop(*servers):
    for live, web, ingest in servers:
        live.stopping = True
        web.shutdown()
        ingest.shutdown()
        web.server_close()
        ingest.server_close()


def login(port, password=PASSWORD):
    code, h, _ = request(port, "POST", "/login", urllib.parse.urlencode({"password": password, "next": "/"}),
                         {"Content-Type": "application/x-www-form-urlencoded"})
    m = re.match(r"pgp_live=([^;]+)", h.get("set-cookie", ""))
    return code, h, (m.group(1) if m else None)


def post(port, run, ev):
    body = ev if isinstance(ev, (bytes, str)) else json.dumps(ev)
    return request(port, "POST", f"/ingest/{run}", body, {"Content-Type": "application/json"}, record=False)


class Stream:
    """An open /events stream, read in a thread."""

    def __init__(self, port, cookie, last=None):
        self.items, self.sock = [], socket.create_connection(("127.0.0.1", port), timeout=30)
        hdr = f"GET /runs/{RUN}/events HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nCookie: pgp_live={cookie}\r\n"
        if last is not None:
            hdr += f"Last-Event-ID: {last}\r\n"
        self.sock.sendall((hdr + "\r\n").encode())
        threading.Thread(target=self.read, daemon=True).start()

    def read(self):
        buf, cur = b"", {}
        try:
            while True:
                chunk = self.sock.recv(65536)
                if not chunk:
                    return
                buf += chunk
                while b"\n" in buf:
                    line, buf = buf.split(b"\n", 1)
                    line = line.decode().rstrip("\r")
                    if line.startswith("id: "):
                        cur["id"] = int(line[4:])
                    elif line.startswith("data: "):
                        cur["data"] = json.loads(line[6:])
                    elif line == "" and "data" in cur:
                        self.items.append((cur["id"], cur["data"]))
                        cur = {}
        except OSError:
            return

    def wait(self, n, timeout):
        end = time.time() + timeout
        while time.time() < end and len(self.items) < n:
            time.sleep(0.02)
        return len(self.items) >= n

    def close(self):
        try:
            self.sock.close()
        except OSError:
            pass


def asset_problems(texts):
    """(file, problem) for each outside URL, inline script or on* attribute in assets/live."""
    out = []
    for name, t in texts.items():
        for u in re.findall(r"(?:https?:)?//[A-Za-z0-9.-]+\.[a-z]{2,}[^\s'\"]*", t):
            if u != "http://www.w3.org/2000/svg":  # the SVG namespace name, never fetched
                out.append((name, u))
        if name.endswith(".html"):
            for tag in re.findall(r"<script\b[^>]*>", t, re.I):
                if "src=" not in tag:
                    out.append((name, tag))
            out += [(name, a) for a in re.findall(r"\son[a-z]+\s*=", t, re.I)]
        if name.endswith(".js") and re.search(r"\beval\s*\(|new\s+Function\s*\(", t):
            out.append((name, "eval"))
    return out


# --- the browser over the DevTools pipe ------------------------------------------------------

def find_chrome():
    for c in (os.environ.get("CHROME"), "google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
              "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"):
        if c and (shutil.which(c) or os.path.isfile(c)):
            return shutil.which(c) or c
    return None


class Chrome:
    """Headless Chrome driven through --remote-debugging-pipe (fd 3 in, fd 4 out): no package needed."""

    def __init__(self, chrome, profile):
        to_r, self.to_w = os.pipe()
        self.from_r, from_w = os.pipe()

        def fds():
            os.dup2(to_r, 3)
            os.dup2(from_w, 4)
        self.p = subprocess.Popen([chrome, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-first-run",
                                   "--no-default-browser-check", f"--user-data-dir={profile}",
                                   "--remote-debugging-pipe", "about:blank"],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, preexec_fn=fds, pass_fds=(3, 4))
        os.close(to_r)
        os.close(from_w)
        self.n, self.buf = 0, b""
        t = self.call("Target.createTarget", {"url": "about:blank"})["targetId"]
        self.session = self.call("Target.attachToTarget", {"targetId": t, "flatten": True})["sessionId"]

    def call(self, method, params=None, session=None, timeout=30):
        self.n += 1
        msg = {"id": self.n, "method": method, "params": params or {}}
        if session:
            msg["sessionId"] = session
        os.write(self.to_w, json.dumps(msg).encode() + b"\0")
        end = time.time() + timeout
        while time.time() < end:
            while b"\0" in self.buf:
                raw, self.buf = self.buf.split(b"\0", 1)
                m = json.loads(raw)
                if m.get("id") == self.n:
                    if "error" in m:
                        raise RuntimeError(f"{method}: {m['error']}")
                    return m.get("result", {})
            r, _, _ = select.select([self.from_r], [], [], 0.5)
            if r:
                chunk = os.read(self.from_r, 1 << 16)
                if not chunk:
                    raise RuntimeError("Chrome closed the pipe")
                self.buf += chunk
        raise RuntimeError(f"{method}: no answer in {timeout} s")

    def eval(self, expr):
        r = self.call("Runtime.evaluate", {"expression": expr, "returnByValue": True}, self.session)
        return r.get("result", {}).get("value")

    def go(self, url):
        self.call("Page.navigate", {"url": url}, self.session)

    def wait(self, expr, timeout=20):
        end = time.time() + timeout
        while time.time() < end:
            try:
                v = self.eval(expr)
            except RuntimeError:
                v = None
            if v:
                return v
            time.sleep(0.2)
        return None

    def close(self):
        try:
            self.call("Browser.close", timeout=5)
        except (RuntimeError, OSError):
            pass
        try:
            self.p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            self.p.kill()


class ReplayHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Content-Security-Policy", pgp_live.CSP)
        super().end_headers()

    def log_message(self, *a):
        pass


def browser_checks(port, work):
    chrome = find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the page in")
        if ci:
            FAILS.append("no Chrome in CI")
        return
    b = Chrome(chrome, os.path.join(work, "chrome"))
    replay_srv = None
    try:
        base = f"http://127.0.0.1:{port}"
        b.go(f"{base}/runs/{RUN}/")
        check("browser: a page without a session goes to the login form",
              b.wait("location.pathname === '/login' && !!document.getElementById('password')"),
              b.eval("location.href"))
        b.eval(f"document.getElementById('password').value = {json.dumps(PASSWORD)}; "
               "document.querySelector('form').submit(); true")
        ok = b.wait(f"location.pathname === '/runs/{RUN}/' && document.documentElement.dataset.pgpStatus === 'ready'")
        check("browser: the login form leads back to the run page, which loads", ok,
              b.eval("location.href + ' ' + document.documentElement.outerHTML.slice(0, 600)"))
        n = b.eval("document.querySelectorAll('#map g.box').length")
        check("browser: the step map draws one box per step shown", isinstance(n, int) and n >= 4, n)
        check("browser: a finished step is painted done",
              b.eval("document.querySelector('g.box[data-step=\"06\"]').dataset.state") == "done",
              b.eval("document.querySelector('g.box[data-step=\"06\"]') && document.querySelector('g.box[data-step=\"06\"]').dataset.state"))
        before = b.eval("document.querySelector('g.box[data-step=\"07\"]').dataset.state")
        post(INGEST_PORT, RUN, {"event": "process_started", "utcTime": "2026-01-02T10:07:00Z",
                                "trace": {"task_id": 20, "name": f"PGX:PHARMCAT_PREPROCESS ({RUN})",
                                          "process": "PGX:PHARMCAT_PREPROCESS", "tag": RUN, "status": "RUNNING",
                                          "start": 1767348420000}})
        after = b.wait("document.querySelector('g.box[data-step=\"07\"]').dataset.state === 'running' && 'running'", 10)
        check("browser: a box repaints after a posted message (07 waiting -> running)",
              before == "waiting" and after == "running", (before, after))
        b.eval("document.querySelector('g.box[data-step=\"11\"]').dispatchEvent(new MouseEvent('click', {bubbles: true})); true")
        detail = b.wait("document.getElementById('detail').textContent.indexOf('Step 11') >= 0 && "
                        "document.getElementById('detail').textContent")
        check("browser: a click on a box opens its detail with its tasks", detail and "ROH" in detail, detail)
        blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
        check("browser: the policy blocked nothing on the live page", not blocked, blocked)

        # replay: the same page and script on frames.json, served as static files
        rd = os.path.join(work, "replay")
        os.makedirs(rd)
        for f in ("live.js", "live.css"):
            shutil.copy(os.path.join(ASSETS, f), rd)
        shutil.copy(os.path.join(FIX, "frames.json"), rd)
        with open(os.path.join(ASSETS, "index.html")) as f:
            page = f.read().replace("{{assets}}", ".").replace("{{run}}", RUN).replace("{{mode}}", "replay") \
                .replace("{{plugins}}", "")
        with open(os.path.join(rd, "index.html"), "w") as f:
            f.write(page)
        handler = type("H", (ReplayHandler,), {"__init__": lambda self, *a, **k: SimpleHTTPRequestHandler.__init__(
            self, *a, directory=rd, **k)})
        replay_srv = ThreadingHTTPServer(("127.0.0.1", 0), handler)
        threading.Thread(target=replay_srv.serve_forever, daemon=True).start()
        b.go(f"http://127.0.0.1:{replay_srv.server_address[1]}/index.html")
        first = b.wait("document.documentElement.dataset.pgpFrame === '0' && "
                       "document.querySelector('g.box[data-step=\"11\"]').dataset.state")
        b.eval("document.getElementById('replay-play').click(); true")
        last = b.wait("document.documentElement.dataset.pgpFrame === '2' && "
                      "document.querySelector('g.box[data-step=\"11\"]').dataset.state", 10)
        found = b.eval("document.getElementById('found').textContent")
        check("browser: replay plays frames.json (11 waiting at frame 0, done at the last, its headline listed)",
              first == "waiting" and last == "done" and "Runs of homozygosity" in (found or ""), (first, last, found))
        blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
        check("browser: the policy blocked nothing in replay", not blocked, blocked)
    finally:
        b.close()
        if replay_srv:
            replay_srv.shutdown()
            replay_srv.server_close()


INGEST_PORT = None


def main():
    global INGEST_PORT
    work = tempfile.mkdtemp(prefix="pgp-live-")
    servers = []
    err = io.StringIO()
    try:
        root = make_root(work)
        state_dir = os.path.join(work, "state")
        pwfile = os.path.join(work, "password")

        # 1. the password
        r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "pgp_live.py"), "hash-password"],
                           input=PASSWORD + "\n", capture_output=True, text=True, timeout=60)
        h = r.stdout.strip()
        check("hash-password prints a pbkdf2_sha256 hash of 600000 rounds, not the password",
              r.returncode == 0 and h.startswith("pbkdf2_sha256$600000$") and PASSWORD not in r.stdout + r.stderr,
              (r.returncode, r.stderr))
        check("the hash verifies the password and no other",
              pgp_live.verify_password(PASSWORD, h) and not pgp_live.verify_password(PASSWORD + "x", h))
        with open(pwfile, "w") as f:
            f.write(h + "\n")
        env = {k: v for k, v in os.environ.items() if k != "PGP_LIVE_PASSWORD_HASH"}
        r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "pgp_live.py"), "--root", root, "--state-dir",
                            state_dir, "--bind", "0.0.0.0", "--port", "0", "--ingest-port", "0"],
                           capture_output=True, text=True, timeout=60, env=env)
        check("a non-loopback --bind without a password hash exits 2", r.returncode == 2 and "password" in r.stderr,
              (r.returncode, r.stderr))

        app = pgp_live.App()
        argv = ["--root", root, "--state-dir", state_dir, "--port", "0", "--ingest-port", "0", "--password-file",
                pwfile, "--allowed-host", "live.example.org", "--plugin-dir", os.path.join(FIX, "plugins"),
                "--fixture-word", "planted-flag", "--quiet"]
        with contextlib.redirect_stderr(err):
            srv = serve(argv, app)
        servers.append(srv)
        live, web, ingest = srv
        port, INGEST_PORT = web.server_address[1], ingest.server_address[1]
        for line in open(os.path.join(FIX, "events.jsonl")):
            code, _, _ = post(INGEST_PORT, RUN, line.strip())
            if code != 200:
                check("the fixture's messages are taken", False, code)

        # 2. login
        code, hd, _ = request(port, "GET", f"/runs/{RUN}/")
        code2, _, _ = request(port, "GET", f"/runs/{RUN}/state.json")
        check("no cookie: a page redirects to /login, data is 401",
              code == 303 and hd.get("location", "").startswith("/login") and code2 == 401, (code, hd, code2))
        with contextlib.redirect_stderr(err):
            code, _, body = login(port, "wrong password")
        log = open(os.path.join(state_dir, "access.log")).read()
        check("a wrong password is 401 and writes a login-fail row", code == 401 and "\tlogin-fail\t" in log, (code, log))
        with contextlib.redirect_stderr(err):
            code, hd, cookie = login(port)
        sc = hd.get("set-cookie", "")
        check("the right password sets an HttpOnly, SameSite=Strict session cookie",
              code == 303 and cookie and "HttpOnly" in sc and "SameSite=Strict" in sc, (code, sc))
        auth = {"Cookie": f"pgp_live={cookie}"}
        code, _, _ = request(port, "GET", f"/runs/{RUN}/state.json", headers=auth)
        check("a valid session reads state.json", code == 200, code)
        code, _, _ = request(port, "GET", "/", headers={"Cookie": "pgp_live=made-up"})
        check("a made-up cookie is not a session", code == 303, code)
        _, _, old = login(port)
        sid = __import__("hashlib").sha256(old.encode()).hexdigest()
        live.sessions[sid]["created"] -= live.cfg.session_hours * 3600 + 1
        code, _, _ = request(port, "GET", f"/runs/{RUN}/state.json", headers={"Cookie": f"pgp_live={old}"})
        check("an expired session is 401", code == 401, code)
        live.failures.clear()
        codes = [login(port, f"wrong {i}")[0] for i in range(6)]
        check("five wrong passwords are 401, the sixth try is 429", codes == [401] * 5 + [429], codes)
        code, _, _ = login(port)
        check("the right password is refused too while blocked", code == 429, code)
        live.failures.clear()

        def leaks():
            files = [os.path.join(state_dir, "access.log"), live.store.path(RUN)]
            return any(PASSWORD in open(p).read() for p in files if os.path.exists(p)) or PASSWORD in err.getvalue()
        check("the plain password is in no log, no stderr and no stored message", not leaks())
        app.access.write("127.0.0.1", "-", "planted", RUN, PASSWORD, 200)
        control("the password written into the access log", not leaks())
        with open(os.path.join(state_dir, "access.log"), "w") as f:
            f.write(pgp_live.AccessLog.HEADER)

        # 3. Host
        code, _, _ = request(port, "GET", "/", headers=dict(auth, Host="evil.example.net"))
        code2, _, _ = request(port, "GET", "/", headers=dict(auth, Host="live.example.org:8443"))
        check("a foreign Host is 403, an --allowed-host name is served", code == 403 and code2 == 200, (code, code2))

        # 4. ingest
        code, _, _ = request(port, "POST", f"/ingest/{RUN}", "{}", auth)
        check("POST /ingest on the browser port is 404", code == 404, code)
        # The server answers 413 from the Content-Length, before the body: send the head and a little.
        sk = socket.create_connection(("127.0.0.1", INGEST_PORT), timeout=10)
        sk.sendall(f"POST /ingest/{RUN} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\n"
                   f"Content-Length: {2 * 1024 * 1024}\r\n\r\n".encode() + b"x" * 1024)
        head = sk.recv(200).decode(errors="replace")
        sk.close()
        check("a 2 MiB message is 413", head.startswith("HTTP/1.0 413") or head.startswith("HTTP/1.1 413"), head)
        code, _, _ = post(INGEST_PORT, RUN, "not json")
        check("a message that is not JSON is 400", code == 400, code)
        code, _, _ = post(INGEST_PORT, "../x", {"event": "started"})
        code2, _, _ = post(INGEST_PORT, "-x", {"event": "started"})
        check("a run named '../x' or '-x' is 400", code == 400 and code2 == 400, (code, code2))
        code, _, _ = request(INGEST_PORT, "POST", f"/ingest/{RUN}", '{"event":"started"}',
                             {"Origin": "https://evil.example.net"}, record=False)
        check("a browser's post (an Origin header) is 403", code == 403, code)
        stored = open(live.store.path(RUN)).read()

        def stored_leaks(text):
            return [x for x in ("HP:0001250", "DEMO-001.vcf.gz", "bcftools roh", "commandLine") if x in text]
        check("metadata.parameters, the command line and the task script are not stored",
              not stored_leaks(stored) and '"outdir":"/data/genome"' in stored, stored_leaks(stored))
        control("the raw message stored", not stored_leaks(open(os.path.join(FIX, "events.jsonl")).read()))

        # 5. the event stream
        steps = run_state.read_steps()
        n0 = len(live.store.events(RUN))
        st = Stream(port, cookie)
        got_all = st.wait(n0, 5)
        t0 = time.time()
        post(INGEST_PORT, RUN, {"event": "process_completed", "utcTime": "2026-01-02T10:06:30Z",
                                "trace": {"task_id": 9, "name": f"CLINICAL:CLINVAR_SCREEN ({RUN})",
                                          "process": "CLINICAL:CLINVAR_SCREEN", "tag": RUN, "status": "COMPLETED",
                                          "exit": 0, "start": 1767348360000, "complete": 1767348390000}})
        live_ok = st.wait(n0 + 1, 2)
        check("the stream sends a posted process_completed within 2 s",
              got_all and live_ok and st.items[-1][1].get("event") == "process_completed" and st.items[-1][1].get("step") == "06",
              (got_all, live_ok, time.time() - t0, st.items[-1:]))
        want = pgp_live.messages(live.store.events(RUN), steps)
        check("what the stream sent equals messages()", st.items == want, (st.items[:3], want[:3]))
        st.close()
        st2 = Stream(port, cookie, last=n0 - 1)
        st2.wait(2, 5)
        check("Last-Event-ID resumes after the id given", [i for i, _ in st2.items] == [n0, n0 + 1],
              [i for i, _ in st2.items])
        st2.close()
        check("no message carries a path or a tag",
              not any(k in json.dumps(want) for k in ("/data/", "workDir", "outdir", '"tag"')), json.dumps(want)[:300])

        # 6. state.json
        _, _, body = request(port, "GET", f"/runs/{RUN}/state.json", headers=auth)
        s = json.loads(body)
        ref = run_state.collect(RUN, os.path.join(root, RUN), events_path=live.store.path(RUN))
        check("state.json's step states are run_state's",
              [(x["step"], x["state"]) for x in s["steps"]] == [(x["step"], x["state"]) for x in ref["steps"]],
              [(x["step"], x["state"]) for x in s["steps"] if x["state"] != "not planned"])
        mine = {x["step"]: x["state"] for x in s["steps"]}
        check("the fixture's steps: 06 done, 07 waiting, 11 done, 12 skipped",
              (mine["06"], mine["07"], mine["11"], mine["12"]) == ("done", "waiting", "done", "skipped"), mine)
        check("a state hook's answer is under plugins", s["plugins"].get("fixture", {}).get("steps") == len(s["steps"]),
              s.get("plugins"))
        import render_report
        tones = [pgp_live.headline("clinvar", {"state": "ok", "title": "ClinVar screen", "data": {"count": c, "hits": []}})
                 for c in (0, 3)]
        check("the ClinVar headline is never 'ok' and carries the report's caveat word for word",
              all(t[1] in ("info", "attention") for t in tones) and tones[1][1] == "attention"
              and all(render_report.CLINVAR_FREQ_NOTE in t[0] for t in tones), tones)
        check("a section with no file has no headline", pgp_live.headline("roh", {"state": "missing"}) is None)

        # 7. the policy
        code, _, _ = request(port, "GET", "/login")
        code, _, _ = request(port, "GET", "/assets/live.css")
        code, _, _ = request(port, "GET", "/nothing/here", headers=auth)
        bad = [p for p, hdr in RESPONSES if "default-src 'none'" not in hdr.get("content-security-policy", "")
               or "frame-ancestors 'none'" not in hdr.get("content-security-policy", "")
               or hdr.get("x-content-type-options") != "nosniff" or hdr.get("cache-control") != "no-store"]
        check(f"every browser response ({len(RESPONSES)}) carries the policy, nosniff and no-store", not bad, bad)
        texts = {f: open(os.path.join(ASSETS, f)).read() for f in os.listdir(ASSETS)
                 if os.path.isfile(os.path.join(ASSETS, f)) and not f.startswith(".")}
        check("assets/live holds no outside URL, inline script, on* attribute or eval", not asset_problems(texts),
              asset_problems(texts))
        planted = dict(texts)
        planted["index.html"] = texts["index.html"].replace("</body>", "<script>alert(1)</script></body>")
        control("an inline script in index.html", not asset_problems(planted))
        planted = dict(texts, **{"live.js": texts["live.js"] + "\nfetch('https://cdn.example.com/x.js');\n"})
        control("an outside URL in live.js", not asset_problems(planted))

        # 8. plugins
        code, _, body = request(port, "GET", f"/fixture/{RUN}", headers=auth)
        code2, _, _ = request(port, "GET", f"/fixture/{RUN}")
        check("a plugin route answers with its own flag behind the session, 401 without",
              code == 200 and json.loads(body) == {"word": "planted-flag", "run": RUN, "dir": True} and code2 == 401,
              (code, body, code2))
        code, _, body = request(port, "GET", f"/runs/{RUN}/extra/note.txt", headers=auth)
        check("a plugin's run file is served", code == 200 and body == b"plugin file\n", (code, body))
        app.routes = [(m, rx, fn, False) for m, rx, fn, _ in app.routes]
        code2, _, _ = request(port, "GET", f"/fixture/{RUN}")
        control("the plugin route registered with auth skipped", code2 == 401)
        app.routes = [(m, rx, fn, True) for m, rx, fn, _ in app.routes]
        adir = os.path.join(work, "assets")
        shutil.copytree(ASSETS, adir)
        for f in ("a.js", "b.js", "notes.txt"):
            open(os.path.join(adir, "plugins", f), "w").write("// plugin\n")
        live.cfg.assets = adir
        _, _, body = request(port, "GET", f"/runs/{RUN}/", headers=auth)
        tags = re.findall(r'<script src="/assets/plugins/([^"]+)"', body.decode())
        check("the run page gets one script tag per plugin file", tags == ["a.js", "b.js"], tags)
        code, _, _ = request(port, "GET", "/assets/plugins/a.js")
        check("a plugin script is served", code == 200, code)
        live.cfg.assets = ASSETS

        # 9. run files
        code, hd, body = request(port, "GET", f"/runs/{RUN}/report.html", headers=auth)
        check("report.html is the run's report, sandboxed", code == 200 and b"clinvar/hits.tsv" in body
              and "sandbox" in hd.get("content-security-policy", ""), (code, hd))
        code, _, body = request(port, "GET", f"/runs/{RUN}/clinvar/hits.tsv", headers=auth)
        check("a file the report links to is served", code == 200 and body.startswith(b"gene\t"), code)
        c1, _, _ = request(port, "GET", f"/runs/{RUN}/clinvar/unlinked.tsv", headers=auth)
        c2, _, _ = request(port, "GET", f"/runs/{RUN}/clinvar/escape.tsv", headers=auth)
        c3, _, _ = request(port, "GET", f"/runs/{RUN}/logs/run_status.tsv", headers=auth)
        c4, _, _ = request(port, "GET", f"/runs/{RUN}/..%2F..%2Foutside.txt", headers=auth)
        check("an unlinked file, a symlink out of the root and a path out of the run are 404",
              (c1, c2, c3, c4) == (404, 404, 404, 404), (c1, c2, c3, c4))
        code, _, body = request(port, "GET", "/", headers=auth)
        check("the index lists the run", code == 200 and f'/runs/{RUN}/'.encode() in body, code)
        log = open(os.path.join(state_dir, "access.log")).read()
        check("the access log has a row per page read, with no query string",
              "\treport\tDEMO-001\t" in log and "\tfile\tDEMO-001\tclinvar/hits.tsv\t200" in log and "?" not in log,
              log[-600:])

        # 10. the browser
        browser_checks(port, work)
    finally:
        stop(*servers)
        shutil.rmtree(work, ignore_errors=True)
    print()
    print("RESULT: ALL PASS" if not FAILS else f"RESULT: {len(FAILS)} FAILED")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
