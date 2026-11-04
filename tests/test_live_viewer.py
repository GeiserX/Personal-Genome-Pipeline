#!/usr/bin/env python3
"""The genome viewer on the live site: bin/owngenome_live_viewer.py, started by
scripts/owngenome-live.sh, on the synthetic sample of tests/test_viewer.py.

Checks:
  1. owngenome-live.sh downloads the igv.js pinned in versions.env (sha256
     checked) and starts the site with the viewer: /runs/<run>/viewer/ without
     a login sends the browser to /login, and the config and the reads answer
     401;
  2. after the login the page carries the viewer's Content-Security-Policy
     and says it is served by the site; the config names the run's BAM, VCF
     and SV VCF; a byte range of the BAM is 206 with exactly those bytes, past
     the end 416; a file the config does not name (a wrong name, a path out
     of data/) is 404; a run that does not exist is 404;
  3. the report the site serves links each position to /runs/<run>/viewer/
     with the site's note, not to the local viewer on 127.0.0.1:8765;
  4. opening the page writes one access-log row, and igv.js's byte ranges
     write none;
  5. a run whose BAM is a symlink out of GENOME_DIR gets no reads track, a
     run whose VCF index is one gets no VCF track and neither file is
     served, and a reference named with a space, '#' and '?' is served from
     the encoded URL its config gives;
  6. a site started without the viewer flags answers the viewer with 404 and
     the reason, and leaves the report's links as they are;
  7. in a headless Chrome (skipped without Chrome, except in CI) the page,
     loaded from the site, starts igv.js (status ready, 3 tracks) and the
     policy blocks nothing.
Checks 2, 3 and 5 are also run against planted mistakes (a viewer that serves
any file of the run folder, a report filter that does nothing, a viewer that
compares the still-quoted name, one that serves an index out of GENOME_DIR)
and must fail.

Run: python3 tests/test_live_viewer.py
"""
import http.client
import importlib.util
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import urllib.parse

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests"))
import owngenome_live  # noqa: E402
import render_report  # noqa: E402
import test_viewer as tv  # noqa: E402

PASSWORD = "viewer test password"
FAILED = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:600]}"))
    if not ok:
        FAILED.append(desc)


def request(port, path, headers=None, method="GET", body=None):
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=20)
    try:
        c.request(method, path, body=body, headers=dict({"Host": f"127.0.0.1:{port}"}, **(headers or {})))
        r = c.getresponse()
        return r.status, {k.lower(): v for k, v in r.getheaders()}, r.read()
    finally:
        c.close()


def login(port):
    code, h, _ = request(port, "/login", {"Content-Type": "application/x-www-form-urlencoded"}, "POST",
                         urllib.parse.urlencode({"password": PASSWORD, "next": "/"}))
    m = re.match(r"owngenome_live=([^;]+)", h.get("set-cookie", ""))
    return {"Cookie": f"owngenome_live={m.group(1)}"} if code == 303 and m else None


def wait_port(port, proc, seconds=120):
    for _ in range(seconds * 2):
        if proc.poll() is not None:
            return False
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=1):
                return True
        except OSError:
            time.sleep(0.5)
    return False


def in_process(argv, plugin_dir=None):
    """The site in this process (plugins from plugin_dir, default bin/), on free ports."""
    argv = argv + ["--port", "0", "--ingest-port", "0", "--quiet"] + (["--plugin-dir", plugin_dir] if plugin_dir else [])
    live, web, ingest = owngenome_live.start(argv)
    for s in (web, ingest):
        threading.Thread(target=s.serve_forever, daemon=True).start()
    return live, web, ingest


def stop(live, web, ingest):
    live.stopping = True
    for s in (web, ingest):
        s.shutdown()
        s.server_close()


def viewer_checks(port, auth, s, label=""):
    """Checks 2: what the viewer serves to a signed-in browser. Returns the config."""
    import json
    code, h, body = request(port, "/runs/S/viewer/", auth)
    csp = h.get("content-security-policy", "")
    check(f"{label}page: 200 with the viewer's policy, served by the site",
          code == 200 and csp == tv.serve.CSP and b"served by the run site" in body
          and b"served from this computer only" not in body, (code, csp, body[:300]))
    code, _, body = request(port, "/runs/S/viewer/igv.min.js", auth)
    c2, _, js = request(port, "/runs/S/viewer/viewer.js", auth)
    check(f"{label}page: igv.min.js and viewer.js are served beside it",
          code == 200 and len(body) > 100000 and c2 == 200 and b"config.json" in js, (code, c2))
    code, _, body = request(port, "/runs/S/viewer/config.json", auth)
    cfg = json.loads(body) if code == 200 else {}
    names = [t["url"].rsplit("/", 1)[1] for t in cfg.get("tracks", [])]
    check(f"{label}config: the run's BAM, VCF and SV VCF, by relative URL",
          names == ["S_sorted.bam", "S.vcf.gz", "S_sv_consensus.vcf.gz"]
          and all(t["url"].startswith("data/") for t in cfg.get("tracks", [])), (code, cfg))
    raw = open(f"{s}/aligned/S_sorted.bam", "rb").read()
    bam = f"/runs/S/viewer/{cfg['tracks'][0]['url']}" if cfg.get("tracks") else "/runs/S/viewer/data/x"
    code, h, body = request(port, bam, dict(auth, Range="bytes=10-109"))
    check(f"{label}reads: a BAM range is 206 with exactly those bytes",
          code == 206 and body == raw[10:110] and h.get("content-range") == f"bytes 10-109/{len(raw)}", (code, h))
    code, h, _ = request(port, bam, dict(auth, Range=f"bytes={len(raw) + 5}-"))
    check(f"{label}reads: a range past the end is 416", code == 416, code)
    bad = [request(port, p, auth)[0] for p in (
        "/runs/S/viewer/data/2/S_report.html", "/runs/S/viewer/data/99/S_sorted.bam",
        "/runs/S/viewer/S_report.html", "/runs/S/viewer/data/2/..%2F..%2FS_report.html",
        "/runs/S/viewer/serve.py", "/runs/NOPE/viewer/config.json")]
    check(f"{label}a file the config does not name, and a run that does not exist, are 404",
          bad == [404] * len(bad), bad)
    return cfg


def report_checks(port, auth, label=""):
    """Checks 3: the served report's links point at the site's viewer."""
    code, _, body = request(port, "/runs/S/report.html", auth)
    page = body.decode(errors="replace")
    check(f"{label}report: positions link to /runs/S/viewer/ with the site's note, none to 127.0.0.1:8765",
          code == 200 and 'href="/runs/S/viewer/#chr1:451"' in page and render_report.VIEWER_URL not in page
          and "this site's genome viewer" in page
          and render_report.VIEWER_NOTE not in page, (code, [x for x in page.splitlines() if "viewer" in x][:4]))


def outside_checks(port, auth, label=""):
    """Checks 5: files out of GENOME_DIR are left out; an odd reference name is served."""
    import json
    code, _, body = request(port, "/runs/T/viewer/config.json", auth)
    tracks = json.loads(body).get("tracks") if code == 200 else None
    check(f"{label}a BAM symlinked out of GENOME_DIR is left out of the config", code == 200 and tracks == [],
          (code, body[:300]))
    code, _, body = request(port, "/runs/U/viewer/config.json", auth)
    cfg = json.loads(body) if code == 200 else {}
    c2, _, _ = request(port, "/runs/U/viewer/data/2/U.vcf.gz", dict(auth, Range="bytes=0-9"))
    c3, _, _ = request(port, "/runs/U/viewer/data/3/U.vcf.gz.tbi", dict(auth, Range="bytes=0-9"))
    check(f"{label}a VCF whose index is symlinked out of GENOME_DIR: no track, neither file served",
          code == 200 and cfg.get("tracks") == [] and c2 == 404 and c3 == 404, (code, cfg, c2, c3))
    url = cfg.get("reference", {}).get("fastaURL", "")
    code, _, body = request(port, f"/runs/U/viewer/{url}", dict(auth, Range="bytes=0-4"))
    check(f"{label}a reference named with a space, '#' or '?' is served from its encoded URL",
          url == "data/0/my%20ref%231%3F.fasta" and code == 206 and body == b">chr1", (url, code, body))


def mutant_dir(work, name, old, new):
    """A plugin folder holding bin/owngenome_live_viewer.py with one planted mistake (beside a link
    to assets/, which the plugin finds from its own folder)."""
    top = os.path.join(work, f"mutant-{name}")
    d = os.path.join(top, "bin")
    os.makedirs(d)
    os.symlink(os.path.join(REPO, "assets"), os.path.join(top, "assets"))
    src = open(os.path.join(REPO, "bin", "owngenome_live_viewer.py")).read()
    if old not in src:
        check(f"mutant '{name}': its line is in owngenome_live_viewer.py", False, old)
    with open(os.path.join(d, "owngenome_live_viewer.py"), "w") as f:
        f.write(src.replace(old, new))
    return d


def main():
    work = tempfile.mkdtemp(prefix="pgp-live-viewer-")
    proc = None
    try:
        gd = f"{work}/genome"
        tv.make_sample(gd)
        s = f"{gd}/S"
        os.makedirs(f"{s}/nextflow")  # the site lists a folder with nextflow/ as a run
        html, _ = tv.findings_report(work)
        with open(f"{s}/S_report.html", "w") as f:
            f.write(html)
        state = f"{work}/state"
        os.makedirs(state)
        with open(f"{state}/password", "w") as f:
            f.write(owngenome_live.hash_password(PASSWORD) + "\n")

        # 1. owngenome-live.sh, with a password
        port, iport = tv.free_port(), tv.free_port()
        env = dict(os.environ, GENOME_DIR=gd, REF_FASTA="reference/synthetic.fasta", OWNGENOME_LIVE_PORT=str(port),
                   OWNGENOME_LIVE_INGEST_PORT=str(iport), OWNGENOME_LIVE_STATE_DIR=state)
        log = open(f"{work}/site.log", "w")
        proc = subprocess.Popen(["bash", f"{REPO}/scripts/owngenome-live.sh", "--quiet"], env=env,
                                stdout=log, stderr=subprocess.STDOUT)
        up = wait_port(port, proc)
        check("owngenome-live.sh: the site starts", up, open(f"{work}/site.log").read())
        if not up:
            raise SystemExit(1)
        pin = tv.igv_pin()
        igv = f"{gd}/tools/igv.js-{pin['IGV_JS_VERSION']}/igv.min.js"
        import hashlib
        check("owngenome-live.sh: the igv.js pinned in versions.env is downloaded, with its sha256",
              os.path.isfile(igv) and hashlib.sha256(open(igv, "rb").read()).hexdigest() == pin["IGV_JS_SHA256"],
              open(f"{work}/site.log").read())
        code, h, _ = request(port, "/runs/S/viewer/")
        c2, _, _ = request(port, "/runs/S/viewer/config.json")
        c3, _, _ = request(port, "/runs/S/viewer/data/2/S_sorted.bam", {"Range": "bytes=0-9"})
        check("no login: the page sends the browser to /login, the config and the reads are 401",
              code == 303 and h.get("location") == "/login?next=/runs/S/viewer/" and c2 == 401 and c3 == 401,
              (code, h.get("location"), c2, c3))
        auth = login(port)
        check("login: the password opens a session", auth is not None)
        auth = auth or {}

        # 2 and 3. the viewer and the served report
        acc = f"{state}/access.log"
        before = open(acc).read().count("\tviewer\t") if os.path.isfile(acc) else 0
        cfg = viewer_checks(port, auth, s)
        report_checks(port, auth)
        code, h, _ = request(port, "/runs/S/viewer", auth)
        check("/runs/S/viewer redirects to /runs/S/viewer/", code == 301 and h.get("location") == "/runs/S/viewer/",
              (code, h.get("location")))

        # 4. the access log
        time.sleep(0.5)
        rows = [r for r in open(acc).read().splitlines() if "\tviewer\t" in r]
        mine = rows[before:]
        check("access log: one row for the page, none for config.json, igv.js or the byte ranges",
              len(mine) == 2 and all("\t200" in r or "\t404" in r for r in mine), mine)
        proc.terminate()
        proc.wait(timeout=10)
        proc = None

        # 5. a BAM, and a VCF index, symlinked out of GENOME_DIR; an odd reference name
        t = f"{gd}/T"
        os.makedirs(f"{t}/aligned")
        os.makedirs(f"{t}/nextflow")
        shutil.copy(f"{s}/aligned/S_sorted.bam", f"{work}/outside.bam")
        shutil.copy(f"{s}/aligned/S_sorted.bam.bai", f"{work}/outside.bam.bai")
        os.symlink(f"{work}/outside.bam", f"{t}/aligned/T_sorted.bam")
        os.symlink(f"{work}/outside.bam.bai", f"{t}/aligned/T_sorted.bam.bai")
        u = f"{gd}/U"
        os.makedirs(f"{u}/vcf")
        os.makedirs(f"{u}/nextflow")
        shutil.copy(f"{s}/vcf/S.vcf.gz", f"{u}/vcf/U.vcf.gz")
        shutil.copy(f"{s}/vcf/S.vcf.gz.tbi", f"{work}/outside.vcf.gz.tbi")
        os.symlink(f"{work}/outside.vcf.gz.tbi", f"{u}/vcf/U.vcf.gz.tbi")
        odd = f"{gd}/reference/my ref#1?.fasta"
        shutil.copy(f"{gd}/reference/synthetic.fasta", odd)
        shutil.copy(f"{gd}/reference/synthetic.fasta.fai", odd + ".fai")
        flags = ["--root", gd, "--state-dir", f"{work}/state2", "--viewer-fasta", odd, "--viewer-igv-js", igv]
        srv = in_process(flags)
        try:
            p = srv[1].server_address[1]
            outside_checks(p, {})

            # 7. headless Chrome, on the site without a password (loopback)
            chrome = tv.find_chrome()
            if not chrome:
                ci = os.environ.get("CI") == "true"
                print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the page in")
                if ci:
                    FAILED.append("no Chrome in CI")
            else:
                r = subprocess.run([chrome, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-first-run",
                                    "--no-default-browser-check", f"--user-data-dir={work}/chrome",
                                    "--virtual-time-budget=20000", "--dump-dom",
                                    f"http://127.0.0.1:{p}/runs/S/viewer/#chr1:451"],
                                   capture_output=True, text=True, timeout=180)
                status = re.search(r'data-pgp-status="([^"]*)"', r.stdout)
                tracks = re.search(r'data-pgp-tracks="([^"]*)"', r.stdout)
                blocked = re.search(r'data-pgp-blocked="([^"]*)"', r.stdout)
                check("browser: the page from the site starts igv.js with the run's 3 tracks",
                      status and status.group(1) == "ready" and tracks and tracks.group(1) == "3",
                      (status and status.group(1), r.stderr[-800:], r.stdout[:600]))
                check("browser: the policy blocked nothing the page needed", not blocked, blocked and blocked.group(1))
        finally:
            stop(*srv)

        # 6. a site without the viewer flags
        srv = in_process(["--root", gd, "--state-dir", f"{work}/state3"])
        try:
            p = srv[1].server_address[1]
            code, _, body = request(p, "/runs/S/viewer/")
            _, _, rep = request(p, "/runs/S/report.html")
            check("without the flags: the viewer is 404 with the reason, the report keeps its local links",
                  code == 404 and b"not set up" in body and render_report.VIEWER_URL.encode() in rep, (code, body))
        finally:
            stop(*srv)

        # Planted mistakes: checks 2, 3 and 5 must fail on them.
        for name, old, new in (
                ("serves any file", "        n = int(dm.group(1)) if dm else len(files)\n",
                 "        if rest and os.path.isfile(os.path.join(d, rest)):\n"
                 "            return serve.send_file(req, head_writer(req), 'application/octet-stream', "
                 "os.path.join(d, rest))\n        n = int(dm.group(1)) if dm else len(files)\n"),
                ("no report filter", "    app.report_filter(point_links)\n", ""),
                ("quoted name", "== urllib.parse.unquote(dm.group(2)):", "== dm.group(2):"),
                ("index out served", "    files = [f if n < 2 or inside(f) else None for n, f in enumerate(files)]\n",
                 "")):
            pd = mutant_dir(work, name.replace(" ", "-"), old, new)
            srv = in_process(flags[:3] + [f"{work}/state-{name.replace(' ', '-')}"] + flags[4:], pd)
            try:
                p = srv[1].server_address[1]
                quiet = []
                real, globals()["check"] = check, (lambda desc, ok, detail="": quiet.append(ok))
                try:
                    viewer_checks(p, {}, s)
                    report_checks(p, {})
                    outside_checks(p, {})
                finally:
                    globals()["check"] = real
                check(f"control: the checks fail on the planted mistake '{name}'", not all(quiet), quiet)
            finally:
                stop(*srv)
    finally:
        if proc and proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
        shutil.rmtree(work, ignore_errors=True)

    print()
    print("RESULT: ALL PASS" if not FAILED else f"RESULT: {len(FAILED)} FAILED")
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
