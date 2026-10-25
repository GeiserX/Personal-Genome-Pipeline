#!/usr/bin/env python3
"""The local genome viewer: scripts/pgp-view.sh and assets/viewer/serve.py,
and the report's links to it (bin/render_report.py), on a synthetic sample.

The sample is made here: a random 3 kb chr1 and 1 kb chr2 reference with its
.fai, 30 reads on chr1 as a coordinate-sorted BAM with its .bai, a bgzipped
small-variant VCF with its .tbi and an SV VCF without an index, each written
byte by byte (no samtools or htslib needed).

Checks:
  1. build_config names the sample's files and the locus, and refuses a locus
     that is not chr[:pos[-end]] of the reference, and a BAM with no index;
  2. pgp-view.sh downloads the pinned igv.js and checks its sha256, finds the
     BAM, VCF and SV VCF in the sample folder, and starts the server, which
     listens on 127.0.0.1 and not on the machine's other address;
  3. the server sends the page with a Content-Security-Policy that keeps it
     on this server, the config with the files and the locus, a byte range of
     the BAM as 206 with exactly those bytes, 416 for a range past the end,
     404 for a path it does not serve, 403 for a foreign Host header;
  4. in a headless Chrome the page starts igv.js on that config (status
     "ready", nothing blocked by the policy) and igv.js reads the BAM and its
     index by range; skipped without Chrome, except in CI;
  5. the HTML report links each ClinVar hit, each mitochondrial variant and
     each secondary-findings variant to the viewer at its position; a file
     name with a space, '#' or '?' is served from its percent-encoded URL.
"""

import gzip
import json
import os
import random
import re
import shutil
import socket
import struct
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
import zlib

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "assets", "viewer"))
sys.path.insert(0, os.path.join(REPO, "bin"))
import serve  # noqa: E402

FAILED = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:600]}"))
    if not ok:
        FAILED.append(desc)


# --- synthetic files ----------------------------------------------------------
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzf_block(data):
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = c.compress(data) + c.flush()
    head = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff" + struct.pack("<H", 6) + b"BC" + struct.pack("<HH", 2, len(body) + 25)
    return head + body + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))


def reg2bin(beg, end):
    end -= 1
    for shift, off in ((14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)):
        if beg >> shift == end >> shift:
            return off + (beg >> shift)
    return 0


def binning_index(n_ref, recs, rec_voff, end_voff):
    """BAI/TBI per-reference part: recs is (ref, beg, end) of every record,
    all in one BGZF block that starts at rec_voff."""
    out = b""
    for r in range(n_ref):
        mine = [x for x in recs if x[0] == r]
        if not mine:
            out += struct.pack("<ii", 0, 0)
            continue
        bins = sorted({reg2bin(b, e) for _, b, e in mine})
        out += struct.pack("<i", len(bins))
        for b in bins:
            out += struct.pack("<Iiqq", b, 1, rec_voff, end_voff)
        last = max(e for _, _, e in mine)
        n_intv = (last - 1 >> 14) + 1
        out += struct.pack("<i", n_intv) + struct.pack("<q", rec_voff) * n_intv
    return out


def write_bam(path, contigs, reads):
    """reads: (ref index, 0-based pos, sequence), sorted; 50M each."""
    text = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{n}\tLN:{len(s)}\n" for n, s in contigs) + "@RG\tID:S\tSM:S\n"
    head = b"BAM\x01" + struct.pack("<i", len(text)) + text.encode() + struct.pack("<i", len(contigs))
    for n, s in contigs:
        head += struct.pack("<i", len(n) + 1) + n.encode() + b"\0" + struct.pack("<i", len(s))
    codes = {c: i for i, c in enumerate("=ACMGRSVTWYHKDBN")}
    body, recs = b"", []
    for i, (ref, pos, seq) in enumerate(reads):
        name = f"r{i}".encode() + b"\0"
        packed = bytes((codes[seq[j]] << 4) | (codes[seq[j + 1]] if j + 1 < len(seq) else 0) for j in range(0, len(seq), 2))
        rec = struct.pack("<iiBBHHHiiii", ref, pos, len(name), 60, reg2bin(pos, pos + len(seq)), 1, 0, len(seq), -1, -1, 0)
        rec += name + struct.pack("<I", len(seq) << 4) + packed + bytes([30] * len(seq)) + b"RGZS\0"
        body += struct.pack("<i", len(rec)) + rec
        recs.append((ref, pos, pos + len(seq)))
    b0, b1 = bgzf_block(head), bgzf_block(body)
    with open(path, "wb") as f:
        f.write(b0 + b1 + BGZF_EOF)
    with open(path + ".bai", "wb") as f:
        f.write(b"BAI\x01" + struct.pack("<i", len(contigs)) + binning_index(len(contigs), recs, len(b0) << 16, (len(b0) + len(b1)) << 16))


def write_vcf(path, contigs, records, index):
    """records: (chrom, pos 1-based, ref, alt, info); bgzipped, with a .tbi."""
    head = "##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={n},length={len(s)}>\n" for n, s in contigs)
    head += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
    head += '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
    head += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    head += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    body = "".join(f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t{i}\tGT\t0/1\n" for c, p, r, a, i in records)
    b0, b1 = bgzf_block(head.encode()), bgzf_block(body.encode())
    with open(path, "wb") as f:
        f.write(b0 + b1 + BGZF_EOF)
    if index:
        names = sorted({c for c, *_ in records}, key=[n for n, _ in contigs].index)
        nm = b"".join(n.encode() + b"\0" for n in names)
        recs = [(names.index(c), p - 1, p - 1 + len(r)) for c, p, r, _, _ in records]
        tbi = b"TBI\x01" + struct.pack("<iiiiiii", len(names), 2, 1, 2, 0, ord("#"), 0) + struct.pack("<i", len(nm)) + nm
        tbi += binning_index(len(names), recs, len(b0) << 16, (len(b0) + len(b1)) << 16)
        with open(path + ".tbi", "wb") as f:
            f.write(bgzf_block(tbi) + BGZF_EOF)


def make_sample(gd):
    rnd = random.Random(7)
    contigs = [("chr1", "".join(rnd.choice("ACGT") for _ in range(3000))),
               ("chr2", "".join(rnd.choice("ACGT") for _ in range(1000)))]
    os.makedirs(f"{gd}/reference")
    fa = f"{gd}/reference/synthetic.fasta"
    with open(fa, "w") as f, open(fa + ".fai", "w") as fai:
        off = 0
        for n, s in contigs:
            off += len(n) + 2
            fai.write(f"{n}\t{len(s)}\t{off}\t60\t61\n")
            lines = [s[i:i + 60] for i in range(0, len(s), 60)]
            f.write(f">{n}\n" + "\n".join(lines) + "\n")
            off += len(s) + len(lines)
    s = f"{gd}/S"
    for d in ("aligned", "vcf", "sv_merged"):
        os.makedirs(f"{s}/{d}")
    reads = [(0, p, contigs[0][1][p:p + 50]) for p in range(400, 700, 10)]
    write_bam(f"{s}/aligned/S_sorted.bam", contigs, reads)
    seq = contigs[0][1]
    write_vcf(f"{s}/vcf/S.vcf.gz", contigs, [("chr1", 451, seq[450], "ACGT".replace(seq[450], "")[0], ".")], True)
    write_vcf(f"{s}/sv_merged/S_sv_consensus.vcf.gz", contigs, [("chr1", 1000, seq[999], "<DEL>", "SVTYPE=DEL;END=1500")], False)
    return fa, contigs


# --- helpers ------------------------------------------------------------------
def free_port():
    with socket.socket() as so:
        so.bind(("127.0.0.1", 0))
        return so.getsockname()[1]


def get(url, headers=None):
    req = urllib.request.Request(url, headers=headers or {})
    try:
        with urllib.request.urlopen(req, timeout=10) as r:
            return r.status, dict(r.headers), r.read()
    except urllib.error.HTTPError as e:
        return e.code, dict(e.headers), e.read()


def other_address():
    """An address of this machine that is not loopback, or None."""
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as so:
            so.connect(("192.0.2.1", 9))
            ip = so.getsockname()[0]
        return None if ip.startswith("127.") else ip
    except OSError:
        return None


def find_chrome():
    for c in (os.environ.get("CHROME"), "google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
              "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"):
        if c and (shutil.which(c) or os.path.isfile(c)):
            return shutil.which(c) or c
    return None


def main():
    work = tempfile.mkdtemp(prefix="pgp-view-")
    proc = None
    try:
        gd = f"{work}/genome"
        fa, contigs = make_sample(gd)
        s = f"{gd}/S"

        # 1. build_config
        cfg, files = serve.build_config("S", fa, f"{s}/aligned/S_sorted.bam", f"{s}/vcf/S.vcf.gz",
                                        f"{s}/sv_merged/S_sv_consensus.vcf.gz", "chr1:451")
        check("config: the bind address is 127.0.0.1", serve.BIND == "127.0.0.1")
        check("config: the locus, the reference and its .fai",
              cfg["locus"] == "chr1:451" and files[0] == fa and files[1] == fa + ".fai"
              and cfg["reference"]["fastaURL"] == "data/0/synthetic.fasta", cfg)
        kinds = [(t["type"], t["format"], os.path.basename(files[int(t["url"].split("/")[1])])) for t in cfg["tracks"]]
        check("config: the BAM, the VCF and the SV VCF of the sample",
              kinds == [("alignment", "bam", "S_sorted.bam"), ("variant", "vcf", "S.vcf.gz"),
                        ("variant", "vcf", "S_sv_consensus.vcf.gz")], kinds)
        check("config: the BAM's .bai and the VCF's .tbi; the SV VCF without one loads whole",
              files[int(cfg["tracks"][0]["indexURL"].split("/")[1])].endswith("S_sorted.bam.bai")
              and cfg["tracks"][1]["indexURL"].endswith("S.vcf.gz.tbi") and cfg["tracks"][2].get("indexed") is False, cfg)
        check("config: no locus opens the first contig", serve.build_config("S", fa)[0]["locus"] == "chr1")
        for bad in ("chr9:100", "chr1:99999", "chr1;alert(1)", "BRCA1 x"):
            try:
                serve.build_config("S", fa, locus=bad)
                check(f"config: locus '{bad}' is refused", False)
            except ValueError:
                check(f"config: locus '{bad}' is refused", True)
        os.rename(f"{s}/aligned/S_sorted.bam.bai", f"{work}/bai")
        try:
            serve.build_config("S", fa, f"{s}/aligned/S_sorted.bam")
            check("config: a BAM without an index is refused", False)
        except ValueError as e:
            check("config: a BAM without an index is refused", "samtools index" in str(e), e)
        os.rename(f"{work}/bai", f"{s}/aligned/S_sorted.bam.bai")

        # 2. pgp-view.sh
        port = free_port()
        env = dict(os.environ, GENOME_DIR=gd, REF_FASTA="reference/synthetic.fasta", PGP_VIEW_PORT=str(port))
        log = open(f"{work}/server.log", "w")
        proc = subprocess.Popen(["bash", f"{REPO}/scripts/pgp-view.sh", "S", "chr1:451"], env=env,
                                stdout=log, stderr=subprocess.STDOUT)
        base = f"http://127.0.0.1:{port}"
        up = False
        for _ in range(240):
            if proc.poll() is not None:
                break
            try:
                with socket.create_connection(("127.0.0.1", port), timeout=1):
                    up = True
                    break
            except OSError:
                time.sleep(0.5)
        out = open(f"{work}/server.log").read()
        check("pgp-view.sh: the server starts and prints its address", up and f"{base}/#chr1:451" in out, out)
        if not up:
            raise SystemExit(1)
        igv = f"{gd}/tools/igv.js-{re.search(r'IGV_JS_VERSION=(.+)', open(f'{REPO}/scripts/pgp-view.sh').read()).group(1)}/igv.min.js"
        check("pgp-view.sh: the pinned igv.js is downloaded into GENOME_DIR/tools", os.path.getsize(igv) > 100000, igv)
        ip = other_address()
        if ip:
            try:
                socket.create_connection((ip, port), timeout=2).close()
                check(f"server: not reachable on {ip}, only on 127.0.0.1", False)
            except OSError:
                check("server: not reachable on the machine's other address, only on 127.0.0.1", True)
        else:
            print("[SKIP] no non-loopback address to test the bind against")

        # 3. what it serves
        code, h, body = get(f"{base}/")
        csp = h.get("Content-Security-Policy", "")
        check("server: the page, with a policy that keeps every request on this server",
              code == 200 and b"igv.min.js" in body and "connect-src 'self'" in csp and "default-src 'none'" in csp, (code, csp))
        code, _, body = get(f"{base}/config.json")
        sc = json.loads(body) if code == 200 else {}
        check("server: the config names the sample's files and the locus",
              sc.get("locus") == "chr1:451" and [t["url"].rsplit("/", 1)[1] for t in sc.get("tracks", [])]
              == ["S_sorted.bam", "S.vcf.gz", "S_sv_consensus.vcf.gz"], sc)
        bam_url = f"{base}/{sc['tracks'][0]['url']}"
        raw = open(f"{s}/aligned/S_sorted.bam", "rb").read()
        code, h, body = get(bam_url, {"Range": "bytes=10-109"})
        check("server: a BAM range comes back as 206 with exactly those bytes",
              code == 206 and body == raw[10:110] and h.get("Content-Range") == f"bytes 10-109/{len(raw)}", (code, h))
        code, h, body = get(bam_url, {"Range": "bytes=-28"})
        check("server: a suffix range (the BGZF EOF block)", code == 206 and body == BGZF_EOF, code)
        code, h, _ = get(bam_url, {"Range": f"bytes={len(raw) + 5}-"})
        check("server: a range past the end is 416", code == 416 and h.get("Content-Range") == f"bytes */{len(raw)}", code)
        code, _, _ = get(f"{base}/data/0/../../reference/synthetic.fasta.fai")
        code2, _, _ = get(f"{base}/serve.py")
        check("server: a path it does not serve is 404", code == 404 and code2 == 404, (code, code2))
        code, _, _ = get(f"{base}/config.json", {"Host": "attacker.example:8765"})
        check("server: a foreign Host header is 403", code == 403, code)
        # a file name with a space, '#' and '?' (a REF_FASTA the user named so)
        odd = f"{work}/odd"
        os.makedirs(odd)
        shutil.copy(fa, f"{odd}/my ref#1?.fasta")
        shutil.copy(fa + ".fai", f"{odd}/my ref#1?.fasta.fai")
        ocfg, ofiles = serve.build_config("S", f"{odd}/my ref#1?.fasta")
        osrv = serve.make_server(ocfg, ofiles, igv, 0, quiet=True)
        threading.Thread(target=osrv.serve_forever, daemon=True).start()
        try:
            obase = f"http://127.0.0.1:{osrv.server_address[1]}"
            code, _, body = get(f"{obase}/{ocfg['reference']['fastaURL']}", {"Range": "bytes=0-4"})
            check("server: a file name with a space, '#' or '?' is served from its encoded URL",
                  ocfg["reference"]["fastaURL"] == "data/0/my%20ref%231%3F.fasta" and code == 206 and body == b">chr1",
                  (ocfg["reference"]["fastaURL"], code, body))
        finally:
            osrv.shutdown()
            osrv.server_close()

        # 4. the page in a headless browser
        chrome = find_chrome()
        if not chrome:
            ci = os.environ.get("CI") == "true"
            print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the page in")
            if ci:
                FAILED.append("no Chrome in CI")
        else:
            before = len(open(f"{work}/server.log").read())
            r = subprocess.run([chrome, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-first-run",
                                "--no-default-browser-check", f"--user-data-dir={work}/chrome",
                                "--virtual-time-budget=20000", "--dump-dom", f"{base}/#chr1:451"],
                               capture_output=True, text=True, timeout=180)
            dom = r.stdout
            status = re.search(r'data-pgp-status="([^"]*)"', dom)
            blocked = re.search(r'data-pgp-blocked="([^"]*)"', dom)
            check("browser: igv.js starts on the config (status ready)", status and status.group(1) == "ready",
                  (status and status.group(1), r.stderr[-800:], dom[:800]))
            check("browser: the policy blocked nothing the page needed", not blocked, blocked and blocked.group(1))
            served = open(f"{work}/server.log").read()[before:]
            check("browser: igv.js read the BAM's index and a range of the BAM",
                  "S_sorted.bam.bai" in served and re.search(r'"GET /data/\d+/S_sorted\.bam HTTP', served), served[-1500:])
            check("browser: every request of the page was served (no 4xx)",
                  not re.search(r'HTTP/1\.[01]" 4\d\d', served), served[-1500:])
            check("browser: igv.js read the reference and the VCF",
                  "synthetic.fasta" in served and "S.vcf.gz" in served, served[-1500:])

        # 5. the report's links
        import collect_summary
        import render_report
        d = f"{work}/report/S"
        os.makedirs(f"{d}/clinvar")
        with open(f"{d}/clinvar/S_clinvar_hits.vcf", "w") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                    "chr1\t451\t2\tT\tG\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Pathogenic;"
                    "CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n")
        html = render_report.html_report(collect_summary.collect("S", d))
        check("report: the ClinVar hit links to the viewer at its position",
              'href="http://127.0.0.1:8765/#chr1:451"' in html and "scripts/pgp-view.sh" in html,
              [l for l in html.splitlines() if "ClinVar" in l][:5])
        summ = collect_summary.collect("S", d)
        for k in ("mito", "clinical"):
            summ["sections"][k]["state"] = "ok"
        summ["sections"]["mito"]["data"]["annotation"] = {
            "source": "synthetic", "variants": [{"pos": 3243, "ref": "A", "alt": "G"}]}
        summ["sections"]["clinical"]["data"]["acmg_sf"] = {
            "version": "v3.3",
            "clinvar_hits": [{"gene": "GENEC", "variant": "chr2:120 A>G", "genotype": "0/1", "significance": "Pathogenic"}],
            "high_impact": [{"gene": "GENED", "variant": "chr1:2000 C>T", "genotype": "0/1", "consequence": "stop_gained"}]}
        html = render_report.html_report(summ)
        want = [f'href="http://127.0.0.1:8765/#{x}"' for x in ("chrM:3244", "chr2:120", "chr1:2000")]
        check("report: the mito and secondary-findings variants link to the viewer at their positions",
              all(w in html for w in want), [l for l in html.splitlines() if "8765" in l][:5])
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
