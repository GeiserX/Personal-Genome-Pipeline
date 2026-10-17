#!/usr/bin/env python3
"""bin/run_delta.py on two synthetic runs of the same sample.

Checks:
  1. no previous run: the delta says so, lists no change, and still writes
     the text, HTML and JSON files;
  2. two identical runs (only the time each was written differs): no change;
  3. a planted ClinVar reclassification and a changed image digest: the
     delta names exactly those two changes, the ClinVar one with the digest
     change as a possible cause only when that image feeds ClinVar, and with
     "the call itself differs" when it does not;
  4. a hit added, one dropped, a re-starred one, a CPIC call, pypgx's CYP2D6,
     a PRS percentile that moved past the band and one that stayed inside it,
     a CPSR count, a section that stopped running, a new ClinVar file date
     and a new versions.env tag: each is listed once, and the reclassified
     hit names the ClinVar file date as the reason;
  5. the HTML escapes what it prints.

Run: python3 tests/test_run_delta.py
"""
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(REPO, "bin", "run_delta.py")

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1


def hit(pos, sig, stars=2, rev="criteria provided, multiple submitters, no conflicts", gt="het", gene="GENEA"):
    return {"stars": stars, "chrom": "chr1", "pos": pos, "ref": "A", "alt": "G", "genotype": gt,
            "gene": gene, "significance": sig, "review_status": rev}


def sec(title, state="ok", data=None):
    return {"title": title, "step": "x", "state": state, "source": None, "file_date": None,
            "note": None, "data": data or {}}


def summary(generated="2026-01-01T00:00:00Z"):
    return {
        "schema_version": 1, "sample": "S", "generated_utc": generated,
        "sections": {
            "clinvar": sec("ClinVar screen", data={"count": 2, "hits": [hit(100, "Pathogenic"),
                                                                      hit(200, "Likely pathogenic")]}),
            "cpsr": sec("Cancer predisposition (CPSR)", data={"classification": {"Benign": 10, "VUS": 3}}),
            "cpic": sec("CPIC drug-gene recommendations", data={"genes": [
                {"gene": "CYP2C19", "diplotype": "*1/*1", "phenotype": "Normal Metabolizer", "status": "normal"}]}),
            "pypgx": sec("pypgx", data={"cyp2d6": "*1/*2"}),
            "prs": sec("Polygenic risk scores", data={"scores": [
                {"pgs_id": "PGS000001", "condition": "Trait one", "percentile": "40.0"},
                {"pgs_id": "PGS000002", "condition": "Trait two", "percentile": "50.0"}]}),
            "roh": sec("Runs of homozygosity"),
        },
        "cyp2d6": {"consensus": None},
    }


MANIFEST = [
    ("run", "written_utc", "2026-01-01 00:00:00 UTC"),
    ("run", "written_by", "run-all.sh"),
    ("run", "git_commit", "abc"),
    ("version", "PYTHON_IMAGE", "python-ref"),
    ("version", "DEEPVARIANT_IMAGE", "deepvariant-ref"),
    ("version", "PHARMCAT_IMAGE", "pharmcat-ref"),
    ("version", "VEP_CACHE_RELEASE", "115"),
    ("image", "python-ref", "digest-aaa"),
    ("image", "deepvariant-ref", "digest-bbb"),
    ("image", "pharmcat-ref", "digest-ccc"),
    ("data", "clinvar_file_date", "20260101 (clinvar/clinvar.vcf.gz)"),
    ("data", "ancestry_panel", "not installed (setup.sh --ancestry-panel)"),
]


def write_run(d, summ, manifest):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "summary.json"), "w") as f:
        json.dump(summ, f)
    with open(os.path.join(d, "run_manifest.tsv"), "w") as f:
        for r in manifest:
            f.write("\t".join(r) + "\n")


def set_row(manifest, sec_, key, value):
    return [(s, k, value if (s, k) == (sec_, key) else v) for s, k, v in manifest]


def delta(w, prev_dir, cur_dir, name):
    out = os.path.join(w, name)
    args = [sys.executable, TOOL, "--sample", "S", "--summary", os.path.join(cur_dir, "summary.json"),
            "--manifest", os.path.join(cur_dir, "run_manifest.tsv"), "--out", out]
    if prev_dir:
        args += ["--previous", prev_dir]
    p = subprocess.run(args, capture_output=True, text=True)
    check(f"{name}: exits 0", p.returncode == 0, p.stderr)
    with open(out + ".json") as f:
        d = json.load(f)
    with open(out + ".txt") as f:
        txt = f.read()
    with open(out + ".html") as f:
        page = f.read()
    return d, txt, page


def main():
    w = tempfile.mkdtemp()
    try:
        prev, cur = os.path.join(w, "history", "2026-01-01T000000Z"), os.path.join(w, "cur")
        write_run(prev, summary(), MANIFEST)

        # 1. no previous run
        write_run(cur, summary("2026-02-01T00:00:00Z"), set_row(MANIFEST, "run", "written_utc", "later"))
        d, txt, page = delta(w, None, cur, "first")
        check("first run: previous is null and no change", d["previous"] is None and d["changes"] == [], d)
        check("first run: the text says there is no previous run", "No previous run of this sample" in txt, txt)
        check("first run: the HTML says so too", "No previous run of this sample" in page, page)

        # 2. identical runs
        d, txt, _ = delta(w, prev, cur, "same")
        check("identical runs: no change", d["changes"] == [], d["changes"])
        check("identical runs: the text says no change", "No change:" in txt, txt)
        check("identical runs: previous is the history folder name", d["previous"] == "2026-01-01T000000Z", d)

        # 3. one reclassification and one new digest
        s = summary("2026-02-01T00:00:00Z")
        s["sections"]["clinvar"]["data"]["hits"][1]["significance"] = "Benign"
        m = set_row(MANIFEST, "image", "deepvariant-ref", "digest-ddd")
        write_run(cur, s, m)
        d, txt, _ = delta(w, prev, cur, "planted")
        kinds = sorted((c["area"], c["kind"], c["item"]) for c in d["changes"])
        check("planted: exactly the two planted changes",
              kinds == [("ClinVar", "reclassified", "chr1:200 A>G (GENEA)"), ("inputs", "input", "DEEPVARIANT_IMAGE")],
              kinds)
        cv = next((c for c in d["changes"] if c["kind"] == "reclassified"), {})
        check("planted: the reclassification reads Likely pathogenic -> Benign",
              cv.get("before") == "Likely pathogenic" and cv.get("after") == "Benign", cv)
        check("planted: the new DeepVariant digest is named as a possible cause",
              "DEEPVARIANT_IMAGE" in cv.get("why", "") and "digest-ddd" in cv.get("why", ""), cv)
        img = next((c for c in d["changes"] if c["item"] == "DEEPVARIANT_IMAGE"), {})
        check("planted: the digest change shows both digests",
              "digest-bbb" in img.get("before", "") and "digest-ddd" in img.get("after", ""), img)
        check("planted: the text lists 2 changes", "2 change(s):" in txt, txt)

        # 3b. the same reclassification with only the PharmCAT digest changed
        m = set_row(MANIFEST, "image", "pharmcat-ref", "digest-eee")
        write_run(cur, s, m)
        d, _, _ = delta(w, prev, cur, "unrelated")
        cv = next((c for c in d["changes"] if c["kind"] == "reclassified"), {})
        check("a PharmCAT digest does not explain a ClinVar change",
              cv.get("why") == "the inputs this delta tracks for it are the same; the call itself differs", cv)

        # 4. one of each
        s = summary("2026-02-01T00:00:00Z")
        hits = s["sections"]["clinvar"]["data"]["hits"]
        hits[0]["stars"], hits[0]["review_status"] = 3, "reviewed by expert panel"
        hits[1]["significance"] = "Pathogenic"
        hits.append(hit(300, "Pathogenic", gene="GENEB"))
        s["sections"]["clinvar"]["data"]["hits"] = [h for h in hits if h["pos"] != 100] + [hit(400, "Pathogenic")]
        s["sections"]["clinvar"]["data"]["hits"].append(dict(hits[0], pos=500))
        s = copy.deepcopy(s)
        s["sections"]["cpic"]["data"]["genes"][0]["diplotype"] = "*1/*2"
        s["sections"]["cpic"]["data"]["genes"][0]["phenotype"] = "Intermediate Metabolizer"
        s["sections"]["pypgx"]["data"]["cyp2d6"] = "*1/*4"
        s["sections"]["prs"]["data"]["scores"][0]["percentile"] = "52.0"
        s["sections"]["prs"]["data"]["scores"][1]["percentile"] = "53.0"
        s["sections"]["cpsr"]["data"]["classification"] = {"Benign": 10, "VUS": 4}
        s["sections"]["roh"]["state"] = "missing"
        m = set_row(MANIFEST, "data", "clinvar_file_date", "20260201 (clinvar/clinvar.vcf.gz)")
        m = set_row(m, "version", "VEP_CACHE_RELEASE", "116")
        write_run(cur, s, m)
        d, txt, _ = delta(w, prev, cur, "many")
        got = sorted((c["area"], c["kind"], c["item"]) for c in d["changes"])
        want = sorted([
            ("inputs", "input", "VEP_CACHE_RELEASE"),
            ("inputs", "input", "clinvar_file_date"),
            ("Runs of homozygosity", "section", "roh"),
            ("ClinVar", "dropped", "chr1:100 A>G (GENEA)"),
            ("ClinVar", "reclassified", "chr1:200 A>G (GENEA)"),
            ("ClinVar", "added", "chr1:300 A>G (GENEB)"),
            ("ClinVar", "added", "chr1:400 A>G (GENEA)"),
            ("ClinVar", "added", "chr1:500 A>G (GENEA)"),
            ("CPSR", "classification count", "VUS"),
            ("PGx (CPIC)", "call", "CYP2C19"),
            ("PGx (pypgx)", "call", "CYP2D6"),
            ("PRS", "percentile", "PGS000001 (Trait one)"),
        ])
        check("each change once, the 3-point PRS move left out", got == want,
              {"extra": sorted(set(got) - set(want)), "missing": sorted(set(want) - set(got))})
        cv = next((c for c in d["changes"] if c["kind"] == "reclassified"), {})
        check("the reclassification names the new ClinVar file date",
              "clinvar_file_date 20260101" in cv.get("why", "") and "20260201" in cv.get("why", ""), cv)
        cp = next((c for c in d["changes"] if c["area"] == "CPSR"), {})
        check("the CPSR count names the new VEP cache release", "VEP_CACHE_RELEASE 115 -> 116" in cp.get("why", ""), cp)

        # 4b. re-starred, on its own
        s = summary("2026-02-01T00:00:00Z")
        s["sections"]["clinvar"]["data"]["hits"][0]["stars"] = 3
        s["sections"]["clinvar"]["data"]["hits"][0]["review_status"] = "reviewed by expert panel"
        write_run(cur, s, MANIFEST)
        d, _, _ = delta(w, prev, cur, "restarred")
        got = [(c["kind"], c["before"], c["after"]) for c in d["changes"]]
        check("a new review status is one re-starred change",
              got == [("re-starred", "2 stars (criteria provided, multiple submitters, no conflicts)",
                       "3 stars (reviewed by expert panel)")], got)

        # 5. HTML escaping
        s = summary("2026-02-01T00:00:00Z")
        s["sections"]["clinvar"]["data"]["hits"][1]["significance"] = "<b>x</b>"
        write_run(cur, s, MANIFEST)
        _, _, page = delta(w, prev, cur, "escape")
        check("the HTML escapes the values", "&lt;b&gt;x&lt;/b&gt;" in page and "<b>x</b>" not in page, page[-600:])
    finally:
        shutil.rmtree(w)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
