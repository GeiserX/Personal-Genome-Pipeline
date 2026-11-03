#!/usr/bin/env python3
"""Step 20's confirmed mtDNA variant list, large deletions (eKLIPse) and mtDNA
copy number, on synthetic inputs.

Checks:
  1. assets/mito_confirmed.tsv: its header, and every row with a ClinVar id
     or a PubMed id (a planted row with neither fails the schema check); the
     first rows are there (m.3243A>G, m.8344A>G, the three LHON variants and
     the common 4977 bp deletion).
  2. a chrM VCF with m.3243A>G at 12% is reported with its level, in the
     summary and in both reports; the same VCF without it reports none.
  3. the copy number: idxstats with chrM at 100 times the autosomal depth
     gives 200 copies per cell (+-5%); a sliced BAM and a BAM with no chrM
     reads are not assessed.
  4. eKLIPse's deletions CSV (decimal commas, quoted fields) is read; the
     common deletion is named within 50 bp of its breakpoints, a deletion
     elsewhere is not, one under 1% is counted but not listed, and a table
     with no deletion lists none.
  5. EKLIPSE_IMAGE is pinned by digest, and its versions.env comment names
     the upstream commit 3606cb2edac983d2623ddc667b49206c3d01373c (the e2e
     case checks the Bioconda build's recipe against it).

Run: python3 tests/test_mito_extended.py
Also: python3 tests/test_mito_extended.py synth-sam ... writes the synthetic
chrM reads tests/e2e/mito-extended-1-deletions.sh aligns nothing for: SAM
records with the soft clips an aligner writes at a deletion's junction.
"""
import argparse
import csv
import gzip
import os
import random
import re
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))

EKLIPSE_COMMIT = "3606cb2edac983d2623ddc667b49206c3d01373c"
FAILS = 0


def check(name, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" ({detail})" if detail and not ok else ""))
    if not ok:
        FAILS += 1


# --- the synthetic reads (synth-sam) ------------------------------------------------

def read_fasta(path):
    with open(path) as f:
        return "".join(ln.strip() for ln in f if not ln.startswith(">")).upper()


def synth_reads(seq, depth, deletion, level, read_len=150, seed=1):
    """[(pos 1-based, cigar, flag, read sequence)] of single-end reads at the
    given depth from molecules of which `level` carry the deletion (first and
    last deleted base, 1-based), or none when deletion is None. A read across
    the junction is written as an aligner writes it: aligned on its longer
    side, the rest soft-clipped (a forward read clipped on the right, a
    reverse read on the left, as eKLIPse reads them); a junction read with a
    side under 30 bases is left out."""
    rnd = random.Random(seed)
    n = int(len(seq) * depth / read_len)
    out = []
    if deletion:
        a, b = deletion
        mol = seq[:a - 1] + seq[b:]
        junction = a - 1  # first base of mol after the junction (0-based)
    for i in range(n):
        if deletion and rnd.random() < level:
            s = rnd.randrange(0, len(mol) - read_len)
            left = junction - s
            r = mol[s:s + read_len]
            if left <= 0:
                out.append((s + 1 + (b - a + 1), f"{read_len}M", rnd.choice((0, 16)), r))
            elif left >= read_len:
                out.append((s + 1, f"{read_len}M", rnd.choice((0, 16)), r))
            elif min(left, read_len - left) < 30:
                continue
            elif left >= read_len - left:
                out.append((s + 1, f"{left}M{read_len - left}S", 0, r))
            else:
                out.append((b + 1, f"{left}S{read_len - left}M", 16, r))
        else:
            s = rnd.randrange(0, len(seq) - read_len)
            out.append((s + 1, f"{read_len}M", rnd.choice((0, 16)), seq[s:s + read_len]))
    return sorted(out)


def synth_sam_main(argv):
    ap = argparse.ArgumentParser(prog="test_mito_extended.py synth-sam")
    ap.add_argument("--fasta", required=True, help="chrM alone (samtools faidx REF chrM)")
    ap.add_argument("--header", required=True, help="the @HD/@SQ lines the reads go under")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--depth", type=float, default=400)
    ap.add_argument("--deletion", help="FIRST-LAST deleted base, e.g. 8470-13446")
    ap.add_argument("--level", type=float, default=0.1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    seq = read_fasta(a.fasta)
    dele = tuple(int(x) for x in a.deletion.split("-")) if a.deletion else None
    reads = synth_reads(seq, a.depth, dele, a.level)
    with open(a.header) as f, open(a.out, "w") as o:
        o.write("".join(ln for ln in f if ln.startswith(("@HD", "@SQ"))))
        o.write(f"@RG\tID:{a.sample}\tSM:{a.sample}\tPL:ILLUMINA\tLB:{a.sample}\n")
        for i, (pos, cigar, flag, r) in enumerate(reads):
            o.write(f"r{i}\t{flag}\tchrM\t{pos}\t60\t{cigar}\t*\t0\t0\t{r}\t{'I' * len(r)}\tRG:Z:{a.sample}\n")
    clipped = sum(1 for _, c, _, _ in reads if "S" in c)
    print(f"{len(reads)} reads, {clipped} soft-clipped at the junction")
    return 0


# --- the checks -------------------------------------------------------------------

CONFIRMED_COLS = ["kind", "pos", "end", "ref", "alt", "hgvs", "gene", "disease", "clinvar_id", "pmids",
                  "heteroplasmy_note"]


def confirmed_errors(path):
    """What is wrong with a mito_confirmed.tsv: a wrong header, a row with no
    ClinVar id and no PubMed id, a malformed id, an unknown kind."""
    errs = []
    with open(path) as f:
        rows = list(csv.DictReader((ln for ln in f if not ln.startswith("#")), delimiter="\t",
                                   quoting=csv.QUOTE_NONE))
    if not rows or list(rows[0].keys()) != CONFIRMED_COLS:
        return [f"header is not {' '.join(CONFIRMED_COLS)}"]
    for r in rows:
        name = r.get("hgvs") or r.get("pos")
        if None in r or None in r.values():
            errs.append(f"{name}: wrong number of columns")
            continue
        cv = r["clinvar_id"] not in ("", ".")
        pm = r["pmids"] not in ("", ".")
        if not cv and not pm:
            errs.append(f"{name}: neither a ClinVar id nor a PubMed id")
        if cv and not r["clinvar_id"].isdigit():
            errs.append(f"{name}: ClinVar id {r['clinvar_id']!r}")
        if pm and not all(x.isdigit() for x in r["pmids"].split(";")):
            errs.append(f"{name}: PubMed ids {r['pmids']!r}")
        if r["kind"] not in ("snv", "deletion"):
            errs.append(f"{name}: kind {r['kind']!r}")
        elif not (r["pos"].isdigit() and r["end"].isdigit() and 1 <= int(r["pos"]) <= int(r["end"]) <= 16569):
            errs.append(f"{name}: positions {r['pos']}-{r['end']}")
    return errs


def test_confirmed_schema(tmp):
    path = os.path.join(REPO, "assets", "mito_confirmed.tsv")
    errs = confirmed_errors(path)
    check("assets/mito_confirmed.tsv: every row has a ClinVar or PubMed id and a valid shape", not errs, "; ".join(errs))
    with open(path) as f:
        text = f.read()
    for want in ("m.3243A>G", "m.8344A>G", "m.11778G>A", "m.3460G>A", "m.14484T>C"):
        check(f"the list has {want}", f"\t{want}\t" in text)
    check("the list has the common 4977 bp deletion (m.8470 to m.13446)",
          re.search(r"^deletion\t8470\t13446\t", text, re.M) is not None)
    bad = os.path.join(tmp, "bad.tsv")
    with open(bad, "w") as f:
        f.write(text + "snv\t3271\t3271\tT\tC\tm.3271T>C\tMT-TL1\tMELAS\t.\t.\tno id at all\n")
    errs = confirmed_errors(bad)
    check("red control: a planted row with neither id fails the schema check",
          any("neither a ClinVar id nor a PubMed id" in e for e in errs), "; ".join(errs) or "no error")


def write_vcf(path, records):
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for pos, ref, alt, af in records:
            f.write(f"chrM\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:AF\t0/1:{af}\n")


def sample_dir(tmp, name, records, deletions=None, idxstats=None):
    import collect_summary
    d = os.path.join(tmp, name)
    os.makedirs(os.path.join(d, "mito"))
    write_vcf(os.path.join(d, "mito", f"{name}_chrM_filtered.vcf.gz"), records)
    if deletions is not None:
        csv_path = os.path.join(tmp, f"{name}_eklipse.csv")
        with open(csv_path, "w") as f:
            f.write("\"Title\";\"5' breakpoint\";\"3' breakpoint\";\"Freq\";\"Freq For\";\"Freq Rev\";\"5' Blast\";"
                    "\"3' Blast\";\"5' Depth\";\"3' Depth\";\"Repetition\"\n")
            for b5, b3, lvl in deletions:
                v = str(lvl).replace(".", ",")
                f.write(f'"{name}";"{b5}";"{b3}";"{v}";"{v}";"{v}";"12";"14";"390.0";"402.0";"none"\n')
        rc = collect_summary.main(["mito-deletions", "--eklipse", csv_path,
                                   "--out", os.path.join(d, "mito", f"{name}_chrM_deletions.tsv")])
        check(f"{name}: mito-deletions reads eKLIPse's CSV", rc == 0)
    if idxstats is not None:
        ip = os.path.join(tmp, f"{name}.idxstats")
        with open(ip, "w") as f:
            f.write(idxstats)
        rc = collect_summary.main(["mito-copies", "--idxstats", ip,
                                   "--out", os.path.join(d, "mito", f"{name}_mtdna_copies.tsv")])
        check(f"{name}: mito-copies reads the idxstats", rc == 0)
    return d


def mito_section(d, name):
    import collect_summary
    _, data = collect_summary.sec_mito(d, name)
    return data


def reports(d, name):
    import collect_summary
    import render_report
    s = collect_summary.collect(name, d)
    return render_report.text_report(s), render_report.html_report(s)


def test_confirmed_snv(tmp):
    d = sample_dir(tmp, "S3243", [(73, "A", "G", "1.00"), (3243, "A", "G", "0.12"), (16519, "T", "C", "0.998")])
    data = mito_section(d, "S3243")
    conf = data.get("confirmed") or []
    check("m.3243A>G at 12% is in the summary's confirmed list", [v.get("hgvs") for v in conf] == ["m.3243A>G"],
          str(conf))
    check("its level is 0.12", bool(conf) and conf[0].get("af") == "0.12", str(conf))
    text, html = reports(d, "S3243")
    check("the text report names m.3243A>G with AF 0.12 and MELAS",
          re.search(r"m\.3243A>G \(MT-TL1\)  AF 0\.12 .*MELAS", text) is not None)
    check("the HTML report counts one confirmed variant and names it",
          re.search(r'Confirmed disease variants</span><span class="value"><span class="badge badge-yellow">1<', html)
          is not None and "m.3243A&gt;G (MT-TL1)  AF 0.12" in html)
    d = sample_dir(tmp, "Snone", [(73, "A", "G", "1.00"), (16519, "T", "C", "0.998")])
    data = mito_section(d, "Snone")
    check("control: without m.3243A>G the confirmed list is empty", data.get("confirmed") == [], str(data.get("confirmed")))
    text, _ = reports(d, "Snone")
    check("control: the text report counts no confirmed variant",
          "Confirmed disease variants (assets/mito_confirmed.tsv): 0" in text)
    import collect_summary
    saved, collect_summary.ASSETS_DIR = collect_summary.ASSETS_DIR, os.path.join(tmp, "no-assets")
    collect_summary._ASSETS.clear()
    try:
        data = mito_section(d, "Snone")
        text, html = reports(d, "Snone")
    finally:
        collect_summary.ASSETS_DIR = saved
        collect_summary._ASSETS.clear()
    check("without the list the summary has no confirmed count, only why",
          "confirmed" not in data and "mito_confirmed.tsv is not in" in (data.get("confirmed_note") or ""), str(data))
    check("without the list neither report shows a count of confirmed variants",
          "Confirmed disease variants (assets" not in text and "Confirmed disease variants</span>" not in html
          and "Confirmed disease variants: not checked" in text)


def idxstats(chrm_reads, auto_reads_per_chrom, chroms=22, length=50_000_000):
    lines = [f"chr{i}\t{length}\t{auto_reads_per_chrom}\t0" for i in range(1, chroms + 1)]
    lines += [f"chrX\t156040895\t{auto_reads_per_chrom}\t0", f"chrM\t16569\t{chrm_reads}\t0", "*\t0\t0\t500"]
    return "\n".join(lines) + "\n"


def test_copies(tmp):
    import collect_summary
    # 10 reads per kb of autosome would be under the floor; 200 per kb is about 30x.
    per_chrom = 200 * 50_000
    chrm = round(100 * 16569 * per_chrom / 50_000_000)
    b = collect_summary.mito_copies(idxstats(chrm, per_chrom))
    got = float(b["mtdna_copies"] or 0)
    check(f"chrM at 100x the autosomal depth reads 200 copies per cell (+-5%): {got}",
          b["status"] == "called" and abs(got - 200) <= 10, str(b))
    b = collect_summary.mito_copies(idxstats(chrm // 2, per_chrom))
    check("half the chrM reads reads half the copies", abs(float(b["mtdna_copies"] or 0) - 100) <= 5, str(b))
    b = collect_summary.mito_copies(idxstats(chrm, 500))
    check("a sliced BAM (few autosomal reads) is not assessed",
          b["status"] == "not assessed" and b["reason"] == "low_autosomal_depth", str(b))
    b = collect_summary.mito_copies(idxstats(0, per_chrom))
    check("a BAM with no chrM reads is not assessed", b["reason"] == "no_chrM_reads", str(b))
    b = collect_summary.mito_copies(idxstats(chrm, per_chrom).replace("chrM", "MT").replace("chr", ""))
    check("contigs named 1..22 and MT work", b["status"] == "called", str(b))
    d = sample_dir(tmp, "Scn", [(73, "A", "G", "1.00")], idxstats=idxstats(chrm, per_chrom))
    text, html = reports(d, "Scn")
    check("the text report gives the copy number",
          re.search(r"mtDNA copy number: 20[0-9]\.[0-9] per cell", text) is not None)
    check("the HTML report gives the copy number",
          re.search(r'mtDNA copies per cell</span><span class="value">20[0-9]\.[0-9]<', html) is not None)


def test_deletions(tmp):
    import collect_summary
    d = sample_dir(tmp, "Sdel", [(73, "A", "G", "1.00")],
                   deletions=[(8469, 13447, 10.2), (8482, 13460, 2.5), (9000, 14000, 3.1), (5000, 7000, 0.4)])
    rows = collect_summary.eklipse_deletions(os.path.join(tmp, "Sdel_eklipse.csv"))
    check("eKLIPse's decimal commas are read", [r["level_pct"] for r in rows] == [10.2, 2.5, 3.1, 0.4],
          str([r["level_pct"] for r in rows]))
    check("the deletion's length is the bases between the breakpoints", rows[0]["length"] == 4977, str(rows[0]))
    dl = mito_section(d, "Sdel").get("deletions") or {}
    listed = dl.get("listed") or []
    named = [(r["bkp5"], r["bkp3"]) for r in listed if r.get("known")]
    check("the common deletion is named at both of its writings across the 13 bp repeat",
          named == [(8469, 13447), (8482, 13460)], str(listed))
    check("a deletion elsewhere (m.9001-13999) is listed but not named",
          any(r["bkp5"] == 9000 and not r.get("known") for r in listed), str(listed))
    check("a deletion under 1% is counted, not listed", dl.get("total") == 4 and len(listed) == 3, str(dl))
    text, html = reports(d, "Sdel")
    check("the HTML report counts four deletions, two of them named",
          '>Large deletions (eKLIPse)</span><span class="value">4<' in html
          and re.search(r'Named in the confirmed list</span><span class="value"><span class="badge badge-yellow">2<', html)
          is not None)
    check("the text report names the common deletion with its level and ClinVar id",
          re.search(r"m\.8470_13446del \(4977 bp\)  level 10\.2%  The common 4977 bp deletion.*ClinVar 430680", text)
          is not None)
    check("the text report says blood often lacks the deletion", "often found in muscle and not in blood" in text)
    d = sample_dir(tmp, "Snodel", [(73, "A", "G", "1.00")], deletions=[])
    dl = mito_section(d, "Snodel").get("deletions") or {}
    check("red control: a table with no deletion lists none", dl.get("total") == 0 and dl.get("listed") == [], str(dl))
    bad = os.path.join(tmp, "bad.csv")
    with open(bad, "w") as f:
        f.write("pos\tref\n")
    try:
        collect_summary.eklipse_deletions(bad)
        ok = False
    except ValueError:
        ok = True
    check("a file that is not eKLIPse's table is refused", ok)
    wrap = os.path.join(tmp, "wrap.csv")
    with open(wrap, "w") as f:
        f.write("\"Title\";\"5' breakpoint\";\"3' breakpoint\";\"Freq\";\"Freq For\";\"Freq Rev\";\"5' Blast\";"
                "\"3' Blast\";\"5' Depth\";\"3' Depth\";\"Repetition\"\n")
        f.write('"S";"16000";"500";"2,0";"2,0";"2,0";"3";"3";"100";"100";"none"\n')
    rows = collect_summary.eklipse_deletions(wrap)
    check("a deletion across the origin (m.16001 to m.499) gets its length, 1068 bp",
          rows[0]["length"] == 1068, str(rows[0]))


def test_image_pin():
    with open(os.path.join(REPO, "versions.env")) as f:
        line = next((ln for ln in f if ln.startswith("EKLIPSE_IMAGE=")), "")
    check("EKLIPSE_IMAGE is pinned by digest", re.match(r'EKLIPSE_IMAGE="[^"@]+:[^"@]+@sha256:[0-9a-f]{64}"', line)
          is not None, line.strip())
    check(f"its comment names the upstream commit {EKLIPSE_COMMIT[:7]}", EKLIPSE_COMMIT in line, line.strip())


def main(argv):
    if argv[:1] == ["synth-sam"]:
        return synth_sam_main(argv[1:])
    tmp = tempfile.mkdtemp()
    try:
        test_confirmed_schema(tmp)
        test_confirmed_snv(tmp)
        test_copies(tmp)
        test_deletions(tmp)
        test_image_pin()
    finally:
        shutil.rmtree(tmp)
    print(f"{FAILS} failed" if FAILS else "all passed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
