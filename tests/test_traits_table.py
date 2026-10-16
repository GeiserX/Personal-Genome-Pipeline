#!/usr/bin/env python3
"""assets/traits.tsv, the curated table of step 42: its schema, its reference
bases and its sources.

  python3 tests/test_traits_table.py                    schema (the unit-tests job)
  python3 tests/test_traits_table.py --fasta REF.fasta  also: each row's ref is the
                                                        FASTA's base, on every contig
                                                        the FASTA has (needs REF.fasta.fai)
  python3 tests/test_traits_table.py --online           also: every PMID resolves in
                                                        PubMed, and every rsID's GRCh38
                                                        position, reference base and
                                                        alt allele are dbSNP's (curl)
  --table FILE checks another table (the e2e case plants a broken copy).

The e2e case traits-panel-1 runs it with --fasta (the fixture's GRCh38
reference) and --online, and on a planted copy with one wrong reference base
and one PubMed id that does not exist, which must fail.

Checks without options:
  1. the table passes bin/traits_report.py's schema check (columns, chr1 to
     chr22, positions, one base each for ref and alt, readings, ancestry
     yes/no, PubMed ids, the two APOE rows, no SNPedia text);
  2. planted problems are each caught: a wrong column count, ref = alt, a
     non-base allele, a duplicate rsID, a PMID that is not a number, a
     missing APOE row;
  3. the FASTA check on a three-base synthetic FASTA: the right base passes,
     a planted wrong base fails, a contig the FASTA lacks is counted as not
     checked.

Run: python3 tests/test_traits_table.py
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import traits_report as tr  # noqa: E402

TABLE = os.path.join(REPO, "assets", "traits.tsv")
FAILS = 0
NC = {f"NC_0000{n:02d}": f"chr{n}" for n in range(1, 23)}


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1


def fasta_bases(fasta, rows):
    """{rsid: base} from FASTA at each row's position, by the .fai; rows on a
    contig the FASTA lacks are left out."""
    fai = {}
    with open(fasta + ".fai") as f:
        for line in f:
            c = line.split("\t")
            fai[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
    out = {}
    with open(fasta, "rb") as fh:
        for r in rows:
            if r["chrom"] not in fai:
                continue
            length, offset, bases, width = fai[r["chrom"]]
            i = int(r["pos"]) - 1
            if i >= length:
                out[r["rsid"]] = "past the end of the contig"
                continue
            fh.seek(offset + (i // bases) * width + i % bases)
            out[r["rsid"]] = fh.read(1).decode().upper()
    return out


def fasta_errors(fasta, rows):
    """(problems, rows checked) of the reference bases against FASTA."""
    got = fasta_bases(fasta, rows)
    errs = [f"{r['rsid']} {r['chrom']}:{r['pos']}: table ref {r['ref']}, FASTA {got[r['rsid']]}"
            for r in rows if r["rsid"] in got and got[r["rsid"]] != r["ref"]]
    return errs, len(got)


def curl_json(url, tries=4):
    last = ""
    for t in range(tries):
        res = subprocess.run(["curl", "-sS", "--fail", "--max-time", "60", url], capture_output=True, text=True)
        if res.returncode == 0:
            try:
                return json.loads(res.stdout)
            except ValueError as e:
                last = f"not JSON: {e}"
        else:
            last = res.stderr.strip()
        time.sleep(5 * (t + 1))
    raise RuntimeError(f"{url}: {last}")


def online_errors(rows):
    """Problems with the sources: PubMed ids that do not resolve, and rsIDs
    whose GRCh38 placement, reference base or alleles differ from dbSNP's."""
    errs = []
    pmids = sorted({p for r in rows for p in r["pmids"].split(",")}, key=int)
    d = curl_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
                  + ",".join(pmids))
    res = d.get("result", {})
    for p in pmids:
        rec = res.get(p) or {}
        if "error" in rec or not rec.get("title"):
            errs.append(f"PMID {p} does not resolve in PubMed ({rec.get('error', 'no record')})")
    print(f"  {len(pmids)} PubMed ids looked up")
    for r in rows:
        time.sleep(0.4)  # NCBI allows three requests a second without a key
        try:
            snp = curl_json(f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{r['rsid'][2:]}")
        except RuntimeError as e:
            errs.append(f"{r['rsid']}: dbSNP lookup failed: {e}")
            continue
        psd = snp.get("primary_snapshot_data")
        if not psd:
            errs.append(f"{r['rsid']}: not a current dbSNP record (merged or withdrawn)")
            continue
        place = None
        for pl in psd.get("placements_with_allele", []):
            traits = pl.get("placement_annot", {}).get("seq_id_traits_by_assembly", [])
            if pl.get("is_ptlp") and traits and traits[0].get("assembly_name", "").startswith("GRCh38"):
                place = pl
                break
        if place is None:
            errs.append(f"{r['rsid']}: no GRCh38 placement in dbSNP")
            continue
        spdi = [a["allele"]["spdi"] for a in place["alleles"]]
        chrom = NC.get(place["seq_id"].split(".")[0], place["seq_id"])
        pos, ref = spdi[0]["position"] + 1, spdi[0]["deleted_sequence"]
        alts = {s["inserted_sequence"] for s in spdi if s["inserted_sequence"] != s["deleted_sequence"]}
        if (chrom, str(pos), ref) != (r["chrom"], r["pos"], r["ref"]):
            errs.append(f"{r['rsid']}: table {r['chrom']}:{r['pos']} {r['ref']}, dbSNP GRCh38 {chrom}:{pos} {ref}")
        elif r["alt"] not in alts:
            errs.append(f"{r['rsid']}: alt {r['alt']} is not one of dbSNP's alleles {sorted(alts)}")
    print(f"  {len(rows)} rsIDs looked up in dbSNP")
    return errs


def schema_controls():
    good = tr.read_table(TABLE)
    errs = tr.check_rows(good)
    check(f"assets/traits.tsv passes the schema check ({len(good)} rows)", not errs, errs)
    check("the table has at least 30 rows", len(good) >= 30, len(good))

    def broken(desc, mutate, want):
        rows = [dict(r) for r in good]
        mutate(rows)
        e = tr.check_rows(rows)
        check(f"control: {desc} is caught", any(want in x for x in e), e)

    broken("ref = alt", lambda rs: rs[0].update(alt=rs[0]["ref"]), "same base")
    broken("a ref that is not a base", lambda rs: rs[0].update(ref="N"), "is not one of A, C, G, T")
    broken("a duplicate rsID", lambda rs: rs[1].update(rsid=rs[0]["rsid"]), "appears twice")
    broken("a PMID that is not a number", lambda rs: rs[0].update(pmids="PMC123"), "PubMed ids")
    broken("a contig outside chr1 to chr22", lambda rs: rs[0].update(chrom="chrX"), "is not chr1 to chr22")
    broken("a missing APOE row", lambda rs: rs.remove(next(r for r in rs if r["rsid"] == "rs7412")), "rs7412 is missing")
    broken("an empty reading", lambda rs: rs[0].update(het=" "), "het is empty")
    broken("SNPedia text", lambda rs: rs[0].update(note="From SNPedia."), "SNPedia")
    work = tempfile.mkdtemp()
    try:
        bad = os.path.join(work, "t.tsv")
        with open(TABLE) as f:
            text = f.read()
        with open(bad, "w") as f:
            f.write(text.replace("\trs671\t", "\trs671\textra\t", 1))
        try:
            tr.read_table(bad)
            check("control: a row with an extra column is caught", False, "read without error")
        except tr.TableError as e:
            check("control: a row with an extra column is caught", "columns" in str(e), e)

        # FASTA check on a synthetic FASTA: chr2 holds the lactase row's base.
        lct = next(r for r in good if r["rsid"] == "rs4988235")
        pos = int(lct["pos"])
        fa = os.path.join(work, "ref.fasta")
        seq = ["A"] * pos
        seq[pos - 1] = lct["ref"]
        width = 60
        body = "".join("".join(seq[i:i + width]) + "\n" for i in range(0, pos, width))
        with open(fa, "w") as f:
            f.write(">chr2\n" + body)
        with open(fa + ".fai", "w") as f:
            f.write(f"chr2\t{pos}\t6\t{width}\t{width + 1}\n")
        e, n = fasta_errors(fa, [lct])
        check("FASTA check: the right base passes", not e and n == 1, (e, n))
        wrong = dict(lct, ref=next(b for b in "ACGT" if b != lct["ref"]))
        e, n = fasta_errors(fa, [wrong])
        check("FASTA check control: a planted wrong reference base fails, naming the row",
              len(e) == 1 and "rs4988235" in e[0], e)
        other = next(r for r in good if r["chrom"] != "chr2")
        e, n = fasta_errors(fa, [other])
        check("FASTA check: a row on a contig the FASTA lacks is not checked", not e and n == 0, (e, n))
    finally:
        shutil.rmtree(work)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", default=TABLE)
    ap.add_argument("--fasta", help="check each row's ref against this FASTA (with its .fai)")
    ap.add_argument("--min-fasta-rows", type=int, default=1,
                    help="fail unless at least this many rows were checked against the FASTA")
    ap.add_argument("--online", action="store_true", help="look the PMIDs up in PubMed and the rsIDs in dbSNP")
    a = ap.parse_args(argv)
    if a.table == TABLE and not a.fasta and not a.online:
        schema_controls()
    else:
        try:
            rows = tr.read_table(a.table)
        except tr.TableError as e:
            check("the table can be read", False, e)
            rows = []
        errs = tr.check_rows(rows)
        check(f"{a.table} passes the schema check", not errs, errs)
        if a.fasta and rows:
            e, n = fasta_errors(a.fasta, rows)
            print(f"  ref checked against {a.fasta}: {n} rows (the others are on contigs it lacks)")
            check(f"every reference base matches the FASTA ({n} rows checked)", not e, e)
            check(f"at least {a.min_fasta_rows} rows were checked against the FASTA", n >= a.min_fasta_rows, n)
        if a.online and rows:
            try:
                e = online_errors(rows)
            except RuntimeError as x:
                e = [f"lookup failed: {x}"]
            check("every PMID resolves and every rsID matches dbSNP's GRCh38 record", not e, e)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
