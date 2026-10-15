#!/usr/bin/env python3
"""bin/clinvar_join_adjacent.awk, the combined-allele joiner of step 06 and the
CLINVAR_SCREEN module, on hand-made records at real GRCh38 chr22 positions.

Checks, with every awk on the machine (the module runs the bcftools image's
awk, mawk):
  1. two records within 10 bp, both het with the same unphased GT, are joined
     into one REF/ALT, the bases between them taken from the reference;
     0/1 and 1/0 are the same het call;
  2. phased records are joined only in the same phase set (FORMAT/PS) with
     the same GT: other phase sets and the other haplotype are not joined;
  3. a het next to a hom, a half call (./1), overlapping records and records
     11 bp apart are not joined; 10 bp apart are;
  4. a single SNV gives no combined record;
  5. finish drops a combined record whose part disagrees with the reference;
  6. the combined records come out in POS order;
  7. the tag pass marks each part of a matched combined allele "part of",
     flags two hits in one gene that can be one allele, and leaves a single
     SNV hit, hits in two genes and hits on two haplotypes unchanged;
  8. with docker (always in CI; elsewhere with PGP_TEST_DOCKER=1), the real
     scripts/06-clinvar-screen.sh on the fixture records in
     tests/fixtures/vcf/combined_allele_*.vcf: the 2-star ClinVar record of the combined allele is a hit, its two 1-star
     parts are hits marked as its parts, and the single SNV hit is unchanged.

The window is 10 bp; set it to 0 in the awk file and checks 1, 2, 3, 6, 7 and
8 fail.

Run: python3 tests/test_clinvar_join_adjacent.py
"""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AWK_FILE = os.path.join(REPO, "bin", "clinvar_join_adjacent.awk")
FIXTURES = os.path.join(REPO, "tests", "fixtures", "vcf")

# GRCh38 chr22 (the same in the NCBI no-ALT analysis set), two short stretches.
REF = {30000000: "AGACTCTGTCGCAAAAAAAAAATTAATGTAACCATTTAGCCATTACCCAGTTCCTTAAAAT",
       30000200: "GCCTGTAATCCCAGCACTTTGGGAGGCTGAG"}

HEAD = ("##fileformat=VCFv4.2\n##contig=<ID=chr22,length=50818468>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)}")
    if not ok:
        FAILS += 1


def ref_at(pos, n):
    for start, seq in REF.items():
        if start <= pos and pos + n <= start + len(seq):
            return seq[pos - start:pos - start + n]
    raise KeyError(f"chr22:{pos} (+{n}) is outside the test's reference stretches")


def rec(pos, ref, alt, gt, ps=None):
    fmt, smp = ("GT", gt) if ps is None else ("GT:PS", f"{gt}:{ps}")
    return f"chr22\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\t{fmt}\t{smp}\n"


def awks():
    found = []
    for name in ("awk", "mawk", "gawk", "original-awk", "nawk"):
        p = shutil.which(name)
        if p and os.path.realpath(p) not in [os.path.realpath(x[1]) for x in found]:
            found.append((name, p))
    if shutil.which("busybox"):
        found.append(("busybox awk", "busybox"))
    return found


def awk_run(awk, args, stdin=""):
    cmd = ([awk] if awk != "busybox" else ["busybox", "awk"]) + args
    p = subprocess.run(cmd, input=stdin, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def fill(vcf):
    """What `bcftools +fill-from-fasta -- -c INFO/PGP_JOIN_SEQ` does: the
    reference bases from POS for the length of REF."""
    out = []
    for line in vcf.splitlines(True):
        if line.startswith("#"):
            out.append(line)
            continue
        f = line.rstrip("\n").split("\t")
        f[7] = f[7].replace("PGP_JOIN_SEQ=.", "PGP_JOIN_SEQ=" + ref_at(int(f[1]), len(f[3])), 1)
        out.append("\t".join(f) + "\n")
    return "".join(out)


def join(awk, records, fill_fn=fill):
    """prepare | fill | finish; the records out, as lists of fields."""
    rc, prep, err = awk_run(awk, ["-v", "mode=prepare", "-f", AWK_FILE], HEAD + "".join(records))
    if rc != 0:
        return None, f"prepare exit {rc}: {err}"
    rc, done, err = awk_run(awk, ["-v", "mode=finish", "-f", AWK_FILE], fill_fn(prep))
    if rc != 0:
        return None, f"finish exit {rc}: {err}"
    return [l.split("\t") for l in done.splitlines() if l and not l.startswith("#")], done


def core(r):
    return (r[1], r[3], r[4], r[8], r[9])


def awk_checks(name, awk, work):
    # 1. unphased het pair, 3 bp apart: AGC>A (a 2 bp deletion) and T>TGGG
    out, raw = join(awk, [rec(30000037, "AGC", "A", "0/1"), rec(30000042, "T", "TGGG", "0/1")])
    check(f"{name}: two het records 3 bp apart become one REF/ALT with the reference bases between them",
          out is not None and [core(r) for r in out] == [("30000037", "AGCCAT", "ACATGGG", "GT", "0/1")], raw)
    check(f"{name}: the combined record names its two parts",
          out is not None and len(out) == 1
          and out[0][7] == "PGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG", raw)
    check(f"{name}: the header keeps PGP_JOINED_FROM and drops the temporary PGP_JOIN_SEQ",
          "ID=PGP_JOINED_FROM" in raw and "ID=PGP_JOIN_SEQ" not in raw, raw)
    out, raw = join(awk, [rec(30000203, "T", "A", "1/0"), rec(30000208, "T", "C", "0/1")])
    check(f"{name}: 1/0 and 0/1 (unphased) are the same het call and are joined",
          out is not None and [core(r) for r in out] == [("30000203", "TGTAAT", "AGTAAC", "GT", "0/1")], raw)

    # 2. phase sets
    out, raw = join(awk, [rec(30000203, "T", "A", "0|1", 100), rec(30000208, "T", "C", "0|1", 100)])
    check(f"{name}: phased, same phase set and GT: joined, with GT and PS",
          out is not None and [core(r) for r in out] == [("30000203", "TGTAAT", "AGTAAC", "GT:PS", "0|1:100")], raw)
    out, raw = join(awk, [rec(30000203, "T", "A", "0|1", 100), rec(30000208, "T", "C", "0|1", 200)])
    check(f"{name}: phased on different phase sets: not joined", out == [], raw)
    out, raw = join(awk, [rec(30000203, "T", "A", "0|1", 100), rec(30000208, "T", "C", "1|0", 100)])
    check(f"{name}: phased on the two haplotypes of one phase set: not joined", out == [], raw)
    out, raw = join(awk, [rec(30000203, "T", "A", "0|1", 100), rec(30000208, "T", "C", "0/1")])
    check(f"{name}: one phased and one unphased record: not joined", out == [], raw)

    # 3. what is not one allele, and the window's edge
    out, raw = join(awk, [rec(30000203, "T", "A", "0/1"), rec(30000208, "T", "C", "1/1")])
    check(f"{name}: a het next to a hom: not joined", out == [], raw)
    out, raw = join(awk, [rec(30000203, "T", "A", "./1"), rec(30000208, "T", "C", "./1")])
    check(f"{name}: half calls (./1): not joined", out == [], raw)
    out, raw = join(awk, [rec(30000203, "TGT", "T", "0/1"), rec(30000204, "G", "A", "0/1")])
    check(f"{name}: overlapping records: not joined", out == [], raw)
    out, raw = join(awk, [rec(30000203, "T", "A", "0/1"), rec(30000213, "G", "T", "0/1")])
    check(f"{name}: 10 bp apart: joined",
          out is not None and [core(r)[:3] for r in out] == [("30000203", "TGTAATCCCAG", "AGTAATCCCAT")], raw)
    out, raw = join(awk, [rec(30000203, "T", "A", "0/1"), rec(30000214, "C", "A", "0/1")])
    check(f"{name}: 11 bp apart: not joined", out == [], raw)

    # 4. a single SNV
    out, raw = join(awk, [rec(30000210, "C", "T", "0/1")])
    check(f"{name}: a single SNV gives no combined record", out == [], raw)

    # 5. a part whose REF is not the reference's base (30000203 is T)
    out, raw = join(awk, [rec(30000203, "C", "A", "0/1"), rec(30000208, "T", "C", "0/1")])
    check(f"{name}: a part that disagrees with the reference: dropped by finish", out == [], raw)

    # 6. POS order: a 12 bp deletion that ends after the two records inside it
    recs = [rec(30000201, "CCTGTAATCCCA", "C", "0/1"), rec(30000205, "T", "G", "0/1"),
            rec(30000208, "T", "C", "0/1"), rec(30000214, "C", "A", "0/1")]
    out, raw = join(awk, recs)
    pos = [int(r[1]) for r in out] if out else []
    check(f"{name}: four combined records, in POS order", pos == [30000201, 30000205, 30000205, 30000208], raw)

    # mode is required
    rc, _, err = awk_run(awk, ["-f", AWK_FILE], HEAD)
    check(f"{name}: no mode: exit 2 with a message", rc == 2 and "mode" in err, f"exit {rc}: {err}")

    # 7. the tag pass on a hits VCF
    hits = os.path.join(work, f"hits-{name.replace(' ', '_')}.vcf")
    cv = ";GENEINFO={g};CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter"
    lines = [
        "chr22\t30000037\t9900001\tAGC\tA\t50\tPASS\t" + cv.format(g="SYNTH1:1")[1:] + "\tGT\t0/1",
        "chr22\t30000038\t9900003\tGCCAT\tCATGGG\t.\tPASS\tPGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG"
        + cv.format(g="SYNTH1:1").replace("single_submitter", "multiple_submitters,_no_conflicts") + "\tGT\t0/1",
        "chr22\t30000042\t9900002\tT\tTGGG\t50\tPASS\t" + cv.format(g="SYNTH1:1")[1:] + "\tGT\t0/1",
        "chr22\t30000210\t9900004\tC\tT\t50\tPASS\t" + cv.format(g="SYNTH2:2")[1:] + "\tGT\t0/1",
        "chr22\t30000300\t9900005\tA\tG\t50\tPASS\t" + cv.format(g="SYNTH3:3")[1:] + "\tGT\t0/1",
        "chr22\t30000305\t9900006\tC\tT\t50\tPASS\t" + cv.format(g="SYNTH3:3")[1:] + "\tGT\t0/1",
        "chr22\t30000400\t9900007\tA\tG\t50\tPASS\t" + cv.format(g="SYNTH4:4")[1:] + "\tGT\t0/1",
        "chr22\t30000404\t9900008\tC\tT\t50\tPASS\t" + cv.format(g="SYNTH5:5")[1:] + "\tGT\t0/1",
        "chr22\t30000500\t9900009\tA\tG\t50\tPASS\t" + cv.format(g="SYNTH6:6")[1:] + "\tGT:PS\t0|1:7",
        "chr22\t30000504\t9900010\tC\tT\t50\tPASS\t" + cv.format(g="SYNTH6:6")[1:] + "\tGT:PS\t1|0:7",
    ]
    with open(hits, "w") as f:
        f.write(HEAD + "\n".join(lines) + "\n")
    rc, out, err = awk_run(awk, ["-v", "mode=tag", "-f", AWK_FILE, hits, hits])
    got = [l for l in out.splitlines() if l and not l.startswith("#")]
    info = {tuple(l.split("\t")[1:2]): l.split("\t")[7] for l in got}
    check(f"{name}: tag: exits 0, header and every record kept",
          rc == 0 and len(got) == len(lines) and out.startswith("##fileformat"), f"exit {rc}: {err}")
    part = "PGP_PART_OF=chr22:30000038:GCCAT:CATGGG"
    check(f"{name}: tag: both parts are marked part of the combined allele",
          info.get(("30000037",), "").endswith(part) and info.get(("30000042",), "").endswith(part), info)
    check(f"{name}: tag: the combined hit itself gets no part mark",
          "PGP_PART_OF" not in info.get(("30000038",), "") and "PGP_ONE_ALLELE" not in info.get(("30000038",), ""),
          info)
    check(f"{name}: tag: the single SNV hit is unchanged", lines[3] in got, info)
    check(f"{name}: tag: two hits in one gene, 5 bp apart, no combined record: flagged both ways",
          info.get(("30000300",), "").endswith("PGP_ONE_ALLELE=chr22:30000305:C:T")
          and info.get(("30000305",), "").endswith("PGP_ONE_ALLELE=chr22:30000300:A:G"), info)
    check(f"{name}: tag: hits in two genes, or on the two haplotypes of a phase set, are unchanged",
          all(x in got for x in lines[6:]), info)


def docker_checks(work):
    """Check 8: the real step 06 with the real bcftools image."""
    env = dict(os.environ)
    versions = {}
    with open(os.path.join(REPO, "versions.env")) as f:
        for line in f:
            if line.startswith("BCFTOOLS_IMAGE="):
                versions["BCFTOOLS_IMAGE"] = line.split("=", 1)[1].split("#")[0].strip().strip('"')
    image = versions["BCFTOOLS_IMAGE"]
    gd = os.path.join(work, "genome")
    for d in ("reference", "clinvar", "S1/vcf"):
        os.makedirs(os.path.join(gd, d), exist_ok=True)
    # chr22 up to the last stretch, N everywhere else, 60 bases a line.
    length = max(s + len(q) for s, q in REF.items()) + 50
    seq = bytearray(b"N" * length)
    for start, q in REF.items():
        seq[start - 1:start - 1 + len(q)] = q.encode()
    fa = os.path.join(gd, "reference", "GRCh38_no_alt_analysis_set.fasta")
    with open(fa, "wb") as f:
        f.write(b">chr22\n")
        for i in range(0, length, 60):
            f.write(bytes(seq[i:i + 60]) + b"\n")
    with open(fa + ".fai", "w") as f:
        f.write(f"chr22\t{length}\t7\t60\t61\n")
    shutil.copy(os.path.join(FIXTURES, "combined_allele_sample.vcf"), os.path.join(gd, "S1", "vcf", "S1.vcf"))
    shutil.copy(os.path.join(FIXTURES, "combined_allele_clinvar.vcf"),
                os.path.join(gd, "clinvar", "clinvar_pathogenic_chr.vcf"))
    uid = f"{os.getuid()}:{os.getgid()}"
    p = subprocess.run(["docker", "run", "--rm", "--user", uid, "-v", f"{gd}:/g", image, "sh", "-c",
                        "set -e; bcftools view -Oz -o /g/S1/vcf/S1.vcf.gz /g/S1/vcf/S1.vcf;"
                        " bcftools index -t /g/S1/vcf/S1.vcf.gz;"
                        " bcftools view -Oz -o /g/clinvar/clinvar_pathogenic_chr.vcf.gz /g/clinvar/clinvar_pathogenic_chr.vcf;"
                        " bcftools index -t /g/clinvar/clinvar_pathogenic_chr.vcf.gz"],
                       capture_output=True, text=True)
    check("step 06: the fixture files compress and index", p.returncode == 0, p.stderr)
    env["GENOME_DIR"] = gd
    p = subprocess.run(["bash", os.path.join(REPO, "scripts", "06-clinvar-screen.sh"), "S1"],
                       capture_output=True, text=True, env=env)
    log = p.stdout + p.stderr
    print("    | " + log.rstrip("\n").replace("\n", "\n    | "))
    check("step 06: exits 0", p.returncode == 0, f"exit {p.returncode}")
    hits_path = os.path.join(gd, "S1", "clinvar", "S1_clinvar_hits.vcf")
    rows = []
    if os.path.isfile(hits_path):
        with open(hits_path) as f:
            rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    by = {(r[1], r[3], r[4]): r for r in rows}
    comb = by.get(("30000038", "GCCAT", "CATGGG"))
    check("step 06: the 2-star ClinVar record of the combined allele is a hit",
          comb is not None and comb[2] == "9900003"
          and "CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts" in comb[7]
          and "PGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG" in comb[7], rows)
    parts = [by.get(("30000037", "AGC", "A")), by.get(("30000042", "T", "TGGG"))]
    check("step 06: both 1-star parts are hits, marked part of the combined allele",
          all(x is not None and "PGP_PART_OF=chr22:30000038:GCCAT:CATGGG" in x[7] for x in parts), rows)
    snv = by.get(("30000210", "C", "T"))
    check("step 06: the single SNV hit is unchanged (its ClinVar fields, no combined-allele tag)",
          snv is not None and snv[2] == "9900004" and "PGP_" not in snv[7], rows)
    check("step 06: four hits in all", len(rows) == 4, rows)
    tsv = os.path.join(gd, "S1", "clinvar", "S1_clinvar_hits.tsv")
    n_tsv = sum(1 for _ in open(tsv)) - 1 if os.path.isfile(tsv) else -1
    check("step 06: the hits TSV has the same four rows", n_tsv == 4, n_tsv)
    check("step 06: the log says which hits are parts of the combined allele",
          "part of a combined allele: chr22:30000038 GCCAT>CATGGG" in log
          and "combined allele of chr22:30000037 AGC>A + chr22:30000042 T>TGGG" in log, "")
    # Files written by the containers belong to the caller (run_in --user).
    subprocess.run(["docker", "run", "--rm", "-v", f"{work}:/w", image, "sh", "-c", "rm -rf /w/genome"],
                   capture_output=True)


def main():
    work = tempfile.mkdtemp()
    try:
        tools = awks()
        check("at least one awk is installed", len(tools) >= 1, "no awk on PATH")
        for name, path in tools:
            awk_checks(name, path, work)
        # Check 8 pulls and runs BCFTOOLS_IMAGE: always in CI, elsewhere only
        # when asked for, so a quick local run stays light.
        if os.environ.get("CI") or os.environ.get("PGP_TEST_DOCKER") == "1":
            if shutil.which("docker"):
                docker_checks(work)
            else:
                check("docker is installed (check 8 runs the real step 06)", False, "no docker on PATH")
        else:
            print("[SKIP] check 8: set PGP_TEST_DOCKER=1 to run the real step 06 here (CI always runs it)")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print()
    print("RESULT: ALL PASS" if FAILS == 0 else f"RESULT: {FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
