#!/usr/bin/env python3
"""bin/pgx_refblock_fill.py: PharmCAT's multi-base positions called 0/0 from
the reference blocks of a gVCF (steps 07 and PHARMCAT_PREPROCESS run it).

The fixture gVCF (tests/fixtures/pharmcat/refblocks.g.vcf) has reference
blocks over nine PharmCAT positions longer than one base, from the subset of
PharmCAT 3.4.0's positions file beside it. The test expands the blocks into
one-base records as `bcftools convert --gvcf2vcf | bcftools view
--trim-alt-alleles -i 'GT!="mis"'` does, runs the script on that stream, and
checks:
  1. five positions are called 0/0, with PharmCAT's ID, REF and ALT, the flag
     PGX_REFBLOCK, the block's genotype (0 on haploid chrX) and the lowest GQ
     over the span: DPYD (two blocks meet inside the span), UGT1A1 *28, a TPMT
     insertion, a RYR1 two-base substitution and a G6PD deletion;
  2. four stay missing: a ./. hole inside the span (CYP2C9, 11 bp), a block
     with GQ 10, a variant record inside the span, a span no block reaches;
  3. each filled record takes the place of the one-base record at its
     position, every other line passes unchanged, and the output stays sorted;
  4. the GQ threshold: 20 by default, a block at exactly --min-gq counts;
  5. a block that stops one base short, a block with a failed FILTER, a ./.
     block with a high GQ, or a deletion that reaches into the span from
     before it, fills nothing;
     a gVCF with two samples or no #CHROM line is an error (exit 2).

Run: python3 tests/test_pgx_refblock_fill.py
"""
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "pgx_refblock_fill.py")
FX = os.path.join(REPO, "tests", "fixtures", "pharmcat")
GVCF = os.path.join(FX, "refblocks.g.vcf")
POSITIONS = os.path.join(FX, "pharmcat_positions-3.4.0-subset.vcf")
FAILS = []

FILLED = {  # (chrom, pos): (REF, ALT, GT, GQ, ID)
    ("chr1", 97450065): ("TG", "T", "0/0", "45", "rs72549303"),
    ("chr2", 233760233): ("CAT", "C,CATAT,CATATAT", "0/0", "50", "rs3064744"),
    ("chr6", 18149032): ("C", "CT", "0/0", "40", "rs759836180"),
    ("chr19", 38572266): ("TC", "CT", "0/0", "35", "rs193922862"),
    ("chrX", 154532055): ("CTCT", "C", "0", "30", "."),
}
LEFT = {
    ("chr10", 94942205): "a ./. hole inside the REF span",
    ("chr10", 94949161): "a block with GQ 10",
    ("chr10", 94949281): "a variant record inside the REF span",
    ("chr6", 18149103): "no block over the span",
}


def check(desc, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS.append(desc)


def expand(gvcf_text):
    """What gvcf2vcf plus `view --trim-alt-alleles -i 'GT!="mis"'` writes:
    one record per base of each block (REF from the reference: N here, the
    script reads only its length), ALT . for a reference call, no ./. records."""
    out = []
    for line in gvcf_text.splitlines(keepends=True):
        if line.startswith("#"):
            out.append(line)
            continue
        c = line.rstrip("\n").split("\t")
        gt = c[9].split(":")[0]
        if "." in re.split(r"[/|]", gt):
            continue
        alts = [a for a in c[4].split(",") if a != "<*>"]
        m = re.search(r"END=(\d+)", c[7])
        if m and not alts:
            for p in range(int(c[1]), int(m.group(1)) + 1):
                out.append("\t".join([c[0], str(p), ".", "N", ".", c[5], c[6], ".", c[8], c[9]]) + "\n")
        else:
            out.append("\t".join(c[:4] + [",".join(alts) or "."] + c[5:]) + "\n")
    return "".join(out)


def run(gvcf_text, *opts, positions=POSITIONS, stdin=None):
    with tempfile.TemporaryDirectory() as tmp:
        g = os.path.join(tmp, "records.vcf")
        with open(g, "w") as f:
            f.write(gvcf_text)
        p = subprocess.run([sys.executable, SCRIPT, *opts, positions, g], text=True, capture_output=True,
                           input=expand(gvcf_text) if stdin is None else stdin)
    return p.returncode, p.stdout, p.stderr


def records(vcf):
    return [l.rstrip("\n").split("\t") for l in vcf.splitlines() if l and not l.startswith("#")]


def fills(vcf):
    return {(r[0], int(r[1])): r for r in records(vcf) if r[7] == "PGX_REFBLOCK"}


def main():
    with open(GVCF) as f:
        gvcf = f.read()
    rc, out, err = run(gvcf)
    print(err.strip())
    check("runs and exits 0", rc == 0, err)
    got = fills(out)
    # 1. filled
    check("five positions called 0/0: DPYD, UGT1A1 *28, TPMT, RYR1, G6PD",
          sorted(got) == sorted(FILLED), sorted(got))
    for key, (ref, alt, gt, gq, vid) in FILLED.items():
        r = got.get(key)
        check(f"{key[0]}:{key[1]}: PharmCAT's ID, REF and ALT, GT {gt}, lowest GQ {gq}",
              r is not None and r[2:5] == [vid, ref, alt] and r[6] == "PASS" and r[8] == "GT:GQ"
              and r[9] == f"{gt}:{gq}", r)
    # 2. left missing
    for key, why in LEFT.items():
        check(f"{key[0]}:{key[1]} stays missing ({why})", key not in got, got.get(key))
    check("the summary on stderr counts 5 of the 9 positions longer than one base",
          "5 of 9 PharmCAT positions longer than one base called 0/0 from reference blocks with GQ >= 20" in err, err)
    # 3. replaced in place, the rest unchanged, sorted
    inp = records(expand(gvcf))
    outr = records(out)
    at = {}
    for r in outr:
        at.setdefault((r[0], int(r[1])), []).append(r)
    check("one record at each filled position (the one-base record is replaced)",
          all(len(at.get(k, [])) == 1 for k in FILLED), {k: len(at.get(k, [])) for k in FILLED})
    kept = [r for r in inp if (r[0], int(r[1])) not in FILLED]
    check("every other record passes unchanged and in order",
          [r for r in outr if r[7] != "PGX_REFBLOCK"] == kept, len(kept))
    check("the bases after a filled position inside its span stay (chr2:233760234-233760235)",
          ("chr2", 233760234) in at and ("chr2", 233760235) in at)
    order = ["chr1", "chr2", "chr6", "chr10", "chr19", "chrX"]
    keys = [(order.index(r[0]), int(r[1])) for r in outr]
    check("the output is sorted", keys == sorted(keys))
    head = [l for l in out.splitlines() if l.startswith("#")]
    check("the header declares PGX_REFBLOCK before #CHROM",
          len(head) >= 2 and head[-1].startswith("#CHROM") and head[-2].startswith("##INFO=<ID=PGX_REFBLOCK,"), head[-3:])
    check("the sample column stays", head and head[-1].split("\t")[9:] == ["REFBLOCKS"], head[-1:])
    # 4. threshold
    _, o, _ = run(gvcf, "--min-gq", "50")
    g50 = fills(o)
    check("--min-gq 50: UGT1A1 (GQ 50) still filled, a block at exactly the threshold counts",
          ("chr2", 233760233) in g50, sorted(g50))
    check("--min-gq 50: DPYD (one block GQ 45) no longer filled", ("chr1", 97450065) not in g50, sorted(g50))
    _, o, _ = run(gvcf, "--min-gq", "51")
    check("--min-gq 51: nothing filled", fills(o) == {}, sorted(fills(o)))
    low = gvcf.replace("0/0:50:30\n", "0/0:19:30\n", 1)  # the first GQ-50 block of a filled position: UGT1A1
    _, o, _ = run(low)
    check("default threshold 20: UGT1A1 with a GQ 19 block is not filled", ("chr2", 233760233) not in fills(o))
    # 5. edges
    hole = gvcf.replace("./.:0:0", "./.:50:0")
    _, o, _ = run(hole)
    check("the ./. hole stays a hole when its GQ is 50 (the genotype decides, not only the GQ)",
          hole != gvcf and ("chr10", 94942205) not in fills(o), sorted(fills(o)))
    # A deletion that starts before the span and reaches into it, with a
    # block over the deleted bases after it (as GATK writes them): the block
    # covers the span, the deletion still overlaps it.
    dele = gvcf.replace(
        "chr2\t233759900\t.\tT\t<*>\t0\t.\tEND=233761000\tGT:GQ:MIN_DP\t0/0:50:30\n",
        "chr2\t233759900\t.\tT\t<*>\t0\t.\tEND=233760230\tGT:GQ:MIN_DP\t0/0:50:30\n"
        "chr2\t233760231\t.\tGCCA\tG,<*>\t30\tPASS\t.\tGT:GQ:DP\t0/1:30:28\n"
        "chr2\t233760232\t.\tC\t<*>\t0\t.\tEND=233761000\tGT:GQ:MIN_DP\t0/0:50:30\n")
    _, o, _ = run(dele)
    check("a deletion from before the span that overlaps it fills nothing, though blocks cover the span",
          dele != gvcf and ("chr2", 233760233) not in fills(o), sorted(fills(o)))
    short = gvcf.replace("END=233761000", "END=233760234")
    _, o, _ = run(short)
    check("a block ending one base before the span ends fills nothing", ("chr2", 233760233) not in fills(o))
    lowq = gvcf.replace("T\t<*>\t0\t.\tEND=233761000", "T\t<*>\t0\tLowQual\tEND=233761000")
    _, o, _ = run(lowq)
    check("a block with a failed FILTER fills nothing", lowq != gvcf and ("chr2", 233760233) not in fills(o))
    two = gvcf.replace("\tREFBLOCKS\n", "\tREFBLOCKS\tOTHER\n")
    rc2, _, e2 = run(two)
    check("a gVCF with two samples: exit 2", rc2 == 2 and "one sample" in e2, (rc2, e2))
    nohead = "".join(l for l in gvcf.splitlines(keepends=True) if not l.startswith("#CHROM"))
    rc3, _, e3 = run(nohead, stdin=expand(gvcf))
    check("a gVCF without its #CHROM line: exit 2", rc3 == 2 and "#CHROM" in e3, (rc3, e3))

    print("\nRESULT:", "ALL PASS" if not FAILS else f"{len(FAILS)} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
