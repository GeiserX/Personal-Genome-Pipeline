#!/usr/bin/env python3
"""bin/constraint_join.awk, the gnomAD constraint loader of steps 23 and 31 and
the Nextflow SLIVAR_PRIORITIZE module, on a six-row constraint table.

Checks:
  1. mis_z is gnomAD v4.1's `mis.z_score` (the old loaders read `mis_z` or
     `missense.z_score`, so mis_z was always '.');
  2. only canonical rows count: a non-canonical Ensembl row listed before the
     canonical one does not win (the old step 31 had no canonical filter);
  3. of a gene's Ensembl and RefSeq canonical rows the Ensembl one wins,
     before (GENEB) or after (GENEC) the RefSeq row;
  4. a table where no gene matches exits 4, a table without mis.z_score exits 3;
  5. every awk on the machine (awk, mawk, gawk, busybox awk) gives the same
     output: the module runs it with the bcftools image's awk, the scripts with
     the host's.

Run: python3 tests/test_constraint_join.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AWK_FILE = os.path.join(REPO, "bin", "constraint_join.awk")

HEADER = ["gene", "gene_id", "transcript", "canonical", "mane_select",
          "lof.oe_ci.upper", "lof.pLI", "mis.z_score", "syn.z_score"]
ROWS = [
    ["GENEA", "ENSG1", "ENST00000000002", "false", "false", "0.95", "0.01", "0.40", "0.2"],
    ["GENEA", "ENSG1", "ENST00000000001", "true", "true", "0.21", "0.99", "3.52", "0.1"],
    ["GENEB", "1001", "NM_000003.4", "true", "false", "0.80", "0.00", "1.10", "0.3"],
    ["GENEB", "ENSG2", "ENST00000000004", "true", "true", "0.52", "0.20", "2.25", "0.4"],
    ["GENEC", "ENSG3", "ENST00000000005", "true", "true", "NA", "NA", "NA", "0.5"],
    ["GENEC", "1002", "NM_000006.2", "true", "false", "0.70", "0.10", "1.50", "0.6"],
]
TABLE = [["CHROM", "POS", "GENE"],
         ["chr1", "100", "GENEA"],
         ["chr2", "200", "GENEB"],
         ["chr3", "300", "GENEC"],
         ["chr4", "400", "GENED"],
         ["chr5", "500", "."]]

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + detail}")
    if not ok:
        FAILS += 1


def write(path, rows):
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(r) + "\n")


def awks():
    found = []
    for name in ("awk", "mawk", "gawk", "original-awk", "nawk"):
        p = shutil.which(name)
        if p and os.path.realpath(p) not in [os.path.realpath(x[1]) for x in found]:
            found.append((name, p))
    if shutil.which("busybox"):
        found.append(("busybox awk", "busybox"))
    return found


def run(awk, constraint, table, *extra):
    cmd = ([awk] if awk != "busybox" else ["busybox", "awk"]) + \
        ["-f", AWK_FILE, "gene_col=GENE", *extra, constraint, table]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def parse(out):
    lines = out.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    return head, {r.split("\t")[2]: dict(zip(head, r.split("\t"))) for r in lines[1:]}


def is_number(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def main():
    work = tempfile.mkdtemp()
    try:
        constraint = os.path.join(work, "constraint.tsv")
        table = os.path.join(work, "table.tsv")
        write(constraint, [HEADER] + ROWS)
        write(table, TABLE)

        tools = awks()
        check("at least one awk is installed", len(tools) >= 1, "no awk on PATH")
        outputs = {}
        for name, path in tools:
            rc, out, err = run(path, constraint, table, "constrained=1")
            outputs[name] = out
            check(f"{name}: exits 0", rc == 0, f"exit {rc}: {err.strip()}")
            if rc != 0:
                continue
            head, rows = parse(out)
            check(f"{name}: header gains LOEUF, pLI, mis_z, CONSTRAINED",
                  head[-4:] == ["LOEUF", "pLI", "mis_z", "CONSTRAINED"], str(head))
            a, b, c, d, dot = (rows.get(g, {}) for g in ("GENEA", "GENEB", "GENEC", "GENED", "."))
            check(f"{name}: GENEA mis_z is numeric (mis.z_score)", is_number(a.get("mis_z", ".")), str(a))
            check(f"{name}: GENEA takes the canonical row (LOEUF 0.21, not the non-canonical 0.95)",
                  a.get("LOEUF") == "0.21" and a.get("mis_z") == "3.52", str(a))
            check(f"{name}: GENEA is CONSTRAINED (LOEUF < 0.35)", a.get("CONSTRAINED") == "YES", str(a))
            check(f"{name}: GENEB takes the Ensembl canonical row over the earlier RefSeq one",
                  b.get("LOEUF") == "0.52" and b.get("pLI") == "0.20" and b.get("mis_z") == "2.25", str(b))
            check(f"{name}: GENEB is not CONSTRAINED", b.get("CONSTRAINED") == "NO", str(b))
            check(f"{name}: GENEC keeps its Ensembl row over the later RefSeq one, NA written as '.'",
                  (c.get("LOEUF"), c.get("pLI"), c.get("mis_z")) == (".", ".", "."), str(c))
            check(f"{name}: GENED (not in the table) gets '.'", d.get("LOEUF") == "." and d.get("mis_z") == ".", str(d))
            check(f"{name}: a row without a gene gets '.'", dot.get("LOEUF") == ".", str(dot))

            rc2, out2, _ = run(path, constraint, table)
            head2 = out2.split("\n", 1)[0].split("\t")
            check(f"{name}: without constrained=1 there is no CONSTRAINED column",
                  rc2 == 0 and head2[-3:] == ["LOEUF", "pLI", "mis_z"], str(head2))

            nomatch = os.path.join(work, "nomatch.tsv")
            write(nomatch, [["CHROM", "POS", "GENE"], ["chr1", "1", "OTHER1"], ["chr1", "2", "OTHER2"]])
            rc3, _, err3 = run(path, constraint, nomatch)
            check(f"{name}: no gene matching the table exits 4", rc3 == 4, f"exit {rc3}: {err3.strip()}")

            old = os.path.join(work, "old_names.tsv")
            write(old, [[h.replace("mis.z_score", "missense.z_score") for h in HEADER]] + ROWS)
            rc4, _, err4 = run(path, old, table)
            check(f"{name}: a table without mis.z_score exits 3", rc4 == 3, f"exit {rc4}: {err4.strip()}")

        distinct = set(outputs.values())
        check(f"every awk gives the same output ({', '.join(outputs)})", len(distinct) == 1,
              "\n".join(f"--- {k}\n{v}" for k, v in outputs.items()))
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
