#!/usr/bin/env python3
"""bin/y_tree_place.py on a synthetic tree, synthetic position tables and
synthetic pileups (no real sample).

The tree:
    (root) -> A (a1, a2/a2alt) -> B (b1, b2, b3) -> C (c1, c2)
                               -> X (x1)        -> Y (y1)
Checks:
  1. every SNP derived along A-B-C: C is the terminal branch at level 3, with
     the counts of each level, and nothing ancestral in its block;
  2. one ancestral SNP inside B (b3) and no read at C: B is the terminal
     branch, and b3 is listed as ancestral inside the terminal block;
  3. no-calls: C's SNPs read at depth 9 (below 10), or by a 50/50 mix of
     bases, are no-calls, not derived; a SNP no table places counts as a
     SNP but not as on the tables;
  4. a negative branch blocks the branches below it: X ancestral, Y derived
     (deeper than B) is not the terminal branch; without X's ancestral call
     Y wins, because it is deeper;
  5. nothing derived anywhere: Hg NA, level 0;
  6. the first table wins for a name both tables place, and a SNP in the
     pseudoautosomal region is never asked for;
  7. the pileup reader skips read starts, ends and indels, and reads '.' and
     ',' as the reference base.

Run: python3 tests/test_y_tree_place.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(REPO, "bin", "y_tree_place.py")
sys.path.insert(0, os.path.join(REPO, "bin"))

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


TREE = {"id": "", "snps": "", "children": [
    {"id": "A", "snps": "a1, a2/a2alt", "children": [
        {"id": "B", "snps": "b1, b2, b3, bnot", "children": [
            {"id": "C", "snps": "c1, c2", "children": []}]},
        {"id": "X", "snps": "x1", "children": [
            {"id": "Y", "snps": "y1", "children": [
                {"id": "Z", "snps": "z1", "children": []}]}]}]}]}

# name -> (position, ancestral, derived); bnot is on no table
FTDNA = {"a1": (3000100, "C", "T"), "a2alt": (3000200, "G", "A"), "b1": (3000300, "A", "G"),
         "b2": (3000400, "T", "C"), "b3": (3000500, "C", "G"), "c1": (3000600, "G", "T"),
         "x1": (3000800, "A", "C"), "y1": (3000900, "T", "A"), "z1": (3001000, "G", "C"),
         "par": (2000000, "A", "G")}
YFULL = {"a1": (3999999, "C", "T"), "c2": (3000700, "A", "T")}


def write_tables(d):
    def table(path, rows):
        with open(path, "w") as f:
            for name, (pos, anc, der) in rows.items():
                f.write(f"chry\t{name}\tR-{name}\t{pos}\t{anc}->{der}\t{anc}\t{der}\n")
    tree = json.loads(json.dumps(TREE))
    tree["children"][0]["snps"] += ", par"
    with open(os.path.join(d, "current_tree.json"), "w") as f:
        json.dump(tree, f)
    table(os.path.join(d, "positions_ftdna_hg38.txt"), FTDNA)
    table(os.path.join(d, "positions_yfull_hg38.txt"), YFULL)
    with open(os.path.join(d, "version.txt"), "w") as f:
        f.write("tree\tYFull YTree test\npositions\ttest tables\n")


def pos_of(name):
    return (FTDNA.get(name) or YFULL[name])


def pileup(path, states, depth=12):
    """states: name -> 'D', 'A', 'mix' or ('D', depth)."""
    with open(path, "w") as f:
        for name, st in states.items():
            dp = depth
            if isinstance(st, tuple):
                st, dp = st
            pos, anc, der = pos_of(name)
            if st == "mix":
                bases = anc * (dp // 2) + der * (dp - dp // 2)
            else:
                bases = (der if st == "D" else anc) * dp
            f.write(f"chrY\t{pos}\tN\t{dp}\t{bases}\t{'I' * dp}\n")


def run(d, states, tag, depth=12, extra=()):
    p = os.path.join(d, f"{tag}.pileup")
    pileup(p, states, depth)
    prefix = os.path.join(d, tag)
    r = subprocess.run([sys.executable, TOOL, "place", "--tree-dir", d, "--pileup", p, "--sample", "S",
                        "--out-prefix", prefix, *extra], capture_output=True, text=True)
    if r.returncode != 0:
        return r, {}, []
    with open(prefix + "_y_tree.txt") as f:
        lines = f.read().splitlines()
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    with open(prefix + "_y_tree_path.tsv") as f:
        pl = f.read().splitlines()
    head = pl[0].split("\t")
    path = [dict(zip(head, x.split("\t"))) for x in pl[1:]]
    return r, row, path


def main():
    d = tempfile.mkdtemp()
    try:
        write_tables(d)
        ALL_A = {"a1": "D", "a2alt": "D"}
        ALL_B = {"b1": "D", "b2": "D", "b3": "D"}

        # 1. deepest branch
        r, row, path = run(d, {**ALL_A, **ALL_B, "c1": "D", "c2": "D", "x1": "A"}, "deep")
        check("1. the tool runs", r.returncode == 0, r.stderr)
        check("1. C is the terminal branch, at level 3", row.get("Hg") == "C" and row.get("Level") == "3", row)
        check("1. the path is A, B, C", [x["Branch"] for x in path] == ["A", "B", "C"], path)
        check("1. B's counts: 4 SNPs, 3 on the tables, 3 derived, 0 ancestral",
              path[1:2] and (path[1]["SNPs"], path[1]["On_tables"], path[1]["Derived"], path[1]["Ancestral"]) == ("4", "3", "3", "0"),
              path)
        check("1. nothing ancestral in the terminal block", row.get("Ancestral_in_terminal") == ".", row)
        check("1. the path totals: 7 derived, 0 ancestral",
              (row.get("Path_derived"), row.get("Path_ancestral")) == ("7", "0"), row)
        check("1. the result names its tree and tables",
              row.get("Tree") == "YFull YTree test" and row.get("Positions") == "test tables", row)
        check("1. the step line names the current tree", "Current tree (YFull YTree test): C (level 3" in r.stdout, r.stdout)

        # 2. one ancestral SNP in the terminal block
        r, row, path = run(d, {**ALL_A, "b1": "D", "b2": "D", "b3": "A", "x1": "A"}, "split")
        check("2. B is the terminal branch with one ancestral SNP in its block",
              row.get("Hg") == "B" and row.get("Ancestral_in_terminal") == "b3", row)
        check("2. B's row lists b3 as ancestral", path[-1:] and path[-1]["Ancestral_SNPs"] == "b3"
              and path[-1]["Derived"] == "2" and path[-1]["Ancestral"] == "1", path)

        # 3. no-calls
        r, row, path = run(d, {**ALL_A, **ALL_B, "c1": ("D", 9), "c2": "mix", "x1": "A"}, "nocall")
        check("3. depth 9 and a 50/50 mix are no-calls: B is the terminal branch",
              row.get("Hg") == "B", row)
        r2, row2, path2 = run(d, {**ALL_A, **ALL_B, "c1": ("D", 9), "c2": "mix", "x1": "A"}, "nocall2",
                              extra=("--min-depth", "5"))
        check("3. with --min-depth 5 the depth-9 read is derived: C, with one no-call",
              row2.get("Hg") == "C" and path2[-1:] and path2[-1]["No_call"] == "1", (row2, path2))
        check("3. a SNP no table places is counted but not on the tables (bnot)",
              path[1:2] and path[1]["SNPs"] == "4" and path[1]["On_tables"] == "3", path)

        # 4. a negative branch blocks the walk below it
        r, row, path = run(d, {**ALL_A, **ALL_B, "x1": "A", "y1": "D", "z1": "D"}, "blocked")
        check("4. X ancestral: Z (deeper, derived) is not the terminal branch, B is",
              row.get("Hg") == "B", row)
        r, row, path = run(d, {**ALL_A, **ALL_B, "y1": "D", "z1": "D"}, "unblocked")
        check("4. X with no call: Z, the deeper branch, wins", row.get("Hg") == "Z"
              and [x["Branch"] for x in path] == ["A", "X", "Y", "Z"], (row, path))

        # 5. nothing derived
        r, row, path = run(d, {"a1": "A", "a2alt": "A"}, "none")
        check("5. no derived call: Hg NA at level 0, an empty path",
              row.get("Hg") == "NA" and row.get("Level") == "0" and path == [], (row, path))
        check("5. the step line says no branch has a derived call", "no branch has a derived call" in r.stdout, r.stdout)

        # 6. tables
        out = os.path.join(d, "pos.txt")
        r = subprocess.run([sys.executable, TOOL, "positions", "--tree-dir", d, out], capture_output=True, text=True)
        pos = open(out).read().split() if r.returncode == 0 else []
        check("6. a1 is asked for at the first table's position, not the second's",
              "3000100" in pos and "3999999" not in pos, pos)
        check("6. c2, on the second table only, is asked for", "3000700" in pos, pos)
        check("6. the pseudoautosomal SNP is never asked for", "2000000" not in pos, pos)
        check("6. every line is chrY<TAB>pos", r.returncode == 0 and all(
            x.split("\t")[0] == "chrY" for x in open(out).read().splitlines()), r.stderr)

        # 7. pileup reader
        import y_tree_place as y
        c = y.base_counts("^I..,$+2ACt-1Ga*", "t")
        check("7. read marks and indels skipped, '.'/',' read as the reference",
              c == {"A": 1, "C": 0, "G": 0, "T": 4}, c)

        # Missing data is an error that says how to install it
        os.remove(os.path.join(d, "positions_yfull_hg38.txt"))
        r = subprocess.run([sys.executable, TOOL, "positions", "--tree-dir", d, out], capture_output=True, text=True)
        check("a missing table stops the tool and names setup.sh --yleaf-data",
              r.returncode != 0 and "setup.sh --yleaf-data" in r.stderr, r.stderr)
    finally:
        shutil.rmtree(d)
    print(f"\n{'FAILED' if FAILS else 'OK'}: {FAILS} failed check(s)")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
