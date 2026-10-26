#!/usr/bin/env python3
"""y_tree_place.py: place a Y chromosome on a current public tree (step 37).

Yleaf 3.2.1, the only Yleaf image, places the sample on the YFull tree it
ships (v10.01). The current YFull YTree is several releases newer, so the
terminal branch Yleaf reports can sit levels above the one the current tree
resolves. This tool walks the current tree itself, from the reads, and step 37
reports both results, each with the tree it comes from.

The tree folder (`setup.sh --yleaf-data` installs it as <yleaf data>/ytree):
  current_tree.json       YFull YTree (YTREE_VERSION, CC BY 4.0): every branch
                          with its id and its defining SNPs by name
                          ("A/B" = one SNP with two names)
  positions_ftdna_hg38.txt, positions_yfull_hg38.txt
                          Yleaf 4.x's GRCh38 SNP tables (GPL-3.0): name,
                          position, ancestral and derived base. A SNP name
                          is looked up in the first table, then the second.
  version.txt             "tree<TAB>..." and "positions<TAB>..." lines

  y_tree_place.py positions --tree-dir DIR OUT
      the "chrY<TAB>pos" list of every tree SNP the tables place outside
      the pseudoautosomal regions, for `samtools mpileup -l`
  y_tree_place.py place --tree-dir DIR --pileup PILEUP --sample S --out-prefix P [--min-depth N]
      reads the pileup (samtools mpileup without -f, as step 37 makes it) and
      writes P_y_tree.txt (one result row) and P_y_tree_path.tsv (one row per
      branch from the root to the terminal one)

A SNP is called when at least --min-depth reads (default 10) cover it and 90%
of them carry one of its two bases: derived (D) or ancestral (A). Anything
else is a no-call (N), and a SNP no table places is "not on the tables".
These are Yleaf's own defaults (-r 10, -b 90).

A branch is positive when it has at least one derived SNP and more derived
than ancestral ones; it is negative when it has an ancestral SNP and no more
derived than ancestral ones. The terminal branch is the deepest positive
branch with no negative branch above it (ties: more derived SNPs on the
path, then the name). Branches with no call at all neither help nor block,
so a branch whose SNPs no table places does not stop the walk below it.
The ancestral SNPs inside the terminal branch's own block are listed: the
sample shares only part of that block, a sign of a sub-branch the tree does
not have yet. Standard library only.
"""
import json
import re
import sys
from pathlib import Path

MIN_DEPTH = 10      # Yleaf's -r default
MAJORITY = 0.9      # Yleaf's -b default (90%)
PAR1_END = 2781479  # GRCh38 chrY PAR1 is 10,001-2,781,479
PAR2_START = 56887903  # GRCh38 chrY PAR2 is 56,887,903-57,217,415
TABLES = ("positions_ftdna_hg38.txt", "positions_yfull_hg38.txt")
TREE = "current_tree.json"
BASES = set("ACGT")


def snp_groups(node):
    """The defining SNPs of a branch: a list of name lists (aliases)."""
    out = []
    for g in (node.get("snps") or "").split(","):
        names = [a.strip() for a in g.split("/") if a.strip()]
        if names:
            out.append(names)
    return out


def walk(tree):
    """(node, parent_index, depth) for every branch, root first."""
    out, stack = [], [(tree, -1, 0)]
    while stack:
        node, parent, depth = stack.pop()
        out.append((node, parent, depth))
        i = len(out) - 1
        for child in reversed(node.get("children") or []):
            stack.append((child, i, depth + 1))
    return out


def load(tree_dir):
    d = Path(tree_dir)
    for need in (d / TREE, *(d / t for t in TABLES)):
        if not need.is_file():
            sys.exit(f"ERROR: {need} not found: install the current Y tree with scripts/setup.sh --yleaf-data <genome_dir>")
    tree = json.loads((d / TREE).read_text())
    nodes = walk(tree)
    wanted = {n for node, _, _ in nodes for g in snp_groups(node) for n in g}
    table = {}
    for t in TABLES:
        seen = {}
        with open(d / t) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 7 or c[1] not in wanted or c[1] in table or c[1] in seen:
                    continue
                if not c[3].isdigit() or c[5] not in BASES or c[6] not in BASES or c[5] == c[6]:
                    continue
                pos = int(c[3])
                if pos <= PAR1_END or pos >= PAR2_START:
                    continue
                seen[c[1]] = (pos, c[5], c[6])
        table.update(seen)
    labels = {}
    vf = d / "version.txt"
    if vf.is_file():
        for line in vf.read_text().splitlines():
            k, _, v = line.partition("\t")
            labels[k.strip()] = v.strip()
    return nodes, table, labels


def place_snp(names, table):
    for n in names:
        if n in table:
            return n, table[n]
    return names[0], None


def positions(tree_dir, out):
    nodes, table, _ = load(tree_dir)
    pos = sorted({table[n][0] for n in table})
    if not pos:
        sys.exit(f"ERROR: no tree SNP is on the position tables in {tree_dir}")
    with open(out, "w") as o:
        for p in pos:
            o.write(f"chrY\t{p}\n")
    print(f"{len(pos)} Y positions of the current tree ({len(nodes)} branches)")


INDEL = re.compile(r"[+-](\d+)")


def base_counts(bases, ref):
    """Counts of A, C, G, T in a samtools mpileup read-base column."""
    counts = dict.fromkeys("ACGT", 0)
    i, n = 0, len(bases)
    ref = ref.upper()
    while i < n:
        ch = bases[i]
        if ch == "^":
            i += 2
            continue
        if ch in "+-":
            m = INDEL.match(bases, i)
            if m:
                i = m.end() + int(m.group(1))
                continue
        if ch in ".,":
            ch = ref
        ch = ch.upper()
        if ch in counts:
            counts[ch] += 1
        i += 1
    return counts


def read_pileup(path):
    calls = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or not c[1].isdigit():
                continue
            calls[int(c[1])] = base_counts(c[4], c[2])
    return calls


def call(snp, pileup, min_depth):
    pos, anc, der = snp
    counts = pileup.get(pos)
    if not counts:
        return "N"
    depth = sum(counts.values())
    if depth < min_depth:
        return "N"
    if counts[der] >= MAJORITY * depth:
        return "D"
    if counts[anc] >= MAJORITY * depth:
        return "A"
    return "N"


def place(tree_dir, pileup_path, sample, prefix, min_depth=MIN_DEPTH):
    nodes, table, labels = load(tree_dir)
    pileup = read_pileup(pileup_path)
    stats = []
    for node, _, _ in nodes:
        s = {"snps": 0, "on_table": 0, "D": 0, "A": 0, "N": 0, "anc": [], "der": []}
        for names in snp_groups(node):
            s["snps"] += 1
            name, snp = place_snp(names, table)
            if snp is None:
                continue
            s["on_table"] += 1
            c = call(snp, pileup, min_depth)
            s[c] += 1
            if c == "A":
                s["anc"].append(name)
            elif c == "D":
                s["der"].append(name)
        stats.append(s)
    positive = [s["D"] >= 1 and s["D"] > s["A"] for s in stats]
    negative = [s["A"] >= 1 and s["D"] <= s["A"] for s in stats]
    # blocked[i]: a negative branch on the path from the root to i, i included
    blocked, path_d = [False] * len(nodes), [0] * len(nodes)
    best = None
    for i, (node, parent, depth) in enumerate(nodes):
        blocked[i] = negative[i] or (parent >= 0 and blocked[parent])
        path_d[i] = stats[i]["D"] + (path_d[parent] if parent >= 0 else 0)
        if positive[i] and not blocked[i]:
            key = (depth, path_d[i], node.get("id") or "")
            if best is None or key[:2] > best[0][:2] or (key[:2] == best[0][:2] and key[2] < best[0][2]):
                best = (key, i)
    path = []
    if best is not None:
        i = best[1]
        while i >= 0:
            if nodes[i][0].get("id"):
                path.append(i)
            i = nodes[i][1]
        path.reverse()
    tree_label = labels.get("tree", "current tree")
    pos_label = labels.get("positions", "position tables")
    called = sum(s["D"] + s["A"] for s in stats)
    term = path[-1] if path else None
    hg = nodes[term][0]["id"] if term is not None else "NA"
    with open(f"{prefix}_y_tree_path.tsv", "w") as o:
        o.write("Level\tBranch\tSNPs\tOn_tables\tDerived\tAncestral\tNo_call\tAncestral_SNPs\n")
        for lvl, i in enumerate(path, 1):
            s = stats[i]
            o.write(f"{lvl}\t{nodes[i][0]['id']}\t{s['snps']}\t{s['on_table']}\t{s['D']}\t{s['A']}\t{s['N']}\t"
                    f"{','.join(s['anc']) or '.'}\n")
    pd = sum(stats[i]["D"] for i in path)
    pa = sum(stats[i]["A"] for i in path)
    pn = sum(stats[i]["N"] for i in path)
    anc_term = ",".join(stats[term]["anc"]) if term is not None and stats[term]["anc"] else "."
    with open(f"{prefix}_y_tree.txt", "w") as o:
        o.write("Sample_name\tHg\tLevel\tPath_derived\tPath_ancestral\tPath_no_call\tAncestral_in_terminal\t"
                "SNPs_called\tMin_depth\tTree\tPositions\n")
        o.write(f"{sample}\t{hg}\t{len(path)}\t{pd}\t{pa}\t{pn}\t{anc_term}\t{called}\t{min_depth}\t"
                f"{tree_label}\t{pos_label}\n")
    if term is None:
        print(f"Current tree ({tree_label}): no branch has a derived call ({called} SNPs called at depth >= {min_depth})")
    else:
        extra = f"; ancestral inside its block: {anc_term}" if anc_term != "." else ""
        print(f"Current tree ({tree_label}): {hg} (level {len(path)}, path SNPs {pd} derived, {pa} ancestral, "
              f"{pn} no call{extra})")
    return hg


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(prog="y_tree_place.py", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("positions")
    p.add_argument("--tree-dir", required=True)
    p.add_argument("out")
    q = sub.add_parser("place")
    for a in ("--tree-dir", "--pileup", "--sample", "--out-prefix"):
        q.add_argument(a, required=True)
    q.add_argument("--min-depth", type=int, default=MIN_DEPTH)
    a = ap.parse_args(argv)
    if a.cmd == "positions":
        positions(a.tree_dir, a.out)
    else:
        if a.min_depth < 1:
            sys.exit("ERROR: --min-depth must be 1 or more")
        place(a.tree_dir, a.pileup, a.sample, a.out_prefix, a.min_depth)


if __name__ == "__main__":
    main(sys.argv[1:])
