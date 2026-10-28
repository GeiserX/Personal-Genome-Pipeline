#!/usr/bin/env bash
# build-sv-rule-regions.sh: rebuild assets/sv_rule_regions.tsv, the regions
# where bin/annotsv_postfilter.py overrides AnnotSV's class 4 or 5.
#
# Usage: scripts/ci/build-sv-rule-regions.sh [OUT]   (default: assets/sv_rule_regions.tsv)
#
# Downloads four pinned files (about 145 MB, each checked against its sha256)
# into a temporary folder, writes OUT and deletes the downloads:
#
#   TRExplorer v2.1's GRCh38 tandem repeat catalog (broadinstitute/trexplorer-catalog)
#   GENCODE 50's basic annotation (the release scripts/lib/common.sh pins)
#   HGNC's complete set of 2026-10-06 (monthly archive)
#   the GRCh38 assembly regions report of NCBI (GCA_000001405.15)
#
# and keeps, all on GRCh38 chr1-22, X and Y:
#
#   tandem_repeat  every catalog repeat of 50 bp or more (the shortest SV)
#                  that overlaps a coding exon of a MANE Select transcript;
#                  name GENE:MOTIF (a motif over 12 bases as its length)
#   ig_tr_locus    the span of the IG and TR genes of GENCODE (gene_type IG_*
#                  or TR_*) of each of IGH, IGK, IGL, TRA (with TRD inside it),
#                  TRB and TRG: per locus the largest group of genes on one
#                  chromosome with gaps under 1 Mb, so the orphon copies on
#                  other chromosomes stay out
#   mhc            the MHC region of the assembly report
#   gene_alias     the GENCODE span of each HGNC pseudogene whose alias or
#                  previous symbol is the approved symbol of a gene with a
#                  GenCC disease curation; name ALIAS|PSEUDOGENE|its HGNC id|
#                  the HGNC id of the gene that holds ALIAS as its symbol
#
# Needs curl, python3 and sha256sum or shasum. It reads about 1 GB of text,
# so run it on a box with room, not in a unit test.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
OUT=${1:-${REPO}/assets/sv_rule_regions.tsv}

TRX_URL=https://github.com/broadinstitute/trexplorer-catalog/releases/download/v2.1/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz
TRX_SHA256=4eb0fbbd36d40a801c2d764cdda247fea8edfae1604087526459064624ec1417
GENCODE_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz
GENCODE_SHA256=5b9ad9aa71f6ab5d55123c51f74723644a895bc0cbaa9cee8fbc11bb9dec0bf7
HGNC_URL=https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-10-06.txt
HGNC_SHA256=8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d
REGIONS_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_regions.txt
REGIONS_SHA256=6b91d2c7961a322936db8fd09b7fb4a9143a590882f28c0a628a4ef26897f6c7

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# get NAME URL SHA256: download URL to WORK/NAME and stop unless it matches.
get() {
  echo "Downloading $2" >&2
  curl -fsSL --retry 3 -o "${WORK}/$1" "$2"
  local got
  got=$(sha256 "${WORK}/$1")
  if [ "$got" != "$3" ]; then
    echo "ERROR: $1 has sha256 ${got}, want $3" >&2
    exit 1
  fi
}
get trx.bed.gz "$TRX_URL" "$TRX_SHA256"
get gencode.gtf.gz "$GENCODE_URL" "$GENCODE_SHA256"
get hgnc.txt "$HGNC_URL" "$HGNC_SHA256"
get regions.txt "$REGIONS_URL" "$REGIONS_SHA256"

python3 - "$WORK" "${OUT}.part" <<'PY'
import bisect
import csv
import gzip
import re
import sys
from collections import defaultdict

work, out = sys.argv[1:3]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ORDER = {c: i for i, c in enumerate(CHROMS)}
MIN_TR = 50
LOCUS_GAP = 1_000_000
ATT = re.compile(r'(\S+) "([^"]*)"')
SRC_TR = "TRExplorer v2.1; GENCODE 50 MANE Select CDS"
SRC_IG = "GENCODE 50 IG and TR genes"
SRC_MHC = "GRCh38 assembly regions (GCA_000001405.15)"
SRC_ALIAS = "HGNC 2026-10-06; GENCODE 50"

rows = []
gene_span = {}                      # HGNC id -> (chrom, start0, end)
ig = defaultdict(list)              # locus -> [(chrom, start0, end)]
cds = defaultdict(list)             # chrom -> [(start0, end, gene)]
with gzip.open(f"{work}/gencode.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] not in ORDER:
            continue
        start, end = int(c[3]) - 1, int(c[4])
        if c[2] == "gene":
            a = dict(ATT.findall(c[8]))
            if a.get("hgnc_id"):
                gene_span[a["hgnc_id"]] = (c[0], start, end)
            name, gtype = a.get("gene_name", ""), a.get("gene_type", "")
            if re.match(r"(IG|TR)_", gtype) and re.match(r"(IG[HKL]|TR[ABDG])", name):
                locus = "TRA" if name.startswith("TRD") else name[:3]
                ig[locus].append((c[0], start, end))
        elif c[2] == "CDS" and 'tag "MANE_Select"' in c[8]:
            cds[c[0]].append((start, end, dict(ATT.findall(c[8]))["gene_name"]))

# tandem_repeat: catalog repeats of MIN_TR bp or more on a MANE Select coding exon
index = {}
for chrom, iv in cds.items():
    iv.sort()
    starts = [s for s, _, _ in iv]
    maxend, m = [], 0
    for _, e, _ in iv:
        m = max(m, e)
        maxend.append(m)
    index[chrom] = (iv, starts, maxend)
with gzip.open(f"{work}/trx.bed.gz", "rt") as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        if c[0] not in index:
            continue
        s, e = int(c[1]), int(c[2])
        if e - s < MIN_TR:
            continue
        iv, starts, maxend = index[c[0]]
        j = bisect.bisect_left(starts, e) - 1
        genes = set()
        while j >= 0 and maxend[j] > s:
            if iv[j][1] > s:
                genes.add(iv[j][2])
            j -= 1
        if genes:
            motif = c[3] if len(c[3]) <= 12 else f"{len(c[3])}bp motif"
            rows.append((c[0], s, e, "tandem_repeat", "class_3", f"{'/'.join(sorted(genes))}:{motif}", SRC_TR))

# ig_tr_locus: per locus the largest group of genes on one chromosome
for locus, genes in ig.items():
    groups = []
    for chrom, s, e in sorted(genes):
        g = groups[-1] if groups else None
        if g and g[0] == chrom and s - g[2] < LOCUS_GAP:
            g[2] = max(g[2], e)
            g[3] += 1
        else:
            groups.append([chrom, s, e, 1])
    chrom, s, e, _ = max(groups, key=lambda g: g[3])
    rows.append((chrom, s, e, "ig_tr_locus", "class_3", locus, SRC_IG))

# mhc: the assembly report's region (1-based, inclusive)
with open(f"{work}/regions.txt") as f:
    mhc = sorted({tuple(l.split("\t")[:4]) for l in f if l.startswith("MHC\t")})
if len(mhc) != 1:
    sys.exit(f"ERROR: want one MHC region in the assembly report, got {mhc}")
_, chrom, s, e = mhc[0]
rows.append((f"chr{chrom}", int(s) - 1, int(e), "mhc", "class_3", "MHC", SRC_MHC))

# gene_alias: a pseudogene that carries a disease gene's symbol as alias or previous symbol
with open(f"{work}/hgnc.txt", newline="") as f:
    hgnc = list(csv.DictReader(f, delimiter="\t"))
approved = {r["symbol"]: r for r in hgnc}
for r in hgnc:
    if r["locus_group"] != "pseudogene" or r["hgnc_id"] not in gene_span:
        continue
    for sym in sorted({x.strip('"') for x in (r["alias_symbol"] + "|" + r["prev_symbol"]).split("|")} - {""}):
        d = approved.get(sym)
        if d and d["hgnc_id"] != r["hgnc_id"] and d["gencc"]:
            chrom, s, e = gene_span[r["hgnc_id"]]
            rows.append((chrom, s, e, "gene_alias", "remap_gene", f"{sym}|{r['symbol']}|{r['hgnc_id']}|{d['hgnc_id']}",
                         SRC_ALIAS))

rows.sort(key=lambda x: (ORDER[x[0]], x[1], x[2], x[3], x[5]))
counts = defaultdict(int)
for x in rows:
    counts[x[3]] += 1
with open(out, "w") as f:
    f.write("# assets/sv_rule_regions.tsv: the GRCh38 regions where bin/annotsv_postfilter.py overrides AnnotSV's class 4\n"
            "# or 5. Written by scripts/ci/build-sv-rule-regions.sh, which says how each rule's rows are chosen;\n"
            "# rebuild it there, do not edit it by hand. start is 0-based and end 1-based, as in a BED file.\n"
            "# Rows: " + ", ".join(f"{k} {counts[k]}" for k in sorted(counts)) + ".\n")
    f.write("#chrom\tstart\tend\trule\taction\tname\tsource\n")
    for x in rows:
        f.write("\t".join(map(str, x)) + "\n")
print("Rows: " + ", ".join(f"{k} {counts[k]}" for k in sorted(counts)), file=sys.stderr)
PY
mv "${OUT}.part" "$OUT"
echo "Wrote ${OUT}" >&2
