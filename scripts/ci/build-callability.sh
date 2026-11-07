#!/usr/bin/env bash
# build-callability.sh: rebuild assets/callability_genes.tsv, the genes where
# short reads often cannot be placed, so a "no finding" there needs a reading
# of its own (step 16b measures each one, the report lists the ones it could
# not read).
#
# Usage: scripts/ci/build-callability.sh [OUT]   (default: assets/callability_genes.tsv)
#        scripts/ci/build-callability.sh --list    the downloads, one "NAME URL SHA256" per line
#
# Downloads seven pinned files (about 85 MB, each checked against its
# sha256) into a temporary folder, or reads them from CALLABILITY_CACHE when
# that folder holds them with the right sha256, writes OUT and deletes the
# downloads:
#
#   GIAB CMRG v1.00's GRCh38 gene coordinates (Wagner 2022, PMID 35132260):
#     the 273 challenging medically relevant genes
#   GIAB genome stratifications v3.6, GRCh38: segmental duplications, the
#     false duplications of GRCh38 (both copies; Behera 2023, PMID 36810122)
#     and the collapsed duplications
#   GENCODE 50's basic annotation with the patches and ALT scaffolds
#   the GRCh38 assembly report of NCBI (GCA_000001405.15), which names the
#     ALT scaffolds
#
# A gene is a gene of GENCODE 50 with an HGNC id, on chr1-22, X or Y unless it
# is ALT-only, and protein-coding unless it is a CMRG gene (a few of those are
# T-cell receptor segments). It is listed with every reason that applies:
#
#   CMRG               a CMRG gene (matched by symbol, or by the protein-coding
#                      gene that covers most of its CMRG span when GENCODE
#                      has renamed it)
#   segdup             half or more of the gene body in a segmental duplication
#   false_duplication  the gene body overlaps a false duplication of GRCh38
#                      (either copy): reads of the real gene get MAPQ 0
#   collapsed          the gene body overlaps a collapsed duplication: two
#                      copies in people, one in GRCh38
#   ALT-only           the gene is on an ALT scaffold and nowhere on chr1-22,
#                      X, Y or M, so the default no-ALT reference cannot hold it
#
# and the step that calls it with a method of its own (owner): parascopy for
# the genes of assets/paralog_loci.tsv, pypgx for the bam genes of
# assets/pgx_genes.tsv that step 32 can call on the no-ALT reference,
# cn_loci:<locus> for a gene a gene row of bin/cn_loci.tsv names (its locus,
# the locus up to its first '_', or a target name), in that order; '.' for
# none. A change to one of those three files changes the owners, so rebuild
# after it.
#
# assets/paralog_loci.tsv names owners only; it is not a source of rows. The
# rows say where short reads fail by the measures above, gene body by gene
# body, and a Parascopy model may cover only part of a gene. Seventeen catalog
# genes are not rows: seven are pseudogenes (CELP, GBAP1, NCF1B, NCF1C,
# OTOAP1, PMS2CL, STRCP1), left out as not protein-coding, and no reason
# above applies to the other ten (ABCC6, CFH, CFHR4, CHRNA7, DEFA3, KRT86,
# NEB, PIK3CA, PLG, VWF). Step 35's own report section reads all of them.
# tests/test_paralog_loci.py holds that list and goes red when it changes.
#
# Needs curl, python3 and sha256sum or shasum. The e2e case
# tests/e2e/callability-1-catalog.sh runs it and fails on any difference with
# the committed file.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
OUT=${1:-${REPO}/assets/callability_genes.tsv}

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
STRAT=${GIAB}/genome-stratifications/v3.6/GRCh38@all
FILES=(
  "cmrg.bed|${GIAB}/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SupplementaryFiles/GRCh38_CMRG_benchmark_gene_coordinates.bed|67ce028e3acc0f4946104b8a6e0b76ce1ba9da8eb8457e3d7185bd1e81cbf203"
  "segdups.bed.gz|${STRAT}/SegmentalDuplications/GRCh38_segdups.bed.gz|09ac029a7d0d681a9f0ed6581a6bc10fc708a834bc3f91c46a9bef504f6055c5"
  "falsedup_correct.bed.gz|${STRAT}/OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz|44d07f373f6841a99fb85bd712e535e3e1364a033fd37cfd42e9e17988277717"
  "falsedup_incorrect.bed.gz|${STRAT}/OtherDifficult/GRCh38_false_duplications_incorrect_copy.bed.gz|d2de1abe40a76cb1d29dfa5ae4c0a9ae05f2301641721f4a78f6b72c3a0b5bcd"
  "collapsed.bed.gz|${STRAT}/OtherDifficult/GRCh38_collapsed_duplication_FP_regions.bed.gz|32b6f04e25afeb12078112a771f76e3f880f3a1d1bb2f454f60c5462644bf8af"
  "gencode.gtf.gz|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.chr_patch_hapl_scaff.basic.annotation.gtf.gz|e9199c3ef075e54f339411a090385dd65ffe6ff59725403230cf1688ce9e8524"
  "report.txt|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt|aa733ce92719f6339b3a4540137d4322894434363314c4b7704ba6691e9eab66"
)

if [ "${1:-}" = --list ]; then
  for f in "${FILES[@]}"; do
    IFS='|' read -r name url sum <<<"$f"
    echo "${name} ${url} ${sum}"
  done
  exit 0
fi

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# get NAME URL SHA256: WORK/NAME from CALLABILITY_CACHE when it matches, else
# downloaded; stop unless it matches.
get() {
  local cached="${CALLABILITY_CACHE:-}/$1" got
  if [ -n "${CALLABILITY_CACHE:-}" ] && [ -f "$cached" ] && [ "$(sha256 "$cached")" = "$3" ]; then
    cp "$cached" "${WORK}/$1"
    return
  fi
  echo "Downloading $2" >&2
  # Tries spread over about five minutes: EBI's FTP server refused every
  # connection for minutes at a time (curl exit 7, 2026-11-06).
  curl -fsSL --retry 6 --retry-all-errors --retry-delay 30 --retry-max-time 300 -o "${WORK}/$1" "$2"
  got=$(sha256 "${WORK}/$1")
  if [ "$got" != "$3" ]; then
    echo "ERROR: $1 has sha256 ${got}, want $3" >&2
    exit 1
  fi
  if [ -n "${CALLABILITY_CACHE:-}" ]; then
    mkdir -p "$CALLABILITY_CACHE" && cp "${WORK}/$1" "$cached"
  fi
}
for f in "${FILES[@]}"; do
  IFS='|' read -r name url sum <<<"$f"
  get "$name" "$url" "$sum"
done

python3 - "$WORK" "$REPO" "${OUT}.part" <<'PY'
import csv
import gzip
import re
import sys
from collections import defaultdict

work, repo, out = sys.argv[1:4]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ORDER = {c: i for i, c in enumerate(CHROMS)}
PRIMARY = set(CHROMS) | {"chrM"}
ATT = re.compile(r'(\S+) "([^"]*)"')
SEGDUP_SHARE = 0.5
REASONS = ("CMRG", "segdup", "false_duplication", "collapsed", "ALT-only")


def bed(name):
    d = defaultdict(list)
    with (gzip.open if name.endswith(".gz") else open)(f"{work}/{name}", "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and not line.startswith(("#", "track", "browser")):
                d[c[0]].append((int(c[1]), int(c[2])))
    for k in d:
        d[k].sort()
    return d


def covered(chrom, s, e, d):
    """Bases of [s, e) inside the intervals of d (merged as it goes)."""
    n, cur = 0, s
    for a, b in d.get(chrom, []):
        if b <= cur:
            continue
        if a >= e:
            break
        a, b = max(a, cur), min(b, e)
        if b > a:
            n += b - a
            cur = b
    return n


# ALT scaffolds: GenBank accession -> UCSC name (GENCODE names them by accession)
alt = {}
with open(f"{work}/report.txt") as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        if not line.startswith("#") and len(c) >= 10 and c[1] == "alt-scaffold":
            alt[c[4]] = c[9]

genes = {}                       # hgnc id -> [symbol, chrom, start0, end] (any gene type)
coding = set()                   # hgnc ids of the protein-coding ones
on_primary = set()               # every gene symbol on chr1-22, X, Y or M
alt_genes = {}                   # hgnc id -> [symbol, UCSC ALT name] (first ALT copy)
with gzip.open(f"{work}/gencode.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.split("\t", 9)
        if c[2] != "gene":
            continue
        a = dict(ATT.findall(c[8]))
        if c[0] in PRIMARY:
            on_primary.add(a.get("gene_name"))
        h = a.get("hgnc_id")
        if not h or a["gene_id"].endswith("_PAR_Y"):
            continue
        if c[0] in ORDER and h not in genes:
            genes[h] = [a["gene_name"], c[0], int(c[3]) - 1, int(c[4])]
            if a.get("gene_type") == "protein_coding":
                coding.add(h)
        elif c[0] in alt and h not in alt_genes and a.get("gene_type") == "protein_coding":
            alt_genes[h] = [a["gene_name"], alt[c[0]]]

reasons = defaultdict(set)
by_symbol = {g[0]: h for h, g in genes.items()}
with open(f"{work}/cmrg.bed") as f:
    for line in f:
        chrom, s, e, sym = line.rstrip("\n").split("\t")[:4]
        s, e = int(s), int(e)
        h = by_symbol.get(sym)
        if h is None:
            # renamed since CMRG v1.00: the protein-coding gene that covers most of its span
            best = max(((min(e, genes[k][3]) - max(s, genes[k][2]), k) for k in coding if genes[k][1] == chrom),
                       default=(0, None))
            if best[0] <= 0:
                sys.exit(f"ERROR: CMRG gene {sym} ({chrom}:{s}-{e}) matches no GENCODE 50 gene")
            h = best[1]
        reasons[h].add("CMRG")

segdups, collapsed = bed("segdups.bed.gz"), bed("collapsed.bed.gz")
falsedup = bed("falsedup_correct.bed.gz")
for k, v in bed("falsedup_incorrect.bed.gz").items():
    falsedup[k] = sorted(falsedup[k] + v)
for h in coding:
    sym, chrom, s, e = genes[h]
    if covered(chrom, s, e, segdups) >= SEGDUP_SHARE * (e - s):
        reasons[h].add("segdup")
    if covered(chrom, s, e, falsedup):
        reasons[h].add("false_duplication")
    if covered(chrom, s, e, collapsed):
        reasons[h].add("collapsed")
for h, (sym, name) in alt_genes.items():
    if h not in genes and sym not in on_primary:
        reasons[h].add("ALT-only")


def table(path):
    with open(path) as f:
        return list(csv.DictReader((x for x in f if x.strip() and not x.startswith("#")), delimiter="\t",
                                   quoting=csv.QUOTE_NONE))


parascopy = {g for r in table(f"{repo}/assets/paralog_loci.tsv") for g in r["genes"].split(",")}
# A gene step 32 skips on the default no-ALT reference (callable_without_alt no) has no
# caller there, and step 32 refuses a BAM with ALT contigs, so it is not pypgx's.
pypgx = {r["gene"] for r in table(f"{repo}/assets/pgx_genes.tsv")
         if r["mode"] == "bam" and r["callable_without_alt"] != "no"}
# A cn_loci gene row reports the genes it names (its locus, as RHD, or LPA for LPA_KIV2,
# or a target, as HBA1 and HBA2), not every gene whose span holds one of its targets.
cn_names = {}
for r in table(f"{repo}/bin/cn_loci.tsv"):
    if r["kind"] == "gene":
        names = {r["locus"], r["locus"].split("_")[0]}
        for part in r["targets"].split(";"):
            names |= {part.split(":")[0], part.split(":")[0].split("_")[0]}
        for n in names:
            cn_names.setdefault(n, r["locus"])


def owner(sym):
    if sym in parascopy:
        return "parascopy"
    if sym in pypgx:
        return "pypgx"
    if sym in cn_names:
        return f"cn_loci:{cn_names[sym]}"
    return "."


rows = []
for h, why in reasons.items():
    if h in genes:
        sym, chrom, s, e = genes[h]
        rows.append((ORDER[chrom], s, sym, h, chrom, s, e, why, owner(sym)))
    else:
        sym, name = alt_genes[h]
        rows.append((len(ORDER), 0, sym, h, name, ".", ".", why, owner(sym)))
rows.sort()
syms = [r[2] for r in rows]
dup = sorted({x for x in syms if syms.count(x) > 1})
if dup:
    sys.exit(f"ERROR: gene symbols listed twice: {', '.join(dup)}")
counts = {k: sum(1 for r in rows if k in r[7]) for k in REASONS}
owned = sum(1 for r in rows if r[8] != ".")
with open(out, "w") as f:
    f.write("# assets/callability_genes.tsv: the genes where short reads are often not placed with confidence, so a\n"
            "# 'no finding' there needs a reading of its own. Step 16b measures each gene on chr1-22, X and Y; the\n"
            "# report lists the ones it could not read. Written by scripts/ci/build-callability.sh, which says how the\n"
            "# rows and the owners are chosen; rebuild it there, do not edit it by hand. start is 0-based and end\n"
            "# 1-based, as in a BED file; an ALT-only gene has its ALT scaffold as chrom and '.' as start and end.\n"
            f"# Genes: {len(rows)} (" + ", ".join(f"{k} {counts[k]}" for k in REASONS) + f"); with an owner: {owned}.\n")
    f.write("gene\thgnc_id\tchrom\tstart\tend\treasons\towner\n")
    for r in rows:
        f.write("\t".join(map(str, (r[2], r[3], r[4], r[5], r[6], ",".join(k for k in REASONS if k in r[7]), r[8])))
                + "\n")
print(f"Genes: {len(rows)} (" + ", ".join(f"{k} {counts[k]}" for k in REASONS) + f"); with an owner: {owned}",
      file=sys.stderr)
PY
mv "${OUT}.part" "$OUT"
echo "Wrote ${OUT}" >&2
