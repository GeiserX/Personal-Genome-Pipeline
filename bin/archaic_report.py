#!/usr/bin/env python3
"""archaic_report.py: Neanderthal and Denisovan segments in one genome, with IBDmix.

IBDmix (Chen et al. 2020, Cell 180:677, doi:10.1016/j.cell.2020.01.012) compares
a modern genome with an archaic one site by site and adds up a log-odds score
(LOD) that a stretch of the modern genome shares a recent ancestor with the
archaic genome, against the chance of the same genotypes in the modern
population. It needs no unadmixed reference population, but it does need the
allele frequency of each site in the modern population, and it takes that from
the modern people in its input. With one genome as the only modern input the
frequency is that person's own genotype and no site scores (find_frequency in
IBDmix's src/Genotype_Reader.cc drops every site below its minor allele count),
so `genotypes` adds PSEUDO stand-ins whose genotypes carry the 1000 Genomes
allele count of the chosen population at each site. IBDmix scores each modern
column from the site's frequency and that column's genotype only, so your rows
are the ones it would give you inside a population sample of that frequency;
the stand-ins' own rows are dropped.

Subcommands:

  table      (setup.sh --archaic) one GRCh37 chromosome: the 1000 Genomes
             phase 3 biallelic SNVs with their allele frequencies, and the
             genotypes of the Altai Neanderthal, the Vindija 33.19 Neanderthal
             and the Altai Denisovan at those sites and wherever one of the
             three carries a non-reference allele, as a sites-only VCF that
             Picard lifts to GRCh38. Each archaic genome is read inside its
             own mask (--mask NAME=BED, the MPI-EVA FilterBed of that
             genome): a call outside it is kept as masked, code 10 plus its
             ALT count. Every call at a site inside an excluded region
             (--exclude BED: segmental duplications, the 5 bp around a modern
             indel) is masked the same way. That is how Chen et al. filter:
             the site stays in IBDmix's input and IBDmix's --mask stops it
             scoring, except that a homozygous mismatch there still counts
             against a segment.
  strict     (setup.sh --archaic) the lifted GRCh38 rows of one chromosome with
             every row outside the 1000 Genomes strict mask marked masked in
             all three archaic columns (the same code 10 to 12), instead of
             dropped, so IBDmix still scores a mismatch there.
  genotypes  (step 43) IBDmix's genotype file of each chromosome: the three
             archaic columns, SAMPLE (yours, from the VCF) and the stand-ins;
             the mask of each archaic genome (chrN.<genome>.mask.bed, the
             positions where its call is masked); and spans.tsv, the span and
             site count of each chromosome.
  report     (step 43) IBDmix's segments for SAMPLE against each archaic genome,
             kept at LOD 4 and 50 kb or longer (the cutoffs of Chen et al.),
             merged into Neanderthal and Denisovan segments, with the share of
             the analysed span, the ten longest segments and the protein-coding
             genes in them.

Codes: an archaic or modern genotype is the number of ALT alleles (0, 1, 2), or
9 when it is unknown, the coding IBDmix's genotype file uses. In the table an
archaic call outside that genome's mask is 10, 11 or 12 (MASKED plus the ALT
count); `genotypes` writes it to IBDmix as 0, 1 or 2 and lists the position in
that genome's mask file (chrN.<genome>.mask.bed), which step 43 passes to
IBDmix with -r.

Used by scripts/setup.sh, scripts/43-archaic.sh and the ARCHAIC_* processes
(bin/ is on the task PATH); tests/test_archaic_report.py checks it. Standard
library only.
"""
import argparse
import bisect
import gzip
import heapq
import os
import re
import sys

ARCHAICS = ("Altai", "Vindija", "Denisova")
NEANDERTHALS = ("Altai", "Vindija")
POPS = ("ALL", "AFR", "AMR", "EAS", "EUR", "SAS")
TABLE_HEADER = ["chrom", "pos", "ref", "alt", "Altai", "Vindija", "Denisova"] + [f"AF_{p}" for p in POPS]
MASKED = 10           # table code: MASKED + ALT count, a call outside the genome's mask
PSEUDO = 100          # stand-in modern people: 200 alleles, about one 1000 Genomes population
MIN_LOD = 4.0         # Chen et al. 2020: LOD 4
MIN_LENGTH = 50000    # and 50 kb
TOP = 10
SEGMENT_HEADER = ["chrom", "start", "end", "length_bp", "ancestry", "matched", "lod", "genes"]
CAVEAT = ("IBDmix finds stretches where your genome matches a Neanderthal or the Denisovan genome more "
          "closely than the allele frequencies of a modern population explain. A stretch counts when either "
          "of your two copies carries it, and the share is of the analysed span of the autosomes. Some of what "
          "it finds in anyone is old ancestry shared with the archaic genomes rather than gene flow from them: "
          "Chen et al. 2020 find Neanderthal sequence in every population they studied, Africans included. "
          "The numbers move with the cutoffs (LOD 4, 50 kb) and with the frequencies used. An estimate, not "
          "a measurement.")


def die(msg, code=1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gt_code(gt):
    """'0/1', '1|1', '0/0' -> 1, 2, 0; a missing, partial or multi-allelic call -> None."""
    a = re.split(r"[/|]", gt)
    if len(a) != 2 or any(x not in ("0", "1") for x in a):
        return None
    return int(a[0]) + int(a[1])


# --- table (setup) ------------------------------------------------------------------

def read_modern(path):
    """1000 Genomes lines: POS REF ALT AF AFR_AF AMR_AF EAS_AF EUR_AF SAS_AF (bcftools query)."""
    with opener(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                die(f"{path}: a line has {len(c)} columns, want 9: {line[:80]!r}")
            if len(c[1]) != 1 or len(c[2]) != 1:
                continue
            yield int(c[0]), 0, (c[1], c[2], c[3:9])


def read_archaic(path, idx):
    """Archaic lines: POS REF ALT GT (bcftools query)."""
    with opener(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                die(f"{path}: a line has {len(c)} columns, want 4: {line[:80]!r}")
            yield int(c[0]), 1 + idx, (c[1], c[2], c[3])


def read_mask(path, chrom):
    """The intervals of one chromosome of a mask BED, sorted and merged.
    A BED is 0-based and half-open, so [start, end) holds the 1-based
    positions start + 1 to end. The chromosome is matched with and without
    the chr prefix, because the archaic masks name it 22 and others chr22."""
    num = str(chrom)[3:] if str(chrom).startswith("chr") else str(chrom)
    want = {num, f"chr{num}"}
    iv = []
    with opener(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            c = line.split()
            if len(c) < 3 or c[0] not in want:
                continue
            s, e = int(c[1]), int(c[2])
            if e <= s:
                continue          # an empty or reversed interval holds no position
            iv.append((s, e))
    iv.sort()
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


class Mask:
    """Whether a 1-based position is inside a mask. Positions must come in
    ascending order, as they do in the table build: the walk over the
    intervals never goes back."""

    def __init__(self, intervals):
        self.iv = intervals
        self.i = 0

    def holds(self, pos):
        p = pos - 1
        while self.i < len(self.iv) and self.iv[self.i][1] <= p:
            self.i += 1
        return self.i < len(self.iv) and self.iv[self.i][0] <= p


def archaic_code(rec, ref, alt):
    """The ALT count of one archaic record at a row with this REF and ALT, or 9."""
    if rec is None:
        return 9
    aref, aalt, gt = rec
    if aref != ref:
        return 9
    code = gt_code(gt)
    if code is None:
        return 9
    if aalt in (".", "<*>") or code == 0:
        return 0 if code == 0 else 9
    return code if aalt == alt else 9


def table_row(pos, modern, arch):
    """The row of one position, or None. modern: (REF, ALT, AFs) or None;
    arch: one (REF, ALT, GT) or None per archaic genome."""
    if modern is not None:
        ref, alt, afs = modern
    else:
        # A site no modern person varies at: the ALT of the first archaic
        # genome with a non-reference call, every population frequency 0.
        cand = [r for r in arch if r is not None and gt_code(r[2]) not in (None, 0)
                and len(r[0]) == 1 and len(r[1]) == 1]
        if not cand:
            return None
        ref, alt, afs = cand[0][0], cand[0][1], ["0"] * len(POPS)
    codes = [archaic_code(r, ref, alt) for r in arch]
    if all(c == 9 for c in codes):
        return None
    if modern is None and not any(c in (1, 2) for c in codes):
        return None
    clean = []
    for v in afs:
        try:
            x = float(v)
        except ValueError:
            x = 0.0
        clean.append(f"{min(max(x, 0.0), 1.0):.6g}")
    return ref, alt, codes, clean


def cmd_table(a):
    names = []
    streams = [read_modern(a.modern)]
    for i, spec in enumerate(a.archaic):
        name, _, path = spec.partition("=")
        if name not in ARCHAICS or not path:
            die(f"--archaic takes NAME=FILE with NAME one of {', '.join(ARCHAICS)}, got {spec!r}", 2)
        names.append(name)
        streams.append(read_archaic(path, i))
    if names != list(ARCHAICS):
        die(f"--archaic must name {', '.join(ARCHAICS)} in that order, got {', '.join(names)}", 2)
    masks = {}
    for spec in a.mask:
        name, _, path = spec.partition("=")
        if name not in ARCHAICS or not path:
            die(f"--mask takes NAME=FILE with NAME one of {', '.join(ARCHAICS)}, got {spec!r}", 2)
        if name in masks:
            die(f"--mask names {name} twice", 2)
        iv = read_mask(path, a.chrom)
        if not iv:
            die(f"{path} holds no interval on chromosome {a.chrom}: a mask that matches nothing "
                "would make every site unknown")
        masks[name] = Mask(iv)
    excluded = []
    for path in a.exclude:
        iv = read_mask(path, a.chrom)
        if not iv:
            die(f"{path} holds no interval on chromosome {a.chrom}: an empty exclusion list is a "
                "download that went wrong, not a chromosome without segmental duplications")
        excluded.append(Mask(iv))
    chrom = a.chrom if a.chrom.startswith("chr") else f"chr{a.chrom}"
    out = open(a.out, "w")
    out.write("##fileformat=VCFv4.2\n")
    out.write(f"##contig=<ID={chrom}>\n")
    for n in ARCHAICS:
        out.write(f'##INFO=<ID={n},Number=1,Type=Integer,Description="ALT alleles in the {n} genome (9: unknown; 10 to 12: 0 to 2 outside its mask)">\n')
    for p in POPS:
        out.write(f'##INFO=<ID=AF_{p},Number=1,Type=Float,Description="1000 Genomes phase 3 ALT frequency ({p})">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    rows = 0
    cur, modern, arch, last = None, None, [None] * len(ARCHAICS), -1

    def flush():
        nonlocal rows
        if cur is None:
            return
        # Every Mask is walked for every position, before any early return,
        # because the walk only goes forward.
        excl = any([m.holds(cur) for m in excluded])
        outside = [excl or (n in masks and not masks[n].holds(cur)) for n in ARCHAICS]
        r = table_row(cur, modern, arch)
        if r is None:
            return
        ref, alt, codes, afs = r
        # A genome's call outside its mask, or at an excluded site, is kept but
        # marked masked (MASKED + code). Step 43 hands IBDmix the call and a
        # mask of those positions, as Chen et al. do: IBDmix scores nothing
        # there except a homozygous mismatch, which still counts against a
        # segment. Dropping the call instead (code 9) let segments run through
        # masked sequence nothing could contradict.
        codes = [c + MASKED if (out and c != 9) else c for c, out in zip(codes, outside)]
        info = ";".join([f"{n}={c}" for n, c in zip(ARCHAICS, codes)] + [f"AF_{p}={v}" for p, v in zip(POPS, afs)])
        out.write(f"{chrom}\t{cur}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")
        rows += 1

    for pos, src, rec in heapq.merge(*streams, key=lambda t: (t[0], t[1])):
        if pos < last:
            die(f"input not sorted by position at {pos}")
        if pos != cur:
            flush()
            cur, modern, arch = pos, None, [None] * len(ARCHAICS)
        last = pos
        if src == 0:
            # A second modern record at one position (a split multi-allelic
            # site) makes the row ambiguous: keep neither.
            modern = rec if modern is None else ("N", "N", ["0"] * len(POPS))
        else:
            arch[src - 1] = rec
    flush()
    out.close()
    print(f"{chrom}: {rows} rows", file=sys.stderr)


def cmd_strict(a):
    """The lifted GRCh38 rows of one chromosome (no header, the 13 table
    columns), once per position, with the calls of every row outside the 1000
    Genomes strict mask marked masked (MASKED plus the ALT count, 9 stays 9).
    The strict mask is the same for the three genomes, so all three are
    marked; step 43 then hands IBDmix those positions as its mask."""
    iv = read_mask(a.strict, a.chrom)
    if not iv:
        die(f"{a.strict} holds no interval on chromosome {a.chrom}: a mask that matches nothing "
            "would mark every site masked")
    starts = [s for s, _ in iv]
    seen = set()
    rows = marked = 0
    with open(a.rows) as f, open(a.out, "w") as out:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) != len(TABLE_HEADER):
                die(f"{a.rows}: a row has {len(c)} columns, want {len(TABLE_HEADER)}: {line[:80]!r}")
            pos = int(c[1])
            if pos in seen:
                continue
            seen.add(pos)
            i = bisect.bisect_right(starts, pos - 1) - 1
            if not (i >= 0 and pos - 1 < iv[i][1]):
                hit = False
                for k in (4, 5, 6):
                    if c[k] in ("0", "1", "2"):
                        c[k] = str(int(c[k]) + MASKED)
                        hit = True
                marked += hit
            out.write("\t".join(c) + "\n")
            rows += 1
    print(f"{a.chrom}: {rows} rows, {marked} outside the strict mask and marked masked", file=sys.stderr)


# setup.sh --archaic stamps its table with the sha256 of everything above this
# section, so an edit below it does not rebuild the table (about 165 GB of downloads).

# --- genotypes (step 43) ----------------------------------------------------------

def pseudo_columns(n):
    """The stand-in columns for each ALT allele count 0..2n: as many 2s and
    one 1 as the count needs, then 0s."""
    cols = []
    for ac in range(2 * n + 1):
        g = ["2"] * (ac // 2) + (["1"] if ac % 2 else [])
        g += ["0"] * (n - len(g))
        cols.append("\t".join(g))
    return cols


def read_vcf(path):
    """(chrom, pos, ref, alt, code) of the biallelic SNVs of a one-sample VCF.
    A record whose FILTER is not PASS or '.' gives code 9 (not read), so a
    DeepVariant RefCall is neither your allele nor the reference."""
    with opener(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                if len(c) != 10:
                    die(f"{path} holds {max(len(c) - 9, 0)} samples; step 43 reads a one-sample VCF "
                        "(keep one with bcftools view -s <sample>).")
                continue
            if len(c) < 10 or len(c[3]) != 1 or c[4] == ".":
                continue
            if "," in c[4] or len(c[4]) != 1:
                yield c[0], int(c[1]), c[3], None, 9
                continue
            fmt = c[8].split(":")
            if "GT" not in fmt:
                continue
            vals = c[9].split(":")
            gt = vals[fmt.index("GT")] if fmt.index("GT") < len(vals) else "."
            code = gt_code(gt)
            if code is None or c[6] not in ("PASS", "."):
                code = 9
            yield c[0], int(c[1]), c[3], c[4], code


CLASH = object()
ARCHAIC_CODES = {str(v) for v in (0, 1, 2, 9, MASKED, MASKED + 1, MASKED + 2)}


def pick_call(prev, call):
    """The call to keep when a position comes twice, as in a VCF split with
    bcftools norm -m-: each is (REF, ALT, code), ALT None for an insertion or a
    multi-allelic record. An SNV wins over such a record, in either order. Of
    one SNV twice, a usable record wins over code 9, and two different usable
    codes are CLASH. Of two different SNVs, the
    one that carries its ALT wins over a 0/0 or unusable one (a split 0/1 at
    A>G,T is A>G 0/1 and A>T 0/0); both 0/0 is the reference; both carrying
    their ALT, or neither readable, is CLASH (unknown, 9)."""
    if prev is None:
        return call
    if prev[1] is None or call[1] is None:
        return call if prev[1] is None and call[1] is not None else prev
    if call[1] == prev[1]:
        if call[2] == 9 or call[2] == prev[2]:
            return prev
        return call if prev[2] == 9 else CLASH
    has_prev, has_call = prev[2] in (1, 2), call[2] in (1, 2)
    if has_prev != has_call:
        return prev if has_prev else call
    if not has_prev and prev[2] == 0 and call[2] == 0:
        return (prev[0], None, 0)
    return CLASH


def read_table(path):
    """Rows of one chromosome of the GRCh38 table, with the header checked."""
    with opener(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head != TABLE_HEADER:
            die(f"{path} does not start with the table header ({' '.join(TABLE_HEADER)})")
        for line in f:
            c = line.rstrip("\n").split("\t")
            yield c


def cmd_genotypes(a):
    if a.pop not in POPS:
        die(f"--pop must be one of {', '.join(POPS)}, got {a.pop!r}", 2)
    af_col = TABLE_HEADER.index(f"AF_{a.pop}")
    tables = {}
    for name in sorted(os.listdir(a.table_dir)):
        m = re.fullmatch(r"(chr(?:[1-9]|1[0-9]|2[0-2]))\.tsv\.gz", name)
        if m:
            tables[m.group(1)] = os.path.join(a.table_dir, name)
    if not tables:
        die(f"{a.table_dir} holds no chrN.tsv.gz: install the archaic data with scripts/setup.sh --archaic")
    os.makedirs(a.out_dir, exist_ok=True)
    pseudo = pseudo_columns(a.pseudo)
    two_n = 2 * a.pseudo
    head = "\t".join(["chrom", "pos", "ref", "alt", *ARCHAICS, "SAMPLE"] + [f"P{i:03d}" for i in range(1, a.pseudo + 1)])
    with open(os.path.join(a.out_dir, "samples.txt"), "w") as f:
        f.write("SAMPLE\n" + "".join(f"P{i:03d}\n" for i in range(1, a.pseudo + 1)))
    spans = []
    vcf = read_vcf(a.vcf)
    pending = next(vcf, None)
    done = set()
    while pending is not None:
        chrom = pending[0]
        if chrom in done:
            die(f"{a.vcf} is not sorted: {chrom} comes back after another chromosome")
        done.add(chrom)
        if chrom not in tables:
            while pending is not None and pending[0] == chrom:
                pending = next(vcf, None)
            continue
        # The sample's calls on this chromosome, merged with the table by position.
        calls, clash = {}, set()
        while pending is not None and pending[0] == chrom:
            pos = pending[1]
            if pos not in clash:
                calls[pos] = pick_call(calls.get(pos), pending[2:])
                if calls[pos] is CLASH:
                    calls[pos] = (pending[2], None, 9)
                    clash.add(pos)
            pending = next(vcf, None)
        n = first = last = dropped = carried = 0
        num = chrom[3:]
        path = os.path.join(a.out_dir, f"{chrom}.gt.gz")
        # Per archaic genome, the runs of written rows whose call is masked:
        # [first, last] positions, closed at the next unmasked row.
        runs = {g: [] for g in ARCHAICS}
        open_run = {g: None for g in ARCHAICS}
        with gzip.open(path, "wt", compresslevel=1) as out:
            out.write(head + "\n")
            for c in read_table(tables[chrom]):
                pos, ref, alt = int(c[1]), c[2], c[3]
                arch = []
                for g, v in zip(ARCHAICS, c[4:7]):
                    if v not in ARCHAIC_CODES:
                        die(f"{tables[chrom]}: {chrom}:{pos} has {g} code {v!r}; the table holds 0, 1, 2, 9 "
                            f"or {MASKED} to {MASKED + 2}")
                    arch.append(int(v))
                call = calls.get(pos)
                if call is None:
                    code = 0
                elif call[0] != ref or (call[1] is not None and call[1] != alt):
                    # Another allele than the row's: IBDmix's generate_gt drops the line.
                    dropped += 1
                    continue
                else:
                    code = call[2]
                if code in (1, 2):
                    carried += 1
                try:
                    ac = int(round(float(c[af_col]) * two_n))
                except ValueError:
                    ac = 0
                ac = min(max(ac, 0), two_n)
                for g, v in zip(ARCHAICS, arch):
                    if v >= MASKED:
                        if open_run[g] is None:
                            open_run[g] = [pos, pos]
                        else:
                            open_run[g][1] = pos
                    elif open_run[g] is not None:
                        runs[g].append(open_run[g])
                        open_run[g] = None
                gts = "\t".join(str(v - MASKED if v >= MASKED else v) for v in arch)
                out.write(f"{num}\t{pos}\t{ref}\t{alt}\t{gts}\t{code}\t{pseudo[ac]}\n")
                if n == 0:
                    first = pos
                last = pos
                n += 1
        if n == 0:
            os.remove(path)
            continue
        # IBDmix's mask (-r): chromosome as in the genotype file, and a line
        # covers the positions start + 1 to end, as in a BED.
        for g in ARCHAICS:
            if open_run[g] is not None:
                runs[g].append(open_run[g])
            with open(os.path.join(a.out_dir, f"{chrom}.{g}.mask.bed"), "w") as f:
                f.writelines(f"{num}\t{s0 - 1}\t{e0}\n" for s0, e0 in runs[g])
        spans.append((chrom, first, last, n, carried, dropped))
    with open(os.path.join(a.out_dir, "spans.tsv"), "w") as f:
        f.write("chrom\tfirst\tlast\tsites\tsites_with_your_alt\tsites_other_allele\n")
        for s in spans:
            f.write("\t".join(map(str, s)) + "\n")
    if not spans:
        die(f"{a.vcf} has no records on a chromosome of the archaic table ({', '.join(sorted(tables))})")
    print(f"{len(spans)} chromosomes, {sum(s[3] for s in spans)} sites", file=sys.stderr)


# --- report (step 43) ------------------------------------------------------------

def read_segments(path, sample="SAMPLE"):
    """IBDmix's output rows for SAMPLE: (chrom, start, end, lod)."""
    segs = []
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head[:5] != ["ID", "chrom", "start", "end", "slod"]:
            die(f"{path} is not IBDmix output (header {' '.join(head[:5])})")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or c[0] != sample:
                continue
            segs.append((c[1], int(c[2]), int(c[3]), float(c[4])))
    return segs


def keep(segs, min_lod=MIN_LOD, min_length=MIN_LENGTH):
    return [s for s in segs if s[3] >= min_lod and s[2] - s[1] >= min_length]


def merge(segs):
    """Overlapping or touching segments of one class as one: (chrom, start, end, lod, members)."""
    out = []
    for chrom, start, end, lod, who in sorted(segs, key=lambda s: (chrom_key(s[0]), s[1], s[2])):
        if out and out[-1][0] == chrom and start <= out[-1][2]:
            c, s0, e0, l0, w0 = out[-1]
            out[-1] = (c, s0, max(e0, end), max(l0, lod), w0 | who)
        else:
            out.append((chrom, start, end, lod, set(who)))
    return out


def chrom_key(c):
    c = str(c)
    c = c[3:] if c.startswith("chr") else c
    return int(c) if c.isdigit() else 99


def overlaps(seg, others):
    return [o for o in others if o[0] == seg[0] and o[1] < seg[2] and seg[1] < o[2]]


def read_genes(path):
    """Protein-coding genes of a GENCODE gene GTF: chrom -> [(start, end, name)]."""
    genes = {}
    if not path:
        return genes
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            if 'gene_type "protein_coding"' not in c[8]:
                continue
            m = re.search(r'gene_name "([^"]+)"', c[8])
            if m:
                genes.setdefault(c[0], []).append((int(c[3]), int(c[4]), m.group(1)))
    return genes


def genes_in(genes, chrom, start, end):
    """The names of the genes overlapping [start, end), once each, in order.
    IBDmix's start is the segment's first site and its end the next site, past
    the segment (its README, --inclusive-end); a GTF gene is [s, e], 1-based."""
    name = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    seen = []
    for s, e, n in sorted(genes.get(name, [])):
        if s < end and e >= start and n not in seen:
            seen.append(n)
    return seen


def classify(per_archaic):
    """Neanderthal: the union of the Altai and Vindija segments. Denisovan: the
    Denisova segments that overlap no Neanderthal one (a stretch the Denisova
    genome and a Neanderthal both match is counted as Neanderthal, the more
    common source outside Oceania). Returns a list of
    (chrom, start, end, ancestry, matched, lod)."""
    nea = merge([(c, s, e, l, {n}) for n in NEANDERTHALS for c, s, e, l in per_archaic.get(n, [])])
    den = merge([(c, s, e, l, {"Denisova"}) for c, s, e, l in per_archaic.get("Denisova", [])])
    rows = []
    for c, s, e, l, who in nea:
        also = {"Denisova"} if overlaps((c, s, e), den) else set()
        rows.append((c, s, e, "Neanderthal", sorted(who) + sorted(also), l))
    for c, s, e, l, who in den:
        if not overlaps((c, s, e), nea):
            rows.append((c, s, e, "Denisovan", ["Denisova"], l))
    rows.sort(key=lambda r: (chrom_key(r[0]), r[1]))
    return rows


def summarise(rows, spans_bp):
    out = {}
    for anc in ("Neanderthal", "Denisovan"):
        r = [x for x in rows if x[3] == anc]
        bp = sum(x[2] - x[1] for x in r)
        key = anc.lower()
        out[f"{key}_segments"] = len(r)
        out[f"{key}_bp"] = bp
        out[f"{key}_percent"] = f"{100.0 * bp / spans_bp:.2f}" if spans_bp else "0.00"
        out[f"{key}_longest_bp"] = max((x[2] - x[1] for x in r), default=0)
    return out


def cmd_report(a):
    spans = []
    with open(a.spans) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head[:4] != ["chrom", "first", "last", "sites"]:
            die(f"{a.spans} is not the spans table of the genotypes step")
        for line in f:
            c = line.rstrip("\n").split("\t")
            spans.append((c[0], int(c[1]), int(c[2]), int(c[3])))
    if not spans:
        die(f"{a.spans} lists no chromosome")
    spans_bp = sum(last - first + 1 for _, first, last, _ in spans)
    per = {n: [] for n in ARCHAICS}
    for chrom, *_ in spans:
        for n in ARCHAICS:
            p = os.path.join(a.ibdmix_dir, f"{chrom}.{n}.txt")
            if not os.path.isfile(p):
                die(f"no IBDmix output {p}")
            per[n] += [(f"chr{c}" if not str(c).startswith("chr") else c, s, e, l)
                       for c, s, e, l in keep(read_segments(p))]
    rows = classify(per)
    genes = read_genes(a.genes)
    summary = summarise(rows, spans_bp)
    os.makedirs(a.out_dir, exist_ok=True)
    seg_path = os.path.join(a.out_dir, f"{a.sample}_archaic_segments.tsv")
    with open(seg_path, "w") as f:
        f.write("\t".join(SEGMENT_HEADER) + "\n")
        for c, s, e, anc, who, lod in rows:
            g = genes_in(genes, c, s, e) if a.genes else []
            f.write(f"{c}\t{s}\t{e}\t{e - s}\t{anc}\t{','.join(who)}\t{lod:.1f}\t{','.join(g) if g else '.'}\n")
    longest = sorted(rows, key=lambda r: (-(r[2] - r[1]), chrom_key(r[0]), r[1]))[:TOP]
    meta = {
        "sample": a.sample,
        "method": "IBDmix",
        "archaic_genomes": "Altai Neanderthal, Vindija 33.19 Neanderthal, Altai Denisovan",
        "population": a.pop,
        "min_lod": f"{MIN_LOD:g}",
        "min_length_bp": MIN_LENGTH,
        "analysed_bp": spans_bp,
        "analysed_chromosomes": len(spans),
        "analysed_sites": sum(s[3] for s in spans),
        "genes": "GENCODE protein-coding" if a.genes else "none (no gene annotation installed)",
        **summary,
    }
    with open(os.path.join(a.out_dir, f"{a.sample}_archaic.tsv"), "w") as f:
        f.write("key\tvalue\n")
        for k, v in meta.items():
            f.write(f"{k}\t{v}\n")
        for i, (c, s, e, anc, who, lod) in enumerate(longest, 1):
            g = genes_in(genes, c, s, e) if a.genes else []
            f.write(f"longest_{i}\t{c}:{s}-{e} {anc} {(e - s) / 1e6:.2f} Mb "
                    f"({'/'.join(who)}); genes: {', '.join(g) if g else 'none'}\n")
    lines = [f"Archaic ancestry of {a.sample} (IBDmix against the Altai and Vindija Neanderthals and the Altai Denisovan)",
             f"  Analysed: {len(spans)} autosomes, {spans_bp / 1e6:.1f} Mb, {meta['analysed_sites']} sites; "
             f"modern frequencies: 1000 Genomes {a.pop}",
             f"  Neanderthal: {summary['neanderthal_percent']}% ({summary['neanderthal_bp'] / 1e6:.1f} Mb in "
             f"{summary['neanderthal_segments']} segments)",
             f"  Denisovan:   {summary['denisovan_percent']}% ({summary['denisovan_bp'] / 1e6:.1f} Mb in "
             f"{summary['denisovan_segments']} segments)"]
    if longest:
        lines.append(f"  The {len(longest)} longest segments:")
        for c, s, e, anc, who, lod in longest:
            g = genes_in(genes, c, s, e) if a.genes else []
            lines.append(f"    {c}:{s}-{e}  {(e - s) / 1e6:.2f} Mb  {anc} ({'/'.join(who)}, LOD {lod:.0f})  "
                         f"genes: {', '.join(g) if g else 'none'}")
    lines.append(f"  {CAVEAT}")
    lines.append(f"  Every segment: {os.path.basename(seg_path)}")
    with open(os.path.join(a.out_dir, f"{a.sample}_archaic.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("table", help="one GRCh37 chromosome of the archaic table (setup.sh --archaic)")
    t.add_argument("--chrom", required=True)
    t.add_argument("--modern", required=True)
    t.add_argument("--archaic", action="append", default=[], help="NAME=FILE, once each for Altai, Vindija, Denisova")
    t.add_argument("--mask", action="append", default=[],
                   help="NAME=BED: read that archaic genome only inside this mask (its MPI-EVA FilterBed)")
    t.add_argument("--exclude", action="append", default=[],
                   help="BED: leave out every site inside these regions (segmental duplications, indel neighbourhoods)")
    t.add_argument("--out", required=True)
    m = sub.add_parser("strict", help="mark the rows outside the GRCh38 strict mask masked (setup.sh --archaic)")
    m.add_argument("--chrom", required=True)
    m.add_argument("--strict", required=True, help="BED: the 1000 Genomes GRCh38 strict mask")
    m.add_argument("--rows", required=True, help="the lifted rows, tab-separated, no header")
    m.add_argument("--out", required=True)
    g = sub.add_parser("genotypes", help="IBDmix's genotype files for one sample (step 43)")
    g.add_argument("--table-dir", required=True)
    g.add_argument("--vcf", required=True)
    g.add_argument("--pop", default="ALL")
    g.add_argument("--pseudo", type=int, default=PSEUDO)
    g.add_argument("--out-dir", required=True)
    r = sub.add_parser("report", help="segments, shares and the ten longest (step 43)")
    r.add_argument("--sample", required=True)
    r.add_argument("--ibdmix-dir", required=True)
    r.add_argument("--spans", required=True)
    r.add_argument("--genes", default="")
    r.add_argument("--pop", default="ALL")
    r.add_argument("--out-dir", required=True)
    a = p.parse_args(argv)
    if a.cmd == "genotypes" and a.pseudo < 1:
        die("--pseudo must be 1 or more", 2)
    {"table": cmd_table, "strict": cmd_strict, "genotypes": cmd_genotypes, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
