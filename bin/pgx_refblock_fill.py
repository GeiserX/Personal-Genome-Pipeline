#!/usr/bin/env python3
"""pgx_refblock_fill.py: call PharmCAT's multi-base positions 0/0 from a gVCF.

Step 07 and PHARMCAT_PREPROCESS expand the reference blocks of a gVCF over
PharmCAT's gene regions with `bcftools convert --gvcf2vcf`, which writes one
record per base, each with a one-base REF. At a position whose REF or ALT is
longer than one base (an insertion or deletion such as UGT1A1 *28, or a
multi-base substitution) PharmCAT's preprocessor sets such a record aside and
lists the position as missing; where PharmCAT's REF is longer than one base
its matcher then discards the record, so the alleles defined there are
uncallable, even when the gVCF has confident reference blocks over the whole
REF span.

This filter reads that expanded VCF on stdin and writes it to stdout with one
change: for each such PharmCAT position where the gVCF shows the sample is
reference, it writes a 0/0 record with PharmCAT's own ID, REF and ALT (INFO
flag PGX_REFBLOCK), in place of the one-base record at the same position. The
gVCF shows it when
  - every base of the REF span lies in reference blocks: records whose ALT is
    only <*>, <X>, <NON_REF> or ., whose genotype is all 0 (0/0, or 0 on a
    haploid chrX), whose FILTER is PASS or ., and whose GQ is MIN_GQ or more;
  - and no other record overlaps the span: no variant record, no ./. block,
    no block below MIN_GQ.
Otherwise the position stays as it was. The record's GT is the block's at the
first base; its GQ is the lowest GQ over the span.

  pgx_refblock_fill.py [--min-gq 20] POSITIONS GVCF < expanded.vcf > out.vcf

POSITIONS is PharmCAT's positions VCF (pharmcat_positions.vcf.bgz in the
image, bgzip or plain). GVCF is the sample's gVCF as plain VCF text, one
sample, at least the records over PharmCAT's regions (bcftools view -R
pharmcat_regions.bed). A summary goes to stderr. Exit 2 on unreadable input.
Standard library only: it runs in the PharmCAT image and under
tests/test_pgx_refblock_fill.py.
"""
import argparse
import bisect
import gzip
import re
import sys

MIN_GQ = 20
BLOCK_ALTS = {"<*>", "<X>", "<NON_REF>", "."}
FLAG = "PGX_REFBLOCK"
FLAG_HEADER = (f'##INFO=<ID={FLAG},Number=0,Type=Flag,Description="0/0 written by '
               'pgx_refblock_fill.py: reference blocks cover the whole REF span">\n')


def open_text(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_positions(path):
    """{chrom: [(start, end, id, ref, alt, gene)]} of the positions whose REF
    or an ALT is longer than one base, sorted by start; and their count."""
    out, n = {}, 0
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            ref, alts = c[3], c[4].split(",")
            if len(ref) == 1 and all(len(a) == 1 for a in alts):
                continue
            gene = next((kv[3:] for kv in (c[7] if len(c) > 7 else "").split(";") if kv.startswith("PX=")), "")
            pos = int(c[1])
            out.setdefault(c[0], []).append((pos, pos + len(ref) - 1, c[2], ref, c[4], gene))
            n += 1
    for v in out.values():
        v.sort()
    return out, n


def is_block(alt, filt, gt, gq, min_gq):
    """A confident reference block (or single-base reference record)."""
    if not set(alt.split(",")) <= BLOCK_ALTS:
        return False
    alleles = re.split(r"[/|]", gt)
    if not alleles or any(a != "0" for a in alleles):
        return False
    if filt not in (".", "PASS"):
        return False
    return gq is not None and gq >= min_gq


def scan_gvcf(path, targets, min_gq):
    """Walk the gVCF. Returns (contig order, {(chrom, start): (gt, min gq)})
    for every target whose span is fully inside confident blocks and touched
    by no other record."""
    order = []
    state = {}  # (chrom, start) -> [covered bases, bad, gt at start, min gq]
    maxlen = {ch: max(t[1] - t[0] for t in ts) for ch, ts in targets.items()}
    starts = {ch: [t[0] for t in ts] for ch, ts in targets.items()}
    samples = None
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                m = re.match(r"##contig=<ID=([^,>]+)", line)
                if m and m.group(1) not in order:
                    order.append(m.group(1))
                continue
            if line.startswith("#"):
                samples = line.rstrip("\n").split("\t")[9:]
                if len(samples) != 1:
                    raise ValueError(f"the gVCF must hold one sample, it holds {len(samples)}")
                continue
            if samples is None:
                raise ValueError("the gVCF has no #CHROM header line")
            c = line.rstrip("\n").split("\t")
            chrom = c[0]
            if chrom not in order:
                order.append(chrom)
            if chrom not in targets:
                continue
            pos = int(c[1])
            m = re.search(r"(?:^|;)END=(\d+)", c[7])
            end = int(m.group(1)) if m else pos + len(c[3]) - 1
            ts, st = targets[chrom], starts[chrom]
            j = bisect.bisect_right(st, end) - 1
            hit = []
            while j >= 0 and st[j] >= pos - maxlen[chrom]:
                if ts[j][1] >= pos:
                    hit.append(ts[j])
                j -= 1
            if not hit:
                continue
            keys = c[8].split(":") if len(c) > 8 else []
            vals = c[9].split(":") if len(c) > 9 else []
            fmt = dict(zip(keys, vals))
            gt = fmt.get("GT", ".")
            try:
                gq = int(fmt["GQ"])
            except (KeyError, ValueError):
                gq = None
            good = is_block(c[4], c[6], gt, gq, min_gq)
            for t in hit:
                s = state.setdefault((chrom, t[0]), [set(), False, None, None])
                if not good:
                    s[1] = True
                    continue
                lo, hi = max(pos, t[0]), min(end, t[1])
                s[0].update(range(lo, hi + 1))
                if lo == t[0]:
                    s[2] = gt
                s[3] = gq if s[3] is None else min(s[3], gq)
    if samples is None:
        raise ValueError("the gVCF has no #CHROM header line")
    ok = {}
    for chrom, ts in targets.items():
        for t in ts:
            s = state.get((chrom, t[0]))
            if s and not s[1] and len(s[0]) == t[1] - t[0] + 1 and s[2]:
                ok[(chrom, t[0])] = (s[2], s[3])
    return order, ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--min-gq", type=int, default=MIN_GQ, help=f"lowest block GQ that counts (default {MIN_GQ})")
    p.add_argument("positions", help="PharmCAT's positions VCF")
    p.add_argument("gvcf", help="the sample's gVCF records over PharmCAT's regions (plain VCF)")
    a = p.parse_args(argv)
    try:
        targets, n = read_positions(a.positions)
        order, ok = scan_gvcf(a.gvcf, targets, a.min_gq)
    except (OSError, ValueError, EOFError) as e:
        print(f"pgx_refblock_fill: ERROR: {e}", file=sys.stderr)
        return 2
    rank = {ch: i for i, ch in enumerate(order)}
    fills = []
    for chrom, ts in targets.items():
        for start, _, vid, ref, alt, gene in ts:
            if (chrom, start) in ok:
                gt, gq = ok[(chrom, start)]
                fills.append(((rank.get(chrom, len(rank)), start),
                              f"{chrom}\t{start}\t{vid}\t{ref}\t{alt}\t.\tPASS\t{FLAG}\tGT:GQ\t{gt}:{gq}\n",
                              f"{gene or '?'} {chrom}:{start}"))
    fills.sort()
    out, i = sys.stdout, 0
    for line in sys.stdin:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                out.write(FLAG_HEADER)
            out.write(line)
            continue
        c = line.split("\t", 5)
        key = (rank.get(c[0], len(rank)), int(c[1]))
        while i < len(fills) and fills[i][0] < key:
            out.write(fills[i][1])
            i += 1
        if i < len(fills) and fills[i][0] == key and len(c[3]) == 1 and c[4] in BLOCK_ALTS:
            out.write(fills[i][1])  # in place of the one-base reference record
            i += 1
            continue
        out.write(line)
    for _, rec, _ in fills[i:]:
        out.write(rec)
    print(f"pgx_refblock_fill: {len(fills)} of {n} PharmCAT positions longer than one base called 0/0 "
          f"from reference blocks with GQ >= {a.min_gq}"
          + (": " + ", ".join(f[2] for f in fills) if fills else "")
          + "; the others stay as the gVCF has them (a variant call, or missing)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
