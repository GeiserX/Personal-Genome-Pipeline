#!/usr/bin/env python3
"""Y-STR repeat counts from HipSTR's calls, in FamilyTreeDNA's marker order (step 45).

Two commands, both run in HIPSTR_IMAGE (its Python, standard library only):

  ystr_report.py bed --panel assets/ystr_panel.tsv --out regions.bed
      The HipSTR region file: every region of the panel's single and parts
      markers, one line each, sorted by position.

  ystr_report.py report --panel assets/ystr_panel.tsv --vcf <id>_y_str.vcf.gz \
      --tsv <id>_y_str.tsv --csv <id>_y_str.csv [--log <id>_y_str.log] [--read-length 150]
      One row per panel marker, in the panel's order. With HipSTR's log, a
      region HipSTR wrote no record for carries the reason the log gives.

How a count is made. HipSTR writes, per region, GB: the base pairs the call
differs from the reference allele. The reference copies of a region are the
panel's (rounded to whole repeats: the region edges of the marker set are
not always on a repeat boundary; halves round up). A region's count is those copies plus
GB / motif length. A marker of several parts is the sum of its parts. A span
marker is one HipSTR region from its first part's start to its last part's
end: HipSTR genotypes each region alone, and on these close-set parts it
wrote no record for one part on real reads, so the step asks for the whole
locus at once; its count is its parts' reference copies plus that region's GB. When
the base pairs left over are not a whole repeat, the count carries them after
a dot, the usual notation of a partial repeat (11.2: eleven repeats and two
bases). The panel's offset, when it gives one, is added last.

Confidence per marker (the worst of its parts):
  high       HipSTR's own haploid filters pass (posterior Q >= 0.9, at most 15 %
             of the reads with a stutter artifact and at most 15 % with an indel
             in the flanks: HipSTR's README) and at least 10 reads were used
  low        a call that fails one of those; the note says which
  no call    HipSTR wrote no record or no genotype for a part
  not callable  multi-copy, absent from the marker set, or a region longer than
             the reads can span (read length minus 10 bases of flank on each side)
"""
import argparse
import csv
import gzip
import sys

MIN_READS = 10        # reads HipSTR used, for a call to count as high confidence
MIN_Q = 0.9           # HipSTR README, filter_haploid_vcf.py --min-call-qual
MAX_STUTTER = 0.15    # --max-call-stutter
MAX_FLANK = 0.15      # --max-call-flank-indel
FLANK = 10            # bases of flank a read needs on each side of the repeat

TSV_COLUMNS = ["order", "panel", "marker", "count", "confidence", "reads", "stutter_reads",
               "flank_indel_reads", "q", "note"]


class PanelError(Exception):
    pass


def read_panel(path):
    """The panel's rows, with each region parsed: (chrom, start, end, period, copies, name)."""
    rows = []
    with open(path, encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f if line.strip() and not line.startswith("#")]
    if not lines:
        raise PanelError(f"{path}: no header")
    header = lines[0].split("\t")
    need = ["order", "panel", "marker", "kind", "regions", "motif", "offset", "note"]
    if header != need:
        raise PanelError(f"{path}: the header is not {' '.join(need)}")
    for n, line in enumerate(lines[1:], 2):
        c = line.split("\t")
        if len(c) != len(need):
            raise PanelError(f"{path}: row {n} has {len(c)} fields, not {len(need)}")
        r = dict(zip(need, c))
        if r["kind"] not in ("single", "parts", "span", "multicopy", "absent"):
            raise PanelError(f"{path}: {r['marker']}: unknown kind {r['kind']}")
        regs = []
        if r["regions"] != ".":
            for part in r["regions"].split(";"):
                try:
                    where, period, copies, name = part.split("/", 3)
                    chrom, span = where.split(":")
                    start, end = (int(x) for x in span.split("-"))
                    regs.append((chrom, start, end, int(period), float(copies), name))
                except ValueError:
                    raise PanelError(f"{path}: {r['marker']}: unreadable region {part}")
        if (r["kind"] in ("single", "parts", "span")) != bool(regs) or (r["kind"] == "single" and len(regs) != 1) \
                or (r["kind"] in ("parts", "span") and len(regs) < 2):
            raise PanelError(f"{path}: {r['marker']}: kind {r['kind']} with {len(regs)} regions")
        if len({p for _c, _s, _e, p, _n, _m in regs}) > 1:
            raise PanelError(f"{path}: {r['marker']}: parts with different motif lengths")
        if r["offset"] != ".":
            try:
                r["offset"] = int(r["offset"])
            except ValueError:
                raise PanelError(f"{path}: {r['marker']}: offset {r['offset']} is not a whole number")
        if r["kind"] == "span":
            if len({c for c, *_x in regs}) > 1 or any(b[1] <= a[2] for a, b in zip(regs, regs[1:])):
                raise PanelError(f"{path}: {r['marker']}: span parts not in order on one chromosome")
            # One region over the whole locus, named after the marker; its copies are the
            # parts' (each rounded, halves up), the repeats the bases between them leave out.
            r["parts"] = regs
            regs = [(regs[0][0], regs[0][1], regs[-1][2], regs[0][3],
                     float(sum(int(x[4] + 0.5) for x in regs)), r["marker"])]
        r["regions"] = regs
        rows.append(r)
    orders = [int(r["order"]) for r in rows]
    if orders != list(range(1, len(rows) + 1)):
        raise PanelError(f"{path}: the order column is not 1 to {len(rows)}")
    return rows


def write_bed(panel, out):
    seen = {}
    for r in panel:
        for chrom, start, end, period, copies, name in r["regions"]:
            seen[name] = (chrom, start, end, period, copies, name)
    with open(out, "w", encoding="utf-8") as f:
        for chrom, start, end, period, copies, name in sorted(seen.values(), key=lambda x: (x[0], x[1], x[2])):
            f.write(f"{chrom}\t{start}\t{end}\t{period}\t{copies:g}\t{name}\n")


def read_calls(path):
    """HipSTR's records by ID: {name: {FORMAT key: value}} of the one sample."""
    calls = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        samples = None
        for line in f:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = c[9:]
                if len(samples) != 1:
                    raise PanelError(f"{path}: {len(samples)} samples; step 45 reads one")
                continue
            if samples is None or len(c) < 10:
                raise PanelError(f"{path}: a record before the header, or one without a sample column")
            calls[c[2]] = dict(zip(c[8].split(":"), c[9].split(":")))
    if samples is None:
        raise PanelError(f"{path}: no #CHROM header line")
    return calls


LOG_REASONS = ("Skipping", "Aborting", "Stutter model training failed", "Haplotype construction failed",
               "Retraining stutter model training failed")


def read_log(path):
    """HipSTR's reason per region: {(chrom, 1-based start, end): first reason line}.

    HipSTR logs "Processing region <chrom> <0-based start> <end>" before each
    region, then, when it gives up on one, a line from LOG_REASONS."""
    why, cur = {}, None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("Processing region "):
                c = line.split()
                try:
                    cur = (c[2], int(c[3]) + 1, int(c[4]))
                except (IndexError, ValueError):
                    cur = None
            elif cur and cur not in why and line.startswith(LOG_REASONS):
                why[cur] = line
    return why


def number(v, kind=float):
    try:
        return kind(v)
    except (TypeError, ValueError):
        return None


def part_call(fmt, period):
    """(bp difference, reads, stutter, flank indels, Q, problems) of one region's record."""
    gt, gb = fmt.get("GT", "."), fmt.get("GB", ".")
    if gt in (".", "./.", ".|.") or gb in (".", "./.", ".|."):
        why = fmt.get("FILTER")
        return None, f"HipSTR made no genotype{' (' + why + ')' if why and why not in ('.', 'PASS') else ''}"
    lens = {number(x, int) for x in gb.replace("|", "/").split("/")}
    if None in lens:
        return None, f"unreadable GB {gb}"
    if len(lens) > 1:
        # Two lengths on a haploid chromosome: a duplicated locus, or noise.
        return None, f"two lengths on the Y chromosome (GB {gb})"
    bp = lens.pop()
    dp = number(fmt.get("DP"), int) or 0
    st = number(fmt.get("DSTUTTER"), int) or 0
    fl = number(fmt.get("DFLANKINDEL"), int) or 0
    q = number(fmt.get("Q"))
    problems = []
    if q is None or q < MIN_Q:
        problems.append(f"posterior {fmt.get('Q', '.')} below {MIN_Q}")
    if dp and st / dp > MAX_STUTTER:
        problems.append(f"{st} of {dp} reads with stutter")
    if dp and fl / dp > MAX_FLANK:
        problems.append(f"{fl} of {dp} reads with a flank indel")
    if dp < MIN_READS:
        problems.append(f"{dp} reads, fewer than {MIN_READS}")
    return (bp, dp, st, fl, q, problems), None


def count_text(ref_copies, bp, period, offset):
    whole, rest = divmod(bp, period)
    n = ref_copies + whole + (offset if isinstance(offset, int) else 0)
    return f"{n}.{rest}" if rest else str(n)


def marker_row(r, calls, read_length, log=None):
    row = {"order": r["order"], "panel": r["panel"], "marker": r["marker"], "count": "", "confidence": "",
           "reads": "", "stutter_reads": "", "flank_indel_reads": "", "q": "", "note": ""}
    notes = [] if r["note"] == "." else [r["note"]]
    if r["kind"] == "multicopy":
        row.update(confidence="not callable",
                   note="not callable from short reads: a multi-copy marker in duplicated Y sequence")
        return row
    if r["kind"] == "absent":
        row.update(confidence="not callable", note="not callable: no region for it in the HipSTR marker set")
        return row
    longest = max(e - s + 1 for _c, s, e, _p, _n, _m in r["regions"])
    if longest + 2 * FLANK > read_length:
        row.update(confidence="not callable",
                   note=f"not callable from short reads: a {longest} bp repeat is longer than {read_length} bp "
                        f"reads can span with {FLANK} bases on each side")
        return row
    period = r["regions"][0][3]
    bp_sum, ref_sum, dps, sts, fls, qs, problems = 0, 0, [], [], [], [], []
    for chrom, start, end, p, copies, name in r["regions"]:
        fmt = calls.get(name)
        if fmt is None:
            if log is None:
                why = f"HipSTR wrote no record for {name}; see HipSTR's log"
            elif (chrom, start, end) in log:
                why = f"HipSTR wrote no record for {name}; its log: {log[(chrom, start, end)]}"
            else:
                why = f"HipSTR wrote no record for {name}, and its log gives no reason"
            row.update(confidence="no call", note="; ".join(notes + [why]))
            return row
        got, why = part_call(fmt, p)
        if got is None:
            row.update(confidence="no call", note="; ".join(notes + [f"{name}: {why}"]))
            return row
        bp, dp, st, fl, q, probs = got
        bp_sum += bp
        ref_sum += int(copies + 0.5)  # halves up
        dps.append(dp)
        sts.append(st)
        fls.append(fl)
        qs.append(q)
        problems += [f"{name}: {x}" if len(r["regions"]) > 1 else x for x in probs]
    row.update(count=count_text(ref_sum, bp_sum, period, r["offset"]),
               confidence="low" if problems else "high",
               reads=str(min(dps)), stutter_reads=str(max(sts)), flank_indel_reads=str(max(fls)),
               q=f"{min(qs):g}" if None not in qs else ".",
               note="; ".join(notes + problems) or ".")
    return row


def report(args):
    panel = read_panel(args.panel)
    calls = read_calls(args.vcf)
    log = read_log(args.log) if args.log else None
    rows = [marker_row(r, calls, args.read_length, log) for r in panel]
    with open(args.tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(TSV_COLUMNS) + "\n")
        for row in rows:
            f.write("\t".join(row[k].replace("\t", " ") for k in TSV_COLUMNS) + "\n")
    with open(args.csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["marker", "value", "confidence"])
        for row in rows:
            w.writerow([row["marker"], row["count"], row["confidence"]])
    tally = {}
    for row in rows:
        tally[row["confidence"]] = tally.get(row["confidence"], 0) + 1
    print("Y-STR markers: " + ", ".join(f"{tally.get(k, 0)} {k}" for k in ("high", "low", "no call", "not callable")))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("bed")
    b.add_argument("--panel", required=True)
    b.add_argument("--out", required=True)
    r = sub.add_parser("report")
    r.add_argument("--panel", required=True)
    r.add_argument("--vcf", required=True)
    r.add_argument("--tsv", required=True)
    r.add_argument("--csv", required=True)
    r.add_argument("--log", help="HipSTR's --log file, for the reason a region has no record")
    r.add_argument("--read-length", type=int, default=150)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "bed":
            write_bed(read_panel(args.panel), args.out)
        else:
            report(args)
    except PanelError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
