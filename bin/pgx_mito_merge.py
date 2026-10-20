#!/usr/bin/env python3
"""pgx_mito_merge.py: an MT-RNR1 outside call for PharmCAT from the chrM calls.

PharmCAT 3.4.0 does not match MT-RNR1 from a VCF: the gene has no allele
definition for its matcher and is read only as an outside call (-po), one
allele, since chrM is not a pair. So the gene was No Result on every sample,
although step 20 (MITO_VARIANTS in Nextflow) calls chrM with Mutect2 at high
depth. This script reads those calls and adds one MT-RNR1 line to the
outside-call file step 07 and PHARMCAT give PharmCAT.

PharmCAT's own MT-RNR1 alleles are read from the jar
(org/pharmgkb/pharmcat/phenotype/MT_RNR1.json): "Reference" and names such as
m.1555A>G. Of the named single-base changes, a chrM record with the same
position, REF and ALT counts when its FILTER is PASS (or .) and its allele
fraction (FORMAT AF, else INFO AF) is MIN_AF or more: homoplasmic. Mutect2
writes 0/1 even at AF 1.0, so the genotype is not read.

  - One homoplasmic named allele: that allele is the call.
  - More than one: the one PharmCAT ranks highest (increased risk, then
    uncertain, then normal) is the call, and the others are listed: PharmCAT
    takes one MT-RNR1 allele.
  - None, and the file has at least one PASS chrM record (so chrM was
    sequenced and called): "Reference".
Nothing is passed, and the gene stays No Result, with a line saying why,
when a named allele is heteroplasmic (AF below MIN_AF), filtered, or has no
allele fraction; when a homoplasmic PASS record at a named position has
another ALT; when a PASS insertion or deletion overlaps a named position or
lies within 10 bases of a named deletion (m.960 and m.961, in a C run whose
calls this script does not map to PharmCAT's names); and when the file has no
PASS chrM record.

  pgx_mito_merge.py --jar pharmcat.jar --mito chrM.vcf.gz [--outside calls.tsv]
                    [--min-af 0.95] --out merged.tsv

--outside is step 36's outside-call file; its lines are copied, except an
MT-RNR1 line (this script's call replaces it). The output is written even
when it is empty. A summary goes to stdout. Exit 2 when the jar has no
MT-RNR1 definitions or the VCF cannot be read. Standard library only: it runs
in the PharmCAT image and under tests/test_pgx_mito_merge.py.
"""
import argparse
import gzip
import json
import re
import sys
import zipfile

MIN_AF = 0.95
GENE = "MT-RNR1"
RESOURCE = "org/pharmgkb/pharmcat/phenotype/MT_RNR1.json"
CHROMS = {"chrM", "MT", "M", "chrMT"}
NEAR_DELETION = 10
SNV_NAME = re.compile(r"^m\.(\d+)([ACGT])>([ACGT])$")
ANY_NAME = re.compile(r"^m\.(\d+)")
RANK = ("increased", "uncertain", "normal")


def say(msg):
    print(f"pgx_mito_merge: {msg}", flush=True)


def die(msg):
    print(f"pgx_mito_merge: ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(2)


def read_alleles(jar):
    """({(pos, ref, alt): (name, function)}, [positions of other named alleles],
    has_reference) from PharmCAT's MT-RNR1 definitions in the jar."""
    try:
        with zipfile.ZipFile(jar) as z:
            data = json.loads(z.read(RESOURCE))
    except (OSError, KeyError, ValueError, zipfile.BadZipFile) as e:
        die(f"cannot read {RESOURCE} from {jar} ({type(e).__name__}: {e})")
    snv, other, has_ref = {}, [], False
    for a in data.get("namedAlleles") or []:
        name = (a or {}).get("name") or ""
        func = (a or {}).get("functionValue") or ""
        if name == "Reference":
            has_ref = True
            continue
        m = SNV_NAME.match(name)
        if m:
            snv[(int(m.group(1)), m.group(2), m.group(3))] = (name, func)
        elif ANY_NAME.match(name):
            other.append(int(ANY_NAME.match(name).group(1)))
    if not snv:
        die(f"{RESOURCE} in {jar} names no MT-RNR1 single-base allele")
    return snv, sorted(set(other)), has_ref


def open_text(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def allele_fractions(info, fmt, sample, n_alt):
    """The AF of each ALT: FORMAT AF of the first sample, else INFO AF; None
    where there is none."""
    keys = fmt.split(":") if fmt else []
    vals = sample.split(":") if sample else []
    raw = None
    if "AF" in keys and keys.index("AF") < len(vals):
        raw = vals[keys.index("AF")]
    else:
        for kv in info.split(";"):
            if kv.startswith("AF="):
                raw = kv[3:]
    out = []
    for i in range(n_alt):
        try:
            out.append(float(raw.split(",")[i]))
        except (AttributeError, IndexError, ValueError):
            out.append(None)
    return out


def read_chrm(path):
    """[(pos, ref, alt, filter, af)] of the chrM records, one per ALT, and the
    count of PASS chrM records."""
    recs, n_pass = [], 0
    try:
        with open_text(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 8:
                    die(f"{path}: a record with {len(c)} columns")
                if c[0] not in CHROMS:
                    continue
                filt = c[6]
                ok = filt in ("PASS", ".")
                n_pass += ok
                alts = c[4].split(",")
                afs = allele_fractions(c[7], c[8] if len(c) > 8 else "", c[9] if len(c) > 9 else "", len(alts))
                for alt, af in zip(alts, afs):
                    if alt in ("*", ".") or alt.startswith("<"):
                        continue
                    recs.append((int(c[1]), c[3], alt, "PASS" if ok else filt, af))
    except (OSError, EOFError, ValueError) as e:
        die(f"cannot read {path} ({type(e).__name__}: {e})")
    return recs, n_pass


def rank(func):
    f = func.lower()
    return next((i for i, w in enumerate(RANK) if f.startswith(w)), len(RANK))


def decide(recs, n_pass, snv, other, has_ref, min_af):
    """(call or None, [lines]) for MT-RNR1."""
    named_pos = {p for p, _, _ in snv}
    calls, held, lines = [], [], []
    for pos, ref, alt, filt, af in recs:
        shown = f"chrM:{pos} {ref}>{alt}" + (f" AF {af:.3f}" if af is not None else " (no AF)")
        if len(ref) == 1 and len(alt) == 1:
            hit = snv.get((pos, ref, alt))
            if hit is None:
                if pos in named_pos and filt == "PASS" and (af is None or af >= min_af):
                    held.append(f"{shown}: a change at an MT-RNR1 position that PharmCAT does not name")
                continue
            name, func = hit
            if filt != "PASS":
                held.append(f"{name} ({shown}) is filtered ({filt})")
            elif af is None:
                held.append(f"{name} ({shown}) has no allele fraction")
            elif af < min_af:
                held.append(f"{name} ({shown}) is heteroplasmic (AF below {min_af:g}); PharmCAT takes "
                            "one homoplasmic allele")
            else:
                calls.append((rank(func), pos, name, func, shown))
            continue
        if filt != "PASS":
            continue
        start = pos
        end = pos + max(len(ref), 2) - 1  # an insertion sits between pos and pos + 1
        near = [p for p in named_pos if start <= p <= end] + \
            [p for p in other if start - NEAR_DELETION <= p <= end + NEAR_DELETION]
        if near:
            held.append(f"an insertion or deletion ({shown}) at or near m.{min(near)}, which this script "
                        "does not map to PharmCAT's names")
    if n_pass == 0:
        lines.append("the chrM call set has no PASS record (chrM not sequenced or not called): "
                     f"no {GENE} call passed, it stays No Result")
        return None, lines
    if held:
        lines += [f"held back: {h}" for h in held]
        lines.append(f"no {GENE} call passed to PharmCAT: it stays No Result; read the records above "
                     "in the chrM call set")
        return None, lines
    if calls:
        calls.sort()
        best = calls[0]
        lines.append(f"{GENE} {best[2]} passed to PharmCAT as an outside call (homoplasmic, {best[4]}; "
                     f"{best[3]})")
        for c in calls[1:]:
            lines.append(f"also homoplasmic, not passed (PharmCAT takes one {GENE} allele): {c[2]} "
                         f"({c[4]}; {c[3]})")
        return best[2], lines
    if not has_ref:
        lines.append(f"no named {GENE} allele, and PharmCAT's definitions have no Reference: no call passed")
        return None, lines
    lines.append(f"{GENE} Reference passed to PharmCAT as an outside call: none of PharmCAT's "
                 f"{len(snv)} named single-base alleles is called, {n_pass} PASS chrM records read")
    return "Reference", lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jar", required=True, help="pharmcat.jar (its MT-RNR1 definitions)")
    ap.add_argument("--mito", required=True, help="chrM VCF of step 20 / MITO_VARIANTS (bgzip or plain)")
    ap.add_argument("--outside", help="step 36's outside-call file, copied into the output")
    ap.add_argument("--min-af", type=float, default=MIN_AF, help=f"homoplasmic from this AF (default {MIN_AF})")
    ap.add_argument("--out", required=True, help="the outside-call file PharmCAT reads")
    a = ap.parse_args(argv)

    snv, other, has_ref = read_alleles(a.jar)
    recs, n_pass = read_chrm(a.mito)
    call, lines = decide(recs, n_pass, snv, other, has_ref, a.min_af)
    for line in lines:
        say(line)

    out = []
    if a.outside:
        try:
            with open(a.outside) as f:
                for line in f:
                    if not line.strip():
                        continue
                    if line.split("\t", 1)[0].strip() == GENE:
                        say(f"dropped the {GENE} line of {a.outside}: the chrM calls decide it")
                        continue
                    out.append(line if line.endswith("\n") else line + "\n")
        except OSError as e:
            die(f"cannot read {a.outside} ({e})")
    if call:
        out.append(f"{GENE}\t{call}\n")
    with open(a.out, "w") as f:
        f.writelines(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
