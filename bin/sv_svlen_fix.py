#!/usr/bin/env python3
"""Give AnnotSV the net length change of a sequence-resolved SV record.

Manta writes a deletion that also inserts a few bases with its REF and ALT
sequences, END at the last deleted base and SVLEN as the deleted length:
REF 137 bp and ALT 2 bp gets SVLEN=-136 and END=POS+136, although the
record shortens the genome by 135 bp. AnnotSV takes the length from SVLEN
and the overlapped coding length from its interval POS-1..END (0-based,
which counts the padding base too), so it reads 137 deleted coding bases,
reports Frameshift "yes" and gives criterion 2E-1 (+0.90, class 4) to an
in-frame change.

For every record whose REF and ALT are both plain sequence (one ALT) of
different lengths, this script sets SVLEN to len(ALT) - len(REF), and for a
net deletion END to POS - 1 + len(REF) - len(ALT): the END AnnotSV itself
writes when it shortens a deletion whose REF is longer than 10 kb, so the
interval it reads is as long as the net change. Symbolic (<DEL>) and
breakend records pass unchanged. A rewritten record keeps its old values in
ORIG_SVLEN and ORIG_END.

Usage: sv_svlen_fix.py IN.vcf[.gz]|- > OUT.vcf
"""
import gzip
import re
import sys

SEQ = re.compile(r"^[ACGTNacgtn]+$")
NEW_INFO = (
    '##INFO=<ID=ORIG_SVLEN,Number=1,Type=Integer,Description="SVLEN before sv_svlen_fix.py set it to the net '
    'length change len(ALT)-len(REF)">',
    '##INFO=<ID=ORIG_END,Number=1,Type=Integer,Description="END before sv_svlen_fix.py set it to '
    'POS-1+len(REF)-len(ALT) for a net deletion, so AnnotSV reads the net length">',
)


def open_text(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def fix(fields):
    """Rewrite one record in place; True when something changed."""
    ref, alt = fields[3], fields[4]
    if "," in alt or not SEQ.match(ref) or not SEQ.match(alt) or len(ref) == len(alt):
        return False
    pos = int(fields[1])
    net = len(alt) - len(ref)
    items = [] if fields[7] in ("", ".") else fields[7].split(";")
    keys = [i.split("=", 1)[0] for i in items]
    info = dict(i.split("=", 1) for i in items if "=" in i)
    svtype = info.get("SVTYPE", "DEL" if net < 0 else "INS")
    changed = []
    if info.get("SVLEN") != str(net):
        changed.append(("SVLEN", str(net)))
    if net < 0 and svtype == "DEL":
        end = str(pos - 1 + len(ref) - len(alt))
        if info.get("END") != end:
            changed.append(("END", end))
    if not changed:
        return False
    for key, value in changed:
        if key in info:
            items.append("ORIG_{0}={1}".format(key, info[key]))
            items[keys.index(key)] = "{0}={1}".format(key, value)
        else:
            items.append("{0}={1}".format(key, value))
    fields[7] = ";".join(items)
    return True


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__.strip().splitlines()[-1] + "\n")
        return 2
    n = 0
    fh = open_text(argv[1])
    out = sys.stdout
    seen = set()
    for line in fh:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                for h in NEW_INFO:
                    if h not in seen:
                        out.write(h + "\n")
            seen.add(line.rstrip("\n"))
            out.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 8 and fix(fields):
            n += 1
        out.write("\t".join(fields) + "\n")
    if fh is not sys.stdin:
        fh.close()
    sys.stderr.write("SV length check: {0} sequence-resolved records set to their net length change\n".format(n))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
