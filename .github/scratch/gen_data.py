#!/usr/bin/env python3
"""Tiny synthetic inputs for the scratch checks. No real genome data.

ref   OUT.fasta NAME:LEN [NAME:LEN ...]
pairs REF.fasta CONTIG N OUT_R1.fastq.gz OUT_R2.fastq.gz
long  REF.fasta CONTIG N OUT.fastq
"""
import gzip
import random
import sys

random.seed(20261001)
COMP = str.maketrans("ACGT", "TGCA")


def rseq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def mutate(seq, every=4000):
    """Second haplotype: a SNV every `every` bp, so callers have something to call."""
    s = list(seq)
    for pos in range(every // 2, len(s), every):
        s[pos] = {"A": "G", "C": "T", "G": "A", "T": "C"}[s[pos]]
    return "".join(s)


def cmd_ref(out, specs):
    with open(out, "w") as fh:
        for spec in specs:
            name, length = spec.split(":")
            seq = rseq(int(length))
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def cmd_pairs(ref, contig, n, out1, out2):
    hap1 = read_fasta(ref)[contig]
    hap2 = mutate(hap1)
    rl = 150
    with gzip.open(out1, "wt") as f1, gzip.open(out2, "wt") as f2:
        for i in range(int(n)):
            hap = hap1 if i % 2 == 0 else hap2
            ins = max(rl + 50, int(random.gauss(400, 30)))
            start = random.randint(0, len(hap) - ins - 1)
            frag = hap[start:start + ins]
            r1, r2 = frag[:rl], frag[-rl:].translate(COMP)[::-1]
            if random.random() < 0.5:
                r1, r2 = r2, r1
            q = "I" * rl
            f1.write(f"@r{i}/1\n{r1}\n+\n{q}\n")
            f2.write(f"@r{i}/2\n{r2}\n+\n{q}\n")


def cmd_long(ref, contig, n, out):
    seq = read_fasta(ref)[contig]
    with open(out, "w") as fh:
        for i in range(int(n)):
            ln = random.randint(3000, 6000)
            start = random.randint(0, len(seq) - ln - 1)
            r = seq[start:start + ln]
            if i % 2:
                r = r.translate(COMP)[::-1]
            fh.write(f"@lr{i}\n{r}\n+\n{'I' * len(r)}\n")


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"ref": lambda a: cmd_ref(a[0], a[1:]),
     "pairs": lambda a: cmd_pairs(*a),
     "long": lambda a: cmd_long(*a)}[cmd](args)
