#!/usr/bin/env python3
"""Write the compound-het phasing fixture: made-up reads on a made-up contig.

  make_fixture.py OUT_DIR [--swap]

writes, for the sample CHPH on contig chrP (80 kb of seeded random sequence):

  reads.sam    read pairs of 100 bp, 500 bp fragments (mates 400 bp apart)
  vep.vcf      the six heterozygous variants with a VEP CSQ, as step 13 writes
               them, for step 31 (all rare, HIGH impact, three genes)
  comphet.vcf  the three pairs as slivar compound-hets names them
               (INFO/slivar_comphet), for bin/comphet_phase.py alone

The three pairs and what the reads show:

  GENEC  SNV at 1000, 2-base deletion at 1003: both on copy 1 (cis); single
         reads cover both
  GENET  SNVs at 5000 and 5400, one on each copy (trans); no single read
         covers both, only a read and its mate do
  GENEF  SNVs at 20000 and 70000, one on each copy; 50 kb apart, so no read
         pair covers both and the phase cannot be read (unknown)

--swap moves GENEC's deletion to copy 2: the same positions, now in trans.
The output is the same for the same arguments; nothing is random at run time.
"""
import os
import random
import sys

SAMPLE = "CHPH"
CONTIG = "chrP"
LENGTH = 80000
READ = 100
FRAG = 500
FIELDS = ("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|CANONICAL|MANE_SELECT|"
          "Existing_variation|MAX_AF|gnomADe_AF|gnomADg_AF|CLIN_SIG")


def reference():
    rng = random.Random(31)
    return "".join(rng.choice("ACGT") for _ in range(LENGTH))


def other(base):
    return "ACGT"["ACGT".index(base) - 1]


def variants(seq, swap=False):
    """[(gene, pos, ref, alt, consequence, copy)]: pos 1-based, copy 0 or 1."""
    snv = lambda pos: (seq[pos - 1], other(seq[pos - 1]))  # noqa: E731
    return [
        ("GENEC", 1000, *snv(1000), "stop_gained", 0),
        ("GENEC", 1003, seq[1002:1005], seq[1002], "frameshift_variant", 1 if swap else 0),
        ("GENET", 5000, *snv(5000), "stop_gained", 0),
        ("GENET", 5400, *snv(5400), "stop_gained", 1),
        ("GENEF", 20000, *snv(20000), "stop_gained", 0),
        ("GENEF", 70000, *snv(70000), "stop_gained", 1),
    ]


def read_on(seq, vs, copy, start):
    """(sequence, cigar) of a READ-base read from 0-based start on one copy."""
    snvs = {v[1] - 1: v[3] for v in vs if v[5] == copy and len(v[2]) == len(v[3])}
    dels = {v[1] + len(v[3]) - 1: len(v[2]) - len(v[3]) for v in vs if v[5] == copy and len(v[2]) > len(v[3])}
    out, ops, p = [], [], start
    while len(out) < READ:
        if p in dels:
            ops.append(("D", dels[p]))
            p += dels[p]
            continue
        out.append(snvs.get(p, seq[p]))
        ops.append(("M", 1))
        p += 1
    cigar, last, n = "", None, 0
    for op, k in ops + [(None, 0)]:
        if op == last:
            n += k
            continue
        if last:
            cigar += f"{n}{last}"
        last, n = op, k
    return "".join(out), cigar


def sam(seq, vs):
    lines = ["@HD\tVN:1.6\tSO:unsorted", f"@SQ\tSN:{CONTIG}\tLN:{LENGTH}", f"@RG\tID:{SAMPLE}\tSM:{SAMPLE}"]
    # Fragment starts (0-based) around each pair: read 1 covers the first
    # position, its mate the second one when they are 400 bp apart.
    starts = {"c": range(920, 991, 7), "t": range(4905, 4996, 10), "f": range(19905, 19996, 10),
              "g": range(69905, 69996, 10)}
    for region, rs in starts.items():
        for copy in (0, 1):
            for s in rs:
                name = f"{region}{copy}_{s}"
                m = s + FRAG - READ
                s1, c1 = read_on(seq, vs, copy, s)
                s2, c2 = read_on(seq, vs, copy, m)
                q = "I" * READ
                lines.append(f"{name}\t99\t{CONTIG}\t{s + 1}\t60\t{c1}\t=\t{m + 1}\t{FRAG}\t{s1}\t{q}\tRG:Z:{SAMPLE}")
                lines.append(f"{name}\t147\t{CONTIG}\t{m + 1}\t60\t{c2}\t=\t{s + 1}\t-{FRAG}\t{s2}\t{q}\tRG:Z:{SAMPLE}")
    return "\n".join(lines) + "\n"


def vcf_head(info):
    return ("##fileformat=VCFv4.2\n"
            f"##contig=<ID={CONTIG},length={LENGTH}>\n"
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            + info
            + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n")


def vep_vcf(vs):
    out = vcf_head('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
                   f'Format: {FIELDS}">\n')
    for gene, pos, ref, alt, cons, _ in vs:
        csq = (f"{alt}|{cons}|HIGH|{gene}|ENSG{gene}|Transcript|ENST{pos}|protein_coding|YES|NM_{pos}.1|"
               f"rs{pos}|0.0001|0.0001|0.0001|")
        out += f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tCSQ={csq}\tGT:GQ:DP:AD\t0/1:50:30:15,15\n"
    return out


def comphet_vcf(vs):
    out = vcf_head('##INFO=<ID=slivar_comphet,Number=.,Type=String,Description="compound hets called by slivar">\n')
    pid = {"GENEC": 1, "GENET": 2, "GENEF": 3}
    for i, (gene, pos, ref, alt, _, _) in enumerate(vs):
        mate = vs[i + 1] if i % 2 == 0 else vs[i - 1]
        tag = f"{SAMPLE}/{gene}/{pid[gene]}/{CONTIG}/{mate[1]}/{mate[2]}/{mate[3]}"
        out += f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tslivar_comphet={tag}\tGT:GQ:DP:AD\t0/1:50:30:15,15\n"
    return out


def write(out_dir, swap=False):
    seq = reference()
    vs = variants(seq, swap)
    os.makedirs(out_dir, exist_ok=True)
    for name, text in (("reads.sam", sam(seq, vs)), ("vep.vcf", vep_vcf(vs)), ("comphet.vcf", comphet_vcf(vs))):
        with open(os.path.join(out_dir, name), "w") as fh:
            fh.write(text)
    return vs


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--swap"]
    if len(args) != 1:
        sys.exit(__doc__.split("\n\n")[1])
    write(args[0], swap="--swap" in sys.argv[1:])
