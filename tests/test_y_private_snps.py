#!/usr/bin/env python3
"""bin/y_private_snps.py on a synthetic VCF and synthetic position tables (no
real sample).

Checks:
  1. a chrY call at a position one of the public SNP tables names is never
     private: the subtraction of the public lists is the whole point of the
     file. Both tables are checked, not only the first;
  2. a call at a position no table names is listed, with its depth, ALT
     fraction, genotype quality and genotype;
  3. the gates: a position in a pseudoautosomal region, a record that did
     not pass the filter, a genotype carrying the reference allele, an
     indel, a second ALT allele, depth below --min-depth, genotype quality
     below --min-gq and an ALT fraction below --min-fraction are all
     dropped, each counted in the file's header;
  4. depth far above the sample's own chrY median is dropped (reads of more
     than one region piling up on one place), and --max-depth-factor 0
     turns that gate off;
  5. a gzipped VCF reads the same as a plain one;
  6. a missing position table stops the tool and names setup.sh --yleaf-data;
  7. a call outside the callable mask (a base not marked "P", or past the
     mask's end) is dropped and counted as outside_callable, whatever its
     quality; --callable-mask points at another mask; a missing mask stops
     the tool.

Run: python3 tests/test_y_private_snps.py
"""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(REPO, "bin", "y_private_snps.py")

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


# name -> position, as the two Yleaf GRCh38 tables place them. The columns are
# the tables' own: chromosome, name, haplogroup, position, mutation, ancestral,
# derived (bin/y_tree_place.py reads the same layout).
FTDNA = {"pub1": 7000100}
YFULL = {"pub2": 7000200}

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chrY,length=57227415>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
"""

# pos, ref, alt, filter, GT:AD:DP:GQ
RECORDS = [
    (7000100, "G", "A", "PASS", "1:0,30:30:50"),    # on the first table
    (7000200, "G", "A", "PASS", "1:0,30:30:50"),    # on the second table
    (7000300, "C", "T", "PASS", "1:0,30:30:50"),    # private
    (7000400, "C", "T", "PASS", "1/1:1,29:30:50"),  # private, diploid call
    (2000000, "C", "T", "PASS", "1:0,30:30:50"),    # PAR1
    (56900000, "C", "T", "PASS", "1:0,30:30:50"),   # PAR2
    (7000500, "C", "T", "RefCall", "1:0,30:30:50"),  # not PASS
    (7000600, "C", "T", "PASS", "0/1:15,15:30:50"),  # carries the reference allele
    (7000700, "C", "CTT", "PASS", "1:0,30:30:50"),  # insertion
    (7000800, "C", "T,G", "PASS", "1:0,20,10:30:50"),  # second ALT allele
    (7000900, "C", "T", "PASS", "1:0,9:9:50"),      # depth below the gate
    (7001000, "C", "T", "PASS", "1:0,30:30:5"),     # genotype quality below the gate
    (7001100, "C", "T", "PASS", "1:6,24:30:50"),    # ALT fraction 0.8
    (7001200, "C", "T", "PASS", "1:0,300:300:50"),  # ten times the median depth
    (7001300, "C", "T", "PASS", "1:0,30:30:50"),    # outside the callable mask
]

# The synthetic mask: every base passes ("P") up to MASK_LEN, except a run of
# "Z" (too many zero mapping-quality reads) around 7001300.
MASK_LEN = 7002000
MASK_Z = (7001250, 7001350)


def write_mask(path, z=MASK_Z):
    seq = "P" * (z[0] - 1) + "Z" * (z[1] - z[0] + 1) + "P" * (MASK_LEN - z[1])
    with gzip.open(path, "wt") as o:
        o.write(">chrY\n")
        for i in range(0, len(seq), 60):
            o.write(seq[i:i + 60] + "\n")


def write_tables(d):
    for name, snps in (("positions_ftdna_hg38.txt", FTDNA), ("positions_yfull_hg38.txt", YFULL)):
        with open(os.path.join(d, name), "w") as o:
            for snp, pos in snps.items():
                o.write(f"chrY\t{snp}\tR\t{pos}\tG->A\tG\tA\n")


def write_vcf(path, records=RECORDS):
    out = gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")
    with out as o:
        o.write(HEADER)
        for pos, ref, alt, filt, sample in records:
            o.write(f"chrY\t{pos}\t.\t{ref}\t{alt}\t50\t{filt}\t.\tGT:AD:DP:GQ\t{sample}\n")


def run(d, vcf, prefix, extra=()):
    r = subprocess.run([sys.executable, TOOL, "--vcf", vcf, "--tree-dir", d,
                        "--sample", "S", "--out-prefix", os.path.join(d, prefix), *extra],
                       capture_output=True, text=True)
    out = os.path.join(d, f"{prefix}_y_private_snps.tsv")
    head, rows = [], []
    if os.path.exists(out):
        for line in open(out).read().splitlines():
            if line.startswith("#"):
                head.append(line)
            elif not line.startswith("Pos\t"):
                rows.append(line.split("\t"))
    return r, head, rows


def main():
    d = tempfile.mkdtemp()
    try:
        write_tables(d)
        write_mask(os.path.join(d, "callable_mask_chrY.fasta.gz"))
        vcf = os.path.join(d, "sample.vcf")
        write_vcf(vcf)
        r, head, rows = run(d, vcf, "S")
        pos = [x[0] for x in rows]
        check("the tool exits 0", r.returncode == 0, r.stderr)
        check("1. a position the first table names (pub1) is not private", "7000100" not in pos, pos)
        check("1. a position the second table names (pub2) is not private", "7000200" not in pos, pos)
        check("1. the header counts the two positions dropped as already public",
              any("on_public_lists=2" in h for h in head), head)
        check("2. the private haploid call is listed", "7000300" in pos, pos)
        check("2. the private diploid ALT/ALT call is listed", "7000400" in pos, pos)
        check("2. its row carries the depth, ALT fraction, quality and genotype",
              ["7000300", "C", "T", "30", "1.000", "50", "1"] in rows, rows)
        check("3. no other position is listed", sorted(pos) == ["7000300", "7000400"], pos)
        for label, count in (("PAR=2", "the two pseudoautosomal calls"), ("not_PASS=1", "the filtered call"),
                             ("genotype=1", "the call with a reference allele"),
                             ("depth=1", "the call below --min-depth"),
                             ("quality=1", "the call below --min-gq"),
                             ("alt_fraction=1", "the call below --min-fraction"),
                             ("excess_depth=1", "the call far above the median depth"),
                             ("outside_callable=1", "the call outside the callable mask")):
            check(f"3. the header counts {count} ({label})", any(label in h for h in head), head)
        check("3. the indel and the second-ALT record never reach the gates (2 of 15 records)",
              any("13 chrY SNV records read" in line for line in r.stdout.splitlines()), r.stdout)
        check("4. --max-depth-factor 0 keeps the deep call",
              "7001200" in [x[0] for x in run(d, vcf, "off", ("--max-depth-factor", "0"))[2]], None)
        gz = os.path.join(d, "s.vcf.gz")
        write_vcf(gz)
        check("5. a gzipped VCF gives the same rows",
              [x[0] for x in run(d, gz, "gz")[2]] == ["7000300", "7000400"], None)
        check("the step line says a private SNP is unconfirmed until another man carries it",
              "unconfirmed until another man" in r.stdout, r.stdout)

        check("7. the call outside the callable mask is not listed", "7001300" not in pos, pos)
        check("7. the header names the mask it used",
              any("callable_mask=callable_mask_chrY.fasta.gz" in h for h in head), head)
        other = os.path.join(d, "other_mask.fa.gz")
        write_mask(other, z=(7000250, 7000350))
        alt = [x[0] for x in run(d, vcf, "alt", ("--callable-mask", other))[2]]
        check("7. --callable-mask reads another mask (7000300 masked there, 7001300 callable)",
              alt == ["7000400", "7001300"], alt)
        os.remove(os.path.join(d, "callable_mask_chrY.fasta.gz"))
        r = subprocess.run([sys.executable, TOOL, "--vcf", vcf, "--tree-dir", d, "--sample", "S",
                            "--out-prefix", os.path.join(d, "nomask")], capture_output=True, text=True)
        check("7. a missing mask stops the tool and names setup.sh --yleaf-data",
              r.returncode != 0 and "setup.sh --yleaf-data" in r.stderr, r.stderr)
        write_mask(os.path.join(d, "callable_mask_chrY.fasta.gz"))

        os.remove(os.path.join(d, "positions_yfull_hg38.txt"))
        r = subprocess.run([sys.executable, TOOL, "--vcf", vcf, "--tree-dir", d, "--sample", "S",
                            "--out-prefix", os.path.join(d, "missing")], capture_output=True, text=True)
        check("6. a missing table stops the tool and names setup.sh --yleaf-data",
              r.returncode != 0 and "setup.sh --yleaf-data" in r.stderr, r.stderr)
    finally:
        shutil.rmtree(d)
    print(f"\n{'FAILED' if FAILS else 'OK'}: {FAILS} failed check(s)")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
