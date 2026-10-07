#!/usr/bin/env python3
"""yleaf_run.py: run Yleaf 3.2.1 (YLEAF_IMAGE) offline, on a pileup made beside it.

The Bioconda image of Yleaf ships no samtools, which Yleaf calls for its BAM
input (`samtools idxstats`, `samtools mpileup`), and Yleaf downloads the whole
hg38 FASTA unless its read-only config names one. So step 37 and the
Y_HAPLOGROUP process run it in three parts:

  yleaf_run.py positions OUT
      (YLEAF_IMAGE) write Yleaf's GRCh38 marker positions as the "chrY<TAB>pos"
      list `samtools mpileup -l` reads
  samtools idxstats BAM > IDXSTATS; samtools mpileup -l OUT -AQ20q1 BAM > PILEUP
      (SAMTOOLS_IMAGE) the two commands Yleaf would run, with its default
      quality threshold of 20
  yleaf_run.py predict --bam BAM --reference FASTA --idxstats IDXSTATS --pileup PILEUP --out DIR
      (YLEAF_IMAGE) Yleaf's own flow, with its reference constant set to FASTA
      (a BAM never needs the sequence), each samtools call served from those
      two files, and its multiprocessing pools replaced by a serial map: a
      failed call raises SystemExit in a pool worker, which kills it and
      leaves Pool.map waiting forever.

DIR/hg_prediction.hg is Yleaf's prediction table (Hg is NA when too few
markers had reads). Standard library plus what Yleaf itself imports.
"""
import multiprocessing
import shutil
import sys
from pathlib import Path

QUALITY = 20   # Yleaf's -q default; the mpileup above uses the same


class SerialPool:
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def map(self, fn, items):
        return list(map(fn, items))


def positions(out):
    from yleaf import yleaf_constants as c
    src = c.DATA_FOLDER / c.HG38 / c.NEW_POSITION_FILE
    seen = set()
    with open(src) as f, open(out, "w") as o:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 3 and cols[3].isdigit() and cols[3] not in seen:
                seen.add(cols[3])
                o.write(f"chrY\t{cols[3]}\n")
    if not seen:
        sys.exit(f"ERROR: no marker positions in {src}")
    print(f"{len(seen)} Y marker positions from {src}")


def predict(argv):
    import argparse
    ap = argparse.ArgumentParser(prog="yleaf_run.py predict")
    for a in ("--bam", "--reference", "--idxstats", "--pileup", "--out"):
        ap.add_argument(a, required=True)
    a = ap.parse_args(argv)
    multiprocessing.Pool = SerialPool
    from yleaf import yleaf_constants
    yleaf_constants.HG38_FULL_GENOME = Path(a.reference)
    from yleaf import Yleaf

    def serve(cmd, stdout_location=None):
        if cmd.startswith("samtools idxstats "):
            with open(a.idxstats) as f:
                shutil.copyfileobj(f, stdout_location)
            stdout_location.flush()
        elif cmd.startswith("samtools mpileup ") and " > " in cmd:
            shutil.copyfile(a.pileup, cmd.rsplit(" > ", 1)[1].strip())
        elif cmd.startswith("samtools index "):
            pass   # the BAM is indexed; Yleaf only asks when it finds no .bai
        else:
            sys.exit(f"ERROR: Yleaf asked for a command this launcher does not serve: {cmd}")

    Yleaf.call_command = serve
    sys.argv = ["Yleaf", "-bam", a.bam, "-o", a.out, "-rg", "hg38", "-force", "-t", "1", "-q", str(QUALITY)]
    Yleaf.main()
    pred = Path(a.out) / "hg_prediction.hg"
    if not pred.is_file() or len(pred.read_text().splitlines()) < 2:
        sys.exit(f"ERROR: Yleaf wrote no prediction ({pred})")


def main(argv):
    if argv[:1] == ["positions"] and len(argv) == 2:
        positions(argv[1])
    elif argv[:1] == ["predict"]:
        predict(argv[1:])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
