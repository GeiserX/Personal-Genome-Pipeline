#!/usr/bin/env python3
"""yleaf_run.py: run Yleaf 4.1.6 (YLEAF_IMAGE) offline, on a pileup made beside it.

Step 37 and the Nextflow processes make one samtools pass per BAM in
SAMTOOLS_IMAGE: Yleaf's pileup, the current tree's pileup (y_tree_place.py)
and the chrY BAM. Yleaf then reads its pileup from a file instead of running
samtools itself, with its data folder (marker positions, haplogroup trees)
set to the copy `setup.sh --yleaf-data` installs from the release archive as
reference/yleaf-<version>/data, the folder that also holds the current tree.
Yleaf also downloads the whole hg38 FASTA when its read-only config names
none, and checks GitHub for a newer release. So it runs in three parts:

  yleaf_run.py positions --data DATA OUT
      (YLEAF_IMAGE) write Yleaf's GRCh38 marker positions as the "chrY<TAB>pos"
      list `samtools mpileup -l` reads
  samtools idxstats BAM > IDXSTATS; samtools mpileup -l OUT -ABQ20q1 BAM > PILEUP
      (SAMTOOLS_IMAGE) the two commands Yleaf would run, with its default
      quality threshold of 20
  yleaf_run.py predict --data DATA --bam BAM --reference FASTA --idxstats IDXSTATS --pileup PILEUP --out DIR
      (YLEAF_IMAGE) Yleaf's own flow on its default tree (YFull v14.01),
      with its data folder set to DATA, its reference constant set to FASTA
      (a BAM never needs the sequence), no update check, each samtools call
      served from those two files (a pileup made with other flags than the
      ones Yleaf asks for stops the run), and its multiprocessing pools
      replaced by a serial map, so every call runs in this process, where
      the served samtools calls are.

DIR/hg_prediction.hg is Yleaf's prediction table (Hg is NA when too few
markers had reads). Standard library plus what Yleaf itself imports.
"""
import multiprocessing
import shutil
import sys
from pathlib import Path

QUALITY = 20   # Yleaf's -q default; the mpileup above uses the same
FLAGS = f"-ABQ{QUALITY}q1"   # what Yleaf 4 gives samtools mpileup at that quality


class SerialPool:
    """multiprocessing.Pool in this process: the initializer runs here, once."""

    def __init__(self, processes=None, initializer=None, initargs=(), **kwargs):
        if initializer is not None:
            initializer(*initargs)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def map(self, fn, items):
        return list(map(fn, items))

    def imap(self, fn, items):
        return map(fn, items)

    imap_unordered = imap


def use_data(data):
    """Point Yleaf's constants at DATA (the yleaf/data folder of its release)."""
    from yleaf import yleaf_constants as c
    # Absolute: Yleaf opens DATA_FOLDER / HG_PREDICTION_FOLDER / tree.json.
    # An absolute HG_PREDICTION_FOLDER replaces DATA_FOLDER in that join; a
    # relative one gives data/data/... (a Nextflow task stages it as "data").
    data = Path(data).absolute()
    for need in (data / c.HG38 / c.NEW_POSITION_FILE, data / "hg_prediction_tables" / c.TREE_FILE):
        if not need.is_file():
            sys.exit(f"ERROR: {need} not found: install Yleaf's data with scripts/setup.sh --yleaf-data <genome_dir>")
    c.DATA_FOLDER = data
    c.HG_PREDICTION_FOLDER = data / "hg_prediction_tables"
    return c


def positions(data, out):
    c = use_data(data)
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
    for a in ("--data", "--bam", "--reference", "--idxstats", "--pileup", "--out"):
        ap.add_argument(a, required=True)
    a = ap.parse_args(argv)
    multiprocessing.Pool = SerialPool
    yleaf_constants = use_data(a.data)
    yleaf_constants.HG38_FULL_GENOME = Path(a.reference)
    if Path(a.pileup).stat().st_size == 0:
        # No read at any marker: Yleaf's pandas read of an empty pileup fails
        # before it predicts, so write the prediction it gives for too few
        # markers (Hg NA) here, with the BAM's mapped reads from idxstats.
        reads = 0
        with open(a.idxstats) as f:
            for line in f:
                cols = line.split("\t")
                if len(cols) > 2 and cols[2].strip().isdigit():
                    reads += int(cols[2])
        Path(a.out).mkdir(parents=True, exist_ok=True)
        (Path(a.out) / "hg_prediction.hg").write_text(
            "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n"
            f"{Path(a.bam).name.rsplit('.', 1)[0]}\tNA\t\t{reads}\t0\tNA\tNA\tNA\tNA\n")
        print("No read at any Y marker: Hg NA (insufficient markers)")
        return
    from yleaf import Yleaf

    def serve(cmd, stdout_location=None):
        if cmd.startswith("samtools idxstats "):
            with open(a.idxstats) as f:
                shutil.copyfileobj(f, stdout_location)
            stdout_location.flush()
        elif cmd.startswith("samtools mpileup ") and " > " in cmd:
            if f" {FLAGS} " not in f" {cmd} ":
                sys.exit(f"ERROR: Yleaf asks for samtools mpileup with other flags than {FLAGS}, the ones the "
                         f"pileup was made with: {cmd}")
            shutil.copyfile(a.pileup, cmd.rsplit(" > ", 1)[1].strip().strip('"'))
        elif cmd.startswith("samtools index "):
            pass   # the BAM is indexed; Yleaf only asks when it finds no .bai
        else:
            sys.exit(f"ERROR: Yleaf asked for a command this launcher does not serve: {cmd}")

    Yleaf.call_command = serve
    sys.argv = ["Yleaf", "-bam", a.bam, "-o", a.out, "-rg", "hg38", "-tree", "yfull", "-force", "-t", "1",
                "-q", str(QUALITY), "--no-update-check"]
    Yleaf.main()
    pred = Path(a.out) / "hg_prediction.hg"
    if not pred.is_file() or len(pred.read_text().splitlines()) < 2:
        # Yleaf 4 logs a failed sample as "Skipping" and goes on: its log says why.
        sys.exit(f"ERROR: Yleaf wrote no prediction ({pred}); its log is in {a.out}")


def main(argv):
    if argv[:2] == ["positions", "--data"] and len(argv) == 4:
        positions(argv[2], argv[3])
    elif argv[:1] == ["predict"]:
        predict(argv[1:])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
