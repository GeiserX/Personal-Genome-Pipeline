# shellcheck shell=sh
# YLEAF_IMAGE row. Yleaf downloads the whole hg38 FASTA on its first run
# unless its config names one, and the image's config is read-only, so step
# 37 points Yleaf's constant at the reference before it starts (the same
# python line as scripts/37-y-haplogroup.sh, which also replaces Yleaf's
# multiprocessing pools by a serial map so a failure stops it instead of
# hanging it). A BAM never needs the FASTA's sequence.
set -e
python3 -c 'import sys, multiprocessing; from pathlib import Path; multiprocessing.Pool = type("SerialPool", (), {"__init__": lambda s, *a, **k: None, "__enter__": lambda s: s, "__exit__": lambda s, *a: False, "map": lambda s, f, xs: list(map(f, xs))}); from yleaf import yleaf_constants; yleaf_constants.HG38_FULL_GENOME = Path(sys.argv[1]); from yleaf import Yleaf; sys.argv = ["Yleaf"] + sys.argv[2:]; Yleaf.main()' /in/ref.fa -bam /in/HG002_slice.bam -o /out/yleaf -rg hg38 -force -t 1
cat /out/yleaf/hg_prediction.hg
