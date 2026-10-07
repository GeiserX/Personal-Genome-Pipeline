# shellcheck shell=sh
# YLEAF_IMAGE row. Yleaf downloads the whole hg38 FASTA on its first run
# unless its config names one, and the image's config is read-only, so step
# 37 points Yleaf's constant at the reference before it starts (the same
# python lines as scripts/37-y-haplogroup.sh). A BAM never needs the FASTA's
# sequence: Yleaf reads its pileup at the Y markers.
set -e
python3 -c 'import sys; from pathlib import Path; from yleaf import yleaf_constants; yleaf_constants.HG38_FULL_GENOME = Path(sys.argv[1]); from yleaf import Yleaf; sys.argv = ["Yleaf"] + sys.argv[2:]; Yleaf.main()' /in/ref.fa -bam /in/HG002_slice.bam -o /out/yleaf -rg hg38 -force -t 2
cat /out/yleaf/hg_prediction.hg
