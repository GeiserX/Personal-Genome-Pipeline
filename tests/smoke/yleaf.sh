# shellcheck shell=sh
# YLEAF_IMAGE row. The image has no samtools, so steps 37 and Y_HAPLOGROUP
# make the pileup in SAMTOOLS_IMAGE and run Yleaf on it with bin/yleaf_run.py
# (the e2e case sv-mito-telomere-steps-3-y-haplogroup runs that on the
# fixture's chrY slice). This row checks what that launcher relies on in the
# image: the GRCh38 marker table, and the functions and calls it replaces.
set -e
python3 - <<'PY'
import inspect
from yleaf import Yleaf, predict_haplogroup, yleaf_constants as c
p = c.DATA_FOLDER / c.HG38 / c.NEW_POSITION_FILE
n = sum(1 for line in open(p) if line.strip())
print(f"{n} marker positions in {p}")
open("/out/positions_count.txt", "w").write(f"{n}\n")
src = inspect.getsource(Yleaf)
missing = [s for s in ("def call_command(", "samtools idxstats", "samtools mpileup", "-AQ{quality_thresh}q1",
                       "multiprocessing.Pool(", "def check_reference(", "HG38_FULL_GENOME") if s not in src]
missing += [s for s in ("multiprocessing.Pool(",) if s not in inspect.getsource(predict_haplogroup)]
print("missing:", missing)
open("/out/launcher_contract.txt", "w").write("ok\n" if not missing else "missing: " + ", ".join(missing) + "\n")
PY
command -v samtools || echo "no samtools in the image (expected)"
