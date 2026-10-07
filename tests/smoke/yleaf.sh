# shellcheck shell=sh
# YLEAF_IMAGE row. The image installs Yleaf's code only: no samtools and no
# data folder. Steps 37 and Y_HAPLOGROUP make the pileup in SAMTOOLS_IMAGE and
# run Yleaf on it with bin/yleaf_run.py and the data setup.sh --yleaf-data
# installs (the e2e case sv-mito-telomere-steps-3-y-haplogroup runs that on
# the fixture's chrY slice). This row fetches that archive as setup.sh does,
# checks its sha256 (YLEAF_DATA_SHA256), and checks what the launcher relies
# on: the marker table, the tree, and the functions and calls it replaces.
set -e
python3 - <<'PY'
import hashlib, inspect, io, os, tarfile, urllib.request
from yleaf import Yleaf, predict_haplogroup, yleaf_constants as c
v, want = os.environ["YLEAF_DATA_VERSION"], os.environ["YLEAF_DATA_SHA256"]
blob = urllib.request.urlopen(f"https://github.com/genid/Yleaf/archive/refs/tags/{v}.tar.gz", timeout=120).read()
got = hashlib.sha256(blob).hexdigest()
print(f"archive sha256 {got}, versions.env {want}")
tar = tarfile.open(fileobj=io.BytesIO(blob))
pos = tar.extractfile(f"Yleaf-{v}/yleaf/data/hg38/{c.NEW_POSITION_FILE}").read().decode()
tree = tar.extractfile(f"Yleaf-{v}/yleaf/data/hg_prediction_tables/{c.TREE_FILE}").read()
n = sum(1 for line in pos.splitlines() if line.strip())
print(f"{n} marker positions, tree {len(tree)} bytes")
open("/out/sha_ok.txt", "w").write("ok\n" if got == want else f"got {got}\n")
open("/out/positions_count.txt", "w").write(f"{n}\n")
src = inspect.getsource(Yleaf)
missing = [s for s in ("def call_command(", "samtools idxstats", "samtools mpileup", "-AQ{quality_thresh}q1",
                       "multiprocessing.Pool(", "def check_reference(", "HG38_FULL_GENOME", "DATA_FOLDER")
           if s not in src]
missing += [s for s in ("multiprocessing.Pool(", "HG_PREDICTION_FOLDER") if s not in inspect.getsource(predict_haplogroup)]
print("missing:", missing)
open("/out/launcher_contract.txt", "w").write("ok\n" if not missing else "missing: " + ", ".join(missing) + "\n")
PY
