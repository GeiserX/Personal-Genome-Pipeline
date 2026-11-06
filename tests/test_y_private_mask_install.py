#!/usr/bin/env python3
"""setup.sh --yleaf-data keeps step 37's chrY callable mask between runs.

The mask bin/y_private_snps.py keeps to lives in the reference folder, in
<yleaf data>/ytree beside the current tree, checked against
YCALLABLE_MASK_SHA256. This test runs the real scripts/setup.sh and
scripts/lib/common.sh from a copy of the repository whose versions.env pins
synthetic tree, table and mask files, with a fake curl that serves them and
logs every URL it is asked for (no network). Yleaf's own tables are already
in place, so only the tree folder is in play:
  1. a first run downloads the three tree files and the mask;
  2. a second run asks for nothing;
  3. a mask that fails its sha256 (one byte added) is downloaded again, and
     nothing else is;
  4. a new tree pin (an older version.txt) downloads the three tree files and
     keeps the mask, with no request for it;
  5. a tree download that fails keeps the checked mask, and no version.txt,
     so no step reads the tree; the next run asks for the tree files only;
  6. a missing mask beside current tables is downloaded on its own.
Before this test, step 3 asked for nothing and left the bad mask in place,
and steps 4 and 5 downloaded the mask again.

Run: python3 tests/test_y_private_mask_install.py
"""
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAILS = 0

FAKE_CURL = r"""#!/usr/bin/env bash
# Serves $SERVE/<URL's base name> to -o FILE (or stdout) and logs the URL to
# $CURL_LOG. A URL matching $CURL_FAIL, or a file not in $SERVE, exits 22.
out="" url=""
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out=$2; shift 2 ;;
    -C|--connect-timeout|--speed-limit|--speed-time|--retry|--max-time) shift 2 ;;
    -*) shift ;;
    *) url=$1; shift ;;
  esac
done
echo "$url" >> "$CURL_LOG"
if [ -n "${CURL_FAIL:-}" ] && [[ "$url" =~ $CURL_FAIL ]]; then exit 22; fi
src="$SERVE/$(basename "$url")"
[ -f "$src" ] || exit 22
if [ -z "$out" ] || [ "$out" = "-" ]; then cat "$src"; else cat "$src" > "$out"; fi
"""


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {detail}"))
    if not ok:
        FAILS += 1


def sha256(path):
    """The file's sha256, or None when it is not there."""
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return None


def main():
    work = tempfile.mkdtemp(prefix="ymask-")
    try:
        root, serve, gd = (os.path.join(work, x) for x in ("repo", "serve", "genome"))
        os.makedirs(os.path.join(root, "scripts", "lib"))
        os.makedirs(serve)
        shutil.copy(os.path.join(REPO, "scripts", "setup.sh"), os.path.join(root, "scripts"))
        shutil.copy(os.path.join(REPO, "scripts", "lib", "common.sh"), os.path.join(root, "scripts", "lib"))
        with open(os.path.join(REPO, "versions.env")) as f:
            env_text = f.read()
        release = re.search(r'^YCALLABLE_MASK_RELEASE="([^"]+)"', env_text, re.M).group(1)
        yleaf = re.search(r'^YLEAF_DATA_VERSION="([^"]+)"', env_text, re.M).group(1)
        # the synthetic files, under the base names of the URLs setup.sh asks for
        files = {"current_tree.json": b'{"id": "synthetic"}\n',
                 "ftdna_positions_hg38.txt": b"A1\tchrY\t100\tA\tG\n",
                 "new_positions.txt": b"B1\tchrY\t200\tC\tT\n",
                 f"{release}.chrY.mask.fasta.gz": gzip.compress(b">chrY\nNNPPPPNN\n", mtime=0)}
        for name, data in files.items():
            with open(os.path.join(serve, name), "wb") as f:
                f.write(data)
        pins = {"YTREE_SHA256": "current_tree.json", "YTREE_FTDNA_SHA256": "ftdna_positions_hg38.txt",
                "YTREE_YFULL_SHA256": "new_positions.txt", "YCALLABLE_MASK_SHA256": f"{release}.chrY.mask.fasta.gz"}
        for var, name in pins.items():
            env_text, n = re.subn(rf'^{var}="[^"]*"', f'{var}="{sha256(os.path.join(serve, name))}"', env_text,
                                  flags=re.M)
            check(f"versions.env has {var}", n == 1)
        with open(os.path.join(root, "versions.env"), "w") as f:
            f.write(env_text)
        data = os.path.join(gd, "reference", f"yleaf-{yleaf}", "data")
        for rel in ("hg38/new_positions.txt", "hg_prediction_tables/tree.json"):
            os.makedirs(os.path.dirname(os.path.join(data, rel)), exist_ok=True)
            with open(os.path.join(data, rel), "w") as f:
                f.write("seeded\n")
        bindir = os.path.join(work, "bin")
        os.makedirs(bindir)
        with open(os.path.join(bindir, "curl"), "w") as f:
            f.write(FAKE_CURL)
        os.chmod(os.path.join(bindir, "curl"), 0o755)
        log = os.path.join(work, "curl.log")
        yt = os.path.join(data, "ytree")
        mask = os.path.join(yt, "callable_mask_chrY.fasta.gz")
        want_mask = sha256(os.path.join(serve, f"{release}.chrY.mask.fasta.gz"))

        def setup(fail=""):
            open(log, "w").close()
            env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"], SERVE=serve, CURL_LOG=log,
                       CURL_FAIL=fail, CONTAINER_ENGINE="true", FETCH_TRIES="1", FETCH_WAIT="0")
            p = subprocess.run(["bash", os.path.join(root, "scripts", "setup.sh"), "--yleaf-data", gd],
                               capture_output=True, text=True, env=env)
            with open(log) as f:
                asked = [os.path.basename(x.strip()) for x in f if x.strip()]
            return p.returncode, asked, p.stdout + p.stderr

        tree3 = ["current_tree.json", "ftdna_positions_hg38.txt", "new_positions.txt"]
        mask_name = f"{release}.chrY.mask.fasta.gz"

        rc, asked, out = setup()
        check("1. first run: the three tree files and the mask are downloaded",
              rc == 0 and asked == tree3 + [mask_name] and os.path.isfile(os.path.join(yt, "version.txt"))
              and sha256(mask) == want_mask, (rc, asked, out[-800:]))
        rc, asked, out = setup()
        check("2. second run: no request at all", rc == 0 and asked == [], (rc, asked, out[-800:]))
        with open(mask, "ab") as f:
            f.write(b"x")
        version_before = open(os.path.join(yt, "version.txt")).read()
        rc, asked, out = setup()
        check("3. a mask that fails its sha256 is downloaded again, and nothing else",
              rc == 0 and asked == [mask_name] and sha256(mask) == want_mask
              and open(os.path.join(yt, "version.txt")).read() == version_before, (rc, asked, out[-800:]))
        with open(os.path.join(yt, "version.txt"), "w") as f:
            f.write("tree\tan older pin\n")
        rc, asked, out = setup()
        check("4. a new tree pin downloads the tree files and keeps the mask, with no request for it",
              rc == 0 and asked == tree3 and sha256(mask) == want_mask
              and open(os.path.join(yt, "version.txt")).read() == version_before, (rc, asked, out[-800:]))
        with open(os.path.join(yt, "version.txt"), "w") as f:
            f.write("tree\tan older pin\n")
        rc, asked, out = setup(fail="current_tree")
        check("5. a failed tree download keeps the checked mask and leaves no version.txt",
              rc == 0 and mask_name not in asked and sha256(mask) == want_mask
              and not os.path.exists(os.path.join(yt, "version.txt"))
              and sorted(os.listdir(yt)) == ["callable_mask_chrY.fasta.gz"], (rc, asked, out[-800:]))
        rc, asked, out = setup()
        check("5. the next run asks for the tree files only", rc == 0 and asked == tree3
              and sha256(mask) == want_mask, (rc, asked, out[-800:]))
        os.remove(mask)
        rc, asked, out = setup()
        check("6. a missing mask beside current tables is downloaded on its own",
              rc == 0 and asked == [mask_name] and sha256(mask) == want_mask, (rc, asked, out[-800:]))
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
