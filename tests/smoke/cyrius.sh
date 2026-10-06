# shellcheck shell=sh
# PYTHON_IMAGE row: install Cyrius as `setup.sh --cyrius` does, every wheel
# checked against its sha256 in scripts/cyrius-constraints.txt (--require-hashes,
# --no-deps, --only-binary), into a directory; then run it from there as step 21
# does, as a module with that directory on PYTHONPATH, on the fixture's Cyrius
# BAM. Cyrius makes no CYP2D6 call on the fixture (Genotype None, as in the
# e2e case), so the row checks that it read depth.
set -e
: "${CYRIUS_VERSION:?versions.env sets it}"
pip install --no-cache-dir --disable-pip-version-check -q --require-hashes --no-deps \
  --only-binary :all: --target /out/cyrius -r /in/cyrius-constraints.txt
PYTHONPATH=/out/cyrius python3 -c 'import importlib.metadata as m; print("cyrius", m.version("cyrius"))' > installed.txt
cat installed.txt
echo /in/HG002_cyrius.bam > manifest.txt
PYTHONPATH=/out/cyrius python3 -m cyrius --manifest manifest.txt --genome 38 --prefix HG002_cyp2d6 --outDir /out --threads 4
cat HG002_cyp2d6.tsv HG002_cyp2d6.json
