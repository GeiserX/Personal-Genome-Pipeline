# shellcheck shell=sh
# PYTHON_IMAGE row: step 21 pip-installs Cyrius into this image, with the
# dependency versions held by scripts/cyrius-constraints.txt, and runs it on
# the fixture's Cyrius BAM. Cyrius makes no CYP2D6 call on the fixture
# (Genotype None, as in the e2e case), so the row checks that it read depth.
set -e
: "${CYRIUS_VERSION:?versions.env sets it}"
pip install --user --no-cache-dir --disable-pip-version-check -q -c /in/cyrius-constraints.txt "cyrius==${CYRIUS_VERSION}"
PATH="${HOME}/.local/bin:${PATH}"
echo /in/HG002_cyrius.bam > manifest.txt
cyrius --manifest manifest.txt --genome 38 --prefix HG002_cyp2d6 --outDir /out --threads 4
cat HG002_cyp2d6.tsv HG002_cyp2d6.json
