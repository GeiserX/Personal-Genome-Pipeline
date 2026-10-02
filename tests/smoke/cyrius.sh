# shellcheck shell=sh
# PYTHON_IMAGE row: step 21 pip-installs Cyrius into this image, with the
# dependency versions held by scripts/cyrius-constraints.txt, and runs it on
# the fixture's Cyrius BAM.
set -e
: "${CYRIUS_VERSION:?versions.env sets it}"
pip install --user --no-cache-dir --disable-pip-version-check -q -c /in/cyrius-constraints.txt "cyrius==${CYRIUS_VERSION}"
PATH="${HOME}/.local/bin:${PATH}"
echo /in/HG002_cyrius.bam > manifest.txt
cyrius --manifest manifest.txt --genome 38 --prefix HG002_cyp2d6 --outDir /out --threads 4
cat HG002_cyp2d6.tsv
