# shellcheck shell=sh
# VERIFYBAMID2_IMAGE row: step 33's VerifyBamID2 call on the fixture slice,
# with the markers of the 1000 Genomes 100k GRCh38 panel that fall inside the
# fixture's regions (verifybamid2_panel.*, 320 of 100,000).
set -e
verifybamid2 --SVDPrefix /smoke/verifybamid2_panel --Reference /in/ref.fa --BamFile /in/HG002_slice.bam \
  --Output vb --NumThread 4 > vb.log 2>&1 || { cat vb.log; exit 1; }
cat vb.log vb.selfSM
