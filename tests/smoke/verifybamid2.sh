# shellcheck shell=sh
# VERIFYBAMID2_IMAGE row: step 33's VerifyBamID2 call on the fixture slice,
# with the markers of the 1000 Genomes 100k GRCh38 panel that fall inside the
# fixture's regions (verifybamid2_panel.*, 320 of 100,000). VerifyBamID2
# refuses fewer than 1,000 markers with reads ("Insufficient Available
# markers"); step 33 then runs it again with --DisableSanityCheck, and so
# does this row: the first call must fail that way, the second must estimate.
set -e
if verifybamid2 --SVDPrefix /smoke/verifybamid2_panel --Reference /in/ref.fa --BamFile /in/HG002_slice.bam \
     --Output checked --NumThread 4 > checked.log 2>&1; then
  echo "the sanity check passed on 320 markers" > checked.result
else
  grep -h 'Insufficient Available markers' checked.log > checked.result || { cat checked.log; exit 1; }
fi
verifybamid2 --SVDPrefix /smoke/verifybamid2_panel --Reference /in/ref.fa --BamFile /in/HG002_slice.bam \
  --Output vb --NumThread 4 --DisableSanityCheck > vb.log 2>&1 || { cat vb.log; exit 1; }
tail -n 5 vb.log
cat checked.result vb.selfSM
