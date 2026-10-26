# shellcheck shell=sh
# VNTYPER_IMAGE row: what step 52 runs, on VNtyper's own public test data
# (tests/smoke/vntyper_testdata.py fetches two BAMs out of its Zenodo archive):
# a sample with a MUC1 VNTR frameshift and one without. The image's samtools
# indexes them, then VNtyper runs with the step's flags from /opt/vntyper,
# the folder its config.json paths are relative to.
set -e
python3 /smoke/vntyper_testdata.py /out/data positive.bam negative.bam
for s in positive negative; do
  samtools index "/out/data/${s}.bam"
  (cd /opt/vntyper && vntyper pipeline --bam "/out/data/${s}.bam" -o "/out/${s}" --sample-name "$s" \
     --threads 2 --fast-mode --reference-assembly hg38 --report-igv off) > "/out/${s}.log" 2>&1 \
    || { tail -n 40 "/out/${s}.log"; exit 1; }
  echo "== ${s}: kestrel/kestrel_result.tsv"; cat "/out/${s}/kestrel/kestrel_result.tsv"
  echo "== ${s}: coverage/coverage_summary.tsv"; cat "/out/${s}/coverage/coverage_summary.tsv"
done
# A BAM with the same header and no read over MUC1, as a genome slice without
# chr1:155 Mb gives: what VNtyper writes then, and its exit code.
samtools view -b -o /out/data/empty.bam /out/data/negative.bam chr1:1-1000 && samtools index /out/data/empty.bam
rc=0
(cd /opt/vntyper && vntyper pipeline --bam /out/data/empty.bam -o /out/empty --sample-name empty \
   --threads 2 --fast-mode --reference-assembly hg38 --report-igv off) > /out/empty.log 2>&1 || rc=$?
echo "$rc" > empty_rc.txt
echo "== empty: exit ${rc}"; tail -n 30 /out/empty.log
ls -R /out/empty 2>/dev/null | head -40
cat /out/empty/kestrel/kestrel_result.tsv /out/empty/coverage/coverage_summary.tsv 2>/dev/null || true
vntyper --version > version.txt 2>&1
ls -R /out/positive | head -60
python3 -c "import json; d = json.load(open('/out/positive/pipeline_summary.json')); print(sorted(d)); print(json.dumps(d, indent=1)[:4000])"
# calls.tsv: sample, Kestrel's confidence, the variant's name, the mean depth
# over the VNTR and VNtyper's coverage verdict.
python3 - <<'PY'
import csv
rows = []
for s in ("positive", "negative"):
    k = [l for l in open(f"/out/{s}/kestrel/kestrel_result.tsv") if not l.startswith("#")]
    r = next(csv.DictReader(k, delimiter="\t"))
    c = next(csv.DictReader(open(f"/out/{s}/coverage/coverage_summary.tsv"), delimiter="\t"))
    rows.append("\t".join([s, r.get("Confidence", ""), r.get("Nomenclature", "") or r.get("Variant", ""),
                           c.get("mean", ""), c.get("coverage_qc", "")]))
open("/out/calls.tsv", "w").write("".join(x + "\n" for x in rows))
print("\n".join(rows))
PY
