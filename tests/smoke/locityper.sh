# shellcheck shell=sh
# shellcheck disable=SC2154  # LOCITYPER_PANEL_URL and _TBI_SHA256 come from versions.env, set by image-smoke.sh
# LOCITYPER_IMAGE row: what setup.sh --locityper-data and step 53 run, on the
# fixture's HG002 reads, checked against HG002's own assembly haplotypes.
#
# The HPRC v2.1 panel (LOCITYPER_PANEL_URL) holds HG002's two assembled
# haplotypes, HG002.1 and HG002.2: they are the truth. The row:
#   1. builds a reference of chr6, chr20 and chr22 from the fixture's, and
#      realigns the fixture's reads over HLA class II, CYP2D6 and chr20:10.0-
#      10.5 Mb to it (Locityper needs the BAM's contigs to be the reference's);
#   2. as setup.sh does: jellyfish k-mer counts, the panel's .tbi checked by
#      LOCITYPER_PANEL_TBI_SHA256, the panel read over https 250 kb either
#      side of each locus, locityper target and augment;
#   3. as step 53 does: preproc on chr20:10.0-10.5 Mb (the fixture has no
#      chr17, where Locityper's own background region is), then genotype;
#   4. aligns each called haplotype to HG002.1 and HG002.2 with locityper
#      align and keeps the better pairing: qv.tsv;
#   5. the red control: every haplotype within QV 20 of either truth
#      haplotype taken out of the panel (HG002's own and the near copies other
#      people carry), the database rebuilt from that panel, genotype again:
#      red_qv.tsv, whose calls cannot reach QV 20.
# The loci are those of assets/locityper_loci.bed (HLA-DRB1, HLA-DQB1,
# CYP2D6); the fixture has no read over FCGR or KIR.
set -e
T=$(nproc)
cd /out

echo "== 1. reference and reads"
samtools faidx /in/ref.fa chr6 chr20 chr22 > mini.fa && samtools faidx mini.fa
samtools view -u -F 0x900 /in/HG002_slice.bam chr6:31500001-33100000 chr20:10000001-10500000 chr22:42000001-42300000 \
  | samtools collate -O -u - collate | samtools fastq -n -1 r1.fq.gz -2 r2.fq.gz -0 /dev/null -s /dev/null -
minimap2 -t "$T" -a -x sr -R '@RG\tID:HG002\tSM:HG002' mini.fa r1.fq.gz r2.fq.gz 2> minimap2.log \
  | samtools sort -@ 2 -o hg.bam - && samtools index hg.bam
rm -f r1.fq.gz r2.fq.gz

echo "== 2. database"
printf 'chr6\t32484347\t32593538\tHLA-DRB1\nchr6\t32652603\t32669077\tHLA-DQB1\nchr22\t42122284\t42152960\tCYP2D6\n' > loci.bed
jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 --threads "$T" --size 400M -o mini.jf mini.fa
wget -q -O panel.tbi "${LOCITYPER_PANEL_URL}.tbi"
sha256sum panel.tbi | cut -d' ' -f1 > panel_tbi_sha256.txt
[ "$(cat panel_tbi_sha256.txt)" = "$LOCITYPER_PANEL_TBI_SHA256" ] || { echo "panel .tbi sha256 $(cat panel_tbi_sha256.txt), want ${LOCITYPER_PANEL_TBI_SHA256}"; exit 1; }
tabix -h "${LOCITYPER_PANEL_URL}##idx##/out/panel.tbi" chr6:32234347-32919077 chr22:41872284-42402960 | bgzip > panel.vcf.gz
tabix -p vcf panel.vcf.gz
locityper target -d db -r mini.fa -j mini.jf -v panel.vcf.gz -g GRCh38 -L loci.bed -@ "$T" > target.log 2>&1 || { tail -n 30 target.log; exit 1; }
locityper augment -d db -@ "$T" > augment.log 2>&1 || { tail -n 30 augment.log; exit 1; }

echo "== 3. preproc and genotype"
locityper preproc -a hg.bam -r mini.fa -j mini.jf -o bg -b chr20:10000001-10500000 -@ "$T" > preproc.log 2>&1 || { tail -n 30 preproc.log; exit 1; }
locityper genotype -a hg.bam -d db -p bg -o gt -@ "$T" > genotype.log 2>&1 || { tail -n 30 genotype.log; exit 1; }
locityper --version > version.txt 2>&1 || true

# qv.py GT_DIR OUT: per locus the call and the QV of each called haplotype to
# the truth haplotype it pairs with (the better of the two pairings); 99 for an
# identical sequence. The sequences come from the full panel (db).
cat > qv.py <<'PY'
import gzip, json, os, subprocess, sys
gt_dir, out = sys.argv[1], sys.argv[2]
rows = []
for locus in ("HLA-DRB1", "HLA-DQB1", "CYP2D6"):
    res = json.load(gzip.open(f"{gt_dir}/loci/{locus}/res.json.gz", "rt"))
    a, b = res["genotype"].split(",")
    fa = f"db/loci/{locus}/haplotypes.fa.gz"
    pairs = sorted({(x, y) for x in (a, b) for y in ("HG002.1", "HG002.2") if x != y})
    q = {(x, y): 99.0 for x in (a, b) for y in ("HG002.1", "HG002.2") if x == y}
    if pairs:
        paf = f"{gt_dir}.{locus}.paf"
        subprocess.run(["locityper", "align", "-i", fa, "-o", paf, "-@", "2", "-p"] + [f"{x},{y}" for x, y in pairs],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for line in open(paf):
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            tags = dict(t.split(":", 2)[::2] for t in c[12:])
            v = 99.0 if tags.get("NM") == "0" else min(99.0, float(tags["qv"]))
            q[(c[0], c[5])] = q[(c[5], c[0])] = v
    p1 = (q.get((a, "HG002.1"), 0.0), q.get((b, "HG002.2"), 0.0))
    p2 = (q.get((a, "HG002.2"), 0.0), q.get((b, "HG002.1"), 0.0))
    best = max(p1, p2, key=min)
    rows.append([locus, res["genotype"], f"{best[0]:.1f}", f"{best[1]:.1f}", f"{min(best):.1f}",
                 f"{res['quality']:.1f}", str(res.get("weight_dist", "nan"))])
with open(out, "w") as f:
    for r in rows:
        f.write("\t".join(r) + "\n")
        print("\t".join(r))
PY
echo "== 4. the calls against HG002's haplotypes (locus, genotype, QV 1, QV 2, lower, quality, weighted distance)"
python3 qv.py gt qv.tsv

echo "== 5. red control: the panel without HG002 and every haplotype within QV 20 of it"
: > red.bed
for l in HLA-DRB1 HLA-DQB1 CYP2D6; do
  locityper align -i "db/loci/${l}/haplotypes.fa.gz" -o "near.${l}.paf" --against HG002.1 HG002.2 --against-div 0.05 -@ "$T" > /dev/null 2>&1
  python3 - "$l" <<'PY'
import gzip, sys
l = sys.argv[1]
near = {"HG002.1", "HG002.2"}
for line in open(f"near.{l}.paf"):
    if line.startswith("#"):
        continue
    c = line.rstrip("\n").split("\t")
    tags = dict(t.split(":", 2)[::2] for t in c[12:])
    if tags.get("NM") == "0" or float(tags.get("qv", "0")) >= 20:
        near.update((c[0], c[5]))
keep, name, n = [], None, 0
for line in gzip.open(f"db/loci/{l}/haplotypes.fa.gz", "rt"):
    if line.startswith(">"):
        name = line[1:].split()[0]
        n += 1
    if name not in near:
        keep.append(line)
open(f"red.{l}.fa", "w").write("".join(keep))
kept = sum(x.startswith(">") for x in keep)
open("red_removed.tsv", "a").write(f"{l}\t{n}\t{n - kept}\n")
print(l, n, "haplotypes,", n - kept, "within QV 20 of HG002 removed")
PY
  awk -v l="$l" -v f="/out/red.${l}.fa" -F'\t' '$4 == l {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" f}' loci.bed >> red.bed
done
grep -h '^>' red.*.fa | grep -c '^>HG002\.' > red_hg002.txt || true
locityper target -d dbred -r mini.fa -j mini.jf -L red.bed -@ "$T" > target_red.log 2>&1 || { tail -n 30 target_red.log; exit 1; }
locityper genotype -a hg.bam -d dbred -p bg -o gtred -@ "$T" > genotype_red.log 2>&1 || { tail -n 30 genotype_red.log; exit 1; }
python3 qv.py gtred red_qv.tsv
rm -rf collate* hg.bam* mini.jf
