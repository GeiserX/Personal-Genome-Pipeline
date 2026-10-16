#!/usr/bin/env bash
# Step 40 on the fixture: the 23andMe-layout export of HG002's gVCF (case 21
# called it), with the CI slice of the site table
# (tests/fixtures/consumer_export/sites.tsv.gz: the 23andMe v5 sites dbSNP
# places inside the fixture's regions).
#   1. step 40 writes both zips; the strict reader passes both;
#   2. every site gets the genotype a second, independent reading of the
#      gVCF gives it (each site looked up among all gVCF records), almost
#      every site is called, and almost none lacks a gVCF record;
#   3. a planted hole (the gVCF record over one called site removed) makes
#      that site '--' and changes no other site;
#   4. round trip: scripts/chip-to-vcf.sh reads the exported 23andMe file and
#      the AncestryDNA file back (bcftools tsv2vcf on GRCh37, UCSC's
#      hg19ToHg38 chain, Picard LiftoverVcf): each lifted record lands on the
#      table's GRCh38 position and has the alleles of case 21's calls there
#      (hom-ref where case 21 called nothing), but for at most 0.1% of the
#      sites, where the list's GRCh37 position and dbSNP's GRCh38 position
#      fall on two copies of a duplicated region. chrM is left out: the hg19
#      chain maps hg19's chrM, not the rCRS positions the file uses;
#   5. red controls: the table shifted by one GRCh38 position stops the
#      export (reference bases differ); shifted by one GRCh37 position, the
#      export runs and the round trip fails;
#   6. the Nextflow module on the same gVCF writes the same files.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
D="${G}/${SAMPLE}"
GVCF="${D}/vcf/${SAMPLE}.g.vcf.gz"
VCF="${D}/vcf/${SAMPLE}.vcf.gz"
check "case 21 wrote a gVCF" vcf_ok "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
SLICE="${REPO}/tests/fixtures/consumer_export/sites.tsv.gz"
SITES="${G}/reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
mkdir -p "$(dirname "$SITES")"
cp "$SLICE" "$SITES"
TOOL="${REPO}/bin/export_consumer.py"
CHAIN_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
CHAIN_SHA256=5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d

# unz ZIP OUT: the one member of ZIP to OUT
unz() { python3 -c 'import sys, zipfile; z = zipfile.ZipFile(sys.argv[1]); open(sys.argv[2], "wb").write(z.read(z.namelist()[0]))' "$1" "$2"; }

# --- 1. step 40 ----------------------------------------------------------------
run_step 40-export-consumer.sh "$SAMPLE" male
check_step_exit 40-export-consumer.sh
OUT="${D}/consumer_export"
Z23="${OUT}/${SAMPLE}_23andme_v5.zip"
ZAN="${OUT}/${SAMPLE}_ancestrydna_layout.zip"
SUM="${OUT}/${SAMPLE}_consumer_export.tsv"
check "the 23andMe zip is written" test -s "$Z23"
check "the AncestryDNA zip is written" test -s "$ZAN"
check "strict reader: the 23andMe layout" python3 "$TOOL" check --layout 23andme "$Z23"
check "strict reader: the AncestryDNA layout" python3 "$TOOL" check --layout ancestrydna "$ZAN"
unz "$Z23" "${CASE_TMP}/export_23andme.txt"
unz "$ZAN" "${CASE_TMP}/export_ancestry.txt"
cat "$SUM"
val() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "$2"; }
N_SITES=$(gzip -dc "$SLICE" | grep -vc '^#')
check_eq "every site of the table is in the file" "$(grep -vc '^#' "${CASE_TMP}/export_23andme.txt")" "$N_SITES"

# --- 2. every site against an independent reading of the gVCF --------------------
# expected_calls GVCF SITES: rsid and the letters the gVCF gives each site
# (male rules), looked up among all records with a bisect, not the step's walk:
# '--' inside a non-reference deletion or longer allele; else the variant
# record that starts at the site; else the reference block over it.
expected_calls() {
  python3 - "$1" "$2" "${REPO}/assets/par_grch38.bed" <<'PY'
import bisect, gzip, re, sys
gvcf, sites, par_bed = sys.argv[1:4]
par = [(l.split()[0], int(l.split()[1]) + 1, int(l.split()[2])) for l in open(par_bed) if l.strip()]
recs = {}
for line in gzip.open(gvcf, "rt"):
    if line.startswith("#"):
        continue
    c = line.rstrip("\n").split("\t")
    pos, ref, alts = int(c[1]), c[3].upper(), c[4].upper().split(",")
    m = re.search(r"(?:^|;)END=(\d+)", c[7])
    end = int(m.group(1)) if m else pos + len(ref) - 1
    keys, vals = c[8].split(":"), c[9].split(":")
    recs.setdefault(c[0], []).append((pos, end, ref, alts, c[6], vals[keys.index("GT")]))
starts = {k: [r[0] for r in v] for k, v in recs.items()}
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
symbolic = lambda a: a in ("<*>", "<NON_REF>", ".")
for line in gzip.open(sites, "rt"):
    if line.startswith("#"):
        continue
    rsid, chrom, pos37, ref37, pos38, ref38, strand = line.rstrip("\n").split("\t")
    p, contig = int(pos38), ("chrM" if chrom == "MT" else "chr" + chrom)
    lst = recs.get(contig, [])
    j = bisect.bisect_right(starts.get(contig, []), p)
    cover = [r for r in lst[max(0, j - 50):j] if r[0] <= p <= r[1]]
    nonref = lambda r: any(i not in ("0", ".", "") for i in re.split(r"[/|]", r[5]))
    inside = [r for r in cover if r[0] < p and not all(symbolic(a) for a in r[3]) and len(r[2]) > 1 and nonref(r)]
    at = [r for r in cover if r[0] == p and not all(symbolic(a) for a in r[3])]
    blk = [r for r in cover if all(symbolic(a) for a in r[3])]
    letters = None
    if inside:
        pass
    elif at:
        r = at[0]
        idx = re.split(r"[/|]", r[5])
        if "." not in idx and not (r[4] not in ("PASS", ".") and any(i != "0" for i in idx)) and r[2][0] == ref38:
            bases = [([r[2]] + r[3])[int(i)][0] for i in idx]
            letters = bases if all(b in "ACGT" for b in bases) else None
    elif blk:
        idx = re.split(r"[/|]", blk[-1][5])
        letters = [ref38] * len(idx) if all(i == "0" for i in idx) else None
    if letters is not None:
        haploid = chrom == "MT" or (chrom in ("X", "Y") and not any(c == contig and s <= p <= e for c, s, e in par))
        if haploid and len(letters) == 2:
            letters = letters[:1] if letters[0] == letters[1] else None
    if letters is not None and strand == "-":
        letters = [comp[b] for b in letters]
    print(f"{rsid}\t{''.join(sorted(letters)) if letters else '--'}")
PY
}
expected_calls "$GVCF" "$SITES" > "${CASE_TMP}/expected.tsv"
grep -v '^#' "${CASE_TMP}/export_23andme.txt" | cut -f1,4 > "${CASE_TMP}/exported.tsv"
DIFF=$(join -t $'\t' <(LC_ALL=C sort "${CASE_TMP}/expected.tsv") <(LC_ALL=C sort "${CASE_TMP}/exported.tsv") | awk -F'\t' '$2 != $3' | wc -l | tr -d ' ')
head -n 5 < <(join -t $'\t' <(LC_ALL=C sort "${CASE_TMP}/expected.tsv") <(LC_ALL=C sort "${CASE_TMP}/exported.tsv") | awk -F'\t' '$2 != $3')
check_eq "sites whose exported genotype differs from the independent reading of the gVCF" "$DIFF" 0
CALLED=$(val called "$SUM")
check_ge "sites called, in percent of ${N_SITES}" "$(( ${CALLED:-0} * 100 / N_SITES ))" 90
NOREC=$(val 'no_call: no gVCF record' "$SUM")
check "sites with no gVCF record (${NOREC:-0}) are under 1% of ${N_SITES}" test "$(( ${NOREC:-0} * 100 ))" -lt "$N_SITES"
check_ge "sites homozygous reference from the gVCF's reference blocks" "$(val from_reference_blocks "$SUM")" 1000
check_eq "reference bases that differ from the table" "$(val reference_mismatch "$SUM")" 0
check_ge "reference bases checked against the table" "$(val reference_checked "$SUM")" 100

# --- 3. a planted hole -----------------------------------------------------------
H=HG002CXH
rm -rf "${G:?}/${H}"
mkdir -p "${G}/${H}/vcf"
HOLE=$(python3 - "$GVCF" "${G}/${H}/vcf/${H}.g.vcf.gz" "${CASE_TMP}/exported.tsv" "$SITES" <<'PY'
import gzip, re, sys
src, dst, exported, sites = sys.argv[1:5]
called = {l.split("\t")[0] for l in open(exported) if l.rstrip("\n").split("\t")[1] != "--"}
site = None
for l in gzip.open(sites, "rt"):
    c = l.rstrip("\n").split("\t")
    if not l.startswith("#") and c[1] == "20" and c[0] in called:
        site = c
        break
p = int(site[4])
dropped = 0
with gzip.open(dst, "wt") as out:
    for line in gzip.open(src, "rt"):
        c = line.split("\t", 8)
        if not line.startswith("#") and c[0] == "chr20":
            m = re.search(r"(?:^|;)END=(\d+)", c[7])
            end = int(m.group(1)) if m else int(c[1]) + len(c[3]) - 1
            if int(c[1]) <= p <= end:
                dropped += 1
                continue
        out.write(line)
print(f"{site[0]} {dropped}")
PY
)
read -r HOLE_ID HOLE_RECS <<< "$HOLE"
echo "planted hole: ${HOLE_ID} (${HOLE_RECS} gVCF record(s) removed)"
run_step 40-export-consumer.sh "$H" male
check_step_exit "40-export-consumer.sh (planted hole)"
unz "${G}/${H}/consumer_export/${H}_23andme_v5.zip" "${CASE_TMP}/hole.txt"
check_eq "the planted hole's site ${HOLE_ID:-none} is --" "$(awk -F'\t' -v r="${HOLE_ID:-none}" '$1 == r {print $4}' "${CASE_TMP}/hole.txt")" "--"
# Only sites the removed record covered may change.
expected_calls "${G}/${H}/vcf/${H}.g.vcf.gz" "$SITES" > "${CASE_TMP}/hole_expected.tsv"
HDIFF=$(join -t $'\t' <(LC_ALL=C sort "${CASE_TMP}/hole_expected.tsv") <(grep -v '^#' "${CASE_TMP}/hole.txt" | cut -f1,4 | LC_ALL=C sort) \
  | awk -F'\t' '$2 != $3' | wc -l | tr -d ' ')
check_eq "planted hole: sites that differ from the independent reading of the holed gVCF" "$HDIFF" 0
CHANGED=$(join -t $'\t' <(LC_ALL=C sort "${CASE_TMP}/exported.tsv") <(grep -v '^#' "${CASE_TMP}/hole.txt" | cut -f1,4 | LC_ALL=C sort) \
  | awk -F'\t' '$2 != $3' | wc -l | tr -d ' ')
check "planted hole: 1 to 50 sites changed, all of them to -- (${CHANGED})" test "$CHANGED" -ge 1 -a "$CHANGED" -le 50
rm -rf "${G:?}/${H}"

# --- 4. round trip through chip-to-vcf.sh ------------------------------------------
# chip-to-vcf.sh reads GRCh37: a FASTA with the table's GRCh37 base at each
# site (N elsewhere; those bases agree with dbSNP's GRCh37 records, see
# tests/fixtures/consumer_export/README.md) and UCSC's chain.
mkdir -p "${G}/reference_hg19" "${G}/liftover"
python3 - "$SITES" "${G}/reference_hg19/human_g1k_v37.fasta" <<'PY'
import gzip, sys
sites, fa = sys.argv[1:3]
order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
bases = {}
for l in gzip.open(sites, "rt"):
    if l.startswith("#"):
        continue
    c = l.rstrip("\n").split("\t")
    bases.setdefault(c[1], {})[int(c[2])] = c[3]
off = 0
with open(fa, "w") as f, open(fa + ".fai", "w") as fai:
    for chrom in [c for c in order if c in bases]:
        n = max(bases[chrom]) + 1
        head = f">{chrom}\n"
        f.write(head)
        off += len(head)
        fai.write(f"{chrom}\t{n}\t{off}\t60\t61\n")
        seq = bytearray(b"N" * n)
        for p, b in bases[chrom].items():
            seq[p - 1] = ord(b)
        for i in range(0, n, 60):
            f.write(seq[i:i + 60].decode() + "\n")
        off += n + (n + 59) // 60
PY
if curl -fsSL --retry 3 -o "${G}/liftover/hg19ToHg38.over.chain.gz" "$CHAIN_URL"; then
  check_eq "UCSC's hg19ToHg38 chain has the pinned sha256" \
    "$(sha256sum "${G}/liftover/hg19ToHg38.over.chain.gz" | cut -d' ' -f1)" "$CHAIN_SHA256"
else
  fail "could not download ${CHAIN_URL}"
fi

# roundtrip NAME TXT FORMAT: chip-to-vcf.sh on TXT as sample NAME; prints
# "lifted moved compared nonref mismatches" against case 21's calls.
roundtrip() {
  local n=$1
  rm -rf "${G:?}/${n}"
  mkdir -p "${G}/${n}/raw"
  cp "$2" "${G}/${n}/raw/${n}_raw.txt"
  run_step chip-to-vcf.sh "$n" "$3" >&2
  bcf query -f '%ID\t%CHROM\t%POS[\t%GT]\t%REF\t%ALT\n' "${n}/vcf/${n}.vcf.gz" > "${CASE_TMP}/${n}_lifted.tsv" 2>/dev/null
  bcf query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT]\n' "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" > "${CASE_TMP}/calls.tsv" 2>/dev/null
  python3 - "${CASE_TMP}/${n}_lifted.tsv" "${CASE_TMP}/calls.tsv" "$SITES" "$2" <<'PY'
import gzip, re, sys
lifted, calls, sites, txt = sys.argv[1:5]
table = {}
for l in gzip.open(sites, "rt"):
    if not l.startswith("#"):
        c = l.rstrip("\n").split("\t")
        table[c[0]] = c
dv = {}
for l in open(calls):
    chrom, pos, ref, alt, flt, gt = l.rstrip("\n").split("\t")
    if flt == "PASS":
        al = [ref] + alt.split(",")
        dv[(chrom, int(pos))] = {al[int(i)][0] for i in re.split(r"[/|]", gt) if i not in (".", "")}
called = 0
for l in open(txt):
    c = l.rstrip("\n").split("\t")
    if l.startswith("#") or c[0] == "rsid" or c[0] not in table or table[c[0]][1] == "MT":
        continue
    if c[-1] not in ("--", "0") and "0" not in c[3:]:
        called += 1
got = moved = compared = nonref = bad = 0
for l in open(lifted):
    rid, chrom, pos, gt, ref, alt = l.rstrip("\n").split("\t")
    t = table.get(rid)
    if t is None or t[1] == "MT":
        continue
    idx = [i for i in re.split(r"[/|]", gt) if i not in (".", "")]
    if not idx:
        continue
    got += 1
    if int(pos) != int(t[4]):
        moved += 1
        if moved <= 5:
            print(f"moved {rid}: GRCh37 {t[1]}:{t[2]} lifts to {chrom}:{pos}, dbSNP's GRCh38 position is {t[4]}", file=sys.stderr)
        continue
    al = [ref] + alt.split(",")
    have = {al[int(i)][0] for i in idx}
    want = dv.get((chrom, int(pos)), {t[5]})
    compared += 1
    nonref += want != {t[5]}
    if have != want:
        bad += 1
        if bad <= 5:
            print(f"mismatch {rid} {chrom}:{pos} lifted {sorted(have)} calls {sorted(want)}", file=sys.stderr)
print(called, got, moved, compared, nonref, bad)
PY
}

# rt_checks LABEL "called lifted moved compared nonref mismatches"
rt_checks() {
  local called got moved compared nonref bad
  read -r called got moved compared nonref bad <<< "$2"
  echo "$1: ${called:-?} called sites off chrM, ${got:-?} lifted, ${moved:-?} moved, ${compared:-?} compared (${nonref:-?} non-reference), ${bad:-?} mismatches"
  check_ge "$1: lifted records, in percent of the called sites off chrM" "$(( ${got:-0} * 100 / (${called:-0} > 0 ? called : 1) ))" 95
  # The list's GRCh37 position and dbSNP's GRCh38 position can fall on two
  # copies of a duplicated region (one site in SMN1/SMN2 on the slice).
  check "$1: lifted records on another GRCh38 position than the table's (${moved:-?}) are at most 0.1% of ${got:-?}" \
    test "$(( ${moved:-1000000} * 1000 ))" -le "${got:-0}"
  check_ge "$1: non-reference sites compared" "${nonref:-0}" 100
  check_eq "$1: genotypes that differ from case 21's calls" "${bad:-x}" 0
}
R23=$(roundtrip HG002RT "${CASE_TMP}/export_23andme.txt" 23andme)
rt_checks "round trip, 23andMe layout" "$R23"
RAN=$(roundtrip HG002RA "${CASE_TMP}/export_ancestry.txt" ancestrydna)
rt_checks "round trip, AncestryDNA layout" "$RAN"
printf '#### Step 40 on the fixture\n\n%s of %s sites called; round trip (called, lifted, moved, compared, non-reference, mismatches): 23andMe %s, AncestryDNA %s\n\n' \
  "${CALLED:-?}" "$N_SITES" "$R23" "$RAN" >> "$E2E_NOTES"

# --- 5. red controls: the table shifted by one position ------------------------------
S38="${G}/reference/consumer/shift38.tsv.gz"
S37="${G}/reference/consumer/shift37.tsv.gz"
gzip -dc "$SITES" | awk -F'\t' -v OFS='\t' '/^#/ {print; next} {$5 = $5 + 1; print}' | gzip -c > "$S38"
gzip -dc "$SITES" | awk -F'\t' -v OFS='\t' '/^#/ {print; next} {$3 = $3 + 1; print}' | gzip -c > "$S37"
X38=$(docker run --rm --user "$(id -u):$(id -g)" -v "${G}:/genome" -v "${REPO}:/repo:ro" "$PYTHON_IMAGE" python3 /repo/bin/export_consumer.py export \
  --gvcf "/genome/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" --sites /genome/reference/consumer/shift38.tsv.gz --sample S38 --sex male \
  --par /repo/assets/par_grch38.bed --out-dir /genome/reference/consumer/out38 2>&1; echo "exit=$?")
echo "$X38" | tail -n 3
check "control: the table shifted by one GRCh38 position stops the export" has 'exit=1$' "$X38"
check "control: it names the reference bases" has 'another reference base than the site table' "$X38"
check "control: it writes nothing" test ! -e "${G}/reference/consumer/out38/S38_23andme_v5.zip"
docker run --rm --user "$(id -u):$(id -g)" -v "${G}:/genome" -v "${REPO}:/repo:ro" "$PYTHON_IMAGE" python3 /repo/bin/export_consumer.py export \
  --gvcf "/genome/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" --sites /genome/reference/consumer/shift37.tsv.gz --sample S37 --sex male \
  --par /repo/assets/par_grch38.bed --out-dir /genome/reference/consumer/out37
unz "${G}/reference/consumer/out37/S37_23andme_v5.zip" "${CASE_TMP}/shift37.txt"
RS=$(roundtrip HG002RS "${CASE_TMP}/shift37.txt" 23andme)
read -r _ s_got _ _ _ s_bad <<< "$RS"
echo "control, GRCh37 positions shifted by one: ${RS}"
check "control: with GRCh37 positions shifted by one the round trip fails (lifted ${s_got:-?} of ${CALLED:-?}, ${s_bad:-?} mismatches)" \
  test "$(( ${s_got:-0} * 100 / (${CALLED:-0} > 0 ? CALLED : 1) ))" -lt 95 -o "${s_bad:-0}" -gt 0
rm -rf "${G:?}/reference/consumer/out38" "${G:?}/reference/consumer/out37" "${G:?}/HG002RT" "${G:?}/HG002RA" "${G:?}/HG002RS" "${G:?}/reference_hg19" "${G:?}/liftover" "$S38" "$S37"

# --- 6. the Nextflow module on the same gVCF ------------------------------------------
if command -v nextflow >/dev/null; then
  NF_OUT="${G}/nf-consumer" WORK="${E2E_WORK}/nf-work-consumer" SHEET="${CASE_TMP}/samplesheet.csv"
  rm -rf "$NF_OUT" "$WORK"
  printf 'sample,vcf,vcf_index,gvcf,gvcf_index,sex\n%s,%s,%s,%s,%s,male\n' "$SAMPLE" "$VCF" "${VCF}.tbi" "$GVCF" "${GVCF}.tbi" > "$SHEET"
  (
    cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "$WORK" \
      --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
      --tools consumer_export --consumer_sites "$SITES" --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB
  ) 2>&1 | tee "$STEP_LOG"
  check_eq "nextflow run (consumer_export) exits 0" "${PIPESTATUS[0]}" 0
  NZ="${NF_OUT}/${SAMPLE}/consumer_export/${SAMPLE}_23andme_v5.zip"
  NA="${NF_OUT}/${SAMPLE}/consumer_export/${SAMPLE}_ancestrydna_layout.zip"
  if [ -s "$NZ" ] && [ -s "$NA" ]; then
    unz "$NZ" "${CASE_TMP}/nf_23andme.txt"
    unz "$NA" "${CASE_TMP}/nf_ancestry.txt"
    check "Nextflow and step 40 write the same 23andMe file (the date line aside)" \
      cmp -s <(grep -v '^# Written' "${CASE_TMP}/nf_23andme.txt") <(grep -v '^# Written' "${CASE_TMP}/export_23andme.txt")
    check "Nextflow and step 40 write the same AncestryDNA file (the date line aside)" \
      cmp -s <(grep -v 'Written 20' "${CASE_TMP}/nf_ancestry.txt") <(grep -v 'Written 20' "${CASE_TMP}/export_ancestry.txt")
    check "Nextflow and step 40 write the same summary" cmp -s "${NF_OUT}/${SAMPLE}/consumer_export/${SAMPLE}_consumer_export.tsv" "$SUM"
  else
    fail "the Nextflow run wrote no consumer_export zips under ${NF_OUT}/${SAMPLE}/consumer_export"
  fi
  rm -rf "$NF_OUT" "$WORK"
else
  fail "nextflow is not on PATH"
fi

rm -rf "${OUT:?}" "${SITES:?}"
finish
