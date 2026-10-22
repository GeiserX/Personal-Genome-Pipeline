#!/usr/bin/env bash
# Step 35's HBA1/HBA2 count from read depth (bin/hba_copy_number.py) on
# synthetic reads, with a planted one-gene deletion. The fixture has no reads
# at the alpha-globin cluster, so the case writes its own: 150 bp reads cut
# from the fixture reference every 5 bp over chr16:95,000-235,000 (30x),
# with the read group of the GIAB slice BAM. Two samples:
#   HG002hba4  the reads as cut: four alpha-globin genes
#   HG002hba3  every second read over HBA2 left out (one HBA2 copy lost, as
#              -a4.2): three
# Each is merged with the chr5 and chr20 reads of the GIAB slice BAM, so the
# Parascopy part of step 35 runs as in pgx-4-parascopy.sh (its stand-in
# reference with an N-filled chr11, background windows on chr20).
# Then the HBA block of modules/local/parascopy/main.nf, cut from the file,
# runs in the same image on the same BAM: its table must equal step 35's.
# Synthetic reads place every read uniquely and carry no GC bias, so this
# shows the counting and the wiring, not the accuracy on a real genome.
# Last, the same block runs on HG002's real Illumina reads at the cluster,
# streamed from GIAB's 60x GRCh38 BAM (the fixture's source) and sampled to
# about 30x: Paraphase on HG002's HiFi reads calls two HBA1 and two HBA2
# (paralogs-paraphase-1), so four is expected. One genome at one depth.
. "$(dirname "$0")/lib.sh"

echo "+ scripts/setup.sh --parascopy-data"
"${REPO}/scripts/setup.sh" --parascopy-data "$GENOME_DIR"

REF="${GENOME_DIR}/reference/parascopy_hba_ref.fasta"
python3 - "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" "$REF" <<'PY'
import shutil, sys
with open(sys.argv[2], "wb") as out, open(sys.argv[1], "rb") as src:
    shutil.copyfileobj(src, out)
    out.write(b">chr11\n")
    left = 135086622            # GRCh38 chr11, as the homology table lists it
    while left:
        n = min(left, 60)
        out.write(b"N" * n + b"\n")
        left -= n
PY
sam faidx "reference/$(basename "$REF")"

W="hba_e2e"
mkdir -p "${GENOME_DIR}/${W}"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam" "${GENOME_DIR}/${W}/slice.bam"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam.bai" "${GENOME_DIR}/${W}/slice.bam.bai"
sam view -b -o "${W}/slice_chr5_chr20.bam" "${W}/slice.bam" chr5 chr20
sam view -H "${W}/slice.bam" > "${GENOME_DIR}/${W}/header.sam"
sam faidx reference/GRCh38_no_alt_analysis_set.fasta chr16:95001-235000 > "${GENOME_DIR}/${W}/chr16.fa"

for P in HG002hba4 HG002hba3; do
  python3 - "${GENOME_DIR}/${W}/chr16.fa" "${GENOME_DIR}/${W}/header.sam" "$P" "${GENOME_DIR}/${W}/${P}_syn.sam" <<'PY'
import sys
fa, header, sample, out = sys.argv[1:]
seq = "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()
start0 = 95000                      # 0-based start of chr16:95001-235000
rg = next((l.split("\t")[1][3:] for l in open(header) if l.startswith("@RG")), None)
hba2 = (172875, 173710)             # bin/hba_copy_number.py's HBA2, 0-based
L, step, n = 150, 5, 0
with open(out, "w") as f:
    f.write(open(header).read())
    for i in range(0, len(seq) - L + 1, step):
        s = start0 + i
        if sample.endswith("3") and s < hba2[1] and s + L > hba2[0]:
            n += 1
            if n % 2:
                continue
        f.write(f"syn{s}\t0\tchr16\t{s + 1}\t60\t{L}M\t*\t0\t0\t{seq[i:i + L]}\t{'I' * L}"
                + (f"\tRG:Z:{rg}" if rg else "") + "\n")
PY
  mkdir -p "${GENOME_DIR}/${P}/aligned"
  sam sort -o "${W}/${P}_syn.bam" "${W}/${P}_syn.sam"
  sam merge -f "${P}/aligned/${P}_sorted.bam" "${W}/slice_chr5_chr20.bam" "${W}/${P}_syn.bam"
  sam index "${P}/aligned/${P}_sorted.bam"
done

BED="${GENOME_DIR}/${W}/parascopy_windows.bed"
awk 'BEGIN { for (s = 10050000; s < 10450000; s += 100) printf "chr20\t%d\t%d\n", s, s + 100 }' > "$BED"

# The module's HBA block, as written in the .nf file, made a plain shell
# script: Groovy's \\ at a line end is the shell's \, ${bam} and ${meta.id}
# become the sample's BAM and name.
MOD="${REPO}/modules/local/parascopy/main.nf"
BLOCK="${CASE_TMP}/module_hba.sh"
sed -n '/# HBA1\/HBA2 from depth (as scripts\/35-paralogs.sh)/,/rm -f hba_regions.bed hba_bedcov.tsv/p' "$MOD" \
  | sed 's/\\\\$/\\/' > "$BLOCK"
check_ge "the module's HBA block was cut from main.nf" "$(wc -l < "$BLOCK" | tr -d ' ')" 8

hba_value() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "$1" 2>/dev/null; }

for P in HG002hba4 HG002hba3; do
  want=${P#HG002hba}
  REF_FASTA="$REF" PARASCOPY_DEPTH_BED="$BED" run_step 35-paralogs.sh "$P"
  check_step_exit 35-paralogs.sh
  T="${GENOME_DIR}/${P}/paralogs/${P}_hba_copy_number.tsv"
  cat "$T" 2>/dev/null
  check_eq "${P}: step 35's HBA table says called" "$(hba_value "$T" status)" called
  check_eq "${P}: alpha-globin genes from depth" "$(hba_value "$T" alpha_globin_genes)" "$want"
  check "${P}: the SMN table is still written" test -s "${GENOME_DIR}/${P}/paralogs/${P}_smn_copy_number.tsv"

  # The module's block in the same image, on the same BAM
  M="${GENOME_DIR}/${W}/${P}_module"
  rm -rf "$M"; mkdir -p "$M"
  sed -e "s|\${bam}|/genome/${P}/aligned/${P}_sorted.bam|g" -e "s|\${meta.id}|${P}|g" "$BLOCK" > "${M}/block.sh"
  docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w "/genome/${W}/${P}_module" \
    "$PARASCOPY_IMAGE" sh -c 'PATH="/pgp-bin:$PATH"; sh -e block.sh' 2>&1 | tail -n 3
  check "${P}: the module's HBA block writes the same table as step 35" cmp -s "$T" "${M}/${P}_hba_copy_number.tsv"
  echo "hba: ${P} $(hba_value "$T" estimate) ($(hba_value "$T" alpha_globin_genes) genes)" >> "$E2E_NOTES"
done

# --- HG002's real Illumina reads at the alpha-globin cluster ------------------
GIAB_BAM="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
sam_https() {  # the image's libcurl has no CA bundle, so lend it the host's
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w "/genome/${W}" \
    "$SAMTOOLS_IMAGE" samtools "$@"
}
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && sleep $((i * 15)); done; return 1; }
R=HG002real
retry sam_https view -b -s 42.5 -o "${R}.bam" "$GIAB_BAM" chr16:95001-235000
rm -f "${GENOME_DIR}/${W}/$(basename "$GIAB_BAM").bai"
sam index "${W}/${R}.bam"
check_ge "HG002 Illumina reads at the cluster" "$(sam view -c "${W}/${R}.bam" 2>/dev/null)" 10000
M="${GENOME_DIR}/${W}/${R}_module"
mkdir -p "$M"
sed -e "s|\${bam}|/genome/${W}/${R}.bam|g" -e "s|\${meta.id}|${R}|g" "$BLOCK" > "${M}/block.sh"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w "/genome/${W}/${R}_module" \
  "$PARASCOPY_IMAGE" sh -c 'PATH="/pgp-bin:$PATH"; sh -e block.sh' 2>&1 | tail -n 3
cat "${M}/${R}_hba_copy_number.tsv" 2>/dev/null
check_eq "HG002 (real reads, about 30x): alpha-globin genes from depth" \
  "$(hba_value "${M}/${R}_hba_copy_number.tsv" alpha_globin_genes)" 4
echo "hba: HG002 real reads, estimate $(hba_value "${M}/${R}_hba_copy_number.tsv" estimate);" \
  "HBA2 $(hba_value "${M}/${R}_hba_copy_number.tsv" hba2_depth)x, HBA1 $(hba_value "${M}/${R}_hba_copy_number.tsv" hba1_depth)x," \
  "flanks $(hba_value "${M}/${R}_hba_copy_number.tsv" flank_left_depth)x and $(hba_value "${M}/${R}_hba_copy_number.tsv" flank_right_depth)x" >> "$E2E_NOTES"

rm -f "$REF" "${REF}.fai"
rm -rf "${GENOME_DIR:?}/${W}"
finish
