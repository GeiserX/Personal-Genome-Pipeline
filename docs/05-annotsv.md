# Step 5: Structural Variant Annotation (AnnotSV)

## What This Does
Classifies every structural variant from Manta that passes duphold's depth filter (step 15; when that file is missing or older than Manta's calls, the script annotates Manta's calls and says so) using ACMG guidelines (class 1-5), adding gene overlap, population frequency, and clinical significance.

## Why
Raw Manta output contains thousands of SVs with no clinical interpretation. AnnotSV tells you which ones matter by cross-referencing known pathogenic SVs, gene databases, and population data.

## Tool
- **AnnotSV** — ACMG-compliant structural variant annotation and classification

## Docker Image
- `ANNOTSV_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Annotation Data
The image holds AnnotSV's code only. Its annotation data (genes, known pathogenic SVs, population frequencies) is a separate 5.3 GB download that unpacks to about 20 GB, and AnnotSV exits with an error without it. The server is slow (about 0.8 MB/s measured from a GitHub runner), so the download can take 1-2 hours; `curl -C -` resumes it. `./scripts/setup.sh` downloads and unpacks it into `${GENOME_DIR}/annotsv_annotations/`. To do it by hand:

```bash
curl -fL -C - -o ${GENOME_DIR}/Annotations_Human_3.5.tar.gz \
  https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz
mkdir -p ${GENOME_DIR}/annotsv_annotations
tar -xzf ${GENOME_DIR}/Annotations_Human_3.5.tar.gz -C ${GENOME_DIR}/annotsv_annotations
```

The script checks for `${GENOME_DIR}/annotsv_annotations/Annotations_Human/Genes/GRCh38` and stops with a pointer to `setup.sh` when it is missing. `run-all.sh` reports step 5 as skipped in that case.

### gnomAD-SV frequencies and MANE Select transcripts

`setup.sh` adds two files to the bundle. Both are pinned in `versions.env` with a sha256, and the run manifest records the pins through its copy of `versions.env`. `./scripts/setup.sh --annotsv-data ${GENOME_DIR}` adds them to an existing install without running the rest of setup.

- **gnomAD-SV common SVs** (`GNOMAD_SV_VERSION`, `GNOMAD_SV_SHA256`). Without population frequencies, a deletion that many people carry is scored as if nobody did, and can reach class 4. `setup.sh` downloads gnomAD-SV's sites file (about 1.7 GB) and keeps the records AnnotSV uses for its benign set: FILTER PASS, type DEL, DUP, INS or INV, AN at least 1,000, AF at least 0.1% and at least 5 homozygous carriers. It writes them as `Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad.v4.0.sv.chrAll.from-v4.1-sites.vcf.gz`, the file name and layout AnnotSV 3.5's gnomAD-SV v4 reader expects, and deletes the download. The v4.1 sites file has no controls-and-biobanks subset, so the frequencies of all genomes are written under the `controls_and_biobanks_*` keys that reader takes. `setup.sh` then runs AnnotSV once, as you, on a one-deletion VCF: AnnotSV merges these records into its `benign_*_SV_GRCh38.sorted.bed` files and deletes the file. Only after that merge does `gnomad-sv.installed` beside it record the release. Sample runs find nothing left to merge, so parallel Nextflow `ANNOTSV` tasks never rewrite those shared files at once, and tasks that run as root never take them over.
- **MANE Select transcripts** (`MANE_RELEASE`, `MANE_SUMMARY_SHA256`). For each gene AnnotSV reports one transcript in `Tx`, `Frameshift`, `Location`, `Location2` and `Exon_count`, by default the one with the longest overlapped coding sequence, which is often a predicted `XM_` transcript. The ranking reads `Frameshift` and `Location`, so the transcript can change the class. `setup.sh` takes the MANE Select RefSeq ids (`NM_`) from NCBI's MANE summary into `annotsv_annotations/mane_select_refseq.txt`, and the step passes that file as `-txFile`, so AnnotSV reports the MANE Select transcript of a gene whenever the SV overlaps it. The exon numbers then match ClinVar and the literature.

What the frequencies change in the ranking, by AnnotSV's rules:

- A loss gets criterion 2F/4O (score -1.00) when it lies entirely inside a gnomAD-SV deletion of at least 1% frequency (`-benignAF`, default 0.01); `B_loss_source`, `B_loss_coord` and `B_loss_AFmax` name that record. A call that reaches a few bases past the gnomAD record does not get it.
- 2F counts only when no stronger criterion of section 2 applies. A deletion that breaks the reading frame of a haploinsufficient gene (2E-1, +0.90) keeps its class even when it is common. With MANE Select as the transcript, an intronic deletion is no longer read as a frameshift of an `XM_` exon, so the benign criterion can apply.

Step 05 and the Nextflow `ANNOTSV` stop with a pointer to `--annotsv-data` when `mane_select_refseq.txt` or `gnomad-sv.installed` is missing, and `validate-setup.sh` fails on either.

## Command
```bash
./scripts/05-annotsv.sh your_sample
```

What the script runs:

```bash
source versions.env   # from the repository root
# The input is step 15's filtered calls; without them, Manta's
# /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz
docker run --rm --user "$(id -u):$(id -g)" \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  "${ANNOTSV_IMAGE}" \
  AnnotSV \
    -SVinputFile /genome/${SAMPLE}/duphold/${SAMPLE}_sv_filtered.vcf.gz \
    -outputFile /genome/${SAMPLE}/annotsv/${SAMPLE}_sv_annotated.tsv \
    -genomeBuild GRCh38 \
    -annotationMode both \
    -annotationsDir /genome/annotsv_annotations \
    -txFile /genome/annotsv_annotations/mane_select_refseq.txt
```

To annotate another SV VCF (for example Sniffles2 output), set `SV_VCF` to its host path; the file must be inside `${GENOME_DIR}`.

## Output
- `${GENOME_DIR}/${SAMPLE}/annotsv/${SAMPLE}_sv_annotated.tsv` — main annotated output (one row per SV, with ACMG class)
- Columns include: SV type, coordinates, overlapping genes and the MANE Select transcript (`Tx`), the benign-set frequency (`B_loss_AFmax`, `B_gain_AFmax`), ACMG classification, ClinVar hits

## ACMG Classification
| Class | Meaning | Action |
|---|---|---|
| 1 | Benign | Ignore |
| 2 | Likely benign | Ignore |
| 3 | Variant of uncertain significance (VUS) | Review if in known disease gene |
| 4 | Likely pathogenic | Investigate — check gene, inheritance, phenotype |
| 5 | Pathogenic | Investigate — known disease-causing SV |

## Important Notes
- **Class 4-5 = pathogenic/likely pathogenic** — these require manual review
- SVs >5MB in short-read WGS are usually artifacts from segmental duplications — do not trust large calls blindly
- Most SVs will be class 2-3 (benign/VUS) — this is normal for a healthy genome
- Input: step 15's `duphold/<sample>_sv_filtered.vcf.gz`, Manta's calls after duphold's depth filter, the file the Nextflow `ANNOTSV` reads. Without it (or when Manta's `diploidSV.vcf.gz` is newer), the script annotates Manta's calls and says so. `SV_VCF=<path>` picks another file, for example Sniffles2's. Do not point it at Manta's `candidateSV.vcf.gz`: those are unfiltered candidates, not calls
- The annotation databases are not in the Docker image: they live in `${GENOME_DIR}/annotsv_annotations/` (see Annotation Data above)
