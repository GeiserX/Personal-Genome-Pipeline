# Why run locally?

## Cost Comparison

| Approach | Cost | What You Get | Data Privacy |
|---|---|---|---|
| **This pipeline** | $0 (free, open source) | Every step in the [pipeline overview](pipeline-overview.md#what-a-default-run-covers) | No step uploads your data (see below) |
| Clinical WGS interpretation | $500-5,000 | 1-page report, selected genes only | Lab retains your data |
| Nebula/Dante report | $0-200 (included/add-on) | Web dashboard, limited depth | Data on company servers |
| 23andMe Health | $229 | ~10 health reports from array data | Data shared with research partners |
| Genetic counselor consultation | $200-500/hour | Expert interpretation of specific findings | HIPAA-protected |

**The pipeline is complementary, not a replacement.** Use it for comprehensive self-analysis, then bring specific findings to a genetic counselor or physician for clinical interpretation.

## Privacy and Security

Your genome is the most permanent piece of personal data you have. Unlike a password, you cannot change it if it leaks.

**This pipeline keeps your data local:**
- No pipeline script sends your reads, alignments or variants anywhere, and the pipeline has no telemetry, analytics or tracking of its own. Every step runs its containers with networking turned off (`--network none`); only `setup.sh` gives a container network, to install data and the opt-in Cyrius.
- The one place the scripts send your data off the machine is your own choice: [step 14](14-imputation-prep.md) and the [chip data guide](chip-data-guide.md#optional-imputation) prepare files for an imputation server, and uploading them is a separate step you take or skip.
- [Reference setup](00-reference-setup.md) downloads the databases and pulls the images once. After that, a run still makes the network calls below. None of these calls carries sample data.

### Network calls during a run

Read from the scripts on 2026-10-05. The step containers have no network (`run_in` in `scripts/lib/common.sh`, and `--network none` for every Nextflow task), so these are downloads the scripts make themselves, outside a container.

| Where | What it fetches | When |
|---|---|---|
| Any step | its Docker image, from Docker Hub or quay.io | the first time an image is used, if `setup.sh` did not pull it |
| Step 13 (VEP) | the 26 GB VEP cache from Ensembl | when the script is run by hand without the cache (`run-all.sh` skips the step instead) |
| `setup.sh --cyrius` (step 21, opt-in) | `cyrius==1.1.1` and its dependencies from PyPI, every wheel checked against its sha256 in `scripts/cyrius-constraints.txt` | once, when you ask for it; step 21 itself downloads nothing |
| Step 25 (PRS) | PGS scoring files from the PGS Catalog (EBI) | the script's first run, then cached in `prs_scores/` (`run-all.sh` downloads nothing: it scores the files already there and skips the step without them) |
| `run-all.sh` | the Nextflow release `NXF_VER` names (`NEXTFLOW_VERSION` in `versions.env`), from GitHub | the first run, when that release is not installed yet |
| Step 26 (ancestry, opt-in) | 1000 Genomes sites and population labels | first run, then cached |
| Step 4b (GRIDSS, opt-in) | the ENCODE blacklist BED from GitHub | first run, then cached |

`setup.sh` installs the IPD-IMGT/HLA release step 8 types against (pinned, `HLA_DB_RELEASE`), and with `--parascopy-data` and `--kir-data` the data of the opt-in step 35 and of KIR typing; the steps read them and download nothing. MultiQC (step 28) runs with `--no-version-check`.

### Opening the reports

A report that loads a file from the internet tells that server your IP address and when you opened it, but nothing from the report. The [e2e test](testing.md) lists every external address in the HTML reports it produces. On 2026-10-02 it found:

- **Step 24 HTML report:** loads nothing; its only external address is a plain link to the project.
- **PharmCAT report (step 7):** loads a Font Awesome stylesheet from `maxcdn.bootstrapcdn.com` and two scripts from `oss.maxcdn.com`.
- **indexcov report (step 16):** loads Chart.js from `cdnjs.cloudflare.com` and jQuery from `code.jquery.com`. Offline, the page opens without its plots.
- **Nextflow report and timeline:** load nothing.

Not checked: the CPSR (step 17) and MultiQC (step 28) reports, which the test cannot produce.

**Recommendations for securing your data:**
- Store genomic data on an encrypted filesystem (LUKS on Linux, FileVault on macOS, BitLocker on Windows)
- Never upload raw FASTQ/BAM/VCF files to unencrypted cloud storage
- If using a NAS, enable encryption at rest
- Be cautious with VCF files in particular — they are small enough to accidentally email or upload
- Consider who has physical access to the machine where your data is stored
- If you delete your data, use `shred` (Linux) or secure erase — standard file deletion leaves data recoverable

**GDPR note:** If you are in the EU, your genomic data is classified as "special category personal data" under GDPR Article 9. Processing it locally for personal use is lawful. Sharing it with third parties (including cloud services) may require explicit consent and appropriate safeguards.

