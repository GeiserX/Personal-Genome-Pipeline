# Genome Viewer: See the Reads Behind a Finding

## What it does

`scripts/pgp-view.sh` opens a sample in [igv.js](https://github.com/igvteam/igv.js), the web version of the Integrative Genomics Viewer, served from your own computer. One page shows:

- the reference genome (`REF_FASTA`),
- the sample's reads (the BAM, or the CRAM from step 34),
- its small variants (the VCF from step 3),
- its structural variants (step 22's consensus, else Manta's VCF from step 4).

It opens at the locus you give. The HTML report (step 24) links to the viewer under the ClinVar hits, the secondary-findings genes and the mitochondrial variants. With the viewer running, a click on one of those links opens the reads at that position.

## Why

A variant call is a summary of the reads. Before you act on a finding, look at the reads under it: how many there are, how many carry the change, whether they are low-quality or mapped at MAPQ 0 (common in genes with a close copy elsewhere, such as PMS2 or CYP2D6), and whether they all start at the same place (a sign of PCR duplicates). Several report notes say "check the reads". This is where you do that.

## Tool

- **igv.js** 3.8.10 (`IGV_JS_VERSION` in `scripts/pgp-view.sh`, checked against a pinned sha256) draws the reads and variants in the browser.
- **`assets/viewer/serve.py`** (Python standard library) serves the page and the sample's files on `127.0.0.1`.

## Docker Image

None. The viewer runs on the host's `python3`, so it needs no image and no port publishing.

## Command

```bash
GENOME_DIR=/path/to/data ./scripts/pgp-view.sh <sample> [locus]
```

`locus` is `chr17`, `chr17:43045000` or `chr17:43044000-43046000`. Without one the viewer opens the first chromosome of the reference. Zoom in to a few kilobases or less to see the reads. The script prints the address to open:

```
Viewer for your_name: http://127.0.0.1:8765/#chr17:43045000
Press Ctrl-C to stop it.
```

Open it in a browser. To jump somewhere else, type a locus in igv.js's search box, or change the part after `#` in the address. Gene names do not work in the search box: looking a name up would need igv.org (see Privacy below).

| Variable | Default | Meaning |
|---|---|---|
| `PGP_VIEW_PORT` | `8765` | The port. The report's links use 8765, so keep it unless the port is taken. |
| `READS` | `<sample>/aligned/<sample>_sorted.bam`, else `..._sorted.cram` | The reads; `none` leaves them out. `ALIGN_DIR` picks another folder, as in the steps. |
| `VCF` | `<sample>/vcf/<sample>.vcf.gz` | The small-variant VCF (`VCF_DIR` as in the steps). |
| `SV_VCF` | `<sample>/sv_merged/<sample>_sv_consensus.vcf.gz`, else `<sample>/manta/results/variants/diploidSV.vcf.gz` | The structural-variant VCF. |
| `REF_FASTA` | as in every step | The reference. Its `.fai` must sit beside it. |

The script leaves out a file it does not find, and says so. A BAM or CRAM needs its index (`.bai`, `.crai`) beside it. A VCF with a `.tbi` is read in pieces. One without a `.tbi` is read whole, which is fine for an SV VCF.

The command needs `python3` on the host, but no Docker. The first run downloads igv.js (version 3.8.10, the `IGV_JS_VERSION` line in the script) into `GENOME_DIR/tools/igv.js-3.8.10/` and checks its sha256. After that it needs no network.

## Privacy

Nothing is uploaded:

- The server (`assets/viewer/serve.py`) listens on `127.0.0.1` only. Other machines on your network cannot reach it.
- It serves only the page, igv.js and the files listed above. It sends nothing else from `GENOME_DIR`.
- Every response carries a Content-Security-Policy that lets the page talk to this server and nowhere else. igv.js has features that call outside services: the genome list and gene-name search on igv.org, and BLAT. The browser blocks those calls. The viewer turns off the genome list, which igv.js would otherwise fetch at start.
- It refuses a request whose `Host` header is not `127.0.0.1` or `localhost`. That stops a web page from reaching the server through a domain name that points at your own machine.

The only download is igv.js itself, once, from jsDelivr's copy of the npm release. It is checked against the pinned sha256 before use.

## On a remote machine

If the data and the pipeline are on another machine (a server, a NAS), run the viewer there and forward the port over ssh from the computer you browse on:

```bash
# on your computer: forward local port 8765 to the remote machine's 127.0.0.1:8765
ssh -N -L 8765:127.0.0.1:8765 user@remote-machine

# on the remote machine (a second ssh session)
GENOME_DIR=/path/to/data ./scripts/pgp-view.sh your_name chr17:43045000
```

Then open `http://127.0.0.1:8765/` on your computer. The report's links work the same way. The reads travel inside the ssh connection, and only the ranges igv.js asks for. Do not expose the port any other way: the server has no login, so anyone who can reach it can read the sample.

## Output

Nothing is written apart from igv.js in `GENOME_DIR/tools/` on the first run. Stop the viewer with Ctrl-C.

## Runtime estimate

The server starts at once. How fast a view loads depends on the disk and the read depth: igv.js reads only the part of the BAM or CRAM under the view.

## Tested by

`tests/test_viewer.py`, in CI's unit tests, on a synthetic sample it writes itself (a random 3 kb reference, 30 reads, one variant, one deletion). It checks that:

- the page's config names the sample's files and the locus;
- the server listens on 127.0.0.1 and not on the machine's other address;
- it answers byte ranges of the BAM exactly, and refuses unknown paths and a foreign `Host`;
- in headless Chrome, igv.js starts on the config, the policy blocks nothing it needs, and it reads the BAM, its index, the reference and the VCF;
- the report links a ClinVar hit to the viewer at its position.

## Notes

- **The report links assume port 8765.** With another `PGP_VIEW_PORT`, type the locus into the viewer yourself.
- **CRAM** is decoded in the browser against the same reference. It must be the FASTA the CRAM was written with.
- **Gene names** are not searchable (see Privacy). Use the coordinates the report shows.
- igv.js is MIT-licensed. Reference: Robinson JT, Thorvaldsdottir H, Turner D, Mesirov JP. igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer (IGV). *Bioinformatics* 39(1), btac830 (2023). https://doi.org/10.1093/bioinformatics/btac830
