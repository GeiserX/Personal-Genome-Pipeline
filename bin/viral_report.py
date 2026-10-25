#!/usr/bin/env python3
"""viral_report.py: the viral DNA in a sample's reads that do not map to the human genome.

  viral_report.py --sample S --viral-fasta viral_genomes.fa
                  --viral-coverage S_viral_coverage.tsv --viral-depth S_viral_depth.tsv
                  --human-coverage S_human_coverage.tsv --out S_viral.tsv

scripts/47-viral.sh and the VIRAL_READS process realign the reads the human
alignment left unmapped (and the reads on chrEBV, the EBV genome GRCh38's
analysis set carries) to the viral genomes setup.sh --viral-data installs, then
write three samtools files this script reads:
  --viral-coverage  samtools coverage of that viral BAM: reads per virus
  --viral-depth     samtools depth -aa of the same BAM: the depth at every base
                    of every viral genome, zero included
  --human-coverage  samtools coverage of the sample's BAM: the human depth
                    (chr1 to chr22, over the bases with at least one read)
Both coverages count duplicates and leave out unmapped, secondary and QC-failed
reads (--ff UNMAP,SECONDARY,QCFAIL), so the two depths are measured alike.

One row per viral genome, in the order of the FASTA. Columns: virus, accession,
length, reads, breadth (the fraction of the genome with at least one read),
mean_depth, mode_depth (the most common depth over the genome's bases),
human_depth, relative_depth (mode_depth / human_depth), copies_per_cell
(2 x mean_depth / human_depth: a diploid cell holds two copies of each human
autosome), call, note.

The call:
  inherited integrated HHV-6 likely
      HHV-6A or HHV-6B with relative_depth >= 0.30. Liu et al. 2020 (PLoS Genet
      16:e1008915, doi:10.1371/journal.pgen.1008915) took the mode of the
      per-base depth over the viral genome against the WGS depth, and called
      integration at 0.30 or more: one integrated copy per cell gives about
      0.5, and viraemia gives 10 to 1000 times less. When both species pass,
      the one with the higher relative depth is called and the other row is
      read as reads shared by the two genomes.
  lab control (phiX)
      Illumina's phiX spike-in, added in the lab; never from the sample.
  viral DNA found     at least MIN_READS reads over at least FOUND_BREADTH of the genome
  a few reads         at least MIN_READS reads, over less than FOUND_BREADTH
  not found           fewer than MIN_READS reads
MIN_READS and FOUND_BREADTH are this pipeline's choices, not published values.
None of these is a diagnosis of infection.

Used by scripts/47-viral.sh and the VIRAL_REPORT process (bin/ is on the task
PATH); tests/test_viral_report.py checks it. Standard library only.
"""
import argparse
import csv
import re
import sys

MIN_READS = 10             # reads before a virus counts as seen at all
FOUND_BREADTH = 0.10       # the fraction of a genome those reads must cover to say "viral DNA found"
ICI_RELATIVE_DEPTH = 0.30  # Liu et al. 2020: mode viral depth / WGS depth for integrated HHV-6
HHV6 = {"NC_001664": "HHV-6A", "NC_000898": "HHV-6B"}
PHIX = {"NC_001422"}
RNA_VIRUSES = {"NC_004102"}  # HCV: an RNA virus with no DNA stage
HEADER = ["virus", "accession", "length", "reads", "breadth", "mean_depth", "mode_depth", "human_depth",
          "relative_depth", "copies_per_cell", "call", "note"]
AUTOSOME = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2])$")

CALL_ICI = "inherited integrated HHV-6 likely"
CALL_CONTROL = "lab control (phiX)"
CALL_FOUND = "viral DNA found"
CALL_FEW = "a few reads"
CALL_NONE = "not found"


def die(msg):
    sys.exit(f"viral_report.py: {msg}")


def read_fasta_titles(path):
    """[(accession, title)] in file order: the first word of each '>' line and the rest."""
    out = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                acc, _, title = line[1:].strip().partition(" ")
                out.append((acc, title.strip() or acc))
    if not out:
        die(f"{path} has no sequence")
    return out


def read_coverage(path):
    """samtools coverage rows by contig name: numreads, length, covbases, meandepth."""
    rows = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        header = None
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = [c.lstrip("#") for c in cols]
                continue
            if header is None or len(cols) != len(header):
                die(f"{path} is not samtools coverage output")
            r = dict(zip(header, cols))
            rows[r["rname"]] = {"reads": int(r["numreads"]), "length": int(r["endpos"]) - int(r["startpos"]) + 1,
                                "covbases": int(r["covbases"]), "meandepth": float(r["meandepth"])}
    return rows


def human_depth(cov):
    """Mean depth over chr1 to chr22's bases with at least one read, or None.
    Bases with no read (assembly gaps, regions a targeted or sliced BAM lacks)
    are left out, so the value is the depth where reads are."""
    total = covered = 0.0
    for name, r in cov.items():
        if AUTOSOME.match(name) and r["covbases"] > 0:
            total += r["meandepth"] * r["length"]
            covered += r["covbases"]
    return total / covered if covered else None


def read_depth(path, lengths):
    """{contig: (breadth, mean, mode)} from samtools depth -aa (contig, pos, depth)."""
    hist, seen = {}, {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for n, line in enumerate(f, 1):
            cols = line.split("\t")
            if len(cols) != 3:
                die(f"{path} line {n} is not samtools depth output (contig, position, depth)")
            h = hist.setdefault(cols[0], {})
            d = int(cols[2])
            h[d] = h.get(d, 0) + 1
            seen[cols[0]] = seen.get(cols[0], 0) + 1
    out = {}
    for name, length in lengths.items():
        h = hist.get(name) or {0: length}
        # No position at all is a genome with no read (samtools depth may leave
        # it out); some positions but not all is a depth file cut short.
        if seen.get(name, 0) not in (0, length):
            die(f"{path} has {seen.get(name, 0)} positions for {name}, its length is {length} (run samtools depth -aa)")
        breadth = 1.0 - h.get(0, 0) / length
        mean = sum(d * c for d, c in h.items()) / length
        mode = min(h, key=lambda d: (-h[d], d))   # the most common depth; a tie goes to the lower one
        out[name] = (breadth, mean, mode)
    return out


def short_name(acc, title):
    if acc.split(".")[0] in HHV6:
        return HHV6[acc.split(".")[0]]
    t = re.sub(r",? (complete (virion )?genome|genome|DNA)\b.*$", "", title)
    return t.strip() or acc


def call_row(acc, reads, breadth, relative):
    """The call and its note for one viral genome (before the HHV-6 species pick)."""
    base = acc.split(".")[0]
    if base in PHIX:
        if reads:
            return CALL_CONTROL, "Illumina adds phiX to the run as a control: these reads were added in the lab."
        return CALL_NONE, "Illumina's phiX control: none of its reads reached this sample's data."
    if base in HHV6 and relative is not None and relative >= ICI_RELATIVE_DEPTH:
        return CALL_ICI, (f"Depth over the genome at {relative:.2f} of the human depth (Liu et al. 2020 call "
                          f"integration at {ICI_RELATIVE_DEPTH:.2f} or more; one copy per cell gives about 0.5). "
                          "Inherited integrated HHV-6 is in every cell and is passed on to about half of the "
                          "children; it can make a blood test read as an active HHV-6 infection.")
    if reads < MIN_READS:
        return CALL_NONE, ""
    note = ""
    if base in HHV6:
        note = ("Below the integration threshold of Liu et al. 2020: not the one-copy-per-cell pattern of "
                "inherited HHV-6." if relative is not None else
                "No human depth (no read on chr1 to chr22): integration is not assessed.")
    elif base in RNA_VIRUSES:
        note = ("An RNA virus with no DNA stage: DNA reads point to a contamination or an integration "
                "(Moustafa et al. 2017).")
    return (CALL_FOUND if breadth >= FOUND_BREADTH else CALL_FEW), note


def build_rows(titles, vcov, vdepth, hdepth):
    rows = []
    for acc, title in titles:
        if acc not in vcov:
            die(f"the viral coverage has no row for {acc}: was the BAM aligned to this FASTA?")
        reads = vcov[acc]["reads"]
        breadth, mean, mode = vdepth[acc]
        relative = mode / hdepth if hdepth else None
        copies = 2 * mean / hdepth if hdepth else None
        call, note = call_row(acc, reads, breadth, relative)
        rows.append({"virus": short_name(acc, title), "accession": acc, "length": vcov[acc]["length"],
                     "reads": reads, "breadth": f"{breadth:.4f}", "mean_depth": f"{mean:.3f}", "mode_depth": mode,
                     "human_depth": f"{hdepth:.2f}" if hdepth else ".",
                     "relative_depth": f"{relative:.3f}" if relative is not None else ".",
                     "copies_per_cell": f"{copies:.4f}" if copies is not None else ".",
                     "call": call, "note": note, "_relative": relative or 0.0})
    # Both HHV-6 species pass: the one with the higher relative depth is the
    # integrated genome; the other's depth is reads shared by the two.
    ici = [r for r in rows if r["call"] == CALL_ICI]
    if len(ici) > 1:
        keep = max(ici, key=lambda r: r["_relative"])
        for r in ici:
            if r is not keep:
                r["call"] = CALL_FOUND
                r["note"] = (f"Reads shared with {keep['virus']}, which has the higher depth and is the integrated "
                             f"species called.")
    for r in rows:
        del r["_relative"]
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--viral-fasta", required=True)
    ap.add_argument("--viral-coverage", required=True)
    ap.add_argument("--viral-depth", required=True)
    ap.add_argument("--human-coverage", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    titles = read_fasta_titles(a.viral_fasta)
    vcov = read_coverage(a.viral_coverage)
    vdepth = read_depth(a.viral_depth, {acc: vcov[acc]["length"] for acc, _ in titles if acc in vcov})
    rows = build_rows(titles, vcov, vdepth, human_depth(read_coverage(a.human_coverage)))
    with open(a.out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    found = [f"{r['virus']} ({r['call']})" for r in rows if r["call"] not in (CALL_NONE,)]
    print(f"{a.sample}: {len(rows)} viral genomes; " + ("; ".join(found) if found else "no viral DNA found"))


if __name__ == "__main__":
    main()
