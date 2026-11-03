#!/usr/bin/env python3
"""mei_report.py: the Alu, L1 and SVA insertions MEGAnE calls, and the known disease insertions.

  mei_report.py --sample S --vcf S_megane.vcf --sites assets/mei_sites.tsv
                --depth S_mei_depth.tsv [--annotsv S_mei_annotsv.tsv]
                --out-calls S_mei.tsv --out-sites S_mei_sites.tsv
  mei_report.py --sample S --vcf S_megane.vcf --annotsv-input S_mei_annotsv_input.vcf

scripts/50-mei.sh and the MEI process run MEGAnE (Kojima et al. 2023, Nat Genet,
PMID 37169872) on the sample's BAM. Its VCF has one record per insertion it
found, ALT <INS:ME>, the repeat class in INFO SVTYPE (SINE/Alu, LINE/L1,
Retroposon/SVA, and older classes too), MEPRED=PASS when the inserted sequence
was identified as that element, the subfamily in INFO MEI (or MEI_lbp and
MEI_rbp, one per breakpoint), and the copy number in FORMAT CN (0, 1 or 2).

--out-calls: one row per insertion of Alu, L1 or SVA with MEPRED=PASS and at
least one copy. Columns: chrom, pos, end (MEGAnE's right breakpoint), element,
subfamily, copies, filter (MEGAnE's FILTER: PASS, or its flags such as LC, low
confidence), confidence (CN_conf: H or L, how sure the copy number is). Other
classes (DNA transposons, old LTRs) are left out: none of them is active in
human genomes, so a call of one is read as an artefact.

With --annotsv (AnnotSV's table of the same calls, which scripts/50-mei.sh
and the MEI_ANNOTSV process write from --annotsv-input), four more columns:
genes (the genes the insertion falls in, AnnotSV's Gene_name), location (each
gene with where in it the insertion falls, AnnotSV's Location: exon3,
intron2-intron2, ...), benign_af (AnnotSV's B_ins_AFmax, the highest allele
frequency of a benign insertion at that place in its benign set, gnomAD-SV
among them) and acmg_class (AnnotSV's ACMG_class). AnnotSV 3.5 scores
deletions and duplications only, so its class of an insertion is NA. Without
--annotsv, or for a call AnnotSV did not annotate, the four are ".".

--annotsv-input OUT writes the calls the table lists as the VCF AnnotSV reads
(ALT <INS:ME>, SVTYPE=INS, END=POS, and an ID, chrom_pos_element, that the
--annotsv table is joined on) and stops; --sites, --depth and the --out
options are not needed then.

--out-sites: every row of the catalog (assets/mei_sites.tsv) with three
columns added:
  status      seen        a listed insertion of the same element inside the
                          row's window (start to end)
              not seen    none, with at least MIN_DEPTH reads of depth over
                          the window
              not callable  the BAM has no reads of that contig, or the mean
                          depth over the window is under MIN_DEPTH: too few
                          reads to tell
  mean_depth  the mean depth over the window (samtools bedcov -c: summed depth
              / window length), "." when the contig is not in the BAM
  call        the insertion seen: position, subfamily, copies, filter; else "."
--depth is samtools bedcov -c over the catalog windows on the BAM's contigs:
chrom, start (0-based), end, gene, summed depth, reads. MIN_DEPTH is this
pipeline's choice, not a published value. "Not seen" is not proof of absence:
an SVA is about 2 kb of repeat, and short reads find its edges, not its body.

Used by scripts/50-mei.sh and the MEI_REPORT process (bin/ is on the task
PATH); tests/test_mei_report.py checks it. Standard library only.
"""
import argparse
import csv
import sys

MIN_DEPTH = 10  # mean depth over a catalog window before "not seen" is said
ELEMENTS = {"SINE/Alu": "Alu", "LINE/L1": "L1", "Retroposon/SVA": "SVA"}
CALL_HEADER = ["chrom", "pos", "end", "element", "subfamily", "copies", "filter", "confidence",
               "genes", "location", "benign_af", "acmg_class"]
ANNOTSV_COLUMNS = ("ID", "Annotation_mode", "Gene_name", "ACMG_class")
SITE_COLUMNS = ["gene", "chrom", "start", "end", "site", "element", "inheritance", "condition", "clinvar", "pmid"]


def die(msg):
    sys.exit(f"mei_report.py: {msg}")


def info_dict(text):
    return dict(x.split("=", 1) if "=" in x else (x, "") for x in text.split(";") if x)


def subfamily(info):
    """The first subfamily MEGAnE names: from MEI, else from either breakpoint."""
    for key in ("MEI", "MEI_lbp", "MEI_rbp"):
        name = info.get(key, "").split("|")[0].split(",")[0]
        if name and name not in ("NA", "pA"):
            return name
    return "."


def read_calls(path, contigs=None):
    """The calls the table lists; the ##contig lines go to contigs when it is a list."""
    calls, header = [], False
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if line.startswith("##"):
                if contigs is not None and line.startswith("##contig="):
                    contigs.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                header = True
                continue
            if not line.strip():
                continue
            if not header:
                die(f"{path}: a record before the #CHROM line (line {n}): not a VCF")
            f_ = line.rstrip("\n").split("\t")
            if len(f_) < 10:
                die(f"{path} line {n}: {len(f_)} columns, a VCF with one sample has 10")
            info = info_dict(f_[7])
            if "SVTYPE" not in info or "MEPRED" not in info:
                die(f"{path} line {n}: no SVTYPE or MEPRED in INFO: not MEGAnE's VCF")
            element = ELEMENTS.get(info["SVTYPE"])
            fmt = dict(zip(f_[8].split(":"), f_[9].split(":")))
            cn = fmt.get("CN", ".")
            if element is None or info["MEPRED"] != "PASS" or not cn.isdigit() or int(cn) < 1:
                continue
            calls.append({"chrom": f_[0], "pos": int(f_[1]), "end": info.get("MEI_rpos", f_[1]),
                          "element": element, "subfamily": subfamily(info), "copies": cn,
                          "filter": f_[6], "confidence": info.get("CN_conf", ".")})
    if not header:
        die(f"{path}: no #CHROM line: not a VCF")
    seen = {}
    for c in calls:
        key = f"{c['chrom']}_{c['pos']}_{c['element']}"
        seen[key] = seen.get(key, 0) + 1
        c["id"] = key if seen[key] == 1 else f"{key}_{seen[key]}"
    return calls


def write_annotsv_input(path, sample, calls, contigs):
    """The calls as the VCF AnnotSV reads: one <INS:ME> record each, END=POS,
    no length (MEGAnE gives none; AnnotSV then skips its minimum-size check)."""
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.writelines(c + "\n" for c in contigs)
        out.write('##ALT=<ID=INS:ME,Description="Insertion of mobile element">\n'
                  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n'
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                  f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for c in calls:
            gt = "0/1" if c["copies"] == "1" else "1/1"
            out.write(f"{c['chrom']}\t{c['pos']}\t{c['id']}\tN\t<INS:ME>\t.\t{c['filter']}\t"
                      f"SVTYPE=INS;END={c['pos']}\tGT\t{gt}\n")


def read_annotsv(path):
    """AnnotSV's table of the calls, by ID: genes, location, benign_af and
    acmg_class from its full row, the location of each gene from its split
    rows. AnnotSV writes an empty file when it had nothing to annotate."""
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
        head = rows[0].keys() if rows else []
    if not rows:
        return {}
    if any(c not in head for c in ANNOTSV_COLUMNS):
        die(f"{path}: not AnnotSV's table (a column of {', '.join(ANNOTSV_COLUMNS)} is missing)")
    out = {}
    for r in rows:
        a = out.setdefault(r["ID"], {"genes": ".", "location": [], "benign_af": ".", "acmg_class": "."})
        if r["Annotation_mode"] == "full":
            a["genes"] = ",".join(g for g in r["Gene_name"].split(";") if g) or "."
            a["benign_af"] = (r.get("B_ins_AFmax") or "").strip() or "."
            a["acmg_class"] = r["ACMG_class"].strip() or "."
        elif r["Annotation_mode"] == "split" and r["Gene_name"] and (r.get("Location") or "").strip():
            a["location"].append(f"{r['Gene_name']}:{r['Location'].strip()}")
    for a in out.values():
        a["location"] = ",".join(a["location"]) or "."
    return out


def read_sites(path):
    rows = [r for r in open(path) if r.strip() and not r.startswith("#")]
    reader = list(csv.DictReader(rows, delimiter="\t"))
    if not reader or any(c not in reader[0] for c in SITE_COLUMNS):
        die(f"{path}: no rows, or a column of {', '.join(SITE_COLUMNS)} is missing")
    return reader


def read_depth(path):
    """samtools bedcov -c over the 4-column BED: chrom, start, end, gene, sum, reads."""
    depth = {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if not line.strip():
                continue
            f_ = line.rstrip("\n").split("\t")
            if len(f_) < 6 or not f_[1].isdigit() or not f_[2].isdigit():
                die(f"{path} line {n}: not samtools bedcov -c output of a 4-column BED")
            length = int(f_[2]) - int(f_[1])
            depth[(f_[3], f_[0])] = float(f_[4]) / length if length > 0 else 0.0
    return depth


def chrom_rank(chrom):
    """chr1 to chr22, then chrX, chrY and the rest, as a karyotype is read."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return int(c) if c.isdigit() else {"X": 23, "Y": 24}.get(c, 25)


def site_status(site, calls, depth):
    key = (site["gene"], site["chrom"])
    if key not in depth:
        return "not callable", ".", "."
    mean = depth[key]
    lo, hi = int(site["start"]), int(site["end"])
    hits = [c for c in calls if c["chrom"] == site["chrom"] and c["element"] == site["element"] and lo <= c["pos"] <= hi]
    if hits:
        c = hits[0]
        call = f"{c['chrom']}:{c['pos']} {c['subfamily']} {c['copies']} cop{'y' if c['copies'] == '1' else 'ies'} {c['filter']}"
        return "seen", f"{mean:.1f}", call
    if mean < MIN_DEPTH:
        return "not callable", f"{mean:.1f}", "."
    return "not seen", f"{mean:.1f}", "."


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sites")
    ap.add_argument("--depth")
    ap.add_argument("--annotsv", help="AnnotSV's table of the --annotsv-input VCF")
    ap.add_argument("--out-calls")
    ap.add_argument("--out-sites")
    ap.add_argument("--annotsv-input", help="write the calls as AnnotSV's input VCF here, and stop")
    a = ap.parse_args(argv)
    contigs = []
    calls = read_calls(a.vcf, contigs)
    if a.annotsv_input:
        write_annotsv_input(a.annotsv_input, a.sample, calls, contigs)
        print(f"{a.sample}: {len(calls)} insertions written for AnnotSV")
        return
    for opt in ("sites", "depth", "out_calls", "out_sites"):
        if not getattr(a, opt):
            ap.error(f"--{opt.replace('_', '-')} is required (or give --annotsv-input)")
    ann = read_annotsv(a.annotsv) if a.annotsv else {}
    for c in calls:
        c.update(ann.get(c["id"], {"genes": ".", "location": ".", "benign_af": ".", "acmg_class": "."}))
    sites = read_sites(a.sites)
    depth = read_depth(a.depth)
    calls.sort(key=lambda c: (chrom_rank(c["chrom"]), c["chrom"], c["pos"]))
    with open(a.out_calls, "w") as out:
        out.write("\t".join(CALL_HEADER) + "\n")
        for c in calls:
            out.write("\t".join(str(c[k]) for k in CALL_HEADER) + "\n")
    with open(a.out_sites, "w") as out:
        out.write("\t".join(SITE_COLUMNS + ["status", "mean_depth", "call"]) + "\n")
        for s in sites:
            out.write("\t".join([s[k] for k in SITE_COLUMNS] + list(site_status(s, calls, depth))) + "\n")
    by = {}
    for c in calls:
        by[c["element"]] = by.get(c["element"], 0) + 1
    print(f"{a.sample}: {len(calls)} insertions ({', '.join(f'{k} {v}' for k, v in sorted(by.items())) or 'none'})")
    if a.annotsv:
        print(f"{a.sample}: AnnotSV annotated {sum(c['id'] in ann for c in calls)} of them, "
              f"{sum(c['genes'] != '.' for c in calls)} inside a gene")


if __name__ == "__main__":
    main()
