#!/usr/bin/env python3
"""bloodtype_report.py: the red cell and platelet antigens one genome predicts.

  bloodtype_report.py --sample S --vcf S.vcf.gz --markers assets/bloodtype_markers.tsv
                      [--depth S_bloodtype_depth.tsv] [--cn-loci bin/cn_loci.tsv]
                      --out S_bloodtype.tsv --evidence S_bloodtype_evidence.tsv

reads the sample's genotype at each defining variant of assets/bloodtype_markers.tsv
and writes one row per blood group system: ABO (with A1 and A2), RhD, RhCE (C/c,
E/e), Kell (K/k, Kp, Js), Duffy (with the GATA-box variant), Kidd, MNS (M/N, S/s),
Lutheran, Diego, Dombrock, Colton, and the platelet antigens HPA-1 to -5 and -15.
Columns: system, phenotype, genotype, confidence, note. It is a prediction from
the DNA, not a serological typing. The evidence file has one row per site and per
RhD depth region: what the VCF and the reads said.

--depth is samtools bedcov -Q 20 -j over the table's sites and the RHD row of
bin/cn_loci.tsv (scripts/41-bloodtype.sh and the BLOODTYPE_DEPTH process make
it from the BAM; `copy_number.py bed --locus RHD` prints that row's regions,
and bin/copy_number.py turns their depth into the copy number). It gives:
  - RhD. The common RhD-negative allele (RHD*01N.01) deletes the whole RHD gene,
    which a VCF does not show. RHD depth over the depth of its two flanks gives
    the copy number (0, 1, 2); exons 3 to 7 without RHD reads in a gene that is
    there is the pattern of an RHD-RHCE hybrid, and RhD is then not determined.
    RHD*08N.01 (RHD pseudogene, c.487-19_504dup with c.609G>A) is an inactive
    copy and is subtracted. Without --depth RhD is not determined.
  - A check of the sites the VCF has no record for. A caller writes no record
    where the sample matches the reference, and none where it had no reads. With
    the depth a missing site counts as reference only with MIN_DEPTH reads at
    MAPQ >= 20, and is not determined below that; without it, a missing site is
    read as reference and the confidence says so (medium).

Confidence: high (every site called or covered by reads), medium (a site was
read as reference with no depth to confirm it), low (a site in a region where
short reads are unreliable, named in the note, or markers that disagree), not
determined. Phase is not read: two variants of one system are placed on the
chromosomes the common alleles carry them on, and the note says so.

Used by scripts/41-bloodtype.sh and the BLOODTYPE process (bin/ is on the task
PATH); tests/test_bloodtype_report.py checks it. Standard library only.
"""
import argparse
import csv
import gzip
import re
import sys

import copy_number as cn

MIN_DEPTH = 8          # reads (MAPQ >= 20) before a site missing from the VCF counts as reference
RHD_HYBRID_EXONS = (3, 4, 5, 6, 7)   # the exons an RHD-CE(3-7)-D hybrid takes from RHCE
TABLE_HEADER = ["system", "phenotype", "genotype", "confidence", "note"]
EVIDENCE_HEADER = ["system", "name", "variant", "rsid", "chrom", "pos", "ref", "alt",
                   "status", "alt_copies", "depth", "note"]
LEVELS = ["high", "medium", "low", "not determined"]
NOT_DETERMINED = "not determined"

# How each two-allele site is written: its label and the phenotype template, with
# {antigen} replaced by + or -. The antigen names are the table's ref_means and
# alt_means; tests/test_bloodtype_report.py checks that every site has a line.
PAIRS = {
    "RHCE_307": ("C/c", "C{C}c{c}"),
    "RHCE_676": ("E/e", "E{E}e{e}"),
    "KEL_578": ("K/k", "K{K}k{k}"),
    "KEL_841": ("Kp", "Kp(a{Kpa}b{Kpb})"),
    "KEL_1790": ("Js", "Js(a{Jsa}b{Jsb})"),
    "JK_838": ("Jk", "Jk(a{Jka}b{Jkb})"),
    "GYPA_59": ("M/N", "M{M}N{N}"),
    "GYPB_143": ("S/s", "S{S}s{s}"),
    "LU_230": ("Lu", "Lu(a{Lua}b{Lub})"),
    "DI_2561": ("Di", "Di(a{Dia}b{Dib})"),
    "DO_793": ("Do", "Do(a{Doa}b{Dob})"),
    "CO_134": ("Co", "Co(a{Coa}b{Cob})"),
    "HPA1_176": ("HPA-1", "HPA-1(a{1a}b{1b})"),
    "HPA2_482": ("HPA-2", "HPA-2(a{2a}b{2b})"),
    "HPA3_2621": ("HPA-3", "HPA-3(a{3a}b{3b})"),
    "HPA4_506": ("HPA-4", "HPA-4(a{4a}b{4b})"),
    "HPA5_1600": ("HPA-5", "HPA-5(a{5a}b{5b})"),
    "HPA15_2108": ("HPA-15", "HPA-15(a{15a}b{15b})"),
}
# The systems in report order. ABO, RhD and Duffy have rules of their own; the
# others are their two-allele sites.
SYSTEMS = [
    ("ABO", None), ("RhD", None), ("RhCE", ["RHCE_307", "RHCE_676"]),
    ("Kell", ["KEL_578", "KEL_841", "KEL_1790"]), ("Duffy", None), ("Kidd", ["JK_838"]),
    ("MNS", ["GYPA_59", "GYPB_143"]), ("Lutheran", ["LU_230"]), ("Diego", ["DI_2561"]),
    ("Dombrock", ["DO_793"]), ("Colton", ["CO_134"]),
    ("HPA-1", ["HPA1_176"]), ("HPA-2", ["HPA2_482"]), ("HPA-3", ["HPA3_2621"]),
    ("HPA-4", ["HPA4_506"]), ("HPA-5", ["HPA5_1600"]), ("HPA-15", ["HPA15_2108"]),
]
SPECIAL_SITES = {"ABO_261", "ABO_802", "ABO_803", "ABO_796", "ABO_1061",
                 "RHD_PSI_DUP", "RHD_PSI_609", "FY_125", "FY_GATA", "FY_265"}


class InputError(Exception):
    pass


# --- inputs ------------------------------------------------------------------

def norm_chrom(c):
    if c.startswith("chr"):
        return c
    return "chrM" if c in ("MT", "M") else "chr" + c


def rhd_regions(cn_loci=None):
    """The RhD depth regions (the RHD row of bin/cn_loci.tsv: flanks, gene,
    exons), by name, shaped as marker rows."""
    row = cn.load_catalog(cn_loci)["RHD"]
    return {n: {"chrom": row["chrom"], "start": s, "end": e, "kind": "region", "system": "RhD", "name": n,
                "hgvs": label, "ref": "", "alt": "", "rsid": ""}
            for n, s, e, label in row["_controls"] + row["_targets"] + row["_details"]}


def read_markers(path, cn_loci=None):
    """Sites of the marker table and the RhD regions of the catalog, by
    name, in table order."""
    sites, regions = {}, rhd_regions(cn_loci)
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("#") and ln.strip()]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if not rows or "kind" not in rows[0]:
        raise InputError(f"{path}: no marker rows (header with a 'kind' column expected)")
    for r in rows:
        r["start"], r["end"] = int(r["start"]), int(r["end"])
        if r["kind"] == "site":
            r["pos"] = r["end"]
            sites[r["name"]] = r
        else:
            raise InputError(f"{path}: row {r['name']} has kind '{r['kind']}' (site; the RhD depth regions "
                             "are the RHD row of bin/cn_loci.tsv)")
    return sites, regions


def open_text(path):
    with open(path, "rb") as f:
        magic = f.read(2)
    return gzip.open(path, "rt", errors="replace") if magic == b"\x1f\x8b" else open(path, errors="replace")


class Site:
    """What the VCF and the reads say at one site: status, alt copies (None when not
    known), read depth, note."""

    def __init__(self, row):
        self.row, self.status, self.n, self.depth, self.note = row, "absent", None, None, ""

    @property
    def caveat(self):
        return self.row.get("caveat") or ""


def read_vcf(path, sites):
    """Genotype the sites from the sample's VCF (the first sample column)."""
    by_pos = {}
    for name, r in sites.items():
        by_pos.setdefault((r["chrom"], r["pos"]), []).append(name)
    out = {name: Site(r) for name, r in sites.items()}
    nsamples = None
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                nsamples = len(line.rstrip("\n").split("\t")) - 9
                if nsamples < 1:
                    raise InputError(f"{path} has no sample column: a genotype is needed at every site")
                if nsamples > 1:
                    raise InputError(f"{path} has {nsamples} samples: this step reads one person's VCF")
                continue
            c = line.split("\t", 2)
            if len(c) < 3:
                continue
            try:
                key = (norm_chrom(c[0]), int(c[1]))
            except ValueError:
                continue
            if key not in by_pos:
                continue
            f9 = line.rstrip("\n").split("\t")
            if len(f9) < 10:
                continue
            ref, alts, filt = f9[3], f9[4].split(","), f9[6]
            fmt, smp = f9[8].split(":"), f9[9].split(":")
            gt = smp[fmt.index("GT")] if "GT" in fmt and fmt.index("GT") < len(smp) else "."
            alleles = re.split(r"[/|]", gt)
            for name in by_pos[key]:
                site = out[name]
                if site.status == "called":
                    continue   # a second record at the position: the first one counts
                genotype_site(site, ref, alts, filt, alleles)
    if nsamples is None:
        raise InputError(f"{path} has no #CHROM header line: not a VCF")
    return out


def genotype_site(site, ref, alts, filt, alleles):
    r = site.row
    want = r["alt"] if r["ref"] == ref else None
    idx = str(alts.index(want) + 1) if want in alts else None
    if "." in alleles or len(alleles) != 2:
        site.status, site.note = "no call", f"genotype {'/'.join(alleles)}"
        return
    # A record whose REF is not the marker's REF says nothing about the marker's ALT,
    # unless it is a reference call over a span that starts with the marker's REF.
    if ref != r["ref"] and not (ref.startswith(r["ref"]) and all(a == "0" for a in alleles)):
        site.status, site.note = "no call", f"REF {ref} at this position, the marker's REF is {r['ref']}"
        return
    if filt not in ("PASS", ".") and any(a != "0" for a in alleles):
        site.status, site.note = "no call", f"FILTER {filt}"
        return
    site.status, site.n, site.note = "called", sum(1 for a in alleles if idx is not None and a == idx), ""
    others = [a for a in alleles if a not in ("0", idx)]
    if others:
        site.note = f"another allele at this position ({ref}>{','.join(alts)}, genotype {'/'.join(alleles)})"


def read_depth(path):
    """samtools bedcov output over the marker sites (name in column 6) and the
    RHD regions (a 4-column BED, name in column 4): mean depth by name."""
    depth = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or c[0].startswith("#"):
                continue
            try:
                start, end, total = int(c[1]), int(c[2]), float(c[-1])
            except ValueError:
                continue
            depth[c[5] if len(c) >= 7 else c[3]] = total / max(end - start, 1)
    return depth


def apply_depth(sites, depth):
    """A site with no VCF record: reference when the reads cover it, else not known."""
    for s in sites.values():
        d = depth.get(s.row["name"]) if depth is not None else None
        s.depth = d
        if s.status != "absent":
            continue
        if depth is None:
            s.status, s.n = "assumed", 0
            s.note = "no VCF record; no BAM depth to confirm the reference"
        elif d is None:
            s.status, s.note = "no call", "no VCF record and no depth row (the BAM has no such contig)"
        elif d >= MIN_DEPTH:
            s.status, s.n = "covered", 0
        else:
            s.status, s.note = "low depth", f"no VCF record and {d:.0f} reads at MAPQ 20 (fewer than {MIN_DEPTH})"


# --- systems -----------------------------------------------------------------

def level_of(used):
    """Confidence from the sites a call used: low for a caveat or a stray allele,
    medium for a site read as reference without depth, else high."""
    lv = "high"
    for s in used:
        if s.caveat or s.note.startswith("another allele"):
            return "low"
        if s.status == "assumed":
            lv = "medium"
    return lv


def notes_of(used, extra=()):
    out = list(extra)
    for s in used:
        if s.caveat and s.caveat not in out:
            out.append(s.caveat)
        if s.note.startswith("another allele"):
            out.append(f"{s.row['hgvs']}: {s.note}")
    if any(s.status == "assumed" for s in used):
        out.append("sites with no VCF record read as reference (no BAM depth to confirm them)")
    return "; ".join(out)


def unknown(system, sites_missing, extra=""):
    why = "; ".join(f"{s.row['hgvs']}: {s.status}" + (f" ({s.note})" if s.note else "") for s in sites_missing)
    return {"system": system, "phenotype": NOT_DETERMINED, "genotype": ".", "confidence": NOT_DETERMINED,
            "note": "; ".join(x for x in (extra, why) if x)}


def sign(present):
    return "+" if present else "-"


def pair_part(name, s):
    """(phenotype part, genotype part) of one two-allele site, or None."""
    if s.n is None:
        return None
    r = s.row
    tmpl = PAIRS[name][1]
    pheno = tmpl.replace("{" + r["alt_means"] + "}", sign(s.n > 0)).replace("{" + r["ref_means"] + "}", sign(s.n < 2))
    alleles = sorted([r["alt_means"]] * s.n + [r["ref_means"]] * (2 - s.n))
    return pheno, "/".join(alleles)


def system_pairs(system, names, S):
    parts, genos, used, missing = [], [], [], []
    for name in names:
        s = S[name]
        p = pair_part(name, s)
        if p is None:
            parts.append(f"{PAIRS[name][0]} {NOT_DETERMINED}")
            missing.append(s)
            continue
        parts.append(p[0])
        genos.append(p[1])
        used.append(s)
    if not used:
        return unknown(system, missing)
    conf = level_of(used) if not missing else "low"
    extra = [f"{s.row['hgvs']}: {s.status}" + (f" ({s.note})" if s.note else "") for s in missing]
    return {"system": system, "phenotype": " ".join(parts), "genotype": " ".join(genos),
            "confidence": conf, "note": notes_of(used, extra)}


def system_abo(S):
    need = [S["ABO_261"], S["ABO_802"], S["ABO_803"]]
    missing = [s for s in need if s.n is None]
    if missing:
        return unknown("ABO", missing)
    n_func, n_o2, n_b = (s.n for s in need)
    non_o = n_func - n_o2
    if non_o < 0 or n_b > non_o:
        return unknown("ABO", [], f"the markers do not add up to two ABO alleles (c.261 insertion x{n_func}, "
                                  f"c.802G>A x{n_o2}, c.803G>C x{n_b})")
    n_a = non_o - n_b
    used, extra = list(need), []
    alleles = []
    a_label = ""
    if n_a:
        s1061 = S["ABO_1061"]
        if s1061.n is None:
            alleles += ["A"] * n_a
            a_label = "A"
            extra.append(f"A1 or A2 not determined ({s1061.row['hgvs']}: {s1061.status})")
        else:
            # A1 or A2 only: c.1061delC's caveat goes to the note, not to the
            # confidence of the group
            if s1061.caveat:
                extra.append(f"A1 or A2 rests on {s1061.row['hgvs']}, {s1061.caveat}")
            if s1061.status == "assumed":
                used.append(s1061)
            n_a2 = min(s1061.n, n_a)
            alleles += ["A1"] * (n_a - n_a2) + ["A2"] * n_a2
            a_label = "A1" if n_a > n_a2 else "A2"
    alleles += ["B"] * n_b + ["O.01"] * (2 - n_func) + ["O.02"] * n_o2
    s796 = S["ABO_796"]
    low = False
    if s796.n is not None and s796.n != n_b:
        low = True
        extra.append(f"c.796C>A (x{s796.n}) and c.803G>C (x{n_b}) disagree: a rare A or B allele")
    pheno = (a_label + ("B" if n_b else "")) or "O"
    conf = "low" if low else level_of(used)
    if n_b and n_a:
        extra.append("A and B on different chromosomes assumed (no phase)")
    return {"system": "ABO", "phenotype": pheno, "genotype": "/".join(alleles), "confidence": conf,
            "note": notes_of(used, extra)}


def system_duffy(S):
    s125, sg, s265 = S["FY_125"], S["FY_GATA"], S["FY_265"]
    missing = [s for s in (s125, sg) if s.n is None]
    if missing:
        return unknown("Duffy", missing)
    nb, ng = s125.n, sg.n
    sil_b = min(ng, nb)
    sil_a = ng - sil_b
    a_expr, b_expr = (2 - nb) - sil_a, nb - sil_b
    used, extra = [s125, sg], []
    weak = 0
    if s265.n is not None:
        used.append(s265)
        weak = min(s265.n, b_expr)
    b_sign = "+" if b_expr - weak > 0 else ("+w" if weak else "-")
    alleles = (["FY*A"] * (2 - nb - sil_a) + ["FY*A GATA"] * sil_a + ["FY*B"] * (b_expr - weak)
               + ["FY*X"] * weak + ["FY*B GATA"] * sil_b)
    if ng:
        extra.append("c.-67T>C (GATA box) stops Duffy expression on red cells only; with no phase it is "
                     "placed on the Fyb chromosome first, where it usually lies (FY*02N.01)")
    return {"system": "Duffy", "phenotype": f"Fy(a{sign(a_expr > 0)}b{b_sign})", "genotype": "/".join(alleles),
            "confidence": level_of(used), "note": notes_of(used, extra)}


def rhd_copies(depth, regions, cn_loci=None):
    """(copies of the gene, {exon: copies}, flank depth) from bedcov means, or None
    (bin/copy_number.py on the RHD row; cn.BrokenControl when a flank has no read
    while the other has the depth the row needs)."""
    b = cn.assess(cn.load_catalog(cn_loci)["RHD"], {n: depth[n] for n in regions if n in depth})
    if b["reason"] == "no_data":
        return None
    if b["reason"] == "low_control_depth":
        return None, {}, b["control_depth"]
    exons = {int(n[len("RHD_EX"):]): c for n, c in b["details"].items()}
    return b["estimate"], exons, b["control_depth"]


def system_rhd(S, depth, regions, cn_loci=None):
    if depth is None:
        return {"system": "RhD", "phenotype": NOT_DETERMINED, "genotype": ".", "confidence": NOT_DETERMINED,
                "note": "needs the BAM: the common RhD-negative allele (RHD*01N.01) deletes the whole RHD gene, "
                        "which a VCF does not show"}
    try:
        est = rhd_copies(depth, regions, cn_loci)
    except cn.BrokenControl as e:
        return {"system": "RhD", "phenotype": NOT_DETERMINED, "genotype": ".", "confidence": NOT_DETERMINED,
                "note": str(e)}
    if est is None:
        return {"system": "RhD", "phenotype": NOT_DETERMINED, "genotype": ".", "confidence": NOT_DETERMINED,
                "note": "no read depth over RHD and its flanks (the BAM has no "
                        f"{regions['RHD_GENE']['chrom']} contig)"}
    if est[0] is None:
        return {"system": "RhD", "phenotype": NOT_DETERMINED, "genotype": ".", "confidence": NOT_DETERMINED,
                "note": f"too few reads beside RHD ({est[2]:.1f}x) for a copy number (at least {rhd_min_flank(cn_loci)}x)"}
    copies, exons, _ = est
    call = cn.cut_copies(copies, cn.load_catalog(cn_loci)["RHD"]["_call"][1])
    sd, s609 = S["RHD_PSI_DUP"], S["RHD_PSI_609"]
    used, extra = [], []
    n_psi = 0
    if call == 0:
        pass   # no RHD left: the pseudogene markers have nothing to say
    elif sd.n is not None and s609.n is not None:
        used += [sd, s609]
        n_psi = min(sd.n, s609.n, call)
    else:
        extra.append("RHD pseudogene (RHD*08N.01) not checked: "
                     + "; ".join(f"{s.row['hgvs']}: {s.status}" for s in (sd, s609) if s.n is None))
    hybrid = [e for e in RHD_HYBRID_EXONS if call >= 1 and exons.get(e, 1.0) < 0.5]
    exon_txt = ", ".join(f"exon {e} {exons[e]:.2f}" for e in sorted(exons))
    geno = (["RHD"] * (call - n_psi) + ["RHD*08N.01"] * n_psi + ["RHD deleted"] * (2 - call))
    basis = f"{'/'.join(geno)} (RHD copies by read depth {copies:.2f})"
    if hybrid:
        return {"system": "RhD", "phenotype": NOT_DETERMINED, "genotype": basis, "confidence": NOT_DETERMINED,
                "note": f"RHD exon{'s' if len(hybrid) > 1 else ''} {', '.join(map(str, hybrid))} without RHD reads "
                        f"while the gene has {copies:.2f} copies: the pattern of an RHD-RHCE hybrid gene "
                        f"(such as RHD*03N.01) or a partial deletion. Copies per exon: {exon_txt}"}
    margin = min(abs(copies - 0.5), abs(copies - 1.5))
    conf = level_of(used) if used else ("high" if call == 0 else "medium")
    if margin < 0.2:
        conf = "low"
        extra.append(f"read depth {copies:.2f} copies is near a copy-number boundary")
    if copies > 2.5:
        extra.append(f"read depth gives {copies:.2f} RHD copies, more than two")
    extra.append("weak and partial D (RHD missense alleles) are not assessed")
    return {"system": "RhD", "phenotype": "D+" if call - n_psi > 0 else "D-", "genotype": basis,
            "confidence": conf, "note": notes_of(used, extra)}


def rhd_min_flank(cn_loci=None):
    """The RHD row's minimum flank depth, written as the old constant was (10)."""
    v = cn.load_catalog(cn_loci)["RHD"]["_min_control"]
    return int(v) if v == int(v) else v


def predict(S, depth, regions, cn_loci=None):
    rows = []
    for system, names in SYSTEMS:
        if system == "ABO":
            rows.append(system_abo(S))
        elif system == "RhD":
            rows.append(system_rhd(S, depth, regions, cn_loci))
        elif system == "Duffy":
            rows.append(system_duffy(S))
        else:
            rows.append(system_pairs(system, names, S))
    return rows


# --- output ------------------------------------------------------------------

def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(h, "") for h in header])


def evidence_rows(S, depth, regions, cn_loci=None):
    rows = []
    for s in S.values():
        r = s.row
        rows.append({"system": r["system"], "name": r["name"], "variant": r["hgvs"], "rsid": r["rsid"],
                     "chrom": r["chrom"], "pos": r["pos"], "ref": r["ref"], "alt": r["alt"], "status": s.status,
                     "alt_copies": "." if s.n is None else s.n,
                     "depth": "." if s.depth is None else f"{s.depth:.1f}", "note": s.note})
    try:
        est = rhd_copies(depth, regions, cn_loci) if depth is not None else None
    except cn.BrokenControl:
        est = None
    flank = est[2] if est else None
    for name, r in regions.items():
        d = (depth or {}).get(name)
        cpy = f"{2 * d / flank:.2f}" if d is not None and flank and not name.startswith("RHD_FLANK") else "."
        rows.append({"system": r["system"], "name": name, "variant": r["hgvs"], "chrom": r["chrom"],
                     "pos": f"{r['start'] + 1}-{r['end']}", "status": "depth" if d is not None else "no depth",
                     "alt_copies": cpy, "depth": "." if d is None else f"{d:.1f}",
                     "note": "copies = 2 x depth / flank depth" if cpy != "." else ""})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--depth", help="samtools bedcov -Q 20 -j over the marker sites and the RHD regions, from the BAM")
    ap.add_argument("--cn-loci", help="the copy-number catalog (default: bin/cn_loci.tsv beside this script)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--evidence", required=True)
    a = ap.parse_args(argv)
    try:
        sites, regions = read_markers(a.markers, a.cn_loci)
        S = read_vcf(a.vcf, sites)
        depth = read_depth(a.depth) if a.depth else None
    except (InputError, OSError, cn.CatalogError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    apply_depth(S, depth)
    rows = predict(S, depth, regions, a.cn_loci)
    write_tsv(a.out, TABLE_HEADER, rows)
    write_tsv(a.evidence, EVIDENCE_HEADER, evidence_rows(S, depth, regions, a.cn_loci))
    print(f"Blood groups of {a.sample}, predicted from the genome (not a serological typing):")
    for r in rows:
        print(f"  {r['system']:<9} {r['phenotype']:<34} {r['confidence']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
