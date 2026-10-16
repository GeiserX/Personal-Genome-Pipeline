#!/usr/bin/env python3
"""paraphase_report.py: step 35b's two ends around Paraphase.

  paraphase_report.py check-bam BAM

reads the BAM's header and exits 0 only when every read group says PL:PACBIO.
Paraphase is written for PacBio HiFi reads (its own help: "BAM or CRAM file
with aligned HiFi reads"), and its licence (the PacBio Software License
Agreement, from Paraphase 3.3.0 on) allows it only on data generated on a
PacBio instrument. A short-read or Oxford Nanopore BAM stops here, before
Paraphase reads it, with exit 3 and the reason.

  paraphase_report.py table --sample S --json S.paraphase.json
                            --vcf-dir S_paraphase_vcfs
                            --out-cn S_paraphase_copy_number.tsv
                            --out-variants S_paraphase_variants.tsv

turns Paraphase's JSON into one row per gene (copy number, call status, what
the region-specific fields say) and its per-region VCFs into one row per
variant a haplotype carries, with the reads behind it. The genes of
GENE_ROWS get their own rows; every other region Paraphase ran gets one row:
the gene's own copies where Paraphase labels haplotypes gene1 (the gene) and
gene2 (its pseudogene), as for GBA and OTOA, else the gene family's copies.

parascopy_smn1() and compare_smn1() are the rule bin/collect_summary.py uses
to set SMN1 from Paraphase beside SMN1 from Parascopy (step 35), so both
reports judge agreement the same way.

Used by scripts/35b-paraphase.sh and the PARAPHASE_CHECK and PARAPHASE_REPORT
processes (bin/ is on the task PATH); tests/test_paraphase_report.py checks
it. Standard library only.
"""
import argparse
import gzip
import json
import os
import re
import struct
import sys
import zlib

CN_HEADER = ["region", "gene", "copy_number", "status", "detail", "total_cn", "region_depth", "genome_depth"]
VARIANT_HEADER = ["region", "chrom", "pos", "ref", "alt", "filter", "haplotypes", "support"]

# Region-specific copy-number fields, per Paraphase region (data/38/config.yaml
# names). A field that is null is a no-call: Paraphase saw one haplotype and
# the depth did not settle one copy or two.
GENE_ROWS = {
    "smn1": [("SMN1", "smn1_cn"), ("SMN2", "smn2_cn"), ("SMN exon 7-8 deletion", "smn_del78_cn")],
    "pms2": [("PMS2", "gene_cn")],
    "strc": [("STRC", "gene_cn")],
    "ncf1": [("NCF1", "gene_cn")],
}

# HBA structural variants Paraphase names from a hybrid haplotype (sv_called).
HBA_SV = {
    "3p7del": "3.7 kb deletion (-a3.7): HBA2 and HBA1 fused into one gene, one alpha-globin gene lost",
    "4p2del": "4.2 kb deletion (-a4.2): HBA2 lost",
    "3p7dup": "3.7 kb duplication (aaa anti-3.7): one extra alpha-globin gene",
    "4p2dup": "4.2 kb duplication (aaa anti-4.2): one extra HBA2",
}
HBA_LARGE = ("an allele with both genes lost (--) is a large deletion such as --SEA, --MED or --FIL; "
             "its breakpoints lie outside the region Paraphase reads, so which one it is cannot be told")

# SMN1 exon 7, c.840 (the base that tells SMN1 from SMN2), GRCh38 1-based.
SMN1_EXON7 = ("chr5", 70951946)


class Unreadable(Exception):
    pass


# --- check-bam -----------------------------------------------------------------

def bam_header(path):
    """The SAM header text of a BAM (BGZF is gzip with many members)."""
    try:
        with gzip.open(path, "rb") as f:
            if f.read(4) != b"BAM\x01":
                raise Unreadable(f"{path} is not a BAM file")
            (l_text,) = struct.unpack("<i", f.read(4))
            return f.read(l_text).decode("ascii", "replace")
    except (OSError, EOFError, struct.error, zlib.error) as e:
        raise Unreadable(f"{path}: {e}") from e


def platforms(header):
    """[(read group ID, PL)] of the header's @RG lines; PL '' when absent."""
    out = []
    for line in header.split("\n"):
        if not line.startswith("@RG"):
            continue
        tags = dict(t.split(":", 1) for t in line.split("\t")[1:] if ":" in t)
        out.append((tags.get("ID", "?"), tags.get("PL", "")))
    return out


def check_bam(path):
    """(ok, message): ok when the BAM has read groups and all are PL:PACBIO."""
    rgs = platforms(bam_header(path))
    if not rgs:
        return False, (f"{path} has no @RG line, so nothing says the reads are PacBio HiFi. Paraphase needs "
                       "a BAM of PacBio HiFi reads with PL:PACBIO in each read group (step 02b with "
                       "PLATFORM=hifi writes one, as pbmm2 does).")
    other = sorted({pl or "none" for _, pl in rgs if pl.upper() != "PACBIO"})
    if other:
        return False, (f"{path} has read groups with PL {', '.join(other)}, not PACBIO. Paraphase is written for "
                       "PacBio HiFi reads, and its licence (PacBio Software License Agreement) allows it only on "
                       "data from a PacBio instrument. For short reads, step 35 (Parascopy) estimates SMN1 and "
                       "SMN2 copy number.")
    return True, f"{path}: {len(rgs)} read group(s), all PL:PACBIO"


# --- table -----------------------------------------------------------------------

def num(v):
    """An int for an integer value, else the value; None stays None."""
    if isinstance(v, float) and v.is_integer():
        return int(v)
    return v


def text(v):
    return "." if v is None or v == "" else str(num(v))


def depth_text(region):
    d = region.get("region_depth")
    if isinstance(d, dict):
        d = d.get("median")
    elif isinstance(d, list) and d:
        d = d[0]
    return text(d)


def hap_names(region):
    """Names of the final haplotypes, each two-copy haplotype counted twice."""
    finals = region.get("final_haplotypes") or {}
    names = list(finals.values()) if isinstance(finals, dict) else list(finals)
    return names + [h for h in (region.get("two_copy_haplotypes") or []) if h in names]


def hba_rows(region):
    rsi = region.get("region_specific_info") or {}
    names = hap_names(region)
    hba1 = sum(1 for h in names if "hba1hap" in h)
    hba2 = sum(1 for h in names if "hba2hap" in h)
    genotype = rsi.get("genotype")
    svs = rsi.get("sv_called") or {}
    named = []
    for hap, val in sorted(svs.items()):
        kind = str(val).split(",")[0]
        named.append(HBA_SV.get(kind, kind) + f" ({hap})")
    detail = []
    if genotype:
        detail.append(f"genotype {genotype} ({genotype.count('a')} alpha-globin genes)")
        if "--" in genotype:
            detail.append(HBA_LARGE)
    else:
        detail.append("no genotype: Paraphase could not resolve the alleles")
    detail += named
    if not named and genotype:
        detail.append("no 3.7 kb or 4.2 kb deletion or duplication found")
    status = "called" if genotype else "no call"
    return [("HBA1", hba1 if names else None, status, "; ".join(detail)),
            ("HBA2", hba2 if names else None, status, "; ".join(detail))]


def rccx_rows(region):
    rsi = region.get("region_specific_info") or {}
    alleles = rsi.get("annotated_alleles") or []
    if isinstance(alleles, dict):
        alleles = list(alleles.values())
    detail = ["RCCX module copies (each carries CYP21A2 or its pseudogene CYP21A1P)"]
    if alleles:
        detail.append("CYP21A2 alleles: " + " / ".join(str(a) for a in alleles))
    if rsi.get("phasing_success") is False:
        detail.append("the copies could not be phased into two alleles")
    total = region.get("total_cn")
    return [("CYP21A2 (RCCX)", total, "called" if total is not None else "no call", "; ".join(detail))]


def generic_rows(name, region):
    rsi = region.get("region_specific_info") or {}
    genes = region.get("genes_in_region") or [name]
    if isinstance(genes, str):
        genes = genes.split(",")
    gene = ",".join(genes)
    detail = []
    fusions = rsi.get("fusions_called")
    if fusions:
        detail.append("fusions: " + json.dumps(fusions, sort_keys=True))
    names = hap_names(region)
    if "gene_cn" in rsi:
        cn = rsi.get("gene_cn")
        detail.insert(0, "copies of the gene, not of its pseudogenes")
    elif any("_gene1hap" in h or "_gene2hap" in h for h in names):
        # Paraphase labels each haplotype as the region's first gene (gene1) or
        # its second gene or pseudogene (gene2): GBA1/GBAP1, OTOA/OTOAP1. A
        # haplotype it cannot place is "unknownhap" and counts for neither.
        gene = genes[0]
        other = genes[1] if len(genes) > 1 else "pseudogene"
        cn = sum("_gene1hap" in h for h in names)
        unknown = sum("_unknownhap" in h for h in names)
        detail.insert(0, f"copies of {gene} alone (Paraphase's gene1 haplotypes); {other} copies "
                      f"{sum('_gene2hap' in h for h in names)} (gene2 haplotypes)"
                      + (f"; {unknown} haplotype(s) not placed on either" if unknown else ""))
    else:
        cn = region.get("total_cn")
        detail.insert(0, "copies of the whole gene family (gene and paralogs together)")
    return [(gene, cn, "called" if cn is not None else "no call", "; ".join(detail))]


def gene_rows(name, region):
    if name == "hba":
        return hba_rows(region)
    if name == "rccx":
        return rccx_rows(region)
    if name in GENE_ROWS:
        rsi = region.get("region_specific_info") or {}
        rows = []
        for gene, field in GENE_ROWS[name]:
            cn = rsi.get(field)
            detail = ("one haplotype seen and the depth does not settle one copy or two"
                      if cn is None and field in rsi else "")
            rows.append((gene, cn, "called" if cn is not None else "no call", detail))
        return rows
    return generic_rows(name, region)


def load_json(path):
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, ValueError) as e:
        raise Unreadable(f"{os.path.basename(path)}: {e}") from e
    if not isinstance(data, dict):
        raise Unreadable(f"{os.path.basename(path)} is not a JSON object of regions")
    return data


def cn_table(data):
    """Rows (CN_HEADER order) for every region of Paraphase's JSON."""
    out = []
    for name in sorted(data, key=lambda k: (k not in GENE_ROWS and k not in ("hba", "rccx"), k.lower())):
        region = data[name]
        if not isinstance(region, dict):
            continue
        for gene, cn, status, detail in gene_rows(name, region):
            out.append([name, gene, text(cn), status, detail or ".", text(region.get("total_cn")),
                        depth_text(region), text(region.get("genome_depth"))])
    return out


def variant_rows(region, path):
    """Rows (VARIANT_HEADER order) of one Paraphase VCF: every record with a
    haplotype that carries the ALT allele."""
    out, samples = [], []
    try:
        with open(path, errors="replace") as f:
            for line in f:
                if line.startswith("##"):
                    continue
                c = line.rstrip("\n").split("\t")
                if line.startswith("#CHROM"):
                    samples = c[9:]
                    continue
                if len(c) < 10:
                    continue
                keys = c[8].split(":")
                carriers, support = [], []
                for hap, val in zip(samples, c[9:]):
                    f_ = dict(zip(keys, val.split(":")))
                    if f_.get("GT") not in ("1", "1/1", "1|1"):
                        continue
                    carriers.append(hap)
                    ad = (f_.get("AD") or "").split(",")
                    alt_reads = ad[1] if len(ad) > 1 else "."
                    support.append(f"{hap}:{alt_reads}/{f_.get('DP', '.')}")
                if carriers:
                    out.append([region, c[0], c[1], c[3], c[4], c[6], ",".join(carriers), ",".join(support)])
    except OSError as e:
        raise Unreadable(f"{os.path.basename(path)}: {e}") from e
    return out


def variants_table(vcf_dir, sample):
    out = []
    if not vcf_dir or not os.path.isdir(vcf_dir):
        return out
    prefix = f"{sample}_"
    for name in sorted(os.listdir(vcf_dir)):
        if not name.endswith(".vcf"):
            continue
        region = name[len(prefix):-4] if name.startswith(prefix) else name[:-4]
        out += variant_rows(region, os.path.join(vcf_dir, name))
    return out


def write_tsv(path, header, rows):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp, path)


# --- SMN1 from both tools ------------------------------------------------------------

def parascopy_smn1(path):
    """SMN1 copy number from step 35's table (<s>_smn_copy_number.tsv): the
    paralog-specific value of the row, or of the homologous region, that holds
    SMN1 exon 7. {'cn': int or None, 'text': str}; cn is None unless the
    filter is PASS and the quality at least 20 (Parascopy's threshold)."""
    chrom, pos = SMN1_EXON7
    with open(path, errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            regions = []
            try:
                regions.append((r["chrom"], int(r["start"]) + 1, int(r["end"])))
            except (KeyError, ValueError):
                continue
            for h in (r.get("homologous_regions") or "*").split(","):
                m = re.match(r"^([^:]+):(\d+)-(\d+)", h.strip())
                if m:
                    regions.append((m.group(1), int(m.group(2)), int(m.group(3))))
            idx = next((i for i, (c, s, e) in enumerate(regions) if c == chrom and s <= pos <= e), None)
            if idx is None:
                continue
            cns = (r.get("psCN") or "").split(",")
            quals = (r.get("psCN_qual") or "").split(",")
            cn_txt = cns[idx] if idx < len(cns) else "?"
            q_txt = quals[idx] if idx < len(quals) else "*"
            filt = r.get("psCN_filter") or "."
            try:
                q = float(q_txt)
            except ValueError:
                q = None
            ok = filt == "PASS" and q is not None and q >= 20 and cn_txt.isdigit()
            return {"cn": int(cn_txt) if ok else None,
                    "text": f"{cn_txt} (quality {q_txt}, filter {filt})" + ("" if ok else ", not resolved")}
    return {"cn": None, "text": "no row covers SMN1 exon 7"}


def compare_smn1(paraphase_cn, parascopy):
    """'agree', 'disagree', or None when either tool has no SMN1 call."""
    if paraphase_cn is None or parascopy is None or parascopy.get("cn") is None:
        return None
    return "agree" if int(paraphase_cn) == parascopy["cn"] else "disagree"


# --- main -------------------------------------------------------------------------

def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check-bam", help="exit 0 only for a BAM whose read groups are all PL:PACBIO")
    c.add_argument("bam")
    t = sub.add_parser("table", help="Paraphase's JSON and VCFs as two tables")
    t.add_argument("--sample", required=True)
    t.add_argument("--json", required=True)
    t.add_argument("--vcf-dir")
    t.add_argument("--out-cn", required=True)
    t.add_argument("--out-variants", required=True)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "check-bam":
            ok, msg = check_bam(a.bam)
            print(msg if ok else f"ERROR: {msg}", file=sys.stdout if ok else sys.stderr)
            return 0 if ok else 3
        rows = cn_table(load_json(a.json))
        variants = variants_table(a.vcf_dir, a.sample)
    except Unreadable as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    if not rows:
        print(f"ERROR: {a.json} has no region", file=sys.stderr)
        return 1
    write_tsv(a.out_cn, CN_HEADER, rows)
    write_tsv(a.out_variants, VARIANT_HEADER, variants)
    for r in rows:
        if r[0] in GENE_ROWS or r[0] in ("hba", "rccx", "GBA", "OTOA"):
            print(f"  {r[1]}: {r[2]} ({r[3]})")
    print(f"  {len(rows)} gene rows, {len(variants)} variant rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
