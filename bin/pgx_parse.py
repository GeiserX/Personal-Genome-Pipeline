#!/usr/bin/env python3
"""pgx_parse.py: the one PharmCAT report.json reader of the pipeline.

Used by scripts/27-cpic-lookup.sh, the CPIC_LOOKUP Nextflow module (bin/ is on
the task PATH), bin/collect_summary.py and tests/test_cpic_parser.py. Standard
library only, so it runs in the plain python image.

  pgx_parse.py cpic-report --sample S --report R.json --outdir DIR
                           [--pypgx S_pypgx_summary.tsv --comparison OUT.tsv]

writes DIR/S_phenotypes.tsv and DIR/S_cpic_recommendations.txt, and with
--pypgx the PharmCAT/pypgx comparison TSV. It exits 1 when the report cannot be
read or yields no gene: a report that parses to nothing is a format change,
never an all-clear.

PharmCAT 3.x writes `genes` as a flat {gene: GeneReport} map and `drugs` as
{guidance source: {drug: DrugReport}} holding only the annotations that match
the sample. 2.x nested `genes` by source, and older versions used a list; all
three shapes are read.
"""
import argparse
import csv
import json
import os
import re
import sys
from datetime import date

# Fallback only: used for a gene when the report itself names no drug for it
# (a report without a `drugs` section or `relatedDrugs`). Source:
# https://cpicpgx.org/guidelines/
STATIC_DRUGS = {
    "ABCG2": "rosuvastatin",
    "CACNA1S": "desflurane,enflurane,halothane,isoflurane,methoxyflurane,sevoflurane,succinylcholine",
    "CYP2B6": "efavirenz,sertraline",
    "CYP2C19": "clopidogrel,voriconazole,escitalopram,citalopram,sertraline,amitriptyline,clomipramine,doxepin,imipramine,trimipramine,lansoprazole,omeprazole,pantoprazole,dexlansoprazole",
    "CYP2C9": "warfarin,phenytoin,flurbiprofen,celecoxib,ibuprofen,lornoxicam,meloxicam,piroxicam,tenoxicam,siponimod",
    "CYP2D6": "codeine,tramadol,hydrocodone,oxycodone,amitriptyline,clomipramine,desipramine,doxepin,imipramine,nortriptyline,trimipramine,fluvoxamine,paroxetine,atomoxetine,ondansetron,tropisetron,tamoxifen,eliglustat",
    "CYP3A5": "tacrolimus",
    "CYP4F2": "warfarin",
    "DPYD": "fluorouracil,capecitabine,tegafur",
    "G6PD": "rasburicase,dapsone,chloroquine,primaquine,nitrofurantoin,methylene-blue",
    "HLA-A": "carbamazepine,oxcarbazepine",
    "HLA-B": "abacavir,carbamazepine,oxcarbazepine,phenytoin,allopurinol",
    "IFNL3": "peginterferon-alfa-2a,peginterferon-alfa-2b",
    "MT-RNR1": "aminoglycosides",
    "NUDT15": "azathioprine,mercaptopurine,thioguanine",
    "RYR1": "desflurane,enflurane,halothane,isoflurane,methoxyflurane,sevoflurane,succinylcholine",
    "SLCO1B1": "simvastatin,atorvastatin,rosuvastatin,pravastatin,pitavastatin,fluvastatin,lovastatin",
    "TPMT": "azathioprine,mercaptopurine,thioguanine",
    "UGT1A1": "atazanavir,belinostat,irinotecan",
    "VKORC1": "warfarin",
}

# Short names for PharmCAT's guidance sources, in the order they are listed.
SOURCES = [
    ("CPIC", "CPIC"),
    ("DPWG", "DPWG"),
    ("FDA Label", "FDA label"),
    ("FDA PGx", "FDA association"),
]

NO_CALL_NAMES = {"", "unknown", "none", "?", "n/a"}
NORMAL_WORDS = {"normal", "typical", "extensive"}


class GeneCall:
    """One gene of the report: what PharmCAT called and how to read it.

    labels are all the diplotypes PharmCAT lists for the gene. More than one
    means the data cannot tell them apart (positions missing from the VCF, for
    example); when their phenotypes differ the gene is 'ambiguous' and its
    first diplotype must not be read as the result."""

    def __init__(self, gene, diplotype, phenotype, called, labels=None, phenotypes=None):
        self.gene = gene
        self.diplotype = diplotype
        self.phenotype = phenotype
        self.called = called
        self.labels = labels or [diplotype]
        self.phenotypes = phenotypes or [phenotype]

    @property
    def ambiguous(self):
        return self.called and len(set(self.phenotypes)) > 1

    @property
    def status(self):
        """'not called', 'ambiguous', 'normal' or 'non-normal' (listed with its drugs)."""
        if not self.called:
            return "not called"
        if self.ambiguous:
            return "ambiguous"
        if is_normal(self.phenotype):
            return "normal"
        return "non-normal"

    def row(self):
        return (self.gene, self.diplotype, self.phenotype, self.status)


def is_normal(phenotype):
    """Normal function, or an HLA result whose every allele test is negative."""
    parts = [p.strip().lower() for p in phenotype.split(";") if p.strip()]
    if not parts:
        return False
    if all(p.endswith(" negative") for p in parts):
        return True
    return all(NORMAL_WORDS & set(p.replace("-", " ").split()) for p in parts)


def _one_diplotype(dip):
    """(label, phenotype, called) of one diplotype object."""
    alleles = [dip.get(a) for a in ("allele1", "allele2")]
    names = [(a or {}).get("name", "") for a in alleles if a]
    label = dip.get("label") or "/".join(n or "?" for n in names) or "?"
    if names:
        called = any((n or "").strip().lower() not in NO_CALL_NAMES for n in names)
    else:
        called = not label.lower().startswith("unknown")
    phenos = [p for p in (dip.get("phenotypes") or []) if p]
    if phenos:
        phenotype = "; ".join(phenos)
    elif dip.get("phenotype"):
        phenotype = dip["phenotype"]
    else:
        phenotype = "no phenotype assigned" if called else "No Result"
    if phenotype.strip().lower() in ("no result", "n/a") and not called:
        phenotype = "No Result"
    return label, phenotype, called


def _diplotype(name, g):
    """GeneCall from a gene object with a diplotype array, else None."""
    if not isinstance(g, dict):
        return None
    dips = [d for d in (g.get("sourceDiplotypes") or g.get("recommendationDiplotypes") or []) if isinstance(d, dict)]
    if not dips:
        return None
    parsed = [_one_diplotype(d) for d in dips]
    label, phenotype, called = parsed[0]
    labels = [p[0] for p in parsed]
    phenotypes = list(dict.fromkeys(p[1] for p in parsed))
    if len(parsed) > 1:
        label = f"{label} (1 of {len(parsed)} possible)"
        if len(phenotypes) > 1:
            shown = ", ".join(phenotypes[:5]) + (", ..." if len(phenotypes) > 5 else "")
            phenotype = f"ambiguous: one of {shown}"
    return GeneCall(name, label, phenotype, called, labels, phenotypes)


def parse_genes(data):
    """(status, [GeneCall]) with status OK, PARSE_EMPTY or UNKNOWN_FORMAT."""
    if not isinstance(data, dict):
        return "UNKNOWN_FORMAT", []
    calls = []
    genes = data.get("genes")
    if isinstance(genes, dict):
        for key, val in genes.items():
            if not isinstance(val, dict):
                continue
            if "sourceDiplotypes" in val or "recommendationDiplotypes" in val:
                r = _diplotype(key, val)                 # flat (3.x): key is the gene
                if r:
                    calls.append(r)
            else:
                for gene_name, g in val.items():         # nested (2.x): key is the source
                    r = _diplotype(gene_name, g)
                    if r:
                        calls.append(r)
    elif isinstance(genes, list):
        for entry in genes:
            if not isinstance(entry, dict):
                continue
            name = entry.get("geneSymbol", entry.get("gene", "Unknown"))
            r = _diplotype(name, entry)
            if r:
                calls.append(r)
            else:
                dl = entry.get("diplotype", "N/A")
                ph = entry.get("phenotype", "N/A")
                if dl != "N/A" or ph != "N/A":
                    calls.append(GeneCall(name, dl, ph, dl not in ("N/A", "Unknown/Unknown")))
    elif isinstance(data.get("geneResults"), list):
        for gr in data["geneResults"]:
            dl = gr.get("diplotype", "N/A")
            calls.append(GeneCall(gr.get("gene", "Unknown"), dl, gr.get("phenotype", "N/A"),
                                  dl not in ("N/A", "Unknown/Unknown")))
    else:
        return "UNKNOWN_FORMAT", []
    seen, out = set(), []
    for c in calls:
        if c.gene in seen:
            continue
        seen.add(c.gene)
        out.append(c)
    if not out:
        return "PARSE_EMPTY", []
    return "OK", out


def _source_name(key):
    for prefix, short in SOURCES:
        if str(key).startswith(prefix):
            return short
    return str(key)


def _source_rank(short):
    names = [s for _, s in SOURCES]
    return names.index(short) if short in names else len(names)


def drug_guidance(data):
    """{gene: [(drug, source, classification, recommendation, labels, phenotype)]}
    from the report's own `drugs` section. PharmCAT lists an annotation for
    every diplotype the sample may have; labels and phenotype say which
    diplotypes of the gene each one is for (see matching())."""
    out = {}
    drugs = data.get("drugs") if isinstance(data, dict) else None
    if not isinstance(drugs, dict):
        return out
    for src_key, by_drug in drugs.items():
        if not isinstance(by_drug, dict):
            continue
        src = _source_name(src_key)
        for drug_name, rep in by_drug.items():
            if not isinstance(rep, dict):
                continue
            for gl in rep.get("guidelines") or []:
                for ann in (gl or {}).get("annotations") or []:
                    if not isinstance(ann, dict):
                        continue
                    phen = ann.get("phenotypes") or {}
                    labels = {}
                    for gt in ann.get("genotypes") or []:
                        for dp in (gt or {}).get("diplotypes") or []:
                            if isinstance(dp, dict) and dp.get("gene"):
                                labels.setdefault(dp["gene"], set()).add(dp.get("label") or "")
                    rec = " ".join(str(ann.get("drugRecommendation") or "").split())
                    cls = ann.get("classification") or ""
                    for g in set(phen) | set(labels):
                        entry = (rep.get("name") or drug_name, src, cls, rec,
                                 frozenset(labels.get(g, ())), str(phen.get(g) or ""))
                        if entry not in out.setdefault(g, []):
                            out[g].append(entry)
    for g in out:
        out[g].sort(key=lambda e: (_source_rank(e[1]), e[0].lower()))
    return out


def related_drugs(data):
    """{gene: [drug]} from each gene's `relatedDrugs` (PharmCAT's gene-drug links)."""
    out = {}
    genes = data.get("genes") if isinstance(data, dict) else None
    if not isinstance(genes, dict):
        return out
    for name, g in genes.items():
        if isinstance(g, dict) and isinstance(g.get("relatedDrugs"), list):
            names = [d.get("name") for d in g["relatedDrugs"] if isinstance(d, dict) and d.get("name")]
            if names:
                out[name] = names
    return out


def matching(call, entries):
    """The guidance entries for the diplotype PharmCAT called: those naming the
    called diplotype, else those for the called phenotype."""
    by_label = [e for e in entries if call.labels[0] in e[4]]
    if by_label:
        return by_label
    return [e for e in entries if not e[4] and e[5] and e[5] in call.phenotype]


def drug_lines(call, guidance, related):
    """The lines that list a gene's medications, best source first: CPIC's
    matched recommendation per drug, then the drugs other sources name."""
    gene = call.gene
    entries = matching(call, guidance.get(gene) or [])
    if entries:
        lines = ["    Drugs with guidance in this PharmCAT report:"]
        shown = set()
        for drug, src, cls, rec, _, _ in entries:
            if src != "CPIC" or drug in shown:
                continue
            shown.add(drug)
            rec = re.sub(r"<[^>]+>", " ", rec)
            rec = " ".join(rec.split())
            head = f"      - {drug} [CPIC{', ' + cls if cls and cls != 'Unspecified' else ''}]"
            if rec:
                head += ": " + (rec if len(rec) <= 200 else rec[:197] + "...")
            lines.append(head)
        others = sorted({e[0] for e in entries if e[1] != "CPIC" and e[0] not in shown})
        if others:
            srcs = sorted({e[1] for e in entries if e[1] != "CPIC"}, key=_source_rank)
            lines.append(f"      also named by {', '.join(srcs)}: " + ", ".join(others))
        return lines
    if related.get(gene):
        return [f"    Drugs PharmCAT links to {gene} (no matched guidance text in the report): "
                + ", ".join(related[gene])]
    if STATIC_DRUGS.get(gene):
        return ["    Drugs (pipeline fallback table; the report names none): "
                + STATIC_DRUGS[gene].replace(",", ", ")]
    return [f"    NOTE: {gene} has a non-normal phenotype but is not in the drug table; "
            "see the PharmCAT HTML report."]


# --- pypgx ---------------------------------------------------------------------

def read_pypgx(path):
    """{gene: (diplotype, phenotype)} from step 32's / PYPGX's summary TSV."""
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene = (row.get("Gene") or "").strip()
            if gene:
                out[gene] = ((row.get("Diplotype") or "").strip(), (row.get("Phenotype") or "").strip())
    return out


def pypgx_called(diplotype):
    return diplotype not in ("", "FAILED", "N/A", "Not called")


def _norm(diplotype):
    # PharmCAT names VKORC1 alleles 'rs9923231 variant (T)' where pypgx writes
    # 'rs9923231', and either tool may put the alleles in either order.
    return sorted(re.sub(r" (variant|reference) \([ACGT]+\)", "", a).strip() for a in diplotype.split("/"))


def compare(calls, pypgx):
    """Rows (gene, PharmCAT, pypgx, match, called_by) for genes either tool called."""
    pc = {c.gene: c for c in calls}
    rows = []
    for gene in sorted(set(pc) | set(pypgx)):
        c = pc.get(gene)
        pc_dip = c.diplotype if c else "Not called"
        pc_ok = bool(c and c.called and not c.ambiguous)
        if c and c.ambiguous:
            pc_dip = f"ambiguous ({len(c.labels)} possible diplotypes)"
        pg_dip = pypgx.get(gene, ("Not called", ""))[0] or "Not called"
        pg_ok = pypgx_called(pg_dip)
        if not pc_ok and not pg_ok:
            continue
        if not pc_ok:
            match = by = "pypgx only"
        elif not pg_ok:
            match = by = "PharmCAT only"
        elif _norm(pc_dip) == _norm(pg_dip):
            match, by = "Yes", "both"
        else:
            match, by = "No", "both"
        rows.append((gene, pc_dip, pg_dip, match, by))
    return rows


def pypgx_warnings(calls, pypgx, guidance, related):
    """Lines for genes PharmCAT reports on but could not call while pypgx did:
    the medications section is silent on them, so say so."""
    lines = []
    for c in calls:
        if (c.called and not c.ambiguous) or c.gene not in pypgx:
            continue
        dip, phen = pypgx[c.gene]
        if not pypgx_called(dip):
            continue
        what = "no result" if not c.called else f"no single result ({len(c.labels)} possible diplotypes)"
        if is_normal(phen):
            # pypgx's call is normal: worth knowing, nothing to review.
            lines.append(f"  NOTE: PharmCAT has {what} for {c.gene}; pypgx (step 32) called {dip} ({phen}).")
            continue
        lines.append(f"  WARNING: PharmCAT has {what} for {c.gene}, but pypgx (step 32) called "
                     f"{dip} ({phen or 'no phenotype'}). The medications above do not cover {c.gene}.")
        # Every drug PharmCAT links to the gene, not only the ones whose
        # guidance happened to match another gene of this sample.
        names = related.get(c.gene) \
            or (STATIC_DRUGS.get(c.gene, "").split(",") if STATIC_DRUGS.get(c.gene) else []) \
            or [e[0] for e in guidance.get(c.gene, [])]
        if names:
            lines.append(f"    Drugs affected by {c.gene}: " + ", ".join(dict.fromkeys(names)))
        lines.append("    Review them with the pypgx call and docs/32-pypgx.md.")
    return lines


# --- report --------------------------------------------------------------------

def write_failure(outdir, sample, why):
    rec = os.path.join(outdir, f"{sample}_cpic_recommendations.txt")
    with open(rec, "w") as out:
        out.write("Pharmacogenomic Drug Recommendations\n=====================================\n")
        out.write(f"Sample: {sample}\nDate: {date.today().isoformat()}\n\n")
        out.write("!! PARSING FAILED -- RESULTS CANNOT BE TRUSTED !!\n")
        out.write(f"{why}\n")
        out.write("This is NOT a clean result. Read the PharmCAT HTML report directly,\n")
        out.write("and report this as a pipeline bug (likely a PharmCAT output-format change).\n")
    with open(os.path.join(outdir, f"{sample}_phenotypes.tsv"), "w") as f:
        f.write("Gene\tDiplotype\tPhenotype\tStatus\n")


def cpic_report(args):
    os.makedirs(args.outdir, exist_ok=True)
    try:
        with open(args.report) as f:
            data = json.load(f)
    except (OSError, ValueError) as e:
        write_failure(args.outdir, args.sample, f"Could not read {args.report}: {e}")
        print(f"ERROR: could not read PharmCAT report {args.report}: {e}", file=sys.stderr)
        return 1
    status, calls = parse_genes(data)
    if status != "OK":
        why = ("The report has no gene results in any known layout."
               if status == "UNKNOWN_FORMAT" else
               "The report's gene section was recognised but no gene diplotype could be read.")
        write_failure(args.outdir, args.sample, why)
        print(f"ERROR: {status}: no gene parsed from {args.report}; refusing to write an all-clear report",
              file=sys.stderr)
        return 1
    guidance = drug_guidance(data)
    related = related_drugs(data)
    pypgx = {}
    if args.pypgx:
        pypgx = read_pypgx(args.pypgx)

    with open(os.path.join(args.outdir, f"{args.sample}_phenotypes.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene", "Diplotype", "Phenotype", "Status"])
        w.writerows(c.row() for c in calls)

    lines = [
        "Pharmacogenomic Drug Recommendations",
        "=====================================",
        f"Sample: {args.sample}",
        f"Date: {date.today().isoformat()}",
        f"Source: PharmCAT {data.get('pharmcatVersion', '')} report ({os.path.basename(args.report)})".rstrip(),
        "",
        "Gene Results:",
        "-" * 72,
        f"{'Gene':<12} {'Diplotype':<30} Phenotype",
    ]
    lines += [f"{c.gene:<12} {c.diplotype:<30} {c.phenotype}" for c in calls]
    lines += ["", "Affected Medications:", "-" * 72, ""]
    listed = [c for c in calls if c.status == "non-normal"]
    for c in listed:
        lines.append(f"  {c.gene} -- {c.phenotype}:")
        lines.append(f"    Diplotype: {c.diplotype}")
        lines += drug_lines(c, guidance, related)
        lines.append("    Action: Consult CPIC guidelines at https://cpicpgx.org/guidelines/")
        lines.append("")
    if not listed:
        lines += ["  No gene with a non-normal phenotype.", ""]

    ambiguous = [c for c in calls if c.status == "ambiguous"]
    if ambiguous:
        lines += ["Genes With More Than One Possible Result:", "-" * 72, ""]
        for c in ambiguous:
            lines.append(f"  {c.gene} -- {len(c.labels)} possible diplotypes with different phenotypes "
                         f"({', '.join(c.phenotypes[:5])}{', ...' if len(c.phenotypes) > 5 else ''}).")
        lines += ["  The data cannot tell them apart (often positions missing from the VCF), so no",
                  "  drug guidance is given for them here. See the PharmCAT HTML report.", ""]

    lines += ["Uncallable Genes:", "-" * 72, ""]
    uncalled = [c for c in calls if c.status == "not called"]
    for c in uncalled:
        lines.append(f"  {c.gene} -- {c.phenotype} (not callable from available data)")
    if not uncalled:
        lines.append("  None -- all genes were successfully called.")
    lines.append("")

    if args.pypgx:
        lines += ["PharmCAT and pypgx:", "-" * 72, ""]
        warn = pypgx_warnings(calls, pypgx, guidance, related)
        lines += warn if warn else ["  No gene that PharmCAT could not call has a pypgx call."]
        lines.append("")

    lines += [
        "NOTE: Uncallable genes may lack coverage, have complex structural",
        "variants, or require data not present in your VCF. Their absence",
        "from the recommendations section does NOT mean normal function.",
        "",
        "-" * 72,
        "DISCLAIMER: These recommendations are based on CPIC clinical",
        "guidelines. Always consult your healthcare provider before",
        "making any medication changes. This is NOT medical advice.",
        "-" * 72,
    ]
    with open(os.path.join(args.outdir, f"{args.sample}_cpic_recommendations.txt"), "w") as out:
        out.write("\n".join(lines) + "\n")

    if args.pypgx and args.comparison:
        rows = compare(calls, pypgx)
        with open(args.comparison, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["Gene", "PharmCAT_diplotype", "pypgx_diplotype", "Match", "Called_by"])
            w.writerows(rows)
        same = sum(1 for r in rows if r[3] == "Yes")
        print(f"Comparison written: {args.comparison} (concordant {same}, other {len(rows) - same})")

    print(f"Genes parsed: {len(calls)}; non-normal: {len(listed)}; ambiguous: {len(ambiguous)}; "
          f"not called: {len(uncalled)}")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cpic-report", help="phenotypes TSV, recommendations text, pypgx comparison")
    c.add_argument("--sample", required=True)
    c.add_argument("--report", required=True, help="PharmCAT report.json")
    c.add_argument("--outdir", required=True)
    c.add_argument("--pypgx", help="pypgx summary TSV (step 32 / PYPGX)")
    c.add_argument("--comparison", help="where to write the PharmCAT/pypgx comparison TSV")
    args = p.parse_args(argv)
    if args.cmd == "cpic-report":
        return cpic_report(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
