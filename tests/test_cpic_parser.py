#!/usr/bin/env python3
"""Regression test for the PharmCAT->CPIC parser logic shared by
scripts/27-cpic-lookup.sh and modules/local/cpic_lookup/main.nf.

Guards the two clinical-safety properties:
  1. The parser auto-detects BOTH the flat 3.x `genes -> {gene -> data}` shape and
     the nested 2.15.x `genes -> {source -> {gene -> data}}` shape, so the script
     and the Nextflow module agree no matter which PharmCAT emits.
  2. A recognized PharmCAT report that yields zero genes FAILS LOUD (PARSE_EMPTY)
     instead of reporting a clean "all genes were successfully called" — a silent
     false-clear is the worst failure mode for a pharmacogenomics report.

This mirrors the embedded parser logic. A fixture test against a real PharmCAT
3.2.0 report.json (captured from an end-to-end run) is tracked as follow-up work.
Run: python3 tests/test_cpic_parser.py
"""
import sys


def parse_gene(name, g):
    if not isinstance(g, dict):
        return None
    dips = g.get('sourceDiplotypes') or g.get('recommendationDiplotypes') or []
    if not dips:
        return None
    dip = dips[0]
    a1 = (dip.get('allele1') or {}).get('name', '?')
    a2 = (dip.get('allele2') or {}).get('name', '?')
    diplotype = dip.get('label') or (a1 + '/' + a2)
    phenos = dip.get('phenotypes') or []
    phenotype = phenos[0] if phenos else dip.get('phenotype', 'N/A')
    return (name, diplotype, phenotype)


def extract(data):
    """Returns (status, results). status is OK | PARSE_EMPTY | UNKNOWN_FORMAT."""
    results = []
    genes = data.get('genes')
    if isinstance(genes, dict):
        for key, val in genes.items():
            if not isinstance(val, dict):
                continue
            if 'sourceDiplotypes' in val or 'recommendationDiplotypes' in val:
                r = parse_gene(key, val)            # flat: key is the gene
                if r:
                    results.append(r)
            else:
                for gene_name, g in val.items():    # nested: key is the source
                    r = parse_gene(gene_name, g)
                    if r:
                        results.append(r)
    elif isinstance(genes, list):
        for entry in genes:
            name = entry.get('geneSymbol', entry.get('gene', 'Unknown'))
            r = parse_gene(name, entry)
            if r:
                results.append(r)
            else:
                dl = entry.get('diplotype', 'N/A')
                ph = entry.get('phenotype', 'N/A')
                if dl != 'N/A' or ph != 'N/A':
                    results.append((name, dl, ph))
    elif isinstance(data.get('geneResults'), list):
        for gr in data['geneResults']:
            results.append((gr.get('gene', 'Unknown'),
                            gr.get('diplotype', 'N/A'),
                            gr.get('phenotype', 'N/A')))
    else:
        return ('UNKNOWN_FORMAT', [])
    seen, dedup = set(), []
    for g, d, p in results:
        if g in seen:
            continue
        seen.add(g)
        dedup.append((g, d, p))
    if not dedup:
        return ('PARSE_EMPTY', [])
    return ('OK', dedup)


def main():
    dip = {'allele1': {'name': '*1'}, 'allele2': {'name': '*2'},
           'phenotypes': ['Rapid Metabolizer']}
    cases = [
        ('flat 3.x',
         {'genes': {'CYP2C19': {'sourceDiplotypes': [dip]},
                    'CYP2D6': {'sourceDiplotypes': [dip]}}}, 'OK', 2),
        ('nested 2.15.x',
         {'genes': {'CPIC': {'CYP2C19': {'sourceDiplotypes': [dip]},
                             'CYP2D6': {'sourceDiplotypes': [dip]}}}}, 'OK', 2),
        ('recommendationDiplotypes',
         {'genes': {'DPYD': {'recommendationDiplotypes':
                             [{'label': 'c.1905+1G>A/Reference',
                               'phenotypes': ['Intermediate Metabolizer']}]}}}, 'OK', 1),
        ('list format',
         {'genes': [{'gene': 'TPMT', 'sourceDiplotypes': [dip]}]}, 'OK', 1),
        ('recognized-but-zero (the false-clear trap)',
         {'genes': {'CYP2C19': {'name': 'CYP2C19', 'someField': 1}}}, 'PARSE_EMPTY', 0),
        ('no gene keys at all',
         {'metadata': {'foo': 'bar'}}, 'UNKNOWN_FORMAT', 0),
    ]
    failures = 0
    for label, data, want_status, want_n in cases:
        status, res = extract(data)
        ok = status == want_status and len(res) == want_n
        print(f"[{'PASS' if ok else 'FAIL'}] {label}: status={status} n={len(res)} "
              f"(want {want_status}/{want_n})")
        if not ok:
            failures += 1
    # The single most important invariant: a zero-gene extraction is NEVER 'OK'.
    assert extract({'genes': {'X': {'y': 1}}})[0] != 'OK', "zero genes must not be OK"
    print("\nRESULT:", "ALL PASS" if failures == 0 else f"{failures} FAILED")
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
