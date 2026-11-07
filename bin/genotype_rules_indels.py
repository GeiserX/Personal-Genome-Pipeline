#!/usr/bin/env python3
"""genotype_rules_indels.py: the APOL1 G2 deletion, read beside step 42's APOL1 rule.

Step 42's rule table (assets/genotype_rules.tsv, bin/traits_report.py) reads
single-base sites only, so its APOL1 rule reads G1 (two missense changes) and
says nothing about G2 (rs71785313), a 6-base deletion of N388 and Y389 a few
bases away (Genovese et al. 2010, PMID 20647424). G1 and G2 never share a copy,
and two risk copies (G1/G1, G1/G2 or G2/G2) are the high-risk genotype. This
module reads G2 from the sample's VCF and puts the two together:

  read_g2(records, kind)   the G2 copies (0, 1, 2) from the VCF records that
                           overlap G2, with a status:
                             called          a G2 call, or reference blocks
                                             (gVCF) or reference calls over
                                             the whole G2 span: 0 copies
                             no_call_record  a VCF without reference blocks and
                                             no record over G2: no G2 call, but
                                             a reference base cannot be told
                                             from a position with no reads
                             not_covered     a gVCF with no genotype over part
                                             of the G2 span
                             not_read        a no-call, a filtered G2 call, an
                                             allele other than G2 inside its
                                             span, a REF that is not the
                                             GRCh38 base, or more than one
                                             record calling G2 (one copy
                                             written twice cannot be told
                                             from two copies)
  combine(rule_line, g2)   the G1 + G2 result and its reading, from step 42's
                           APOL1 line (G0/G0, G0/G1, G1/G1) and read_g2's answer.

A deletion in a repeat can be written at several positions. A record is G2
when applying its allele to GRCh38_CONTEXT gives the same sequence as G2, so a
left-aligned and a right-shifted record both count. The pinned dbSNP record,
DBSNP, is what NCBI's Variation API gave for rs71785313 on 2026-11-06 (the VCF
fields of NC_000022.11:36265995:ATAATTATAA:ATAA); tests/e2e/genotype-rules-2-g2.sh
compares it with dbSNP's VCF of the build versions.env pins, and
tests/test_genotype_rules_indels.py checks the marker against it and against
GRCh38_CONTEXT. Standard library only.
"""
import re

# rsID, chrom, VCF POS, REF, ALT of rs71785313 in dbSNP's GRCh38 VCF (left-aligned).
DBSNP = ("rs71785313", "chr22", 36265995, "AATAATT", "A")
G2 = {"rule": "APOL1", "name": "G2", "rsid": "rs71785313", "chrom": "chr22", "pos": 36265995,
      "ref": "AATAATT", "alt": "A"}
# GRCh38 chr22:36265975-36266025, the plus strand (UCSC hg38, 2026-11-06). dbSNP's
# SPDI for G2 spans 36265996-36266005 (ATAATTATAA becomes ATAA): any 6-base
# deletion inside that span is the same sequence.
CONTEXT_START = 36265975
GRCh38_CONTEXT = "AGAAGCTAAACATTCTCAACAATAATTATAAGATTCTGCAGGCGGACCAAG"
CONTEXT_END = CONTEXT_START + len(GRCh38_CONTEXT) - 1
# The G2 span: the anchor base of the left-aligned record to the end of the
# ambiguous span. A genotype over every base of it reads G2.
SPAN = (36265995, 36266005)
G1_COPIES = {"G0/G0": 0, "G0/G1": 1, "G1/G1": 2}
SYMBOLIC = re.compile(r"^(<.*>|\*|\.)$")


def apply(pos, ref, alt):
    """GRCh38_CONTEXT with ref at pos replaced by alt, or None when ref is not
    the GRCh38 bases there or the record does not lie inside the context."""
    i = pos - CONTEXT_START
    if i < 0 or pos + len(ref) - 1 > CONTEXT_END or GRCh38_CONTEXT[i:i + len(ref)] != ref:
        return None
    return GRCh38_CONTEXT[:i] + alt + GRCh38_CONTEXT[i + len(ref):]


def marker_errors(m=G2, dbsnp=DBSNP):
    """The problems of the G2 marker, as messages: against the pinned dbSNP
    record (a swapped pair named as such) and against GRCh38_CONTEXT."""
    errs = []
    rec = (m["rsid"], m["chrom"], m["pos"], m["ref"], m["alt"])
    if rec != tuple(dbsnp):
        if (m["ref"], m["alt"]) == (dbsnp[4], dbsnp[3]):
            errs.append(f"{m['rsid']}: ref and alt are swapped (marker {m['ref']}>{m['alt']}, "
                        f"dbSNP REF {dbsnp[3]})")
        else:
            errs.append(f"{m['rsid']}: the marker is {m['chrom']}:{m['pos']} {m['ref']}>{m['alt']}, "
                        f"dbSNP {dbsnp[1]}:{dbsnp[2]} {dbsnp[3]}>{dbsnp[4]}")
    hap = apply(m["pos"], m["ref"], m["alt"])
    if hap is None:
        errs.append(f"{m['rsid']}: REF {m['ref']} is not the GRCh38 bases at {m['chrom']}:{m['pos']}")
    elif len(GRCh38_CONTEXT) - len(hap) != 6:
        errs.append(f"{m['rsid']}: the marker is not a 6-base deletion")
    return errs


G2_HAPLOTYPE = apply(G2["pos"], G2["ref"], G2["alt"])


def record_end(r):
    """The last base a VCF record covers: INFO END for a reference block, else POS + len(REF) - 1."""
    m = re.search(r"(?:^|;)END=(\d+)", r[7]) if len(r) > 7 else None
    return int(m.group(1)) if m else int(r[1]) + len(r[3]) - 1


def overlaps(r):
    """Whether a split VCF data line is on chr22 and covers part of the G2 span."""
    return (len(r) > 4 and r[0] == G2["chrom"] and r[1].isdigit()
            and int(r[1]) <= SPAN[1] and record_end(r) >= SPAN[0])


def read_g2(records, kind):
    """G2 from the VCF records that overlap SPAN (other records are skipped).
    kind is 'gvcf' (reference blocks: a base with no record has no reads) or
    'vcf' (variant calls only). Returns {status, copies, gq, text}; copies is
    None unless status is called or no_call_record."""
    copies, covered, gqs, problems, seen, g2_records = 0, set(), [], [], False, 0
    for r in records:
        if not overlaps(r):
            continue
        if len(r) < 10:
            problems.append("a record over G2 has no sample column")
            continue
        seen = True
        pos, ref = int(r[1]), r[3]
        alleles = [ref] + r[4].split(",")
        fmt = dict(zip(r[8].split(":"), r[9].split(":")))
        gt = re.split(r"[/|]", fmt.get("GT", "."))
        if any(not g.isdigit() for g in gt):
            problems.append(f"a no-call at {r[0]}:{pos}")
            continue
        idx = [int(g) for g in gt]
        if any(i >= len(alleles) for i in idx):
            problems.append(f"a genotype at {r[0]}:{pos} names an allele the record lacks")
            continue
        if all(i == 0 for i in idx):
            covered.update(range(pos, record_end(r) + 1))
            gqs.append(fmt.get("GQ", "."))
            continue
        if r[6] not in ("PASS", "."):
            problems.append(f"a call at {r[0]}:{pos} that failed the caller's filter ({r[6]})")
            continue
        n = 0
        for i in idx:
            if i == 0:
                continue
            a = alleles[i]
            inside = pos >= CONTEXT_START and record_end(r) <= CONTEXT_END
            if inside and apply(pos, ref, ref) is None:
                problems.append(f"REF {ref} at {r[0]}:{pos} is not the GRCh38 base (is the VCF on GRCh38?)")
                break
            hap = None if SYMBOLIC.match(a) or not inside else apply(pos, ref, a)
            if hap == G2_HAPLOTYPE:
                n += 1
            else:
                problems.append(f"an allele other than G2 inside its span ({r[0]}:{pos} {ref}>{a})")
                break
        else:
            if n:
                g2_records += 1
                if g2_records > 1:
                    # Two records of the one deletion (say a left-aligned and a shifted copy of
                    # the same call) cannot be told from two G2 copies.
                    problems.append("more than one record calls the G2 deletion")
                    break
            copies += n
            covered.update(range(pos, record_end(r) + 1))
            gqs.append(fmt.get("GQ", "."))
    gq = ",".join(gqs) if gqs else "."
    if problems:
        return {"status": "not_read", "copies": None, "gq": gq,
                "text": f"G2 (rs71785313) not read: {problems[0]}."}
    if copies:
        return {"status": "called", "copies": copies, "gq": gq,
                "text": f"G2 (rs71785313): {copies} cop{'ies' if copies > 1 else 'y'} of the deletion called."}
    if all(p in covered for p in range(SPAN[0], SPAN[1] + 1)):
        return {"status": "called", "copies": 0, "gq": gq,
                "text": "G2 (rs71785313): the reference over the whole deletion span, no G2 copy."}
    if kind == "vcf" and not seen:
        return {"status": "no_call_record", "copies": 0, "gq": ".",
                "text": "G2 (rs71785313): no G2 deletion among the VCF's calls."}
    return {"status": "not_covered", "copies": None, "gq": gq,
            "text": "G2 (rs71785313) not read: no genotype over part of the deletion span."}


TWO_COPIES = ("Two APOL1 risk copies ({label}): the high-risk genotype of Genovese 2010, which went with "
              "focal segmental glomerulosclerosis and hypertension-attributed kidney failure in African "
              "Americans. A risk factor, not a diagnosis.")


def combine(line, g2):
    """{result, reading} of G1 and G2 together, from step 42's APOL1 line and read_g2's answer."""
    g1 = G1_COPIES.get(line.get("result")) if line.get("status") == "called" else None
    if g1 is None and g2["status"] == "called" and g2["copies"] == 2:
        # Two G2 copies fill both copies of the gene: G2/G2 whatever the G1 sites read.
        return {"result": "G2/G2", "reading": TWO_COPIES.format(label="G2/G2") + " The G1 sites were not read, "
                "but two G2 copies leave no copy for G1."}
    if g1 is None:
        return {"result": ".", "reading": "Not read: the G1 sites of the APOL1 rule were not read."}
    if g2["copies"] is None:
        return {"result": ".", "reading": "The G1 result above stands on its own."}
    c = g2["copies"]
    if g1 + c > 2:
        return {"result": ".", "reading": f"{g1} G1 and {c} G2 copies do not fit two copies of the gene "
                                          "(G1 and G2 never share a copy): not read. Check the calls."}
    label = "/".join(["G0"] * (2 - g1 - c) + ["G1"] * g1 + ["G2"] * c)
    if g1 + c == 2:
        reading = TWO_COPIES.format(label=label)
    elif g1 + c == 1:
        reading = (f"One APOL1 risk copy ({label}): one copy did not raise the risk in Genovese 2010 "
                   "(a recessive pattern).")
    else:
        reading = "No G1 and no G2 copy: no APOL1 kidney risk copy."
    if g2["status"] == "no_call_record":
        label += " (no G2 call)"
        reading += (" No G2 deletion was called, but this VCF has no reference blocks, so a G2 copy at a "
                    "position with no reads would not show.")
    return {"result": label, "reading": reading}
