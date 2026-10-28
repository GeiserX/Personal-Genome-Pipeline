#!/usr/bin/env python3
"""bin/strchive_convert.py and the motif count on its catalog.

Reads the two STRchive files versions.env pins (STRCHIVE_VERSION and the
two sha256), downloaded once into a cache directory (STRCHIVE_DIR, default
<tmp>/pgp-strchive-<version>) and checked against their sha256, and checks:

  1. the files: a copy with one byte flipped fails the sha256 check (the
     control that the check can fail);
  2. the converter: one entry per STRchive ExpansionHunter locus (82), no
     name twice, every ExpansionHunter LocusId traced to one entry by its
     STRchive id; the report's key loci keep their names (HTT, FMR1,
     C9ORF72, ATXN1, DMPK, RFC1); the thresholds of HTT, FMR1 and the FAME
     locus SAMD12; the loci without thresholds and why; and a planted extra
     ExpansionHunter locus with no STRchive record stops it (exit 2, and
     neither the catalog nor the ranges table is written);
  3. bin/rfc1_motif.py on the catalog's motifs: every RFC1 case of
     tests/test_rfc1_motif.py gives the same verdict with the catalog's RFC1
     motifs and locus as with the built-in ones; at SAMD12 (FAME1) reads
     with a TTTCA expansion give pathogenic_motif and reads with a TTTTA
     expansion benign_motif, reads of the reference allele are not counted,
     and the control: the built-in RFC1 motifs call the same TTTCA reads
     uncertain, so the catalog's motifs are what label them;
  4. Stranger (STRANGER_IMAGE, through docker) on the catalog: HTT at its
     benign maximum, one above it and at its pathogenic minimum gives
     normal, pre_mutation and full_mutation; SAMD12 at the reference count
     normal and one above pre_mutation; MUC1 no status; and the control:
     Stranger's bundled GRCh38 catalog does not call HTT 27 pre_mutation, so
     the check tells the two catalogs apart.

Run: python3 tests/test_str_catalog.py
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests"))
import rfc1_motif as rm  # noqa: E402
import strchive_convert as sc  # noqa: E402

FAILS = []


def check(desc, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {desc}" + ("" if ok else f": {detail}"))
    if not ok:
        FAILS.append(desc)


def env(name):
    for line in open(os.path.join(REPO, "versions.env")):
        if line.startswith(name + "="):
            return line.split("=", 1)[1].split("#")[0].strip().strip('"')
    raise SystemExit(f"FAIL no {name} line in versions.env")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def files():
    """(release, EH catalog path, loci path), downloaded and checked."""
    rel = env("STRCHIVE_VERSION")
    d = os.environ.get("STRCHIVE_DIR") or os.path.join(tempfile.gettempdir(), f"pgp-strchive-{rel}")
    os.makedirs(d, exist_ok=True)
    out = []
    for name, var in ((f"STRchive-disease-loci-{rel}.hg38.expansionhunter.json", "STRCHIVE_EH_SHA256"),
                      (f"STRchive-loci-{rel}.json", "STRCHIVE_LOCI_SHA256")):
        path = os.path.join(d, name)
        if not (os.path.exists(path) and sha256(path) == env(var)):
            url = f"https://github.com/dashnowlab/STRchive/releases/download/{rel}/{name}"
            for i in range(3):
                try:
                    with urllib.request.urlopen(url, timeout=60) as r, open(path + ".part", "wb") as f:
                        shutil.copyfileobj(r, f)
                    break
                except OSError as e:
                    print(f"download {i + 1}/3 of {url} failed: {e}")
                    time.sleep(5 * (i + 1))
            os.replace(path + ".part", path)
        if sha256(path) != env(var):
            raise SystemExit(f"FAIL {name} does not have the sha256 {var} pins")
        out.append(path)
    return rel, out[0], out[1]


def build_catalog():
    """The converted catalog (list of entries) and its path, in a temp file; the
    ranges table (--ranges) is str_ranges.tsv beside it."""
    rel, eh, loci = files()
    out = os.path.join(tempfile.mkdtemp(), "catalog.json")
    rc = sc.main(["--eh", eh, "--loci", loci, "--release", rel, "--out", out,
                  "--ranges", os.path.join(os.path.dirname(out), "str_ranges.tsv")])
    if rc != 0:
        raise SystemExit("FAIL the converter failed on the pinned files")
    return json.load(open(out)), out


def part_files(rel, eh, loci):
    with tempfile.TemporaryDirectory() as tmp:
        bad = os.path.join(tmp, "flipped.json")
        data = bytearray(open(eh, "rb").read())
        data[len(data) // 2] ^= 1
        open(bad, "wb").write(data)
        check("control: a copy with one byte flipped fails the sha256 check",
              sha256(bad) != env("STRCHIVE_EH_SHA256"))


def part_convert(rel, eh, loci, cat, cat_path):
    src = json.load(open(eh))
    check("STRchive's ExpansionHunter catalog has 82 loci", len(src) == 82, len(src))
    check("the catalog has one entry per ExpansionHunter locus", len(cat) == len(src), len(cat))
    by_src = {}
    for e in cat:
        by_src.setdefault(e["SourceId"], []).append(e["LocusId"])
    check("every ExpansionHunter LocusId has exactly one entry",
          sorted(by_src) == sorted(x["LocusId"] for x in src) and all(len(v) == 1 for v in by_src.values()),
          {k: v for k, v in by_src.items() if len(v) != 1})
    names = [v for e in cat for v in (e.get("VariantId") if isinstance(e.get("VariantId"), list) else [e["LocusId"]])]
    check("no variant name twice", len(names) == len(set(names)), sorted(n for n in names if names.count(n) > 1))
    by = {e["LocusId"]: e for e in cat}
    check("the report's key loci and RFC1 keep their names",
          {"HTT", "FMR1", "C9ORF72", "ATXN1", "DMPK", "RFC1"} <= set(by), sorted(by))
    check("one gene with two loci keeps STRchive's ids", {"EIEE1_ARX", "PRTS_ARX", "HFG_HOXA13-I"} <= set(by))
    htt = by["HTT"]
    check("HTT: normal to 26, pathogenic from 36 (STRchive), the CAG repeat is REPID HTT",
          (htt.get("NormalMax"), htt.get("PathologicMin"), htt.get("VariantId")) == (26, 36, ["HTT", "HTT_CCG"]), htt)
    check("FMR1: normal to 44, pathogenic from 201", (by["FMR1"].get("NormalMax"), by["FMR1"].get("PathologicMin"))
          == (44, 201), by["FMR1"])
    s12 = by["SAMD12"]
    check("SAMD12 (FAME1): no benign maximum in STRchive, so normal is the reference count (1 TGAAA), pathogenic "
          "from 105; TTTTA benign, TTTCA pathogenic",
          (s12.get("NormalMax"), s12.get("PathologicMin"), s12.get("BenignMotifs"), s12.get("PathogenicMotifs"),
           s12.get("VariantId")) == (1, 105, ["TAAAA"], ["TGAAA"], ["SAMD12_TAAAA", "SAMD12"]), s12)
    check("RFC1: AAAGG is counted benign, CAGGG pathogenic",
          (by["RFC1"].get("BenignMotifs"), by["RFC1"].get("PathogenicMotifs"))
          == (["AAAAG", "AAAGG"], ["AAGGG", "ACAGG", "CAGGG"]), by["RFC1"])
    check("motif loci: RFC1, the FAME loci, SCA31 and SCA37",
          {"RFC1", "SAMD12", "STARD7", "YEATS2", "MARCHF6", "RAPGEF2", "TNRC6A", "RAI1", "BEAN1", "DAB1"}
          <= {e["LocusId"] for e in cat if "PathogenicMotifs" in e})
    none = {e["LocusId"]: e["NoThresholds"] for e in cat if "NoThresholds" in e}
    check("no thresholds: the four repeats longer than two reads, two without a pathogenic size, one contraction",
          sorted(none) == ["CEL", "MIR7-2", "MUC1", "PLIN4", "POLG", "TYMS", "VWA1"], none)
    check("a locus without thresholds has no NormalMax, so Stranger gives it no status",
          not any("NormalMax" in by[k] or "PathologicMin" in by[k] for k in none))
    check("every threshold pair is normal below pathogenic",
          all(e["NormalMax"] < e["PathologicMin"] for e in cat if "NormalMax" in e))
    check("the release is in SourceDisplay", by["HTT"].get("SourceDisplay") == f"STRchive-{rel}", by["HTT"])
    with tempfile.TemporaryDirectory() as tmp:
        planted = os.path.join(tmp, "eh.json")
        json.dump(src + [{"LocusId": "PLANTED_GENE", "LocusStructure": "(CAG)*",
                          "ReferenceRegion": "chr1:1000-1030", "VariantType": "Repeat"}], open(planted, "w"))
        out, rng = os.path.join(tmp, "out.json"), os.path.join(tmp, "ranges.tsv")
        rc = sc.main(["--eh", planted, "--loci", loci, "--release", rel, "--out", out, "--ranges", rng])
        check("control: a planted ExpansionHunter locus with no STRchive record stops the converter (exit 2, "
              "no catalog and no ranges table)", rc == 2 and not os.path.exists(out) and not os.path.exists(rng), rc)


# --- 3. the motif count on the catalog --------------------------------------------------
class Read:
    """The pysam.AlignedSegment attributes rfc1_motif reads."""

    def __init__(self, name, contig, start, cigar, seq, unmapped=False):
        self.query_name, self.reference_start, self.query_sequence = name, start, seq
        self.reference_id = self.next_reference_id = 0
        self.next_reference_start = start
        self.cigartuples = rm.parse_cigar(cigar) if cigar else None
        self.is_read1, self.is_unmapped, self.is_reverse = True, unmapped, False
        self.is_duplicate = self.is_secondary = self.is_supplementary = self.is_qcfail = False
        self.is_paired, self.mate_is_unmapped = True, False


def count(spec, reads):
    chrom, start, end = spec["span"]
    units = rm.min_units(start, end, spec["period"])
    counts = rm.count_motifs(rm.read_segments(reads, start - rm.FLANK, end + rm.FLANK, None, True),
                             units, spec["period"], spec["known"])
    return counts, rm.verdict(counts, units, spec)


def part_motifs(cat_path):
    import test_rfc1_motif as t
    specs = rm.catalog_specs(cat_path)
    sp = specs["RFC1"]
    check("RFC1 from the catalog: chr4, every repeat of the entry (AAAAG and AAGGG), unit 5, CANVAS, recessive",
          (sp["span"], sp["period"], sp["disease"], sp["recessive"]) == (("chr4", 39348424, 39348483), 5, "CANVAS", True),
          sp)
    check("RFC1: the catalog's locus asks for the same 12 units as the built-in one",
          rm.min_units(*sp["span"][1:], sp["period"]) == rm.min_units(t.START, t.END) == 12)
    for kind in ("benign", "pathogenic", "mixed", "acagg", "unknown", "few"):
        at, elsewhere = t.case(kind)
        mate_of = lambda r: elsewhere.get(r.query_name)  # noqa: E731
        lo, hi = t.START - rm.FLANK, t.END + rm.FLANK
        units = rm.min_units(t.START, t.END)
        plain = rm.verdict(rm.count_motifs(rm.read_segments(at, lo, hi, mate_of, True), units), units)
        _, s0, e0 = sp["span"]
        u2 = rm.min_units(s0, e0, 5)
        cat = rm.verdict(rm.count_motifs(rm.read_segments(at, s0 - rm.FLANK, e0 + rm.FLANK, mate_of, True), u2,
                                         5, sp["known"]), u2, sp)
        check(f"RFC1 {kind}: the catalog's motifs give the built-in verdict ({plain[4]}), shares and note",
              cat == plain, (cat, plain))

    fame = specs["SAMD12"]
    chrom, start, end = fame["span"]
    check("SAMD12 from the catalog: chr8, the TTTTA and TTTCA repeats, unit 5",
          (chrom, start, end, fame["period"]) == ("chr8", 118366812, 118366918, 5), fame)
    left = "GATTACACCTGGTCAGTGCA" * 3
    ref_allele = left + "TAAAA" * 20 + "TGAAA" + left[:30]   # the reference: (TAAAA)20 TGAAA

    def reads(unit, n, tag):
        return [Read(f"{tag}{i}", chrom, start, "150M", unit * 30) for i in range(n)]

    ref_reads = [Read(f"r{i}", chrom, start - 60, f"{len(ref_allele)}M", ref_allele) for i in range(6)]
    counts, v = count(fame, reads("TTTCA", 5, "c") + ref_reads)
    check("SAMD12: five reads of a TTTCA expansion are pathogenic_motif (TGAAA), the six reference reads not counted",
          (counts, v[0], v[4]) == ({"TGAAA": 5}, "TGAAA", "pathogenic_motif"), (counts, v))
    check("SAMD12: the note does not ask for a second copy (FAME1 is dominant)", "both copies" not in v[5], v[5])
    counts, v = count(fame, reads("TTTTA", 5, "a") + ref_reads)
    check("SAMD12: five reads of a TTTTA expansion are benign_motif (TAAAA)",
          (counts, v[0], v[4]) == ({"TAAAA": 5}, "TAAAA", "benign_motif"), (counts, v))
    counts, v = count(fame, ref_reads)
    check("SAMD12: the reference allele alone is not resolved (20 units, the floor is 22)",
          (counts, v[4], rm.min_units(start, end)) == ({}, "not_resolved", 22), (counts, v))
    plain = dict(fame, known=rm.KNOWN)
    _, v = count(plain, reads("TTTCA", 5, "c"))
    check("control: with the built-in RFC1 motifs the same TTTCA reads are uncertain_motif", v[4] == "uncertain_motif", v)


# --- 4. Stranger on the catalog -----------------------------------------------------------
FMT = "GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC"


def stranger(tmp, catalog, records):
    """STR_STATUS per record id from STRANGER_IMAGE on a VCF of (id, chrom, pos, end, count) records."""
    head = [line for line in open(os.path.join(REPO, "tests", "smoke", "stranger_input.vcf.in")) if line.startswith("##")]
    alts = sorted({n for *_, n in records})
    head += [f'##ALT=<ID=STR{a},Description="Allele comprised of {a} repeat units">\n' for a in alts]
    head += [f"##contig=<ID={c},length=250000000>\n" for c in sorted({r[1] for r in records})]
    with open(os.path.join(tmp, "in.vcf"), "w") as f:
        f.writelines(head)
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for i, (rid, chrom, pos, end, n) in enumerate(records):
            f.write(f"{chrom}\t{pos}\t{rid}.{i}\tA\t<STR{n}>\t.\tPASS\tSVTYPE=STR;END={end};REF=1;RL=3;RU=CAG;"
                    f"VARID={rid};REPID={rid}\t{FMT}\t"
                    f"1/1:SPANNING/SPANNING:{n}/{n}:{n}-{n}/{n}-{n}:20/20:0/0:0/0:30.0\n")
    # Without a catalog: the GRCh38 catalog Stranger bundles, as step 9b passes it.
    run = (["stranger", "--repeats-file", f"/w/{os.path.basename(catalog)}", "/w/in.vcf"] if catalog else
           ["sh", "-c", 'exec stranger --repeats-file "$(python -c \'import os, stranger; print(os.path.join('
            'os.path.dirname(stranger.__file__), "resources", "variant_catalog_grch38.json"))\')" /w/in.vcf'])
    cmd = ["docker", "run", "--rm", "--network", "none", "-u", f"{os.getuid()}:{os.getgid()}", "-v", f"{tmp}:/w",
           env("STRANGER_IMAGE")] + run
    print("+ " + " ".join(cmd), flush=True)
    try:
        p = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        return None
    if p.returncode != 0:
        print(p.stderr)
        return None
    out = {}
    for line in p.stdout.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(x.split("=", 1) if "=" in x else (x, "") for x in f[7].split(";"))
        out[f[2]] = info.get("STR_STATUS", "")
    return out


def part_stranger(cat_path):
    recs = [("HTT", "chr4", 3074877, 3074933, 26), ("HTT", "chr4", 3074877, 3074933, 27),
            ("HTT", "chr4", 3074877, 3074933, 36), ("SAMD12", "chr8", 118366913, 118366918, 1),
            ("SAMD12", "chr8", 118366913, 118366918, 2), ("MUC1", "chr1", 155188506, 155192239, 60)]
    with tempfile.TemporaryDirectory() as tmp:
        os.chmod(tmp, 0o755)
        shutil.copy(cat_path, os.path.join(tmp, "catalog.json"))
        got = stranger(tmp, os.path.join(tmp, "catalog.json"), recs)
        if got is None:
            check("Stranger runs on the catalog (docker and STRANGER_IMAGE)", False, "it did not run")
            return
        want = {"HTT.0": "normal", "HTT.1": "pre_mutation", "HTT.2": "full_mutation",
                "SAMD12.3": "normal", "SAMD12.4": "pre_mutation", "MUC1.5": ""}
        check("Stranger on the catalog: HTT 26 normal, 27 pre_mutation, 36 full_mutation; SAMD12 1 normal, 2 "
              "pre_mutation; MUC1 no status", got == want, got)
        bundled = stranger(tmp, None, recs[:3])
        check("control: Stranger's bundled GRCh38 catalog grades HTT 27 otherwise",
              bundled is not None and bundled.get("HTT.1") != "pre_mutation", bundled)


def main():
    rel, eh, loci = files()
    part_files(rel, eh, loci)
    cat, cat_path = build_catalog()
    part_convert(rel, eh, loci, cat, cat_path)
    part_motifs(cat_path)
    part_stranger(cat_path)
    if FAILS:
        print(f"{len(FAILS)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
