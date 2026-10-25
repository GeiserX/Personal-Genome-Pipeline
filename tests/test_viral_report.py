#!/usr/bin/env python3
"""bin/viral_report.py: the viral DNA in the reads that do not map to the human genome (step 47).

Each case writes a synthetic viral FASTA, the samtools coverage and depth -aa
files with chosen depths, and a human coverage file, runs the script, and
checks the rows it writes:

  1. the integrated HHV-6 pattern: the whole HHV-6B genome at half the human
     depth is called inherited integrated HHV-6; just above and just below
     Liu et al. 2020's 0.30 threshold flip the call
  2. low-level reads (viraemia-like depth over the whole genome) are viral DNA
     found but not integration; a few reads over a small part are "a few
     reads"; under 10 reads is not found
  3. phiX is the lab control whatever its depth, and never integration
  4. both HHV-6 species over the threshold: only the deeper one is called
  5. no read on chr1 to chr22: no human depth, integration not assessed
  6. the human depth is over the autosomal bases with reads: chrX, chrY and
     unplaced contigs are left out, and so are bases with no read
  7. damaged input stops the script: a depth file with missing positions,
     a coverage file that is not samtools coverage output; a genome with no
     position at all in the depth file has no read, and is not found
  8. collect_summary.py's viral section and render_report.py's text block
     read the table

Run: python3 tests/test_viral_report.py
"""
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "viral_report.py")
sys.path.insert(0, os.path.join(REPO, "bin"))
import viral_report as vr  # noqa: E402

FAILS = []
COV_HEAD = "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
# accession, title, length: small genomes keep the files small; the rules read fractions, not sizes.
GENOMES = [("NC_001664.4", "Human betaherpesvirus 6A, variant A DNA, complete virion genome, isolate U1102", 2000),
           ("NC_000898.1", "Human herpesvirus 6B, complete genome", 2000),
           ("NC_007605.1", "Human gammaherpesvirus 4, complete genome", 1500),
           ("NC_004102.1", "Hepatitis C virus genotype 1, complete genome", 900),
           ("NC_001422.1", "Escherichia phage phiX174, complete genome", 500)]


def check(name, got, want):
    if got == want:
        print(f"PASS {name}")
    else:
        print(f"FAIL {name}: got {got!r}, want {want!r}")
        FAILS.append(name)


def write_inputs(d, depths, reads, human=None):
    """depths: {acc: function(pos) -> depth}; reads: {acc: n}; human: [(rname, length, covbases, meandepth)]."""
    with open(os.path.join(d, "v.fa"), "w") as f:
        for acc, title, length in GENOMES:
            f.write(f">{acc} {title}\n" + "A" * length + "\n\n")
    with open(os.path.join(d, "vcov.tsv"), "w") as f, open(os.path.join(d, "vdepth.tsv"), "w") as g:
        f.write(COV_HEAD)
        for acc, _, length in GENOMES:
            fn = depths.get(acc, lambda p: 0)
            dep = [fn(p) for p in range(1, length + 1)]
            cov = sum(1 for x in dep if x)
            f.write(f"{acc}\t1\t{length}\t{reads.get(acc, 0)}\t{cov}\t{100 * cov / length:.4g}\t"
                    f"{sum(dep) / length:.6g}\t30\t60\n")
            g.writelines(f"{acc}\t{p}\t{x}\n" for p, x in enumerate(dep, 1))
    if human is None:
        human = [("chr1", 100000, 100000, 30.0), ("chr22", 50000, 50000, 30.0)]
    with open(os.path.join(d, "hcov.tsv"), "w") as f:
        f.write(COV_HEAD)
        for rname, length, covbases, mean in human:
            f.write(f"{rname}\t1\t{length}\t1000\t{covbases}\t{100 * covbases / length:.4g}\t{mean:.6g}\t30\t60\n")


def run(depths=None, reads=None, human=None, raw=None):
    """Rows by accession, or the error text when the script fails."""
    with tempfile.TemporaryDirectory() as d:
        write_inputs(d, depths or {}, reads or {}, human)
        if raw:
            raw(d)
        out = os.path.join(d, "out.tsv")
        p = subprocess.run([sys.executable, SCRIPT, "--sample", "S", "--viral-fasta", f"{d}/v.fa",
                            "--viral-coverage", f"{d}/vcov.tsv", "--viral-depth", f"{d}/vdepth.tsv",
                            "--human-coverage", f"{d}/hcov.tsv", "--out", out], capture_output=True, text=True)
        if p.returncode:
            return p.stderr
        with open(out) as f:
            lines = [ln.rstrip("\n").split("\t") for ln in f]
        head = lines[0]
        check("the header is the documented one", head, vr.HEADER)
        return {r[1]: dict(zip(head, r)) for r in lines[1:]}


def flat(x):
    return lambda p: x


A, B, EBV, HCV, PHIX = (g[0] for g in GENOMES)

# 1. integration: whole genome at half the human depth (30x human, 15x viral)
ici = run({B: flat(15)}, {B: 400})
check("1 HHV-6B at half the human depth: inherited integrated HHV-6", ici[B]["call"], vr.CALL_ICI)
check("1 relative depth 0.5", ici[B]["relative_depth"], "0.500")
check("1 about one copy per cell", ici[B]["copies_per_cell"], "1.0000")
check("1 full breadth", ici[B]["breadth"], "1.0000")
check("1 the short name", ici[B]["virus"], "HHV-6B")
check("1 HHV-6A with no read is not found", ici[A]["call"], vr.CALL_NONE)
check("1 rows in FASTA order", list(ici), [g[0] for g in GENOMES])
check("1 the note cites the threshold", "Liu et al. 2020" in ici[B]["note"], True)
above = run({B: flat(9)}, {B: 300})
check("1 relative 0.30 is integration", (above[B]["relative_depth"], above[B]["call"]), ("0.300", vr.CALL_ICI))
below = run({B: flat(8)}, {B: 300})
check("1 relative 0.27 is not", (below[B]["relative_depth"], below[B]["call"]), ("0.267", vr.CALL_FOUND))
check("1 and its note says it is below the threshold", "Below the integration threshold" in below[B]["note"], True)
# The mode, not the mean: a telomere-like pile-up on a short stretch does not lift a genome with no reads elsewhere.
pile = run({B: lambda p: 400 if p <= 100 else 0}, {B: 300})
check("1 a 100-base pile-up at 400x has mode 0: not integration", (pile[B]["mode_depth"], pile[B]["call"]),
      ("0", vr.CALL_FEW))
check("1 and its mean is high", pile[B]["mean_depth"], "20.000")

# 2. low level, few reads, none
low = run({EBV: lambda p: 1 if p % 2 else 0}, {EBV: 40})
check("2 EBV over half the genome: viral DNA found", low[EBV]["call"], vr.CALL_FOUND)
check("2 its mode depth is 0 or 1, relative far below integration", float(low[EBV]["relative_depth"]) < 0.1, True)
few = run({EBV: lambda p: 1 if p <= 100 else 0}, {EBV: 12})
check("2 12 reads over 7% of EBV: a few reads", few[EBV]["call"], vr.CALL_FEW)
nine = run({EBV: lambda p: 1 if p <= 1000 else 0}, {EBV: 9})
check("2 nine reads: not found", nine[EBV]["call"], vr.CALL_NONE)
hcv = run({HCV: flat(2)}, {HCV: 50})
check("2 HCV DNA is flagged as an RNA virus", "RNA virus with no DNA stage" in hcv[HCV]["note"], True)
check("2 HHV-6 rule never applies to EBV at integration depth",
      run({EBV: flat(15)}, {EBV: 500})[EBV]["call"], vr.CALL_FOUND)

# 3. phiX
px = run({PHIX: flat(30)}, {PHIX: 200})
check("3 phiX at human depth is the lab control", px[PHIX]["call"], vr.CALL_CONTROL)
check("3 one phiX read is still the lab control", run({PHIX: lambda p: 1 if p < 5 else 0}, {PHIX: 1})[PHIX]["call"],
      vr.CALL_CONTROL)
check("3 no phiX read: not found", run()[PHIX]["call"], vr.CALL_NONE)

# 4. both species pass: the deeper one wins
both = run({A: flat(10), B: flat(15)}, {A: 300, B: 400})
check("4 HHV-6B (deeper) is the integrated call", both[B]["call"], vr.CALL_ICI)
check("4 HHV-6A is read as shared reads", (both[A]["call"], "Reads shared with HHV-6B" in both[A]["note"]),
      (vr.CALL_FOUND, True))

# 5. no human depth
nohum = run({B: flat(15)}, {B: 400}, human=[("chrX", 1000, 1000, 30.0)])
check("5 no autosome: human depth '.'", nohum[B]["human_depth"], ".")
check("5 and integration is not assessed", (nohum[B]["call"], "integration is not assessed" in nohum[B]["note"]),
      (vr.CALL_FOUND, True))

# 6. the human depth
hd = vr.human_depth({"chr1": {"length": 1000, "covbases": 500, "meandepth": 10.0, "reads": 0},
                     "22": {"length": 1000, "covbases": 1000, "meandepth": 40.0, "reads": 0},
                     "chrX": {"length": 1000, "covbases": 1000, "meandepth": 99.0, "reads": 0},
                     "chrUn_KI270302v1": {"length": 1000, "covbases": 1000, "meandepth": 99.0, "reads": 0},
                     "chrEBV": {"length": 1000, "covbases": 1000, "meandepth": 99.0, "reads": 0}})
check("6 human depth over autosomal bases with reads: (10*1000 + 40*1000) / 1500", round(hd, 4), round(50000 / 1500, 4))

# 7. damaged input


def cut_depth(d):
    with open(f"{d}/vdepth.tsv") as f:
        lines = f.readlines()
    with open(f"{d}/vdepth.tsv", "w") as f:
        f.writelines(lines[:-10])


def bad_cov(d):
    with open(f"{d}/hcov.tsv", "w") as f:
        f.write("chr1\t1\t2\n")


def drop_genome(d):
    with open(f"{d}/vdepth.tsv") as f:
        lines = [ln for ln in f if not ln.startswith(HCV + "\t")]
    with open(f"{d}/vdepth.tsv", "w") as f:
        f.writelines(lines)


gone = run({B: flat(15)}, {B: 400}, raw=drop_genome)
check("7 a genome with no read that depth left out entirely is not found",
      (gone[HCV]["call"], gone[HCV]["breadth"], gone[B]["call"]), (vr.CALL_NONE, "0.0000", vr.CALL_ICI))
err = run(raw=cut_depth)
check("7 a depth file with positions missing stops the script", isinstance(err, str) and "samtools depth -aa" in err, True)
err = run(raw=bad_cov)
check("7 a coverage file that is not samtools coverage output stops it",
      isinstance(err, str) and "not samtools coverage output" in err, True)

# 8. the summary and the text report
import collect_summary as cs  # noqa: E402
import render_report as rr  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    write_inputs(d, {B: flat(15), EBV: flat(1), PHIX: flat(3)}, {B: 400, EBV: 30, PHIX: 20})
    os.makedirs(f"{d}/S/viral")
    vr.main(["--sample", "S", "--viral-fasta", f"{d}/v.fa", "--viral-coverage", f"{d}/vcov.tsv",
             "--viral-depth", f"{d}/vdepth.tsv", "--human-coverage", f"{d}/hcov.tsv",
             "--out", f"{d}/S/viral/S_viral.tsv"])
    path, data = cs.sec_viral(f"{d}/S", "S")
    check("8 the summary finds the table", path.endswith("viral/S_viral.tsv"), True)
    check("8 the summary's HHV-6 line", data["hhv6"], "inherited integrated HHV-6B likely")
    check("8 the summary lists what was found", [x["virus"] for x in data["found"]],
          ["HHV-6B", "Human gammaherpesvirus 4", "Escherichia phage phiX174"])
    check("8 no table: no section", cs.sec_viral(f"{d}/none", "S"), (None, {}))
    lines = rr.viral_lines(data)
    check("8 the text block opens with the HHV-6 line", lines[0], "HHV-6: inherited integrated HHV-6B likely")
    check("8 and ends with the note", lines[-1], rr.VIRAL_NOTE)
    html = "".join(rr.viral_body(data))
    check("8 the HTML card names HHV-6B", "HHV-6B" in html, True)
    with open(f"{d}/S/viral/S_viral.tsv", "w") as f:
        f.write("virus\tcall\n")
    try:
        cs.sec_viral(f"{d}/S", "S")
        check("8 a table without its columns is unreadable", "read", "Unreadable")
    except cs.Unreadable:
        check("8 a table without its columns is unreadable", "Unreadable", "Unreadable")

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
