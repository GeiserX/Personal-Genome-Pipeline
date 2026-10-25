#!/usr/bin/env python3
"""Both reports point to where a haplogroup can be explored further (steps 12 and 37).

Each case writes the haplogroup files of one made-up sample S into a bash
output folder, collects the summary and renders both reports:

  1. a male sample (Yleaf table with R-M269, haplogrep3 with H5a7): the Y and
     the mtDNA links are built from those names, in the text and the HTML
     report, with the files written for the services
  2. Yleaf's paragroup form R-M269*(xR-L23): the links use R-M269
  3. a female sample (no Yleaf table): no Y links; the mtDNA links stay
  4. no call (Hg NA, haplogroup '.'): no links at all
  5. a branch name with characters a URL cannot carry raw is escaped

The branch names are the fixture's and public examples, not anyone's result.

Run: python3 tests/test_haplogroup_explore.py
"""
import os
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402
import render_report  # noqa: E402

FAILS = []
Y_HEAD = "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n"


def check(desc, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {desc}" + ("" if ok else f": {detail}"))
    if not ok:
        FAILS.append(desc)


def put(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def render(d):
    summ = collect_summary.collect("S", d)
    return summ, render_report.text_report(summ), render_report.html_report(summ)


def sample(tmp, name, y_hg=None, mt_hg="H5a7", files=True):
    d = os.path.join(tmp, name)
    put(f"{d}/mito/S_haplogroup.txt", f'"SampleID"\t"Haplogroup"\t"Rank"\n"S"\t"{mt_hg}"\t"0.98"\n')
    if files:
        put(f"{d}/mito/S_chrM.fasta", ">S_chrM\nGATCACAGGT\n")
    if y_hg is not None:
        put(f"{d}/y_haplogroup/S_y_haplogroup.txt", Y_HEAD + f"S_sorted\t{y_hg}\tM269\t100\t42\t0.97\t1\t1\t0.97\n")
        if files:
            put(f"{d}/y_haplogroup/S_chrY.bam", "BAM")
    return d


def main():
    with tempfile.TemporaryDirectory() as tmp:
        # 1. male
        summ, txt, html = render(sample(tmp, "male", y_hg="R-M269"))
        for url in ("https://discover.familytreedna.com/y-dna/R-M269", "https://www.yfull.com/tree/R-M269/",
                    "https://discover.familytreedna.com/mtdna/H5a7", "https://www.yfull.com/mtree/H5a7/"):
            check(f"male: {url} in the text report", url in txt, txt)
            check(f"male: {url} as a link in the HTML report", f'<a href="{url}">' in html, "")
        check("male: the docs section for the routes is linked",
              "docs/37-y-haplogroup.md#explore-your-branch-further" in txt
              and "docs/12-mito-haplogroup.md#explore-your-branch-further" in txt, txt)
        check("male: the chrY BAM is named as the file to share",
              "Files for those services: y_haplogroup/S_chrY.bam" in txt and "y_haplogroup/S_chrY.bam" in html, txt)
        check("male: the mtDNA FASTA is named as the file to share",
              "Files for those services: mito/S_chrM.fasta" in txt and "mito/S_chrM.fasta" in html, txt)
        check("male: the reports say the pipeline uploads nothing", "Nothing is uploaded by this pipeline" in txt, txt)

        # 2. paragroup
        summ, txt, html = render(sample(tmp, "para", y_hg="R-M269*(xR-L23)"))
        check("paragroup: the link is built from R-M269",
              "https://www.yfull.com/tree/R-M269/" in txt and "y-dna/R-M269*" not in txt and "%2A" not in txt, txt)

        # 3. female: no Yleaf table, so no Y section and no Y link
        summ, txt, html = render(sample(tmp, "female"))
        check("female: no Y link", "y-dna/" not in txt and "yfull.com/tree/" not in txt
              and "y-dna/" not in html and "yfull.com/tree/" not in html, txt)
        check("female: the mtDNA links stay", "https://discover.familytreedna.com/mtdna/H5a7" in txt
              and "https://www.yfull.com/mtree/H5a7/" in html, txt)

        # 4. no call
        summ, txt, html = render(sample(tmp, "nocall", y_hg="NA", mt_hg=""))
        check("no call: no explore links", "discover.familytreedna.com" not in txt and "yfull.com" not in txt
              and "discover.familytreedna.com" not in html, txt)
        check("no call: the Y section still says insufficient markers", "insufficient markers" in txt, txt)

        # 5. escaping (a made-up PhyloTree-style name)
        summ, txt, html = render(sample(tmp, "esc", mt_hg="L3e2b+152", files=False))
        check("escaped: '+' is encoded in the mtDNA links", "mtdna/L3e2b%2B152" in txt and "mtree/L3e2b%2B152/" in txt, txt)
        check("no files written: no 'Files for those services' line", "Files for those services" not in txt, txt)

    print(f"\n{len(FAILS)} failed" if FAILS else "\nall passed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
