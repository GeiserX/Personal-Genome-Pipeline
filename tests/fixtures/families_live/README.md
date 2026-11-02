# families_live

A synthetic `GENOME_DIR` for `tests/test_pedigree_svg.py` and `tests/test_pgp_live_family.py`. Every name and number is invented; no file holds a genotype.

- `family_a/`: a family folder as steps 38 and 39 leave it. `members.tsv` has three generations: a sequenced paternal grandfather and grandmother (`gf_a`, `gm_a`), their son `father_a`, his partner `mother_a`, her mother `gm_m` (named, not sequenced), and two children, `child_a` and `child_b`. The other files stand in for step 38's and step 39's outputs, with made-up values. `notes.txt` is a file no step writes: the site must not serve it.
- `family_b/`: a couple with only `pedigree.tsv` (step 39 run, step 38 not), so the pedigree is drawn from it.
- `child_a/`: a sample run folder (`logs/run_status.tsv`), the run page that links to its family.
