# What 1.0 Keeps Stable

OwnGenome 1.0 makes one promise to anyone who scripts around it: the command options, the config keys and the names and formats of the result files do not break in a 1.x release. A break waits for 2.0. This page says what the promise covers, what counts as a break, and how CI holds it.

The full list is one file, [`assets/output_contract.tsv`](https://github.com/GeiserX/OwnGenome/blob/main/assets/output_contract.tsv). Every row there is either `stable`, `experimental` or `internal`.

## What is covered

- **Nextflow params:** every param in `nextflow_schema.json` (`--input`, `--outdir`, `--tools`, `--clinvar`, ...), by name and meaning.
- **`--tools` names:** every name `main.nf` accepts (`pharmcat`, `clinvar`, `vcfanno`, ...).
- **`run-all.sh`:** its two arguments (`<sample> <male|female>`) and every switch in its `SWITCHES` list (`TOOLS`, `THREADS`, `SKIP_TRIM`, `INTERVALS`, `REF_FASTA`, ...).
- **Result files:** every file the Nextflow pipeline publishes, which is also what `run-all.sh` writes, at its path under `GENOME_DIR/<sample>/` (or `--outdir`). For each one the contract gives its step, its format and the fields a reader can count on: the columns of a table, the INFO and FORMAT fields of a VCF header, the top-level keys of a JSON file (`summary.json`'s keys are the ones its schema, `tests/schema/summary.schema.json`, requires).
- **Single-step scripts:** the files listed with leg `bash` or `both`. A row is `both` when the step's script (`scripts/<step>-*.sh`) writes the same path under `GENOME_DIR`. Some single steps write to other folders than the pipeline does (step 7 writes PharmCAT's reports to `vcf/`, the pipeline to `pharmcat/`; step 16b writes to `mosdepth/`, the pipeline to `coverage/`). Both names are kept. A single-step output that is not in the file is not part of the promise yet.

## What counts as a break

A 1.x release never does any of these to a `stable` row:

- remove or rename a file, or move it to another folder;
- remove or rename a column of a table, an INFO or FORMAT field of a VCF, or a top-level key of a JSON file;
- remove or rename a param, a `--tools` name or a `run-all.sh` switch, or change what it means.

These are not breaks, and can come in any 1.x release:

- a new file, column, field, key, param, tool or switch;
- a new row in a table, or a different value where the science changed: a newer ClinVar gives other hits, a better caller gives other calls. The contract is about the shape of the results, not about their content;
- a change to an `experimental` or `internal` row.

When a 1.x change touches an output on purpose, its pull request says so, and we add a column or a file rather than rename one.

## Stable, experimental and internal

- **`stable`:** kept until 2.0.
- **`experimental`:** an opt-in engine that is new in 1.0 and has not been run by many people yet: Y haplogroups and Y-STRs, Parascopy and Paraphase, blood groups, Lp(a), VNtyper, viral DNA, CHIP, MEI, MoChA, Exomiser, archaic and ancient ancestry, fine-grained ancestry, sample identity (somalier and VerifyBamID2), the CRAM archive, the family and couple steps, KIR and OptiType. Their files, params, tools and switches are listed, so you can see what exists, but a 1.x release may change them. Each change is named in the release notes. An experimental engine becomes stable in a 1.x release by changing its rows; it never goes the other way.
- **`internal`:** a task file that the pipeline does not publish (a temporary VCF, a work index). It is listed so the list stays complete, and it promises nothing.

The engines that run by default, and every step that existed before 1.0, are `stable`.

## How CI holds it

Two checks read the contract, both with [`scripts/ci/check-output-contract.py`](https://github.com/GeiserX/OwnGenome/blob/main/scripts/ci/check-output-contract.py):

- **On every pull request, with no run** (`tests/test_output_contract.py`, in the `unit-tests` job of the Guard workflow): every param of `nextflow_schema.json`, every `--tools` name of `main.nf` and every switch of `run-all.sh` has a row, and every row still names one. Every file row of a Nextflow process still matches a file that process writes and publishes, in the folder of a `publishDir` that keeps it (a `saveAs` that returns `null` for the file does not count), and every file a module publishes has a row. A published row with leg `nextflow` fails when its step's script writes the same path too, since its leg is then `both`; the check reads the script with its shell variables filled in, so a path it cannot follow is not flagged, and a `both` row it does not find is not refused. Every default in `nextflow_schema.json` must equal the value in `nextflow.config`. A removed or renamed stable item fails as `BREAKING`. A new one with no row fails too, with the message to add it, so the contract is always the full list. The test also plants each kind of break on a copy of the code (a renamed output, a removed param, a removed switch, a new tool and a new output with no row, a schema default that differs from `nextflow.config`, a row in the folder a `saveAs` drops the file from, and a `nextflow` row a step script also writes) and requires each to fail.
- **On the HG002 end-to-end run** (`tests/e2e/output-contract-1-hg002.sh`): every stable file whose `e2e` column names a leg must be in the results that leg wrote, under its contracted name, with its required fields. The bash leg is the single-step scripts the base cases run; the Nextflow leg is case `60-nextflow`. Rows whose step the fixture run does not reach have `e2e` set to `-` and are held by the static check only. The two files `run-all.sh` writes itself, `run_manifest.tsv` and `logs/run_status.tsv`, are listed but no check holds them yet: neither the fixture run nor a module writes them. The same case copies the Nextflow results, drops the `clnsig` column of the ClinVar hits table and renames the ROH summary, and requires the check to fail on both. It also prints every file each leg wrote with its fields, which is what to read when the contract is written again.

## Changing the contract

- **A new output, param, tool or switch:** add its row in the same pull request. Make it `experimental` when it belongs to a new opt-in engine, `stable` otherwise.
- **A break that cannot wait for 2.0:** keep the old name beside the new one for the rest of 1.x, with both rows `stable`, say on the step's doc page that the old name is deprecated, and say it in the pull request and the release notes.
- **Promoting an engine:** change its rows from `experimental` to `stable` and name it in the release notes.

[Testing](testing.md) has how CI runs; [Pipeline overview](pipeline-overview.md) has what each step is for.
