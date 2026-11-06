# Nextflow Parameters

Every parameter of the Nextflow pipeline (`main.nf`), grouped as in [`nextflow_schema.json`](https://github.com/GeiserX/OwnGenome/blob/main/nextflow_schema.json), with its default and the description the schema gives. `--pedigree` is read by `main.nf` but is not in the schema yet; its row is written here. Give one as `--<name> <value>` to `nextflow run`, or after the sex to `run-all.sh`, which hands the options after the sex to `nextflow run` unchanged. `run-all.sh` already passes the databases it finds under `GENOME_DIR` ([Full run](getting-started.md#full-run)). How a run is built: [Nextflow Execution](nextflow.md).

A default of none means the parameter is unset; each description says what a run does without it.

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Default | What it does |
|---|---|---|
| `--input` | none | Path to comma-separated samplesheet: sample, then fastq_1,fastq_2 (aligned and called), or bam,bam_index (called), or vcf,vcf_index (a VCF with one sample) with an optional bam,bam_index; sex (male or female) is required on every row that is called and when expansion_hunter is enabled. A row can also give cram,crai instead of the BAM, and gvcf,gvcf_index beside the VCF; every column: [Samplesheet Format](nextflow.md#samplesheet-format). |
| `--pedigree` | none | The family workflow (step 38, and step 39 for two samples the pedigree marks as partners): a tab-separated pedigree (sample, sex, father, mother, partner; `assets/pedigree.example.tsv`) or a families file with a leading family column (`assets/families.example.tsv`), naming samplesheet rows. Needs `--somalier_sites` and at least two of the named samples. [Family workflow](nextflow.md) and [Multi-sample comparison](multi-sample.md#family-analysis-step-38-opt-in) say what it runs. |
| `--outdir` | `./results` | The output directory where results will be saved. |
| `--allow_unfiltered` | `false` | Accept a VCF that has no FILTER=PASS record by treating FILTER '.' as PASS. Without it such a VCF stops the run, because every PASS-only step would report zero hits. |
| `--publish_dir_mode` | `copy` | Method used to save pipeline results to output directory. One of: `symlink`, `rellink`, `link`, `copy`, `copyNoFollow`, `move`. |

## Reads to VCF

Alignment and variant calling of the samplesheet rows that start from FASTQ or from a BAM without a VCF.

| Parameter | Default | What it does |
|---|---|---|
| `--skip_trim` | `false` | Align the raw reads, without fastp trimming (SKIP_TRIM=true in scripts/01b-fastp-qc.sh). |
| `--minimap2_index` | none | minimap2 index built with -x sr. Default: &lt;reference base&gt;.sr.mmi beside the FASTA, the file scripts/02-alignment.sh builds; when neither exists the run builds one. |
| `--intervals` | none | Space-separated regions DeepVariant calls (its --regions, INTERVALS in scripts/03-deepvariant.sh). Unset calls the whole genome. |
| `--sex_check` | `fail` | What to do when the sex inferred from a BAM differs from the samplesheet's, by indexcov (every BAM, before any BAM step) or by somalier (with sample_qc): fail stops the run; warn logs it and goes on with the declared sex. One of: `fail`, `warn`. |

## Reference data

Paths to reference files required by the pipeline.

| Parameter | Default | What it does |
|---|---|---|
| `--reference` | none | Path to the GRCh38 reference FASTA the BAMs were aligned to: the no-ALT analysis set that scripts/setup.sh installs as reference/GRCh38_no_alt_analysis_set.fasta. |
| `--clinvar` | none | Path to ClinVar pathogenic VCF (bgzipped). Optional — skip ClinVar screen if not provided. |
| `--clinvar_index` | none | Path to ClinVar pathogenic VCF tabix index (.tbi). |

## Tool selection

Control which analysis tools to run.

| Parameter | Default | What it does |
|---|---|---|
| `--tools` | `pharmcat,cpic,vcfanno,roh,prs,mito_haplogroup,telomere_hunter,mosdepth,mito_variants,html_report,multiqc` | Comma-separated list of tools to run; an unknown name stops the run. DB-heavy tools (vep, slivar, clinical_filter, cpsr, clinvar, expansion_hunter, sample_qc) and extras (stranger, ancestry (needs --ancestry_ref and prs), ancestry_fine (needs --ancestry_fine_panel), ancient (needs --ancient_panel; samples with a gVCF), pypgx, SV callers, cram_archive, cyrius, parascopy, paraphase (PacBio HiFi BAMs only), y_haplogroup (male samples), y_str (Y-STR repeat counts, male samples), consumer_export (samples with a gVCF), bloodtype (blood group antigens; RhD needs the BAM), lpa (the LPA KIV-2 repeat count; needs the BAM), vntyper (MUC1 VNTR frameshifts behind ADTKD-MUC1; needs the BAM), viral (viral DNA in the unmapped reads; needs --viral_ref and a BAM), exomiser (runs when --hpo is set), chip (clonal haematopoiesis from the BAM; needs --chip_gff), archaic (needs --archaic_data), mei (Alu, L1 and SVA insertions with MEGAnE; needs --megane_kmers and a BAM), mocha (mosaic chromosomal alterations; needs --mocha_data), locityper (a second opinion at HLA class II, KIR, CYP2D6 and FCGR; samples with a BAM; needs --locityper_data)) are opt-in. prs and vcfanno are skipped with a warning until --pgs_scoring or a score file is set. With pharmcat, PharmCAT also gets the HLA types of hla_typing and a CYP2D6 call pypgx and cyrius agree on (PGX_CONSENSUS). |

## Annotation databases

Paths to annotation databases. Some tools (VEP, CPSR, ExpansionHunter) require their databases and will fail at startup if enabled without them. Others produce degraded output with warnings.

| Parameter | Default | What it does |
|---|---|---|
| `--vep_cache` | none | VEP offline cache directory (release 116). |
| `--vep_cache_cpsr` | none | VEP cache for CPSR (release 115, PCGR_VEP_CACHE_RELEASE in versions.env; separate from the VEP step's release 116 cache). |
| `--pcgr_data` | none | PCGR/CPSR reference data bundle directory. |
| `--pypgx_bundle` | none | pypgx-bundle reference data directory for star allele calling. |
| `--cadd_snv` | none | CADD whole-genome pre-scored SNV file (tsv.gz). |
| `--cadd_snv_index` | none | CADD SNV tabix index (.tbi). |
| `--cadd_indel` | none | CADD whole-genome pre-scored indel file (tsv.gz). |
| `--cadd_indel_index` | none | CADD indel tabix index (.tbi). |
| `--spliceai_snv` | none | SpliceAI pre-scored SNV VCF (vcf.gz). |
| `--spliceai_snv_index` | none | SpliceAI SNV tabix index (.tbi). |
| `--spliceai_indel` | none | SpliceAI pre-scored indel VCF (vcf.gz). |
| `--spliceai_indel_index` | none | SpliceAI indel tabix index (.tbi). |
| `--revel` | none | REVEL pathogenicity scores file. |
| `--revel_index` | none | REVEL tabix index (.tbi). |
| `--alphamissense` | none | AlphaMissense pathogenicity scores file. |
| `--alphamissense_index` | none | AlphaMissense tabix index (.tbi). |
| `--gnomad_constraint` | none | gnomAD gene constraint TSV (LOEUF, pLI scores). |

## Clinical & population genetics

Paths to clinical and population genetics reference data.

| Parameter | Default | What it does |
|---|---|---|
| `--pgs_scoring` | none | Folder of GRCh38-harmonised PGS Catalog scoring files (setup.sh downloads the scores of assets/pgs_scores.tsv into &lt;genome_dir&gt;/prs_scores). prs scores each with pgsc_calc. |
| `--ancestry_ref` | none | pgsc_calc's ancestry reference panel (pgsc_1000G_v1.tar.zst, installed by setup.sh --ancestry-panel with its _GRCh38_sites.tsv beside it). With it prs projects the sample onto the panel, reports each score as a percentile among the most similar population, and writes the ancestry table of step 26; without it the scores are raw sums. |
| `--ancestry_fine_panel` | none | Step 26b's reference: the folder setup.sh --ancestry-fine-panel writes (reference/ancestry_fine/hgdp1kgp_v2, gnomAD's phased HGDP + 1000 Genomes panel with ADMIXTURE's group frequencies). With ancestry_fine in --tools, ADMIXTURE projects the sample on the panel's groups and RFMix paints its chromosomes. With --pedigree, the family workflow also phases the joint call against it with Beagle and lists the DNA segments each pair of members shares (hap-ibd). ADMIXTURE is free for academic use and RFMix for academic research use only. |
| `--ancient_panel` | none | Step 44's reference: the folder setup.sh --ancient-panel writes (reference/ancient/aadr_v66.p1, the Allen Ancient DNA Resource's 1240K panel with smartpca's present-day West Eurasian axes). With ancient in --tools, smartpca projects the sample onto those axes and the report lists the nearest ancient groups. Samples with a gVCF only. |
| `--pgsc_calc` | none | The pgsc_calc checkout setup.sh makes (&lt;genome_dir&gt;/tools/pgsc_calc-&lt;release&gt;). With it prs runs pgsc_calc offline; without, the PRS task fetches GitHub's archive of the release versions.env pins, checked against its sha256. |
| `--prs_keep_ambiguous` | `false` | prs: keep each score's strand-ambiguous rows (A/T and C/G pairs) when its other rows agree with the reference on the forward strand, and pass --keep_ambiguous to pgsc_calc; a score that fails that check keeps dropping them. Rows with a gVCF (WGS) only: a row with a VCF alone stops the PRS task. Off by default: pgsc_calc drops every ambiguous row. The summary's Ambiguous column says what happened to them per score. |
| `--expansion_catalog` | none | ExpansionHunter variant catalog JSON. |
| `--hla_dat` | none | Pre-downloaded IPD-IMGT/HLA hla.dat file for T1K HLA typing. |
| `--hla_genes` | none | Gene annotation (GTF) T1K takes each HLA gene's GRCh38 coordinates from: the gene lines of GENCODE's basic annotation that scripts/setup.sh installs as reference/gencode.v&lt;release&gt;.basic.genes.gtf. Required with hla_typing. |
| `--gencode_exons` | none | The exons of GENCODE's basic annotation as a BED, which scripts/setup.sh installs as reference/gencode.v&lt;release&gt;.basic.exons.bed. With clinvar, CLINVAR_CDNA_CHECK reads it for the cDNA read check of the ClinVar hits (without it every hit is not checked); with --pedigree, the couple carrier screen uses deletions over an exon. |
| `--mito_clinvar` | none | The whole ClinVar VCF with chr names, its .tbi beside it, which scripts/setup.sh installs as clinvar/clinvar_chr.vcf.gz. With mito_variants, MITO_VARIANTS_ANNOTATE gives each chrM call its ClinVar record, as step 20 does; without it the ClinVar columns of &lt;id&gt;_chrM_annotated.tsv are '.'. |
| `--gnomad_mito` | none | gnomAD v3.1's chrM sites VCF (GNOMAD_MITO_URL in versions.env), which scripts/setup.sh installs under gnomad_mito/. With mito_variants, MITO_VARIANTS_ANNOTATE gives each chrM call its gnomAD frequencies, as step 20 does; without it those columns are '.'. |
| `--kir` | `false` | With hla_typing: also type the KIR genes, in a second T1K pass with its kir-wgs preset. Needs --kir_dat. |
| `--kir_dat` | none | IPD-KIR's kir.dat of release KIR_DB_RELEASE (versions.env), which `scripts/setup.sh --kir-data` installs as kir/IPD-KIR_&lt;release&gt;/kir.dat. Required with --kir. |
| `--hla_qc` | `false` | With hla_typing: also type HLA-A, HLA-B and HLA-C with OptiType (OPTITYPE_IMAGE), a second caller, and write its two-field concordance with T1K. HLA-A and HLA-B then reach PharmCAT only when the two callers agree. |
| `--cyrius_install` | none | The Cyrius install `scripts/setup.sh --cyrius` makes (tools/cyrius-&lt;version&gt;), hash-locked by scripts/cyrius-constraints.txt. Required with cyrius. Cyrius is under the PolyForm Strict licence 1.0.0: non-commercial use only. |
| `--parascopy_data` | none | Parascopy's GRCh38 homology table and 1000 Genomes models, which `scripts/setup.sh --parascopy-data` installs as reference/parascopy-&lt;version&gt;. Required with parascopy. |
| `--yleaf_data` | none | Yleaf's data folder (GRCh38 marker positions and haplogroup tree), which `scripts/setup.sh --yleaf-data` installs as reference/yleaf-&lt;version&gt;/data: the Yleaf image ships the code only. Its ytree/ subfolder, from the same command, holds the current YFull YTree the step also places the sample on. Required with y_haplogroup. |
| `--hpo` | none | The symptoms of one person as Human Phenotype Ontology ids, comma separated (HP:0001250,HP:0001263; look them up at https://hpo.jax.org). Setting it runs exomiser (step 48), which ranks the genes and variants that could explain them; it needs --exomiser_data. |
| `--exomiser_data` | none | Exomiser's data folder, which `scripts/setup.sh --exomiser-data` installs as reference/exomiser-&lt;release&gt; (its &lt;release&gt;_hg38 and &lt;release&gt;_phenotype folders). Required with --hpo. |
| `--family_name` | none | With --pedigree: the family to run when the file is a families file (a pedigree with a leading family column, assets/families.example.tsv) that holds more than one family, and the family's folder under --outdir. Without it the folder is family/. |
| `--exomiser_proband` | none | The sample the --hpo terms describe. Needed when the samplesheet has more than one sample. |
| `--exomiser_assembly` | `hg38` | The assembly of the Exomiser data. Keep hg38: hg19 is there for Exomiser's cut-down test release, which CI runs. One of: `hg38`, `hg19`. |
| `--medications` | none | Opt-in: a plain text list of the drugs taken, one per line (the generic name first; '#' starts a comment). With cpic, a strong or moderate CYP2D6 inhibitor on it gives CYP2D6 a phenoconverted phenotype by CPIC's activity-score rule beside the genotype one, and a CYP2C19 inhibitor gets a note (assets/pgx_inhibitors.tsv). Never set by default. |
| `--archaic_data` | none | Step 43's table: the folder setup.sh --archaic builds (reference/archaic/&lt;ARCHAIC_DATA_VERSION&gt;): the Altai and Vindija Neanderthal and Altai Denisovan genotypes with 1000 Genomes frequencies, lifted to GRCh38. Required with archaic. |
| `--archaic_pop` | `ALL` | archaic: the 1000 Genomes group whose allele frequencies IBDmix scores your genome against. One of: `ALL`, `AFR`, `AMR`, `EAS`, `EUR`, `SAS`. |
| `--consumer_sites` | none | The 23andMe v5 site table (GRCh37 and GRCh38 positions), which `scripts/setup.sh --consumer-sites` builds as reference/consumer/23andme_v5_sites_dbsnp&lt;build&gt;.tsv.gz. Required with consumer_export. |
| `--megane_kmers` | none | MEGAnE's k-mer set of the reference for step 50, the .mk file with its .mi count file beside it, which `scripts/setup.sh --megane-kmers` builds as reference/megane/&lt;reference name&gt;.mk. Required with mei. |
| `--mocha_data` | none | MoChA's phasing panel and genetic map for step 51 (the 1000 Genomes high-coverage panel, MOCHA_PANEL_VERSION in versions.env), which `scripts/setup.sh --mocha-data` installs as reference/mocha-&lt;version&gt;. Required with mocha. |
| `--locityper_data` | none | Locityper's database for step 53 (counts.jf, db/ and loci.bed, built from the HPRC pangenome panel, LOCITYPER_PANEL_VERSION in versions.env, and the reference), which `scripts/setup.sh --locityper-data` installs as reference/locityper-&lt;version&gt;. Required with locityper. |
| `--locityper_bg_region` | none | The background region locityper preproc measures read depth and insert size on, as chrom:start-end. Default: Locityper's own for GRCh38 (chr17:72062001-76562000). Set it only when that region is not in the reference. |
| `--viral_ref` | none | The viral genomes of step 47 (NCBI RefSeq, VIRAL_REF_ACCESSIONS in versions.env), which `scripts/setup.sh --viral-data` installs as reference/viral-&lt;version&gt;/viral_genomes.fa. Required with viral. |
| `--chip_gff` | none | Ensembl's GRCh38 gene annotation, which `scripts/setup.sh --chip-data` installs as reference/ensembl/Homo_sapiens.GRCh38.&lt;release&gt;.chr.gff3.gz: chip takes each CHIP gene's MANE Select transcript from it. Required with chip. |
| `--chip_germline_resource` | none | chip: gnomAD's AF-only germline resource for Mutect2 (somatic/af-only-gnomad.hg38.vcf.gz, .tbi beside it; optional). Without it inherited variants are left out by their VAF only. |
| `--chip_pon` | none | chip: a Mutect2 panel of normals (somatic/1000g_pon.hg38.vcf.gz, .tbi beside it; optional). |
| `--parascopy_population` | `EUR` | parascopy: the 1000 Genomes population whose model parameters are used. One of: `AFR`, `AMR`, `EAS`, `EUR`, `SAS`. |
| `--parascopy_depth_bed` | none | parascopy: background windows (BED, all of one size) for a BAM that covers only part of the genome; default Parascopy's own GRCh38 windows. |
| `--paraphase_regions` | none | paraphase: the Paraphase regions to run, comma separated (for example smn1,hba,rccx,GBA,pms2,strc,ncf1,OTOA); empty runs all of its GRCh38 regions. |
| `--slivar_bin` | none | Deprecated and ignored: slivar now runs from the pinned image (SLIVAR_IMAGE in versions.env). Accepted with a warning for this release only. |
| `--cnvpytor_resources` | none | Directory with the pinned CNVpytor GC and mask .pytor files (see docs/00-reference-setup.md). Required by cnvpytor: the 1.3.2 image ships without them. |
| `--annotsv_annotations` | none | AnnotSV annotation directory (contains Annotations_Human/), built once with AnnotSV's INSTALL_annotations.sh. Required by annotsv: the AnnotSV image ships no annotation data. |
| `--delly_exclude` | none | Optional Delly exclude map, passed as `delly sr -x` (for example Delly's excludeTemplates/human.hg38.excl.tsv). |
| `--manta_call_regions` | none | Optional bgzipped BED of the regions Manta calls on, passed as `configManta.py --callRegions`; its .tbi must sit beside it. For example chr1-22, X and Y, which keeps Manta off the ALT, decoy and unplaced contigs. |
| `--cytoband` | none | Optional UCSC GRCh38 chromosome bands (cytoBand.hg38.txt, chr1-22, X and Y; scripts/setup.sh installs it), passed as `telomerehunter -b`. Without it TelomereHunter uses its own hg19 bands. |
| `--somalier_sites` | none | somalier's GRCh38 sites VCF (sites.hg38.vcf.gz; scripts/setup.sh installs it under reference/somalier/). Required by sample_qc. |
| `--verifybamid2_panel` | none | Folder holding VerifyBamID2's marker panel, the .UD, .mu and .bed files of one prefix (scripts/setup.sh installs the 1000 Genomes 100k GRCh38 panel under reference/verifybamid2/). Required by sample_qc. |
| `--freemix_warn` | `0.03` | VerifyBamID2's FREEMIX (estimated share of reads from another person) above which sample_qc reports possible contamination. A warning only; it never stops the run. |

## Resource limits

Set maximum resource limits for the pipeline.

| Parameter | Default | What it does |
|---|---|---|
| `--max_cpus` | `16` | Maximum number of CPUs that can be requested for any single job. |
| `--max_memory` | `64.GB` | Maximum amount of memory that can be requested for any single job. |
| `--max_time` | `48.h` | Maximum amount of time that can be requested for any single job. |
