# What the Pipeline Does Not Report

A report lists what was found. It cannot list everything that was not looked at. This page does: the results we leave out on purpose, and the things short reads cannot see at all. When a result on this page is missing from your report, that silence is not a normal result. It means nobody looked, or nobody could.

Each item names its source. A PMID opens at `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`.

## Not reported on purpose

Some of these are results other services sell; the rest are tools we chose not to run. We leave each out because the evidence does not support a reading for one person, or because another step already owns the locus.

- **Genotype-based diet panels.** Professional reviews found the evidence too weak to base a person's diet on their genotype (Camp and Trujillo, J Acad Nutr Diet 2014, PMID 24439821; Keathley et al., BMJ Open 2022, PMID 35193914).
- **MTHFR C677T as health advice.** Step 42 shows it as an enzyme note only. Medical geneticists advise against MTHFR testing for thrombophilia, pregnancy loss or relatives of a carrier (Hickey et al., Genet Med 2013, PMID 23288205).
- **Longevity or biological-age scores beyond APOE.** In 1 million parent lifespans, the lifespan variants act mostly through dementia, heart disease and lung cancer, and a score built from them moves the average lifespan by about five years between the top and the bottom tenth: a group average, not a forecast for one person (Timmers et al., eLife 2019, PMID 30642433). Step 42 reports APOE, folded away behind a warning.
- **Educational attainment, cognition and personality scores.** Much of what such scores predict shrinks within families, where upbringing and population structure are held fixed, so they say little about one person (Howe et al., Nat Genet 2022, PMID 35534559).
- **ACTN3 as training advice.** Step 42 reports the R577X genotype as muscle fibre type. A consensus statement says no genetic test should guide training or talent choice (Webborn et al., Br J Sports Med 2015, PMID 26582191).
- **Coffee and CYP1A2 heart risk.** The link came from one case-control study (Cornelis et al., JAMA 2006, PMID 16522833). A prospective study of up to 347,077 people found no such interaction (Zhou and Hyppönen, Am J Clin Nutr 2019, PMID 30838377). Step 42 reports caffeine breakdown and coffee intake, with no heart reading.
- **Health readings for a mitochondrial haplogroup.** Most reported links come from small case-control studies that mtDNA population structure misleads (Salas and Elson, J Genet Genomics 2015, PMID 25953355). Step 12 reports the haplogroup as maternal ancestry, and checks it for a second person's DNA.
- **AMY1 copy number as obesity risk.** In about 4,500 people whose copies were counted directly, AMY1 copy number did not go with obesity (Usher et al., Nat Genet 2015, PMID 26098870). An AMY1 copy number is a count, with no reading.
- **CCL3L1 copy number as HIV susceptibility.** The original link (Gonzalez et al., Science 2005, PMID 15637236) did not replicate in larger cohorts (Urban et al., Nat Med 2009, PMID 19812560). Step 35 gives the copy number with no reading ([copy number loci](cn-loci.md)).
- **Chronotype from one SNP.** Being a morning or evening person is spread over 351 loci, each with a tiny effect (Jones et al., Nat Commun 2019, PMID 30696823).
- **qpAdm ancestry fits on one person.** qpAdm needs chosen source and outgroup populations, and its fits depend on that choice (Harney et al., Genetics 2021, PMID 33772284). Step 44 reports a PCA placement among ancient groups instead, a different method.
- **ExpansionHunter Denovo on one sample.** It finds new expansions by comparing many genomes against each other, as an outlier test (Dolzhenko et al., Genome Biol 2020, PMID 32345345). One genome has nothing to compare with. Step 9 genotypes the known disease repeats.
- **SMNCopyNumberCaller, Gauchian, Aldy and StellarPGx.** Each locus has one owner: Parascopy (step 35) for SMN1 and SMN2, pypgx (step 32) for the star alleles. Gauchian's GBA1 calls were also found unreliable (Tayebi et al., Commun Biol 2025, PMID 40346301).
- **ClassifyCNV, peddy, Hail, arcasHLA, GATK gCNV and HLA\*LA.** Other steps cover the same ground: AnnotSV classifies CNVs (step 5), somalier checks identity and relatedness (step 33; Pedersen et al., Genome Med 2020, PMID 32664994), T1K types HLA (step 8; Song et al., Genome Res 2023, PMID 37169596) and CNVpytor calls CNVs from depth (step 18).
- **Exomiser without symptoms.** Exomiser ranks genes by how well they explain a list of HPO terms (Smedley et al., Nat Protoc 2015, PMID 26562621). Without terms the ranking means nothing, so step 48 runs only when you give them.

## Not visible on short reads

These are limits of the data, not of the pipeline. A clean report says nothing about them. Each needs another kind of test.

- **The F8 intron 22 inversion** (about 45% of severe haemophilia A). It joins two copies of a near-identical 9.5 kb repeat that lie about 500 kb apart, and a read of 150 bases falls inside one copy (Lakich et al., Nat Genet 1993, PMID 8275087).
- **The size of the FSHD1 D4Z4 array.** The array is a row of 3.3 kb units, up to about 100 of them, with a near-identical copy on chromosome 10. Short reads cannot count the units or tell the two arrays apart (Lemmers et al., Science 2010, PMID 20724583).
- **The exact size and the interruptions of a repeat allele longer than the DNA fragment.** No read spans it, so its length is an estimate and its sequence inside is unknown (Dolzhenko et al., Genome Res 2017, PMID 28887402).
- **Methylation (5mC) and imprinting epimutations.** Sequencing by synthesis reads 5mC as a plain C, and PCR copies do not carry the mark, so a standard WGS library shows no methylation. It needs bisulfite or enzymatic conversion, or long reads (Frommer et al., PNAS 1992, PMID 1542678). A deletion of 1 kb or more over an imprinting centre can still be called from depth, as any deletion can.
- **Heterodisomic uniparental disomy from one person.** Two different chromosomes from one parent look like a normal pair. Only the parents' genotypes show it (Nakka et al., Am J Hum Genet 2019, PMID 31607426). Step 11 flags isodisomy, which shows as a long run of homozygosity.
- **Proving an SMN1 "2+0" carrier.** Two SMN1 copies on one chromosome and none on the other give the same count as one on each. A copy number cannot tell them apart (Feng et al., Genet Med 2017, PMID 28125085). Step 35 says so with its count.
- **The order and identity of copies in near-identical arrays.** Short reads count copies, but cannot say which copy sits where or which variant belongs to which copy (Vollger et al., Science 2022, PMID 35357917).
- **Multi-copy and long Y-STRs.** Markers such as DYS385a/b lie in duplicated, palindromic Y sequence, and some repeats are longer than a read can span (Ballantyne et al., Am J Hum Genet 2010, PMID 20817138). Step 45 marks them not callable.
- **Variants in paralog segments with no paralog-specific sites.** Where two copies are identical over a read's length, a read cannot be placed on one copy, so a variant there cannot be given to a gene (Ebbert et al., Genome Biol 2019, PMID 31104630). Parascopy, in step 35, tells the copies apart only at the sites where they differ (Prodanov and Bansal, Nat Commun 2022, PMID 35680869; [paralogs](35-paralogs.md)).
- **Small blood clones.** At 30x a variant at a few percent VAF shows in one or two reads, too few to call. Step 49's detection limit at 30x is about 14% VAF, against the 2% that defines CHIP (Steensma et al., Blood 2015, PMID 25931582). Step 51 reports no mosaic chromosomal alteration under a cell fraction of 5%, and its read-sampling floor at the sample's depth can be higher (Loh et al., Nature 2018, PMID 29995854). This is a sampling limit, not a read-length one: deeper sequencing lowers it.

## Not read on the default reference, or not in this genome

Short reads can see these, but the reference or the alignment of one genome can hide them. The first comes back with another reference; the second is listed per genome.

- **U2AF1, KCNE1, CBS and CRYAA on chr21.** GRCh38 holds a false copy of part of chr21's long arm on its short arm, so reads from these genes get mapping quality 0 and no caller sees them (Behera et al., Genome Biol 2023, PMID 36810122). `setup.sh --masked-reference` builds a reference with the false copies masked, and a sample aligned to it reads them again ([the masked reference](00-reference-setup.md#the-masked-reference-optional)).
- **Genes inside segmental duplications and on the false duplications.** Step 16b lists, per genome, the hard genes its reads could not cover, so a "no finding" there is marked as not read (Ebbert et al., Genome Biol 2019, PMID 31104630; [gene callability](16b-mosdepth.md#gene-callability)).
