snake_orthofinder (snakefile):
inputs: genomes and gffs of the 13 aculeate species (see supplementary table for file locations; P. rugosus used P. barbatus as a proxy)
	1.identifies longest isoform using gene annotation files and genomes for each of the 13 aculeate species (AGAT)
	2.identify orthologs across the 13 species with the longest isoform output (OrthoFinder)
snake_alignment (snakefile):
inputs: longest isoforms of the single copy orthologs found in OrthoFinder
	1. aligns the dna sequence for each single copy orthoglog from orthofinder output (MACSE)
	2. procure go terms for each single copy ortholog (interproscan)
snake_busted (snakefile):
inputs: dna alignments of the single copy orthologs and ortholog gene tree outputs from orthoFinder
snakefile 1. detect signature of positive selection on single copy orthologs for each single copy ortholog (HyPhy: BUSTED)

trait_dat.csv: predictors of recombination rate in aculeates and their recombination rate

pgls_analysis.R:
inputs: trait_dat.csv
	phylogenetic least squares of each combination of predictor of genome-wid recombination rate (average number of crossovers per chromosome)
		calculates AIC and leave-one-out cross validation r-squared values for each model
rer_converge.R:
inputs: trait_dat.csv, all single copy ortholog multiple sequence alignments
	runs rerconverge to look for evolutionary rate relationships between orthologs and predictors of GWRR and GWRR

fdr_rerconverge_busted.R:
inputs: rerconverge outputs for each single copy ortholog from rer_converge.R, busted results of each single copy ortholog
	finds genes which are significant for BUSTED and rerconverge after BH p<value adjustment for both rerconverge and BUSTED

gene_set_enrichment_analysis.R:
inputs: rerconverge outputs for each single copy ortholog
	gene set enrichment analyses of genes with a positive evolutioary rate relationship for each predictor of GWRR and GWRR,
 	using the rerconverge output p-value to rank gene importance

