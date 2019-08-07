# codonDT_snakemake
Snakemake pipeline to infer single-codon and codon-pair dwell times as well as gene flux from ribosome profiling data using generalized linear model (GLM) with negative binomial noise (Gobet et al., BioRxiv, 2019). ! Pipeline under development for generalization to any kind of ribosome profiling datasets !

## Pipeline

1. Download the genome, CDS and gtf files from ENSEMBL according to the species defined in the configuration file
2. Build genome index using STAR 
3. Download SRA run files (SRR token) specified in the sample spreadsheet from GEO database
4. Convert SRA to Fastq files
5. Merge Fastq files run from the sample according to the the sample spreadsheet definition
6. STAR alignement to genome with inline adapter clipping (defined in the configuration file)
7. Bam files indexing 
8. Read counts and CDS positions are retrieved
9. Parse the downloaded CDS file to be used as a reference in the GLM fit
10. Load the CDS parsed and read count files to generate the matrix for the fit
11. Select gene and positions. Make the fit the glm4 function.
12. Compute coefficients p-value and rescale the coefficients according to our convention (see method section in the paper)

## Authors
Cédric Gobet
