#!/usr/bin/env python

###########################################################################
## Author: Gobet Cédric
## Email: cedric.gobet@epfl.ch
## Date: 25/03/2019
## 
## The purpose of this pipeline is to infer gene flux and codon dwell times
# from ribosome profiling data.
## The pipeline is still in active development.
###########################################################################


import pandas as pd
import glob 
import numpy as np

#### Load configuration and sample sheet ####
configfile: "config.yaml"

df = pd.read_csv(config["samples"],sep="\t")
df.SAMPLES = df.SAMPLES + "_" + df.Type
SRR = df.SRR
spec = config['species']
types = df.Type

df = df.set_index('SRR') 

GTF_URL = config[spec]['gtf']
CDS_URL = config[spec]['cds']
GENOME_URL = config[spec]['genome']

homedir=config['homedir']
workdir: config['workdir'] 

SAMPLES=df.SAMPLES.unique()

### Function definition ###
def get_sample_sra_rel(wildcards):
    return "Data/Raw/" + df[df['SAMPLES']== wildcards.sample].index + ".fastq"

def get_rna_from_ribo(wildcards):
    sample_name = str(wildcards.sample).split('_RIBO')[0]
    sample_rna = sample_name + "_RNA"

    if df[df['SAMPLES'] == sample_rna].size == 0:
        return "Data/Counting/" + str(wildcards.sample) + "_ncount.RData"
    else:
        return "Data/Fit/" + sample_rna + "_fit_" + str(wildcards.pair) + ".RData"

##--------------------------------------##
##  Target rule                         ##
##--------------------------------------##
rule all:
     input:
        expand("Data/Fit/{sample}_plot_{pair}.pdf", sample=SAMPLES , pair=["24:25","25:26","23:24"])

##--------------------------------------##
##  Download gtf from Ensembl           ##
##--------------------------------------##

rule download_ensembl_gtf:
    output: 
        gtf = 'references/' + spec + '/ensembl.gtf'
    params: 
        url = GTF_URL
    shell: "wget -O {output.gtf}.gz {params.url}; gunzip {output.gtf}.gz"

##--------------------------------------##
##  Download cds from Ensembl           ##
##--------------------------------------##

rule download_ensembl_cds:
    output: 
        cds = 'references/' + spec + '/ensembl.cds.fa'
    params: 
        url = CDS_URL
    shell: "wget -O {output.cds}.gz {params.url}; gunzip {output.cds}.gz"

##--------------------------------------##
##  Download dna from Ensembl           ##
##--------------------------------------##

rule download_ensembl_genome:
    output: 
        genome = 'references/' + spec + '/ensembl.genome.fa'
    params: 
        url = GENOME_URL
    shell: "wget -O {output.genome}.gz {params.url}; gunzip {output.genome}.gz"

##--------------------------------------##
##  Generate STAR genome index          ##
##--------------------------------------##

rule run_index_star:
    input: fa = rules.download_ensembl_genome.output.genome,
       gtf = rules.download_ensembl_gtf.output.gtf
    output: genome = directory('references/' + spec + '/genome')
    shell: "module load gcc/7.4.0; module load star/2.7.0e; mkdir {output.genome}; STAR --runMode genomeGenerate --runThreadN 12 --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeDir {output.genome}"

##----------------------------------------------------------------##
## Download sra files from GEO  with fastq-dump and output fastq  ##
##----------------------------------------------------------------##

rule downloadSRA:
    output:
         "Data/Raw/{srr}.fastq"
    shell: "module add sra-toolkit/2.9.6; fasterq-dump {wildcards.srr} -O Data/Raw/"


##--------------------------------------##
## Trimm Fastq in case of UMI           ##
##--------------------------------------##

#rule trimfastq:
 #   input:
  #      ""
   # output:
#	""
 #   params: adapter = config["adapter"]
  #  shell: "fastx_clipper -Q33 -a {params.adapter}  -l 25  -i {input} -o {output}";


##--------------------------------------##
## Merge fastq from the same sample     ##
##--------------------------------------##
rule mergefastq:
    input: 
        get_sample_sra_rel
    output:
        "Data/Raw/{sample}.fastq.merge"
    shell: "cat {input} > {output}"

##--------------------------------------##
##  STAR alignment to the genome        ##
##--------------------------------------##


rule runstar:
    input:
        fastq = "Data/Raw/{sample}.fastq.merge" , genome= rules.run_index_star.output.genome
    output:
        "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam"
    params: star_params = "--outSAMtype BAM SortedByCoordinate --seedSearchStartLmax 15 --limitBAMsortRAM 61000000000",
            star_adapter = config["adapter"]
    shell: "module load gcc/7.4.0;  module load star/2.7.0e; STAR --genomeDir {input.genome} {params.star_params} --clip3pAdapterSeq {params.star_adapter} --outFileNamePrefix Data/Mapping/{wildcards.sample} --readFilesIn {input.fastq} --runThreadN 12"

##--------------------------------------##
##  BAM files indexing                  ##
##--------------------------------------##

rule samindex:
    input:
        "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam"
    output:
        "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam.bai"
    shell: "module load gcc/7.4.0; module load samtools/1.9; samtools index {input}"


##--------------------------------------##
##  Fragment size distribution          ##
##--------------------------------------##

rule sizedistrib:
    input:
        bam ="Data/Mapping/{sample}Aligned.sortedByCoord.out.bam", bam_index="Data/Mapping/{sample}Aligned.sortedByCoord.out.bam.bai"
    output:
        "Data/Mapping/{sample}_fragment_size.txt"
    shell: "module load gcc/7.4.0; module load samtools/1.9; perl {homedir}Script/SizeDistrib.pl {input.bam} > {output}"


##--------------------------------------##
##  Plot size distribution              ##
##--------------------------------------##

rule plotsizedistrib:
    input:
    	"Data/Mapping/{sample}_fragment_size.txt"
    output:
        "Data/Mapping/{sample}_fragment_size.pdf"
    shell: "module add gcc/7.4.0;  module add openblas/0.3.6-openmp; module add r/3.6.0; Rscript {homedir}Script/PlotSizeDistrib.R {input} {output}"


##--------------------------------------##
## Read counting and CDS position       ##
##--------------------------------------##

rule countreads:
    input:
        bam ="Data/Mapping/{sample}Aligned.sortedByCoord.out.bam", bam_index="Data/Mapping/{sample}Aligned.sortedByCoord.out.bam.bai", cds=rules.download_ensembl_cds.output.cds
    output:
        "Data/Counting/{sample}.count" 
    params: L1 = config["L1"],
	    L2 = config["L2"],
	    STRAND = config["library"],
    shell: "module load gcc/7.4.0; module load samtools/1.9; perl {homedir}Script/CountingFullSeq.pl {input.bam} {params.L1} {params.L2} {params.STRAND} {input.cds} {output}"

##--------------------------------------##
## Parse CDS for ref. in the fit        ##
##--------------------------------------##

rule parsecds:
    input:
        cds = rules.download_ensembl_cds.output.cds
    output:
        parse_cds = 'references/' + spec + '/ensembl.cds.parse.fa'
    shell: "module load gcc/7.4.0; module load samtools/1.9; perl {homedir}Script/CdsAllSeq.pl {input} >  {output}"

##--------------------------------------##
## Load count file and gen. matrix      ##
##--------------------------------------##

rule loaddata:
    input:
        count_2="Data/Counting/{sample}.count", 
        parse_cds=rules.parsecds.output.parse_cds
    output:
        "Data/Counting/{sample}_ncount.RData"
    shell: "module add gcc/7.4.0;  module add openblas/0.3.6-openmp; module add r/3.6.0; Rscript {homedir}Script/LoadAndGenData.R {input.count_2} {output} {input.parse_cds} "

##--------------------------------------##
## Make GLM fit (DT and flux)           ##
##--------------------------------------##

rule makefit:
    input:
        "Data/Counting/{sample}_ncount.RData"
    output:
        fit="Data/Fit/{sample}_fit_{pair}.RData", pred="Data/Fit/{sample}_fit_{pair}.RData.pred", cor="Data/Fit/{sample}_fit_{pair}.RData.cor"
    params: codon = '1:40',
	    mode = 'simple'
    wildcard_constraints: sample=".*RNA.*"
    shell: "module add gcc/7.4.0; module add openblas/0.3.6-openmp; module add r/3.6.0; Rscript {homedir}Script/MakeFit.R {input} NULL {output.fit} {params.codon} {wildcards.pair} {params.mode}"


##------------------------------------------------------------------##
## Make GLM fit (DT and flux) with RNA-seq prediciton as an offset  ##
##------------------------------------------------------------------##

rule makefit_combined:
    input:
        ncount="Data/Counting/{sample}_ncount.RData",
	rnafit= get_rna_from_ribo
    output:
        fit="Data/Fit/{sample}_fit_{pair}.RData", pred="Data/Fit/{sample}_fit_{pair}.RData.pred", cor="Data/Fit/{sample}_fit_{pair}.RData.cor"
    params: codon = '1:40',
            mode = 'combined'
    wildcard_constraints: sample=".*RIBO.*"
    shell: "module add gcc/7.4.0; module add openblas/0.3.6-openmp; module add r/3.6.0; Rscript {homedir}Script/MakeFit.R {input.ncount} {input.rnafit} {output.fit} {params.codon} {wildcards.pair} {params.mode}"


##--------------------------------------##
## Compute p-value and rescale coef.    ##
##--------------------------------------##

rule coepval:
    input:
        "Data/Fit/{sample}_fit_{pair}.RData"
    output:
        "Data/Fit/{sample}_coe_pval_{pair}.RData"
    shell: "module add gcc/7.4.0; module add openblas/0.3.6-openmp; module add r/3.6.0; Rscript {homedir}Script/CoeAndPVal.R {input} {output}"


##-----------------------------------------------##
## Plot single, pair heatmaps and scatterplot    ##
##-----------------------------------------------##

rule heatmap:
    input:
        coe="Data/Fit/{sample}_coe_pval_{pair}.RData", size_2="Data/Mapping/{sample}_fragment_size.pdf"
    output:
        "Data/Fit/{sample}_plot_{pair}.pdf"
    shell: "module add gcc/7.4.0; module add openblas/0.3.6-openmp; module add r/3.6.0; Rscript {homedir}Script/PlotHeatmaps.R {input.coe} {output}"


