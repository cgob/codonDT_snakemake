#!/usr/bin/env python

###########################################################################
## Author: Gobet Cédric
## Email: cedric.gobet@epfl.ch
## Date: 03/05/2021
## 
## The purpose of the Ribo-DT is to infer gene flux and codon dwell times
# from ribosome profiling data.
## The pipeline is under continuous development.
###########################################################################


import pandas as pd
import glob 
import numpy as np
import re

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

refdir=config['refdir']
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


def get_AsiteRNA_from_ribo(wildcards):
    sample_name = str(wildcards.sample)
    return "Data/A_site_offset/" + re.sub("RNA", "RIBO", sample_name) + "_A_site_pos_inferred.tsv"

### define pair interaction to compute (24 = E site, 25 = P site, 26 = A site)

pair_pos = ['24:25','25:26','24:26']

##--------------------------------------##
##  Target rule                         ##
##--------------------------------------##
rule all:
     input:
        expand("Data/Fit/{sample}_plot_{pair}.pdf", sample=SAMPLES , pair= pair_pos), "Data/Tables/summary_flux.tsv", "Data/Tables/summary_single_DT.tsv", "Data/Tables/summary_pair_DT.tsv"
##--------------------------------------##
##  Download gtf from Ensembl           ##
##--------------------------------------##

rule download_ensembl_gtf:
    output: 
        gtf = refdir + spec + '/ensembl.gtf'
    params: 
        url = GTF_URL
    shell: "wget -O {output.gtf}.gz {params.url}; gunzip {output.gtf}.gz"

##--------------------------------------##
##  Download cds from Ensembl           ##
##--------------------------------------##

rule download_ensembl_cds:
    output: 
        cds = refdir + spec + '/ensembl.cds.fa'
    params: 
        url = CDS_URL
    shell: "wget -O {output.cds}.gz {params.url}; gunzip {output.cds}.gz"

##--------------------------------------##
##  Download dna from Ensembl           ##
##--------------------------------------##

rule download_ensembl_genome:
    output: 
        genome = refdir + spec + '/ensembl.genome.fa'
    params: 
        url = GENOME_URL
    shell: "wget -O {output.genome}.gz {params.url}; gunzip {output.genome}.gz"

##--------------------------------------##
##  Generate STAR genome index          ##
##--------------------------------------##

rule run_index_star:
    input: fa = rules.download_ensembl_genome.output.genome,
       gtf = rules.download_ensembl_gtf.output.gtf
    output: genome = directory(refdir + spec + '/genome')
    shell: "mkdir {output.genome}; STAR --runMode genomeGenerate --runThreadN 12 --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeDir {output.genome}"

##----------------------------------------------------------------##
## Download sra files from GEO  with fastq-dump and output fastq  ##
##----------------------------------------------------------------##

rule downloadSRA:
    output:
         "Data/Raw/{srr}.fastq"
    shell: "fasterq-dump {wildcards.srr} -O Data/Raw/"


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
    shell: "STAR --genomeDir {input.genome} {params.star_params} --clip3pAdapterSeq {params.star_adapter} --outFileNamePrefix Data/Mapping/{wildcards.sample} --readFilesIn {input.fastq} --runThreadN 12"

##--------------------------------------##
##  BAM files indexing                  ##
##--------------------------------------##

rule samindex:
    input:
        "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam"
    output:
        "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam.bai"
    shell: "samtools index {input}"


##--------------------------------------##
##  Fragment size distribution          ##
##--------------------------------------##

rule sizedistrib:
    input:
        bam ="Data/Mapping/{sample}Aligned.sortedByCoord.out.bam", bam_index="Data/Mapping/{sample}Aligned.sortedByCoord.out.bam.bai"
    output:
        "Data/Mapping/{sample}_fragment_size.txt"
    shell: "perl {homedir}Script/SizeDistrib.pl {input.bam} > {output}"


##--------------------------------------##
##  Plot size distribution              ##
##--------------------------------------##

rule plotsizedistrib:
    input:
    	"Data/Mapping/{sample}_fragment_size.txt"
    output:
        "Data/Mapping/{sample}_fragment_size.pdf"
    shell: "Rscript {homedir}Script/PlotSizeDistrib.R {input} {output}"

##--------------------------------------##
##  Compute pile_up plot at start codon ##
##--------------------------------------##

rule pileupplot:
    input:
        bam ="Data/Mapping/{sample}Aligned.sortedByCoord.out.bam", 
	bam_index="Data/Mapping/{sample}Aligned.sortedByCoord.out.bam.bai",
        gtf = rules.download_ensembl_gtf.output.gtf
    output:
        "Data/A_site_offset/{sample}_A_site_pos.tsv"
    params:
        strand  = config["library"],
        A_site_end = config["A_site_end"]
    wildcard_constraints: sample=".*RIBO.*"   
    shell: "perl {homedir}Script/compute_profile_all.pl {input.gtf} {input.bam} {params.strand} {params.A_site_end} > {output}"

##--------------------------------------------##
##  Compute A site position from pile_up plot ##
##--------------------------------------------##

rule findAsite:
    input:
        A_site="Data/A_site_offset/{sample}_A_site_pos.tsv"
    output:
        tsv="Data/A_site_offset/{sample}_A_site_pos_inferred.tsv",
        pdf="Data/A_site_offset/{sample}_A_site_pos_inferred.pdf"
    params:
        L1 = config["L1"],
        L2 = config["L2"],
        A_site_end = config["A_site_end"]
    wildcard_constraints: sample=".*RIBO.*"   
    shell: "Rscript {homedir}Script/find_A_pos.R {input.A_site} {params.L1} {params.L2} {params.A_site_end} {output.tsv} {output.pdf}"


##--------------------------------------##
## Read counting and CDS position       ##
##--------------------------------------##

rule countreads:
    input:
        bam = "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam", 
	bam_index = "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam.bai",
	cds = rules.download_ensembl_cds.output.cds,
	A_site_pos = get_AsiteRNA_from_ribo
    output:
        "Data/Counting/{sample}.count" 
    params: 
        L1 = config["L1"],
        L2 = config["L2"],
        STRAND = config["library"],
        A_site_end = config["A_site_end"]
    shell: "perl {homedir}Script/CountingFullSeq_Apos.pl {input.bam} {params.L1} {params.L2} {params.STRAND} {params.A_site_end} {input.cds} {input.A_site_pos} {output}"

##--------------------------------------##
## Parse CDS for ref. in the fit        ##
##--------------------------------------##

rule parsecds:
    input:
        cds = rules.download_ensembl_cds.output.cds
    output:
        parse_cds = refdir + spec + '/ensembl.cds.parse.fa'
    shell: "perl {homedir}Script/CdsAllSeq.pl {input} >  {output}"

##--------------------------------------##
## Load count file and gen. matrix      ##
##--------------------------------------##

rule loaddata:
    input:
        count_2="Data/Counting/{sample}.count", 
        parse_cds=rules.parsecds.output.parse_cds
    output:
        "Data/Counting/{sample}_ncount.RData"
    params: filter_1 = config['filter_1']
    shell: "Rscript {homedir}Script/LoadAndGenData.R {input.count_2} {output} {input.parse_cds} {params.filter_1}"

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
    shell: "Rscript {homedir}Script/MakeFit.R {input} NULL {output.fit} {params.codon} {wildcards.pair} {params.mode}"


##------------------------------------------------------------------##
## Make GLM fit (DT and flux) with RNA-seq prediciton as an offset  ##
##------------------------------------------------------------------##

rule makefit_combined:
    input:
        ncount="Data/Counting/{sample}_ncount.RData",
	rnafit= get_rna_from_ribo
    output:
        fit="Data/Fit/{sample}_fit_{pair}.RData"
    params: codon = '1:40',
            mode = 'combined'
    wildcard_constraints: sample=".*RIBO.*"
    shell: "Rscript {homedir}Script/MakeFit.R {input.ncount} {input.rnafit} {output.fit} {params.codon} {wildcards.pair} {params.mode}"


##--------------------------------------##
## Compute p-value and rescale coef.    ##
##--------------------------------------##

rule coepval:
    input:
        "Data/Fit/{sample}_fit_{pair}.RData"
    output:
        "Data/Fit/{sample}_coe_pval_{pair}.RData"
    shell: "Rscript {homedir}Script/CoeAndPVal.R {input} {output}"


##-----------------------------------------------##
## Plot single, pair heatmaps and scatterplot    ##
##-----------------------------------------------##

rule heatmap:
    input:
        coe="Data/Fit/{sample}_coe_pval_{pair}.RData", size_2="Data/Mapping/{sample}_fragment_size.pdf"
    output:
        "Data/Fit/{sample}_plot_{pair}.pdf"
    params: pval_thresh = config['filter_2']
    shell: "Rscript {homedir}Script/PlotHeatmaps.R {input.coe} {output} {params.pval_thresh}"

##---------------------------------------------------------------------##
## Output table with coefficients, standard error, t-value and p-value ##
##---------------------------------------------------------------------##

rule output_table:
    input:
        expand("Data/Fit/{sample}_coe_pval_{pair}.RData", sample=SAMPLES , pair= pair_pos)
    output:
        "Data/Tables/summary_flux.tsv", "Data/Tables/summary_single_DT.tsv", "Data/Tables/summary_pair_DT.tsv"
    shell: "Rscript {homedir}Script/OutputTable.R {input}"
