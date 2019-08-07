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


#### Load configuration and sample sheet ####
configfile: "config.yaml"

df = pd.read_csv(config["samples"],sep="\t")
SAMPLES=df.SAMPLES.unique()
SRR=df["SRR"]
spec = config['species']

df = df.set_index('SRR') 
GTF_URL=config[spec]['gtf']
CDS_URL=config[spec]['cds']
GENOME_URL=config[spec]['genome']

homedir=config['homedir']
workdir: config['workdir'] 

### Function definition ###
def get_sample_sra_rel(wildcards):
    return "Data/Raw/" + df[df['SAMPLES']== wildcards.sample].index + ".fastq.gz"


##--------------------------------------##
##  Target rule                         ##
##--------------------------------------##
rule all:
     input:
        expand("Data/Fit/{sample}_coe_pval_{pair}.RData", sample=SAMPLES, pair=["24:25","25:26","23:24"])

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

##--------------------------------------##
## Download sra files from GEO          ##
##--------------------------------------##

rule downloadSRA:
    output: 
         "Data/Raw/{srr}.sra"
    params:
         sra = lambda wildcards: str(wildcards.srr)[0:6]
    shell: "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{params.sra}/{wildcards.srr}/{wildcards.srr}.sra -O Data/Raw/{wildcards.srr}.sra"

##--------------------------------------##
## Convert SRA to Fastq                 ##
##--------------------------------------##

rule sra2fastq:
    input:
        "Data/Raw/{srr}.sra"
    output:
        "Data/Raw/{srr}.fastq.gz"
    shell: "module load sra-toolkit/2.9.6; fastq-dump.2.9.6 --gzip {input} -O Data/Raw/"


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
        "Data/Raw/{sample}.fastq.gz.merge"
    shell: "cat {input} > {output}"

##--------------------------------------##
##  STAR alignment to the genome        ##
##--------------------------------------##


rule runstar:
    input:
        fastq = "Data/Raw/{sample}.fastq.gz.merge" , genome= rules.run_index_star.output.genome
    output:
        "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam"
    params: star_params = "--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --seedSearchStartLmax 15 --limitBAMsortRAM 61000000000",
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
        count="Data/Counting/{sample}.count", 
        parse_cds=rules.parsecds.output.parse_cds
    output:
        "Data/Counting/{sample}_ncount.RData"
    shell: "module add gcc/7.4.0; module add r/3.6.0; Rscript {homedir}Script/LoadAndGenData.R {input.count} {output} {input.parse_cds} "

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
    shell: "module add gcc/7.4.0; module add r/3.6.0; Rscript {homedir}Script/MakeFit.R {input} {output.fit} {params.codon} {wildcards.pair} {params.mode}"

##--------------------------------------##
## Compute p-value and rescale coef.    ##
##--------------------------------------##

rule coepval:
    input:
        "Data/Fit/{sample}_fit_{pair}.RData"
    output:
        "Data/Fit/{sample}_coe_pval_{pair}.RData"
    shell: "module add gcc/7.4.0; module add r/3.6.0; Rscript {homedir}Script/CoeAndPVal.R {input} {output}"

