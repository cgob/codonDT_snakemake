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
import os
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
RIBO_SAMPLES=[x for x in SAMPLES if "RIBO" in x]

#### UMI / PCR-deduplication settings ####
# When enabled, raw fastq go through adapter trimming + UMI deduplication
# (Script/DedupUMI.pl) before STAR. When disabled the merged fastq is mapped
# directly, so the SRA route keeps working unchanged.
UMI = config.get('umi', {}) or {}
UMI_ENABLED = bool(UMI.get('enabled', False))

#### Parallelisation of the read counting ####
# countreads is the pipeline bottleneck: one samtools view per transcript.
# It is sharded by contig, which is exact rather than approximate (see the
# comment in Script/CountingFullSeq_Apos.pl). Shard membership is computed by
# Script/make_shards.py; only the shard COUNT has to be known at parse time.
N_SHARDS = int(config.get('count_shards', 25))
SHARDS = [str(i) for i in range(N_SHARDS)]

#### A-site initiation-peak search window ####
# Metagene coordinates (0 = last nt of the AUG). Sign must match A_site_end;
# find_A_pos.R errors out if it does not. Falls back to the historical
# monosome defaults when the key is absent from config.yaml.
A_SITE_WINDOW = config.get('A_site_window')
if not A_SITE_WINDOW:
    A_SITE_WINDOW = [-20, -10] if config['A_site_end'] == '5p' else [10, 20]
A_SITE_WINDOW = [int(A_SITE_WINDOW[0]), int(A_SITE_WINDOW[1])]

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


def get_merged_fastq(wildcards):
    # Pre-merged libraries are often kept gzipped (and may be symlinks to a
    # read-only location, so they cannot be decompressed in place). cutadapt
    # reads .gz directly, so hand the compressed file straight to dedup_umi.
    merged = "Data/Raw/" + wildcards.sample + ".fastq.merge"
    if not os.path.exists(merged) and os.path.exists(merged + ".gz"):
        return merged + ".gz"
    return merged


def get_star_fastq(wildcards):
    if UMI_ENABLED:
        return "Data/Raw/" + wildcards.sample + ".fastq.dedup"
    return "Data/Raw/" + wildcards.sample + ".fastq.merge"


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
        expand("Data/Fit/{sample}_plot_{pair}.pdf", sample=SAMPLES , pair= pair_pos), "Data/Tables/summary_flux.tsv", "Data/Tables/summary_single_DT.tsv", "Data/Tables/summary_pair_DT.tsv", expand("Data/A_site_offset/{sample}_A_site_profiles.pdf", sample=RIBO_SAMPLES)
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

##----------------------------------------------------##
##  Adapter trimming and UMI PCR-deduplication         ##
##----------------------------------------------------##
## Only part of the DAG when config['umi']['enabled'] is true.
## Read layout: 5'-[UMI left][insert][UMI right][3' adapter]-3'

rule dedup_umi:
    input:
        get_merged_fastq
    output:
        "Data/Raw/{sample}.fastq.dedup"
    params:
        adapter = UMI.get('adapter', ''),
        umi_left = UMI.get('left', 0),
        umi_right = UMI.get('right', 0),
        min_insert = UMI.get('min_insert', 10)
    log:
        "Data/Raw/{sample}.dedup.log"
    shell: "perl {homedir}Script/DedupUMI.pl {input} {params.adapter} {params.umi_left} {params.umi_right} {params.min_insert} {output} 2> {log}"

##--------------------------------------##
##  STAR alignment to the genome        ##
##--------------------------------------##


rule runstar:
    input:
        fastq = get_star_fastq , genome= rules.run_index_star.output.genome
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
        A_site_end = config["A_site_end"],
        win_lo = A_SITE_WINDOW[0],
        win_hi = A_SITE_WINDOW[1]
    wildcard_constraints: sample=".*RIBO.*"   
    shell: "Rscript {homedir}Script/find_A_pos.R {input.A_site} {params.L1} {params.L2} {params.A_site_end} {output.tsv} {output.pdf} {params.win_lo} {params.win_hi}"


##------------------------------------------------------##
##  A-site profile diagnostic, all lengths, wide window ##
##------------------------------------------------------##
## Independent of L1/L2 and of A_site_window: it plots every length present
## in the pile-up over a wide range, with no peak selection. That is what is
## needed to CHOOSE the size window and the search window in the first place -
## find_A_pos.R only ever plots L1:L2, so a population outside the current
## thresholds is invisible there.

rule plot_A_site_profiles:
    input:
        A_site = "Data/A_site_offset/{sample}_A_site_pos.tsv"
    output:
        pdf = "Data/A_site_offset/{sample}_A_site_profiles.pdf"
    params:
        xlo = -100,
        xhi = 60,
        win_lo = A_SITE_WINDOW[0],
        win_hi = A_SITE_WINDOW[1],
        A_site_end = config["A_site_end"]
    wildcard_constraints: sample=".*RIBO.*"
    shell: "Rscript {homedir}Script/plot_A_site_profiles.R {input.A_site} {output.pdf} {params.xlo} {params.xhi} {params.win_lo} {params.win_hi} {params.A_site_end}"


##--------------------------------------##
## Read counting and CDS position       ##
##--------------------------------------##

rule make_shards:
    input:
        cds = rules.download_ensembl_cds.output.cds
    output:
        expand("Data/Counting/shards/{shard}.txt", shard=SHARDS)
    params:
        n = N_SHARDS,
        outdir = "Data/Counting/shards"
    shell: "python3 {homedir}Script/make_shards.py {input.cds} {params.n} {params.outdir}"

rule countreads_shard:
    input:
        bam = "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam",
	bam_index = "Data/Mapping/{sample}Aligned.sortedByCoord.out.bam.bai",
	cds = rules.download_ensembl_cds.output.cds,
	A_site_pos = get_AsiteRNA_from_ribo,
	shard = "Data/Counting/shards/{shard}.txt"
    output:
        temp("Data/Counting/{sample}.shard{shard}.count")
    params:
        L1 = config["L1"],
        L2 = config["L2"],
        STRAND = config["library"],
        A_site_end = config["A_site_end"]
    wildcard_constraints: sample="[^.]+", shard="[0-9]+"
    shell: "perl {homedir}Script/CountingFullSeq_Apos.pl {input.bam} {params.L1} {params.L2} {params.STRAND} {params.A_site_end} {input.cds} {input.A_site_pos} {output} {input.shard}"

rule countreads:
    input:
        expand("Data/Counting/{{sample}}.shard{shard}.count", shard=SHARDS)
    output:
        "Data/Counting/{sample}.count"
    wildcard_constraints: sample="[^.]+"
    shell: "cat {input} > {output}"

##--------------------------------------##
## Parse CDS for ref. in the fit        ##
##--------------------------------------##

rule parsecds:
    input:
        cds = rules.download_ensembl_cds.output.cds
    output:
        parse_cds = refdir + spec + '/ensembl.cds.parse.fa'
    shell: "perl {homedir}Script/CdsAllSeq.pl {input} >  {output}"

##--------------------------------------------------##
## Cache of every in-frame 40-codon CDS window       ##
##--------------------------------------------------##
## Its own rule on purpose: LoadAndGenData.R used to build this behind an
## "if (!file.exists)" guard, so the parallel loaddata jobs all raced to write
## the same file and readers died on a half-written RData.

rule prep_cds_rdata:
    input:
        parse_cds = rules.parsecds.output.parse_cds
    output:
        rdata = refdir + spec + '/ensembl.cds.parse.fa.RData'
    shell: "Rscript {homedir}Script/PrepCdsRData.R {input.parse_cds}"

##--------------------------------------##
## Load count file and gen. matrix      ##
##--------------------------------------##

rule loaddata:
    input:
        count_2="Data/Counting/{sample}.count", 
        parse_cds=rules.parsecds.output.parse_cds,
        cds_rdata=rules.prep_cds_rdata.output.rdata
    output:
        "Data/Counting/{sample}_ncount.RData"
    params: filter_1 = config['filter_1']
    wildcard_constraints: sample="[^.]+"
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
