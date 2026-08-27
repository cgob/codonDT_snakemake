#!/usr/bin/env Rscript
## Build the cached reference of every in-frame 40-codon CDS window.
##
## This used to live inside LoadAndGenData.R behind an `if (!file.exists(...))`
## guard. With one loaddata job per sample they all start together, all see the
## cache missing, and all write the same path - readers then hit a half-written
## file and die with "ReadItem : type inconnu 0". It is a rule of its own now,
## so snakemake builds it exactly once before any loaddata runs.
##
## Usage: Rscript PrepCdsRData.R <ensembl.cds.parse.fa>
##        writes <ensembl.cds.parse.fa>.RData

args = commandArgs(trailingOnly=TRUE)
cds_path = args[1]

require("methods")
require("MatrixModels")
require("Matrix")
require('MASS')
require('data.table')
require('stringi')

  cod= fread(cds_path,colClasses=c('character','character'),stringsAsFactors=F,header=F) # Load CDS 
  cod=as.data.frame(cod)
  cod=cod[sapply(cod$V2,nchar)==120,] # Remove sequence smaller than 120 nucleotides
  N=grep('N',cod$V2)  # Remove sequences containing "N"
   if(length(N)!=0){
   cod=cod[-N,]
   }
  
  DD=matrix(NA,nrow=nrow(cod),ncol=40) #Initialize matrix in which sequences are splitted in codons
  
  for(i in 1:nrow(cod)){
    b=unlist(stri_extract_all_regex(cod$V2[i], '.{1,3}')) # Split sequences in 3-nucleotides codons
    DD[i,]=b
  }
  
  ncount=data.frame(cod,DD)
  aa=unlist(sapply(1:40, function(x) grep("TAG|TGA|TAA|N",ncount[,x+2]))) # remove sequences with stop codons and N
  if(length(aa)!=0){
    ncount=ncount[-aa,]
  }
  
  #Remove non-unique 30-mers sequences
  ncount$seq=apply(ncount[,22:32],1,function(x) paste0(as.character(unlist(x)), collapse=""))
  ss.1=split(1:nrow(ncount),ncount$seq)
  ss.2=sapply(ss.1,length)
  ss.3=names(ss.2[ss.2>1])
  ncount=ncount[!ncount$seq%in%ss.3,]
  rownames(ncount)=paste0(ncount$V1,ncount$V2)
  ncount=ncount[,-match(c("V2","seq"),names(ncount))]
  save(ncount,file=paste0(cds_path,".RData"))
