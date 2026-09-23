#!/usr/bin/env Rscript
######################


args = commandArgs(trailingOnly=TRUE)
LoadAndGenData = function(input_file,output_file,cds_path, filter_1){
  ########  
  require("methods")
  require("MatrixModels")
  require("Matrix")
  require('MASS') 
  require('data.table')
  require('stringi') 
  #######
  filter_1 = as.numeric(filter_1)
  ############ Reference of every CDS position: built once by PrepCdsRData.R
  ############ (its own snakemake rule, so parallel loaddata jobs cannot race).
  if(!file.exists(paste0(cds_path,".RData"))){
    stop("missing ", paste0(cds_path, ".RData"), " - run PrepCdsRData.R first")
  }
############################################################################################################  

############################ Load reference CDS and count RFP data. Merging of the two data.frame #############

  load(paste0(cds_path,".RData"))
  rfp=fread(input_file,colClasses=c('character','character','numeric','numeric','character'),stringsAsFactors=F)  
  rfp=as.data.frame(rfp)
  rfp=rfp[sapply(rfp$V2,nchar)==120,]
  rownames(rfp)=paste0(rfp$V1,rfp$V2)
  rfp=rfp[rownames(rfp)%in%rownames(ncount),]
  ncount$rfp =0
  ncount$wig=NA
  ncount[rownames(rfp),'rfp']=rfp[,3]
  ncount[rownames(rfp),'wig']=rfp[,4]
  ncount[rownames(rfp),'length']=rfp[,5]
  rm(rfp)
  gc()
  rownames(ncount)=1:nrow(ncount)
############################################################################################################  

##################### Selection of genes and CDS positions stastifying our criterion  ################################################
  ss=split(1:nrow(ncount),ncount$V1)
  gg=sapply(ss,function(x) sum(ncount$rfp[x]))
  
  names(ncount)[1]='gene'
  names(ncount)[42]='count.rfp'
  
  ss=split(1:nrow(ncount),ncount$gene)
  gg=sapply(ss,function(x) length(ncount$count.rfp[x][ncount$count.rfp[x]!=0])/length(ncount$count.rfp[x]))
  ncount=ncount[ncount$gene%in%names(gg[gg>0.05]),] # Genes with less than 5% positions covered are removed
  
  ncount.nn=subset(ncount,count.rfp!=0)
  ss=split(1:nrow(ncount.nn),ncount.nn$gene)
  gg=sapply(ss,length)
  ncount=ncount[ncount$gene%in%names(gg[gg>5]),] # Genes with less than 5 positions are removed
  
  ss=split(1:nrow(ncount),ncount$gene)
  gg=sapply(ss,function(x) sum(ncount$count.rfp[x]))
  ncount=ncount[ncount$gene%in%names(gg[gg> filter_1]),] # Genes with less than 100 reads are removed
  rm(ncount.nn)
  
  for(i in 1:40){
    ncount[,i+1]=as.factor(as.character(unlist(ncount[,i+1]))) # Preparation of the data.frame for the fit with the right number of factors.
  }
  
  save(ncount,file=output_file)
}
###################################################################################################################################

LoadAndGenData(args[1], args[2], args[3], args[4] )

