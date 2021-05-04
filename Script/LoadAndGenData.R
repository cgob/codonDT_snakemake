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
  ############ Generate RData file with every cds position if the file doesn't exist ########################
  if(!file.exists(paste0(cds_path,".RData"))){
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

