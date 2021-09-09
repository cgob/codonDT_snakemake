#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
require(reshape2)
require(Biostrings)
require(data.table)

##########################
change_cc2aa=function(cod){
  sub = as.character(AMINO_ACID_CODE[as.character(GENETIC_CODE[cod])])
  sub[is.na(sub)]='Stop'
  cod = paste(sub,cod,sep="->")
}

##########################
colr=c("#87597C", "#DF9158", "#D6ACE3", "#8CE55E", "#65E7D1", "#D042E4", "#E4D94E", "#DFD5DF", "#7C84DC", "#DB6ECB", "#6E9574", "#8CD5E0", "#78A9D4", "#7AE296", "#DED3B2", "#D3D57F", "#C0EAC9", "#7753E2", "#E64B7D", "#DF9DA7")
colr=data.frame(colr,stringsAsFactors = F)
rownames(colr)=AMINO_ACID_CODE[1:20]

cod.aa=data.frame(codon=names(GENETIC_CODE),amino_acid=AMINO_ACID_CODE[GENETIC_CODE],combi=paste(AMINO_ACID_CODE[GENETIC_CODE],names(GENETIC_CODE),sep="->"),stringsAsFactors=FALSE)
cod.aa=cod.aa[!is.na(cod.aa$amino_acid),]
cod.aa$color=colr[cod.aa$amino_acid,1]
rownames(cod.aa)=cod.aa$combi
cod.aa$one.letter=names(AMINO_ACID_CODE)[match(cod.aa$amino_acid,AMINO_ACID_CODE)]
cod.aa$comb.2=paste0(cod.aa$one.letter,":",cod.aa$codon)
ref.mat=cbind(c(-23:-1,c('E','P','A'),3:16),paste0('X',1:40))
##########

 
load_files=function(files){
  DD=list()
  a=1
  for(i in files){
    tempi=data.frame(t(get(load(i))))
    rownames(tempi)=c('coef','std.error','t.value','pvalue')
    tempi['coef',] =tempi['coef',]/log(2)
    tempi['pvalue',which(tempi['pvalue',] == 0)]=1e-320
    DD[[a]]=(tempi)
    a=a+1
  }
  df.genes=rbindlist(DD,fill=TRUE)
  df.genes=t(df.genes)
  df.genes=data.frame(df.genes)
  
  nami=paste(rep(files,each=4),rep(c('coef','std.error','t.value','pvalue'),length(files)),sep="_")
  names(df.genes)= nami
  names(df.genes)=gsub("_coe_pval","",names(df.genes))
  names(df.genes)=gsub(".RData","",names(df.genes))
  colnames(df.genes)=gsub("25:26","P:A",colnames(df.genes))
  colnames(df.genes)=gsub("24:25","E:P",colnames(df.genes))
  colnames(df.genes)=gsub("24:26","E:A",colnames(df.genes))
  df.genes
}


df=load_files(args)
df.PA=df[grep('^......$',rownames(df)),]
df.PA=df.PA[-grep('^X',rownames(df.PA)),]
df.PA=round(df.PA,4)
df.PA$site.1 =  change_cc2aa(substring(rownames(df.PA),1,3))
df.PA$site.2 = change_cc2aa(substring(rownames(df.PA),4,6))

df.g= df[grep('gene',rownames(df)),]
rownames(df.g)=gsub('gene','',rownames(df.g))
df.g=round(df.g,4)

df.g$gene=rownames(df.g)

df.X=df[grep('^X',rownames(df)),]
df.X=round(df.X,4)

df.X$codon=rep(change_cc2aa(sapply(strsplit(rownames(df.X)[grep('X24',rownames(df.X))],"X24"),"[[",2)),40)
df.X$pos=rep(ref.mat[,1],each=61)

write.table(df.g,file="Data/Tables/summary_flux.tsv",sep="\t",quote=F,row.names=F)
write.table(df.X,file="Data/Tables/summary_single_DT.tsv",sep="\t",quote=F,row.names=F)
write.table(df.PA,file="Data/Tables/summary_pair_DT.tsv",sep="\t",quote=F,row.names=F)
