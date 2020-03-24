require(plyr)
require(Biostrings)
require(RColorBrewer)
require(ggplot2)
require(gplots)
require(data.table)
require(gtools)
require(grid)
require(gridExtra)

#########################

args = commandArgs(trailingOnly=TRUE)
files=args[1]
pdf.out=args[2]

##########################

extract_x=function(x,d,i){
  sapply(strsplit(as.character(x),d),"[[",i)
  
}

##########################

change_cc2aa=function(cod){
  sub = as.character(AMINO_ACID_CODE[as.character(GENETIC_CODE[cod])])
  sub[is.na(sub)]='Stop'
  cod = paste(sub,cod,sep="->")
}

##########################

plot_heatmap.single=function(dff.c,margin,titi,colv){
  dff.c=dff.c[order(extract_x(rownames(dff.c),"->",1),extract_x(rownames(dff.c),"->",2)),]
  dff.c.2=dff.c
  rownames(dff.c.2)=cod.aa[rownames(dff.c),'comb.2']
  
  heatmap.2(t(as.matrix(dff.c.2)) ,
            Rowv=FALSE,
            Colv=colv,
            dendrogram="col",
            cexRow = 0.6,
            cexCol=.4,
            trace="none",
            col=colorRampPalette(brewer.pal(11,"RdBu"))(100),
            scale='none',key=TRUE,
            ColSideColors=as(cod.aa[rownames(dff.c),'color'],"character"),
            margins=margin,
            density.info='none',
            main=titi,
            na.color = 'black')
}

#############################################

plot_heatmap.pair=function(dff.c,titi,site.1,site.2,clust.1,clust.2){
  
  dff.c.2=dff.c
  rownames(dff.c.2)=cod.aa[rownames(dff.c),'comb.2']
  colnames(dff.c.2)=cod.aa[colnames(dff.c),'comb.2']
  
  heatmap.2(dff.c.2,col=colorRampPalette(brewer.pal(11,"RdBu"))(100),
            Rowv=clust.1,
            Colv=clust.2,
            dendrogram="both",
            cexRow = 0.6,
            cexCol=0.6,
            trace="none",
            scale='none',
            key=TRUE,
            sepcolor='white',
            main=titi,
            ColSideColors=as(cod.aa[rownames(dff.c),'color'],"character"),
            RowSideColors=as(cod.aa[colnames(dff.c),'color'],"character"),
            xlab=site.1,
            ylab= site.2,
            density.info='none', 
            notecol="black",
            na.color = 'black')
} 

#############################################

generate_data=function(fit){
  
  df.genes=fit
  df.genes.PA=df.genes[grep('^......$',rownames(df.genes)),]
  df.genes.PA=df.genes.PA[-grep('^X',rownames(df.genes.PA)),]
  df.genes.g= df.genes[grep('gene',rownames(df.genes)),]
  df.genes.X=df.genes[grep('^X',rownames(df.genes)),]
  
  PA.p=df.genes.PA[,grep('pval',names(df.genes.PA)),drop=F]
  PA.c=df.genes.PA[,grep('coef',names(df.genes.PA)),drop=F]/log(2)
  SI.p=df.genes.X[,grep('pval',names(df.genes.X)),drop=F]
  SI.c=df.genes.X[,grep('coef',names(df.genes.X)),drop=F]/log(2)
  G.p= df.genes.g[,grep('pval',names(df.genes.g)),drop=F]
  G.c= df.genes.g[,grep('coef',names(df.genes.g)),drop=F]/log(2)
  
  PA.p[PA.p==0]=1e-320
  SI.p[SI.p==0]=1e-320
  G.p[G.p==0]=1e-320
  
  return(list(SI.c=SI.c,SI.p=SI.p,PA.c=PA.c,PA.p=PA.p,G.c=G.c,G.p=G.p))
}


################################################


load_files=function(files,nb){
  
  tempi=data.frame((get(load(files))))
  tempi=tempi[,c(1,4)]
  colnames(tempi)=c('coef','pvalue')
  tempi
}


plot_heatmaps=function(fit,i,site.1,site.2){
  
  
  SI.c=fit$SI.c;SI.p=fit$SI.p;PA.c=fit$PA.c;PA.p=fit$PA.p;G.c=fit$G.c;G.p=fit$G.p
  
  DF.SI.p=matrix(SI.p[,1],61,nrow(SI.p)/61)
  DF.SI.c=matrix(SI.c[,1],61,nrow(SI.c)/61)
  
  rownames(DF.SI.c)=change_cc2aa(sapply(strsplit(rownames(SI.c)[grep('X24',rownames(SI.c))],"X24"),"[[",2))
  colnames(DF.SI.c)=c(seq(-23,-1,1),c('E','P','A'),3:16)
  
  plot_heatmap.single(DF.SI.c,c(5,7),i,T)
  plot_heatmap.single(DF.SI.c,c(5,7),i,F)
  plot_heatmap.single(DF.SI.c[,grep("E|P|A",colnames(DF.SI.c))],c(32,5),i,T)
  plot_heatmap.single(DF.SI.c[,grep("E|P|A",colnames(DF.SI.c))],c(32,5),i,F)
  
  DD.s.m=melt(DF.SI.c)
  DD.s.m$amino =cod.aa[as.character(DD.s.m[,1]),'amino_acid']
  DD.s.m$codon =cod.aa[as.character(DD.s.m[,1]),'codon']
  DD.s.m$fit=extract_x(DD.s.m$Var2,"_",1)
  names(DD.s.m)[1:2]=c('mixed','site')
  
  g1=ggplot(data=subset(DD.s.m,site%in%c("E","P","A") & value!=0),aes(x=site,y=value,fill=amino,group=interaction(site,amino))) + 
    geom_point(position = position_dodge(width=0.8),shape=21,size=3.5,col='white')+
    ylab("dwell times (log2, centered)") + 
    theme_bw() + 
    xlab("") +
    scale_fill_manual(values=as.character(cod.aa$color[match(levels(as.factor(DD.s.m$amino)),cod.aa$amino_acid)])) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(size=17),axis.text.y = element_text(size=17),axis.title = element_text(size=15),legend.text =  element_text(size=15)) 
  
  print(g1)
  
  P.site=  substring(rownames(PA.c),1,3)
  A.site=  substring(rownames(PA.c),4,6)
  
  DF.PA.c=data.frame(A=A.site,P=P.site,z=(PA.c[,1]))
  
  DF.PA.c=reshape2::acast(DF.PA.c,P~A,value.var="z")
  colnames(DF.PA.c)=change_cc2aa(colnames(DF.PA.c))
  rownames(DF.PA.c)=change_cc2aa(rownames(DF.PA.c))
  
  DF.PA.p=data.frame(A=A.site,P=P.site,z=PA.p[,1])
  DF.PA.p=reshape2::acast(DF.PA.p,P~A,value.var="z")
  DF.SI.c[DF.SI.p>=0.05]=0
  DF.PA.c[DF.PA.p>=0.05]=0

  if(!all(DF.PA.c==0)){  
  plot_heatmap.pair(DF.PA.c,i,"site 1","site 2",TRUE,TRUE)
  }
}
################################################

colr=c("#87597C", "#DF9158", "#D6ACE3", "#8CE55E", "#65E7D1", "#D042E4", "#E4D94E", "#DFD5DF", "#7C84DC", "#DB6ECB", "#6E9574", "#8CD5E0", "#78A9D4", "#7AE296", "#DED3B2", "#D3D57F", "#C0EAC9", "#7753E2", "#E64B7D", "#DF9DA7")
colr=data.frame(colr,stringsAsFactors = F)
rownames(colr)=AMINO_ACID_CODE[1:20]

cod.aa=data.frame(codon=names(GENETIC_CODE),amino_acid=AMINO_ACID_CODE[GENETIC_CODE],combi=paste(AMINO_ACID_CODE[GENETIC_CODE],names(GENETIC_CODE),sep="->"),stringsAsFactors=FALSE)
cod.aa=cod.aa[!is.na(cod.aa$amino_acid),]
cod.aa$color=colr[cod.aa$amino_acid,1]
rownames(cod.aa)=cod.aa$combi
cod.aa$one.letter=names(AMINO_ACID_CODE)[match(cod.aa$amino_acid,AMINO_ACID_CODE)]
cod.aa$comb.2=paste0(cod.aa$one.letter,":",cod.aa$codon)

##############################################

fit.1=load_files(files)
DAT.gen=generate_data(fit.1)

pdf(pdf.out)
plot_heatmaps(DAT.gen,'test','1','2')
dev.off()
  

