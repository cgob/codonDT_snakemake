#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
########  
require("methods")
require("MatrixModels")
require("Matrix")
require('MASS') 
require('data.table')
require('stringi') 
#######

MakeFit=function(input_file,input_file_2, output_file,fit.param,fit.inter,modi){

###### Evalutation of input parameters for the relative positions of the single and codon pair covariates ###### 
 load(input_file)
 fit.param=eval(parse(text=fit.param))
 fit.param=paste0('X',fit.param)
 sp=unlist(strsplit(fit.inter,split=":"))
 fit.inter=paste0('X',sp)
 pairi=paste0(fit.inter,collapse="")
 
 if(input_file == input_file_2){
  modi = "simple"
 }
 
 ncount[,'count.rfp']=as.numeric(as.character(unlist(ncount[,'count.rfp']))) 
 ncount[,pairi]=paste0(ncount[,fit.inter[1]],ncount[,fit.inter[2]])

#################################################################################################################
sum.pair=tapply(ncount$count.rfp,ncount[,pairi],sum)
nb.pair=tapply(ncount$count.rfp,ncount[,pairi],length)
low.pair=names(which(sum.pair < 10 | nb.pair < 5))
ss=split(1:nrow(ncount),ncount$gene)
mean.gene =sapply(ss,function(x) mean(ncount$count.rfp[x],na.rm=T))
ncount$mean.gene=mean.gene[ncount$gene]
ncount[ncount[,pairi]%in%low.pair,'count.rfp']=ncount[ncount[,pairi]%in%low.pair,'mean.gene']
 
##### Compute overdisperion parameter by linear regression of pairs reads counts mean and variance ############# 
 ncount$pagene=paste0(ncount$gene,ncount[,pairi])
 ss=split(1:nrow(ncount),ncount$pagene)
 ss.l=sapply(ss,length)
 ss=ss[ss.l>1]
 m.raw=sapply(ss, function(x) mean(ncount$count.rfp[x]))
 v.raw=r=sapply(ss, function(x) var(ncount$count.rfp[x]))
 alpha=as.numeric(coefficients(lm(formula = v.raw - m.raw ~ 0 + I(m.raw^2))))
 theta=1/alpha
###############################################################################################################

##### Output count per gene ############
 ss=split(1:nrow(ncount),ncount$gene)
 gg=sapply(ss,function(x) sum(ncount$count.rfp[x]))
 DFF.2=data.frame(COUNT=gg)
########################################
 
######### Fit with GLM for without RNA-Seq offset ###############
# ncount=subset(ncount,ncount$gene%in%rownames(subset(DFF.2,COUNT>300)))
 ncount[grep('^AAA|AAA$',ncount[,pairi]),pairi]='AACAAC' # Use all the pairs with AAA + AAC:AAC as ref. 
if(modi=='simple'){
 formu=as.formula(paste("count.rfp ~ 0 + gene "  ,paste(" + ",fit.param,collapse="",sep=""),"+",paste0(fit.inter[1],fit.inter[2]))) # Use single and pairs input parameters as covariates
 fit.rna=glm4(formu, data=ncount,family=negative.binomial(theta),MXITER=400,offset=rep(log(sum(DFF.2$COUNT)),nrow(ncount)),doFit=T, sparse=T, verbose=F) # Fit with theta as dispersion and the library size as offset
 save(fit.rna,file=output_file) # Save the fit
 

#################################################################

###################### Prediction and correlation between data and fit for simple fit  ###############################
mm=sparse.model.matrix(data=ncount,formu)
ite=intersect(colnames(mm),names(coefficients(fit.rna)))
predi.1=mm[,ite]%*%coefficients(fit.rna)[ite] 
ncount$predi=exp(as.numeric(predi.1) + log(sum(DFF.2$COUNT)))
ncount$res.p=resid(fit.rna,type='pearson')
ncount$res.d=resid(fit.rna,type='deviance')
ncount$res.w=resid(fit.rna,type='working')
ncount$res.r=resid(fit.rna,type='response')
ncount$diff=ncount$count.rfp - ncount$predi
save(ncount,file=paste0(output_file,".pred"))
ss=split(1:nrow(ncount),ncount$gene)
ss.r=sapply(ss,function(x) cor(log(1+ncount$predi[x]),log(1+ncount$count.rfp[x])))
save(ss.r,file=paste0(output_file,".cor"))
#######################################################################################################
}
################ Fit RNA-Seq first (if available) and use the prediction as an offset for the ribosome profiling fit
 if(modi=='combined'){
   fit.rna=get(load(input_file_2))
   formu=as.formula(paste("count.rfp ~ 0 + gene ", paste(" + ",fit.param,collapse="",sep=""),"+",paste0(fit.inter[1],fit.inter[2])))
   ncount=ncount[paste0('gene',ncount$gene)%in%names(coefficients(fit.rna)),]
   mm=sparse.model.matrix(data=ncount,formu)
   ite=intersect(colnames(mm),names(coefficients(fit.rna)))
   
   predi.1=mm[,ite]%*%coefficients(fit.rna)[ite] 
   ncount$predi=as.numeric(predi.1)
   rm(fit.rna)
   gc() 
   fit.rfp.rna=glm4(formu , data=ncount,family=negative.binomial(theta),offset=ncount$predi + rep(log(sum(DFF.2$COUNT)),nrow(ncount)), doFit=T, sparse=T, verbose=F)
   save(fit.rfp.rna,file=output_file)
 } 
##################################################################################################################### 
}
MakeFit(args[1],args[2],args[3],args[4],args[5], args[6])
