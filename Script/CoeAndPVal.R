#!/usr/bin/env Rscript

########  
require("methods")
require("MatrixModels")
require("Matrix")
require('MASS') 
require('data.table')
require('stringr') 
#######
args = commandArgs(trailingOnly=TRUE)

CoeAndPval =function(input_path,output_path){ 
shifter <- function(x, n = 0) {
    if (n == 0) x else c(tail(x, -n), head(x, n))
}

fit =  get(load(file=input_path))
rm(fit.rna)
rm(fit.rfp.rna)
gc()
coe=coef(fit)
coe.g=coe[grep('gene',names(coe))]
coe.s=coe[grep('^X',names(coe))]
coe.s=coe.s[-grep('X.+X',names(coe.s))]
#### ADD MISSING PARAMETERS =0 AND REORDER
paramo=names(fit@pred@X@contrasts)
fit.param=paramo[2:(length(paramo)-1)]
fit.inter=unlist(str_extract_all(paramo[length(paramo)], "X.."))
coe.s[paste0(fit.param,"AAA")]=0
coe.s=coe.s[order(as.numeric(unlist(str_extract_all(names(coe.s), "[0-9]+"))),unlist(str_extract_all(names(coe.s), "...$")))]
coe.p=coe[grep('X.+X',names(coe))]
names(coe.p)=gsub("^X\\d+X\\d+","",names(coe.p))
codon=unique(gsub("^X\\d+","",names(coe.s)))

pairi=expand.grid(codon,codon)
PA=data.frame(rep(0,nrow(pairi)),row.names=paste0(pairi[,1],pairi[,2]))
PA[names(coe.p),1]=coe.p
coe.p=t(PA[,1])
names(coe.p)=rownames(PA)
coe.p=coe.p[order(names(coe.p))]

####### dimension single, pairs, flux parameters 

n.g=length(coe.g) 
n.s=length(coe.s)
n.p=length(coe.p)
n.c.2=length(fit.param)
N=n.g + n.s + n.p
n.c=length(codon)
########## Substract colMeans (A site) pairs and add to single A site
#FLUX PART MATRIX
X.3.flux=cbind(diag(n.g),Matrix(0,n.g,n.s+n.p,sparse=T))
 
#SINGLE MATRIX
y=rep(c(1,rep(0,n.c-1)),n.c)
X.s=1/n.c*rbind(t(sapply(0:(n.c-1),function(x) shifter(y,-x))))
X.s.s=Matrix(0,nrow=n.s,ncol=n.p,sparse=TRUE)
X.s.s[grep(paste0("^",fit.inter[2],"...$"),names(coe.s)),]=X.s
X.3.s=cbind(Matrix(0,nrow=n.s,ncol=n.g,sparse=T),cbind(diag(n.s),X.s.s))

#PAIRS MATRIX
y=c(c(n.c-1,rep(0,n.c-1)),rep(c(-1,rep(0,n.c-1)),n.c-1))
X.3=1/n.c*rbind(y,t(sapply(1:(n.p-1),function(x) shifter(y,-x))))
X.3.p=cbind(Matrix(0,nrow=n.p,ncol=n.g+n.s,sparse=T),X.3)
X.3=rbind(X.3.flux,X.3.s,X.3.p)

##################Substract rowMeans (P site) pairs and add to single P site  
#FLUX PART MATRIX
 
X.2.flux=cbind(diag(n.g),Matrix(0,n.g,n.s+n.p,sparse=T))

#SINGLE MATRIX
y= rep(1,n.c)
X.s=1/n.c*t(kronecker(diag(n.c), y))
X.s.s=Matrix(0,nrow=n.s,ncol=n.p,sparse=T)
X.s.s[grep(paste0("^",fit.inter[1],"...$"),names(coe.s)),]=X.s
X.2.s=cbind(Matrix(0,nrow=n.s,ncol=n.g,sparse=T),diag(n.s),X.s.s)

#PAIRS MATRIX
y=c(n.c-1,rep(-1,n.c-1))
y.1=t(sapply(0:(n.c-1),function(x) shifter(y,-x)))
X.2=1/n.c*kronecker(diag(n.c), y.1)
X.2.p= cbind(Matrix(0,nrow=n.p,ncol=n.g+n.s,sparse=T),X.2)
X.2=rbind(X.2.flux,X.2.s,X.2.p)
#
############### Substract rowMeans individual site ###3
X.1=Matrix(diag(N),sparse=T)
y=c(n.c-1,rep(-1,n.c-1))
y.1=t(sapply(0:(n.c-1),function(x) shifter(y,-x)))
X.1.1=1/n.c*kronecker(diag(n.c.2), y.1)
X.1[(n.g+1):(n.g+n.s),(n.g+1):(n.g+n.s)]=X.1.1

######Verification ##########
coe.f=c(coe.g,coe.s,coe.p)

coe.mat=X.1%*%X.2%*%X.3%*%coe.f

# Direct normalization
PA.m=matrix(coe.p,n.c,n.c,byrow=T,dimnames=list(codon,codon))
PA.m.2=sweep(PA.m,2,colMeans(PA.m),FUN='-')
PA.m.3=sweep(PA.m.2,1,rowMeans(PA.m.2),FUN='-')

s.m=matrix(coe.s,nrow=n.c.2,ncol=n.c,byrow=T)
s.m[grep(paste0("^",fit.inter[2],"$"),fit.param),]=s.m[grep(paste0("^",fit.inter[2],"$"),fit.param),]+colMeans(PA.m)
s.m[grep(paste0("^",fit.inter[1],"$"),fit.param),]=s.m[grep(paste0("^",fit.inter[1],"$"),fit.param),]+rowMeans(PA.m.2)
mmm=rowMeans(s.m)
s.m=sweep(s.m,1,mmm,FUN="-")
coe.g=coe.g+sum(mmm)
coe.n=c(coe.g,as.numeric(t(s.m)),as.numeric(t(PA.m.3)))

# COMPARISION
#all.equal(as.numeric(coe.mat),as.numeric(coe.n))

######### Add missing coefficient to the covariance matrix and order it correctly##########
cov=Matrix::chol2inv(fit@pred@fac)
names(coe)=gsub("^X\\d+X\\d+","",names(coe))
colnames(cov) = names(coe)
rownames(cov) = names(coe)
na.coef=names(coe.f[which(is.na(match(names(coe.f),colnames(cov)))==T)])
cov.1=Matrix(0,ncol=length(na.coef),nrow=nrow(cov)+length(na.coef),sparse=T)
rownames(cov.1)=c(rownames(cov),na.coef)
colnames(cov.1)=na.coef
cov.2=Matrix(0,nrow=length(na.coef),ncol=ncol(cov),sparse=T)

rownames(cov.2)=na.coef
colnames(cov.2)=colnames(cov)
cov=rbind(cov,cov.2)
cov=cbind(cov,cov.1)
ii=match(names(coe.f),colnames(cov))
cov=cov[,ii] 
cov=cov[ii,]
 
cov.norm=X.3%*%cov%*%t(X.3)
cov.norm=X.2%*%cov.norm%*%t(X.2)
cov.norm=X.1%*%cov.norm%*%t(X.1)
#
#########Comptue the p-value##########
k=length(fit@pred@X@Dimnames[[2]])
n=length(fit@pred@X@Dimnames[[1]])
df.r=n - k
dispersion = sum(residuals(fit,type='pearson')^2)/df.r
s.err=sqrt(dispersion*diag(cov.norm))
tvalue <- coe.n/s.err
pvalue <- 2 * pt(-abs(tvalue), df.r)
coef.table <- cbind(coe.n, s.err, tvalue, pvalue)
dn <- c("Estimate", "Std. Error")
dimnames(coef.table) <- list(names(coe.n), c(dn, "t value", "Pr(>|t|)"))
rownames(coef.table)=names(coe.f)
save(coef.table, file=output_path)
}

CoeAndPval(args[1],args[2])
