#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
dat = read.table(args[1], sep="\t", stringsAsFactors = F )

pdf(args[2])

 par(pty='s',mfrow=c(1,3))
 plot(dat$V1,dat$V2,type='l',xlab='Size',ylab="Nb. Reads")
 plot(dat$V1[20:40],dat$V2[20:40],type='l',xlab='Size',ylab="Nb. Reads")
 plot(dat$V1[50:70],dat$V2[50:70],type='l',xlab='Size',ylab="Nb. Reads")

dev.off()
