#!/usr/bin/env Rscript
library(data.table)
args = commandArgs(trailingOnly=TRUE)

start_pos = args[1]
L_1 = args[2]
L_2 = args[3]
A_site_end = args[4]
out_file_tsv = args[5]
out_file_pdf = args[6]

table_pos = fread(start_pos, sep = "\t", stringsAsFactors = F)[, -203]

length_pos = as.numeric(gsub("L:","",table_pos$V1))
table_pos=table_pos[,-1]
colnames(table_pos)=as.character(c(-100:0,1:100))

tp = rowSums(table_pos) > 10
table_pos = subset(table_pos, tp)
length_pos = length_pos[tp]

sum_pos = sweep(table_pos,
                MARGIN = 1,
                rowSums(table_pos),
                FUN = '/')

ss.pos = split(1:nrow(sum_pos), length_pos)
sum_pos.l = sapply(ss.pos, function(x)
  colMeans(sum_pos[x, ]))
l=as.character(L_1:L_2)
l=l[l%in%colnames(sum_pos.l)]

pdf(out_file_pdf)
par(mfrow=c(2,2),pty='s')
for(pos in l){
  plot(1:201, sum_pos.l[,pos], xlim=c(70,130), main=pos, xaxt='n', lty='blank', cex=0); abline(v=101, col='grey')
  axis(at=1:201, labels=rownames(sum_pos.l), side=1, cex.axis=0.4)
  lines(seq(1,nrow(sum_pos.l),3), sum_pos.l[seq(1,nrow(sum_pos.l),3), pos], col='darkred', type='h')
  lines(seq(2,nrow(sum_pos.l),3), sum_pos.l[seq(2,nrow(sum_pos.l),3), pos], col='darkblue', type='h')
  lines(seq(3,nrow(sum_pos.l),3), sum_pos.l[seq(3,nrow(sum_pos.l),3), pos], col='darkgreen', type='h')
}
dev.off()

if(A_site_end == '5p'){
pos_am=25:29
}else{
pos_am=38:41
}

A_pos=t(sapply(l,function(x) list(a=names(which.max(sum_pos.l[seq(2,nrow(sum_pos.l),3)[pos_am], x])),
                                                  b=names(which.max(sum_pos.l[seq(3,nrow(sum_pos.l),3)[pos_am], x])),
                                                  c=names(which.max(sum_pos.l[seq(1,nrow(sum_pos.l),3)[pos_am], x])))))

rn=rownames(A_pos)
A_pos=apply(A_pos,2,as.numeric)
rownames(A_pos)=rn
A_pos=abs(A_pos)

write.table(A_pos, file=out_file_tsv, sep="\t", quote=F,col.names = F)
