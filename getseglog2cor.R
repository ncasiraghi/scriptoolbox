#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
  message("\nERROR:\t3 arguments required")
  message("\nUSAGE:\tRscript getseglog2corr.R <outdir> <allelicImbalanceTable.txt> <ploidyTable.txt> \n")
  quit()
}

outdir <- args[1]
ait.table <- args[2]
plt.table <- args[3]

setwd(outdir)

ait = read.delim(ait.table,as.is = T,stringsAsFactors = F)
plt = read.delim(plt.table,as.is = T,stringsAsFactors = F)

segcorr <- ait[,c("sample","chr","start","end","log2","log2.corr")]

toolname <- "FACETS"
col.raw <- "#f46d43"
col.cor <- "#66bd63"

pdf('log2corr_distributions.pdf', h = 219, w = 297, paper='A4r')
par(pty="s",mfrow=c(1,2))
hist(ait$log2,main=NULL,xlim=c(-2,2),xlab = paste("log2","\n(",toolname,")"),breaks = 50,col=col.raw,border=col.raw,las=1)
mtext(paste0("pooled samples (n=",length(unique(ait$sample)),")"),line = 2)
hist(ait$log2.corr,main=NULL,xlim=c(-2,2),xlab = paste("log2.corrected","\n(",toolname,"+ CLONET)"),breaks = 50,col=col.cor,border=col.cor,las=1)
mtext(paste0("pooled samples (n=",length(unique(ait$sample)),")"),line = 2)

for(n in unique(ait$sample)){
  dd <- ait[which(ait$sample==n),]
  par(pty="s",mfrow=c(1,2))
  hist(dd$log2,main=NULL,xlim=c(-2,2),xlab = paste("log2","\n(",toolname,")"),breaks = 50,col=col.raw,border=col.raw,las=1)
  mtext(n,line = 2)
  hist(dd$log2.corr,main=NULL,xlim=c(-2,2),xlab = paste("log2.corrected","\n(",toolname,"+ CLONET)"),breaks = 50,col=col.cor,border=col.cor,las=1)
  mtext(n,line = 2)
  ploidy <- round(plt$ploidy[which(plt$sample==n)],2)
  purity <- round(1-unique(ait$adm[which(ait$sample==n)]),2)
  mtext(paste("ploidy:",ploidy,"| purity:",purity),line = 1)
}
dev.off()

write.table(x = segcorr,file = "output_log2corr.seg",quote = F,row.names = F,col.names = T,sep = "\t")


