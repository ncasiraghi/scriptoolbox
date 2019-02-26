#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
  message("\nERROR:\t3 arguments required")
  message("\nUSAGE:\tRscript getinformativesnps.R <CLONET.sif.file> <aseq.folder.normal> <aseq.folder.tumor> \n")
  quit()
}

sif.file <- args[1]
aseq.folder.normal <- args[2]
aseq.folder.tumor <- args[3]

aseq.normals <- list.files(aseq.folder.normal,pattern = "\\.PILEUP\\.ASEQ$",full.names = T)
aseq.tumors <- list.files(aseq.folder.tumor,pattern = "\\.PILEUP\\.ASEQ$",full.names = T)
sif <- read.delim(sif.file,as.is = T,stringsAsFactors = F)

for(n in unique(sif$Normal.Bam.Name)){
  setwd(aseq.folder.normal)
  message(basename(n))
  aseq.file <- grep(aseq.normals,pattern = gsub(basename(n),pattern = ".bam",replacement = ".PILEUP.ASEQ"),value = T)
  snps.normal <- gsub(basename(aseq.file),pattern = ".PILEUP.ASEQ",replacement = ".snps")
  cmd <- paste("head -n1",aseq.file,">",paste0(snps.normal,".tmp"))
  system(cmd)
  cmd <- paste("awk -F'\t' '{ if($16 >= 0.2 && $16 <= 0.8) { print } }'",aseq.file,">>",paste0(snps.normal,".tmp"))
  system(cmd)
  cmd <- paste("cut -f 1-3,5-10,16,17",paste0(snps.normal,".tmp"),">",snps.normal)
  system(cmd)
  file.remove(paste0(snps.normal,".tmp"))
  rs <- unique(fread(input = snps.normal,sep = "\t",header = T,stringsAsFactors = F,verbose = F,select = c(3),data.table = F))
  outfile <- file.path(aseq.folder.normal,"snplist.txt")
  cat(rs$dbsnp,sep = "\n",file = outfile,append = F)
  # extract snps from tumors
  setwd(aseq.folder.tumor)
  tumors.list <- unique(sif$Tumor.Bam.Name[which(sif$Normal.Bam.Name==n)])
  for(m in tumors.list){
    message(basename(m))
    aseq.file <- grep(aseq.tumors,pattern = gsub(basename(m),pattern = ".bam",replacement = ".PILEUP.ASEQ"),value = T)
    snps.tumor <- gsub(basename(aseq.file),pattern = ".PILEUP.ASEQ",replacement = ".snps")
    cmd <- paste("head -n1",aseq.file,">",paste0(snps.tumor,".tmp"))
    system(cmd)
    cmd <- paste("grep -Fwf",outfile,aseq.file,">>",paste0(snps.tumor,".tmp"))
    system(cmd)
    cmd <- paste("cut -f 1-3,5-10,16,17",paste0(snps.tumor,".tmp"),">",snps.tumor)
    system(cmd)
    file.remove(paste0(snps.tumor,".tmp"))
  }
  file.remove(file.path(aseq.folder.normal,"snplist.txt"))
}


