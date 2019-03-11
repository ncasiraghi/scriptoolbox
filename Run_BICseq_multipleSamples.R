#!/usr/bin/env Rscript
# get Rpackage from:
# http://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq/BICseq_1.1.11.tar.gz

library(BICseq)

# Source BICseq configuration file
message('[',Sys.time(),']\tLoading configuration file ...')
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\n\tError!\n\tUsage: Run_BICseq_multipleSamples.R BICseq_Configure_multipleSamples.R\n")
  quit()
}
configuration_file <- args[1]
source(configuration_file)
# Load sample info file
SampleInfoFile = read.table(sif,header=T,sep="\t",stringsAsFactors = F)
for(id in 1:nrow(SampleInfoFile)){
  message('[',Sys.time(),']\tProcessing sample case:\t',basename(SampleInfoFile[id,1]))
  CASE = SampleInfoFile[id,1]
  CTRL = SampleInfoFile[id,2]
  ## BICseq segmentation
  message('[',Sys.time(),']\tCompute bicseq object ...')
  bicseq <- BICseq(sample = CASE, reference = CTRL, seqNames = seqNames)
  message('[',Sys.time(),']\tCompute segs ...')
  segs <- getBICseg(object = bicseq, bin = bin, lambda = lambda, winSize = winSize, quant = quant, mult = mult)
  message('[',Sys.time(),']\tCompute log2 copy ratios ...')
  bins <- BICseq:::getRatios(bin(segs), what = "bin")
  seg.summary <- BICseq:::getSummary(object = segs)
  message('[',Sys.time(),']\tWriting seg file ...')
  seg.summary <- as.data.frame(append(seg.summary, list(sample = gsub(x = basename(CASE),pattern = ".bam",replacement = "")), after = 0))
  seg.names <- c("sample","chrom","start","end","probes","logR")
  seg.summary <- seg.summary[, !names(seg.summary) %in% c("referenceReads","log10.pvalue")]
  output.name.seg <- paste(basename(CASE),bin,lambda,winSize,"seg",sep='.')
  write.table(x = seg.summary, file = output.name.seg, quote = F,row.names = F, col.names = seg.names, sep ="\t")
  message('[',Sys.time(),']\tSample completed.')
}
message('[',Sys.time(),']\tDone.')