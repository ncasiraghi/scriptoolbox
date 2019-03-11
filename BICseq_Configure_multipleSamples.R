# set working/output directory
outdir="/elaborazioni/sharedCO/UK_Het_Project/BICseq"

# Sample Info File
sif='/elaborazioni/sharedCO/UK_Het_Project/BICseq/sif_multipleSamples.csv'

## Select Human Genome Assembly
# hg19
seqNames <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
# b37
#seqNames <- c(1:22, "X", "Y")

## BICseq Parameters

# the initial genomic bin size in base pair for read counts
bin = 1000
# the penalty of the BIC
lambda = 50
# the window size for outlier identification
winSize = 200
# the probability of the read count quantile
quant = 0.95
# a positive number; a genomic position s is considered as an outlier if it has more than mult*quantile number of aligned read
mult = 1