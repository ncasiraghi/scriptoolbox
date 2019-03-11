#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2){
  message("\nERROR:\t2 arguments required")
  message("\nUSAGE:\tRscript getcns.R <folder with .cnr outputs from cnvkit.py> <outdir> \n")
  quit()
}

datadir = args[1]
outdir = args[2]

# get cns and segs

setwd(outdir)
cnrlist <- list.files(datadir,pattern = "\\.cnr$",full.names = T)
cnvkit <- "python /icgc/dkfzlsdf/analysis/B260/users/n790i/tools/cnvkit/cnvkit.py"

for(cnr in cnrlist){
  message(basename(cnr))
  cbs <- gsub(basename(cnr),pattern = "\\.cnr",replacement = ".cbs.cns")
  # cbs
  #cmd <- paste(cnvkit,"segment",cnr,"-o -p 10",cbs)
  cmd <- paste(cnvkit,"segment -t 1e-6 -p 10",cnr,"-o",cbs)
  system(cmd)
  # flasso
  if(F){
    flasso <- gsub(basename(cnr),pattern = "\\.cnr",replacement = ".flasso.cns")
    cmd <- paste(cnvkit,"segment -m flasso",cnr,"-o",flasso)
    system(cmd)
  }
}

# generate seg file (cbs method)
out.seg <- "allsamples.cbs.seg"
cmd <- paste(cnvkit,"export seg *.cbs.cns -o",out.seg)
system(cmd)

if(FALSE){
  # generate seg file (flasso method)
  out.seg <- "allsamples.flasso.seg"
  cmd <- paste(cnvkit,"export seg *.flasso.cns -o",out.seg)
  system(cmd)
}


