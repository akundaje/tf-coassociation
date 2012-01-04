rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments

# Print usage
print.usage <- function(){
  cat("Rscript batch.convert.mtrx2Rdata.R [mtrxDir] [useRelaxed]\n")
  cat("Converts all .mtrx file in directories in [mtrxDir]\n")
  cat("   [mtrxDir]: directory containing .mtrx association files\n")
  cat("   [useRelaxed]: (OPTIONAL) Set to F if you want to set relaxed peak values(<=0) to 0\n")
}

if (nargs < 1) {
  print.usage()
  q(save="no",status=1)
}

matrix.dir <- args[[1]] # mtrx data directory
use.relaxed <- T # Default work directory

if (nargs > 1) {
  use.relaxed <- as.logical(args[[2]])
}

source("assoc.matrix.utils.R")

for (each.dir in list.dirs(matrix.dir)) {
  cat("Converting directory ",each.dir,"\n")
  batch.read.assoc.file.to.Rdata(each.dir,use.relaxed=use.relaxed)
}
