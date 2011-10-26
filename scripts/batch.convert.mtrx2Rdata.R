rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments

# Print usage
print.usage <- function(){
  cat("Rscript batch.convert.mtrx2Rdata.R [mtrxDir] [workingDir]\n")
  cat("Converts all .mtrx file in directories in [mtrxDir]")
  cat("   [mtrxDir]: directory containing .mtrx association files\n")
  cat("   [workingDir]: (OPTIONAL) working directory\n")
}

if (nargs < 1) {
  print.usage()
  q(save="no",status=1)
}

matrix.dir <- args[[1]] # mtrx data directory
work.dir <- tempfile() # Default work directory

if (nargs > 1) {
  work.dir <- args[[2]]
}

if (! file.exists(work.dir)) {
  dir.create(work.dir)
}

source("assoc.matrix.utils.R")
platform <- "linux"
rfhome <- initialize.rulefit( work.dir=work.dir, rf.package.path=Sys.getenv("RULEFITBASE") )

for (each.dir in list.dirs(matrix.dir)) {
  cat("Converting directory ",each.dir,"\n")
  batch.read.assoc.file.to.Rdata(each.dir)
}
