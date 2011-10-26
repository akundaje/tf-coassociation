rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments
rm.target=F  # Set to true if you want to remove target variable from predictors

# Print usage
print.usage <- function(){
  cat("Rscript crossval.model.randneg.R [ruleFitResultRdataFile] [replaceFlag]\n")
  cat("Compute cross-validation results for a dataset\n")
  cat("   [ruleFitResultRdataFile]: Rdata file containing rulefit.results directory\n")
  cat("   [replaceFlag]: (OPTIONAL) T/F (If set to F then if rulefit.results has a variable called cv the run is aborted) Default: T\n")
}

if (nargs < 1) {
  print.usage()
  q(save="no",status=1)
}

rulefit.results.file <- args[[1]] # Directory containing rulefit.results Rdata files from several random negative set runs
if (! file.exists(rulefit.results.file)) {
  cat("Input File ", rulefit.results.file,"does not exist\n")
  q(save="no",status=1)
}

replace.flag <- T
if (nargs > 1) {
  replace.flag <- as.logical(args[[2]])
}

work.dir <- tempfile() # Default work directory
if (nargs > 3) {
  work.dir <- args[[4]]
}
if (! file.exists(work.dir)) {
  dir.create(work.dir)
}

source("assoc.matrix.utils.R")
platform <- "linux"
rfhome <- initialize.rulefit( work.dir=work.dir, rf.package.path=Sys.getenv("RULEFITBASE") )

load(rulefit.results.file) # Loads rulefit.results

# Check if CV results already exist
if ( ("cv" %in% names(rulefit.results)) && (!replace.flag) ) {
  q(save="no",status=0)
}

rulefit.results <- run.cv.rulefit(rulefit.results)

save(list="rulefit.results",file=rulefit.results.file)