rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments
rm.target=T  # Set to true if you want to remove target variable from predictors
use.null=T # Set to T if you want to compute null models for feature interactions

# Print usage
print.usage <- function(){
  cat("Rscript run.model.posneg.R [posRdata] [negRdata] [outputFile] [workingDir]\n")
  cat("Combines a positive and negative set, learns a rulefit model\n")
  cat("   [posRdata]: .Rdata file containing association matrix\n")
  cat("   [negRdata]: .Rdata file containing association matrix\n")
  cat("   [outputFile]: Path to output File\n")
  cat("   [replaceFlag]: (OPTIONAL) If set to a F, then if output file exits, the run will quit and not replace it, DEFAULT:F\n")
  cat("   [workingDir]: (OPTIONAL) working directory, DEFAULT: tempfile()\n")
}

if (nargs < 3) {
  print.usage()
  q(save="no",status=1)
}

pos.Rdata.file <- args[[1]] # Association matrix Rdata file (Contains list called assoc.data)
if (! file.exists(pos.Rdata.file)) {
  cat("Positive Input matrix Rdata file ",pos.Rdata.file,"does not exist\n")
  q(save="no",status=1)
}

neg.Rdata.file <- args[[2]] # Association matrix Rdata file (Contains list called assoc.data)
if (! file.exists(neg.Rdata.file)) {
  cat("Negative Input matrix Rdata file ",neg.Rdata.file,"does not exist\n")
  q(save="no",status=1)
}

output.file <- args[[3]] # Output File name

replace.flag <- F # do not replace existing file
if (nargs > 3) {
  replace.flag <- as.logical(args[[4]])
}

work.dir <- tempfile() # Default work directory
if (nargs > 4) {
  work.dir <- args[[5]]
}
if (! file.exists(work.dir)) {
  dir.create(work.dir)
}

source("assoc.matrix.utils.R")
platform <- "linux"
rfhome <- initialize.rulefit( work.dir=work.dir, rf.package.path=Sys.getenv("RULEFITBASE") )

output.dir <- get.file.parts(output.file)$path
if (! file.exists(output.dir)) {
  dir.create(output.dir, recursive=T)
}

if (! replace.flag && file.exists(output.file)) {
  cat("Output file exists: ",output.file,". Use replace.flag=T to replace the file\n")
  q(save="no",status=0)
}

# ================
# Run Rulefit
# ================
# $assoc.matrix(DATAFRAME): dataFrame that is the association matrix (rows: TF sites, cols: partner TFs)
#    rownames are peakIDs. Rows are sorted by the numeric part of PeakIDs. First column is target
# $target.name(STRING): name of target
# $assoc.mtrx.file <- path of association .mtrx file
# $assoc.R.file <- path of association .Rdata file
# $expr.val <- expression values (if not valid it is assigned NA)

# Compute optimal rulefit model and variable importance
#   rulefit.results$rfmod
#   rulefit.results$dataset
#   rulefit.results$vi
#   rulefit.results$int.strength
#   rulefit.results$pair.interactions

cat("Computing model and variable importance ...\n")
rulefit.results <- learn.posneg.rulefit.model(pos.Rdata.file, neg.Rdata.file, rm.target=rm.target)
rulefit.results <- get.var.imp( rulefit.results,class=1 )
save(list="rulefit.results",file=output.file)

# # Compute interaction strength
cat("Computing global interaction strength ...\n")
rulefit.results <- get.int.strength( rulefit.results, use.null=use.null )
save(list="rulefit.results",file=output.file)

# Compute pairwise interactions
cat("Computing all pairwise interactions\n")
rulefit.results <- get.all.partner.pair.interactions(rulefit.results, use.import=T, use.null=use.null)
save(list="rulefit.results",file=output.file)

# Compute cross-validation error
rulefit.results <- run.cv.rulefit(rulefit.results)
save(list="rulefit.results",file=output.file)