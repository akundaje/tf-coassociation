rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments
# rm.target=F  # Set to true if you want to remove target variable from predictors
# trim.target=T # Set to false if you want to use all rows in the co-binding matrix (if T, then rows with target.TF binding values < 0 are removed)
# append.null=F # Set to T if you want to use randomized features as extra null model features

# Print usage
print.usage <- function(){
  cat("Rscript run.model.randneg.R [mtrxRdata] [outputFile] [workingDir]\n")
  cat("Samples a random negative set, learns a rulefit model, computes variable importance, interaction strength and all pairwise interactions\n")
  cat("   [mtrxRdata]: .Rdata file containing association matrix\n")
  cat("   [outputFile]: Path to output File\n")
  cat("   [replaceFlag]: (OPTIONAL) If set to a 0, then if output file exits, the run will quit and not replace it, DEFAULT:0\n")
  cat("   [workingDir]: (OPTIONAL) working directory, DEFAULT: tempfile()\n")
  cat("   [rm.target]: (OPTIONAL) Set to T if you want to remove target variable, DEFAULT: F\n")
  cat("   [trim.target]: (OPTIONAL) Set to T if you want to remove rows for which target TF has values < 0, DEFAULT: F()\n")
  cat("   [append.null]: (OPTIONAL) Set to T if you want to append matrix with extra null randomized features, DEFAULT: F()\n")
}

if (nargs < 2) {
  print.usage()
  q(save="no",status=1)
}

mtrx.rdata.file <- args[[1]] # Association matrix Rdata file (Contains list called assoc.data)
if (! file.exists(mtrx.rdata.file)) {
  cat("Input matrix Rdata file ",mtrx.rdata.file,"does not exist\n")
  q(save="no",status=1)
}

output.file <- args[[2]] # Output file

replace.flag <- 0 # do not replace existing file
if (nargs > 2) {
  replace.flag <- as.logical(args[[3]])
}

work.dir <- tempfile() # Default work directory
if (nargs > 3) {
  work.dir <- args[[4]]
}
if (! file.exists(work.dir)) {
  dir.create(work.dir)
}

rm.target <- F
if (nargs > 4) { # rm.target
  rm.target <- as.logical(args[[5]])
}

trim.target <- F
if (nargs > 5) { # trim.target
  trim.target <- as.logical(args[[6]])
}

append.null <- F
if (nargs > 6) {
  append.null <- as.logical(args[[7]])
}

source("assoc.matrix.utils.R")
platform <- "linux"
rfhome <- initialize.rulefit( work.dir=work.dir, rf.package.path=Sys.getenv("RULEFITBASE") )

output.dir <- get.file.parts(output.file)$path
if (! file.exists(output.dir)) {
  dir.create(output.dir, recursive=T)
}

#output.file <- file.path(output.dir, sprintf("%s.randneg.rfresults.Rdata", mtrx.rdata.file.name))

if (! replace.flag & file.exists(output.file)) {
  cat("Output file exists: ",output.file,". Use replace.flag=1 to replace the file\n")
  q(save="no",status=0)
}

# ================
# Run Rulefit
# ================
# This loads a list called assoc.data
# $assoc.matrix(DATAFRAME): dataFrame that is the association matrix (rows: TF sites, cols: partner TFs)
#    rownames are peakIDs. Rows are sorted by the numeric part of PeakIDs. First column is target
# $target.name(STRING): name of target
# $assoc.mtrx.file <- path of association .mtrx file
# $assoc.R.file <- path of association .Rdata file
# $expr.val <- expression values (if not valid it is assigned NA)
load(mtrx.rdata.file) 

# Compute optimal rulefit model and variable importance
#   rulefit.results$rfmod
#   rulefit.results$dataset
#   rulefit.results$vi
#   rulefit.results$int.strength
#   rulefit.results$pair.interactions
cat("Computing model and variable importance ...\n")
#rulefit.results <- get.average.randneg.model(assoc.data,rm.target=rm.target)
rulefit.results <- get.var.imp( sample.randneg.rulefit.model(assoc.data=assoc.data,
                                                             rm.target=rm.target,
                                                             trim.target=trim.target,
                                                             append.null=append.null) )
save(list="rulefit.results",file=output.file)

# Compute interaction strength
cat("Computing global interaction strength ...\n")
rulefit.results <- get.int.strength( rulefit.results )
save(list="rulefit.results",file=output.file)

# Compute pairwise interactions
cat("Computing all pairwise interactions\n")
rulefit.results <- get.all.partner.pair.interactions(rulefit.results, use.import=T)
save(list="rulefit.results",file=output.file)

# Compute cv
cat("Computing cross-validation accuracy\n")
rulefit.results <- run.cv.rulefit(rulefit.results)
save(list="rulefit.results",file=output.file)
