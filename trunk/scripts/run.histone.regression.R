rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments

# Print usage
print.usage <- function(){
  cat("Rscript run.histone.regression.R [histonePredictorMatrixFile] [outputDir]\n")
  cat("Creates a regression model linking histone marks to TF binding strength\n")
  cat("   [histonePredictorMatrixFile]: File containing histone predictors corresponding to TF peaks\n")
  cat("   [outputDir]: Directory where you want output files to go\n")
}

# Check arguments
if (nargs < 2) {
  print.usage()
  q(save="no",status=1)
}

work.dir <- tempfile() # Default work directory
if (! file.exists(work.dir)) {
  dir.create(work.dir)
}
# setwd('/media/fusion10/work/encode/learning/combinatorics/src/scripts/')
source('assoc.matrix.utils.R')
platform="linux"
# rfhome <- initialize.rulefit(rf.package.path='/media/fusion10/work/encode/learning/combinatorics/src/code')
rfhome <- initialize.rulefit( work.dir=work.dir, rf.package.path=Sys.getenv("RULEFITBASE") )

assoc.data <- args[[1]]
if (! file.exists(assoc.data)) {
  cat("Association input file ", assoc.data," does not exist\n")
  q(save="no",status=1)
}

# Create output directory if it doesn't exist
output.dir <- args[[2]]
if (! file.exists(output.dir)) { dir.create(output.dir, recursive=T) }

# Filter predictors to remove redundant predictors

# Learn association model
rulefit.results <- learn.histone.to.tf.regression.model(assoc.data,inverted=T)

# Create subdirectory named after target name in output.dir
output.subdir <- sprintf("%s/%s",output.dir,rulefit.results$dataset$target.name)
if (! file.exists(output.subdir)) { dir.create(output.subdir, recursive=T) }

# Get variable importance and create plot
rulefit.results <- get.var.imp(rulefit.results,class=0)
save(list="rulefit.results",
     file=sprintf("%s/%s/%s.Rdata",output.dir,rulefit.results$dataset$target.name,rulefit.results$dataset$target.name))              

# Get interaction strength and plot it
rulefit.results <- get.int.strength(rulefit.results)
save(list="rulefit.results",
     file=sprintf("%s/%s/%s.Rdata",output.dir,rulefit.results$dataset$target.name,rulefit.results$dataset$target.name))              

# Get all pairwise partner interactions and plot pairwise interaction matrix
rulefit.results <- get.all.partner.pair.interactions(rulefit.results,use.import=T)

# Compute cross-validation results
rulefit.results <- run.cv.rulefit(rulefit.results)
cat("Cross-validation R-square=",rulefit.results$cv$rsquare,"\n")

save(list="rulefit.results",
     file=sprintf("%s/%s/%s.Rdata",output.dir,rulefit.results$dataset$target.name,rulefit.results$dataset$target.name))              

plot.importance(rulefit.results,
                output.dir=output.dir,
                output.filename=sprintf("predictor.importance.%s.png",rulefit.results$dataset$target.name) )

plot.heatmap(data=rulefit.results$pair.interactions,
             replace.diag=T,
             replace.na=F,
             filt.thresh=1e-7,
             clust.method="ward",
             to.file=sprintf("%s/%s/cond.pairwise.int.matrix.%s.png",output.dir,rulefit.results$dataset$target.name,rulefit.results$dataset$target.name),
             row.title="Histone Marks",
             col.title="Histone Marks",
             title.name=rulefit.results$dataset$target.name)
# Also plot all pairwise interaction barplot
plot.pairwise(rulefit.results, output.dir=output.dir)

# Plot single variable partial dependence plots
