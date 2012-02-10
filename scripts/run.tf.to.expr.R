rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments
source('assoc.matrix.utils.R')
two.stage <- T
low.expr.thresh <- 0
rm.low.expr <- F
min.feature.support <- 3
tree.size <- 6
test.reps <- 3
model.type <- "both"

# Print usage
print.usage <- function(){
  cat("Rscript run.tf.to.expr.R [tfAssocRdataFile] [exprRdataFile] [output.dir] [rna.type] [comp.type] [tail.type]\n")
  cat("Fits TF to expression regression model\n")
  cat("   [tfAssocRdataFile]: input Rdata file containing gene centric TF association matrix\n")
  cat("   [exprRdataFile]: input Rdata file containing different types of expression data for all genes\n")
  cat("   [output.Dir]: Output directory where you want runs and results to be saved\n")
  cat("   [rna.type]: (OPTIONAL) Type of RNA eg. rnaseq, cage, ditag. If not specified all types will be used one by one\n")
  cat("   [comp.type]: (OPTIONAL) Compartment of RNA eg. wc, nuc, cy. If not specified all types will be used one by one\n")
  cat("   [tail.type]: (OPTIONAL) Tail type of RNA eg. plus, minus. If not specified plus will be used\n")
}

# Process arguments
if (nargs < 3) {
  print.usage()
  q(save="no",status=1)
}

assoc.data.file <- args[[1]] # TF association data file
if (! file.exists(assoc.data.file)) {
  cat("TF association file ", assoc.data.file,"does not exist\n")
  q(save="no",status=1)
}

expr.file <- args[[2]] # Expression data file
if (! file.exists(expr.file)) {
  cat("Expression file ", expr.file,"does not exist\n")
  q(save="no",status=1)
}

output.dir <- args[[3]]
if (! file.exists(output.dir)) {
  dir.create(output.dir,recursive=T)
}

rna.type <- c("cage","rnaseq")
if (nargs > 3) {
  rna.type <- args[[4]]
}

compartment.type <- c("wc","cy","nuc")
if (nargs > 4) {
  compartment.type <- args[[5]]
}

tail.type <- "plus"
if (nargs > 5) {
  tail.type <- args[[6]]
}

# Initialize rulefit
platform <- "linux"
work.dir <- tempfile() # Default work directory
if (! file.exists(work.dir)) {
  dir.create(work.dir)
}
rfhome <- initialize.rulefit( work.dir=work.dir, rf.package.path=Sys.getenv("RULEFITBASE") )

for (ir in rna.type) {
  for (ic in compartment.type) {
    for (it in tail.type) {            
      if (two.stage) {
        curr.output.file <- file.path(output.dir, paste(ir,ic,it,"2stage","Rdata",sep=".") )
        # Create two-stage model
        cat("============\nComputing two-stage model\n============\n")
        two.stage.rulefit.results <- learn.tf.to.expr.rulefit.model(assoc.data=assoc.data.file,
                                                          expr=expr.file,
                                                          filter.expr=c(ir,ic,it),                                                          
                                                          two.stage.model=T,
                                                          low.expr.thresh=low.expr.thresh,
                                                          min.feature.support=min.feature.support,
                                                          rm.low.expr=rm.low.expr,
                                                          tree.size=tree.size,
                                                          test.reps=test.reps,
                                                          model.type=model.type)        
  
        # Create two-stage randomized model
        cat("============\nComputing two-stage randomized model\n============\n")
        rand.two.stage.rulefit.results <- learn.tf.to.expr.rulefit.model(assoc.data=assoc.data.file,
                                                          expr=expr.file,
                                                          filter.expr=c(ir,ic,it),                                                          
                                                          two.stage.model=T,
                                                          low.expr.thresh=low.expr.thresh,
                                                          min.feature.support=min.feature.support,                                                                         
                                                          rm.low.expr=rm.low.expr,                                                                         
                                                          randomize=1,
                                                          tree.size=tree.size,
                                                          test.reps=test.reps,
                                                          model.type=model.type)        
        
        # Save the results
        save(list=c("two.stage.rulefit.results",
                    "rand.two.stage.rulefit.results",
                    "assoc.data.file",
                    "expr.file",
                    "ir",
                    "ic",
                    "it"),
             file=curr.output.file)
        
      } else {
        curr.output.file <- file.path(output.dir, paste(ir,ic,it,"Rdata",sep=".") )
        # Create regression model
        cat("============\nComputing regression model\n============\n")
        rulefit.results <- learn.tf.to.expr.rulefit.model(assoc.data=assoc.data.file,
                                                          expr=expr.file,
                                                          filter.expr=c(ir,ic,it),                                                          
                                                          two.stage.model=F,
                                                          low.expr.thresh=low.expr.thresh,
                                                          min.feature.support=min.feature.support,
                                                          rm.low.expr=rm.low.expr,
                                                          tree.size=tree.size,
                                                          test.reps=test.reps,
                                                          model.type=model.type)        
        
        # Create randomized regression model
        cat("============\nComputing randomized regression model\n============\n")
        rand.rulefit.results <- learn.tf.to.expr.rulefit.model(assoc.data=assoc.data.file,
                                                          expr=expr.file,
                                                          filter.expr=c(ir,ic,it),                                                          
                                                          two.stage.model=F,
                                                          low.expr.thresh=low.expr.thresh,
                                                          min.feature.support=min.feature.support,
                                                          rm.low.expr=rm.low.expr,                                                                       
                                                          randomize=1,
                                                          tree.size=tree.size,
                                                          test.reps=test.reps,
                                                          model.type=model.type)                                                                       
        
        # Save the results
        save(list=c("rulefit.results",                  
                    "rand.rulefit.results",
                    "assoc.data.file",
                    "expr.file",
                    "ir",
                    "ic",
                    "it"),
             file=curr.output.file)        
        }            
    }
  }
}