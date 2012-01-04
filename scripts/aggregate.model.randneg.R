rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments
rm.target=F  # Set to true if you want to remove target variable from predictors
source('assoc.matrix.utils.R')

# Print usage
print.usage <- function(){
  cat("Rscript aggregate.model.randneg.R [inputDir]\n")
  cat("Aggregates variable importance, interaction strength and pairwise interactions over several random runs within a directory\n")
  cat("   [inputDir]: directory containing .Rdata files with results from different random negative sets\n")
}

if (nargs < 1) {
  print.usage()
  q(save="no",status=1)
}

input.dir <- args[[1]] # Directory containing rulefit.results Rdata files from several random negative set runs
if (! file.exists(input.dir)) {
  cat("Input Directory ", input.dir,"does not exist\n")
  q(save="no",status=1)
}

all.Rdata.files <- list.files(path=input.dir, pattern=".*Rdata$", full.names=T) # Get names of Rdata files
all.Rdata.files <- all.Rdata.files[! grepl("average",all.Rdata.files) ]
n.Files <- length(all.Rdata.files)
if (n.Files == 0) {
  q(save="no",status=0)
}
out.file <- gsub('iter[0-9]+', 'average', all.Rdata.files[[1]])

# RULEFIT RESULTS LIST (rulefit.results)
#  $rfmod: Rulefit model object
#  $assoc.classf.data: Rulefit classification data frame
#  $vi: variable importance DATA FRAME
#  $int.strength: interaction strength DATA FRAME
#  $pair.interactions: pairwise interactions DATA FRAME MATRIX

# Process first file
load(all.Rdata.files[[1]])
av.imp <- data.frame()
av.int.strength <- data.frame()
av.pair.interactions <- list()
av.cv <- data.frame()
pair.names <- rownames(rulefit.results$pair.interactions)
for ( i in pair.names ) {
  av.pair.interactions[[i]] <- data.frame()
}
target.name <- rulefit.results$dataset$target.name
rulefit.results$cv$lo <- NULL

for ( curr.file in all.Rdata.files ) {
  cat(curr.file,"\n")
  # load file
  load(curr.file)
  # skip if all vi are 0
  if ( all(rulefit.results$vi ==0) ) {
    next
  }
  av.imp <- rbind(av.imp, rulefit.results$vi)
  av.int.strength <- rbind(av.int.strength, rulefit.results$int.strength)
  if ( "cv" %in% names(rulefit.results) ) {
   av.cv <- rbind(av.cv, as.data.frame(rulefit.results$cv)) 
  }  
  for ( j in pair.names ) {
    av.pair.interactions[[j]] <- rbind(av.pair.interactions[[j]], rulefit.results$pair.interactions[j,])
  }
}

rulefit.results <- list(vi=av.imp,
                        int.strength=av.int.strength,
                        pair.interactions=av.pair.interactions,
                        cv=av.cv,
                        target.name=target.name)
save(list="rulefit.results",file=out.file)

# Compute average statistics from aggregated results
rulefit.results <- get.average.vi( rulefit.results )
rulefit.results <- get.average.int.strength( rulefit.results )
rulefit.results <- get.average.pairwise( rulefit.results )
rulefit.results <- get.average.cv( rulefit.results )
save(list="rulefit.results",file=out.file)
