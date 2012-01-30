rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments
source('assoc.matrix.utils.R')

# Print usage
print.usage <- function(){
  cat("Rscript global.posneg.factor.importance.R [inputDir] [outFilePrefix] [outputDir]\n")
  cat("Uses posneg rulefit results from several TFs to plot global importance matrix\n")
  cat("   [inputDir]: directory containing .Rdata files with aggregated rulefit results\n")
  cat("   [outFilePrefix]: Output file name prefix\n")
  cat("   [outputDir]: (OPTIONAL) output directory, DEFAULT: same as input directory\n")
}

# Process arguments
if (nargs < 2) {
  print.usage()
  q(save="no",status=1)
}

input.dir <- args[[1]] # Directory containing rulefit.results Rdata files from several random negative set runs
if (! file.exists(input.dir)) {
  cat("Input Directory ", input.dir,"does not exist\n")
  q(save="no",status=1)
}

output.prefix <- args[[2]]

output.dir <- input.dir
if (nargs > 2) {
  output.dir <- args[[3]]
  if (! file.exists(output.dir)) { dir.create(output.dir, recursive=T) }  
}

# Create output file names
out.Rdata.file <- file.path(output.dir, sprintf("%s.allTF.factor.importance.Rdata",output.prefix))
out.vi.file <- file.path(output.dir, sprintf("%s.allTF.factor.importance.matrix.xls",output.prefix))
out.vi.plot <- file.path(output.dir, sprintf("%s.allTF.factor.importance.pdf",output.prefix))

out.accuracy.plot <- file.path(output.dir, sprintf("%s.allTF.accuracy.pdf",output.prefix))
out.auc.plot <- file.path(output.dir, sprintf("%s.allTF.auc.pdf",output.prefix))
  
# Get list of Rdata files
all.Rdata.files <- list.files(path=input.dir, pattern=".*Rdata$", full.names=T) # Get names of Rdata files
n.Files <- length(all.Rdata.files)
if (n.Files == 0) {
  q(save="no",status=0)
}

all.vi.matrix <- data.frame()
all.accuracy <- data.frame(target.name="target",accuracy=0,AUC=0,stringsAsFactors=F)

count.vi <- 1
count.acc <- 1

for (each.file in all.Rdata.files) {
  cat(get.file.parts(each.file)$fullname,"\n")
  load(each.file) # Will load rulefit results
  curr.target <- rulefit.results$dataset$target.name
  
  # Update accuracy and AUC
  all.accuracy[count.acc,"target.name"] <- rulefit.results$dataset$target.name
  all.accuracy[count.acc,"accuracy"] <- 1 - rulefit.results$cv$errave
  all.accuracy[count.acc,"AUC"] <- 1 - rulefit.results$cv$omAUC
  count.acc <- count.acc + 1
  
  # Get all vi matrix and plot it
  curr.vi <- rulefit.results$vi
  rownames(curr.vi) <- curr.target
  if (nrow(curr.vi)>0) {
    curr.vi[1,curr.target] <- 100 #Add target variable
    if (length(all.vi.matrix)==0) { # Check if all.vi.matrix is empty
      all.vi.matrix <- curr.vi    
    } else {
      all.vi.matrix <- rbind(all.vi.matrix, curr.vi) # row bind vi values
    }
    count.vi <- count.vi + 1
  }
}

# Save the output Rdata
global.rulefit.results <- list(output.prefix=output.prefix,
                               all.vi.matrix=all.vi.matrix,
                               all.accuracy=all.accuracy)
save(list="global.rulefit.results",file=out.Rdata.file)

# Write the 3 matrices
write.table(as.matrix(all.vi.matrix), file=out.vi.file, quote=F, sep="\t", col.names=NA, na="")

# Plot the 3 matrices
all.vi.matrix <- filter.cols(filter.rows(all.vi.matrix))
rownames(all.vi.matrix) <- standardize.name(rownames(all.vi.matrix))
colnames(all.vi.matrix) <- standardize.name(colnames(all.vi.matrix))
if (nrow(all.vi.matrix) > 0) {
  plot.heatmap(data=all.vi.matrix,
               to.file=out.vi.plot,                               
               row.title="Target TF context",
               col.title="TFs",
               title.name="Conditional Partner TF importance",
               filt.thresh=NA,
               pseudo.count=1e-30,
               logval=F,
               replace.diag=F,
               replace.na=F,
               num.breaks=255,
               clust.method="ward",
               dist.metric="pearson",
               #break.lowerbound=30,
               #break.upperbound=80,                 
               break.type="linear")
}
