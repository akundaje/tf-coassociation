rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments

# Print usage
print.usage <- function(){
  cat("Rscript plot.biclusters.R [assocData] [outputFile]\n")
  cat("Plots co-binding map for a particular target TF domain\n")
  cat("   [assocData]: .Rdata file containing association matrix\n")
  cat("   [outputFile]: Path to output File (OPTIONAL) Default: inputDir/inputFile.pdf\n")
  cat("   [replaceFlag]: (OPTIONAL) If set to a 0, then if output file exits, the run will quit and not replace it, DEFAULT:0\n")  
}

if (nargs < 1) {
  print.usage()
  q(save="no",status=1)
}

assoc.data.file <- args[[1]] # Association matrix Rdata file (Contains list called assoc.data)
if (! file.exists(assoc.data.file)) {
  cat("Association Data file ",assoc.data.file,"does not exist\n")
  q(save="no",status=1)
}

output.file <- sprintf("%s.pdf",assoc.data.file)
if (nargs > 1) {
  output.file <- args[[2]] # Output File name
}

replace.flag <- F # do not replace existing file
if (nargs > 2) {
  replace.flag <- as.logical(args[[3]])
}

source("assoc.matrix.utils.R")

output.dir <- get.file.parts(output.file)$path
if (! file.exists(output.dir)) {
  dir.create(output.dir, recursive=T)
}

if ( (!replace.flag) & file.exists(output.file)) {
  cat("Output file exists: ",output.file,". Use replace.flag=1 to replace the file\n")
  q(save="no",status=0)
}

# Load association data (assoc.data)
load(assoc.data.file)
#colnames(assoc.data$assoc.matrix) <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "", toupper(colnames(assoc.data$assoc.matrix)), ignore.case=T)
#assoc.data$assoc.matrix <- filter.cols(assoc.data$assoc.matrix)
assoc.data$assoc.matrix <- filter.cols(assoc.data$assoc.matrix)
colnames(assoc.data$assoc.matrix) <- standardize.name(colnames(assoc.data$assoc.matrix))
#assoc.data$assoc.matrix <- scale(assoc.data$assoc.matrix)
clust.results <- plot.heatmap(t(assoc.data$assoc.matrix),
             to.file=output.file, 
             col.title="Target Sites", 
             row.title="TF peak strength", 
             title.name=assoc.data$target.name, 
             filt.thresh=NA, 
             replace.diag=F, 
             break.type="linear",
             clust.method=c("average","ward"),
             dist.metric=c("pearson","euclidean"),
             break.lowerbound=1e-4, # For TF-centric matrices
             break.upperbound=0.5,  # For TF-centric matrices
             #break.lowerbound=0.5, # For TF-centric matrices
             #break.upperbound=1.5,  # For TF-centric matrices
             row.optimal.order=T,
             col.optimal.order=F,                              
             scale="none")

clust.results$assoc.data.file <- assoc.data.file
clust.results$figure.file <- output.file
save(list="clust.results",file=paste(output.file,".clust.Rdata",sep=""))
