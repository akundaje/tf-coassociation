rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments

# Print usage
print.usage <- function(){
  cat("Rscript get.allTF.average.coassociation.R [inputDir] [outFilePrefix] [outputDir]\n")
  cat("Uses aggregate rulefit results from all focus factor contexts to obtain average and variance of interaction scores over all pairs of TFs\n")
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

source('assoc.matrix.utils.R')
require(reshape)

# Get list of Rdata files
all.Rdata.files <- list.files(path=input.dir, pattern=".*average.*Rdata$", full.names=T) # Get names of Rdata files
n.Files <- length(all.Rdata.files)
if (n.Files == 0) {
  q(save="no",status=0)
}

global.pairwise.int.table <- data.frame()

# Collapse scores across all focus domains
for (each.file in all.Rdata.files) {

  cat(get.file.parts(each.file)$fullname,"\n")
  load(each.file) # Will load rulefit results
  curr.target <- rulefit.results$target.name
  
  # Get all target-based pairwise interactions and plot it
  curr.pairwise <- melt.array( as.matrix(rulefit.results$mean.pairwise.int.matrix), varnames=c("p1","p2") )
  curr.pairwise$target <- curr.target
  
  # Add to table
  if (length(global.pairwise.int.table)==0) {
    global.pairwise.int.table <- curr.pairwise
  } else {
    global.pairwise.int.table <- rbind(global.pairwise.int.table, curr.pairwise)
  }

}

# For each TF compute partner scores across all target contexts
partners.over.all.focus <- list()
cat("Accumulating partners of each TF over all focus factor contexts ...\n")
for (p in unique(global.pairwise.int.table$p1)) {
  cat(p,"\n")
  temp.matrix <- as.data.frame( cast( global.pairwise.int.table[(global.pairwise.int.table$p1==p),] , 
                                      formula=target ~ p2 )
                                )
  rownames(temp.matrix) <- as.character(temp.matrix$target)
  temp.matrix$target <- NULL
  partners.over.all.focus[[p]] <- temp.matrix

  # Plot the matrix
  temp.output.dir <- file.path(output.dir, p)
  if (!file.exists(temp.output.dir)) {
    dir.create(temp.output.dir,recursive=T)
  }
  
  output.filename <- file.path( output.dir , p, sprintf("partners.over.all.focus.domains.matrix.%s.pdf",p) )  
  
  if (nrow(temp.matrix) > 0) {
    temp.matrix <- filter.cols(filter.rows(temp.matrix))
    rownames(temp.matrix) <- standardize.name(rownames(temp.matrix))
    colnames(temp.matrix) <- standardize.name(colnames(temp.matrix))
    tryCatch(
      clust.results <- plot.heatmap(data=temp.matrix, 
                                    use.as.dist=F,
                                    show.dendro="both",
                                    symm.cluster=F,
                                    to.file=output.filename,
                                    row.title="Focus TF domains",
                                    col.title=sprintf("Partners of %s",p),
                                    title.name=sprintf("Partners of %s over all focus factor domains",p),
                                    filt.thresh=1e-7,
                                    subtract.filt.thresh=T,
                                    #pseudo.count=1e-30,
                                    logval=T,
                                    replace.diag=T,
                                    replace.na=T,
                                    num.breaks=255,
                                    clust.method="average",
                                    dist.metric="euclidean",
                                    break.type="linear",
                                    break.lowerbound=10^4.5,
                                    #break.upperbound=5
                                    ),                   
      error = function(e) cat("ERROR ",p,"\n"))
  }
}

# Compute average of all pairs over all contexts
temp.table <- global.pairwise.int.table
temp.table$value <- log10(temp.table$value) - log10(1e-7)
temp.table$value[is.na(temp.table$value)] <- NA
temp.table$value[which(temp.table$value<=0)] <- NA

median.global.pairwise.matrix <- as.data.frame(cast(temp.table, formula=p1 ~ p2, fun.aggregate=function(x){median(x,na.rm=T)}))
rownames(median.global.pairwise.matrix) <- median.global.pairwise.matrix$p1
median.global.pairwise.matrix$p1 <- NULL
median.global.pairwise.matrix[!is.finite(as.matrix(median.global.pairwise.matrix))] <- NA

mean.global.pairwise.matrix <- as.data.frame(cast(temp.table, formula=p1 ~ p2, fun.aggregate=function(x){mean(x,na.rm=T)}))
rownames(mean.global.pairwise.matrix) <- mean.global.pairwise.matrix$p1
mean.global.pairwise.matrix$p1 <- NULL
mean.global.pairwise.matrix[!is.finite(as.matrix(mean.global.pairwise.matrix))] <- NA

mad.global.pairwise.matrix <- as.data.frame(cast(temp.table, formula=p1 ~ p2, fun.aggregate=function(x){mad(x,na.rm=T)}))
rownames(mad.global.pairwise.matrix) <- mad.global.pairwise.matrix$p1
mad.global.pairwise.matrix$p1 <- NULL
mad.global.pairwise.matrix[!is.finite(as.matrix(mad.global.pairwise.matrix))] <- NA

std.global.pairwise.matrix <- as.data.frame(cast(temp.table, formula=p1 ~ p2, fun.aggregate=function(x){sd(x,na.rm=T)}))
rownames(std.global.pairwise.matrix) <- std.global.pairwise.matrix$p1
std.global.pairwise.matrix$p1 <- NULL
std.global.pairwise.matrix[!is.finite(as.matrix(std.global.pairwise.matrix))] <- NA

max.global.pairwise.matrix <- as.data.frame(cast(temp.table, formula=p1 ~ p2, fun.aggregate=function(x){max(x,na.rm=T)}))
rownames(max.global.pairwise.matrix) <- max.global.pairwise.matrix$p1
max.global.pairwise.matrix$p1 <- NULL
max.global.pairwise.matrix[!is.finite(as.matrix(max.global.pairwise.matrix))] <- NA

robust.max.global.pairwise.matrix <- as.data.frame(cast(temp.table, formula=p1 ~ p2, fun.aggregate=function(x){quantile(x,probs=0.75,na.rm=T)}))
rownames(robust.max.global.pairwise.matrix) <- robust.max.global.pairwise.matrix$p1
robust.max.global.pairwise.matrix$p1 <- NULL
robust.max.global.pairwise.matrix[!is.finite(as.matrix(robust.max.global.pairwise.matrix))] <- NA

# Create output file names
out.Rdata.file <- file.path(output.dir, sprintf("%s.allTF.average.coassociation.Rdata",output.prefix))
cat("Saving data ..\n")
save(list=c("global.pairwise.int.table",
            "median.global.pairwise.matrix",
            "mean.global.pairwise.matrix",
            "mad.global.pairwise.matrix",
            "std.global.pairwise.matrix",
            "max.global.pairwise.matrix",
            "robust.max.global.pairwise.matrix",
            "partners.over.all.focus",
            "input.dir",
            "output.dir",
            "output.prefix"),
     file=out.Rdata.file)



# Plot the median matrix
cat("Plotting median matrix ..\n")
if (nrow(median.global.pairwise.matrix) > 0) {
  out.median.plot <- file.path(output.dir, sprintf("%s.allTF.median.coassociation.matrix.pdf",output.prefix))
  median.global.pairwise.matrix <- filter.cols(filter.rows(median.global.pairwise.matrix))
  rownames(median.global.pairwise.matrix) <- standardize.name(rownames(median.global.pairwise.matrix))
  colnames(median.global.pairwise.matrix) <- standardize.name(colnames(median.global.pairwise.matrix))
  median.global.pairwise.matrix <- as.matrix(median.global.pairwise.matrix)
  median.global.pairwise.matrix <- median.global.pairwise.matrix[ rownames(median.global.pairwise.matrix), rownames(median.global.pairwise.matrix)]
  diag(median.global.pairwise.matrix) <- NA
  median.global.pairwise.matrix[median.global.pairwise.matrix==0] <- NA
  median.global.pairwise.matrix <- as.data.frame(median.global.pairwise.matrix)
  
  clust.results <- plot.heatmap(data=median.global.pairwise.matrix, 
                                use.as.dist=F,
                                show.dendro="both",
                                symm.cluster=T,
                                to.file=out.median.plot,
                                row.title="TFs",
                                col.title="TFs",
                                title.name="Median pairwise scores across all contexts",
                                filt.thresh=NA,
                                subtract.filt.thresh=F,
                                pseudo.count=1e-30,
                                logval=F,
                                replace.diag=T,
                                replace.na=T,
                                num.breaks=255,
                                clust.method="ward",
                                dist.metric="euclidean",
                                break.type="linear",
                                break.lowerbound=5,
                                #break.upperbound=5
                                )                 
}


# Plot the mean matrix
cat("Plotting mean matrix ..\n")
if (nrow(mean.global.pairwise.matrix) > 0) {
  out.mean.plot <- file.path(output.dir, sprintf("%s.allTF.mean.coassociation.matrix.pdf",output.prefix))
  mean.global.pairwise.matrix <- filter.cols(filter.rows(mean.global.pairwise.matrix))
  rownames(mean.global.pairwise.matrix) <- standardize.name(rownames(mean.global.pairwise.matrix))
  colnames(mean.global.pairwise.matrix) <- standardize.name(colnames(mean.global.pairwise.matrix))
  mean.global.pairwise.matrix <- as.matrix(mean.global.pairwise.matrix)
  mean.global.pairwise.matrix <- mean.global.pairwise.matrix[ rownames(mean.global.pairwise.matrix), rownames(mean.global.pairwise.matrix)]
  diag(mean.global.pairwise.matrix) <- NA
  mean.global.pairwise.matrix[mean.global.pairwise.matrix==0] <- NA
  mean.global.pairwise.matrix <- as.data.frame(mean.global.pairwise.matrix)
  
  clust.results <- plot.heatmap(data=mean.global.pairwise.matrix, 
                                use.as.dist=F,
                                show.dendro="both",
                                symm.cluster=T,
                                to.file=out.mean.plot,
                                row.title="TFs",
                                col.title="TFs",
                                title.name="Mean pairwise scores across all contexts",
                                filt.thresh=NA,
                                subtract.filt.thresh=F,
                                pseudo.count=1e-30,
                                logval=F,
                                replace.diag=T,
                                replace.na=T,
                                num.breaks=255,
                                clust.method="ward",
                                dist.metric="pearson",
                                break.type="linear",
                                break.lowerbound=5,
                                #break.upperbound=5
                                )                 
}

# Plot the max matrix
cat("Plotting max matrix ..\n")
if (nrow(max.global.pairwise.matrix) > 0) {
  out.max.plot <- file.path(output.dir, sprintf("%s.allTF.max.coassociation.matrix.pdf",output.prefix))
  max.global.pairwise.matrix <- filter.cols(filter.rows(max.global.pairwise.matrix))
  rownames(max.global.pairwise.matrix) <- standardize.name(rownames(max.global.pairwise.matrix))
  colnames(max.global.pairwise.matrix) <- standardize.name(colnames(max.global.pairwise.matrix))
  max.global.pairwise.matrix <- as.matrix(max.global.pairwise.matrix)
  max.global.pairwise.matrix <- max.global.pairwise.matrix[ rownames(max.global.pairwise.matrix), rownames(max.global.pairwise.matrix)]
  diag(max.global.pairwise.matrix) <- NA
  max.global.pairwise.matrix[max.global.pairwise.matrix==0] <- NA
  max.global.pairwise.matrix <- as.data.frame(max.global.pairwise.matrix)
  
  clust.results <- plot.heatmap(data=max.global.pairwise.matrix, 
                                use.as.dist=F,
                                show.dendro="both",
                                symm.cluster=T,
                                to.file=out.max.plot,
                                row.title="TFs",
                                col.title="TFs",
                                title.name="Max pairwise scores across all contexts",
                                filt.thresh=NA,
                                subtract.filt.thresh=F,
                                pseudo.count=1e-30,
                                logval=F,
                                replace.diag=T,
                                replace.na=T,
                                num.breaks=255,
                                clust.method="ward",
                                dist.metric="pearson",
                                break.type="quantile",
                                break.lowerbound=5.5,
                                #break.upperbound=5.5
                                )                 
}

# Plot the robust max matrix
cat("Plotting robust max matrix ..\n")
if (nrow(robust.max.global.pairwise.matrix) > 0) {
  out.robust.max.plot <- file.path(output.dir, sprintf("%s.allTF.robust.max.coassociation.matrix.pdf",output.prefix))
  robust.max.global.pairwise.matrix <- filter.cols(filter.rows(robust.max.global.pairwise.matrix))
  rownames(robust.max.global.pairwise.matrix) <- standardize.name(rownames(robust.max.global.pairwise.matrix))
  colnames(robust.max.global.pairwise.matrix) <- standardize.name(colnames(robust.max.global.pairwise.matrix))
  robust.max.global.pairwise.matrix <- as.matrix(robust.max.global.pairwise.matrix)
  robust.max.global.pairwise.matrix <- robust.max.global.pairwise.matrix[ rownames(robust.max.global.pairwise.matrix), rownames(robust.max.global.pairwise.matrix)]
  diag(robust.max.global.pairwise.matrix) <- NA
  robust.max.global.pairwise.matrix[robust.max.global.pairwise.matrix==0] <- NA
  robust.max.global.pairwise.matrix <- as.data.frame(robust.max.global.pairwise.matrix)
  
  clust.results <- plot.heatmap(data=robust.max.global.pairwise.matrix, 
                                use.as.dist=F,
                                show.dendro="both",
                                symm.cluster=T,
                                to.file=out.robust.max.plot,
                                row.title="TFs",
                                col.title="TFs",
                                title.name="75 prcntile pairwise scores across all contexts",
                                filt.thresh=NA,
                                subtract.filt.thresh=F,
                                pseudo.count=1e-30,
                                logval=F,
                                replace.diag=T,
                                replace.na=T,
                                num.breaks=255,
                                clust.method="ward",
                                dist.metric="pearson",
                                break.type="linear",
                                break.lowerbound=5,
                                #break.upperbound=5
                                )                 
}


# Plot the mad matrix
cat("Plotting MAD matrix ..\n")
if (nrow(mad.global.pairwise.matrix) > 0) {
  out.mad.plot <- file.path(output.dir, sprintf("%s.allTF.mad.coassociation.matrix.pdf",output.prefix))
  mad.global.pairwise.matrix <- filter.cols(filter.rows(mad.global.pairwise.matrix))
  rownames(mad.global.pairwise.matrix) <- standardize.name(rownames(mad.global.pairwise.matrix))
  colnames(mad.global.pairwise.matrix) <- standardize.name(colnames(mad.global.pairwise.matrix))
  mad.global.pairwise.matrix <- as.matrix(mad.global.pairwise.matrix)
  mad.global.pairwise.matrix[mad.global.pairwise.matrix==0] <- NA
  #   mad.global.pairwise.matrix <- mad.global.pairwise.matrix[ rownames(mad.global.pairwise.matrix), rownames(mad.global.pairwise.matrix)]
  #   diag(mad.global.pairwise.matrix) <- NA
  mad.global.pairwise.matrix <- as.data.frame(mad.global.pairwise.matrix)
  
  clust.results <- plot.heatmap(data=mad.global.pairwise.matrix, 
                                use.as.dist=F,
                                show.dendro="both",
                                symm.cluster=T,
                                to.file=out.mad.plot,
                                row.title="TFs",
                                col.title="TFs",
                                title.name="Median Abs. Dev. of pairwise scores across all contexts",
                                filt.thresh=NA,
                                subtract.filt.thresh=F,
                                pseudo.count=1e-30,
                                logval=F,
                                replace.diag=F,
                                replace.na=T,
                                num.breaks=255,
                                clust.method="ward",
                                dist.metric="euclidean",
                                break.type="linear",
                                #break.lowerbound=0.5,
                                #break.upperbound=5
                                )                 
}

# Plot the std matrix
cat("Plotting STD matrix ..\n")
if (nrow(std.global.pairwise.matrix) > 0) {
  out.std.plot <- file.path(output.dir, sprintf("%s.allTF.std.coassociation.matrix.pdf",output.prefix))
  std.global.pairwise.matrix <- filter.cols(filter.rows(std.global.pairwise.matrix))
  rownames(std.global.pairwise.matrix) <- standardize.name(rownames(std.global.pairwise.matrix))
  colnames(std.global.pairwise.matrix) <- standardize.name(colnames(std.global.pairwise.matrix))
  std.global.pairwise.matrix <- as.matrix(std.global.pairwise.matrix)
  std.global.pairwise.matrix[std.global.pairwise.matrix==0] <- NA
  #   std.global.pairwise.matrix <- std.global.pairwise.matrix[ rownames(std.global.pairwise.matrix), rownames(std.global.pairwise.matrix)]
  #   diag(std.global.pairwise.matrix) <- NA
  std.global.pairwise.matrix <- as.data.frame(std.global.pairwise.matrix)
  
  clust.results <- plot.heatmap(data=std.global.pairwise.matrix, 
                                use.as.dist=F,
                                show.dendro="both",
                                symm.cluster=F,
                                to.file=out.std.plot,
                                row.title="TFs",
                                col.title="TFs",
                                title.name="Std Dev. of pairwise scores across all contexts",
                                filt.thresh=NA,
                                subtract.filt.thresh=F,
                                pseudo.count=1e-30,
                                logval=F,
                                replace.diag=F,
                                replace.na=T,
                                num.breaks=255,
                                clust.method="ward",
                                dist.metric="euclidean",
                                break.type="linear",
                                #break.lowerbound=1,
                                #break.upperbound=5
                                )                 
}
