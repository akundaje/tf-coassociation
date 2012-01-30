rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments
source('assoc.matrix.utils.R')

# Print usage
print.usage <- function(){
  cat("Rscript get.allTF.ruleFit.results.R [inputDir] [outFilePrefix] [outputDir]\n")
  cat("Uses aggregate rulefit results from several randneg runs to plot global interaction figures and data files\n")
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
out.Rdata.file <- file.path(output.dir, sprintf("%s.allTF.globalVals.Rdata",output.prefix))

out.vi.file <- file.path(output.dir, sprintf("%s.allTF.factor.importance.matrix.xls",output.prefix))
out.vi.stats.file <- file.path(output.dir, sprintf("%s.allTF.factor.importance.stats.xls",output.prefix))
out.vi.plot <- file.path(output.dir, sprintf("%s.allTF.factor.importance.pdf",output.prefix))

out.scaled.vi.plot <- file.path(output.dir, sprintf("%s.allTF.scaled.factor.importance.pdf",output.prefix))

out.int.file <- file.path(output.dir, sprintf("%s.allTF.int.strength.matrix.xls",output.prefix))
out.int.stats.file <- file.path(output.dir, sprintf("%s.allTF.int.strength.stats.xls",output.prefix))
out.int.plot <- file.path(output.dir, sprintf("%s.allTF.int.strength.pdf",output.prefix))

out.pair.file <- file.path(output.dir, sprintf("%s.allTF.target.pair.interact.matrix.xls",output.prefix))
out.pair.stats.file <- file.path(output.dir, sprintf("%s.allTF.target.pair.interact.stats.xls",output.prefix))
out.pair.plot <- file.path(output.dir, sprintf("%s.allTF.target.pair.interact.pdf",output.prefix))

out.accuracy.plot <- file.path(output.dir, sprintf("%s.allTF.accuracy.pdf",output.prefix))
out.auc.plot <- file.path(output.dir, sprintf("%s.allTF.auc.pdf",output.prefix))
  
# Get list of Rdata files
all.Rdata.files <- list.files(path=input.dir, pattern=".*average.*Rdata$", full.names=T) # Get names of Rdata files
n.Files <- length(all.Rdata.files)
if (n.Files == 0) {
  q(save="no",status=0)
}

all.vi.matrix <- data.frame()
all.int.strength.matrix <- data.frame()
global.pairwise.int.matrix <- data.frame()
global.pairwise.int.matrix.scaled <- data.frame()
all.accuracy <- data.frame(target.name="target",accuracy=0,accuracy.lqr=0,accuracy.hqr=0,AUC=0,AUC.lqr=0,AUC.hqr=0,stringsAsFactors=F)

count.vi <- 1
count.int <-1
count.pair <-1
count.acc <- 1

for (each.file in all.Rdata.files) {
  cat(get.file.parts(each.file)$fullname,"\n")
  load(each.file) # Will load rulefit results
  curr.target <- rulefit.results$target.name
  
  # Update accuracy and AUC
  all.accuracy[count.acc,"target.name"] <- rulefit.results$target.name
  all.accuracy[count.acc,"accuracy"] <- 1 - rulefit.results$mean.cv["errave","mean.val"]
  all.accuracy[count.acc,"accuracy.lqr"] <- 1 - rulefit.results$mean.cv["errave","hqr"]
  all.accuracy[count.acc,"accuracy.hqr"] <- 1 - rulefit.results$mean.cv["errave","lqr"]
  all.accuracy[count.acc,"AUC"] <- 1 - rulefit.results$mean.cv["omAUC","mean.val"]
  all.accuracy[count.acc,"AUC.lqr"] <- 1 - rulefit.results$mean.cv["omAUC","hqr"]
  all.accuracy[count.acc,"AUC.hqr"] <- 1 - rulefit.results$mean.cv["omAUC","lqr"]
  count.acc <- count.acc + 1
  
  # Get all vi matrix and plot it
  curr.vi <- data.frame(rulefit.results$mean.vi$mean.val)
  rownames(curr.vi) <- rulefit.results$mean.vi$tf.name
  colnames(curr.vi) <- curr.target
  if (nrow(curr.vi)>0) {
    if (length(all.vi.matrix)==0) { # Check if all.vi.matrix is empty
      all.vi.matrix <- curr.vi    
    } else {
      curr.vi <- as.data.frame(curr.vi[rownames(all.vi.matrix),])
      rownames(curr.vi) <- rownames(all.vi.matrix)
      colnames(curr.vi) <- curr.target      
      all.vi.matrix <- cbind(all.vi.matrix, curr.vi) # column bind vi values
    }
    # Write to out.vi.stats.file
    curr.vi <- rulefit.results$mean.vi
    curr.vi$target.tf.context <-  curr.target # Add new column with target name
    rearranged.colnames <- c("target.tf.context", "tf.name", "mean.val", "std.val", "lqr", "hqr")
    curr.vi <- curr.vi[, rearranged.colnames]
    colnames(curr.vi) <- c("target.tf.context", "tf.name", "median.tf.importance", "std", "lower.quartile", "upper.quartile")
    if (count.vi==1) {
      write.table(curr.vi, file=out.vi.stats.file, quote=F, sep="\t", col.names=T, row.names=F, na="")
    } else {
      write.table(curr.vi, file=out.vi.stats.file, quote=F, sep="\t", col.names=F, append=T, row.names=F, na="")
    }
    count.vi <- count.vi + 1
  }
  
  # Get all int.strength matrix and plot it
  curr.int <- data.frame(rulefit.results$mean.int.strength$mean.val)
  rownames(curr.int) <- rulefit.results$mean.int.strength$tf.name
  colnames(curr.int) <- curr.target
  if (nrow(curr.int)>0) {
    if (length(all.int.strength.matrix)==0) {
      all.int.strength.matrix <- curr.int
    } else {
      curr.int <- as.data.frame(curr.int[rownames(all.int.strength.matrix),])
      rownames(curr.int) <- rownames(all.int.strength.matrix)
      colnames(curr.int) <- curr.target
      all.int.strength.matrix <- cbind(all.int.strength.matrix, curr.int)
    }
    # Write to out.int.stats.file
    curr.int <- rulefit.results$mean.int.strength
    curr.int$target.tf.context <-  curr.target # Add new column with target name
    rearranged.colnames <- c("target.tf.context", "tf.name", "mean.val", "std.val", "lqr", "hqr")
    curr.int <- curr.int[, rearranged.colnames]
    colnames(curr.int) <- c("target.tf.context", "tf.name", "median.int.strength", "std", "lower.quartile", "upper.quartile")
    if (count.int==1) {
      write.table(curr.int, file=out.int.stats.file, quote=F, sep="\t", col.names=T, row.names=F, na="")
    } else {
      write.table(curr.int, file=out.int.stats.file, quote=F, sep="\t", col.names=F, append=T, row.names=F, na="")
    }
    count.int <- count.int + 1
  }
  
  # Get all target-based pairwise interactions and plot it
  curr.pairwise <- rulefit.results$mean.pairwise.int.matrix[curr.target, ]
  if (! all(is.na(curr.pairwise))) {
    if (length(global.pairwise.int.matrix)==0){
      global.pairwise.int.matrix <- curr.pairwise
    } else {
      curr.pairwise <- curr.pairwise[,colnames(global.pairwise.int.matrix)]
      global.pairwise.int.matrix <- rbind(global.pairwise.int.matrix, curr.pairwise)
    }
    # Write to out.pair.stats.file
    curr.pairwise <- rulefit.results$aggregate.pairwise.interactions[[curr.target]]
    curr.pairwise$target.tf.context <-  curr.target # Add new column with target name
    rearranged.colnames <- c("target.tf.context", "tf.name", "mean.val", "std.val", "lqr", "hqr")
    curr.pairwise <- curr.pairwise[, rearranged.colnames]
    colnames(curr.pairwise) <- c("target.tf.context", "tf.name", "median.pairwise.int.strength", "std", "lower.quartile", "upper.quartile")
    if (count.pair==1) {
      write.table(curr.pairwise, file=out.pair.stats.file, quote=F, sep="\t", col.names=T, row.names=F, na="")
    } else {
      write.table(curr.pairwise, file=out.pair.stats.file, quote=F, sep="\t", col.names=F, append=T, row.names=F, na="")
    }
    count.pair <- count.pair + 1
  }
}

# Transpose all.vi.matrix and all.int.strength.matrix
all.vi.matrix <- as.data.frame(t(all.vi.matrix))
all.int.strength.matrix <- as.data.frame(t(all.int.strength.matrix))

# Save the output Rdata
global.rulefit.results <- list(output.prefix=output.prefix,
                               all.vi.matrix=all.vi.matrix,
                               all.int.strength.matrix=all.int.strength.matrix,
                               global.pairwise.int.matrix=global.pairwise.int.matrix,
                               all.accuracy=all.accuracy)
save(list="global.rulefit.results",file=out.Rdata.file)

# Write the 3 matrices
write.table(as.matrix(all.vi.matrix), file=out.vi.file, quote=F, sep="\t", col.names=NA, na="")
write.table(as.matrix(all.int.strength.matrix), file=out.int.file, quote=F, sep="\t", col.names=NA, na="")
write.table(as.matrix(global.pairwise.int.matrix), file=out.pair.file, quote=F, sep="\t", col.names=NA, na="")

# Plot the 3 matrices
# rownames(all.vi.matrix) <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "", toupper( gsub("K562b|Hepg2b", "B-",rownames(all.vi.matrix)) ), ignore.case=T)
# colnames(all.vi.matrix) <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "", toupper( gsub("K562b|Hepg2b", "B-",colnames(all.vi.matrix)) ), ignore.case=T)
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
               #dist.metric="spearman",
               #break.lowerbound=30,
               #break.upperbound=80,                 
               break.type="linear")
}

# rownames(all.int.strength.matrix) <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "", toupper( gsub("K562b|Hepg2b", "B-",rownames(all.int.strength.matrix)) ), ignore.case=T)
# colnames(all.int.strength.matrix) <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "", toupper( gsub("K562b|Hepg2b", "B-",colnames(all.int.strength.matrix)) ), ignore.case=T)
all.int.strength.matrix <- filter.cols(filter.rows(all.int.strength.matrix))
rownames(all.int.strength.matrix) <- standardize.name(rownames(all.int.strength.matrix))
colnames(all.int.strength.matrix) <- standardize.name(colnames(all.int.strength.matrix))
if (nrow(all.int.strength.matrix) > 0) {
  plot.heatmap(data=all.int.strength.matrix,
               to.file=out.int.plot,
               row.title="Target TF context",
               col.title="TFs",
               title.name="Conditional interaction strength",
               filt.thresh=1e-7,
               pseudo.count=0,
               logval=T,
               replace.diag=F,
               replace.na=F,
               num.breaks=255,
               clust.method="ward",
               break.lowerbound=6,
               break.type="linear")
}

# rownames(global.pairwise.int.matrix) <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "", toupper( gsub("K562b|Hepg2b", "B-",rownames(global.pairwise.int.matrix)) ), ignore.case=T)
# colnames(global.pairwise.int.matrix) <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "", toupper( gsub("K562b|Hepg2b", "B-",colnames(global.pairwise.int.matrix)) ), ignore.case=T)
global.pairwise.int.matrix <- filter.cols(filter.rows(global.pairwise.int.matrix))
rownames(global.pairwise.int.matrix) <- standardize.name(rownames(global.pairwise.int.matrix))
colnames(global.pairwise.int.matrix) <- standardize.name(colnames(global.pairwise.int.matrix))

if (nrow(global.pairwise.int.matrix) > 0) {
#   clust.results <- plot.heatmap(data=global.pairwise.int.matrix,
#                show.dendro="none",
#                symm.cluster=T,
#                to.file=out.pair.plot,
#                row.title="Target TFs",
#                col.title="Partners of target TF",
#                title.name="Global pairwise interactions",
#                filt.thresh=1e-7,
#                pseudo.count=1e-30,
#                logval=F,
#                replace.diag=T,
#                replace.na=T,
#                num.breaks=255,
#                clust.method="ward",
#                #dist.metric="spearman",
#                break.type="quantile",
#                break.lowerbound=1e-3,
#                #break.upperbound=1e-3
#                )               
  clust.results <- plot.heatmap(data=global.pairwise.int.matrix, 
                                use.as.dist=F,
                                show.dendro="both",
                                symm.cluster=T,
                                to.file=out.pair.plot,
                                row.title="Target TFs",
                                col.title="Partners of target TF",
                                title.name="Global pairwise interactions",
                                filt.thresh=1e-7,
                                pseudo.count=0,
                                logval=T,
                                replace.diag=T,
                                replace.na=T,
                                num.breaks=255,
                                clust.method="ward",
                                dist.metric="spearman",
                                break.type="linear",
                                break.lowerbound=4.5,
                                #break.upperbound=1e-3
                                )                 
}

# plot.heatmap(data=temp.matrix,row.cluster=F,col.cluster=F,
#                show.dendro="none",
#                #symm.cluster=T,
#                to.file="pair.matrix.png",#NULL,#out.pair.plot,
#                row.title="Target TFs",
#                col.title="Partners of target TF",
#                title.name="Global pairwise interactions",
#                filt.thresh=1e-7,
#                pseudo.count=1e-30,
#                logval=F,
#                replace.diag=T,
#                replace.na=T,
#                num.breaks=255,
#                clust.method="ward",
#                #dist.metric="spearman",
#                break.type="quantile",
#                break.lowerbound=1e-3,
#                #break.upperbound=1e-3
#                )
# Plot accuracy
require(ggplot2)
rownames(all.accuracy) <- all.accuracy$target.name
all.accuracy <- filter.rows(all.accuracy)
all.accuracy$target.name <- standardize.name(all.accuracy$target.name)

axes.format <- opts(plot.title = theme_text(size=12,vjust=1),                    
                    axis.text.x = theme_text(size=10,colour="black"),
                    axis.text.y = theme_text(size=10,colour="black",hjust=1),
                    axis.title.x = theme_text(size=12),
                    axis.title.y = theme_text(size=12,angle=90),
                    legend.title = theme_text(size=10,hjust=0),
                    legend.text = theme_text(size=10)                      
                    )

p <- ggplot(all.accuracy)
p <- p + 
  geom_bar( aes(x=reorder(target.name,accuracy),y=accuracy,fill=accuracy)) +
  geom_errorbar( aes(x=reorder(target.name,accuracy), ymax=accuracy.hqr, ymin=accuracy.lqr) )

axes.labels <- labs(x = "TF", y = "Accuracy") # axes labels
plot.title <- "Accuracy"
p <- p + axes.labels + axes.format + opts(title=plot.title) + coord_flip()
if (nrow(all.accuracy) > 50) {
  p <- p + opts(axis.text.y = theme_text(size=7,colour="black",hjust=1))
}

if (tolower(get.file.parts(out.accuracy.plot)$ext) == "png") {
  ggsave(file=out.accuracy.plot, plot=p, width=6, height=10, dpi=600)  
} else {
  ggsave(file=out.accuracy.plot, plot=p, width=6, height=10)  
}

# Plot AUC
p <- ggplot(all.accuracy)
p <- p + 
  geom_bar( aes(x=reorder(target.name,AUC),y=AUC,fill=AUC)) +
  geom_errorbar( aes(x=reorder(target.name,AUC), ymax=AUC.hqr, ymin=AUC.lqr) )

axes.labels <- labs(x = "TF", y = "AUC") # axes labels
plot.title <- "AUC"
p <- p + axes.labels + axes.format + opts(title=plot.title) + coord_flip()
if (nrow(all.accuracy) > 50) {
  p <- p + opts(axis.text.y = theme_text(size=7,colour="black",hjust=1))
}

if (tolower(get.file.parts(out.auc.plot)$ext) == "png") {
  ggsave(file=out.auc.plot, plot=p, width=6, height=10, dpi=600)  
} else {
  ggsave(file=out.auc.plot, plot=p, width=6, height=10)  
}
                    