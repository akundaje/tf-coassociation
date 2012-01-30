rm(list=ls())
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
setwd('/media/fusion10/work/encode/learning/combinatorics/src/scripts/')
source('assoc.matrix.utils.R')
platform = "linux"
rfhome <- initialize.rulefit(rf.package.path='/media/fusion10/work/encode/learning/combinatorics/src/code')
#load('/media/fusion10/work/encode/learning/combinatorics/data/GM12878/positiveSets/OA/SigMtrx_GM12878ATF3.OA.mtrx.Rdata')
#load('/media/fusion10/work/encode/learning/combinatorics/data/GM12878/positiveSets/OA/SigMtrx_GM12878POL3.OA.mtrx.Rdata')

# for (cell in c("GM12878","K562","Hepg2","Hela-S3","H1hesc")) {
#   #for (j in c("OA","OF","OT","OO")) {
#    for (j in c("OA")) {
#     for (i in list.files(path=file.path("/media/fusion10/work/encode/learning/combinatorics/results/randneg",cell,j),pattern=".*Rdata$",full.names=T)) {
#       tryCatch( plot.pairwise(i,output.dir=file.path("/media/fusion10/work/encode/learning/combinatorics/results/figures/randneg",cell,j)), error = function(e) cat("ERROR\n"))
# #      plot.int.strength(i,output.dir=file.path("/media/fusion10/work/encode/learning/combinatorics/results/figures/randneg",cell,j))
# #      plot.importance(i,output.dir=file.path("/media/fusion10/work/encode/learning/combinatorics/results/figures/randneg",cell,j))
#     }
#   }
# }

# for (cell in c("Hela-S3")) {
#   for (j in c("OA_average_relaxed")) {   
#     for (i in list.files(path=file.path("/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/randneg_iter",cell,j),pattern=".*Rdata$",full.names=T)) {
#       cat(get.file.parts(i)$fullname,"\n")
#       rulefit.results <- get.average.pairwise( get.average.int.strength( get.average.vi(i) ) )
#       save(list="rulefit.results", file=i)
#     }
#   }
# }

# cell <- args[[1]]
# # for (cell in c("GM12878","K562","Hepg2","Hela-S3","H1hesc")) {
#   for (j in c("OA_average_relaxed")) {
#     for (i in list.files(path=file.path("/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/randneg_iter",cell,j),pattern=".*Rdata$",full.names=T)) {
#       cat(get.file.parts(i)$fullname,"\n")
#       tryCatch( plot.average.pairwise.matrix(rulefit.results=i,output.dir=file.path("/media/fusion10/work/encode/learning/combinatorics/figures/TFCentric/randneg_iter",cell,j)), error = function(e) cat("ERROR\n"))
#       tryCatch( write.table.average.pairwise.matrix(rulefit.results=i,output.dir=file.path("/media/fusion10/work/encode/learning/combinatorics/figures/TFCentric/randneg_iter",cell,j)), error = function(e) cat("ERROR\n"))
#     }
#   }
# # }
# 
#cell <- args[[1]]
# odir <- "/media/fusion10/work/encode/learning/combinatorics/figures/TFCentric/randneg_iter"
# for (cell in c("GM12878","Hepg2","Hela-S3","H1hesc","K562")) {
#   for (j in c("OA_average")) {
#     for (i in list.files(path=file.path("/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/randneg_iter",cell,j),pattern=".*Rdata$",full.names=T)) {
#       cat(get.file.parts(i)$fullname,"\n")
#       tryCatch( plot.average.vi(rulefit.results=i,output.dir=file.path(odir,cell,j)), error = function(e) cat("ERROR\n"))
#       tryCatch( plot.average.int.strength(rulefit.results=i,output.dir=file.path(odir,cell,j)), error = function(e) cat("ERROR\n"))
#       #tryCatch( plot.average.pairwise(rulefit.results=i,output.dir=file.path(odir,cell,j)), error = function(e) cat("ERROR\n"))
#       tryCatch( plot.average.pairwise.matrix(rulefit.results=i,output.dir=file.path(odir,cell,j)), error = function(e) cat("ERROR\n"))
#       tryCatch( write.table.average.pairwise.matrix(rulefit.results=i,output.dir=file.path(odir,cell,j)), error = function(e) cat("ERROR\n"))      
#     }
#   }
# }

# for (i in list.files(path="/media/fusion10/work/encode/learning/combinatorics/results/GeneCentric/randneg_iter_average",pattern=".*Rdata$",full.names=T)) {
#   cat(get.file.parts(i)$fullname,"\n")
#   rulefit.results <- get.average.pairwise( get.average.int.strength( get.average.vi(i) ) )  
#   save(list="rulefit.results", file=i)
# }
# 
# for (i in list.files(path="/media/fusion10/work/encode/learning/combinatorics/results/GeneCentric/randneg_iter_average",pattern=".*Rdata$",full.names=T)) {
#   tryCatch( plot.average.vi(rulefit.results=i,output.dir='/media/fusion10/work/encode/learning/combinatorics/figures/GeneCentric/randneg_iter_average'), error = function(e) cat("ERROR\n"))
#   tryCatch( plot.average.int.strength(rulefit.results=i,output.dir='/media/fusion10/work/encode/learning/combinatorics/figures/GeneCentric/randneg_iter_average'), error = function(e) cat("ERROR\n"))
#   tryCatch( plot.average.pairwise(rulefit.results=i,output.dir='/media/fusion10/work/encode/learning/combinatorics/figures/GeneCentric/randneg_iter_average'), error = function(e) cat("ERROR\n"))
#   tryCatch( plot.average.pairwise.matrix(rulefit.results=i,output.dir='/media/fusion10/work/encode/learning/combinatorics/figures/GeneCentric/randneg_iter_average'), error = function(e) cat("ERROR\n"))
#   tryCatch( write.table.average.pairwise.matrix(rulefit.results=i,output.dir='/media/fusion10/work/encode/learning/combinatorics/figures/GeneCentric/randneg_iter_average'), error = function(e) cat("ERROR\n"))
# }

# input.dir <- "/media/fusion10/work/encode/learning/combinatorics/results/Histones/magnitude_average"
# for (i in list.files(path=input.dir,pattern=".*Rdata$",full.names=T)) {
#   rulefit.results <- get.average.pairwise( get.average.int.strength( get.average.vi(i) ) )
#   save(list="rulefit.results", file=i)
# }

# 
# for (cell in c("OA_average",
#                "OA_max_relaxed_average",
#                "OA_max_relaxed_tree6_average",
#                "OA_max_relaxed_trim_average",
#                "OA_max_relaxed_trim_tree10_average",
#                "OA_max_relaxed_trim_tree6_average",
#                "OA_max_relaxed_trim_tree8_average",
#                "OA_relaxed_average",
#                "OA_tree10_average",
#                "OA_tree6_average",
#                "OA_tree8_average",
#                "OF_average",
#                "OF_max_relaxed_average",
#                "OF_max_relaxed_tree6_average",
#                "OF_max_relaxed_trim_average",
#                "OF_max_relaxed_trim_tree6_average",
#                "OF_tree6_average",
#                "OO_average",
#                "OT_average"
#                )) {
#   input.dir <- sprintf("/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/randneg_iter/K562/average/%s",cell)
#   output.dir <- sprintf("/media/fusion10/work/encode/learning/combinatorics/figures/TFCentric/randneg_iter/K562/%s",cell)
# 
# #   input.dir <- "/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/randneg_iter/K562/average/OA_tree8_average"
# #   output.dir <- "/media/fusion10/work/encode/learning/combinatorics/figures/TFCentric/randneg_iter/K562/OA_tree8_average"
#   dir.create(path=output.dir,recursive=T)
#   for (i in list.files(path=input.dir,pattern=".*Rdata$",full.names=T)) {
#   #   tryCatch( plot.average.vi(rulefit.results=i,output.dir=output.dir,thresh=5), error = function(e) cat("ERROR\n"))
#   #   tryCatch( plot.average.int.strength(rulefit.results=i,output.dir=output.dir,thresh=1e-3), error = function(e) cat("ERROR\n"))
#     #tryCatch( plot.average.pairwise(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#     tryCatch( plot.average.pairwise.matrix(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#   #   tryCatch( write.table.average.pairwise.matrix(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#   }
# }

# for (cell in c("randneg_iter","null_append_randneg_iter")) {
# # for (cell in c("K562")) {  
# #   input.dir <- sprintf("/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/randneg_iter/%s/average/null_append_OA_max_tree6_average",cell)
# #   output.dir <- sprintf("/media/fusion10/work/encode/learning/combinatorics/figures/TFCentric/randneg_iter/%s/null_append_OA_max_tree6_average",cell)
#   input.dir <- sprintf("/media/fusion10/work/encode/learning/combinatorics/results/GeneCentric/%s/average/",cell)
#   output.dir <- sprintf("/media/fusion10/work/encode/learning/combinatorics/figures/GeneCentric/%s",cell)  
#   dir.create(path=output.dir,recursive=T)
#   for (i in list.files(path=input.dir,pattern=".*Rdata$",full.names=T)) {
#     tryCatch( plot.average.vi(rulefit.results=i,output.dir=output.dir,thresh=5), error = function(e) cat("ERROR\n"))
#     tryCatch( plot.average.int.strength(rulefit.results=i,output.dir=output.dir,thresh=1e-3), error = function(e) cat("ERROR\n"))
#     #tryCatch( plot.average.pairwise(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#     tryCatch( plot.average.pairwise.matrix(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#     tryCatch( write.table.average.pairwise.matrix(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#   }
# }

for (cell in c("GM12878","K562","Hela-S3","Hepg2","H1hesc")) {
  input.dir <- sprintf("/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/posneg_rand/%s",cell)
  output.dir <- sprintf("/media/fusion10/work/encode/learning/combinatorics/figures/TFCentric/posneg_rand/%s",cell)  
  dir.create(path=output.dir,recursive=T)
  for (i in list.files(path=input.dir,pattern=".*Rdata$",full.names=T)) {
    tryCatch( plot.importance(rulefit.results=i,output.dir=output.dir,filt.thresh=1e-5), error = function(e) cat("ERROR\n"))
    tryCatch( plot.pairwise.matrix(rulefit.results=i,output.dir=output.dir,use.null=T), error = function(e) cat("ERROR\n"))
  }
}
# input.dir <- "/media/fusion10/work/encode/learning/combinatorics/results/PWMs/all"
# for (i in list.files(path=input.dir,pattern=".*Rdata$",full.names=T)) {
#   rulefit.results <- get.average.pairwise( get.average.int.strength( get.average.vi(i) ) )
#   save(list="rulefit.results", file=i)
# }

# input.dir <- "/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/temp"
# for (i in list.files(path=input.dir,pattern=".*Rdata$",full.names=T)) {
#   rulefit.results <- get.average.pairwise( get.average.int.strength( get.average.vi(i) ) )
#   save(list="rulefit.results", file=i)
# }

# input.dir <- "/media/fusion10/work/encode/learning/combinatorics/results/TFCentric/temp"
# output.dir <- "/media/fusion10/work/encode/learning/combinatorics/figures/TFCentric/temp"
# for (i in list.files(path=input.dir,pattern=".*OA.*Rdata$",full.names=T)) {
#   tryCatch( plot.average.vi(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#   tryCatch( plot.average.int.strength(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#   tryCatch( plot.average.pairwise(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#   tryCatch( plot.average.pairwise.matrix(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
#   tryCatch( write.table.average.pairwise.matrix(rulefit.results=i,output.dir=output.dir), error = function(e) cat("ERROR\n"))
# }
