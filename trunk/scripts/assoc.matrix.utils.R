# =======================================
# IMPORTANT
# =======================================
# YOU MUST RUN THE FOLLOWING BEFORE USING ANY OF THE OTHER FUNCTIONS IN THIS PACKAGE
# > rfhome <- initialize.rulefit(work.dir,rf.package.path)
# > platform <- "linux"

# ASSOCIATION FILE SPECIFICATIONS:
#  Assumes that association file name is of form [Prefix]_SigMtrx_[TargetName].[overlapType].mtrx
#  Requires header with Partner TFnames and [TargetName]
#  First column MUST be 'PeakID' representing peakids of the form [TargetName]_Pk_[peakId]_[overlapType]

# ASSOCIATION LIST (assoc.data)
# $assoc.matrix(DATAFRAME): dataFrame that is the association matrix (rows: TF sites, cols: partner TFs)
#    rownames are peakIDs. Rows are sorted by the numeric part of PeakIDs. First column is target
# $target.name(STRING): name of target
# $assoc.mtrx.file <- path of association .mtrx file
# $assoc.R.file <- path of association .Rdata file
# $expr.val <- expression values (if not valid it is assigned NA)

# CLASSIFICATION DATAFRAME
#  $x.vals : feature matrix
#  $y.vals : labels
#  $target.name : focus TF
#  $rm.target: if set to False, target.name is still part of x.vals

# RULEFIT RESULTS LIST (rulefit.results)
#  $rfmod: Rulefit model object
#  $dataset: Rulefit classification data frame
#  $vi: variable importance DATA FRAME
#  $int.strength: interaction strength DATA FRAME
#  $pair.interactions: pairwise interactions DATA FRAME MATRIX

# AGGREGATE RULEFIT RESULTS LIST (rulefit.results)
#  $vi: variable importance DATA FRAME (niter X TFs)
#  $mean.vi (DATA FRAME) rows are TFs, columns below
#     $mean.val : median values
#     $std.val : std deviations
#     $lqr: lower quartile
#     $hqr: upper quartile
#     $val.names: names of partners
#  $int.strength: interaction strength DATA FRAME (niter X TFs)
#  $mean.int.strength: (DATA FRAME) rows are TFs, columns below
#     $mean.val : median values
#     $std.val : std deviations
#     $lqr: lower quartile
#     $hqr: upper quartile
#     $val.names: names of partners
#  $pair.interactions: pairwise interactions LIST (TFs). Each element is a data frame (niter X TFs)
#  $aggregate.pairwise.interactions (a LIST containing data.frames for each TF) which has
#     $mean.val : median values
#     $std.val : std deviations
#     $lqr: lower quartile
#     $hqr: upper quartile
#     $val.names: names of partners
#  $mean.pairwise.int.matrix (data.frame) TFS X TFS  

# =================================================================================================================================
# =================================================================================================================================
# GENERIC UTILITIES
# =================================================================================================================================
# =================================================================================================================================

initialize.rulefit <- function( work.dir="/media/fusion10/work/encode/learning/combinatorics/src/rulefit" , rf.package.path=NA){
  # ===================================
  # Initializes a rulefit directory to start analyses
  # ===================================
  # work.dir: working directory where rulefit files will reside
  # rf.package.path: path to rulefit package
  
  if (is.na(rf.package.path)) {
    rf.package.path <- as.character( Sys.getenv("RULEFITBASE") )
  }
  
  if (! file.exists(work.dir)) {
    stop( sprintf("Rulefit working directory: %s does not exist", work.dir) )
  }

  if (! file.exists(rf.package.path)) {
    stop( sprintf("Rulefit code base directory: %s does not exist", rf.package.path) )
  }
  
  platform = "linux"
  if (! file.exists(work.dir)) {
    system(paste("mkdir",work.dir))    # Create working directory if it doesn't exist
  }
  system( sprintf( "rm -rf %s/*", work.dir))
  system( sprintf( "cp -r %s/* %s" , rf.package.path , work.dir ) ) # copy rulefit package to working directory
  rfhome <- work.dir
  source(file.path(rfhome,"rulefit.r"))
  library(akima,lib.loc=rfhome)
  return(rfhome)
}

restore.rf.model <- function( rulefit.results ){
  # ===================================
  # Restore a rulefit run
  # Returns rulefit.results
  # ===================================  
  # rulefit.results: Rdata file containing the rulefit.result OR rulefit.results object
  
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  rfrestore( rulefit.results$rfmod , x=rulefit.results$dataset$x.vals , y=rulefit.results$dataset$y.vals, wt=rep(1,nrow(rulefit.results$dataset$x.vals)) )
  return(rulefit.results)
}

get.file.parts <- function(file.fullpath) {
  # ===================================
  # Function will take a file name with path and split the file name into 
  # Returns
  # $path: full path (excluding the file part)
  # $fullname: path is removed
  # $name: path and extension is removed
  # $ext: extension only
  # ===================================    
  
  if (! is.character(file.fullpath)) {
  	stop('File name must be a string')
	}
	
  if ( grepl( paste(.Platform$file.sep,'$',sep=""), file.fullpath ) ) {
    return(list(path=gsub(paste(.Platform$file.sep,'$',sep="") , '' , file.fullpath),
                fullname='',
                name='',
                ext=''))
  }
  
	file.parts <- strsplit(as.character(file.fullpath), .Platform$file.sep, fixed=TRUE)[[1]] # split on file separator
	
	if (length(file.parts) == 0) { # if empty file name
		return(list(path='',
						fullname='',
						name='',
						ext='')
		)
	} else {
		if (length(file.parts) == 1) { # if no path then just the file name itself
			file.path <- '.'
			file.fullname <- file.parts
		} else {
			file.path <- paste(file.parts[1:(length(file.parts)-1)], collapse=.Platform$file.sep) # 1:last-1 token is path
			file.fullname <- file.parts[length(file.parts)] # last token is filename
		}        
		file.fullname.parts <- strsplit(file.fullname,'.',fixed=TRUE)[[1]] # split on .
		if (length(file.fullname.parts) == 1) { # if no extension
			file.ext <- ''
			file.name <- file.fullname.parts
		} else {
			file.ext <- paste('.', file.fullname.parts[length(file.fullname.parts)], sep="") # add the . to the last token
			file.name <- paste(file.fullname.parts[1:(length(file.fullname.parts)-1)], collapse=".")
		}
		return(list(path=file.path,
						fullname=file.fullname,
						name=file.name,
						ext=file.ext))
	}         	
} # end: get.file.parts()


filter.cols <- function(data,rm.treatments=F) {
  # ===================================
  # Remove specific columns from a data frame
  # ===================================  
  rem.cols <- c("GM12878OCT2m",
                "GM12878STAT1",
                "GM12878GCN5",
                "GM12878SPT20",
                "GM12878NFKBTNFa",
                "K562BRF2",
                "K562HEY1",
                "K562NR4A1",
                "K562XRCC4",
                "K562POL2S2",
                "HepG2HEY1",
                "HeLaS3GCN5",
                "HeLaS3SPT20",
                "Helas3CmycIFNg6h",
                "Helas3CmycIFNg30",
                "Helas3Cmyc",
                "expr.val")
  if (rm.treatments) {
    rem.cols <- c(rem.cols,
                  colnames(data)[grep("IFNa|IFNg",colnames(data),ignore.case=T)])
  }
  rm.idx <- match(rem.cols, colnames(data))
  rm.idx <- rm.idx[! is.na(rm.idx)]
  if (length(rm.idx) > 0) {
    data <- data[, -rm.idx ]
  }  
  return(data)
}


filter.rows <- function(data,rm.treatments=F) {
  # ===================================
  # Remove specific rows and columns from a data frame
  # ===================================  
  rem.rows <- c("GM12878OCT2m",
                "GM12878STAT1",
                "GM12878GCN5",
                "GM12878SPT20",
                "GM12878NFKBTNFa",
                "K562BRF2",
                "K562HEY1",
                "K562NR4A1",
                "K562XRCC4",
                "K562POL2S2",
                "HepG2HEY1",
                "HeLaS3GCN5",
                "HeLaS3SPT20",
                "Helas3CmycIFNg6h",
                "Helas3CmycIFNg30",
                "Helas3Cmyc",                
                "expr.val")
  if (rm.treatments) {
    rem.cols <- c(rem.cols,
                  colnames(data)[grep("IFNa|IFNg",colnames(data),ignore.case=T)])
  }
  rm.idx <- match(rem.rows, rownames(data))
  rm.idx <- rm.idx[! is.na(rm.idx)]
  if (length(rm.idx) > 0) {
    data <- data[ -rm.idx , ]  
  }  
  return(data)
}

standardize.name <- function(temp.names, conversion.file="name.conversion.tab") {
  # Converts arbitrary TF names to standardized names
  # temp.names: array of TF names
  # conversion.file: 3 columns, Col1(cell.line), Col2(arbit.name), Col3(std.name)
  
  conv.table <- read.table(file="name.conversion.tab",header=T,sep="\t",stringsAsFactors=F)
  conv.idx <- match(tolower(temp.names), tolower(conv.table$arbit.name)) # match names
  not.found <- which(is.na(conv.idx)) # find names that are not found in conv.table
  
  conv.names <- conv.table$std.name[conv.idx]
  conv.names[not.found] <- as.character(temp.names[not.found]) # Use original names for names not found in conv.table
  return(conv.names)
}

get.normal.score <- function(x) {
  # Converts values in x to normal scores (larger values get larger normal scores)
  y <- qnorm( ( rank(x, na.last="keep") - 0.375 ) / ( sum(!is.na(x)) + .25) ) 
  return(y)
}

plot.peaks.to.nearest.gene.distribution <- function(peak.distance.bed.file, proximal.cutoff=5000, distal.cutoff=50000){
  # =================================================
  # plots distribution of distance to nearest gene for all peaks in a dataset
  # peak.distance.bed.file: BED file containing peak coordinates and coordinates of nearest gene (peak_chr,peak_start,peak_stop,peak_id,tss_chr,tss_start,tss_stop,gene_name,distance)  
  # proximal.cutoff: distance cutoff for proximal
  # distal.cutoff: distance cutoff for distal
  # =================================================
  
  require(ggplot2)
#   peak.distance.bed.file <- "/media/fusion10/work/encode/learning/combinatorics/download/NearestGene/FirstClsstGene_KHHHG/K562/K562bGATA2.dist.bed.gz"
#   proximal.cutoff=5000 # prox.cutoff: distance cutoff for proximal 
#   distal.cutoff=50000 # distal.cutoff: distance cutoff for distal
  
  output.file <- gsub(pattern="\\.[^/]+$",replacement=".dist2tss.ecdf.png",x=peak.distance.bed.file)
  target.name <- gsub(pattern=".*/",replacement="",x=peak.distance.bed.file)
  distance.table <- read.table(file=peak.distance.bed.file,
                               header=T,
                               sep="\t",
                               #col.names=c("peak.chr","peak.start","peak.stop","peak.id","tss.chr","tss.start","tss.stop","tss.id","dist")
                               )
  
  distance.table$dist <- abs(distance.table$dist) + 1 # Convert 0s to 1
  distance.cutoff <- data.frame(prox=proximal.cutoff, dist=distal.cutoff)
  
  n.proximal <- sum(distance.table$dist <= proximal.cutoff+1,na.rm=T) # number of proximal sites
  prox.label <- sprintf("%d peaks",n.proximal)                               
  n.distal <- sum(distance.table$dist > distal.cutoff+1,na.rm=T) # number of distal sites
  distal.label <- sprintf("%d peaks",n.distal)
  n.sites <- nrow(distance.table)
  title.label <- sprintf("%s\n%d peaks", target.name, n.sites)

  temp.ecdf <- ecdf(distance.table$dist)
  dist.ecdf <- data.frame( x=knots(temp.ecdf), y=temp.ecdf(knots(temp.ecdf)) )
                               
  p <- ggplot(dist.ecdf) +
    geom_area(aes(x=x,y=y), color="orange", fill="orange", alpha=0.5) +
    geom_vline(data=distance.cutoff, aes(xintercept=prox), linetype=2, color="red") + 
    #geom_hline(aes(yintercept=temp.ecdf(proximal.cutoff)), linetype=2, color="red") +
    geom_vline(data=distance.cutoff, aes(xintercept=dist), linetype=2, color="blue") +
    #geom_hline(aes(yintercept=temp.ecdf(distal.cutoff)), linetype=2, color="blue") +
    #geom_text(aes(x=proximal.cutoff,y=0.2*n.sites,label=prox.label),color="red",hjust=1, vjust=0) +
    #geom_text(aes(x=distal.cutoff,y=0.7*n.sites,label=distal.label),color="blue",hjust=0, vjust=0) +    
    scale_x_log10() + 
    xlab("Distance from nearest TSS") + 
    ylab("Empirical CDF") + 
    theme_bw() +
    opts(title=title.label,
         axis.text.x = theme_text(size=14,colour="black"),
         axis.text.y = theme_text(size=14,colour="black",hjust=1),
         axis.title.x = theme_text(size=20),
         axis.title.y = theme_text(size=20,angle=90, vjust=0.3))
                               
  ggsave(filename=output.file,plot=p,width=5,height=5,dpi=600)
}

# =================================================================================================================================
# =================================================================================================================================
# TF CENTRIC UTILITIES
# =================================================================================================================================
# =================================================================================================================================

read.assoc.file <- function( assoc.file, std.thresh=NA, use.relaxed=T ) {
  # ===================================
  # Parses and reads an association table (needs to have headers for each column)
  #   First column MUST be 'PeakID' representing peakids of the form [TargetName]_Pk_[peakId]_[overlapType]
  #   sorts the rows by peakidx 
  #   remove low std. dev. columns/partners
  # Returns 
  # $assoc.matrix: dataFrame that is the association matrix (rows: TF sites, cols: partner TFs)
  # $target.name: name of target
  # ===================================  
  #assoc.file: association file ( Assumes that file name is of form [Prefix]_SigMtrx_[TargetName].[overlapType].mtrx )
  #std.thresh: columns with stddev. < str.thresh are removed from analysis
  #use.relaxed: If set to F, then all values <= 0 in the association matrix is set to 0
  
  assoc.matrix <- read.table( assoc.file, header=TRUE )
  target.name <- gsub( '(^.+SigMtrx_)|(\\.[^/]*mtrx(\\.gz)?$)' , '' , assoc.file ) # Remove prefix and suffix (.mtrx)
  
  rownames(assoc.matrix) <- assoc.matrix$PeakID # Peak names are in column named PeakID
  assoc.matrix$PeakID <- gsub( '(^.*Pk_)|(_[^_]+$)' , '' , assoc.matrix$PeakID ) # Convert PeakIds to numbers
  assoc.matrix <- assoc.matrix[ order( as.numeric( assoc.matrix$PeakID ) ) , ] # reorder rows by PeakId
  assoc.matrix <- assoc.matrix[ , -1 ] # Remove PeakId column  
  
  # Extract target column
  target.idx <- ( colnames(assoc.matrix) %in% target.name ) # Col index containing target.name
  target <- assoc.matrix[ , target.idx ] # target column
  assoc.matrix <- assoc.matrix[ , !target.idx ] # remove target column    

  # Processed relaxed norm peaks
  if (!use.relaxed) {
    assoc.matrix[(assoc.matrix <= 0)] <- 0
  }
  
  # Remove columns that have very low std
  if (! is.na( std.thresh )) {
    col.std <- apply( data.matrix(assoc.matrix) , 2 , sd ) # compute std for each column
    remove.col <- which( col.std < std.thresh ) # Find columns that have std < threshold
    if (length(remove.col) > 0) {
      assoc.matrix <- assoc.matrix[ , -remove.col]
    }
  }
  
  # Set target column as first column in partners matrix
  assoc.matrix <- cbind(target , assoc.matrix)
  colnames(assoc.matrix)[1] <- target.name
  
  # Reorder PeakIDs if they dont match target ranks
  peak.Ids <- as.numeric( gsub( '(^.*Pk_)|(_[^_]+$)' , '' , rownames(assoc.matrix) ) ) # Convert PeakIds to numbers
  order.target <- order(assoc.matrix[,target.name])
  n.rows <- length(peak.Ids)
  if ( (peak.Ids[1] == order.target[1]) & (peak.Ids[n.rows] == order.target[n.rows] ) ) {
    cat("Reversing labels\n")
    rownames(assoc.matrix) <- rev(rownames(assoc.matrix))
    peak.Ids <- as.numeric( gsub( '(^.*Pk_)|(_[^_]+$)' , '' , rownames(assoc.matrix) ) ) # Convert PeakIds to numbers
    assoc.matrix <- assoc.matrix[ order( peak.Ids ) , ] # reorder rows by PeakId
  }  
  
  # Remove bad columns
  assoc.matrix <- filter.cols(data=assoc.matrix)
  
  return( list(
    assoc.matrix=assoc.matrix,
    target.name=target.name ) )
} # end: read.assoc.file

read.expr.file <- function(expr.file){
  # ===================================
  # Read expression file (Col1: peakIdname [tab] Col2: ExprValue)
  # Returns
  #   $expr.val: Expression values (the rownames are peakId names)
  # ===================================  
  # expr.file: path to expression file (Col1: peakIdname [tab] Col2: ExprValue)
  
  expr.data <- read.table( file=expr.file , header=FALSE , row.names=1 )
  colnames(expr.data)[1] <- "expr.val"
  return(expr.data)
}

assoc.file.to.Rdata <- function( assoc.file, expr.file=NULL, output.dir=NULL, use.relaxed=T ) {
  # ===================================
  # Converts .mtrx file to a R object and saves in an Rdata file
  # ===================================  
  # assoc.file: association text file
  # expr.file: expression file (Col1: peakId, Col2: expression value)
  # output.dir: directory to store corresponding Rdata files
  # use.relaxed: If set to F, then all values <= 0 in the association matrix is set to 0
  if( is.null(output.dir) ) {
    output.dir <- get.file.parts(assoc.file)$path # if output.dir is not set make it equal to assoc.dir
  }

  assoc.data <- read.assoc.file(assoc.file,use.relaxed=use.relaxed)
  output.file <- file.path( output.dir , paste( get.file.parts(gsub("\\.gz$","",assoc.file))$fullname , '.Rdata' , sep="" ) )
  assoc.data$assoc.mtrx.file <- assoc.file
  assoc.data$assoc.R.file <- output.file
  
  if(! is.null(expr.file)) {
    expr.data <- read.expr.file(expr.file)
    assoc.data$assoc.matrix <- merge( assoc.data$assoc.matrix , expr.data , by="row.names" , all.x=TRUE , all.y=FALSE ) # Perform a join of assoc.data and expression data on peakIds
    rownames(assoc.data$assoc.matrix) <- assoc.data$assoc.matrix[,1] # reassign row names (since they get tranferred as col 1 after the join!)
    assoc.data$assoc.matrix <- assoc.data$assoc.matrix[,-1] # remove Column1
  } else {
    assoc.data$assoc.matrix$expr.val <- NA
  }

  save(list="assoc.data",file=output.file)
  
}

batch.read.assoc.file.to.Rdata <- function( assoc.dir , expr.file=NULL , output.dir=NULL, use.relaxed=T) {
  # ===================================
  # Reads all .mtrx files in a directory, 
  # converts them to R data frame and stores them
  # as .mtrx.Rdata files
  # ===================================  
  # assoc.dir: directory containing association files
  # expr.file: expression file (Col1: peakId, Col2: expression value)
  # output.dir: directory to store corresponding Rdata files
  # use.relaxed: If set to F, then all values <= 0 in the association matrix is set to 0
    
  # Search for all .mtrx files in assoc.dir
  assoc.file.paths <- dir( path=assoc.dir , pattern="\\.mtrx(\\.gz)?$" , full.names=TRUE , recursive=TRUE ) 
  
  for ( each.file in assoc.file.paths ) {
    cat("Processing file " , each.file , "\n")
    try( assoc.file.to.Rdata( each.file, expr.file, output.dir, use.relaxed ) , silent=T )
  }
}

truncate.OF.using.OA <- function( oa.dir, of.dir ) {
  # Use peak ids from OA file to truncate OF files
  oa.names <- list.files(oa.dir,pattern="*.mtrx.Rdata",full.names=T)
  of.names <- file.path( of.dir, gsub("\\.OA\\.",".OF.",basename(oa.names)) )
  count <- 0
  for (i in oa.names) {
    count <- count + 1
    if ( !file.exists(of.names[count]) ) {
      next
    } else {
      cat(sprintf("Processing %s\n",basename(i)))
      load(i)
      cat(dim(assoc.data$assoc.matrix),"\n")
      oa.peak.ids <- gsub("_OA","_OF",rownames(assoc.data$assoc.matrix))
      load(of.names[count])
      of.peak.ids <- rownames(assoc.data$assoc.matrix)
      assoc.data$assoc.matrix <- assoc.data$assoc.matrix[(of.peak.ids %in% oa.peak.ids) , ]
      cat(dim(assoc.data$assoc.matrix),"\n")
      save(list="assoc.data",file=of.names[count])
    }
  }
}

update.Rdata.with.expr.zscr <- function(assoc.Rdata.file, expr.zscr.file) {
  # ===================================
  # Update Rdata file with expression z-scores
  # ===================================  
  # assoc.Rdata.file: Rdata file containing assoc.data association data
  # expr.file: expression file (Col1: peakId, Col2: expression value)
  
  load(assoc.file) # loads assoc.data
  expr.data <- read.expr.file(expr.zscr.file)
    
  assoc.data$assoc.matrix <- merge( assoc.data$assoc.matrix , expr.data , by="row.names" ,all.x=TRUE , all.y=FALSE)
  rownames( assoc.data$assoc.matrix ) <- assoc.data$assoc.matrix[,1]
  assoc.data$assoc.matrix <- assoc.data$assoc.matrix[,-1]
  save( list="assoc.data" , file=assoc.file )
}

batch.update.Rdata.with.expr.zscr <- function(assoc.dir, expr.zscr.file) {
  # ===================================
  # Update all Rdata files in a directory with expression z-scores
  # ===================================  
  # assoc.dir: directory containing association files
  # expr.file: expression file (Col1: peakId, Col2: expression value)
  
 assoc.file.paths <- dir( path=assoc.dir , pattern="\\.mtrx\\.Rdata$" , full.names=TRUE , recursive=TRUE )
 for (each.file in assoc.file.paths) {
   cat("Processing file " , get.file.parts(each.file)$fullname , "\n")
   update.Rdata.with.expr.zscr( each.file, expr.zscr.file ) 
 }
}

randomize.assoc.matrix <- function(assoc.data, num.rand=1, rand.dim=2, change.row.names=T){
  # ===================================
  # Create randomized association matrix
  # Each column is randomized individually
  # rownames are assigned '_random' suffix
  # num.rand: number of random instantiations of the matrix that should be row bound to give the final random matrix
  # rand.dim: 0 (randomize rows and columns), 1 (randomize rows), 2 (randomize columns)
  # ===================================  
  # assoc.data$assoc.matrix
  # assoc.data$target.name
  
  # Load assoc.data if it is the name of a file
  if (is.character(assoc.data)) {
    load(assoc.data)    
  }
  
  if (rand.dim==0) {
    
    random.assoc.matrix <- as.data.frame( apply( assoc.data$assoc.matrix , 2 , sample ) ) # apply the 'sample' function to each column of assoc.matrix
    random.assoc.matrix <- as.data.frame( apply( random.assoc.matrix , 1 , sample ) ) # apply the 'sample' function to each row of assoc.matrix
    if (change.row.names) {
      rownames(random.assoc.matrix) <- paste( rownames(random.assoc.matrix) , 'random' , sep="_" ) # Add '_random' suffix to each row name
    } else {
      rownames(random.assoc.matrix) <- rownames(assoc.data$assoc.matrix)
    }
    
    colnames(random.assoc.matrix) <- colnames(assoc.data$assoc.matrix)
    
    if (num.rand > 1) {
      for (i in c(2:num.rand)) {
        temp.matrix <- as.data.frame( apply( assoc.data$assoc.matrix , 2 , sample ) )
        temp.matrix <- as.data.frame( apply( temp.matrix , 1 , sample ) )
        rownames(temp.matrix) <- paste( rownames(temp.matrix) , 'random' , sep="_" ) # Add '_random' suffix to each row name        
        colnames(temp.matrix) <- colnames(assoc.data$assoc.matrix)
        random.assoc.matrix <- rbind( random.assoc.matrix, temp.matrix  )
      }
    }
    
  } else {
    
    random.assoc.matrix <- apply( assoc.data$assoc.matrix , rand.dim , sample ) # apply the 'sample' function to each row/column of assoc.matrix
    if (rand.dim == 1) {
      random.assoc.matrix <- as.data.frame(t(random.assoc.matrix))
    } else {
      random.assoc.matrix <- as.data.frame(random.assoc.matrix)
      rownames(random.assoc.matrix) <- rownames(assoc.data$assoc.matrix)
    }
    if (change.row.names) {
      rownames(random.assoc.matrix) <- paste( rownames(random.assoc.matrix) , 'random' , sep="_" ) # Add '_random' suffix to each row name      
    } else {
      rownames(random.assoc.matrix) <- rownames(assoc.data$assoc.matrix)
    }    
    colnames(random.assoc.matrix) <- colnames(assoc.data$assoc.matrix)
    if (num.rand > 1) {
      for (i in c(2:num.rand)) {
        temp.matrix <- apply( assoc.data$assoc.matrix , rand.dim , sample )
        if (rand.dim == 1) {
          temp.matrix <- as.data.frame(t(temp.matrix))
        } else {
          temp.matrix <- as.data.frame(temp.matrix)
        }        
        rownames(temp.matrix) <- paste( rownames(temp.matrix) , 'random' , sep="_" ) # Add '_random' suffix to each row name
        colnames(temp.matrix) <- colnames(assoc.data$assoc.matrix)
        random.assoc.matrix <- rbind( random.assoc.matrix, temp.matrix  )
      }
    }
    
  }  
  # random.assoc.matrix[ , assoc.data$target.name ] <- -1 # Set target column to -1
  return(random.assoc.matrix)  
}

make.assoc.classf.rand.dataset <- function(assoc.data, rm.target=F, num.rand=1, trim.target=T, append.null=F, null.mode=2, null.replace=F){
  # ===================================
  # Create classification dataset based on random negative set
  # Returns
  #  $x.vals : feature matrix
  #  $y.vals : labels
  #  $target.name : focus TF
  #  $rm.target: if set of False, then target.name is part of x.vals
  # ===================================  
  # assoc.data$assoc.matrix
  # assoc.data$target.name
  # rm.target: if set to TRUE, target column is removed
  # num.rand: number of random instantiations of the matrix that should be row bound to give the final random matrix
  # trim.target: T/F (If set to T then all rows with target TF values < 0 are removed)
  # append.null: T/F (If set to T then randomized versions of each feature (col) of the association matrix is added as an extra feature)
  # null.mode: 0/1/2 (0: randomize rows and columns, 1: randomize rows independently, 2: randomize columns independently)
  # null.replace: T/F Set to T if you want to replace the true matrix with a randomized one
  
  if (is.character(assoc.data)) {
    load(assoc.data)    
  }
  
  assoc.data$assoc.matrix$expr.val <- NULL # Remove expression column
  
  # Append null
  if (append.null | null.replace) {
    if (! null.replace) {
      null.features <- randomize.assoc.matrix(assoc.data=assoc.data, rand.dim=null.mode)
      rownames(null.features) <- rownames(assoc.data$assoc.matrix)
      colnames(null.features) <- paste( colnames(null.features) , 'random' , sep="_" )
      assoc.data$assoc.matrix <- cbind(assoc.data$assoc.matrix, null.features)
    } else {
      null.features <- randomize.assoc.matrix(assoc.data=assoc.data, rand.dim=null.mode)
      rownames(null.features) <- rownames(assoc.data$assoc.matrix)
      colnames(null.features) <- colnames(assoc.data$assoc.matrix)
      assoc.data$assoc.matrix <- null.features
    }
  }
  
  # If trim.target remove rows for which target TF has negative values
  if (trim.target) {
    idx <- which(assoc.data$assoc.matrix[ , assoc.data$target.name ] >= 0)
    assoc.data$assoc.matrix <- assoc.data$assoc.matrix[ idx, ]
  }
  
  random.assoc.matrix <- randomize.assoc.matrix(assoc.data,num.rand=1) # Create random negative set  
  x.vals <- as.data.frame( rbind( assoc.data$assoc.matrix , random.assoc.matrix ) ) # put negative set below positive set  
  
  n.pos <- dim(assoc.data$assoc.matrix)[[1]] # number of positive examples
  y.vals <- x.vals[ , 1 ] # Create labels
  y.vals[c(1:n.pos)] <- 1
  y.vals[c( (n.pos+1) : (2*n.pos) )] <- -1

  x.vals <- filter.cols(x.vals) # Remove unwanted columns
  if (rm.target) {
    x.vals[ , assoc.data$target.name ] <- NULL # Remove target column
  }    
  
  return( list(
    x.vals=x.vals,
    y.vals=y.vals,
    rm.target=rm.target,
    target.name=assoc.data$target.name ) )
}

make.assoc.classf.posneg.dataset <- function(pos.assoc.data , neg.assoc.data, rm.target=T) {
  # ===================================
  # Create classification dataset based on user-defined negative set
  # Returns
  #  $x.vals : feature matrix
  #  $y.vals : labels
  #  $target.name : focus TF
  #  $rm.target: if set of False, then target.name is part of x.vals
  # ===================================  
  # pos.assoc.data: positive association data (list)
  # neg.assoc.data: negative association data (list)
  # rm.target: if set of False, then target.name is part of x.vals
  
  # assoc.data$assoc.matrix
  # assoc.data$target.name
  
  if (is.character(pos.assoc.data)) {
    load(pos.assoc.data)
    pos.assoc.data <- assoc.data
  }

  if (is.character(neg.assoc.data)) {
    load(neg.assoc.data)
    neg.assoc.data <- assoc.data
  }
  
  # Check that pos and neg set have same colnames and if not use intersection
  pos.names <- colnames(pos.assoc.data$assoc.matrix)
  neg.names <- colnames(neg.assoc.data$assoc.matrix)
  common.names <- intersect(pos.names, neg.names)
  if ( (length(pos.names) != length(common.names)) | (length(neg.names) != length(common.names)) ) {
    warning("All Positive and negative set partner TFs do not match. Using intersection")    
  }
  pos.assoc.data$assoc.matrix <- pos.assoc.data$assoc.matrix[,common.names]
  neg.assoc.data$assoc.matrix <- neg.assoc.data$assoc.matrix[,common.names]
  x.vals <- as.data.frame( rbind( pos.assoc.data$assoc.matrix, neg.assoc.data$assoc.matrix ) ) # put negative set below positive set
    
  # Set labels
  n.pos <- dim(pos.assoc.data$assoc.matrix)[[1]] # number of positive examples
  n.neg <- dim(neg.assoc.data$assoc.matrix)[[1]] # number of negative examples
  y.vals <- as.numeric(x.vals[,1])
  y.vals[c(1 : n.pos)] <- 1
  y.vals[c( (n.pos+1) : (n.pos+n.neg) )] <- -1
  
  x.vals <- filter.cols(x.vals) # Remove unwanted columns
  if (rm.target) {
    x.vals[ , pos.assoc.data$target.name ] <- NULL # Remove target column  
  }
  
  x.vals$expr.val <- NULL # Remove expression column
  
  return( list(
    x.vals=x.vals,
    y.vals=y.vals,
    rm.target=rm.target,
    target.name=pos.assoc.data$target.name ) )
}

make.tf.centric.tf.to.expr.dataset <- function(assoc.data,
                                    peak.distance.expr.bed.file,                                    
                                    rm.zero.expr=NA){
  # ===================================
  # Create Gene Centric expression dataset
  # ===================================  
  # assoc.data$assoc.matrix (For relaxed peak thresholds, binding values range from 0 to 2. For non-relaxed, range is 0 to 1)
  # assoc.data$target.name
  
  # Load TF data
  if (is.character(assoc.data)) {
    load(assoc.data)    
  }

  # Remove TF peaks for which there is no coassociated TF data
  x.vals <- filter.cols(assoc.data$assoc.matrix)
  x.vals <- x.vals[ (apply(x.vals,1,function(x) sum(x,na.rm=T))>0) , ]
  
  # Load expression data
  distance.table <- read.table(file=peak.distance.expr.bed.file,
                               header=T,
                               sep="\t",                               
                             )
  
  y.vals <- log2(distance.table$tss.cage+1)
  names(y.vals) <- as.character(distance.table$peak.id)
  
  # Remove zero valued expression data if required
  if (!is.na(rm.zero.expr)) {
    y.vals <- y.vals[y.vals>rm.zero.expr]
  }
  
  # Match gene names for expr and tf data
  common.gene.names <- intersect(names(y.vals),rownames(x.vals))
  x.vals <- x.vals[ match(common.gene.names,rownames(x.vals)) , ]
  y.vals <- y.vals[ match(common.gene.names,names(y.vals)) ]
      
  return(list(x.vals=x.vals,y.vals=y.vals,target.name=assoc.data$target.name))
}

run.rulefit <- function(assoc.classf.data, mode="class", corr.penalty=3, model.type="both", tree.size=6, test.reps=0){
  # ===================================
  # Run Rulefit
  # Returns rulefit object
  # ===================================
  # assoc.classf.data$x.vals : feature matrix
  # assoc.classf.data$y.vals : labels
  # assoc.classf$target.name : focus TF
  # assoc.class$rm.target: if set to TRUE then remove target variable
  # mode: "class" or "regress"
  
  rfmod <- rulefit(x=assoc.classf.data$x.vals ,
                   model.type=model.type,
                   #sparse=1,
                   inter.supp=corr.penalty,
                   xmiss=9e30, 
                   y=assoc.classf.data$y.vals , 
                   rfmode=mode,
                   #mod.sel=2,
                   max.rules=2000,
                   tree.size=tree.size,
                   test.reps=test.reps,
                   quiet=T
                   )  
  return(rfmod)
}

compute.rsquare <- function(y.true, y.pred) {
  # Computes coefficient of determination R-square
  r.square <- 1 - sum( (y.true - y.pred)^2 ) / sum( (y.true - mean(y.true))^2 )
  return(r.square)
}

run.cv.rulefit <- function(rulefit.results, nfold=10) {
  # ===================================
  # Run 10 fold cross-validation on rulefit results
  # Returns rulefit object
  # ===================================
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
    #  $rfmod: Rulefit model object
    #  $dataset: Rulefit classification data frame
    #  $vi: variable importance DATA FRAME
    #  $int.strength: interaction strength DATA FRAME
    #  $pair.interactions: pairwise interactions DATA FRAME MATRIX

  # Load rulefit.results if input is a data list
  rulefit.results <- restore.rf.model( rulefit.results )  
  rulefit.results$cv = rfxval (nfold=nfold, quiet=T)
  #rulefit.results$cv$lo <- NULL
  if ( any(grepl( pattern="rmse", x=names(rulefit.results$cv), fixed=T)) ) {
    rulefit.results$cv$rsquare <- 1 - ( rulefit.results$cv$rmse^2 / 
      mean((rulefit.results$dataset$y.vals - mean(rulefit.results$dataset$y.vals,na.rm=T))^2 ,na.rm=T) )
  }
  return(rulefit.results)
}

make.barplot <- function( vals , labels=NULL , sort.flag=TRUE , top.N=NULL , y.label=NULL , x.label=NULL, title.name="BARPLOT" , to.file=NULL ){
  # ===================================
  # Make horizontal bar plot
  # ===================================
  # vals: values for each bar height (numeric) or data frame
  # labels: optional labels for each bar
  # sort.flag: resort vals in ascending order
  # top.N: Only display top.N bars (sort.flag operates on vals before top.N)
  # title.name: title for bar plot
  # to.file: name of pdf output file to print figure to
  
  if ( (!is.data.frame(vals)) && (!is.matrix(vals)) ) { # If numeric data
    vals <- data.frame( vals=vals , names=c(1:length(vals)) ) # Convert to data frame with 2 columns    
  } else {
    vals <- as.data.frame(vals)
    colnames(vals) <- "vals" # Rename the column of the data frame to vals
    vals$names <- rownames(vals) # add an extra name column
  }
  
  if (!is.null(labels)) {
    vals$names <- labels # Add labels if necessary
  }

  if (sort.flag == TRUE) {
    resort.idx <- order(vals$vals)
    vals <- vals[resort.idx,]
  }
  
  if (! is.null(top.N)) {
    get.idx <- c( (nrow(vals)-top.N+1) : nrow(vals) )
    vals <- vals[get.idx,]    
  }    
  
  library(ggplot2)
  axes.format <- opts(plot.title = theme_text(size=12,vjust=1),                    
                      axis.text.x = theme_text(size=10,colour="black"),
                      axis.text.y = theme_text(size=10,colour="black",hjust=1),
                      axis.title.x = theme_text(size=12),
                      axis.title.y = theme_text(size=12,angle=90),
                      legend.title = theme_text(size=10,hjust=0),
                      legend.text = theme_text(size=10)                      
                      )
                      
  if (length(vals) == 0) {return()}
  
  if (sort.flag) {                      
    p1 <- ggplot(vals) + geom_bar( aes( x=reorder(names,vals) , y=vals ), stat="identity", fill=I("grey30") ) 
  } else {
    p1 <- ggplot(vals) + geom_bar( aes( x=names , y=vals ), stat="identity", fill=I("blue") )
  }
  axes.labels <- labs(x = y.label, y = x.label) # axes labels
  p1 <- p1 + axes.labels + axes.format + opts(title=title.name) + coord_flip()

  if (nrow(vals) > 50) {
    p1 <- p1 + opts(axis.text.y = theme_text(size=7,colour="black",hjust=1))
  }

  if (! is.null(to.file) ) {
    file.ext <- get.file.parts(to.file)$ext
    if (tolower(file.ext) == '.png') {      
      ggsave(file=to.file, plot=p1, width=4, height=10, dpi=600)      
    } else {
      ggsave(file=to.file, plot=p1, width=4, height=10)      
    }
  } else {
    p1
  }
}

get.var.imp <- function(rulefit.results, class=1){
  # ===================================
  # Computes variable importance from a rulefit model
  # Returns
  #   rulefit.results$rfmod
  #   rulefit.results$dataset
  #   rulefit.results$vi
  #   rulefit.results$int.strength
  #   rulefit.results$pair.interactions
  # ===================================  
  # rulefit.results$rfmod
  # rulefit.results$dataset
  # rulefit.results$vi
  # rulefit.results$int.strength
  # rulefit.results$pair.interactions
  # class: 1/0/-1, 1: positive class, -1: negative class, 0: all examples, or a specific set of examples
  
  if (is.character(rulefit.results)) {
    load(rulefit.results)
  }

  if (length(class)==1) {
    if (class==1){
      select.idx <- (rulefit.results$dataset$y.vals == 1)
    } else if (class == -1) {
      select.idx <- (rulefit.results$dataset$y.vals == -1)
    } else {
      select.idx <- is.finite(rulefit.results$dataset$y.vals)
    }
  } else {
    select.idx <- class
  }

  partner.names <- colnames(rulefit.results$dataset$x.vals)
  n.partners <- length(partner.names)
  vi <- varimp( impord=FALSE , plot=FALSE , x=rulefit.results$dataset$x.vals[ select.idx , ] )
  vi <- t(vi$imp)
  # Check if vi is all NaNs (happens when model fails due to very few examples)
  if (all(!is.finite(vi))) {
    vi <- as.data.frame( matrix(0, nrow=1, ncol=n.partners) )
  } else {
    vi <- as.data.frame(vi)
  }
  colnames(vi) <- partner.names
  rulefit.results$vi <- vi
  
  return(rulefit.results)
}

get.null.models <- function(rulefit.results, ntimes=10) {
  # Computes null models for a rulefit model
  # rulefit.results$rfmod
  # rulefit.results$dataset
  # rulefit.results$vi
  # rulefit.results$int.strength
  # rulefit.results$pair.interactions
  # Adds rulefit.results$null.models
  if (is.character(rulefit.results)) {
    load(rulefit.results)
  }
  
  #rulefit.results <- restore.rf.model(rulefit.results)
  cat("Computing null models ...\n");
  
  if ( any( names(rulefit.results) == "null.models" ) ) {
    rulefit.results$null.models <- intnull(ntimes, null.mods=rulefit.results$null.models, quiet=T)
  } else {
    rulefit.results$null.models <- intnull(ntimes, quiet=T)
  }    
  return(rulefit.results)
}

get.int.strength <- function( rulefit.results , plot=FALSE, use.null=F) {
  # ===================================
  # Add interaction strengths to rulefit results
  # Returns
  #   rulefit.results$rfmod
  #   rulefit.results$dataset
  #   rulefit.results$vi
  #   rulefit.results$int.strength
  #   rulefit.results$pair.interactions
  # ===================================
  # rulefit.results$rfmod
  # rulefit.results$dataset
  # rulefit.results$vi
  # rulefit.results$int.strength
  # rulefit.results$int.strength.null.ave (OPTIONAL)
  # rulefit.results$int.strength.null.std (OPTIONAL)
  # rulefit.results$pair.interactions
  # rulefit.results$pair.interactions.null.mean (OPTIONAL)
  # rulefit.results$pair.interactions.null.std (OPTIONAL)
  # rulefit.results$null.models (OPTIONAL)
  
  # plot: TRUE/FALSE/filename, filename: will save interaction strength plot to file
  
  # Load file is rulefit.results is a file name
  if (is.character(rulefit.results)) {
    load(rulefit.results)
  }
  
  partner.names <- colnames( rulefit.results$dataset$x.vals )
    
  if (use.null) {
    
    # Check if null models are computed. If not compute them
    if ( any( names(rulefit.results) == "null.models" ) ) {
      if (all(is.na(rulefit.results$null.models))) {
        rulefit.results <- get.null.models(rulefit.results)
      }
    } else {
      rulefit.results <- get.null.models(rulefit.results)
    }
    
    temp.int <- interact( partner.names, null.mods <- rulefit.results$null.models, plot=F)
    
    rulefit.results$int.strength <- as.data.frame( t(temp.int$int) )
    colnames(rulefit.results$int.strength) <- partner.names        
    
    rulefit.results$int.strength.null.mean <- as.data.frame( t(temp.int$nullave) )
    colnames(rulefit.results$int.strength.null.mean) <- partner.names
    
    rulefit.results$int.strength.null.std <- as.data.frame( t(temp.int$nullstd) )
    colnames(rulefit.results$int.strength.null.std) <- partner.names            
    
  } else {
    int.strength <- interact( partner.names , plot=F ) 
    rulefit.results$int.strength <- as.data.frame(t(int.strength))
    colnames(rulefit.results$int.strength) <- partner.names    
  }  
  
  if (is.logical(plot)) {
    if (plot) {
      title.name <- paste("PartnerTF Interaction strength:", rulefit.results$dataset$target.name)
      make.barplot( int.strength , partner.names , title.name=title.name )
    }
  } else {
      title.name <- paste("PartnerTF Interaction strength:", rulefit.results$dataset$target.name)
      make.barplot( int.strength , partner.names , title.name=title.name , to.file=plot)    
  }
    
  return(rulefit.results)
}

get.partner.pair.interactions <- function( rulefit.results,
                                           var.rank=1,
                                           var.idx=NULL,
                                           plot=FALSE,
                                           use.import=T,
                                           int.thresh=1e-7,
                                           pred.optim=50,
                                           use.null=F ){
  # ===================================
  # Get pairwise factor interactions for a particular partner
  # ===================================  
  # rulefit.results$rfmod
  # rulefit.results$dataset
  # rulefit.results$vi
  # rulefit.results$int.strength  
  # rulefit.results$int.strength.null.ave (OPTIONAL)
  # rulefit.results$int.strength.null.std (OPTIONAL)
  # rulefit.results$pair.interactions
  # rulefit.results$pair.interactions.null.mean (OPTIONAL)
  # rulefit.results$pair.interactions.null.std (OPTIONAL)
  # rulefit.results$null.models (OPTIONAL)
  
  # var.rank: rank (based on interaction strength) of partner TF to use to get interactions (ONLY used if var.idx=NULL)
  # var.idx: index or name of partner TF to use to get interactions
  # plot: TRUE/FALSE/filename, filename: will save interaction strength plot to file
  # use.import: T/F scales the interaction strengths by variable importance
  # pred.optim: if number of predictors is greater than pred.optim then only compute pairwise interaction scores for predictors with interaction potential > int.thresh
  # int.thresh: interaction threshold to use to consider pairwise interactions
  # use.null: if set to T then null models will be used to compute null interaction scores

  if (is.character(rulefit.results)) {
    load(rulefit.results)
  }
  
  partner.names <- colnames( rulefit.results$dataset$x.vals )
  
  # Initialize pair.interactions if necessary
  if ( ! any( names(rulefit.results) == "pair.interactions" ) ) {
    rulefit.results$pair.interactions <- data.frame(matrix( data=NA , nrow=length(partner.names) , ncol=length(partner.names) ) )
    rownames(rulefit.results$pair.interactions) <- partner.names
    colnames(rulefit.results$pair.interactions) <- partner.names      
  }

  # Compute interaction strengths if int.strength is all NA
  if ( all( is.na(rulefit.results$int.strength) ) ) {
    rulefit.results <- get.int.strength(rulefit.results, use.null=use.null)    
  }
  
  # Initialize rulefit.results$pair.interactions.null.mean and rulefit.results$pair.interactions.null.std if required  
  if (use.null) {    
    # Compute null models if required
    if ( any( names(rulefit.results) == "null.models" ) ) {
      if (all(is.na(rulefit.results$null.models))) {
        rulefit.results <- get.null.models(rulefit.results)
      }
    } else {
      rulefit.results <- get.null.models(rulefit.results)
    }
    
    # Initialize rulefit.results$pair.interactions.null.mean
    if ( ! any( names(rulefit.results) == "pair.interactions.null.mean" ) ) {
      rulefit.results$pair.interactions.null.mean <- rulefit.results$pair.interactions
    }
    
    # Initialize rulefit.results$pair.interactions.null.std
    if ( ! any( names(rulefit.results) == "pair.interactions.null.std" ) ) {
      rulefit.results$pair.interactions.null.std <- rulefit.results$pair.interactions
    }    
  }
  
  opt.order <- order( rulefit.results$int.strength , decreasing=TRUE ) # sort partners by decreasing interaction strength
  target.name <- rulefit.results$dataset$target.name
  
  # Get the predictor whose interactions you want to get
  if (is.null(var.idx)){
    var.idx <- opt.order[var.rank]
  }
  
  if (! is.numeric(var.idx)){  
    var.idx <- which(partner.names %in% var.idx)
  }
  
  # name of target partner TF
  var.name <- partner.names[var.idx]
    
  # Get a filtered set of predictors to compare to
  if ( ( length(partner.names) > pred.optim) && (! is.null(int.thresh) ) ) {
    valid.other.idx <- which(rulefit.results$int.strength >= int.thresh)    
  } else {
    valid.other.idx <- c(1:length(partner.names))    
  }
    
  other.idx <- setdiff(valid.other.idx,var.idx) # All other TFs
  
  if (use.null) {
    temp.int2var <- twovarint(var.idx, other.idx, plot=FALSE , import=use.import, null.mods=rulefit.results$null.models)
    int2var <- temp.int2var$int
  } else {
    int2var <- twovarint(var.idx, other.idx, plot=FALSE , import=use.import) 
  }  

  if ( ( length(partner.names) > pred.optim) && (! is.null(int.thresh) ) ) {
    topN <- sum(int2var >= int.thresh)
  } else {
    topN <- NULL
  }
    
  title.name=paste("Pairwise Interactions (Ui =",use.import,") of",var.name,"given",target.name)
  if (is.logical(plot)) {
    if (plot){      
      make.barplot( int2var , partner.names[other.idx] , top.N=topN , title.name=title.name )
    }
  } else {    
    make.barplot( int2var , partner.names[other.idx] , top.N=topN , title.name=title.name , to.file=plot)
  }
  
  rulefit.results$pair.interactions[ var.idx , other.idx ] <- int2var
  if (use.null) {
    rulefit.results$pair.interactions.null.mean[ var.idx , other.idx ] <- temp.int2var$nullave
    rulefit.results$pair.interactions.null.std[ var.idx , other.idx ] <- temp.int2var$nullstd     
  }
  
  return(rulefit.results)
}

get.all.partner.pair.interactions <- function(rulefit.results, use.import=T, int.thresh=1e-7, pred.optim=50, use.null=F) {
  # ===================================
  # Computes all pairwise interactions
  # ===================================  
  # rulefit.results$rfmod
  # rulefit.results$dataset
  # rulefit.results$vi
  # rulefit.results$int.strength
  # rulefit.results$int.strength.null.ave (OPTIONAL)
  # rulefit.results$int.strength.null.std (OPTIONAL)
  # rulefit.results$pair.interactions
  # rulefit.results$pair.interactions.null.mean (OPTIONAL)
  # rulefit.results$pair.interactions.null.std (OPTIONAL)
  # rulefit.results$null.models (OPTIONAL)

  # use.import: T/F scales the interaction strengths by variable importance
  # pred.optim: if number of predictors is greater than pred.optim then only compute pairwise interaction scores for predictors with interaction potential > int.thresh
  # int.thresh: interaction threshold to use to consider pairwise interactions
  # use.null: if set to T, then null models are used to compute null values of interaction strengths

  if (is.character(rulefit.results)) {
    load(rulefit.results)
  }
  
  num.partners <- length(rulefit.results$int.strength)
  if ( (num.partners > pred.optim) && (! is.null(int.thresh) ) ) {
    valid.idx <- which(rulefit.results$int.strength >= int.thresh)
  } else {
    valid.idx <- c(1:num.partners)
  }
  cat("Computing pairwise interactions for ",length(valid.idx), " of ", length(rulefit.results$int.strength), " predictors\n")
  for (vidx in valid.idx) {
    cat("\t",vidx,"..\n")
    rulefit.results <- get.partner.pair.interactions( rulefit.results, 
                                                      var.idx=vidx, 
                                                      plot=F, 
                                                      use.import=T, 
                                                      int.thresh=int.thresh, 
                                                      pred.optim=pred.optim,
                                                      use.null=use.null )
  }
  return(rulefit.results)
}

sample.randneg.rulefit.model <- function(assoc.data , rm.target=F, trim.target=T, append.null=F, null.mode=2, null.replace=F){
  # ===================================
  # Sample a rulefit model
  # (1) Creates a random negative set and returns a rulefit model for it
  # Returns
  #   $rfmod: rulefit model
  #   $dataset: sampled dataset
  #   $vi: variable importance (place holder data.frame of n.cols #partners)
  #   $int.strength: interaction strengths (placeholder data.frame of length #partners)
  #   $pair.interactions: pairwise interactions (placeholder data.frame of size #partners X #partners)
  # ===================================  
  # assoc.data$assoc.matrix
  # assoc.data$target.name
  # rm.target: if set to TRUE then target TF is not used in constructing the model
  # trim.target: T/F (If set to T then all rows with target TF values < 0 are removed)
  # append.null: T/F (If set to T then randomized versions of each feature (col) of the association matrix is added as an extra feature)
  # null.mode: 0/1/2 (0: randomize rows and columns, 1: randomize rows independently, 2: randomize columns independently)
  # null.replace: T/F Set to T if you want to replace the true matrix with a randomized one  
  
  assoc.classf.data <- make.assoc.classf.rand.dataset(assoc.data, 
                                                      rm.target=rm.target,
                                                      trim.target=trim.target,
                                                      append.null=append.null,
                                                      null.mode=null.mode,
                                                      null.replace=null.replace)
  ntrue <- (assoc.classf.data$y.vals == 1)
  rfmod <- run.rulefit(assoc.classf.data)
  
  # Create place holder for variable importance
  partner.names <- colnames(assoc.classf.data$x.vals)
  n.partners <- length(partner.names)
  vi <- as.data.frame( matrix( data=NA, nrow=1, ncol=n.partners) )
  colnames(vi) <- partner.names
  
  # Create place holder for interaction strengths
  int.strength <- data.frame(matrix( data=NA , nrow=1 , ncol=n.partners ))
  colnames(int.strength) <- partner.names
  
  # Create place holder for pairwise interactions
  pair.interactions <- data.frame(matrix( data=NA , nrow=n.partners , ncol=n.partners ))
  rownames(pair.interactions) <- partner.names
  colnames(pair.interactions) <- partner.names
  
  return( list(
    rfmod=rfmod,
    dataset=assoc.classf.data,
    vi=vi,
    int.strength=int.strength,
    pair.interactions=pair.interactions) )
}

# get.average.randneg.model <- function( assoc.data , iter=50 , plot=FALSE, rm.target=F ) {
#   # DEPRECATED FUNCTION
#   # ===================================
#   # (1) Create multiple random negative sets
#   # (2) Computes average factor importance over all sets
#   # (3) Returns the dataset and model whose factor importance is most correlated with average factor importance    
#   # Returns a list with variables
#   #  $rfmod: Rulefit model object
#   #  $assoc.classf.data: Rulefit classification data frame
#   #  $vi: variable importance DATA FRAME
#   #  $int.strength: interaction strength DATA FRAME
#   #  $pair.interactions: pairwise interactions DATA FRAME MATRIX
#   # ===================================  
#   # assoc.data$assoc.matrix
#   # assoc.data$target.name
#   # iter: number of negative sets to average over
#   
#   # Place holders for each sampled dataset and model
#   sampled.results <- list()
# 
#   # Run first iteration
#   curr.model <- get.var.imp( sample.randneg.rulefit.model( assoc.data, rm.target) )
#   final.rulefit.model <- curr.model # Initialize final selected model
#   curr.model$int.strength <- NULL # Remove int.strength
#   curr.model$pair.interactions <- NULL # Remove pair.interactions
#   avi <- curr.model$vi
#   curr.model$vi <- curr.model$vi / max(curr.model$vi + 1e-30) # Divide by max to get relative variable importance
#   sampled.results[[1]] <- curr.model
#   
#   for (i in 2:iter){
#     cat('Iteration ',i,'... \n')
#     curr.model <- get.var.imp( sample.randneg.rulefit.model( assoc.data, rm.target) )
#     curr.model$int.strength <- NULL # Remove int.strength
#     curr.model$pair.interactions <- NULL # Remove pair.interactions
#     avi <- avi + curr.model$vi
#     curr.model$vi <- curr.model$vi / max(curr.model$vi + 1e-30) # Divide by max to get relative variable importance
#     sampled.results[[i]] <- curr.model    
#   }
#   
#   avi <- avi / max(avi + 1e-30) # Divide by max to get relative variable importance
#   
#   # Create plot in necessary
#   if (is.logical(plot)) {
#     if (plot) {
#       title.name <- paste("Average PartnerTF Importance wrt", assoc.data$target.name)
#       make.barplot( as.numeric(avi) , colnames(avi) , title.name=title.name)
#     }    
#   } else {
#     title.name <- paste("Average PartnerTF Importance wrt", avi$target.name)
#     make.barplot( as.numeric(avi) , colnames(avi) , title.name=title.name , to.file=plot )
#   }
#   
#   # Get closest model and dataset  
#   max.corr <- 0
#   model.idx <- 1
#   for (i in 1:iter) {
#     # Check if std.dev of importance vectors are 0
#     if ( ( sd(as.numeric(avi)) == 0 ) | ( sd(as.numeric(sampled.results[[i]]$vi))==0 ) ) {
#       curr.corr <- 0
#     } else {
#       curr.corr <- abs(cor( as.numeric(avi) , as.numeric(sampled.results[[i]]$vi) , method="spearman" ))
#     }
#     if (curr.corr >= max.corr) {
#       max.corr <- curr.corr
#       model.idx <- i
#     }    
#   }
#   final.rulefit.model$rfmod <- sampled.results[[model.idx]]$rfmod
#   final.rulefit.model$dataset <- sampled.results[[model.idx]]$dataset
#   final.rulefit.model$vi <- avi
#   
#   return(final.rulefit.model)
# }

learn.posneg.rulefit.model <- function(pos.assoc.data , neg.assoc.data, rm.target=T){
  # ===================================
  # Sample a rulefit model
  # (1) Creates a random negative set and returns a rulefit model for it
  # Returns
  #   $rfmod: rulefit model
  #   $dataset: sampled dataset
  #   $vi: variable importance (place holder data.frame of n.cols #partners)
  #   $int.strength: interaction strengths (placeholder data.frame of length #partners)
  #   $pair.interactions: pairwise interactions (placeholder data.frame of size #partners X #partners)
  # ===================================  
  # assoc.data$assoc.matrix
  # assoc.data$target.name
  # rm.target: if set to TRUE then target TF is not used in constructing the model
    
  assoc.classf.data <- make.assoc.classf.posneg.dataset(pos.assoc.data, neg.assoc.data, rm.target)
  ntrue <- (assoc.classf.data$y.vals == 1)
  rfmod <- run.rulefit(assoc.classf.data,corr.penalty=1,tree.size=4)
  
  # Create place holder for variable importance
  partner.names <- colnames(assoc.classf.data$x.vals)
  n.partners <- length(partner.names)
  vi <- as.data.frame( matrix( data=NA, nrow=1, ncol=n.partners) )
  colnames(vi) <- partner.names
  
  # Create place holder for interaction strengths
  int.strength <- data.frame(matrix( data=NA , nrow=1 , ncol=n.partners ))
  colnames(int.strength) <- partner.names
  # int.strength expected null
  int.strength.null.mean <- int.strength
  # int.strength std. null
  int.strength.null.std <- int.strength
  
  # Create place holder for pairwise interactions
  pair.interactions <- data.frame(matrix( data=NA , nrow=n.partners , ncol=n.partners ))
  rownames(pair.interactions) <- partner.names
  colnames(pair.interactions) <- partner.names
  pair.interactions.null.mean <- pair.interactions
  pair.interactions.null.std <- pair.interactions
  
  return( list(
    rfmod=rfmod,
    dataset=assoc.classf.data,
    vi=vi,
    int.strength=int.strength,
    int.strength.null.mean=int.strength.null.mean,
    int.strength.null.std=int.strength.null.std,
    pair.interactions=pair.interactions,
    pair.interactions.null.mean=pair.interactions.null.mean,
    pair.interactions.null.std=pair.interactions.null.std) )
}
  
plot.heatmap <- function(data,
                         use.as.dist=F,
                         to.file=NULL, 
                         row.title="rows", 
                         col.title="cols", 
                         title.name=NULL, 
                         filt.thresh=1e-7, 
                         pseudo.count=1e-30, 
                         logval=F, 
                         replace.diag=T, 
                         replace.na=T,
                         num.breaks=255, 
                         break.type="quantile",
                         break.lowerbound=filt.thresh,
                         break.upperbound=NA,
                         clust.method="average", 
                         dist.metric="euclidean", 
                         scale="none", 
                         row.cluster=T, 
                         col.cluster=T,
                         symm.cluster=F,
                         show.dendro="both") {
  # ===================================
  # Plot clustered heatmap of associations
  # data: any data frame (Rows: are binding sites, Cols: partner TFs)
  # use.as.dist: use data directly as a similarity matrix (data must be symmetric matrix/data frame)
  # to.file: png/pdf file that you want to save the figure to (default: no saving)
  # row.title: axis title for rows
  # col.title: axis title for columns
  # title.name: plot title
  # filt.thresh: used to remove rows and cols with all values < filt.thresh
  # pseudo.count: uniform random numbers scaled by pseudo.count are added to the matrix to avoid 0 std for constant columns
  # logval: T/F . If set to T, then the matrix is log transformed before clustering. (filt.thresh will also be log transformed)
  # replace.diag: If set to T, then matrix diagonal values are replaced by maximum value in the matrix
  # replace.na: If set of T, then NA values are replaced by minimum value in the matrix
  # num.breaks: number of breaks in colors (The breaks correspond to uniformly sampled quantiles from the distribution of values in the matrix, excluding all values below filt.thresh)
  # break.type: type of color breaks, quantile: means the colors are adjusted to uniformly spaced quantiles, linear: colors are placed on the linear scale
  # break.lowerbound: For values below break.lowerbound are ignored and set to the lowest color
  # break.upperbound: For values above break.upperbound are ignored and set to the highest color
  # clust.method: linkage method e.g. "complete/average/ward/single"
  # dist.metric="euclidean/pearson/spearman/binary/manhattan"  
  # scale: "row", "col", "none" whether to standardize rows/columns or none
  # row.cluster=T : T or F to cluster rows OR a numeric vector with the desired row order OR a dendrogram object
  # col.cluster=T : T or F to cluster columns OR a numeric vector with the desired row order OR a dendrogram object
  # symm.cluster=F : if set to T then column clustering is set equal to row clustering
  # show.dendro="both"  : which dendrograms to show "row", "column", "both" or "none"
  # NOTE: This is currently horrendously slow for large number of rows
  # TODO Improvements:
  # (1) Don't cluster rows/cols if there are too many  
  # ===================================
  library(fastcluster)
  library(gplots)
  library(fBasics)
  
  # Remove columns with all NAs
  #data <- filter.cols(data)
  if (!symm.cluster) {
    na.idx <- apply( data , 2 , function(x) all(is.na(x)) )
    data <- data[ , !na.idx]
    # Remove rows with all NAs
    na.idx <- apply( data , 1 , function(x) all(is.na(x)) )
    data <- data[!na.idx, ]  
  
    # Remove columns with all very small values
    if (! is.na(filt.thresh) ) {
      na.idx <- apply( data , 2 , function(x) all((x<filt.thresh),na.rm=T) )
      data <- data[ , !na.idx]
      # Remove rows with all very small values
      na.idx <- apply( data , 1 , function(x) all((x<filt.thresh),na.rm=T) )
      data <- data[!na.idx, ]
    }
  }
  data.size <- dim(data)
  if (data.size[[1]] < 1) {return()}
  
  # Add a small random number to each value to avoid problems with clustering
  clean.data <- as.matrix(data) + (pseudo.count * matrix( data=runif(prod(data.size)), data.size[1], data.size[2]) )
  if (logval) {
    clean.data <- log10(clean.data)
    clean.data[is.infinite(clean.data)] <- NA    
    if (!is.na(filt.thresh)) { 
      filt.thresh <- log10(filt.thresh)
      clean.data <- clean.data - filt.thresh
      clean.data[which(clean.data < 0)] <- 0
      }
  }
  breaks.data <- as.vector(clean.data)
  min.val <- min(breaks.data,na.rm=T)
  max.val <- max(breaks.data,na.rm=T)  
  
  # Replace diagonal values with max if required
  if (replace.diag) { 
    for (r in rownames(clean.data)) {
      if (r %in% colnames(clean.data)) {
        if (is.na(clean.data[r,r])) { clean.data[r,r] <- max.val }
      }
    }
  }
  
  # Replace NAs with minimum value
  if (replace.na) { clean.data[is.na(clean.data)] <- min.val } # set NAs to minimum value
  
  # Generate color breaks
  # number of parts to split the color map into (max 3 parts min.val:break.lowerbound , lowerbound:upperbound, upperbound:max.val)
  if (!is.na(break.lowerbound)) {
    if (break.lowerbound < min.val) {break.lowerbound <- NA}
  }
  if (!is.na(break.upperbound)) {
    if (break.upperbound > max.val) {break.upperbound <- NA}
  }
  
  n.breaks <- num.breaks
  
  if (is.na(break.lowerbound) && is.na(break.upperbound)) { # Full scale
    
    if (break.type == "quantile") {
      breaks.vals <- quantile(breaks.data,prob=seq(0,1,length.out=n.breaks),na.rm=T)
      breaks.vals <- c(min.val,breaks.vals,max.val)
    } else if (break.type == "linear") {
      breaks.vals <- seq(min.val,max.val,length.out=n.breaks)
    }
      
  } else if (is.na(break.lowerbound) && !is.na(break.upperbound)) {  # min:upper , upper:max    
    
    n.1.2 <- round(2*n.breaks/3) # number of break values
    n.3 <- round(n.breaks/3)
    if (break.type == "quantile") {
      # min:upper
      breaks.vals.1.2 <- quantile( breaks.data[breaks.data<break.upperbound], prob=seq(0,1,length.out=n.1.2), na.rm=T)
      breaks.vals.1.2 <- c(min.val, breaks.vals.1.2, break.upperbound)
      #upper:max      
      breaks.vals.3 <- quantile( breaks.data[breaks.data>=break.upperbound], prob=seq(0,1,length.out=n.3), na.rm=T)
      breaks.vals.3 <- c(breaks.vals.3, max.val)      
      breaks.vals <- c(breaks.vals.1.2, breaks.vals.3)
    } else if (break.type == "linear") {
      breaks.vals <- c( seq(min.val, break.upperbound, length.out=n.1.2),
                        seq(break.upperbound, max.val, length.out=n.3))
    }
    
  } else if (!is.na(break.lowerbound) && is.na(break.upperbound)) {
    
    n.1 <- round(n.breaks/3) # number of break values
    n.2.3 <- round(2*n.breaks/3)
    if (break.type == "quantile") {
      # min:lower
      breaks.vals.1 <- quantile( breaks.data[breaks.data<break.lowerbound], prob=seq(0,1,length.out=n.1), na.rm=T)
      breaks.vals.1 <- c(min.val, breaks.vals.1, break.lowerbound)
      #lower:max      
      breaks.vals.2.3 <- quantile( breaks.data[breaks.data>=break.lowerbound], prob=seq(0,1,length.out=n.2.3), na.rm=T)
      breaks.vals.2.3 <- c(breaks.vals.2.3, max.val)      
      breaks.vals <- c(breaks.vals.1, breaks.vals.2.3)
    } else if (break.type == "linear") {
      breaks.vals <- c( seq(min.val, break.lowerbound, length.out=n.1),
                        seq(break.lowerbound, max.val, length.out=n.2.3))
    }    
    
  } else {
    
    n.1 <- round(n.breaks/3) # number of break values
    n.2 <- round(n.breaks/3) 
    n.3 <- round(n.breaks/3) 
    if (break.type == "quantile") {
      # min:lower
      breaks.vals.1 <- quantile( breaks.data[breaks.data<break.lowerbound], prob=seq(0,1,length.out=n.1), na.rm=T)
      breaks.vals.1 <- c(min.val, breaks.vals.1, break.lowerbound)
      # lower:upper
      breaks.vals.2 <- quantile( breaks.data[ (breaks.data>=break.lowerbound) & (breaks.data<break.upperbound)], 
                                 prob=seq(0,1,length.out=n.2), na.rm=T)
      breaks.vals.1 <- c(breaks.vals.2, break.upperbound)      
      # upper:max      
      breaks.vals.3 <- quantile( breaks.data[breaks.data>=break.upperbound], prob=seq(0,1,length.out=n.3), na.rm=T)
      breaks.vals.3 <- c(breaks.vals.3, max.val)      
      breaks.vals <- c(breaks.vals.1, breaks.vals.2, breaks.vals.3)
    } else if (break.type == "linear") {
      breaks.vals <- c( seq(min.val, break.lowerbound, length.out=n.1),
                        seq(break.lowerbound, break.upperbound, length.out=n.2),
                        seq(break.upperbound, max.val, length.out=n.3))
    }    
    
  }

  all.colors <- seqPalette( (length(breaks.vals)-1) , "YlOrRd" )

  # Old code  
#   temp.min.val <- min.val
#   temp.max.val <- max.val
#   if (! is.na(break.lowerbound)) {  
#     breaks.data <- breaks.data[breaks.data>break.lowerbound] # Remove low values
#     temp.min.val <- break.lowerbound
#   }
#   if (! is.na(break.upperbound)) {  
#     breaks.data <- breaks.data[breaks.data<break.upperbound] # Remove low values
#     temp.max.val <- break.upperbound
#   }   
#     
#   n.breaks <- num.breaks
#   if (break.type == "quantile") {
#     breaks.vals <- quantile(breaks.data,prob=seq(0,1,length.out=n.breaks),na.rm=T)
#     breaks.vals <- c(temp.min.val,breaks.vals,temp.max.val)
#   } else if (break.type == "linear") {
#     breaks.vals <- seq(temp.min.val,temp.max.val,length.out=n.breaks)
#   }
#   
#   #all.colors <- heat.colors(length(breaks.vals)-3)
#   #all.colors <- heatPalette(length(breaks.vals)-3)
#   #all.colors <- rev( divPalette( (length(breaks.vals)-3) , "RdBu" ) )  
#   #all.colors <- rev(focusPalette( (length(breaks.vals)-3) , "redfocus" ))
#   #all.colors <- rampPalette( (length(breaks.vals)-3) , "blue2red" )
#   #all.colors <- rev( redgreen((length(breaks.vals)-3)) )
#   all.colors <- seqPalette( (length(breaks.vals)-3) , "YlOrRd" )
#   all.colors <- ( c(all.colors[1], all.colors, all.colors[length(all.colors)]) )
  
  # Select image size
  if ( !is.null(to.file) ) {
    if (max(data.size) > 150) {
      plot.width <- 15
      plot.height <- 15      
    } else {
      plot.width <- 11
      plot.height <- 11      
    }
    file.ext <- get.file.parts(to.file)$ext
    if (tolower(file.ext) == '.png') {
      png(filename=to.file, width=plot.width, height=plot.height, units="in", res=600)      
    } else {
      pdf(file=to.file,width=plot.width,height=plot.height)
    }    
  }
  
  # Ajust font sizes
  #cex.val <- 1
  cex.val <- 0.9
  if (max(data.size) > 50) {
    cex.val <- 0.7
  }
  if (max(data.size) > 70) {
    cex.val <- 0.6
  }
  if (max(data.size) < 20) {
    cex.val <- 1
  }
  
  # Decide whether to show row or column names
  lab.row <- NULL
  row.sep <- c(1:nrow(clean.data))     
  lab.col <- NULL
  col.sep <- c(1:ncol(clean.data))

  if (nrow(clean.data) > 200) {
    lab.row <- NA
    row.sep <- NULL
    col.sep <- NULL
  }
  if (ncol(clean.data) > 200) {
    lab.col <- NA
    row.sep <- NULL
    col.sep <- NULL
  }
  
  # Compute clustering
#   orig.clean.data <- clean.data
#   clean.data[clean.data >= filt.thresh] <- 0.5
#   clean.data[clean.data >= 2*filt.thresh] <- 1
#   clean.data[clean.data < filt.thresh] <- 0
#   clean.data <- clean.data + (pseudo.count * matrix( data=runif(prod(clean.data)), nrow(clean.data), ncol(clean.data) ) )
  
  row.cluster.results <- T
  col.cluster.results <- T
  if (grepl(pattern="pearson|spearman",x=dist.metric)) {
    if (is.logical(row.cluster)) {
      if (row.cluster) {
        if (use.as.dist) {
          row.cluster.results <- hclust( as.dist( -clean.data ), method=clust.method )  
        } else {        
          row.cluster.results <- hclust( as.dist( 1 - cor( t(clean.data),method=dist.metric,use="na.or.complete" )^2),method=clust.method )
        }
        row.cluster <- as.dendrogram(row.cluster.results)
      }  
    }
    if (is.logical(col.cluster)) {
      if (col.cluster) {
        if (use.as.dist) {
          col.cluster.results <- hclust( as.dist( -t(clean.data)), method=clust.method )
        } else {        
          col.cluster.results <- hclust( as.dist( 1 - cor( clean.data,method=dist.metric,use="na.or.complete" )^2),method=clust.method )
        }
        col.cluster <- as.dendrogram(col.cluster.results)        
      }  
    }
  } else {
    if (is.logical(row.cluster)) {
      if (row.cluster) {        
        if (use.as.dist) {
          row.cluster.results <- hclust( as.dist( -clean.data ), method=clust.method )  
        } else {
          row.cluster.results <- hclust( dist( clean.data, method=dist.metric ), method=clust.method )  
        }        
        row.cluster <- as.dendrogram(row.cluster.results)
      }  
    }
    if (is.logical(col.cluster)) {
      if (col.cluster) {
        if (use.as.dist) {
          col.cluster.results <- hclust( as.dist( -t(clean.data)), method=clust.method )
        } else {
          col.cluster.results <- hclust( dist( t(clean.data), method=dist.metric ), method=clust.method ) 
        }        
        col.cluster <- as.dendrogram(col.cluster.results)        
      }  
    }    
  }
  
#   clean.data <- orig.clean.data
  # Check if user wants to cluster rows and columns symmetrically
  if ( symm.cluster && (nrow(clean.data) == ncol(clean.data)) ) {
    if (all(rownames(clean.data) %in% colnames(clean.data))) {
      m.idx <- match(rownames(clean.data),colnames(clean.data))
      clean.data <- clean.data[,m.idx]
      col.cluster <- row.cluster
      col.cluster.results <- row.cluster.results      
    }
  }

  # Plot heat map
  if (scale == "none") {
    heatmap.2( clean.data,
             Rowv = row.cluster, 
             Colv = col.cluster,
             dendrogram=show.dendro,
             hclustfun = function(x) hclust(x,method=clust.method),
             cexRow = cex.val,
             cexCol = cex.val,
             scale=scale,
             margins = c(9,9),
             #col = cm.colors(256),
             #col = gray( seq(0,1,length.out=(length(breaks.vals)-1)) ),
             #col = heat.colors(length(breaks.vals)-1),
             col = all.colors,
             breaks = breaks.vals,
             density.info="none",
             trace="none",
             keysize=0.8,
             colsep=col.sep,
             rowsep=row.sep,
             sepcolor="grey",
             sepwidth=c(0.01,0.01),
             na.color="white",
             xlab=col.title,
             ylab=row.title,
             labRow=lab.row,
             labCol=lab.col,   
             main=title.name,
             las = 2)
  } else {
    heatmap.2( clean.data,
             Rowv = row.cluster, 
             Colv = col.cluster,
             dendrogram=show.dendro,   
             hclustfun = function(x) hclust(x,method=clust.method),
             cexRow = cex.val,
             cexCol = cex.val,
             scale=scale,
             margins = c(9,9),
             #col = cm.colors(256),
             #col = gray( seq(0,1,length.out=(length(breaks.vals)-1)) ),
             #col = heat.colors(length(breaks.vals)-1),
             col = all.colors,             
             density.info="none",
             trace="none",
             keysize=0.8,
             colsep=col.sep,
             rowsep=row.sep,
             sepcolor="grey",
             sepwidth=c(0.01,0.01),
             na.color="white",
             xlab=col.title,
             ylab=row.title,
             labRow=lab.row,
             labCol=lab.col,   
             main=title.name,
             las = 2)
      
    }
             
  if (!is.null(to.file)) { dev.off() }
  
  invisible(list(row.cluster=row.cluster.results,
                 col.cluster=col.cluster.results))
}

plot.importance <- function(rulefit.results, output.dir=NULL, output.filename=NULL, ext="pdf", filt.thresh=5){
  # ===================================
  # Plots variable importance
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # output.filename: OPTIONAL file name (no path)
  # ext: OPTIONAL plot type (png/pdf)
  
  # Load rulefit.results if input is a data list
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  target.name <- rulefit.results$dataset$target.name
  plot.title <- sprintf("TF Importance | %s",target.name)
  
  if (!is.null(output.dir)) {
    output.dir <- file.path(output.dir,target.name) 
    # Create output directory if it doesnt exist
    if (!file.exists(output.dir)){
      dir.create(output.dir,recursive=T)
    }
    if (is.null(output.filename)) {
      output.filename <- file.path( output.dir , sprintf("factor.importance.%s.%s",target.name,ext) )
    } else {
      output.filename <- file.path( output.dir , get.file.parts(output.filename)$fullname )
    }
  }
  
  if (! is.null(filt.thresh) ) {
    rulefit.results$vi <- rulefit.results$vi[, (rulefit.results$vi >= filt.thresh) ]
  }
  rulefit.results$vi <- filter.cols(rulefit.results$vi)
  colnames(rulefit.results$vi) <- standardize.name(colnames(rulefit.results$vi))
  
  make.barplot( as.numeric(rulefit.results$vi), labels=colnames(rulefit.results$vi), title.name=plot.title , to.file=output.filename )
}

plot.int.strength <- function(rulefit.results, output.dir=NULL, output.filename=NULL, ext="png", filt.thresh=1e-7){
  # ===================================
  # Plots interaction strength
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # output.filename: OPTIONAL file name (no path)
  # ext: OPTIONAL plot type (png/pdf)
  
  # Load rulefit.results if input is a data list
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  target.name <- rulefit.results$dataset$target.name
  plot.title <- sprintf("Interaction Strength | %s",target.name)
  
  if ( !is.null(output.dir) ) {
    output.dir <- file.path(output.dir,target.name) 
    # Create output directory if it doesnt exist
    if (!file.exists(output.dir)){
      dir.create(output.dir,recursive=T)
    }  
  
    if (is.null(output.filename)) {
      output.filename <- file.path( output.dir , sprintf("int.strength.%s.%s",target.name,ext) )
    } else {
      output.filename <- file.path( output.dir , get.file.parts(output.filename)$fullname )
    }
  }

  if (! is.null(filt.thresh) ) {
    rulefit.results$int.strength <- rulefit.results$int.strength[, (rulefit.results$int.strength >= filt.thresh) ]
  }
  
  make.barplot( as.numeric(rulefit.results$int.strength), labels=toupper(colnames(rulefit.results$int.strength)), title.name=plot.title , to.file=output.filename )
}

plot.pairwise <- function(rulefit.results, output.dir=NULL, ext="png", filt.thresh=1e-7){
  # ===================================
  # Plots Pairwise interaction strength for each partner TFs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # ext: OPTIONAL output file type (png/pdf)
  # filt.thresh: only consider interaction strengths > filt.thresh  
  
  # Load rulefit.results if input is a data list
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$dataset$target.name
  output.dir <- file.path(output.dir,target.name) 
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir, recursive=T)
  }
  
  n.tfs <- dim(rulefit.results$pair.interactions)[[1]] # number of TFs
  
  for (i in c(1:n.tfs)) {
    curr.vector <- rulefit.results$pair.interactions[i,]
    curr.vector <- curr.vector[,-i]
    curr.tf <- rownames(rulefit.results$pair.interactions)[[i]]
    curr.vector <- curr.vector[ , (!is.na(curr.vector)) ]
    topN <- sum(curr.vector>filt.thresh) # Get number of variables with significant interaction scores
    if (topN == 0) {next}
    cat("Plotting TF ",curr.tf," given ",target.name,"\n")
    plot.title <- sprintf("Pair interactions of %s | %s", curr.tf, target.name)
    output.filename <- file.path( output.dir , sprintf("pair.interact.%s.given.%s.%s", curr.tf, target.name, ext) )
    make.barplot( vals=as.numeric(curr.vector), labels=toupper(colnames(curr.vector)), title.name=plot.title , to.file=output.filename, top.N=topN )    
  }
}

plot.pairwise.matrix <- function(rulefit.results, output.dir, output.filename=NA, ext="pdf", filter.thresh=1e-7, use.null=F){
  # ===================================
  # Plots interaction strength
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # output.filename: OPTIONAL file name (no path)
  # ext: OPTIONAL plot type (png/pdf)
  # filter.thresh: Threshold used to filter and normalize scores
  # use.null: if set to T, then if null scores are available they will be subtracted from the true scores

  # Load rulefit.results if input is a character vector
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$dataset$target.name
  output.dir <- file.path(output.dir,target.name) 
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir, recursive=T)
  }

  plot.title <- sprintf("Conditional Pairwise Interactions | %s",target.name)

  if (is.na(output.filename)) {
    output.filename <- file.path( output.dir , sprintf("cond.pairwise.int.matrix.%s.%s",target.name,ext) )
  } else {
    output.filename <- file.path( output.dir , get.file.parts(output.filename)$fullname )
  }
  
  val.data <- rulefit.results$pair.interactions
  if (use.null) {
    if (any(names(rulefit.results)=="pair.interactions.null.mean")) {
      rulefit.results$pair.interactions.null.mean[is.na(rulefit.results$pair.interactions.null.mean)] <- 0
      val.data <- val.data - rulefit.results$pair.interactions.null.mean
      val.data[val.data < 0] <- 0
    }
  }
  val.data <- filter.cols( filter.rows(val.data) )
  rownames(val.data) <- standardize.name(rownames(val.data))
  colnames(val.data) <- standardize.name(colnames(val.data))
  
#   plot.heatmap( data=val.data,
#                 show.dendro="none",
#                 to.file=output.filename,
#                 row.title="Transcription Factors",
#                 col.title="Transcription Factors",
#                 title.name="",
#                 filt.thresh=filter.thresh,
#                 pseudo.count=1e-30,
#                 logval=F,
#                 replace.diag=T,
#                 replace.na=T,
#                 num.breaks=255,
#                 #clust.method="ward",
#                 clust.method="single",
#                 break.lowerbound=1e-3,
#                 break.type="quantile")
  
  plot.heatmap( data=val.data,
                use.as.dist=T,
                show.dendro="none",
                to.file=output.filename,
                row.title="Transcription Factors",
                col.title="Transcription Factors",
                title.name="",
                filt.thresh=filter.thresh,
                pseudo.count=0,
                logval=T,
                replace.diag=T,
                replace.na=T,
                num.breaks=255,
                #clust.method="complete",
                clust.method="ward",
                break.lowerbound=4.5,                  
                break.type="linear")  
}

plot.singleplot <- function(rulefit.results, output.dir=NULL, ext="png", filt.thresh=1e-7){
  # ===================================
  # Plots Pairwise interaction strength for each partner TFs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # ext: OPTIONAL output file type (png/pdf)
  # filt.thresh: only consider interaction strengths > filt.thresh  
  
  # Load rulefit.results if input is a data list
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$dataset$target.name
  output.dir <- file.path(output.dir,target.name) 
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir, recursive=T)
  }
  
  n.tfs <- dim(rulefit.results$pair.interactions)[[1]] # number of TFs
  
  for (i in c(1:n.tfs)) {
#     curr.vector <- rulefit.results$pair.interactions[i,]
#     curr.vector <- curr.vector[,-i]
#     curr.tf <- rownames(rulefit.results$pair.interactions)[[i]]
#     curr.vector <- curr.vector[ (!is.na(curr.vector)) , ]
#     topN <- sum(curr.vector>filt.thresh) # Get number of variables with significant interaction scores
#     if (topN == 0) {next}
#     cat("Plotting TF ",curr.tf," given ",target.name,"\n")
#     plot.title <- sprintf("Pair interactions of %s | %s", curr.tf, target.name)
#     output.filename <- file.path( output.dir , sprintf("pair.interact.%s.given.%s.%s", curr.tf, target.name, ext) )
#     make.barplot( as.numeric(curr.vector), labels=toupper(colnames(curr.vector)), title.name=plot.title , to.file=output.filename, top.N=top.N )    
  }
}

# ###########################################################################
# ###########################################################################
# THE FOLLOWING FUNCTIONS OPERATE ONLY AFTER aggregrate.model.randneg.R 
# is run on the multiple random negative set runs
# ###########################################################################
# ###########################################################################
get.average.vi <- function(rulefit.results) {
  # ===================================
  # Get aveage statistics of variable importance from several random runs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # RETURNS rulefit.results with an extra field mean.vi (data.frame) which has
  #     $mean.val : median values
  #     $std.val : std deviations
  #     $lqr: lower quartile
  #     $hqr: upper quartile
  #     $val.names: names of partners
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # Load rulefit.results if input is a data list
  
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name

  # Get mean and std of vi
  val.names <- colnames(rulefit.results$vi)
  mean.val <- apply(rulefit.results$vi, 2, function(x) median(x,na.rm=T)) # get mean of each column
  std.val <- apply(rulefit.results$vi, 2, function(x) sd(x,na.rm=T)) # get std of each column
  lqr <- apply(rulefit.results$vi, 2, function(x) quantile(x,0.25,na.rm=T)) # get lqr of each column
  hqr <- apply(rulefit.results$vi, 2, function(x) quantile(x,0.75,na.rm=T)) # get hqr of each column
  val.data <- data.frame( mean.val=as.vector(mean.val), std.val=as.vector(std.val), lqr=as.vector(lqr), hqr=as.vector(hqr), tf.name = val.names ) # Create data frame
  
  rulefit.results$mean.vi <- val.data
  return(rulefit.results)  
}

get.average.cv <- function(rulefit.results) {
  # ===================================
  # Plots variable importance from several random runs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # RETURNS rulefit.results with an extra field mean.cv (data.frame) which has
  #     $mean.val : median values
  #     $std.val : std deviations
  #     $lqr: lower quartile
  #     $hqr: upper quartile
  #     $cv.names: names of various cross validation parameters
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # Load rulefit.results if input is a data list
  
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  rulefit.results$cv$lo <- NULL
  # Get mean, std, lqr and hqr of all cross-validation metrics
  cv.names <- colnames(rulefit.results$cv)
  mean.val <- apply(rulefit.results$cv, 2, function(x) median(x,na.rm=T)) # get mean of each column
  std.val <- apply(rulefit.results$cv, 2, function(x) sd(x,na.rm=T)) # get std of each column
  lqr <- apply(rulefit.results$cv, 2, function(x) quantile(x,0.25,na.rm=T)) # get lqr of each column
  hqr <- apply(rulefit.results$cv, 2, function(x) quantile(x,0.75,na.rm=T)) # get hqr of each column
  val.data <- data.frame( mean.val=as.vector(mean.val), std.val=as.vector(std.val), lqr=as.vector(lqr), hqr=as.vector(hqr), cv.names = cv.names ) # Create data frame
  rownames(val.data) <- cv.names
  
  rulefit.results$mean.cv <- val.data
  return(rulefit.results)  
}

plot.average.vi <- function(rulefit.results, output.dir, output.filename=NA, ext="png", thresh=1) {
  # ===================================
  # Plots variable importance from several random runs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing aggregated rulefit.results list OR the aggregated rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # output.filename: OPTIONAL file name (no path)
  # ext: OPTIONAL plot type (png/pdf)
  # thresh: OPTIONAL threshold to use to filter variables whose mean importance is < thresh
  # Load rulefit.results if input is a data list
  
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  output.dir <- file.path(output.dir,target.name) 
  
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir,recursive=T)
  }
  plot.title <- sprintf("TF Importance | %s",target.name)

  if (is.na(output.filename)) {
    output.filename <- file.path( output.dir , sprintf("factor.importance.%s.%s",target.name,ext) )
  } else {
    output.filename <- file.path( output.dir , get.file.parts(output.filename)$fullname )
  }

  library(ggplot2)
  axes.format <- opts(plot.title = theme_text(size=12,vjust=1),                    
                      axis.text.x = theme_text(size=10,colour="black"),
                      axis.text.y = theme_text(size=10,colour="black",hjust=1),
                      axis.title.x = theme_text(size=12),
                      axis.title.y = theme_text(size=12,angle=90),
                      legend.title = theme_text(size=10,hjust=0),
                      legend.text = theme_text(size=10)                      
                      )
  
  # Get mean and std of vi
  val.data <-rulefit.results$mean.vi
  val.data <- val.data[val.data$mean.val >= thresh, ]
  if (length(val.data$tf.name) == 0) {return()}
  
  # Remove unwanted tfs
  rownames(val.data) <- val.data$tf.name
  val.data <- filter.rows(val.data)
  #val.data$tf.name <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "",toupper( gsub("K562b|Hepg2b", "B-", val.data$tf.name) ), ignore.case=T)
  val.data$tf.name <- standardize.name(val.data$tf.name)
  
  p1 <- ggplot(val.data) + 
    #geom_bar( aes( x=reorder(tf.name,mean.val) , y=mean.val, fill=mean.val) ) + 
    geom_bar( aes( x=reorder(tf.name,mean.val) , y=mean.val), fill="red3", alpha=0.8 ) +
    geom_errorbar( aes( x=reorder(tf.name,mean.val), ymax=hqr, ymin=lqr) )
  axes.labels <- labs(x = "TF", y = "Relative variable importance") # axes labels
  #axes.labels <- labs(x = "", y = "") # axes labels
  p1 <- p1 + 
    axes.labels + 
    axes.format + 
    #scale_fill_gradient("VarImp") + 
    opts(title=plot.title) + 
    coord_flip()

  if (nrow(val.data) > 50) {
    p1 <- p1 + opts(axis.text.y = theme_text(size=7,colour="black",hjust=1))
  }

  if (tolower(ext) == "png") {
    ggsave(file=output.filename, plot=p1, width=6, height=10, dpi=600)  
  } else {
    ggsave(file=output.filename, plot=p1, width=6, height=10)  
  }
}

get.average.int.strength <- function(rulefit.results) {
  # ===================================
  # Plots average interaction strength from several aggregated random runs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # RETURNS rulefit.results with an extra field mean.int.strength (data.frame) which has
  #     $mean.val : median values
  #     $std.val : std deviations
  #     $lqr: lower quartile
  #     $hqr: upper quartile
  #     $val.names: names of partners  
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # Load rulefit.results if input is a data list
  
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  
  # Get mean and std of vi
  val.names <- colnames(rulefit.results$int.strength)
  mean.val <- apply(rulefit.results$int.strength, 2, function(x) median(x,na.rm=T)) # get mean of each column
  std.val <- apply(rulefit.results$int.strength, 2, function(x) sd(x,na.rm=T)) # get std of each column
  lqr <- apply(rulefit.results$int.strength, 2, function(x) quantile(x,0.25,na.rm=T)) # get lqr of each column
  hqr <- apply(rulefit.results$int.strength, 2, function(x) quantile(x,0.75,na.rm=T)) # get lqr of each column
  val.data <- data.frame( mean.val=as.vector(mean.val), std.val=as.vector(std.val), lqr=as.vector(lqr), hqr=as.vector(hqr), tf.name = val.names )
  rulefit.results$mean.int.strength <- val.data
  return(rulefit.results)
}

plot.average.int.strength <- function(rulefit.results, output.dir, output.filename=NA, ext="png", thresh=1e-7) {
  # ===================================
  # Plots average interaction strength from several random runs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # output.filename: OPTIONAL file name (no path)
  # ext: OPTIONAL plot type (png/pdf)
  # thresh: OPTIONAL used to filter variables whose mean int.strength is < thresh
  # Load rulefit.results if input is a data list
  
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  output.dir <- file.path(output.dir,target.name) 
  
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir,recursive=T)
  }
  
  plot.title <- sprintf("Interaction Strength | %s",target.name)

  if (is.na(output.filename)) {
    output.filename <- file.path( output.dir , sprintf("int.strength.%s.%s",target.name,ext) )
  } else {
    output.filename <- file.path( output.dir , get.file.parts(output.filename)$fullname )
  }
  
  library(ggplot2)
  axes.format <- opts(plot.title = theme_text(size=12,vjust=1),                    
                      axis.text.x = theme_text(size=10,colour="grey30"),
                      axis.text.y = theme_text(size=10,colour="grey30",hjust=1),
                      axis.title.x = theme_text(size=12),
                      axis.title.y = theme_text(size=12,angle=90),
                      legend.title = theme_text(size=10,hjust=0),
                      legend.text = theme_text(size=10)                      
                      )
  
  # Get mean and std of vi
  val.data <- rulefit.results$mean.int.strength
  val.data <- val.data[val.data$mean.val >= thresh, ]
  if (length(val.data$tf.name) == 0) {return()}
    
  # Remove unwanted tfs
  rownames(val.data) <- val.data$tf.name
  val.data <- filter.rows(val.data)
  #val.data$tf.name <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "",toupper( gsub("K562b|Hepg2b", "B-", val.data$tf.name) ), ignore.case=T)
  val.data$tf.name <- standardize.name(val.data$tf.name)    
  
  p1 <- ggplot(val.data) + 
    geom_bar( aes( x=reorder(tf.name,mean.val) , y=mean.val, fill=mean.val) ) + 
    geom_errorbar( aes( x=reorder(tf.name,mean.val), ymax=hqr, ymin=lqr) )                      
  axes.labels <- labs(x = "TF", y = "Interaction Strength")
  p1 <- p1 + axes.labels + axes.format + scale_fill_gradient("IntStrength") + opts(title=plot.title) + coord_flip()
                      
  if (nrow(val.data) > 50) {
    p1 <- p1 + opts(axis.text.y = theme_text(size=7,colour="grey30",hjust=1))
  }

  if (tolower(ext) == "png") {
    ggsave(file=output.filename, plot=p1, width=6, height=10, dpi=600)  
  } else {
    ggsave(file=output.filename, plot=p1, width=6, height=10)  
  }
}
    
get.average.pairwise <- function(rulefit.results) {
  # ===================================
  # Plots Pairwise interaction strength for each partner TF from several aggregated random runs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # RETURNS rulefit.results with an two extra fields 
  # aggregate.pairwise.interactions (a LIST containing data.frames for each TF) which has
  #     $mean.val : median values
  #     $std.val : std deviations
  #     $lqr: lower quartile
  #     $hqr: upper quartile
  #     $val.names: names of partners
  # mean.pairwise.int.matrix (data.frame) TFS X TFS  
  # ===================================  
  # rulefit.results: Rdata file name containing rulefit.results list OR the rulefit.results LIST
  # Load rulefit.results if input is a data list
  
  # Load rulefit.results if input is a data list
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  partner.names <- names(rulefit.results$pair.interactions) # names of each LHS partner
  n.tfs <- length(partner.names) # number of TFs
  
  mean.pairwise.int.matrix <- as.data.frame(matrix(NA, nrow=n.tfs, ncol=n.tfs)) # Place holder for mean pairwise interaction matrix
  rownames(mean.pairwise.int.matrix) <- partner.names
  colnames(mean.pairwise.int.matrix) <- partner.names
  
  aggregate.pairwise.interactions <- list() # each element of the list has statistics for conditional pairwise interactions of each partner TF
  
  for (i in c(1:n.tfs)) {
    lhs.tf.name <- partner.names[[i]] # current partner    
    curr.matrix <- rulefit.results$pair.interactions[[i]]
    val.names <- colnames(curr.matrix)
    mean.val <- apply(curr.matrix, 2, function(x) median(x,na.rm=T)) # get mean of each column    
    std.val <- apply(curr.matrix, 2, function(x) sd(x,na.rm=T)) # get std of each column
    lqr <- apply(curr.matrix, 2, function(x) quantile(x,0.25,na.rm=T)) # get lqr of each column
    hqr <- apply(curr.matrix, 2, function(x) quantile(x,0.75,na.rm=T)) # get hqr of each column
    val.data <- data.frame( mean.val=as.vector(mean.val), std.val=as.vector(std.val), lqr=as.vector(lqr), hqr=as.vector(hqr), tf.name = val.names )
    mean.pairwise.int.matrix[lhs.tf.name, names(mean.val)] <- mean.val
    aggregate.pairwise.interactions[[lhs.tf.name]] <- val.data
  }
  
  rulefit.results$mean.pairwise.int.matrix <- mean.pairwise.int.matrix
  rulefit.results$aggregate.pairwise.interactions <- aggregate.pairwise.interactions
  return(rulefit.results)
}
    
plot.average.pairwise <- function(rulefit.results, output.dir, ext="png", filter.thresh=1e-7) {
  # ===================================
  # Plots average pairwise interactions from several random runs
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing aggregated rulefit.results list OR the aggregated rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # ext: OPTIONAL plot type (png/pdf)
  # Load rulefit.results if input is a data list
  
  # Load rulefit.results if input is a data list
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  output.dir <- file.path(output.dir,target.name) 
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir, recursive=T)
  }
  
  library(ggplot2)
  axes.format <- opts(plot.title = theme_text(size=12,vjust=1),                    
                      axis.text.x = theme_text(size=10,colour="grey30"),
                      axis.text.y = theme_text(size=10,colour="grey30",hjust=1),
                      axis.title.x = theme_text(size=12),
                      axis.title.y = theme_text(size=12,angle=90),
                      legend.title = theme_text(size=10,hjust=0),
                      legend.text = theme_text(size=10)                      
                      )
  
  partner.names <- names(rulefit.results$aggregate.pairwise.interactions) # names of each LHS partner
  n.tfs <- length(partner.names) # number of TFs
  
  for (i in c(1:n.tfs)) {
    lhs.tf.name <- partner.names[[i]] # current partner
    plot.title <- sprintf("Pair interactions of %s | %s", lhs.tf.name, target.name)
    cat("Plotting TF ", lhs.tf.name, " given ", target.name, "\n")
    output.filename <- file.path( output.dir , sprintf("pair.interact.%s.given.%s.%s", lhs.tf.name, target.name, ext) )
    
    val.data <- rulefit.results$aggregate.pairwise.interactions[[i]]
    rownames(val.data) <- val.data$tf.name
    val.data <- filter.rows(val.data)
    #val.data$tf.name <- gsub("GM12878|K562|HelaS3|Hepg2|H1hesc", "",toupper( gsub("K562b|Hepg2b", "B-", val.data$tf.name) ), ignore.case=T)
    val.data$tf.name <- standardize.name(val.data$tf.name)                          

    if ( max(val.data$mean.val,na.rm=T) < filter.thresh) { # skip if max value is < threshold
      next
    }    
    val.data <- droplevels( val.data[ !is.na(val.data$mean.val) , ] )
    val.data <- droplevels( val.data[ (val.data$mean.val >= filter.thresh), ] )
    
    if (length(val.data) == 0) {next}
    
    p1 <- ggplot(val.data) + 
      geom_bar( aes( x=reorder(tf.name,mean.val) , y=mean.val, fill=mean.val) ) + 
      geom_errorbar( aes( x=reorder(tf.name,mean.val), ymax=hqr, ymin=lqr ) )                  
    axes.labels <- labs(x = "TF", y = "Pairwise Interaction Strength")
    p1 <- p1 + axes.labels + axes.format + scale_fill_gradient("IntStrength") + opts(title=plot.title) + coord_flip()
                      
    if (nrow(val.data) > 50) {
      p1 <- p1 + opts(axis.text.y = theme_text(size=7,colour="grey30",hjust=1))
    }
    if (nrow(val.data) < 10) {
      p1 <- p1 + opts(axis.text.y = theme_text(size=12,colour="grey30",hjust=1))
    }

    if (tolower(ext) == "png") {
      ggsave(file=output.filename, plot=p1, width=6, height=10, dpi=600)  
    } else {
      ggsave(file=output.filename, plot=p1, width=6, height=10)  
    }    
  }  
}
  
plot.average.pairwise.matrix <- function(rulefit.results, output.dir, output.filename=NA, ext="pdf", filter.thresh=1e-7) {
  # ===================================
  # Plots heatmap of average/aggregated conditional pairwise interactions
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing aggregated rulefit.results list OR the aggregated rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # output.filename: OPTIONAL file name (no path)
  # ext: OPTIONAL plot type (png/pdf)
  # filter.thresh: threshold to filter interactions
  # Load rulefit.results if input is a data list
  
  # Load rulefit.results if input is a character vector
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  output.dir <- file.path(output.dir,target.name) 
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir, recursive=T)
  }

  plot.title <- sprintf("Conditional Pairwise Interactions | %s",target.name)

  if (is.na(output.filename)) {
    output.filename <- file.path( output.dir , sprintf("cond.pairwise.int.matrix.%s.%s",target.name,ext) )
  } else {
    output.filename <- file.path( output.dir , get.file.parts(output.filename)$fullname )
  }
  
  val.data <- rulefit.results$mean.pairwise.int.matrix
  val.data <- filter.cols( filter.rows(val.data) )
  rownames(val.data) <- standardize.name(rownames(val.data))
  colnames(val.data) <- standardize.name(colnames(val.data))
  
#   plot.heatmap( data=val.data,
#                 show.dendro="none",
#                 to.file=output.filename,
#                 row.title="Transcription Factors",
#                 col.title="Transcription Factors",
#                 title.name="",
#                 filt.thresh=filter.thresh,
#                 pseudo.count=1e-30,
#                 logval=F,
#                 replace.diag=T,
#                 replace.na=T,
#                 num.breaks=255,
#                 #clust.method="ward",
#                 clust.method="single",
#                 break.lowerbound=1e-3,
#                 break.type="quantile")
  
  plot.heatmap( data=val.data,
                use.as.dist=T,
                show.dendro="none",
                to.file=output.filename,
                row.title="Transcription Factors",
                col.title="Transcription Factors",
                title.name="",
                filt.thresh=filter.thresh,
                pseudo.count=0,
                logval=T,
                replace.diag=T,
                replace.na=T,
                num.breaks=255,
                #clust.method="complete",
                clust.method="ward",
                break.lowerbound=4.5,                  
                break.type="linear")  
}

plot.score.dist.average.pairwise.matrix <- function(rulefit.results, output.dir, output.filename=NA, ext="pdf", filter.thresh=1e-7, feature.type="split") {
  # ===================================
  # Plots heatmap of average/aggregated conditional pairwise interactions
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing aggregated rulefit.results list OR the aggregated rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # output.filename: OPTIONAL file name (no path)
  # ext: OPTIONAL plot type (png/pdf)
  # filter.thresh: Threshold used
  # feature.type: all/null/main/split all: use all features, null: use null features, main: use true features, split: plot main and null on the same plot
  # Load rulefit.results if input is a data list
  
  # Load rulefit.results if input is a character vector
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  output.dir <- file.path(output.dir,target.name) 
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir, recursive=T)
  }

  plot.title <- sprintf("Pairwise Interactions Score distributions| %s",target.name)

  if (is.na(output.filename)) {
    output.filename <- file.path( output.dir , sprintf("score.dist.pairwise.int.matrix.%s.%s.%s",feature.type,target.name,ext) )
  } else {
    output.filename <- file.path( output.dir , get.file.parts(output.filename)$fullname )
  }
  
  val.data <- rulefit.results$mean.pairwise.int.matrix
  val.data <- filter.cols( filter.rows(val.data) )
  
  # Split into main and null
  main.row.idx <- grep(".*random.*", rownames(val.data), invert=T)
  null.row.idx <- grep(".*random.*", rownames(val.data))
  main.col.idx <- grep(".*random.*", colnames(val.data), invert=T)
  null.col.idx <- grep(".*random.*", colnames(val.data))
  
  # Replace 0s and NAs with smallest non-zero score
  all.vals <- as.vector( as.matrix(val.data) )
  all.vals[which(all.vals==0)] <- NA
  min.score <- min(all.vals,na.rm=T)
  all.vals[is.na(all.vals)] <- min.score
  
  # Replace 0s and NAs with smallest non-zero score
  main.vals <- as.vector( as.matrix( val.data[main.row.idx,main.col.idx] ) )
  main.vals[which(main.vals==0)] <- NA
  main.vals[is.na(main.vals)] <- min.score
  main.vals <- log10(main.vals)
  
  # Get null vals as the set of all scores corresponding to any interaction involving a null feature
  null.vals <- vector()
  null.exist <- F
  if (length(null.row.idx) > 0) {
    null.exist <- T
    null.vals <- as.vector( as.matrix( val.data[null.row.idx,main.col.idx] ) )
  }
  if (length(null.col.idx) > 0) {
    null.exist <- T
    null.vals <- c(null.vals, as.vector( as.matrix( val.data[main.row.idx,null.col.idx] ) ) )
  }
  if ( (length(null.row.idx) > 0) & (length(null.col.idx) > 0) ) {
    null.exist <- T
    null.vals <- c( null.vals, as.vector( as.matrix( val.data[null.row.idx,null.col.idx] ) ) )
  }
  if (null.exist) {
    # Replace 0s and NAs with smallest non-zero score
    null.vals[which(null.vals==0)] <- NA
    null.vals[is.na(null.vals)] <- min.score
    null.vals <- log10(null.vals)
  }
  
  # Convert to data frames for plotting using ggplot
  main.vals.data.frame <- data.frame(type="true", log.scores=main.vals - log10(filter.thresh))
  all.vals.data.frame <- main.vals.data.frame
  if (null.exist) {
    null.vals.data.frame <- data.frame(type="null", log.scores=null.vals - log10(filter.thresh))
    all.vals.data.frame <- rbind(all.vals.data.frame, null.vals.data.frame)    
  }
  
  # Plot figure
  require(ggplot2)
  if ( any(feature.type == c("all","split")) ) {
    p <- ggplot(all.vals.data.frame) 
  } else if (feature.type == "main") {
    p <- ggplot(main.vals.data.frame)
  } else {
    p <- ggplot(null.vals.data.frame)
  }
  
  if (feature.type == "all") {
    p <- p + geom_density(aes(x=log.scores)) + opts(title=target.name)
  } else {
    p <- p + geom_density(aes(x=log.scores,color=type,fill=type),alpha=I(0.3), size=I(0.2)) + opts(title=target.name)
  }
  
  # Obtain thresholds based on quantiles of interaction scores of null features
  thresholds <- NA
  if (null.exist) {
    thresholds <- quantile( null.vals, 1-c(0.05,0.01,0.005,0.001,0.0005,0.0001) ) - log10(filter.thresh)
  }
  
  ggsave(filename=output.filename, plot=p, dpi=300, width=5,height=5)
  
  invisible(list(log.main.vals=main.vals,
                 log.null.vals=null.vals,
                 log.all.vals.df=all.vals.data.frame,
                 thresholds=thresholds))
}

write.table.average.pairwise.matrix <- function(rulefit.results, output.dir, filter.thresh=1e-7) {
  # ===================================
  # Write average/aggregated conditional pairwise interaction matrix to a tab delimited file
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing aggregated rulefit.results list OR the aggregated rulefit.results LIST
  # output.dir: directory where you want to save all figure
  # output.filename: OPTIONAL file name (no path)
  # filter.thresh: threshold to filter interactions
  # Load rulefit.results if input is a data list
  
  # Load rulefit.results if input is a character vector
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  output.dir <- file.path(output.dir,target.name) 
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir, recursive=T)
  }

  # Write pairwise interactions as matrix
  output.filename <- file.path( output.dir , sprintf("cond.pairwise.int.matrix.%s.xls",target.name) )
  val.data <- as.matrix(rulefit.results$mean.pairwise.int.matrix)
  write.table(val.data, file=output.filename, quote=F, sep="\t", col.names=NA)
  
  # Write pairwise interactions as adjacency list
  output.filename <- file.path( output.dir , sprintf("cond.pairwise.int.stats.%s.xls",target.name) )
  partner.names <- names(rulefit.results$aggregate.pairwise.interactions) # names of all partners
  count = 1
  for (p1 in partner.names) {
    curr.pairwise <- rulefit.results$aggregate.pairwise.interactions[[p1]]
    if (nrow(curr.pairwise)==0) {next}
    curr.pairwise$target.tf.context <-  target.name # Add new column with target name
    curr.pairwise$tf1 <-  p1 # Add new column with tf1
    rearranged.colnames <- c("target.tf.context", "tf1","tf.name", "mean.val", "std.val", "lqr", "hqr")
    curr.pairwise <- curr.pairwise[, rearranged.colnames]
    colnames(curr.pairwise) <- c("target.tf.context", "tf1", "tf2", "median.pairwise.int.strength", "std", "lower.quartile", "upper.quartile")
    if (count==1) {
      write.table(curr.pairwise, file=output.filename, quote=F, sep="\t", col.names=T, row.names=F, na="")
    } else {
      write.table(curr.pairwise, file=output.filename, quote=F, sep="\t", col.names=F, append=T, row.names=F, na="")
    }
    count=count+1
  }
}

write.symmetric.average.pairwise.matrix <- function(rulefit.results, output.dir, sparsify=T, filt.thresh=1e-7) {
  # ===================================
  # Write average/aggregated conditional pairwise interaction matrix to a tab delimited file
  # Takes as input rulefit.results (list) OR
  # Rdata file name (string) that contains rulefit.results 
  # ===================================  
  # rulefit.results: Rdata file name containing aggregated rulefit.results list OR the aggregated rulefit.results LIST
  # output.dir: directory where you want to save all figure

  # Load rulefit.results if input is a character vector
  if (is.character(rulefit.results)) {
    load(rulefit.results)    
  }
  
  target.name <- rulefit.results$target.name
  #output.dir <- file.path(output.dir,target.name) 
  # Create output directory if it doesnt exist
  if (!file.exists(output.dir)){
    dir.create(output.dir, recursive=T)
  }

  # Write pairwise interactions as matrix
  output.filename <- file.path( output.dir , sprintf("cond.pairwise.int.symm.log.matrix.%s.tab",target.name) )
  val.data <- rulefit.results$mean.pairwise.int.matrix
  # Remove cols with very small values
  na.idx <- apply( val.data , 2 , function(x) all((x<filt.thresh),na.rm=T) )
  val.data <- val.data[ , !na.idx]
  # Remove rows with all very small values
  na.idx <- apply( val.data , 1 , function(x) all((x<filt.thresh),na.rm=T) )
  val.data <- val.data[!na.idx, ]
  # Remove blacklisted TFs
  val.data <- filter.rows(filter.cols(val.data)) 
  # Make symmetric
  val.data <- (val.data + t(val.data))/2
  # Remove lower triangle
  val.data[lower.tri(val.data)] <- NA
  # Standardize names
  rownames(val.data) <- standardize.name(rownames(val.data))
  colnames(val.data) <- standardize.name(colnames(val.data))
  # Convert to log10
  val.data <- log10(as.matrix(val.data))
  val.data[is.infinite(val.data)] <- NA
  if (sparsify) {
    val.data <- val.data - log10(filt.thresh)
  } else {
    val.data <- val.data - min(val.data,na.rm=T)  
  }  
  val.data[which(val.data <= 0)] <- NA
  #val.data[val.data < 0] <- 0
  # Remove rows and cols with all NAs
#   na.idx <- apply( val.data , 2 , function(x) all(is.na(x)) )
#   val.data <- val.data[ , !na.idx]
  
  # Remove rows with all NAs
  na.idx <- apply( val.data , 1 , function(x) all(is.na(x)) )
  val.data <- val.data[!na.idx, ]
  val.data <- val.data[ , !na.idx]
  write.table(val.data, file=output.filename, quote=F, sep="\t", col.names=NA, na="-")  
}

compute.proximal.distal.diff.importance <- function(input.dir, 
                                                    peak.distance.file, 
                                                    output.filename=NULL,
                                                    proximal.cutoff=5000,                                                    
                                                    distal.cutoff=10000,
                                                    rev.peak.ids=T) {
  # Computes differential relative importance for all factors comparing proximal vs. distal sites of the target factor
  # input.dir: directory containing multiple rulefit results from randomized negative sets
  # peak.distance.file: Contains peak distance to nearest TSS information
  # output.filename: output figure file name
  # proximal.cutoff: distance threshold to use as "proximal" class
  # distal.cutoff: distance threshold to use for "distal" class
  # rev.peak.ids: Reverse peak ids in the distance file (This is to be used when peak labeled rank1 in distance file refers to weakest peak, whereas peak labeled rank1 in association data is the strongest peak)
  
  # Check that input directory exists
  if (! file.exists(input.dir)) {
    stop("Input Directory ", input.dir," does not exist\n")
  }
  
  # Check that peak.distance.file exists
  if (! file.exists(peak.distance.file)) {
    stop("Peak2TSS distance file ", peak.distance.file," does not exist\n")
  }
  
  # Autogenerate outputfile name
  if (is.null(output.filename)) {
    output.stub <- gsub(pattern="\\.+$",replacement="",x=get.file.parts(peak.distance.file)$name)
    output.filename <- file.path(input.dir,
                             sprintf("%s.prox.%d.dist.%d.png", output.stub, proximal.cutoff, distal.cutoff))
  }
    
  # Read and parse distance file, get proximal and distal peak ids
  distance.table <- read.table(file=peak.distance.file,
                               header=T,
                               sep="\t",
                               #col.names=c("peak.chr","peak.start","peak.stop","peak.id","tss.chr","tss.start","tss.stop","tss.id","dist"),
                               stringsAsFactors=F)
    
  proximal.peak.ids <- distance.table$peak.id[distance.table$dist <= proximal.cutoff]
  distal.peak.ids <- distance.table$peak.id[distance.table$dist > distal.cutoff]  
    
  # Get list of Rdata files in directory
  all.Rdata.files <- list.files(path=input.dir, pattern=".*Rdata$", full.names=T) # Get names of Rdata files
  n.Files <- length(all.Rdata.files)
  if (n.Files == 0) {
    stop("No Rdata files found in input directory", input.dir, "\n")
  }

  # Reverse peak ids if necessary  
  if (rev.peak.ids) {
    load(all.Rdata.files[1])
    assoc.data.peak.ids <- rownames(rulefit.results$dataset$x.vals)[rulefit.results$dataset$x.vals == 1]
    assoc.data.peak.ranks <- as.numeric( gsub( '(^.*Pk_)|(_[^_]+$)' , '' , assoc.data.peak.ids ) ) # Convert PeakIds to numbers
    assoc.data.peak.ids <- assoc.data.peak.ids[order(assoc.data.peak.ranks)] # Order the assoc.data peak ids
    dtable.peak.ids <- rev( assoc.data.peak.ids ) # dtable.ids is reverse of assoc.ids (each index matches up)
    proximal.peak.ids <- assoc.data.peak.ids[dtable.peak.ids %in% proximal.peak.ids] # Find dtable indices that match proximal.peak.ids and then translate to assoc.ids
    distal.peak.ids <- assoc.data.peak.ids[dtable.peak.ids %in% distal.peak.ids] # Find dtable indices that match distal.peak.ids and then translate to assoc.ids
  }  
  
  # Load each Rdata file, compute differential importance and store them
  proximal.vi <- data.frame()
  distal.vi <- data.frame()
  for (each.file in all.Rdata.files) {
    rulefit.results <- restore.rf.model(each.file) # load Rdata.file
    proximal.rulefit.results <- get.var.imp(rulefit.results,class=proximal.peak.ids)
    proximal.vi <- rbind(proximal.vi, proximal.rulefit.results$vi)
    distal.rulefit.results <- get.var.imp(rulefit.results,class=distal.peak.ids)
    distal.vi <- rbind(distal.vi, distal.rulefit.results$vi)
  }
  
  diff.vi <- distal.vi - proximal.vi
  median.diff.vi <- apply(diff.vi,2,function(x) median(x,na.rm=T))  
  lqr.diff.vi <- apply(diff.vi,2,function(x) quantile(x,0.25,na.rm=T))
  hqr.diff.vi <- apply(diff.vi,2,function(x) quantile(x,0.75,na.rm=T))
  val.data <- data.frame(mean.val=median.diff.vi,lqr=lqr.diff.vi,hqr=hqr.diff.vi,tf.name=names(median.diff.vi), color.val=(median.diff.vi > 0))
  
  rownames(val.data) <- val.data$tf.name
  val.data <- filter.rows(val.data)
  val.data$tf.name <- standardize.name(val.data$tf.name)
  require(ggplot2)
  axes.format <- opts(plot.title = theme_text(size=12,vjust=1),                    
                      axis.text.x = theme_text(size=16,colour="black"),
                      axis.text.y = theme_text(size=10,colour="black",hjust=1),
                      axis.title.x = theme_text(size=12),
                      axis.title.y = theme_text(size=12,angle=90),
                      legend.position="none",
                      legend.title = theme_text(size=10,hjust=0),
                      legend.text = theme_text(size=10)                      
                      )
  
  p1 <- ggplot(val.data) +     
    geom_bar( aes( x=reorder(tf.name,mean.val) , y=mean.val, fill=color.val), alpha=0.8 ) +
    geom_errorbar( aes( x=reorder(tf.name,mean.val), ymax=hqr, ymin=lqr) )
  axes.labels <- labs(x = "TF", y = "Differential importance") # axes labels
  p1 <- p1 + 
    axes.labels + 
    axes.format + 
    #opts(title=plot.title) +
    coord_flip()
  
  if (nrow(val.data) > 50) {
    p1 <- p1 + opts(axis.text.y = theme_text(size=7,colour="black",hjust=1))
  }

  ggsave(file=output.filename, plot=p1, width=6, height=10, dpi=600)    
  
  return(list(proximal.vi=proximal.vi,
              distal.vi=distal.vi,
              vi.data=val.data))
}

compute.expr.high.low.diff.importance <- function(input.dir, 
                                                  peak.distance.expr.file, 
                                                  output.filename=NULL,
                                                  dist.cutoff=NA,
                                                  low.expr.cutoff=1,                                                    
                                                  high.expr.cutoff=4,
                                                  rm.zero.expr=T) {
  # Computes differential relative importance for all factors comparing low vs. high expression genes associated with peaks

  # Check that input directory exists
  if (! file.exists(input.dir)) {
    stop("Input Directory ", input.dir," does not exist\n")
  }
  
  # Check that peak.distance.expr.file exists
  if (! file.exists(peak.distance.expr.file)) {
    stop("Peak2TSS distance+expression file ", peak.distance.expr.file," does not exist\n")
  }
  
  # Autogenerate outputfile name
  if (is.null(output.filename)) {
    output.stub <- gsub(pattern="\\.+$",replacement="",x=get.file.parts(peak.distance.expr.file)$name)
    output.filename <- file.path(input.dir,
                             sprintf("%s.nozeros.%d.dist.%d.expr.low.%d.high.%d.png", output.stub, rm.zero.expr, dist.cutoff,low.expr.cutoff, high.expr.cutoff))
  }
    
  # Read and parse distance file, get low and high peak ids
  distance.table <- read.table(file=peak.distance.expr.file,
                               header=T,
                               sep="\t",
                               #col.names=c("peak.chr","peak.start","peak.stop","peak.id","tss.chr","tss.start","tss.stop","tss.id","dist"),
                               stringsAsFactors=F)
  # Remove zero expression
  if (rm.zero.expr) {
    distance.table <- droplevels(distance.table[distance.table$tss.cage != 0, ])
  }
  # Remove peaks that are beyond the distance cutoff
  if (!is.na(dist.cutoff)) {
    cat(sprintf("Number of peaks passing distance cutoff = %d of %d\n",sum(distance.table$dist <= dist.cutoff),length(distance.table$dist)))
    distance.table <- droplevels(distance.table[distance.table$dist <= dist.cutoff, ])
  }
  
  # Get low expression and high expression peaks
  low.peak.ids <- distance.table$peak.id[distance.table$tss.cage <= low.expr.cutoff]
  cat(sprintf("Number of peaks classified as LOW expression = %d\n",length(low.peak.ids)))
  high.peak.ids <- distance.table$peak.id[distance.table$tss.cage > high.expr.cutoff]
  cat(sprintf("Number of peaks classified as HIGH expression = %d\n",length(high.peak.ids)))
  
  # Get list of Rdata files in directory
  all.Rdata.files <- list.files(path=input.dir, pattern=".*Rdata$", full.names=T) # Get names of Rdata files
  n.Files <- length(all.Rdata.files)
  if (n.Files == 0) {
    stop("No Rdata files found in input directory", input.dir, "\n")
  }
  
  # Load each Rdata file, compute differential importance and store them
  low.vi <- data.frame()
  high.vi <- data.frame()
  for (each.file in all.Rdata.files) {
    rulefit.results <- restore.rf.model(each.file) # load Rdata.file
    low.rulefit.results <- get.var.imp(rulefit.results,class=low.peak.ids)
    low.vi <- rbind(low.vi, low.rulefit.results$vi)
    high.rulefit.results <- get.var.imp(rulefit.results,class=high.peak.ids)
    high.vi <- rbind(high.vi, high.rulefit.results$vi)
  }
  
  diff.vi <- high.vi - low.vi
  median.diff.vi <- apply(diff.vi,2,function(x) median(x,na.rm=T))  
  lqr.diff.vi <- apply(diff.vi,2,function(x) quantile(x,0.25,na.rm=T))
  hqr.diff.vi <- apply(diff.vi,2,function(x) quantile(x,0.75,na.rm=T))
  val.data <- data.frame(mean.val=median.diff.vi,lqr=lqr.diff.vi,hqr=hqr.diff.vi,tf.name=names(median.diff.vi), color.val=(median.diff.vi > 0))
  
  rownames(val.data) <- val.data$tf.name
  val.data <- filter.rows(val.data)
  val.data$tf.name <- standardize.name(val.data$tf.name)
  require(ggplot2)
  axes.format <- opts(plot.title = theme_text(size=12,vjust=1),                    
                      axis.text.x = theme_text(size=16,colour="black"),
                      axis.text.y = theme_text(size=10,colour="black",hjust=1),
                      axis.title.x = theme_text(size=12),
                      axis.title.y = theme_text(size=12,angle=90),
                      legend.position="none",
                      legend.title = theme_text(size=10,hjust=0),
                      legend.text = theme_text(size=10)                      
                      )
  
  p1 <- ggplot(val.data) +     
    geom_bar( aes( x=reorder(tf.name,mean.val) , y=mean.val, fill=color.val), alpha=0.8 ) +
    geom_errorbar( aes( x=reorder(tf.name,mean.val), ymax=hqr, ymin=lqr) )
  axes.labels <- labs(x = "TF", y = "Differential importance") # axes labels
  p1 <- p1 + 
    axes.labels + 
    axes.format + 
    #opts(title=plot.title) +
    coord_flip()
  
  if (nrow(val.data) > 50) {
    p1 <- p1 + opts(axis.text.y = theme_text(size=7,colour="black",hjust=1))
  }

  ggsave(file=output.filename, plot=p1, width=6, height=10, dpi=600)    
  
  return(list(low.vi=low.vi,
              high.vi=high.vi,
              vi.data=val.data))
}

# =================================================================================================================================
# =================================================================================================================================
# GENE CENTRIC FUNCTIONS
# =================================================================================================================================
# =================================================================================================================================
read.gc.assoc.file <- function( gc.assoc.file, std.thresh=NA ) {
  # ===================================
  # Parses and reads a gene-centric association table (needs to have headers for each column)
  #   First column MUST be 'Gene' representing gene names
  #   Optionally remove low std. dev. columns/partners
  # Returns 
  # $assoc.matrix: dataFrame that is the association matrix (rows: genes, cols: TFs)
  # $target.name: name of target
  # ===================================  
  #gc.assoc.file: association file ( Assumes that file name is of form [Prefix]_SigMtrx_[TargetName].[overlapType].mtrx )
  #std.thresh: columns with stddev. < str.thresh are removed from analysis
  
  assoc.matrix <- read.table( gc.assoc.file, header=TRUE )
  target.name <- gsub( '(^.+/)|(Rgn_TFSig_Mtrx\\.mtrx$)|(\\.[^/]*mtrx$)|(FootprintMtrx_)|(PWMMtrx_)' , '' , gc.assoc.file ) # Remove directory and suffix (.mtrx)
  
  rownames(assoc.matrix) <- assoc.matrix$Gene # Gene names are in column named Gene
  assoc.matrix <- assoc.matrix[ , -1 ] # Remove Gene column  
  
  # Remove columns that have very low std
  if (! is.na( std.thresh )) {
    col.std <- apply( data.matrix(assoc.matrix) , 2 , sd ) # compute std for each column
    remove.col <- which( col.std < std.thresh ) # Find columns that have std < threshold
    if (length(remove.col) > 0) {
      assoc.matrix <- assoc.matrix[ , -remove.col]
    }
  }
  
  return( list(
    assoc.matrix=assoc.matrix,
    target.name=target.name ) )
} # end: read.gc.assoc.file

  
gc.assoc.file.to.Rdata <- function( gc.assoc.file, output.dir=NULL ) {
  # ===================================
  # Converts .mtrx file to a R object and saves in an Rdata file
  # ===================================  
  # gc.assoc.file: gene centric association text file
  # output.dir: directory to store corresponding Rdata files
  
  if( is.null(output.dir) ) {
    output.dir <- get.file.parts(gc.assoc.file)$path # if output.dir is not set make it equal to assoc.dir
  }

  assoc.data <- read.gc.assoc.file(gc.assoc.file)
  output.file <- file.path( output.dir , paste( get.file.parts(gc.assoc.file)$fullname , '.Rdata' , sep="" ) )
  assoc.data$assoc.mtrx.file <- gc.assoc.file
  assoc.data$assoc.R.file <- output.file
  
  save(list="assoc.data",file=output.file)  
}

batch.read.gc.assoc.file.to.Rdata <- function( assoc.dir , output.dir=NULL) {
  # ===================================
  # Reads all .mtrx files in a directory, 
  # converts them to R data frame and stores them
  # as .mtrx.Rdata files
  # ===================================  
  # assoc.dir: directory containing association files
  # output.dir: directory to store corresponding Rdata files
    
  # Search for all .mtrx files in assoc.dir
  gc.assoc.file.paths <- dir( path=assoc.dir , pattern="\\.mtrx$" , full.names=TRUE , recursive=TRUE ) 
  
  for ( each.file in gc.assoc.file.paths ) {
    cat("Processing file " , each.file , "\n")
    try( gc.assoc.file.to.Rdata( each.file, output.dir ) , silent=T )
  }
}

consolidate.expression.data <- function( expr.file, norm.type="asinh", pseudocount=1e-5, process.reps="average" ) {
  # Reads in RNA/CAGE tables, takes log2 transform and then averages replicates
  # expr.file: expression data, different columns for different expression data types
  # norm.type: normalization mode
  #       "none" : no transformation
  #       "log": add pseudocount and then take log2
  #       "logmin" : add minimum value as pseudocount and then use log
  #       "sqrt": square root transform
  #       "asinh" : inverse sinh which can be interpretted similar to log but works for 0 and negative values as well
  #       "normscore": convert to normal scores normx = qnorm((rank(x) - 0.375)/(sum(!is.na(x)) + .25))
  # pseudocount: pseudocount to be added to expr if log is used
  # process.reps: how to process replicates
  #               "average" : will average replicates
  #               "indiv": randomly select one of the reps
  
  expr <- read.table( expr.file , header=T , row.names=1 , sep="\t" )  
  expr.colnames <- colnames(expr)
  expr.colnames <- gsub( "rep[1-9]+" , "rep0", expr.colnames )
  unique.expr.colnames <- unique(expr.colnames)
  idx <- match(expr.colnames, unique.expr.colnames)
  
  ncols.f.expr <- length(unique.expr.colnames)
  final.expr <- expr[, c(1:ncols.f.expr)]
  colnames(final.expr) <- unique.expr.colnames
  
  # Aggregate reps if required
  for (i in c(1 : ncols.f.expr)) {
    curr.idx <- which(idx==i)
    if (process.reps=="indiv") {
      curr.idx <- curr.idx[1]
    }
    if (length(curr.idx) > 1) {
      final.expr[ , i ] <- apply( expr[ , curr.idx], 1, mean )  
    } else {
      final.expr[ , i ] <- expr[ , curr.idx]
    }          
  }
  
  # Normalize expression values
  if (norm.type=="logmin") {
    final.expr[final.expr==0] <- NA
    min.val <- min(final.expr,na.rm=T)
    final.expr <- final.expr + min.val
    final.expr[is.na(final.expr)] <- min.val
    final.expr <- log2(final.expr)    
  } else if (norm.type=="log") {
    final.expr <- log2(final.expr + pseudocount)
  } else if (norm.type=="sqrt") {
    final.expr <- sqrt(final.expr)
  } else if (norm.type=="asinh") {
    final.expr <- asinh(final.expr)
  } else if (norm.type=="normscore") {
    r.names <- rownames(final.expr)
    c.names <- colnames(final.expr)
    final.expr <- as.data.frame(
      apply(final.expr,
            2,
            function (x) qnorm((rank(x,ties.method="random") - 0.375)/(sum(!is.na(x)) + .25)) ) )
    rownames(final.expr) <- r.names
    colnames(final.expr) <- c.names
  }
  
  return(final.expr)  
}

# old.make.expr.classf.data <- function(assoc.data,expr.upper=1,expr.lower=1,regress=FALSE){
#   # ===================================
#   # Create Expression classification dataset
#   # ===================================  
#   # assoc.data$assoc.matrix
#   # assoc.data$target.name
#   
#   x.vals <- assoc.data$assoc.matrix
#   y.vals <- x.vals$expr.val
#   x.vals$expr.val <- NULL
#   
#   if (regress==FALSE){
#     grtr.idx <- (y.vals >= expr.upper)
#     lsr.idx <- (y.vals < expr.lower)
#     y.vals[grtr.idx] <- 1
#     y.vals[lsr.idx] <- -1
#   } else {
#     y.vals <- sqrt(y.vals)
#   }
#     
#   return(list(x.vals=x.vals,y.vals=y.vals,target.name=assoc.data$target.name))
# }

make.gene.centric.tf.to.expr.dataset <- function(assoc.data,    # Gene centric association dataset
                                                 expr,    # expression dataset
                                                 filter.expr=c("Cy","plus","Cage"),   # terms to filter columns of expr by (grep)
                                                 rm.zero.expr=NA  # will remove genes whose expression values are < rm.zero.expr
                                                 ){
  # ===================================
  # Create Gene Centric expression dataset
  # ===================================  
  # assoc.data$assoc.matrix (For relaxed peak thresholds, binding values range from 0 to 2. For non-relaxed, range is 0 to 1)
  # assoc.data$target.name
  
  # Load TF data
  if (is.character(assoc.data)) {
    load(assoc.data)    
  }

  # Remove genes for which there is no TF data
  x.vals <- filter.cols(assoc.data$assoc.matrix,rm.treatments=T)
  x.vals <- x.vals[ (apply(x.vals,1,function(x) sum(as.numeric(x>0),na.rm=T))>-1) , ]
  
  # Load expression data
  if (is.character(expr)) {
    load(expr)    
  }
  
  # Get specific expression data type
  expr.types <- colnames(expr)
#   filter.expr <- c("Cy","plus","Cage")
  for (i in filter.expr) {
    expr.types <- expr.types[grep(i,expr.types,ignore.case=T)]
  }
  if (length(expr.types) > 1) {
    cat("Multiple expr types match .. Choosing first match\n")
    expr.types <- expr.types[1]
  }
  cat(expr.types,"\n")
  y.vals <- expr[,expr.types]
  names(y.vals) <- rownames(expr)
  y.vals <- y.vals - min(y.vals) # Make the y values start at 0
  
  # Remove low expression data if required
  if (!is.na(rm.zero.expr)) {
    y.vals <- y.vals[y.vals>rm.zero.expr]
  }
  
  # Match gene names for expr and tf data
  common.gene.names <- intersect(names(y.vals),rownames(x.vals))
  x.vals <- x.vals[ match(common.gene.names,rownames(x.vals)) , ]
  y.vals <- y.vals[ match(common.gene.names,names(y.vals)) ]
      
  return(list(x.vals=x.vals,y.vals=y.vals,target.name=assoc.data$target.name))
}

learn.tf.to.expr.rulefit.model <- function(assoc.data,    # Gene centric association dataset
                                           expr,    # expression dataset
                                           filter.expr=c("Cy","plus","Cage"),   # terms to filter columns of expr by (grep)
                                           rm.zero.expr=NA,  # will remove genes whose expression values are < rm.zero.expr
                                           two.stage.model=T,
                                           randomize=NA # Set to 0 (randomize rows and cols), 1 (randomize rows), 2 (randomize columns)
                                           ){
  # ===================================
  # Sample a rulefit model
  # Returns
  #   $rfmod: rulefit model
  #   $dataset: sampled dataset
  #   $vi: variable importance (place holder data.frame of n.cols #partners)
  #   $int.strength: interaction strengths (placeholder data.frame of length #partners)
  #   $pair.interactions: pairwise interactions (placeholder data.frame of size #partners X #partners)
  # ===================================  
  # assoc.data$assoc.matrix
  # assoc.data$target.name
  tree.size=10
  test.reps=3
  
  if (is.character(assoc.data)) {
    load(assoc.data)
  }
  
  if (!is.na(randomize)) {
    assoc.data$assoc.matrix <- randomize.assoc.matrix(assoc.data,rand.dim=randomize,change.row.names=F)
  }
  
  assoc.classf.data <- make.gene.centric.tf.to.expr.dataset(assoc.data, expr, filter.expr, rm.zero.expr)
  
  if (two.stage.model) {
    class.dataset <- assoc.classf.data
    min.val <- min(class.dataset$y.vals)
    class.dataset$y.vals[class.dataset$y.vals > min.val] <- 1 # Set non-zero values to 1
    class.dataset$y.vals[class.dataset$y.vals <= min.val] <- -1 # Set non-zero values to 1
    rfmod.class <- run.rulefit(class.dataset, mode="class",tree.size=tree.size,test.reps=test.reps) # Learn classification model
    
    regress.dataset <- assoc.classf.data
    keep.idx <- (regress.dataset$y.vals > min.val) # Only keep non-zero values
    regress.dataset$x.vals <- regress.dataset$x.vals[keep.idx,]
    regress.dataset$y.vals <- regress.dataset$y.vals[keep.idx]
    rfmod.regress <- run.rulefit(regress.dataset, mode="regress",tree.size=tree.size,test.reps=test.reps)
  } else {
    rfmod <- run.rulefit(assoc.classf.data, mode="regress",tree.size=tree.size,test.reps=test.reps)
  } 
  
  # Create place holder for variable importance
  partner.names <- colnames(assoc.classf.data$x.vals)
  n.partners <- length(partner.names)
  vi <- as.data.frame( matrix( data=NA, nrow=1, ncol=n.partners) )
  colnames(vi) <- partner.names
  
  # Create place holder for interaction strengths
  int.strength <- data.frame(matrix( data=NA , nrow=1 , ncol=n.partners ))
  colnames(int.strength) <- partner.names
  
  # Create place holder for pairwise interactions
  pair.interactions <- data.frame(matrix( data=NA , nrow=n.partners , ncol=n.partners ))
  rownames(pair.interactions) <- partner.names
  colnames(pair.interactions) <- partner.names
  
  if (two.stage.model) {
    
    # Create classification results
    class.results <- list(rfmod=rfmod.class,
                          dataset=class.dataset,
                          vi=vi,
                          int.strength=int.strength,
                          pair.interactions=pair.interactions)
    class.results <- run.cv.rulefit(class.results)
    
    # Create regression results
    regress.results <- list(rfmod=rfmod.regress,
                            dataset=regress.dataset,
                            vi=vi,
                            int.strength=int.strength,
                            pair.interactions=pair.interactions)
    regress.results <- run.cv.rulefit(regress.results)
    
    # Combine predictions
    regress.results <- restore.rf.model(regress.results)
    y.regress <- rfpred(class.results$dataset$x.vals)    
    y.pred <- as.numeric(class.results$cv$lo > 0) * y.regress
    
    combined.results <- list(dataset=assoc.classf.data,
                             y.pred=y.pred)
    
    rulefit.results <- list(class.results=class.results,
                            regress.results=regress.results,
                            combined.results=combined.results)
    
  } else {
    
    rulefit.results <- list(rfmod=rfmod,
                            dataset=assoc.classf.data,
                            vi=vi,
                            int.strength=int.strength,
                            pair.interactions=pair.interactions)
    rulefit.results <- run.cv.rulefit(rulefit.results)
  }
  
  return( rulefit.results )
}

# =================================================================================================================================
# =================================================================================================================================
# HISTONE MARK FUNCTIONS
# =================================================================================================================================
# =================================================================================================================================
read.histone.assoc.file <- function( histone.assoc.file, std.thresh=NA, col.types="all") {
  # ===================================
  # Parses and reads a gene-centric association table (needs to have headers for each column)
  #   First column MUST be 'Gene' representing gene names
  #   Optionally remove low std. dev. columns/partners
  # Returns 
  # $assoc.matrix: dataFrame that is the association matrix (rows: genes, cols: histone magnitude and shape)
  # $target.name: name of target
  # ===================================  
  #histone.assoc.file: association file
  #std.thresh: columns with stddev. < str.thresh are removed from analysis
  # col.types: "all", "mag", "shape"
  
  assoc.matrix <- read.table( histone.assoc.file, header=TRUE, na.strings=c("NA","NaN","") )
  target.name <- gsub( '(^.+/)|(Rgn_TFSig_Mtrx\\.mtrx$)|(\\.[^/]*mtrx$)|(FootprintMtrx_)|(PWMMtrx_)' , '' , histone.assoc.file ) # Remove directory and suffix (.mtrx)
  
#   rownames(assoc.matrix) <- assoc.matrix$Gene # Gene names are in column named Gene
#   assoc.matrix <- assoc.matrix[ , -1 ] # Remove Gene column  
  
  # Remove columns that have very low std
  if (! is.na( std.thresh )) {
    col.std <- apply( data.matrix(assoc.matrix) , 2 , sd ) # compute std for each column
    remove.col <- which( col.std < std.thresh ) # Find columns that have std < threshold
    if (length(remove.col) > 0) {
      assoc.matrix <- assoc.matrix[ , -remove.col]
    }
  }
  
  # Filter columns based on column type
  if (col.types == "mag") {
    assoc.matrix <- assoc.matrix[, grep( pattern=".*Mag.*", x=colnames(assoc.matrix) ) ]    
  } else if (col.types == "shape") {
    assoc.matrix <- assoc.matrix[, grep( pattern=".*Shape.*", x=colnames(assoc.matrix) ) ]
  }
  
  return( list(
    assoc.matrix=assoc.matrix,
    target.name=target.name ) )
} # end: read.histone.assoc.file

  
histone.assoc.file.to.Rdata <- function( histone.assoc.file, output.dir=NULL, col.types="all" ) {
  # ===================================
  # Converts .mtrx file to a R object and saves in an Rdata file
  # ===================================  
  # histone.assoc.file: gene centric association text file
  # output.dir: directory to store corresponding Rdata files
  
  if( is.null(output.dir) ) {
    output.dir <- get.file.parts(histone.assoc.file)$path # if output.dir is not set make it equal to assoc.dir
  }

  assoc.data <- read.histone.assoc.file(histone.assoc.file, col.types=col.types)
  output.file <- file.path( output.dir , paste( get.file.parts(histone.assoc.file)$fullname , '.', col.types, '.Rdata' , sep="" ) )
  assoc.data$assoc.mtrx.file <- histone.assoc.file
  assoc.data$assoc.R.file <- output.file
  
  save(list="assoc.data",file=output.file)  
}

batch.read.histone.assoc.file.to.Rdata <- function( assoc.dir , output.dir=NULL, col.types="all") {
  # ===================================
  # Reads all .mtrx files in a directory, 
  # converts them to R data frame and stores them
  # as .mtrx.Rdata files
  # ===================================  
  # assoc.dir: directory containing association files
  # output.dir: directory to store corresponding Rdata files
    
  # Search for all .mtrx files in assoc.dir
  histone.assoc.file.paths <- dir( path=assoc.dir , pattern="\\.mtrx$" , full.names=TRUE , recursive=TRUE ) 
  
  for ( each.file in histone.assoc.file.paths ) {
    cat("Processing file " , each.file , "\n")
    try( histone.assoc.file.to.Rdata( each.file, output.dir, col.types=col.types) , silent=T )
  }
}

make.histone.to.tf.regression.dataset <- function(assoc.data, y=NULL, inverted=T, logval=F, filter.dup.var=F) {
  # ===================================
  # Creates an input dataset for rulefit x: histone marks, y: normal scores of TF peak ranks or supplied y values
  # y: if not supplied then y is inferred from the rownames of assoc.data$assoc.matrix
  # inverted: if set to T then lower values of provided or inferred values of 'y' are considered higher ranks
  # ===================================
  # Returns a list with members
  #  $x.vals : feature matrix
  #  $y.vals : labels
  #  $target.name : taken from assoc.data

  if (is.character(assoc.data)) {
    load(assoc.data)    
  }
  
  x.vals <- assoc.data$assoc.matrix
  
  # Filter predictors as required
  x.vals <- filter.cols(x.vals) # Remove unwanted columns
  #x.vals <- x.vals[, -grep(pattern=".*Dnase.*",x=colnames(x.vals))] # Removes DNase predictors
  
  #x.vals[x.vals<2] <- 0
  
  # Convert predictor values to log scale if required
  if (logval) {
    x.vals <- log2(x.vals+0.1)
  }
  
  
  # Convert row names of x to labels y
  if (is.null(y)) {
    y.vals <- rownames(x.vals)
    y.vals <- gsub( pattern=".+_([0-9]+)$" , replacement="\\1" , x=y.vals )
    y.vals <- as.numeric(y.vals)
  } else {
    y.vals <- as.numeric(y)
    if (length(y.vals) != nrow(x.vals)) {
      stop("number of rows of X not matching number of elements in Y\n")
    }
  }  
  
  # Convert y to normal score ranks
  y.vals <- rank(y.vals,na.last="keep")
  if (inverted) {
    y.vals <- rev(y.vals)
  }
  y.vals <- get.normal.score(y.vals)
    
  return( list(
    x.vals=x.vals,
    y.vals=y.vals,
    target.name=assoc.data$target.name ) )  
}
  
learn.histone.to.tf.regression.model <- function(assoc.data, inverted=T){
  # ===================================
  # Learns a rulefit model that regresses histone marks (magnitude &/or shapes to the TF peak ranks)
  #   $rfmod: rulefit model
  #   $dataset: sampled dataset
  #   $vi: variable importance (place holder data.frame of n.cols #partners)
  #   $int.strength: interaction strengths (placeholder data.frame of length #partners)
  #   $pair.interactions: pairwise interactions (placeholder data.frame of size #partners X #partners)
  # ===================================  
  # assoc.data$assoc.matrix
  # assoc.data$target.name
  # rm.target: if set to TRUE then target TF is not used in constructing the model
    
  hist.regress.data <- make.histone.to.tf.regression.dataset(assoc.data,inverted=inverted)
  rfmod <- run.rulefit(hist.regress.data,mode="regress")
  
  # Create place holder for variable importance
  predictor.names <- colnames(hist.regress.data$x.vals)
  n.predictors <- length(predictor.names)
  vi <- as.data.frame( matrix( data=NA, nrow=1, ncol=n.predictors) )
  colnames(vi) <- predictor.names
  
  # Create place holder for interaction strengths
  int.strength <- data.frame(matrix( data=NA , nrow=1 , ncol=n.predictors ))
  colnames(int.strength) <- predictor.names
  
  # Create place holder for pairwise interactions
  pair.interactions <- data.frame(matrix( data=NA , nrow=n.predictors , ncol=n.predictors ))
  rownames(pair.interactions) <- predictor.names
  colnames(pair.interactions) <- predictor.names
  
  return( list(
    rfmod=rfmod,
    dataset=hist.regress.data,
    vi=vi,
    int.strength=int.strength,
    pair.interactions=pair.interactions) )
}

merge.histone.tf.datasets <- function(tf.assoc.data, hist.assoc.data, output.dir=NA) {
  # Merges a tf coassociation dataset for a target TF with a histone coassociation dataset for the target TF
  
  # Load TF data if required
  if (is.character(tf.assoc.data)) {
    load(tf.assoc.data)
    tf.assoc.data <- assoc.data
  }
  
  # Load histone data if required
  if (is.character(hist.assoc.data)) {
    load(hist.assoc.data)
    hist.assoc.data <- assoc.data
  }
  
  # NOTE: TF datasets have rowname indices flipped so sort them opposite to histone data
  
  # Check that number of rows are the same
  if( nrow(hist.assoc.data$assoc.matrix) != nrow(tf.assoc.data$assoc.matrix) ) {
    stop("number of rows of TF dataset not matching number of rows of histone dataset\n")
  }
  
  # Sort rows and then merge
  tf.row.ids <- as.numeric(gsub(pattern=".*_([0-9]+)_?.*",replacement="\\1",x=rownames(tf.assoc.data$assoc.matrix)))
  tf.row.ids <- length(tf.row.ids) - tf.row.ids + 1
  hist.row.ids <- as.numeric(gsub(pattern=".*_([0-9]+)_?.*",replacement="\\1",x=rownames(hist.assoc.data$assoc.matrix)))
  hist.row.perm <- match(hist.row.ids, tf.row.ids)
  merged.assoc.matrix <- cbind(tf.assoc.data$assoc.matrix, hist.assoc.data$assoc.matrix[hist.row.perm,])
  assoc.data <- list(assoc.matrix=merged.assoc.matrix,
                     target.name=tf.assoc.data$target.name,
                     hist.target.name=hist.assoc.data$target.name)
  
  if (!is.na(output.dir)) {
    save(list="assoc.data", file=file.path(output.dir, paste(tf.assoc.data$target.name,".mtrx.tfhist.Rdata",sep="")))
  }
  
  return(assoc.data)              
}

# =======================================
# TODO
# (1) Affinity matrix using random forests, clustering using fuzzy c-means
# (2) Proximal vs distal variable importance
# (3) Positive vs Negative datasets for TFs
# (4) Regression for TFs
# (5) TF to expression models, importance, pairwise interactions, pdf plots
# =======================================

sample.unsupervised.proximity.matrix <- function(assoc.data){
  # ===================================
  # Use a random forest model to learn an affinity matrix
  # Returns
  # $assoc.data
  # $proximity: average proximity matrix  
  # ===================================  
  # assoc.data$assoc.matrix
  # assoc.data$target.name
  
  library(randomForest)
  
  if (is.character(assoc.data)) {
    load(assoc.data)    
  }
  
  # Remove columns with all NAs
  assoc.data$assoc.matrix[, which(apply(is.na(assoc.data$assoc.matrix),2,all)) ] <- NULL
  
  # Learn unsupervised random forest model
  forest <- randomForest(x=assoc.data$assoc.matrix,proximity=T)
  cat("Error rate of random forest = ", forest$err.rate, "\n")
    
  return( list(
    proximity=forest$proximity,
    assoc.data=assoc.data
    ) )
}

aggregate.unsupervised.proximity.matrices <- function(prox.dir) {
  # ===================================
  # Aggregates proximity matrices from multiple runs (*.proximity.Rdata) in a single directory
  # Returns
  # $assoc.data
  # $proximity: average proximity matrix
  # ===================================

  if (! file.exists(prox.dir)) {    
    stop(cat("Input Directory ", prox.dir,"does not exist\n"))
  }
  
  all.Rdata.files <- list.files(path=input.dir, pattern=".*proximity.Rdata$", full.names=T) # Get names of Rdata files
  all.Rdata.files <- all.Rdata.files[! grepl("average",all.Rdata.files) ]
  n.Files <- length(all.Rdata.files)
  if (n.Files == 0) {
    stop(cat("No files matching .*proximity.Rdata in directory\n"))
  }
  
  proximity <- NULL
  count=0
  for ( curr.file in all.Rdata.files ) {
    cat(curr.file,"\n")
    # load file (random.forests.results)
    load(curr.file)
    count <- count + 1
    if (is.null(proximity)) {
      proximity <- random.forest.results$proximity
    } else {
      proximity <- proximity + random.forest.results$proximity
    }
  }
  proximity <- proximity / count
  return( list(
    proximity=proximity,
    assoc.data=random.forest.results$assoc.data))
}

plot.MDS <- function (random.forest.results, labels=NULL, k = 2, eig.cutoff=0.05, make.plot=T, palette = NULL, pch = '.', cex=3, ...) {
  # ===================================
  # Creates a multi-dimensional scaling plot from a proximity matrix from random forests
  # Returns MDS coordinates for each of the points
  # random.forest.results: output of random forest run
  # $assoc.data
  # $proximity
  # labels: list of classification or regression labels
  # k: number of dimensions to consider (if set to NULL, the code will try to optimize k to explain as much variance as possible)
  # eig.cutoff: for auto-tuning 'k', select all dimensions whose eig is within 0.05 of the max eig
  # make.plot: if set to T, MDS plots will be generated
  # ===================================

  n.vars <- ncol(random.forest.results$assoc.data$assoc.matrix) # number of original variables/predictors
  
  if (is.null(k)) {
    # Optimize k
    rf.mds <- stats:::cmdscale(1 - proximity, eig = TRUE, k = 2)
    eig.dist <- rf.mds$eig / max(rf.mds$eig)
    k <- sum(eig.dist >= eig.cutoff)
    if (k > n.vars) { k <- n.vars - 1}
    if (k < 2) { k <- 2}
    rf.mds <- stats:::cmdscale(1 - proximity, eig = TRUE, k = k)
  } else {
    rf.mds <- stats:::cmdscale(1 - proximity, eig = TRUE, k = k)
  }  
  colnames(rf.mds$points) <- paste("Dim", 1:k)
  
  # Disable plotting if k is > 20
  if (k > 20) {
    make.plot <- F
  }
  if (make.plot) {
    op <- par(pty = "s")
    on.exit(par(op))
    if (is.null(labels)) {
      nlevs <- 1
      labels <- c(1:nrow(proximity))*0 + 1
    } else {
      nlevs <- nlevels(labels) 
    }    
    if (is.null(palette)) {
      palette <- if (require(RColorBrewer) && nlevs < 12 && nlevs > 3) 
        brewer.pal(nlevs, "Set1")
      else rainbow(nlevs)
    }
    if (k <= 2) {
      plot(rf.mds$points, col = palette[as.numeric(labels)], pch = pch, cex=cex, ...)
    } else {
      pairs(rf.mds$points, col = palette[as.numeric(labels)], pch = pch, cex=cex, ...)
    }
  }
  invisible(rf.mds)
}
  
# get.proximity.dendrogram <- function(random.forest.results, clust.method="ward", mds.flag=F) {
#   # ===================================
#   # Gets cluster dendrogram for a proximity matrix
#   # random.forest.results: output of random forest run
#   # $assoc.data
#   # $proximity
#   # clust.method: linkage method
#   # mds.flag: if set to T, then first multi-dimensional scaling is performed on the
#   # ===================================
#   library(fastcluster)
#   
#   if ()
# }