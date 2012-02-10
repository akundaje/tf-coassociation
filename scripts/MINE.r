library("rJava")
.jinit(classpath="MINE.jar")

MINE <- function (
		input.filename,
		method=c("master.variable", "all.pairs", "adjacent.pairs"),
		master.variable=NA,
		permute.data=F,
		max.num.boxes.exponent=0.6,
		required.common.vals.fraction=0,
		job.id=NA,
		gc.wait=Inf,
		num.clumps.factor=15,
		debug.level=0) {

	if (gc.wait==Inf) {
		gc.wait <- J("java.lang.Integer")$MAX_VALUE
	} else {
		gc.wait <- as.integer(gc.wait)
	}
	
	if (is.na(job.id)) {
		job.id <- paste(
			"B=n^",
			max.num.boxes.exponent,
			",k=",
			num.clumps.factor,
			sep="")
		if (!is.na(master.variable)) {
			job.id <- paste(job.id, ",mv=", master.variable, sep="")
		}
		if (permute.data) {
			job.id <- paste(job.id, "permuted", sep=",")
		}
	}

	if(missing(input.filename)) {
		stop("you must specify an input filename via the parameter input.filename")
	}

	if(missing(method)) {
		stop("you must specify a value for the parameter 'method'")
	}

	if (permute.data && method=="all.pairs") {
		stop("You cannot specify both permute and all.pairs")
	}
	
	if (debug.level < 0 | debug.level > 4) {
		stop("Invalid debug level (must be between 0 and 4 inclusive)")
	}
	
	debug.stream <- .jnew(
		"java/io/PrintWriter", 
		.jcast(J("java/lang/System")$err, "java/io/OutputStream"))

	debug.out <- .jnew(
		"java/io/BufferedWriter",
		.jcast(debug.stream, "java/io/Writer"))

	debug.level <- as.integer(debug.level)
	dataset <- .jnew(
		"data.Dataset", 
		input.filename, 
		debug.level,
		debug.out)

	analysis <- NULL
	if (method=="all.pairs") {
		analysis <- .jnew(
			"analysis.Analysis",
			dataset,
			J("analysis.Analysis")$AnalysisStyle$allPairs)
	} else if (method=="adjacent.pairs") {
		analysis <- .jnew(
			"analysis.Analysis",
			dataset,
			J("analysis.Analysis")$AnalysisStyle$adjacentPairs)
	} else if (method=="master.variable"){
		if(is.na(master.variable)) {
			stop("You cannot specify method=\"master.variable\" without specifying a value for the parameter 'master.variable'")
		}
		analysis <- .jnew(
			"analysis.Analysis",
			dataset,
			as.integer(master.variable))
	} else {
		stop("the 'method' parameter must be either \"all.pairs\", \"adjacent.pairs\", or \"master.variable\"")
	}
	
	if (!is.null(analysis)) {
		results <- analysis$getSortedResults(
			J("main/BriefResult")$class,
			input.filename, 
			.jfloat(required.common.vals.fraction),
			.jfloat(max.num.boxes.exponent),
			.jfloat(num.clumps.factor),
			gc.wait,
			job.id,
			debug.level,
			debug.out)
		.jcall("main/Analyze", "V", "printResults", results, input.filename, job.id)

		cat("Analysis finished. See file \"")
		cat(input.filename,",",job.id,",Results.csv\"",sep="")
		cat(" for output\n")
	}
}

