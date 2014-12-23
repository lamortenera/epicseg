#this is not the actual complete list of options
#for segmentCLI: the complete list is a merge between
#the two lists 'segmentOptions' and 'reportOptions' 
#which is done in the first lines of 'segmentCLI
getSegmentOptions <- function(){list(
	list(arg="--counts", type="character", required=TRUE, parser=readCounts,
	help="Path to the count matrix"),
	list(arg="--regions", type="character", required=TRUE, parser=readRegions,
	help="Path to the .bed file with the genomic regions associated to 
	the count matrix. Only the first three fields of the bed file
	will be read. (required)"),
	list(arg="--nstates", type="integer",
	help="Number of states to use for training (required if the option '--model'
	is not set)"),
	list(arg="--model", type="character", parser=readModel,
	help="Path to the file with the parameters of the HMM.
	In train mode, the model contains the initial parameters for
	learning, some fields can be left out. In prediction mode
	all fields must be present and the model is used for producing
	the segmentation without training.
	(required if the option '--nstates' is not set)"),
	list(arg="--predict", flag=TRUE,
	help="Toggle the predict mode and untoggle the train mode"),
	list(arg="--collapseInitP", type="logical", meta="true_or_false",
	help="In case a model with multiple initial probabilities is provided,
	should those probabilities be averaged? If you are not sure about what
	this means, don't set this option (the default is FALSE in train mode
	and TRUE in predict mode)"),
	list(arg="--nthreads", type="integer",
	help="Number of threads to be used (default is 1)"),
	list(arg="--maxiter", type="integer", 
	help="Maximum number of iterations in train mode (default is 200)"),
	list(arg="--save_rdata", flag=TRUE,
	help="Save an R data archive with the results of the segmentation. To
	see what they are, inside R, type 'library(epicseg)' and '?segment'"))}

segmentCLI <- function(args, prog){
	#name the options in reportOptions
	reportOptions <- makeOptions(getReportOptions())
	#options for 'reportCLI' to use also in 'segmentCLI'
	computedOptNames <- c("segments", "model", "counts")
	subReportOptNames <- setdiff(names(reportOptions), computedOptNames)
	#merge the two option lists
	segmentOptions <- c(getSegmentOptions(), reportOptions[subReportOptNames])
	#parse the options
	opt <- parseArgs(segmentOptions, args, prog)
	#check required arguments
	if (is.null(opt$nstates) && is.null(opt$model)) stop("provide either 'nstates' or 'model'")
	#split the options for 'segment' from the options for 'report'
	sopt <- opt[setdiff(names(opt), subReportOptNames)]
	ropt <- opt[intersect(names(opt), subReportOptNames)]
	#deal with 'predict' (its name in the 'segment' function is 'train', so it's the opposite)
	if (!is.null(sopt$predict)) sopt$train <- !sopt$predict; sopt$predict <- NULL
	#deal with 'save_rdata'
	save_rdata <- !is.null(sopt$save_rdata)
	sopt$save_rdata <- NULL
	#call 'segment'
	segmentation <- do.call(segment, sopt)
	#set the computed options for 'report'
	ropt$segments <- segmentation$segments
	ropt$model <- segmentation$model
	ropt$counts <- sopt$counts
	if (save_rdata) ropt$rdata <- segmentation
	#call 'report'
	cat("producing report\n")
	do.call(report, ropt)
	
}

#these options are not available from the CLI and can be set in 'segment' with the ellipsis (...)
advancedOpts <- c("tol","verbose","nbtype","init","init.nlev")
#' Learn a model and produce a segmentation
#'
#' @param counts Count matrix matching with the \code{regions} parameter.
#' Each row of the matrix represents a mark and each column a bin
#' resulting from dividing the genomic regions into non-overlapping
#' bins of equal size. The rows of the matrix must be named with the 
#' name of the mark and the names must be unique.
#' @param regions GRanges object containing the genomic regions of interest.
#' 	Each of these regions corresponds to a set of bins and each bin to a column
#' 	of the count matrix. The binsize is automatically derived by comparing
#' 	the columns of the count matrix with the width of the regions.
#' @param nstates Number of states to learn.
#' @param model A list with the parameters that describe
#' the HMM. Missing parameters will be learned, and the provided
#' parameters will be used as initial parameters for the learning
#' algorithm. If \code{train==FALSE} the parameter set must be 
#' complete and no learning will take place.
#' @param train If TRUE, the parameters will be learned, otherwise
#' the provided parameters (with the \code{model} option) will be 
#' used without learning to produce a segmentation.
#' @param collapseInitP In case a model with multiple initial probabilities
#' is provided, should those probabilities be averaged and reduced to
#' one initial probabilities vector? If you are not sure
#' about what this means, don't set this option.
#' @param nthreads number of threads used for learning
#' @param maxiter Maximum number of iterations for learning.
#' @param ... Advanced options for learning. Type
#' \code{epicseg:::advancedOpts} to see which options are allowed,
#' and type \code{?kfoots} to see what the options do.
#' @return A list with the following arguments:
#' 	\item{segments}{The segmentation as a GRanges object.
#' 	The slot \code{names} of this object contains a number from
#' 	1 to nstates saying which state each segment belongs to.}
#' 	\item{model}{A list containing all the parameters of the model.}
#' 	\item{posteriors}{A matrix of size \code{nstates*ncol(counts)} containing the posterior
#'			probability that a given datapoint is generated by the given state}
#' 	\item{states}{An integer vector of length \code{ncol(counts)} saying
#' 	which state each bin is associated to (using the posterior decoding
#' 	algorithm). This vector is used to create the \code{segments} argument.}
#' 	\item{viterbi}{Same as \code{states}, but using the viterbi algorithm.}
#' 	\item{loglik}{the log-likelihood of the whole dataset.}
#' @export
segment <- function(counts, regions, nstates=NULL, model=NULL, train=TRUE, collapseInitP=!train, 
										nthreads=1, maxiter=200, ...){
	#REFORMAT KFOOTS OPTIONS
	kfootsOpts <- list(...)
	if (length(kfootsOpts)>0){
		optNames <- names(kfootsOpts)
		if (is.null(optNames) || any(is.na(optNames))) stop("the advanced options must be named")
		if (!all(optNames %in% advancedOpts)) stop("check 'epicseg::advancedOpts' for the list of allowed advanced options")
	}
	kfootsOpts <- c(list(nthreads=nthreads, framework="HMM"), kfootsOpts)
	#CHECK ARGUMENTS
	validateCounts(counts)
	if (is.null(regions) || !inherits(regions, "GRanges")) stop("'regions' must be a GenomicRanges object")
	#check consistency between regions and counts
	binsize <- checkBinsize(regions, ncol(counts))
	if (is.null(model) && !train) stop("no model provided, train mode necessary")
	if (!is.null(model)){
		#if we can't figure out the number of states from the model
		#we can completely discard it (it is an empty model)
		dims <- validateModel(model, input=train)
		if (is.null(dims$nstates)) model <- NULL
	}
	
	if (!is.null(model)) {
		#figure out the number of states
		if (!is.null(nstates) && nstates != dims$nstates) {
			warning("inconsistent 'nstates' given, using the value provided in the model")
			nstates <- dim$nstates
		}
		if (!is.null(model$marks) && !is.null(model$emisP)) model <- matchModelToCounts(model, counts)
		if (!is.null(model$initP)){
			#adapt the initPs to the current set of observations
			if (!collapseInitP && ncol(model$initP) != length(regions)){
				warning("collapsing the initial probabilities for compatibility")
				collapseInitP <- TRUE
			}
			if (collapseInitP){
				avgInitP <- rowMeans(model$initP)
				model$initP <- matrix(avgInitP, nrow=length(avgInitP), ncol=length(regions))
			}
		}
		nstates <- model$nstates
		if (is.null(model$emisP)) kfootsOpts$k <- nstates
		else kfootsOpts$k <- model$emisP
		kfootsOpts$trans <- model$transP
		kfootsOpts$initP <- model$initP
	} else if (!is.null(nstates)){
		kfootsOpts$k <- nstates
	} else stop("provide either an existing HMM or a desired number of states")
	
	kfootsOpts$maxiter <- maxiter
	#seqlens depend on the genomic regions
	kfootsOpts$seqlens <- width(regions)/binsize
	#pass the count matrix
	kfootsOpts$counts <- counts
	
	#deal with the 'train' option
	if (!train){
		kfootsOpts$maxiter <- 1
		kfootsOpts$verbose <- FALSE
	}
	
	#CALL KFOOTS
	fit <- do.call(kfoots, kfootsOpts)
	
	#REORDER STATES
	fit <- reorderFitStates(fit)
	#MAKE SEGMENTS
	segms <- statesToSegments(fit$clusters, regions)
	#MAKE MODELS
	if (train){
		model <- list(nstates=length(fit$models), marks=rownames(counts), emisP=fit$models, transP=fit$trans, initP=fit$initP)
	}
	
	list(segments=segms, model=model, posteriors=fit$posteriors, states=fit$clusters, viterbi=fit$viterbi$vpath, loglik=fit$loglik)
}


reorderFitStates <- function(fit){
	nstates <- length(fit$models)
	#determine a distance matrix between states
	mps <- t(log(getMeanMatrix(fit$models) + 1))
	dmat <- t(as.matrix(dist(mps, method="euclidean")))
	#get a small weight path through all the states
	path <- smallWeightHamiltonianPath(dmat)
	if (length(path) != nstates) stop("something went wrong..")
	#let's make state 1 the one with more counts
	if (sum(mps[,path[1]]) < sum(mps[,path[nstates]])){
		path <- rev(path)
	}
	
	#let's permute every object now
	fit$models <- fit$models[path, drop=F]
	fit$trans <- fit$trans[path, path, drop=F]
	fit$initP <- fit$initP[path,, drop=F]
	fit$posteriors <- fit$posteriors[path,, drop=F]
	#compute inverse permutation
	io <- integer(nstates); io[path] <- 1:nstates
	fit$clusters <- io[fit$clusters]
	fit$viterbi$vpath <- io[fit$viterbi$vpath]
	
	fit
}

#METHODS TO VALIDATE AND SANITIZE ARGUMENTS
checkBinsize <- function(regions, nbins){
	genomicSpan <- sum(as.numeric(width(regions)))
	if (genomicSpan %% nbins != 0) stop("total genomic span is not a multiple of the number of bins")
	binsize <- genomicSpan/nbins
	if (any((width(regions) %% binsize) != 0)) stop("regions widths are not multiples of binsize")
	binsize
}

#make sure that the rownames of counts 
#are in the same order as the marks in the model
matchModelToCounts <- function(model, counts){
	if (!setequal(model$marks, rownames(counts))) {
		stop("provided models have different names than the rows of the 'counts' matrix")
	}
	mperm <- match(rownames(counts), model$marks)
	#this function reorders the 'marks' and 'emisP' fields of the model
	model <- reoderMarks(model, mperm)
	model
}



#MISC
statesToSegments <- function(states, regions){
	h <- statesToSegments_helper(regions, states)
	GRanges(seqnames=h$chrs, IRanges(start=h$starts, end=h$ends, names=h$states))
}


