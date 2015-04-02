getGetcountsOptions <- function(){list(
	list(arg="--regions", type="character", required=TRUE, parser=readRegions,
	help="Path to the bed file with the genomic regions of interest.
	These regions will be automatically partitioned into smaller, 
	consecutive bins. Only the first three fields in the file matter. 
	If the region lengths are not multiples of the given binsize
	a new bed file will be produced where each coordinate 
	is a multiple of binsize. Use this new file together with
	the count matrix for later analyses."),
	list(arg="--mark", type="character", required=TRUE, vectorial=TRUE, parser=readMarks, meta="label:path",
	help="Mark name and path to a bam file where to extract reads from.
	The bam files must be indexed and the chromosome names must match with 
	those specified in the bed file. 
	This option must be repeated for each mark, for example:
	'-m H3K4me3:/path1/foo1.bam -m H3K36me3:/path2/foo2.bam'"),
	list(arg="--target", type="character", required=TRUE, parser=validatePath,
	help="Full path to the count matrix that will be produced as output.
	To save the matrix as an R archive, provide a path that ends
	in '.Rdata' or '.rda', otherwise it will be saved as a text file.
	Using R archives, writing the file will be a bit slower,
	but reading it will be considerably faster. (required)"),
	list(arg="--binsize", type="integer",
	help="Size of a bin in base pairs. Each given region will be partitioned into
	bins of this size. (default 200)"),
	list(arg="--mapq", type="integer", vectorial=TRUE,
	help="Minimum mapping quality for the reads (see the bam format
	specification for the mapq field). Only reads with the mapq field
	above or equal to the specified value will be considered. 
	Specify this option once for all marks, or specify it as many times 
	as there are marks. (default 0)"),
	list(arg="--pairedend", type="logical", vectorial=TRUE, meta="true_or_false",
	help="Set this option to TRUE or FALSE to activate or deactivate
	the paired-end mode. Only read pairs where both ends
	are mapped will be considered and assigned to the bin where the 
	midpoint of the read pair is located. Specify this option 
	once for all marks, or specify it as many times as there are marks.
	If this flag is set, the 'shift' option will be ignored. 
	(default FALSE)"),
	list(arg="--shift", type="integer", vectorial=TRUE,
	help="Shift the reads in the 5' direction by a fixed number
	of base pairs. The read will be assigned to the bin where the 
	shifted 5' end of the read is located. Specify this option 
	once for all marks, or specify it as many times as there are marks.
	This option will be ignored in paired-end mode. (default 75)"))}

getcountsCLI <- function(args, prog){
	opt <- parseArgs(getGetcountsOptions(), args, prog)
	#remove the target field (it's not an argument of the R interface)
	target <- opt$target; opt$target <- NULL
	#we need to know the binsize in order to know if the regions have to be refined or not
	if (is.null(opt$binsize)) opt$binsize <- formals(getcounts)$binsize
	#check if the regions are all right
	if (any(width(opt$regions) %% opt$binsize != 0)) {
		#regions must be refined
		targetext <- file_ext(target)
		rtarget <- sub(paste0("\\.", targetext, "$"), "_refined_regions.bed", target)
		cat(sep="", 
		"[check regions] The regions must be refined in order to be\n",
		"[check regions] compatible with the chosen binsize.\n",
		"[check regions] Writing refined regions to the file\n",
		"[check regions] '", rtarget, "'\n",
		"[check regions] Use this file for later analyses\n")
		opt$regions <- refineRegions(opt$regions, opt$binsize)
		writeRegions(opt$regions, rtarget)
	}
	#call getcounts
	counts <- do.call(getcounts, opt)
	#write results
	cat(sep="", "writing count matrix to the file '", target, "'\n")
	writeCounts(counts, target)
}

#' Make a count matrix
#'
#' Make a count matrix from a GRanges object and a list of bam files
#' @param regions GRanges object containing the genomic regions of interest.
#' 	Each region will be divided into non-overlapping bins of equal size.
#' 	The reads falling in each bin will be the elements of the final count matrix.
#' @param marks Named list where each item is a path to an indexed bam file.
#' 	The chromosome names in the bam file must match with the sequence names
#' 	in the GRanges object.
#' @param binsize The size of each bin in basepairs. Each region must have
#' 	a width multiple of \code{binsize}, otherwise an error will be thrown.
#' 	Use the function \code{refineRegions} to make sure that your GRanges object
#' 	satisfies this constraint.
#' @param mapqs Threshold on the 'mapq' field in the bam format. Reads
#' 	with a mapping quality strictly lower than 'mapq' will be discarded.
#' @param pairedends Set this option to TRUE or FALSE to activate or deactivate
#' 	the paired-end mode. Only read pairs where both ends
#' 	are mapped will be considered and assigned to the bin where the 
#' 	midpoint of the read pair is located. 
#' 	If this flag is set, the 'shift' option will be ignored. 
#' @param shifts Shift the reads in the 5' direction by a fixed number
#'	of base pairs. The read will be assigned to the bin where the 
#'	shifted 5' end of the read is located. This option is ignored when
#'  the paired-end mode is active. (default: 75)
#' @details The \code{mapqs, pairedends, shifts} parameters can be provided
#' 	once and they will be used for all bam files, or as many times as there
#' 	are bam files, in which case the i-th path on the list will be matched 
#' 	to the i-th element in the vector.
#' @return A count matrix where the rows are the different marks
#' 	labelled with the names of the list \code{marks}, and the columns
#' 	are the different bins in the regions.
#' @export
getcounts <- function(regions, marks, binsize=200, mapqs=0, pairedends=FALSE, shifts=75){
	#ARGUMENT CHECKING
	if (!inherits(regions, "GRanges")) stop("regions must be a GRanges object")
	if (length(binsize) != 1 || binsize <= 0) stop("invalid binsize")
	if (any(width(regions) %% binsize != 0)) stop("region widths must be multiples of binsize")
	if (!is.list(marks)) stop("marks must be a list")
	if (is.null(names(marks))) stop("the 'marks' list must be named")
	if (any(!sapply(marks, file.exists))) stop("invalid paths provided (file doesn't exist)")
	if (!is.logical(pairedends)) stop("'pairedends' must be of type 'logical'")
	
	nmarks <- length(marks)
	if (nmarks==0) stop("no marks provided")
	mapqs <- fixLength(mapqs, nmarks)
	pairedends <- fixLength(pairedends, nmarks)
	shifts <- fixLength(shifts, nmarks)
	
	
	cat("getting counts\n")
	countsList <- lapply(1:nmarks, function(i){
		path <- marks[[i]]
		mark <- names(marks)[i]
		mapq <- mapqs[i]
		paired.end <- pairedends[i]
		shift <- shifts[i]
		if (paired.end) shift <- 0
		paired.end <- ifelse(paired.end, "midpoint", "ignore")
		bprof <- bamProfile(path, regions, paired.end=paired.end, binsize=binsize, shift=shift)
		unlist(as.list(bprof))
	})
	
	cat("pileup done, storing everything in a matrix\n")
	counts <- t(bindCols(countsList))
	rownames(counts) <- names(marks)
	counts
}

fixLength <- function(v, len){
	if (length(v)==1) return(rep(v, len))
	if (length(v)==len) return(v)
	stop(paste0("cannot interpret the given vector as a vector of length ", len,
	":\n provide either ", len, " elements or just 1"))
}

#' Refine regions
#'
#' Refine a GRanges object to make sure that it is compatible with a 
#' binning scheme of a given binsize. There is more than one way of doing it.
#' In the way it is done here, the start coordinates and the end coordinates
#' of the provided regions will become respectively the next number of the form 
#' \code{binsize*k + 1} and the previous number of the form \code{binsize*k},
#' so that the refined regions will always be contained in the original ones.
#' @param regions GRanges object containing the genomic regions of interest.
#' @param binsize The size of each bin in basepairs.
#' @return A GRanges object with the refined regions.
#' @export
refineRegions <- function(regions, binsize){
	starts <- start(regions)-1
	ends <- end(regions)
	newstarts <- binsize*(ceiling(starts/binsize))
	newends <- binsize*(floor(ends/binsize))
	valid <- newends > newstarts
	regions <- regions[valid]
	newstarts <- newstarts[valid]
	newends <- newends[valid]
	start(regions) <- newstarts+1
	end(regions) <- newends
	
	regions
}

readMarks <- function(marksArg){
	marks <- list()
	for (m in strsplit(marksArg, split=":")){
		if (length(m) != 2) stop("invalid mark specification")
		m[1] <- sanitizeFilename(m[1])
		marks[[m[1]]] <- m[2]
	}
	marks
}

validatePath <- function(path){
	if (!file.create(path)) stop("invalid path specified")
	path
}

