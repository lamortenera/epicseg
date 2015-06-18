defSuffix <- "_norm"
triggerOverwrite <- "-"
getNormalizeCountsOptions <- function(){
    list(
    list(arg="--counts", type="character", required=TRUE, vectorial=TRUE, 
    help="Paths to the count matrices such as those produced by `getcounts`"),
    list(arg="--suffix", type="character", parser=sanitizeFilename,
    default=defSuffix,
    help=paste0("Each normalized count matrix will be written in a filename 
    similar to that of the unnormalized matrix. A suffix will be added 
    before the file extension, For example,
    input: `path/to/counts.Rda`, output: `path/to/counts", defSuffix, ".Rda`.
    Set this option to change this behaviour. To overwrite the orignal matrices,
    give the argument `", triggerOverwrite, "`.")),
    list(arg="--nthreads", type="integer", help="Number of threads to use"))
}

normalizecountsCLI <- function(args, prog){
    opt <- parseArgs(getNormalizeCountsOptions(), args, prog)
    #get the real suffix
    if (is.null(opt$suffix)) opt$suffix <- defSuffix
    if (opt$suffix=="-") opt$suffix <- ""
    opt$suffix <- sanitizeFilename(opt$suffix)
    #get the new filenames
    targets <- sapply(opt$counts, function(src){
        fext <- file_ext(src, withDot=TRUE)
        res <- gsub(paste0(fext, "$"), paste0(opt$suffix, fext), src)
        checkWritable(res)
        res
    })
    message("reading matrices\n")
    clist <- lapply(opt$counts, readCounts)
    #check if the normalization function accepts the argument nthreads
    normFun <- formals(normalizecounts)[["normFun"]]
    writeNthreads <- "nthreads" %in% names(formals(normFun)) && 
                     "nthreads" %in% names(opt)
    message("normalizing matrices\n")
    if (writeNthreads) {
        clist2 <- normalizecounts(clist, nthreads=opt$nthreads)
    } else {
        clist2 <- normalizecounts(clist)
    }
    message("writing matrices\n")
    for (i in seq_along(clist2)) writeCounts(clist2[[i]], targets[i])
}

#' Normalize a list of count matrices
#'
#' Given a list of count matrices, corresponding rows across the different count
#' matrices are made comparable with a normalization procedure. 
#' Rows correspond if they have the same row name in all count matrices.
#' @param clist A list of count matrices. Each matrix must have the same set of
#'     row names and the same number of columns.
#' @param normFun A normalization function to use. Such a function takes a
#'     matrix as input argument and returns a matrix with the same dimensions
#'     as output. The goal of this function is to make the columns comparable.
#'     The functions \code{quantileNormalization} and \code{linearNormalization}
#'     are two already implemented, customizable functions.
#' @param nthreads The maximum number of threads. Parallelization is done on the
#'     histone marks using \code{mclapply}. If also \code{normFun} accepts the 
#'     \code{nthreads} argument, parallelization will be nested, with the outer 
#'     level using \code{min(nthreads, nmarks)} threads, and the inner level 
#'     using \code{floor(nthreads/nthreadsOuter)} threads.
#' @param ... options passed to normFun. 
#' @return A list of count matrices similar to the input list, but with the
#'     rows of each matrix normalized.
#' @export
normalizecounts <- function(clist, normFun=quantileNormalization, nthreads=1, ...){
    clist <- validateCList(clist)
    if (length(clist)==1) return(clist)
    #convert clist in a list of matrices for each mark with dims [bin, dataset]
    mlist <- clist2mlist(clist, nthreads=nthreads)
    #if the list was named, set the column names of each element in mlist
    if (!is.null(names(clist))){
        dnames <- list(NULL, names(clist))
        for (i in seq_along(mlist)) setDimnames_unsafe(mlist[[i]], dnames)
    }
    #apply the function to each mark
    nthreadsOuter <- min(nthreads, length(mlist))
    nthreadsInner <- floor(nthreads/nthreadsOuter)
    if ("nthreads" %in% names(formals(normFun))){
        mlist <- safe_mclapply(mc.cores=nthreadsOuter, mlist, normFun, 
            nthreads=nthreadsInner, ...)
    } else mlist <- safe_mclapply(mc.cores=nthreadsOuter, mlist, normFun, ...)
    #convert back to a clist
    setNames(mlist2clist(mlist, nthreads=nthreads), names(clist))
}
