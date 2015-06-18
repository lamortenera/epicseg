getGetcountsOptions <- function(){list(
    list(arg="--regions", type="character", required=TRUE, parser=readRegions,
    help="Path to the BED file with the genomic regions of interest.
    These regions will be automatically partitioned into smaller, 
    consecutive bins. Only the first three fields in the file matter. 
    If the region lengths are not multiples of the given binsize
    a new bed file will be produced where each coordinate 
    is a multiple of binsize. Use this new file together with
    the count matrix for later analyses."),
    list(arg="--mark", type="character", required=TRUE, vectorial=TRUE, meta="label:path",
    help="Mark name and path to a bam file where to extract reads from.
    The bam files must be indexed and the chromosome names must match with 
    those specified in the bed file. Entries with the same mark name will
    be treated as replicates and collapsed into one experiment.
    This option must be repeated for each mark, for example:
    `-m H3K4me3:/path1/foo1.bam -m H3K36me3:/path2/foo2.bam`"),
    list(arg="--target", type="character", required=TRUE, parser=validatePath,
    help="Full path to the count matrix that will be produced as output.
    To save the matrix as an R archive, provide a path that ends
    in `.Rdata` or `.rda`, otherwise it will be saved as a text file.
    Using R archives, writing the file will be a bit slower,
    but reading it will be considerably faster."),
    list(arg="--binsize", type="integer", default=formals(getcounts)$binsize,
    help="Size of a bin in base pairs. Each given region will be partitioned into
    bins of this size."),
    list(arg="--mapq", type="integer", vectorial=TRUE, default=0,
    help="Minimum mapping quality for the reads (see the bam format
    specification for the mapq field). Only reads with the mapq field
    above or equal to the specified value will be considered. 
    Specify this option once for all marks, or specify it as many times 
    as there are marks."),
    list(arg="--pairedend", type="logical", vectorial=TRUE, meta="true_or_false", default=FALSE,
    help="Set this option to TRUE or FALSE to activate or deactivate
    the paired-end mode. Only read pairs where both ends
    are mapped will be considered and assigned to the bin where the 
    midpoint of the read pair is located. Specify this option 
    once for all marks, or specify it as many times as there are marks.
    If this flag is set, the `shift` option will be ignored."),
    list(arg="--shift", type="integer", vectorial=TRUE, default=75,
    help="Shift the reads in the 5' direction by a fixed number
    of base pairs. The read will be assigned to the bin where the 
    shifted 5' end of the read is located. Specify this option 
    once for all marks, or specify it as many times as there are marks.
    This option will be ignored in paired-end mode."),
    list(arg="--nthreads", type="integer", default=formals(getcounts)$nthreads,
    help="Number of threads to use."))}

getcountsCLI <- function(args, prog){
    opt <- parseArgs(getGetcountsOptions(), args, prog)
    target <- opt$target
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
    #make bamtab object
    bamtab <- makeBamtab(opt$mark, shift=opt$shift, mapq=opt$mapq, pairedend=opt$pairedend)
    #call getcounts
    counts <- getcounts(regions=opt$regions, bamtab=bamtab, binsize=opt$binsize, nthreads=opt$nthreads)
    #write results
    cat(sep="", "writing count matrix to the file '", target, "'\n")
    writeCounts(counts, target)
}

#' Make a count matrix
#'
#' Make a count matrix from a GRanges object and a list of bam files
#' @param regions GRanges object containing the genomic regions of interest.
#'     Each region will be divided into non-overlapping bins of equal size.
#'     The reads falling in each bin will be the elements of the final count matrix.
#' @param bamtab Data frame describing how to get the counts from each file. 
#'     The following columns are required: 'mark' (name of the histone mark),
#'     'path' (path to the bam file). The following columns are optional:
#'     'mapq' (minimum mapping quality), 'shift' (shift in the 3' to 5' direction
#'     applied to each read), 'pairedend' (paired end mode or not)
#' @param binsize The size of each bin in basepairs. Each region must have
#'     a width multiple of \code{binsize}, otherwise an error will be thrown.
#'     Use the function \code{refineRegions} to make sure that your GRanges object
#'     satisfies this constraint.
#' @param repFun In case there are replicate experiments (i.e. entries in
#'     bampath with the same 'mark' name), they will be collapsed into 
#'     one experiment using this function. This function takes a matrix of 
#'     counts as input where the columns are different experiments 
#'     and the rows are different bins and returns a vector with as many elements
#'     as the rows of the input matrix.
#' @param nthreads Number of threads. Parallelization is done on the histone
#'     marks using \code{mclapply}.
#' @return A count matrix where the rows are the different marks
#'     labelled with the names of the list \code{marks}, and the columns
#'     are the different bins in the regions.
#' @export
getcounts <- function(regions, bamtab, binsize=200, repFun=defaultRepFun, 
                                                                    nthreads=1){
    #ARGUMENT CHECKING
    if (!inherits(regions, "GRanges")) stop("regions must be a GRanges object")
    if (length(binsize) != 1 || binsize <= 0) stop("invalid binsize")
    if (any(width(regions) %% binsize != 0)) stop("region widths must be multiples of binsize")
    bamtab <- validateBamtab(bamtab)
    
    npaths <- nrow(bamtab)
    if (npaths==0) stop("no marks provided")
       
    message("getting counts")
    countsList <- safe_mclapply(mc.cores=nthreads, 1:npaths, function(i){
        b <- bamtab[i,]
        if (b$pairedend) b$shift <- 0
        b$pairedend <- ifelse(b$pairedend, "midpoint", "ignore")
        
        bprof <- bamProfile(b$path, regions, binsize=binsize, 
            shift=b$shift, mapq=b$mapq, paired.end=b$pairedend)
        unlist(as.list(bprof))
    })
    
    if (any(duplicated(bamtab$mark))){
        message("putting replicate experiments together")
        countsList <- lapply(setNames(nm=unique(bamtab$mark)), function(nm){
            subList <- countsList[bamtab$mark == nm]
            if (length(subList)==1) return( subList[[1]] )
            vmat <- bindCols(subList)
            repFun(vmat)
        })
    } else {
        names(countsList) <- bamtab$mark
    }
    
    message("pileup done, storing everything in a matrix")
    t(bindCols(countsList, nthreads=nthreads))
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



validatePath <- function(path){
    if (!file.create(path)) stop("invalid path specified")
    path
}

