validateCounts <- function(counts){
    if (is.null(counts) || !(is.matrix(counts))) stop("'counts' must be a matrix")
    if (is.null(rownames(counts))) stop("rows of the 'counts' matrix must be named")
    if (anyDuplicated(rownames(counts))) stop("row names of the 'counts' matrix must be unique")
    
}

validateCList <- function(clist){
    sapply(clist, validateCounts)
    if (length(clist) < 2) return(clist)
    marks <- rownames(clist[[1]])
    ncols <- ncol(clist[[1]])
    #make count matrices compatible (i-th row is always the i-th mark)
    for (i in 2:length(clist)){
        newmarks <- rownames(clist[[i]])
        #check compatibility
        if (!setequal(marks, newmarks) || ncol(clist[[i]]) != ncols) {
            stop("non-compatible count matrices") }
        if (any(marks != newmarks)){
            #rearrange rows
            clist[[i]] <- clist[[i]][marks,]
        }
    }
    clist
}

readCounts <- function(path){
    if (isRdata(path)){
        counts <- readRdata(path)
    } else {
        countsList <- read.table(path, header=TRUE, sep="\t", colClasses="integer", row.names=NULL, check.names=FALSE)
        counts <- t(bindCols(countsList))
    }
    validateCounts(counts)
    counts
}



writeCounts <- function(counts, path){
    validateCounts(counts)
    if (isRdata(path)){
        save(counts, file=path)
    } else {
        writeCountsTXT(counts, rownames(counts), path)
    }
}
