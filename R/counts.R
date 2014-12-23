validateCounts <- function(counts){
	if (is.null(counts) || !(is.matrix(counts))) stop("'counts' must be a matrix")
	if (is.null(rownames(counts))) stop("rows of the 'counts' matrix must be named")
	if (anyDuplicated(rownames(counts))) stop("row names of the 'counts' matrix must be unique")
	
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
