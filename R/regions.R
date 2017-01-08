#just makes sure that the second and third fields are numbers
getBEDFormat <- function(bedline) {
    lineArgs <- strsplit(bedline, "\t")[[1]]
    if (length(lineArgs)<3) return(NULL)
    #I couldn't figure out any other way of checking if 
    #a field codes an integer
    matches <- grep("^[[:digit:]]+$", lineArgs[2:3])
    if (length(matches) < 2) return(NULL)
    length(lineArgs)
}

#reads a bed file up to the 6th field, 
#fields after the 6th are discarded
#also the score field (the 5th) is discarded
readRegions <- function(path) {
    stopifnot(length(path)==1)
    #tries to guess the right number of columns and header lines to be skipped
    #should adapt the function bed2GR so that it can parse the metaData
    nLines <- 5
    firstLines <- readLines(path, nLines)
    nFields <- NULL
    for (i in seq_along(firstLines)){
        nFields <- getBEDFormat(firstLines[i])
        if (!is.null(nFields)) break
    }
    if (is.null(nFields)) stop("unable to find a proper bed line in the file")
    #i stores the index of the first proper bed line
    #nFields stores the number of fields
    what <- list(character(), integer(), integer(), character(), NULL, character())[1:(min(nFields, 6))]
    if (nFields > 6){
        for (j in 7:nFields){
            what[[j]] <- NULL
        }
    }
    regions <- scan(path, what, sep="\t", skip=(i-1), flush=TRUE, quiet=TRUE)
    
    gr <- GRanges(seqnames=regions[[1]], IRanges(start=regions[[2]]+1, end=regions[[3]]))
    #set the name of the regions
    if (nFields >= 4) names(gr) <- regions[[4]]
    #set the strand of the regions
    if (nFields >= 6) {
        if (any(!(regions[[6]] %in% c("*",".","+","-")))) {
            warning("invalid strand specification in bed file: ignoring strand")
        } else strand(gr) <- sub("\\.", "*", regions[[6]])
    }

    gr
}

int2str <- function(v) format(v, trim=TRUE, scientific=FALSE)

#writes a bed file up to the 6th field, 
#metadata elements are discarded
#if the strand is present, the score field will be set to 0,
#and if the regions are not named, they will be all named as *
writeRegions <- function(gr, path){
    stopifnot(length(path)==1)
    if (!inherits(gr, "GRanges")) stop("expecting a GRanges object")
    tab <- data.frame(
    chr=      as.character(seqnames(gr)),
    start=    start(gr)-1,
    end=      end(gr))
    
    tab[,2] <- int2str(tab[,2])
    tab[,3] <- int2str(tab[,3])
    
    if (any(strand(gr)!="*")){
        nms <- names(gr)
        if (is.null(nms)) nms <- rep("*", length(gr))
        scores <- rep(0, length(gr))
        tab$name <- nms
        tab$score <- scores
        tab$strand <- as.character(strand(gr))
    } else if (!is.null(names(gr))) tab$name <- names(gr)
    
    write.table(tab, file=path, sep="\t", quote=F, row.names=F, col.names=F)
}
