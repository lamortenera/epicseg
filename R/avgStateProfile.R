avgStateProfile <- function(genes, segm, nstates, before=200, after=200, reflen=round(median(width(genes)))){
    #dimensions of the final matrix
    sfact <- factor(names(segm), levels=1:nstates)
    if (any(is.na(sfact))) stop("names of the GRanges object must be of the kind '1','2', ..'nclasses'")
    npos <- reflen + before + after
    totlen <- nstates*npos
    
    #extend the genes
    gstr <- strand(genes)
    extgenes <- genes
    start(extgenes[gstr!="-"]) <- start(genes[gstr!="-"]) - before
    end(extgenes[gstr!="-"]) <- end(genes[gstr!="-"]) + after
    start(extgenes[gstr=="-"]) <- start(genes[gstr=="-"]) - after
    end(extgenes[gstr=="-"]) <- end(genes[gstr=="-"]) + before
    #segments should be unstranded
    strand(segm) <- "*"
    
    #find all the overlaps
    ov <- findOverlaps(extgenes, segm, type="any", select="all")
    gidx <- queryHits(ov)
    sidx <- subjectHits(ov)
    
    #find start and end index of each overlap
    gcomb <- extgenes[gidx]
    scomb <- segm[sidx]
    posref <- as.logical(strand(gcomb)!="-")
    ovstart <- integer(length(ov))
    ovend <- integer(length(ov))
    ovstart[posref] <- start(scomb)[posref] - start(gcomb)[posref]
    ovstart[!posref] <- end(gcomb)[!posref] - end(scomb)[!posref]
    ovend[posref] <- end(gcomb)[posref] - end(scomb)[posref]
    ovend[!posref] <- start(scomb)[!posref] - start(gcomb)[!posref]
    
    #rescale the coordinates
    gwidth <- width(genes)[gidx]
    int <- ovstart - before
    int[int < 0] <- 0
    int[int > gwidth] <- gwidth[int > gwidth]
    ovstart <- ovstart + round(int*reflen/gwidth) - int
    if (any(ovstart>=npos)) stop("something went wrong in rescaling the coordinates")
    #same with the ends
    int <- ovend - after
    int[int < 0] <- 0
    int[int > gwidth] <- gwidth[int > gwidth]
    ovend <- ovend + round(int*reflen/gwidth) - int
    if (any(ovend>=npos)) stop("something went wrong in rescaling the coordinates")
    
    #for each start position increase the count in the corresponding column of mat
    # the row is determined by the state of the segment
    #indices start at 0
    rows <- as.integer(sfact[sidx])-1
    rnames <- levels(sfact)
    cols <- ovstart; cols[cols<0] <- 0
    matidx <- rows + nstates*cols
    mat1 <- kfoots:::sumAt(rep(1, length(matidx)), matidx, totlen, zeroIdx=T)
    #for each end position decrease the count in the next column in mat
    ovend <- npos - ovend
    cols <- ovend
    matidx <- (rows + nstates*cols)[cols<npos]
    mat2 <- kfoots:::sumAt(rep(1, length(matidx)), matidx, totlen, zeroIdx=T)
    #final matrix
    mat <- matrix(mat1-mat2, nrow=nstates, ncol=npos)
    
    #do cumsum on each row
    mat <- t(apply(mat, 1, cumsum))
    rownames(mat) <- rnames
    colnames(mat) <- 1:npos - before
    
    #sort the row names
    mat <- mat[as.character(1:nstates),,drop=FALSE]
    
    mat
}

plotProfile <- function(mat, colors, legend.pos="top", xlab="offset", 
    ylab="number of regions", main="average state profile", xlim=range(xs), 
    ylim=c(0, max(mat)), legend=T, legend.bg=par("bg"), lwd=par("lwd"),...){
    if (is.null(colnames(mat))) stop("colnames(mat) must be set")
    xs <- as.integer(colnames(mat))
    spacer <- abs(xs[1])
    reflen <- ncol(mat) - 2*spacer
    if (reflen < 0) stop("unable to determine reflen")
    plot(NA, NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...)
    
    if (reflen > 1){
        abline(v=0, lty=2)
        abline(v=reflen, lty=2)
        mtext("start", side=3, at=0)
        mtext("end", side=3, at=reflen)
    }
    
    for (i in 1:nrow(mat)){
        lines(colnames(mat), mat[i,], col=colors[i], lwd=lwd)
    }
    
    if (legend) {
        legend(legend=gsub("_", " ", rownames(mat)), col=colors, lty=1, 
        x=legend.pos, bg=legend.bg, lwd=lwd)
    }
  
}

#clustmeans is a matrix where the columns are clusters and the rows are marks
#the rows must be named.
automaticColoring <- function(clustmeans, markToCol=list(red="H3K4me3", green4="H3K36me3", blue=c("H3K9me3", "H3K27me3"), gold="H3K4me1")){
    nstates <- ncol(clustmeans)
    if (is.null(rownames(clustmeans))) {
        stop("clustmeans matrix must have rownames set")}
    available_marks <- rownames(clustmeans)
    for (marks in markToCol){
        if (!any(marks %in% available_marks)) stop("cannot find matching marks")
    }
    
    #clusters that never occur should have NaNs in the matrix clustmeans
    validClust <- apply(clustmeans, 2, function(v) all(is.finite(v)))
    clustmeans <- clustmeans[,validClust, drop=F]
    n <- ncol(clustmeans)
    if (n < 3) stop("too few valid clusters...")
    
    anchcol <- names(markToCol)
    nanchors <- length(anchcol)
    anchors <- integer(nanchors)
    
    
    for (i in seq_along(markToCol)){
        marks <- markToCol[[i]]
        scores <- colSums(clustmeans[intersect(marks, available_marks), , drop=F])
        #make sure we don't pick a cluster that we already picked
        scores[anchors[anchors>0]] <- -1
        anchors[i] <- which.max(scores)
    }
    
    #compute distance matrix
    #dmat <- kfoots:::KL_dist_mat(clustmeans+1e-6, 1)
    dmat <- as.matrix(dist(t(log(clustmeans+1e-6))))
    #similarity matrix
    smat <- 1/dmat; smat[!is.finite(smat)] <- 1e200
    #weights
    wts <- smat[anchors,]
    #make the weight matrix sparse:
    #two non-zero entries for non-anchor clusters
    #one non-zero entry for anchor clusters
    for (i in 1:n){
        o <- order(wts[,i], decreasing=T)
        if (i %in% anchors){
            wts[o[2:n],i] <- 0
        } else {
            wts[o[3:n],i] <- 0
        }
    }
    #make sure the wts sum up to 1
    wts <- apply(wts, 2, function(v) v/sum(v))
    
    #anchor rgb colors
    anchrgb <- t(col2rgb(anchcol))
    #valid cluster colors
    validrgb <- rgb(round(t(apply(wts, 2, function(wt) colSums(wt*anchrgb)))), maxColorValue=255)
    #colors for cluster that don't occur (this shouldn't matter)
    nonvalidrgb <- rgb(t(col2rgb("white")), alpha=0, maxColorValue=255)
    
    res <- rep(nonvalidrgb, nstates)
    res[validClust] <- validrgb
    
    res
}

#converts R colors to colors for the bed format
col2bedCol <- function(colors){
    rgbs <- col2rgb(colors)
    apply(rgbs, 2, paste, collapse=",")
}

