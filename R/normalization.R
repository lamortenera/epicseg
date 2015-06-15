validateVmat <- function(vmat){
    if (!is.matrix(vmat)) stop("vmat must be a matrix")
    if (!is.numeric(vmat)) stop("'vmat' must be numeric")
    if (length(vmat) <= 0) stop("'vmat' must be non-empty")
}

quantileNormalization <- function(vmat, ref=c("median", "min", "mean"), nthreads=1){
    validateVmat(vmat)
    if (!is.numeric(ref)){
        ref <- match.arg(ref)
        ref <- getRef(vmat, ref, nthreads)
    }
    if (length(ref) != nrow(vmat)) stop("reference vector has the wrong length")
    
    quantileNorm(vmat, ref, nthreads=nthreads)
}

defaultSFFun <- function(vmat, method="RLE", ...){
    edgeR::calcNormFactors(vmat, method=method, ...)
}

linearNormalization <- function(vmat, sfFun=defaultSFFun, ...){
    validateVmat(vmat)
    #compute scaling factors
    sf <- sfFun(vmat, ...)
    
    #scale vectors
    res <- round(vmat*sf[col(vmat)])
    storage.mode(res) <- "integer"
    res
}


defaultRepFun <- function(vmat, normFun=quantileNormalization, nthreads=1, ...){
    normFunOpts <- list(...)
    normFunOpts$vmat <- vmat
    if ("nthreads" %in% names(formals(normFun))) {
        normFunOpts[["nthreads"]] <- nthreads
    }
    nvmat <- do.call(normFun, normFunOpts)
    colSummary(t(nvmat), "median", nthreads=nthreads)
}
