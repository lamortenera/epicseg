context("Normalization routines")

source("utils.R")

test_that("Normalization routines",{
    expect_equal(testSortCounts(c(4,3,2,1)), c(1,2,3,4))
    expect_equal(testSortCounts(100*c(4,3,2,1)), 100*c(1,2,3,4))
    v <- sample(10000)
    expect_equal(testSortCounts(v), sort(v))
    v <- sample(10000, size=100)
    expect_equal(testSortCounts(v), sort(v))
    
    expect_equal(testMeanAndMedian(c(1,2), "mean"), 1)
    expect_equal(testMeanAndMedian(c(1,2), "median"), 1)
    expect_equal(testMeanAndMedian(c(1,101), "mean"), 51)
    expect_equal(testMeanAndMedian(c(1,2,101), "median"), 2)
    
    
    Rmean <- function(v){ as.integer(mean(v)) }
    Rmedian <- function(v) { as.integer(median(v)) }
    Rmin <- min
    fn2f <- list("median"=Rmedian, "mean"=Rmean, "min"=min)
    
    for (i in 1:10){
        v <- sample(10000+i, 100)
        expect_equal(testMeanAndMedian(v,"mean"), Rmean(v))
        expect_equal(testMeanAndMedian(v,"median"), Rmedian(v))
    }
    
    RgetRef <- function(mat, f){
        smat <- apply(mat, 2, sort)
        apply(smat, 1, f)
    }
    
    for (fn in names(fn2f)){
        expect_equal(getRef(matrix(c(1,2,3,4,4,3,2,1), ncol=2), fn), c(1,2,3,4))
    }
    
    for (i in 1:10){
        mat <- matrix(sample(10000), ncol=100)
        for (fn in names(fn2f)){
            expect_equal(getRef(mat, fn), RgetRef(mat, fn2f[[fn]]), label=fn)
            expect_equal(getRef(mat, fn), getRef(mat, fn, nthreads=4))
        }
    }
    
    #this is guaranteed to work only if the numbers in each column are unique
    RquantileNorm <- function(mat, ref){
        apply(mat, 2, function(v){
            o <- order(v)
            ans <- integer(length(v))
            ans[o] <- ref
            ans
        })
    }
    
    expect_equal(quantileNorm(matrix(c(1,4,8,10,4,3,2,1), ncol=2), c(1,3,9,27)), 
        matrix(c(1,3,9,27,27,9,3,1), ncol=2))
    
    #each count is unique
    mat <- matrix(sample(10000), ncol=100)
    ref <- sort(sample(10000,100))
    expect_equal(quantileNorm(mat, ref), RquantileNorm(mat, ref))
    
    
    #repeated counts
    mat <- matrix(sample(100, size=10000, replace=T), ncol=100)
    ref <- sort(sample(10000,100))
    expect_equal(quantileNorm(mat, ref), quantileNorm(mat, ref, nthreads=3))
    #make sure that each column has a different seed for the shuffles
    mat[,c(1,2)] <- 0
    nmat <- quantileNorm(mat, ref)
    expect_true(any(nmat[,1] != nmat[,2]))
    
    #quantileNormalization <- function(vmat, ref=c("median", "min", "mean"))
    mat <- cbind(c(0,2,4,8), c(3,0,6,9), c(0,7,21,14))
    nmat <- quantileNormalization(mat, "min")
    expect_equal(nmat, cbind(c(0,2,4,8), c(2,0,4,8), c(0,2,8,4)))
    nmat <- quantileNormalization(mat, "median")
    expect_equal(nmat, cbind(c(0,3,6,9), c(3,0,6,9), c(0,3,9,6)))
    nmat <- quantileNormalization(mat, "mean")
    expect_equal(nmat, cbind(c(0,4,8,12), c(4,0,8,12), c(0,4,12,8)))
    
    #linearNormalization <- function(vmat, what=c("coverage"), ref=c("median", "min", "mean"))
    #I'm just going to test that it doesn't produce errors
    qfun <- function(p) { function(v) quantile(v, p) }
    for (stat in c(qfun(.5), qfun(1), sum)){
        for (ref in list("median", "min", "mean", 10)){
            linearNormalization(mat, stat=stat, ref=ref)
        }
    }
    
    #test that aggregating 3 tracks is what you expect
    binsize <- 157
    regions <- readRegions(shortcuts$regs.bed)
    regions <- refineRegions(regions, binsize=binsize)
    labs <- c("H3K4me3", "H3K36me3", "H3K9me3")
    paths <- unlist(shortcuts[c("bam1", "bam2", "bam3")])
    mp <- paste(sep=":", c(labs, rep("mixed", 3)), rep(paths,2))
            
    bamtab <- makeBamtab(mp)
    counts <- suppressMessages(getcounts(regions, bamtab, binsize=binsize))
    expect_equal(defaultRepFun(t(counts[labs,])), counts["mixed",])
    
    
    #test clist2mlist and mlist2clist
    nmats <- 4
    nc <- 500
    nr <- 500
    clist <- lapply(1:nmats, function(r) matrix(rpois(nr*nc, lambda=500), nrow=nr))
    R_clist2mlist <- function(clist){
        mlist <- aperm(simplify2array(clist), c(2,3,1))
        lapply(1:dim(mlist)[3], function(i) mlist[,,i])
    }
    expect_equal(R_clist2mlist(clist), clist2mlist(clist))
    for (i in seq_along(clist)) rownames(clist[[i]]) <- as.character(1:nr)
    expect_equal(mlist2clist(clist2mlist(clist)), clist)
    for (i in seq_along(clist)) dimnames(clist[[i]]) <- NULL
    for (nthreads in c(1,3)){
        expect_equal(mlist2clist(clist2mlist(clist, nthreads=nthreads), nthreads=nthreads), clist)
    }
    
    
    
    #test that normalizing the counts is what you expect
    os <- matrix(cbind(c(1,2,3), c(2,1,3), c(3,2,1)))
    
    clist <- lapply(1:ncol(os), function(i){
        bamtab <- makeBamtab(paste(sep=":", labs, paths[os[,i]]))
        suppressMessages(getcounts(regions, bamtab, binsize=binsize))
    })
    clist2 <- normalizecounts(clist)
    
    for (lab in labs){
        oldvmat <- sapply(clist, function(counts) counts[lab,])
        newvmat <- sapply(clist2, function(counts) counts[lab,])
        expect_equal(quantileNormalization(oldvmat), newvmat)
    }
    
    
})
