context("R interface")

test_that("multiple datasets works",{
    #make a list of count matrices
    nmats <- 4; nc <- 500; nr <- 5
    clist <- lapply(1:nmats, function(r) matrix(rpois(nr*nc, lambda=500), nrow=nr))
    marks <- paste0("mark", 1:nr)
    for (i in seq_along(clist)) rownames(clist[[i]]) <- marks
    dsetNames <- paste0("dataset",1:nmats)
    names(clist) <- dsetNames

    #normalize counts
    clist <- normalizecounts(clist)
    
    #make some matching regions
    binsize <- 200
    gr <- GRanges(seqnames="chr1", IRanges(start=1, width=binsize*nc))

    #segment
    s <- suppressMessages(segment(clist, gr, 5, verbose=F))
    
    expect_equal(names(s$segments), dsetNames)
    expect_equal(dimnames(s$posteriors)[[3]], dsetNames)
    expect_equal(dimnames(s$states)[[2]], dsetNames)
    expect_equal(dimnames(s$viterbi)[[2]], dsetNames)
    expect_equal(s$model$marks, marks)
    
    report(s$segments, s$model, outdir=tempdir(), prefix="test")}
    
)
