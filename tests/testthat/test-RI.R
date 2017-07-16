source("utils.R")
context("R interface")

make_clist <- function(nmats, nr, nc){
    clist <- lapply(1:nmats, function(r) matrix(rpois(nr*nc, lambda=500), nrow=nr))
    marks <- paste0("mark", 1:nr)
    for (i in seq_along(clist)) rownames(clist[[i]]) <- marks
    dsetNames <- paste0("dataset",1:nmats)
    names(clist) <- dsetNames
    clist
}

test_that("can bind matrices",{
    # originally this test was designed to make sure
    # that large matrices can be bound, 
    # where the number of elements of the resulting
    # matrix is more than a int32 can represent. Unfortunately
    # I am unable to test that (I don't have the appropriate hardware)
    # so this test is only checking if bindClist works
    # with normal matrices
    nmats <- 100; nc <- 100; nr <- 5
    clist <- make_clist(nmats, nr, nc)
    mat1 <- bindClist(clist, nthreads=4)
    mat2 <- do.call(cbind, clist)
    expect_equal(mat1, mat2) 
})

test_that("multiple datasets works",{
    #make a list of count matrices
    nmats <- 4; nc <- 500; nr <- 5
    clist <- make_clist(nmats, nr, nc)
    dsetNames <- names(clist)
    marks <- rownames(clist[[1]])
    #normalize counts
    clist <- normalizecounts(clist)
    for (method in c("TMM", "RLE")){
        normalizecounts(clist, epicseg:::linearNormalization, method=method)
    }
    
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
    
    report(s$segments, s$model, outdir=tempdir(), prefix="test")
})


test_that("bins as regions throws error",{
    #make a list of count matrices
    nmats <- 4; nc <- 500; nr <- 5
    clist <- make_clist(nmats, nr, nc)

    #make some matching regions
    binsize <- 200
    gr_good <- GRanges(seqnames="chr1", IRanges(start=1, width=binsize*nc))
    gr_bad <- GRanges(seqnames="chr1", IRanges(start=1 + binsize*(0:(nc-1)), width=binsize))

    #everything should work
    suppressMessages(segment(clist, gr_good, 5, verbose=F))
    suppressMessages(segment(clist, gr_good, 5, verbose=F))
    #everything should fail
    expect_error(suppressMessages(segment(clist[[1]], gr_bad, 5, verbose=F)),
        regexp="no region contains consecutive bins")
    expect_error(suppressMessages(segment(clist, gr_bad, 5, verbose=F)),
        regexp="no region contains consecutive bins")
})


test_that("kfoots error handler", {
    # test getBin function
    gr = GRanges(seqnames=c(1,    1,    2,   3), 
            IRanges(start=c(200,  800,  400, 600),
                    end  =c(1000, 1000, 800, 800)-1))
    binsize <- 200
    expect_equal(sum(width(gr)) %% binsize, 0)
    bins <- unlist(tile(gr, width=binsize))
    nbins <- length(bins)
    mybins <- do.call(c, sapply(1:nbins, getBin, regions=gr, binsize=binsize))
    mybins2 <- do.call(c, sapply(nbins + (1:nbins), getBin, regions=gr, binsize=binsize))
    
    f <- function(gr) paste0(seqnames(gr), start(gr), end(gr), sep="@")
    expect_equal(f(mybins), f(bins))
    expect_equal(f(mybins2), f(bins))
    
    
    # make sure an underflow returns an exception with the word 'underflow'
    uflow <- get(load(system.file("extdata/uflow_minimal.Rdata", package="epicseg")))
    model <- list(emisP=uflow$models, transP=uflow$trans, initP=uflow$initP, 
                  marks=rownames(uflow$counts))
    counts <- uflow$counts
    regions <- GRanges(seqnames="foo", IRanges(start=200, width=400))
    expect_error(segment(counts, regions, model=model, verbose=F),
                 regexp="underflow", ignore.case=T)
    # try case of a count list
    clist <- list(dset1=counts, dset2=counts)
    model2 <- model
    nstates <- length(model$emisP)
    unif_init_p <- rep(1/nstates, nstates)
    model2$initP <- cbind(unif_init_p, model$initP)
    expect_error(segment(clist, regions, model=model2, verbose=F),
                 regexp="underflow", ignore.case=T)
})

test_that("avgStateProfile", {
	genes <- readRegions(system.file("extdata/genes.bed", package="epicseg"))
	segms <- readRegions(system.file("extdata/5s_segmentation.bed", package="epicseg"))
	expect_runs(avgStateProfile(genes, segms, 5))
	mid_genes <- genes
	start(genes) <- 0.5*(start(genes) + end(genes))
	end(genes) <- start(genes)
	expect_runs(avgStateProfile(genes, segms, 5))
})
