#do the opposite of "expect_error"
runs <- function(expr){
    res <- try(force(expr), TRUE)
    no_error <- !inherits(res, "try-error")
    if (no_error) {
        return(expectation(TRUE, "code generated an error", 
        "code did not generate an error"))
    }
    else {
        expectation(FALSE, paste0("threw an error:\n", res[1]), "no error thrown")
    }
}
expect_runs <- function(object, info = NULL, label = NULL){
    if (is.null(label)) {
        label <- testthat:::find_expr("object")
    }
    expect_that(object, runs, info = info, label = label)
}

tmpdir <- file.path(tempdir(), "epicseg_test")
dir.create(tmpdir, showW=F)
#these are shortcuts for the command lines
shortcuts <- list(
    init="Rscript epicseg.R",
    tmpdir=tmpdir,
    bam1=system.file("extdata","H3K4me3.bam", package="epicseg"),
    bam2=system.file("extdata","H3K36me3.bam", package="epicseg"),
    bam3=system.file("extdata","H3K9me3.bam", package="epicseg"),
    bed1=system.file("extdata","genes.bed", package="epicseg"),
    regs.bed=system.file("extdata","contigs.bed", package="epicseg"))
