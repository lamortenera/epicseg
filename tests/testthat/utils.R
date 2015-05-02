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

tmpdir <- tempdir()
#these are shortcuts for the command lines
shortcuts <- list(
    init="Rscript epicseg.R",
    tmpdir=tempdir(),
    bam1="/project/epigenome/IMR90/H3K4me3.bam",
    bam2="/project/epigenome/IMR90/H3K36me3.bam",
    bam3="/project/epigenome/IMR90/H3K27me3.bam",
    bed1="/project/ale/home/data/kfoots_paper/data_old/regions/knownGenes.bed",
    bed2="/project/ale/home/data/kfoots_paper/data_old/mansegm/IMR90/mansegm.bed",
    tests="/project/ale/home/data/epicseg_pkg/tests",
    regs.bed="/project/ale/home/data/epicseg_pkg/tests/data/regions.bed")
