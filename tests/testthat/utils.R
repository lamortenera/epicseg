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


#I would be amazed if this worked also in windows
#I also expect that calling this function twice is the same as calling it once
configureSys <- function(){
    #create a temporary directory where all output will be stored
    tmpdir <- file.path(tempdir(), "epicseg_test")
    dir.create(tmpdir, showW=F)
    #make a launcher for epicseg
    launcherPath <- file.path(tmpdir, "epicseg.R")
    epicseg:::getLauncher(launcherPath)
    #grant write permissions
    Sys.chmod(launcherPath, "700")
    #add tmpdir to the PATH variable
    PATHvar <- Sys.getenv("PATH")
    #but first check that it is not already there
    paths <- strsplit(PATHvar, ":")[[1]]
    if (!tmpdir %in% paths) PATHvar <- paste(tmpdir, PATHvar, sep=":")
    #these are environment variables used in the command line
    shortcuts <- list(
        PATH=PATHvar,
        outdir=tmpdir,
        indir =system.file("extdata", package="epicseg"))
    #set the environment variables
    for (name in names(shortcuts)){
        #I didn't find a cleaner way of doing this...
        Rcmd <- paste0("Sys.setenv(\"", name, "\"=\"", shortcuts[[name]], "\")")
        eval(parse(text=Rcmd))
    }
}
