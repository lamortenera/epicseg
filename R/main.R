#' Chromatin segmentation
#'
#' The package provides methods for finding chromatin states in 
#' epigenomic data. There is an R interface as well as a 
#' command line interface.
#'
#' @name epicseg
#' @docType package
#' @author Alessandro Mammana \email{mammana@@molgen.mpg.de}
#' @useDynLib epicseg
#' @import parallel
#' @import RColorBrewer
#' @import Rcpp
#' @import IRanges
#' @import GenomicRanges
#' @import kfoots
#' @import bamsignals
#' @import edgeR
NULL


getLauncher <- function(dest="epicseg.R"){
    RscriptPath <- file.path(Sys.getenv("R_HOME"), "bin", "Rscript")
    shebang <- NULL
    if (!file.exists(RscriptPath)) {
        warning("Rscript executable not found at the expected location")
    } else shebang <- paste0("#!", RscriptPath)
    
    epicsegPath <- path.package("epicseg")
    #if the package was loaded with devtools::load_all it is not properly
    #installed and we need a slightly different launcher
    if (tryCatch("epicseg" %in% devtools::dev_packages(), error=function(e) FALSE)){
        devtoolsPath <- find.package("devtools")
        loadDevtools <- paste0("library(devtools, lib.loc=\"", dirname(devtoolsPath), "\")")
        loadEpicseg <- paste0("devtools::load_all(\"", epicsegPath, "\", quiet=TRUE)")
        loadLibs <- paste(sep="\n", loadDevtools, loadEpicseg)
    } else {
        loadLibs <- paste0("library(epicseg, lib.loc=\"", dirname(path.package("epicseg")), "\")")
    }
    
    writeLines( c(shebang, "cat('loading epicseg\\n')", loadLibs,
    "epicseg:::CLI(args=commandArgs(trailingOnly=TRUE), epicseg:::getProg())"),
    dest)
}

#CLI SUBPROGRAMS MANAGEMENT
print_CLI <- function(prog, subprograms){
    cat(sub("%prog", prog, "usage: %prog [subprogram] [subprogram options]"), fill=TRUE)
    cat("Choose one of the available subprograms", fill=TRUE)
    cat("\n")
    cat("Subprograms:", sep="\n")
    for (i in seq_along(subprograms)){
        progname = names(subprograms)[i]
        progdesc = subprograms[[i]]$desc
        cat("\t")
        cat(progname)
        cat("\n\t")
        cat(progdesc)
        cat("\n\n")
    }
}

#the CLI subprograms are:
#segmentCLI, in the file segment.R
#reportCLI, in the file report.R
#getcountsCLI, in the file getcounts.R
#normalizecountsCLI, in the file normalizecounts.R
getCLIsubprograms <- function(){list(
    getcounts=list(desc="Produce a counts matrix from several bam files", 
    fun=getcountsCLI, cliargs=getGetcountsOptions),
    normalizecounts=list(desc="Normalize several count matrices", 
    fun=normalizecountsCLI, cliargs=getNormalizeCountsOptions),
    segment=list(desc="Produce a segmentation and a report",
    fun=segmentCLI, cliargs=getSegmentOptions),
    report=list(desc="Produce a report for a given segmentation", 
    fun=reportCLI, cliargs=getReportOptions))}
    
CLI <- function(args, prog){
    CLIsubprograms <- getCLIsubprograms()
    if (args[1] %in% names(CLIsubprograms)){
        CLIsubprograms[[args[1]]]$fun(args[-1], paste(prog, args[1]))
    } else {
        print_CLI(prog, CLIsubprograms)
        quit(status=1)
    }
}


#copied and pasted from the getopt function get_Rscript_filename
getProg <- function(){
    prog <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
    if (.Platform$OS.type == "windows") {
        prog <- gsub("\\\\", "\\\\\\\\", prog)
    }
    prog
}


#MISC
#this is used only for testing itself
#the first two arguments are actually discarded
simulateCLI <- function(cmd){
    #remove newlines
    cmd <- gsub("\n", "", cmd)
    #split arguments
    args <- strsplit(cmd, "[ ]+")[[1]]
    if (length(args) < 2) stop("need at least 2 args")
    #redefining commandArgs called in the driver script
    commandArgs <- function(trailingOnly=FALSE){
        sargs <- c()
        if (length(args) > 2) sargs <- args[3:length(args)]
        if (!trailingOnly){
            sargs <- c(paste0("--file=",args[2]), sargs)
        }
        sargs
    }
    CLI(commandArgs(trailingOnly=TRUE), args[2])
}

file_ext <- function(x, withDot=FALSE){ 
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    if (!withDot) pos <- pos + 1
    ifelse(pos > 0L, substring(x, pos), "")
}

isRdata <- function(path){
    tolower(file_ext(path)) %in% c("rda", "rdata")
}

readRdata <- function(path){
    oname <- load(path)
    if (length(oname)!=1) stop("expecting one object in the saved R data")
    get(oname)
}

rename <- function(l, oldname, newname){
    tmp <- l[[oldname]]
    l[[oldname]] <- NULL
    l[[newname]] <- tmp
    l
}

label_sc_path <- function(l_p, unique.labels=FALSE){
    npaths <- length(l_p)
    lp <- strsplit(l_p, split=":")
    lens <- sapply(lp, length)
    if (any(lens != 2)) {
        stop(
        paste0("Invalid input. Expecting input of the form LABEL:PATH, found:\n",
         l_p[which(lens!=2)[1]]))}
    lp <- simplify2array(lp)
    if (unique.labels && anyDuplicated(lp[1,])) {
        ex <- lp[1,which(duplicated(lp[1,]))[1]]
        stop(paste0("The following label appeared more than once: ", ex))
    }
    data.frame(label=lp[1,], path=lp[2,], stringsAsFactors=F)
}

#make sure that errors occurring in mclapply get propagated
propagateErrors <- function(l){
    for (el in l) if (inherits(el, "try-error")) stop(el)
    l
}
safe_mclapply <- function(...) propagateErrors(mclapply(...))
