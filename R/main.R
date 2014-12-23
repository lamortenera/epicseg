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
#' @import RColorBrewer
#' @import Rcpp
#' @import IRanges
#' @import GenomicRanges
#' @import kfoots
#' @import bamsignals
NULL

#writes a file that can be used as a launcher
#to use epicseg through the 
#command line interface. The parameter 'dest' should end in '.R'
getLauncher <- function(dest){
	writeLines(c(
	"cat('loading required libraries\\n')",
	"library(epicseg)",
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
CLI <- function(args, prog){
	CLIsubprograms <- list(
	segment=list(desc="Produce a segmentation and a report", fun=segmentCLI),
	report=list(desc="Produce a report for a given segmentation", fun=reportCLI),
	getcounts=list(desc="Produce a counts matrix from several bam files", fun=getcountsCLI))
	
	if (args[1] %in% names(CLIsubprograms)){
		CLIsubprograms[[args[1]]]$fun(args[2:length(args)], paste(prog, args[1]))
	} else {
		print_CLI(prog, CLIsubprograms)
		quit(status = 1)
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

file_ext <- function(x){ 
	pos <- regexpr("\\.([[:alnum:]]+)$", x)
	ifelse(pos > -1L, substring(x, pos + 1L), "")
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

