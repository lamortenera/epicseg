getReportOptions <- function(){list(
    list(arg="--segments", type="character", required=TRUE, vectorial=TRUE,
    help="Path to the bed file containing the segmentation. The name field
    in the bed file must be a number from 1 to the maximum number of states.
    More than one path can be provided if the bed files describe alternative
    segmentations of the same genomic regions and if the paths are labelled, 
    i.e. of the form `LABEL:PATH`. "),
    list(arg="--model", type="character",  required=TRUE, parser=readModel,
    help="Path to the file with the parameters of the HMM.
    All fields must be present."),
    list(arg="--outdir", type="character", parser=validateOutdir, default=".",
    help="Path to the folder where all results will be saved. If
    the folder doesn't exist it will be created."),
    list(arg="--prefix", type="character", parser=sanitizeFilename, default="",
    help="Prefix to append to all produced files. The prefix should
    not be a path (use the -o option for that)."),
    list(arg="--colors", type="character", parser=readColors,
    help="Path to a text file containing one color per chromatin state.
    Each color must be on a separate line, and they can be specified
    either as RGB colors `R,G,B` where 0 <= R,G,B < 256, or as R colors
    (like `red`, `gold`, `green4`...)."),
    list(arg="--labels", type="character", parser=readLines,
    help="Path to a text file containing one short text label
    per chromatin state. Each label must be on a separate line."),
    list(arg="--annot", type="character", vectorial=TRUE, parser=readAnnotations, meta="label:path",
    help="Genomic annotation to compare with the segmentation. It must
    be specified with a short title (no spaces) and a path to the bed
    file containing the annotation (the strand matters), separated by `:`.
    This option can be repeated, for example:
    `--annot genes:/path1/genes.bed --annot cpg:/path2/cpgislands.bed`."))}

reportCLI <- function(args, prog){
    opt <- parseArgs(getReportOptions(), args, prog)
    #deal with the two scenarios for the segmentation specification
    #unlabelled path to a bed file, acceptable if there is only 1
    if (length(opt$segments) == 1 && !grepl(":", opt$segments)){
        opt$segments <- readRegions(opt$segments)
    } else {#labelled paths to bed files
        lp <- label_sc_path(opt$segments, unique.labels=TRUE)
        opt$segments <- 
            GRangesList(lapply(setNames(lp$path, lp$label), readRegions))
    }
    #call report
    htmlpath <- do.call(report, opt)
    cat("results written to the file: '", htmlpath, "'\n", sep="")
}

#' Produce an html report of a segmentation
#'
#' @param segments GRanges or GRangesList object containing the segmentation. 
#'  If GRanges, the \code{names} slot must be a number from 1 to the maximum 
#'  number of states. Same is true for the unlisted data in a GRangesList.
#'  If a GRangesList, the paths to of each bed file will depend on the names
#'  given by \code{names(segments)}.
#' @param model A list with the complete set of the parameters that describe
#' the HMM, such as that returned by the \code{segment} function.
#' @param outdir Output directory where the report will be created.
#' @param prefix Prefix prepended to all output files.
#' @param colors A character vector assigning one color per state.
#' each color must be a valid R color.
#' @param labels A character vector assigning one name per state
#' @param annots Named list where each item is a GRanges object to be
#' compared with the segmentation.
#' @param rdata R objects that will be saved as Rdata archives and will be
#' part of the report.
#' @param autoColors a list with two fields that controls how
#' the colors are automatically assigned to each state. The value \code{NULL}
#' disables automatic coloring.
#' @details A web page will be created with plots linked to
#'     tables in text format.
#' @return The path to the newly created webpage
#' @export
report <- function(segments, model, outdir=".", 
                        prefix="", colors=NULL, labels=NULL, annots=list(), 
                        rdata=NULL, autoColors=NULL){
    
    #CHECKING ARGUMENTS
    #you should make sure that when report is called from 'segmentCLI'
    #no error causes the program to abort (all assertions done here
    #should always pass)
    validateModel(model, strict=FALSE)
    nstates <- model$nstates
    #'segmlist' will store the GRangesList object, 'segments' the unlisted data
    if (inherits(segments, "GRanges")){
        segmlist <- GRangesList(segments)
    } else if (inherits(segments, "GRangesList")){
        segmlist <- segments
        segments <- unlist(segments, use.names=FALSE)
        #this is not a complete check, but it should be enough for most purposes
        ws <- sum(as(width(segmlist), "CompressedNumericList"))
        if (any(ws != ws[1])) stop("GRanges not on the same genomic regions")
    } else if (is.null(segments)) { stop("'segments' cannot be NULL")
    } else stop("'segments' can be a GRangesList or a GRanges object")
    
    snames <- names(segments)
    if (is.null(snames)) stop("segments must be annotated with the state number")
    inames <- as.integer(snames)
    if (any(is.na(inames) | inames <= 0 | inames > nstates)) stop("invalid segment names")
    if (!all(sapply(annots, inherits, what="GRanges"))) stop("'annots' must be a list of GenomicRanges objects")
    if (length(annots) > 0 && is.null(names(annots))) stop("'annots' elements must be named")
    validateOutdir(outdir)
    #SANITIZING INPUT
    #no errors allowed from here on, at most warnings
    prefix <- sanitizeFilename(prefix)
    if (!is.null(colors) && !is.null(autoColors)) {
        colors <- automaticColoring(getMeanMatrix(model$emisP), autoColors)
    }
    colors <- pickColors(nstates, colors)
    labels <- pickLabels(nstates, labels)
    
    #DOING REPORT
    doc <- c(
        htmlHeader("EpiCSeg report"), 
        "<center>",
        reportRdata(rdata, outdir, prefix),
        reportModel(model, labels, colors, outdir, prefix),
        reportSegmList(segmlist, labels, colors, outdir, prefix),
        reportAnnots(annots, segments, labels, colors, outdir, prefix),
        "</center>",
        htmlFooter())
    
    
    htmlpath <- makePath(outdir, prefix, "report.html")
    writeLines(doc, htmlpath)
    
    htmlpath
}

makePath <- function(outdir, prefix, last){
    file.path(outdir, paste0(prefix, last))
}

reportRdata <- function(rdata, outdir, prefix){
    if (is.null(rdata)) return(NULL)
    path <- makePath(outdir, prefix, "rdata.Rdata")
    save(rdata, file=path)
    htmlLink("additional segmentation data as an R archive", path)
}

reportNBs <- function(nbs, modelPath, outdir, prefix){
    path <- makePath(outdir, prefix, "nbs.png")
    png(path)
    plotNBs(nbs, main="total count across marks", xlab="mean +- standard deviation", ylab=NA)
    dev.off()
    #return html
    htmlImgLink(path, modelPath)
}

reportMeans <- function(means, modelPath, statecolors, outdir, prefix){
    path <- makePath(outdir, prefix, "means.png")
    #plot mean matrix
    myheat(t(means), xlab="mark", ylab="state", zlab="mean count", 
        main="mean counts", dev=path, col=heatpal, rowColors=statecolors)
    #return html
    htmlImgLink(path, modelPath)
}

reportLMeans <- function(lmeans, statecolors, outdir, prefix){
    paths <- makePath(outdir, prefix, c("lmeans.txt", "lmeans.png"))
    #write table
    write.table(lmeans, col.names=T, row.names=T, file=paths[1], quote=F, sep="\t")
    #make plot
    myheat(t(lmeans), xlab="mark", ylab="state", zlab="log(mean count + 1)", 
        main="log of mean counts", dev=paths[2], col=heatpal,rowColors=statecolors)
    #return html
    htmlImgLink(paths[2], paths[1])
}

reportTrans <- function(transP, modelPath, statecolors, outdir, prefix){
    path <- makePath(outdir, prefix, "transP.png")
    #plot transition probabilities
    myheat(transP, xlab="state to", ylab="state from", zlab="probability",
        main="transition\nmatrix", dev=path, col=heatpal, rowColors=statecolors, 
        colColors=statecolors)
    #return html
    htmlImgLink(path, modelPath)
}

reportModel <- function(model, labels, statecolors, outdir, prefix){
    #setting the right labels
    names(model$emisP) <- labels
    rownames(model$transP) <- labels
    colnames(model$transP) <- labels
    #convert parameters from list to matrix format
    nbs <- getNBs(model$emisP)
    means <- getMeanMatrix(model$emisP)
    rownames(means) <- model$marks
    #find a good permutation of the marks for plotting
    #determine a distance matrix between marks
    lmeans <- log(means + 1)
    dmat <- as.matrix(dist(lmeans, method="euclidean"))
    #get a small weight path through all the marks
    marksOrder <- smallWeightHamiltonianPath(dmat)
    #write the model
    modelPath <- makePath(outdir, prefix, "model.txt")
    writeModel(model, modelPath)
    
    nbs_html <- reportNBs(nbs, modelPath, outdir, prefix)
    means_html <- reportMeans(means[marksOrder,,drop=F], modelPath, statecolors, outdir, prefix)
    lmeans_html <- reportLMeans(lmeans[marksOrder,,drop=F], statecolors, outdir, prefix)
    trans_html <- reportTrans(model$transP, modelPath, statecolors, outdir, prefix)
    
    htmlSection("1. Model parameters",
        paste(sep="\n",
        nbs_html,"<br>",
        htmlMatToTable(valign="top", matrix(nrow=1, c(means_html, lmeans_html, trans_html)))
        )
    )
}

reportSegmList <- function(segmlist, labels, colors, outdir, prefix){
    nstates <- length(colors)
    #write colors
    colorsPath <- makePath(outdir, prefix, "colors.txt")
    writeColors(colors, colorsPath)
    #decide paths
    if (length(segmlist)==1 || is.null(names(segmlist))){
        segmPaths <- makePath(outdir, prefix, "segmentation.bed")
    } else {
        if (is.null(names(segmlist))) names(segmlist) <- 1:length(segmlist)
        segmPaths <- makePath(outdir, prefix, paste0("segmentation_", names(segmlist), ".bed"))
    }
    if (any(labels != as.character(1:nstates))) {
        segmPaths <- gsub(".bed$", "_labelled.bed", segmPaths) }
    
    for (i in seq_along(segmlist)){
        segmentsToBed(segmlist[[i]], labels, col2bedCol(colors), segmPaths[i])
    }
    
    #write table
    segments <- unlist(segmlist, use.names=FALSE)
    tabPaths <- makePath(outdir, prefix, c("table.txt", "table.png"))
    inames <- as.integer(names(segments))
    tab <- kfoots:::sumAt(width(segments), inames, size=nstates, zeroIdx=FALSE)
    tab <- tab/sum(tab)
    write.table(matrix(tab, ncol=1), col.names=F, row.names=F, file=tabPaths[1], quote=F, sep="\t")
    png(tabPaths[2])
    barplot(rev(tab), col=rev(colors), main="state frequencies", xlab="fraction of all bins",
                horiz=TRUE, border=NA, names.arg=rev(labels), las=2)
    dev.off()
    
    if (length(segmlist)==1) {
        bedLinks <- htmlLink("segmentation as a .bed file", segmPaths)
    } else {
        bedLinks <- paste0(collapse="<br>\n",
        htmlLink(
        paste0("segmentation as a .bed file (", names(segmlist), ")"), 
        segmPaths))
    }
    htmlSection("2. Segmentation",
        paste(sep="\n",
            htmlLink("state colors", colorsPath),"<br>",
            bedLinks,"<br>",
            htmlImgLink(tabPaths[2], tabPaths[1])
        )
    )
}

reportAnnot <- function(annot, name, segments, labels, colors, outdir, prefix){
    nstates <- length(colors)
    paths <- makePath(outdir, prefix, paste0("annot_", name, c(".txt", ".png")))
    cat("processing annotation: '", name, "'\n", sep="")
    prof <- avgStateProfile(annot, segments, nstates, before=5000, after=5000)
    rownames(prof) <- labels
    plotProfileAndLegend2Dev(prof, colors, dev=paths[2], main=paste0("states vs ", name))
    write.table(prof, file=paths[1], col.names=T, row.names=F, quote=F, sep="\t")
    
    htmlImgLink(paths[2], paths[1])
}

reportAnnots <- function(annots, segments, labels, colors, outdir, prefix){
    #process annotations independently
    annotDocs <- sapply(seq_along(annots), function(i){
        reportAnnot(annots[[i]], names(annots)[i], segments, labels, colors, outdir, prefix)
    })
    if (length(annotDocs)==0) return(NULL)
    
    htmlSection("3. Genomic annotations", htmlMatToTable(matrix(nrow=1, annotDocs)))
}


pickLabels <- function(nstates, labels){
    if (!is.null(labels) ){
        if (length(labels) != nstates) {
            warning("wrong number of state labels provided")
            labels <- NULL
        } else labels <- sanitizeFilename(labels)
    }
    if (is.null(labels)) labels <- as.character(1:nstates)
    labels
}

pickColors <- function(nstates, colors){
    #1. try to use user-defined colors
    if (!is.null(colors) && length(colors) != nstates){
            warning("wrong number of state colors provided")
            colors <- NULL
    }
    #check if the given colors are valid R colors
    tryCatch(col2bedCol(colors),
        error = function(e){
            warning(paste0("invalid colors:\n", e$message))
            colors <- NULL
    })
    #3. fall back to some arbitrary color palette
    if (is.null(colors)){
            set.seed(13)
            colors <- sample(qual.pal(nstates))
            set.seed(NULL)
    }
    colors
}


#average rank of counts distributed according to the given parameters (mu and r) 
#if considered together with all the empirical counts
expRank <- function(counts, mus, rs=Inf, normalize=T){
    if (length(rs)==1) rs <- rep(rs, length(mus))
    if (length(mus) != length(rs)) stop("invalid number of rs provided (one or as many as the mus)")
    #rank when the simulated count has a value equal or less to the maximal empirical count
    #(separately for each count in [0, length(countsTab)-1])
    countsTab <- kfoots:::tabFast(counts)
    ccountsTab <- cumsum(countsTab)
    totc <- length(counts)
    maxc <- length(countsTab)-1
    condRank <- countsTab/2 + 1 + c(0,ccountsTab[1:maxc])
    sapply(seq_along(mus), function(i){
        mu <- mus[i]
        r <- rs[i]
        if (is.finite(mu*r)){
            ps <- dnbinom(0:maxc, mu=mu, size=r)
            p <- pnbinom(maxc, mu=mu, size=r, lower.tail=F)
        } else {
            ps <- dpois(0:maxc, lambda=mu)
            p <- ppois(maxc, lambda=mu, lower.tail=F)
        }
        rnk <- p*(totc+1) + sum(ps*condRank)
        if (normalize) rnk <- (rnk-1)/totc
        rnk
    })
}

#in the argument 'means', rows are marks, columns are states
#the rows must be named
#the argument 'rs' is the r parameter for each state
#the order of the states must match between 'means' and 'rs'
#in the result, rows are marks, columns are states
expRankMatrix <- function(counts, means, rs, normalize=T){
    if (any(rownames(means) != rownames(counts))) stop("mean matrix not matching")
    
    #loop on the marks
    erm <- t(sapply(1:nrow(means), function(i){
        expRank(counts[i,], means[i,], rs, normalize=normalize)
    }))
    
    dimnames(erm) <- dimnames(means)
    erm
}


#METHODS TO VALIDATE AND SANITIZE ARGUMENTS
readAnnotations <- function(annotFields){
    lp <- label_sc_path(annotFields, unique.labels=TRUE)
    lp$label <- sanitizeFilename(lp$label)
    lapply(setNames(lp$path, lp$label), readRegions)
}


#check that an input string can be used
#as a filename 
sanitizeFilename <- function(fname){
    gsub("[^[:alnum:]_-]+", "_", fname)
}

checkWritable <- function(path){
    if (file.access(path, 0)<0){
        #file doesn't exist, see if you can create and delete one
        tryCatch({
            writeLines("checking path", path)
            file.remove(path)
        }, error = function(e){
            stop(paste0("something is wrong with the provided 'outdir' parameter:\n",
            e$message))
        })
    } else {
        #file exists, check that you have write permissions
        if (file.access(path, 2)<0) {
            stop(paste0("no write permissions to the file: ", path))}
    }
}

#this is not really a 'reader',
#it is only checking that oudir is a valid
#output directory by trying to create and delete
#a dummy file in that directory
validateOutdir <- function(outdir){
    if (length(outdir) != 1) stop("expecting a single path")
    path <- file.path(outdir, ".check_path.txt")
    checkWritable(path)
    outdir
}

