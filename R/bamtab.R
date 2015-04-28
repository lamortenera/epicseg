bamtabDefaults <- list(shift=75, mapq=0, pairedend=FALSE)

validateBamtab <- function(bamtab){
    if (!is.data.frame(bamtab)) stop("'bamtab' must be a 'data.frame'")
    reqfields <- c("mark", "path")
    optfields <- names(bamtabDefaults)
    if (!all(reqfields %in% names(bamtab))) stop("missing required fields")
    if (!all(names(bamtab) %in% c(reqfields, optfields))) stop("invalid fields")
    #add defaults
    for (n in names(bamtabDefaults)){
        if (!n %in% names(bamtab)) {
            bamtab[[n]] <- rep(bamtabDefaults[[n]], nrow(bamtab))
        }
    }
    #check path
    if (!is.character(bamtab$path)) stop("invalid path specification")
    if (any(!file.exists(bamtab$path))) stop("BAM file does not exist")
    #check mapq
    if(!is.numeric(bamtab$mapq) || any(bamtab$mapq < 0 | bamtab$mapq > 255)){
        stop("invalid 'mapq'")}
    #check shift
    if (!is.numeric(bamtab$shift)) stop("invalid 'shift'")
    #check pairedend
    if (!is.logical(bamtab$pairedend)) stop("invalid 'pairedend'")

    bamtab
}

makeBamtab <- function(mark_sc_path, shift=NULL, mapq=NULL, pairedend=NULL){
    npaths <- length(mark_sc_path)
    mp <- strsplit(mark_sc_path, split=":")
    if (any(sapply(mp, length) != 2)) stop("invalid mark specification")
    mp <- simplify2array(mp)
    bamtab <- data.frame(mark=mp[1,], path=mp[2,], stringsAsFactors=F)
    for (nm in names(bamtabDefaults)){
        if (!is.null(get(nm))) bamtab[[nm]] <- fixLength(get(nm), npaths)
    }
    bamtab
}

fixLength <- function(v, len){
    if (length(v)==1) return(rep(v, len))
    if (length(v)==len) return(v)
    stop(paste0("cannot interpret the given vector as a vector of length ", len,
    ":\n provide either ", len, " elements or just 1"))
}
