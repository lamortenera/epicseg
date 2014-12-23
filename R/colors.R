defaultAutoColors <- list(markToCol=list(red="H3K4me3", green4="H3K36me3", blue=c("H3K9me3", "H3K27me3"), gold="H3K4me1"), background="gray")

validateAutoColors <- function(autoColors, marks){
	if (length(marks) < length(autoColors$markToCol) + !is.null(autoColors$background)){
		stop("too few marks for the given color scheme")
	}
	for (m2c in autoColors$markToCol){
		if (!any(m2c %in% marks)) stop("marks in the color scheme not in the data")
	}
}

#converts R colors to colors for the bed format
col2bedCol <- function(colors){
	rgbs <- col2rgb(colors)
	apply(rgbs, 2, paste, collapse=",")
}


#extends the palettes from RColorBrewer to an infinite number of colors
#this is thought for sequential palettes, results with other palettes are weird
seq.pal <- function(n, name="RdYlGn"){
	if (!name %in% rownames(brewer.pal.info)) stop("palette not available in RColorBrewer")
	maxcol <- brewer.pal.info[name, "maxcolors"]
	origpal <- brewer.pal(maxcol, name)
	cramp <- colorRamp(origpal)
	rgb(cramp(seq(0, 1, length.out=n)), max=255)
}

#extends the palettes from RColorBrewer to an infinite number of colors
#this is thought for qualitative palettes
qual.pal <- function(n, name="Paired"){
	if (!name %in% rownames(brewer.pal.info)) stop("palette not available in RColorBrewer")
	maxcol <- brewer.pal.info[name, "maxcolors"]
	#RColorBrewer issues a warning in this case (for this palette)
	if (n < 3) return(brewer.pal(3, name)[1:n])
	if (n < maxcol) return(brewer.pal(n, name))
	origpal <- brewer.pal(maxcol, name)
	cramp <- colorRamp(origpal)
	samples <- seq(0, 1, length.out=n)
	#snap floating point nums to multiples 1/(maxcol-1) to get
	#as many original colors as possible
	origsamples <- seq(0, 1, length.out=maxcol)
	for (os in origsamples){
		dists <- abs(os - samples)
		samples[which.min(dists)] <- os
	}
	rgb(cramp(samples), max=255)
}

heatpal <- rev(seq.pal(100))


toRGBVector <- function(v){
	if (length(v)!=3) return(rep(NA,3))
	iv <- as.integer(v)
	if (any(iv < 0 || iv > 255)) return(rep(NA,3))
	iv
}

validateColors <- function(colors){
	invisible(col2rgb(colors))
}

writeColors <- function(colors, path){
	validateColors(colors)
	writeLines(col2bedCol(colors), path)
}

readColors <- function(path){
	txt <- readLines(path)
	#first, check if the colors are in bed format (R,G,B)
	rgbmat <- t(sapply(strsplit(txt, ","), toRGBVector))
	if (!any(is.na(rgbmat))){
		#the colors are in bed format
		rgb(rgbmat, maxColorValue=255)
	} else {
		#check if we can use them as R colors
		validateColors(txt)
		txt
	}
}
