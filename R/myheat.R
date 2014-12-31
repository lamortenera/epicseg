openDevice <- function(dev, Win=7, Hin=7){
	if (is.null(dev)){#this should be X11... dunno what happens in windows
		dev.new(width=Win, height=Hin)
	} else {
		if (grepl(".png$", dev)){
			png(dev, width=72*Win, height=72*Hin)
		} else if (grepl(".pdf$", dev)){
			pdf(dev, width=Win, height=Hin)
		} else stop("the path should end in .png or .pdf")
	}
}

makeBoundingBox <- function(bty, bty.lwd, bty.col, nx, ny){
	if (!bty=="n"){
		usr <- par("usr")
		if (bty=="o"){
			lines(rep(usr[1], 2), usr[3:4], lwd=bty.lwd, col=bty.col) #first vertical
			lines(rep(usr[2], 2), usr[3:4], lwd=bty.lwd, col=bty.col) #last vertical
			lines(usr[1:2], rep(usr[3],2), lwd=bty.lwd, col=bty.col) #first horizontal
			lines(usr[1:2], rep(usr[4],2), lwd=bty.lwd, col=bty.col) #last horizontal
		} else if (bty=="l"){
			lines(rep(usr[1], 2), usr[3:4], lwd=bty.lwd, col=bty.col) #first vertical
			lines(usr[1:2], rep(usr[3],2), lwd=bty.lwd, col=bty.col) #first horizontal
		} else if (bty=="#"){
			xs <- seq(usr[1], usr[2], length.out=(nx+1))
			ys <- seq(usr[3], usr[4], length.out=(ny+1))
			
			xsL <- as.vector(rbind(xs, xs, rep(NA, length(xs))))
			ysL <- as.vector(rbind(ys, ys, rep(NA, length(ys))))
			lines(xsL, rep(c(usr[3:4], NA), nx+1), lwd=bty.lwd, col=bty.col) #all verticals
			lines(rep(c(usr[1:2], NA), ny+1), ysL, lwd=bty.lwd, col=bty.col) #all horizontals
		} else stop("bty must be one of 'o', 'l', and '#'")
	}
}


myheat <- function(mat, dev=NULL, L=0.3, L2lheight=0.5, xlab=NULL, ylab=NULL, main=NULL, 
		col=heatpal,
		zlim=range(mat, na.rm=T), 
		zlab="value",
		ck.antialiasing=F,
		L2lwd=2/0.3, 
		bty="n",
		bty.col="black", symm=F, ...){
	
	nx <- ncol(mat)
	ny <- nrow(mat)
	if (symm){
		if (is.null(colnames(mat)) && !is.null(rownames(mat))) colnames(mat) <- rownames(mat)
		if (is.null(rownames(mat)) && !is.null(colnames(mat))) rownames(mat) <- colnames(mat)
	}
	
	if (is.null(colnames(mat))) {colnames(mat) <- 1:nx}
	if (is.null(rownames(mat))) {rownames(mat) <- 1:ny}
	
	
	#to get some R parameters we need to open and close the device...
	openDevice(dev)
	csi <- par("csi")
	strH <- strwidth("M", units="inches")
	xlablen <- max(strwidth(colnames(mat),units="inches")/csi)
	ylablen <- max(strwidth(rownames(mat),units="inches")/csi)
	dev.off()
	
	#plot size in inches
	pxi <- L*nx
	pyi <- L*ny
	#margin size in inches
	mar <- c(2.1,1.1,2,1) + 2*c(!is.null(xlab), !is.null(ylab), !is.null(main), 0) #constant additional spacing (not due to label lengths)
	mar <- mar + c(xlablen, ylablen, 0, 0)#if there are axes labels, add two more lines of margin
	mai <- L*L2lheight*mar*csi/strH #keeping the same ratio between line height and interline spacing as in standard R plotting
	
	#figure size in inches
	xi <- pxi + sum(mai[c(2,4)])
	yi <- pyi + sum(mai[c(1,3)])
	
	if (!is.null(zlab)){
		#we also plot the color key!
		ck.pxi <- L
		ck.mar <- c(mar[1], 0, mar[3], 3.1)
		ck.mai <- ck.mar*L*L2lheight*csi/strH
		ck.xi <- ck.pxi + sum(ck.mai[c(2,4)])
	} else {
		ck.xi <- 0
	}
	
	
	#opening device
	openDevice(dev, Win=(xi+ck.xi), Hin=yi)
	#backup the original graphical parameters
	opar <- par()
	#set up the new ones
	cex <- L*L2lheight/strwidth("M", units="inches")
	
	if (!is.null(zlab)) layout(matrix(c(1,2), nrow=1), widths=c(xi, ck.xi))
	
	par(mai=mai, cex=cex, las=2, bty="n")
	#plot the heatmap
	if (symm){#delete the cells in the upper diagonal
		n <- ncol(mat)
		for (i in 1:(n-1)){ mat[i, (i+1):n] <- NA }
	}
	image(t(mat[ny:1,]), zlim=zlim, xaxt="n", yaxt="n", main=main, col=col,...)
	
	#xaxis
	xticks <- colnames(mat)
	xpos <- (0:(nx-1))/max(1,(nx-1))
	#yaxis
	yticks <- rownames(mat)
	ypos <- rev((0:(ny-1))/max(1,(ny-1)))
	axis(side=1, at=xpos, labels=xticks, tick=F)
	axis(side=2, at=ypos, labels=yticks, tick=F)
	#xaxis label
	par(las=0)
	if (!is.null(xlab)) mtext(xlab, side=1, line=mar[1]-2.1, cex=par("cex"))
	#yaxislabel
	if (!is.null(ylab)) mtext(ylab, side=2, line=mar[2]-1.1, cex=par("cex"))
	#bounding box
	makeBoundingBox(bty, L*L2lwd, bty.col, nx, ny)
	
	if (!is.null(zlab)){
		#plot the color key!!!
		#plot.new()
		par(mai=ck.mai, cex=cex, las=0, bty="n")
		n <- length(col)
		ck.z <- seq(zlim[1], zlim[2], length.out=n)
		ck.mat <- t(as.matrix(ck.z))
		ck.x <- 0
		ck.y <- ck.z
		image(ck.x, ck.y, ck.mat, xaxt="n", yaxt="n", col=col, xlab=NA, ylab=NA)
		if (ck.antialiasing){
			usr <- par("usr")
			ys <- seq(usr[3], usr[4], length.out=(n+1)) #all the y coordinates between two rectangles
			ys <- ys[2:n] #discard the first and the last
			#interpolate colors
			cramp <- colorRamp(col)
			step <- 1/(n-1)
			colMiddle <- rgb(cramp(seq(from=step/2, to=1-(step/2), length.out=(n-1))), maxColorValue=255)
			#plot lines to glue together adjacent rectangles
			for (i in 1:(n-1)) abline(h=ys[i], lwd=1, col=colMiddle[i])
			#lines(xsL, ysL, col=cols, lwd=1, type="h")
		}
		axis(side=4, lwd=cex, lwd.ticks=cex)
		mtext(zlab, side=4, line=2, cex=par("cex"))
	}
	
	
	#restore previous parameters (some of them cannot be restored and we would get a warning)
	options(warn=-1)
	par(opar)
	options(warn=0)
	#closing device, if it's not X11
	if (!is.null(dev)) dev.off()
}


