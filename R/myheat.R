#like the image function, but with additional xlim and ylim arguments that specify
#(using user coordinates) in which subportion of the plot to draw the image

#plot.new() must have already been called in order to use this function
myimage <- function(mat, xlim, ylim, col, border=NA, lty=par("lty"), 
    lwd=par("lwd"), zlim=range(mat[is.finite(mat)]), ...){
    
    stopifnot(is.matrix(mat))
    stopifnot(ncol(mat)*nrow(mat) > 0)
    cs <- col(mat)
    rs <- row(mat)#nrow(mat) - row(mat) + 1
    xlen <- diff(xlim)
    ylen <- diff(ylim)
    xdelta <- xlen/ncol(mat)
    ydelta <- ylen/nrow(mat)
    
    xall <- xlim[1] + xdelta*(0:ncol(mat))
    yall <- ylim[1] + ydelta*(0:nrow(mat))
    xlow <- xall[cs]
    xhigh <- xall[cs+1]
    ylow <- yall[rs]
    yhigh <- yall[rs+1]
    #in the typical representation of a matrix the first row is at the top
    ylow <- rev(ylow); yhigh <- rev(yhigh)
    
    #a matrix where for each element of the matrix we have 5 coordinates (bottom-left, top-left, top-right, bottom-right, bottom-left)
    xs <- rbind(xlow, xlow, xhigh, xhigh, xlow, NA)
    ys <- rbind(ylow, yhigh, yhigh, ylow, ylow, NA)
    
    #make sure that when mat has values in (1+o):(k+o) (integers) and you provide k colors
    #nothing gets interpolated and value i+o gets color i.
    #the minimum should get the first color, and the maximum the last color
    nc <- length(col)
    if (!missing(zlim) && (any(!is.finite(zlim)) || diff(zlim) <= 0)) 
        stop("invalid z limits")
    zi <- round((nc-1)*(mat - zlim[1L])/diff(zlim)) + 1
    zi[zi <= 0 | zi > nc] <- NA

    colors <- col[zi]
    
    polygon(xs, ys, col=colors, border=NA, ...)
    
    if (!(is.null(border) || is.na(border))){
        #draw vertical lines
        #a matrix where for each line we take start, end and NA
        vxs <- rbind(xall, xall, NA)
        vys <- rbind(rep(ylim[1], length(xall)), rep(ylim[2], length(xall)), NA)
        lines(vxs, vys, lty=lty, col=border)
        
        hxs <- rbind(rep(xlim[1], length(yall)), rep(xlim[2], length(yall)), NA)
        hys <- rbind(yall, yall, NA)
        lines(hxs, hys, lty=lty, lwd=lwd,   col=border)
    }
}


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

#get the margin size in inches
getMAI <- function(mat, targetStrH, dev, xlab, ylab, zlab, main){
    #to get some R parameters we need to open and close the device...
    openDevice(dev)
    
    #the difference between csi and strH is that csi contains also the interline
    #as far as I understood
    csi <- par("csi") #height of (default-sized) characters in inches.
    strH <- strheight("M", units="inches")#line height in inches with cex=1
    cex <- targetStrH/strH
    par(cex=cex)
    xlablen <- max(strwidth(colnames(mat),units="inches"))
    ylablen <- max(strwidth(rownames(mat),units="inches"))
    mar2mai <- mean(par("mai")/par("mar"))
    dev.off()
    #margin size
    #constant additional spacing (not due to label lengths)
    #if there are axes labels, add two more lines of margin
    mar <- c(3.1,2.1,2.1,.5) + 1*c(!is.null(xlab), !is.null(ylab), !is.null(main), (!is.null(zlab))*3) 
    #margins using height of a line of text as a measure
    #margins using inches as a measure
    mai <- mar*mar2mai + c(xlablen, ylablen, 0, 0)
    mai
}

squareCoords <- function(mat, hasRC, hasCC, hasCK){
    #user space limits measured in heatmap squares
    sxlim <- c(0, ncol(mat))
    sylim <- c(0, nrow(mat))
    #icolor key we need 1.5 more squares
    if (hasCK){
        sxlim[2] <- sxlim[2]+1.5
    }
    #row colors, 1.5 more squares 
    if (hasRC){
        sxlim[1] <- sxlim[1]-1.5
    }
    #column colors, 1.5 more squares
    if (hasCC){
        sylim[1] <- sylim[1]-1.5
    }
    list(sxlim=sxlim, sylim=sylim)
}

domyheat <- function(mat, xlab=NULL, ylab=NULL, main=NULL, col=heatpal,
        zlim=range(mat, na.rm=T), zlab="value", rowColors=NULL, colColors=NULL, ...){
    
    nx <- ncol(mat)
    ny <- nrow(mat)
    scoords <- squareCoords(mat, !is.null(rowColors), !is.null(colColors), !is.null(zlab))
    sxlim <- scoords$sxlim
    sylim <- scoords$sylim
    
    plot(NA, NA, xlim=sxlim, ylim=sylim, xaxs="i", yaxs="i", xlab=NA, ylab=NA, main=NA, xaxt="n", yaxt="n", bty="n")
    #plot the heatmap
    myimage(mat, zlim=zlim, xlim=c(0,nx), ylim=c(0,ny), col=col, ...)
    
    #plot the color key
    if (!is.null(zlab)){
        zseq <- seq(zlim[1], zlim[2], length.out=length(col))
        ck.mat <- matrix(ncol=1,zseq)
        myimage(ck.mat, col=col, zlim=zlim, xlim=c(nx+.5, nx+1.5), ylim=c(0,ny))
        zTicks <- mypretty(zseq)
        zAt <- ny*(zTicks-zlim[1])/diff(zlim)
        par(las=0)
        axis(side=4, labels=zTicks, at=zAt, lwd=par("cex"), line=0)
        mtext(zlab, side=4, line=2, cex=par("cex"), at = mean(c(0,ny)))
    }
    #plot the row colors
    if (!is.null(rowColors)){
        rc.mat <- matrix(ncol=1, 1:ny)
        myimage(rc.mat, col=rowColors, xlim=c(-1.5,-.5), ylim=c(0,ny))
    }
    #plot the column colors
    if (!is.null(colColors)){
        cc.mat <- matrix(nrow=1, 1:nx)
        myimage(cc.mat, col=colColors, xlim=c(0,nx), ylim=c(-1.5, -.5))
    }
    
    par(las=2)
    #xaxis
    if (!is.null(colnames(mat))){
        axis(side=1, at=(1:nx-.5), labels=colnames(mat), tick=F)
    }
    #yaxis
    if (!is.null(rownames(mat))){
        yticks <- rev(rownames(mat))
        axis(side=2, at=(1:ny-.5), labels=yticks, tick=F)
    }
    #xaxis label
    par(las=0)
    #main
    if (!is.null(main)){
        mtext(main, side=3, line=1.5, cex=1.2*par("cex"), font=2, at=mean(c(0, nx)))
    }
    if (!is.null(xlab)) {
        mtext(xlab, side=1, cex=par("cex"), 
        line=(par("mar")[1]-2.1), at = mean(c(0,nx))) }
    #yaxislabel
    if (!is.null(ylab)) {
        mtext(ylab, side=2, cex=par("cex"), line=(par("mar")[2]-1.1), 
        at = mean(c(0,ny)))}
    

}

myheat <- function(mat, dev=NULL, L=0.3, L2strH=0.45, xlab=NULL, ylab=NULL, main=NULL, 
        col=heatpal, zlim=range(mat, na.rm=T), zlab="value", 
        rowColors=NULL, colColors=NULL, ...){
    
    nx <- ncol(mat)
    ny <- nrow(mat)
    if (!is.null(rowColors) && length(rowColors)!=ny) stop("invalid rowColors")
    if (!is.null(colColors) && length(colColors)!=nx) stop("invalid colColors")
    ck <- !is.null(zlab)
    
    #line height in inches
    targetStrH <- L*L2strH
    #figure margins in inches so that the labels fit
    mai <- getMAI(mat, targetStrH, dev, xlab, ylab, zlab, main)
    
    #plot limits measured in heatmap squares
    scoords <- squareCoords(mat, !is.null(rowColors), !is.null(colColors), ck)
    sxlim <- scoords$sxlim
    sylim <- scoords$sylim
    
    #plot size in inches
    pxi <- L*diff(sxlim)
    pyi <- L*diff(sylim)
    
    #figure size in inches
    xi <- pxi + sum(mai[c(2,4)])
    yi <- pyi + sum(mai[c(1,3)])
    
    #opening device
    openDevice(dev, Win=xi, Hin=yi)
    cex <- targetStrH/strheight("M", units="inches")
    
    par(mai=mai, cex=cex, las=2)
    #start a plot with user coordinates given by the squares from the bottom-left corner
    
    domyheat(mat, xlab=xlab, ylab=ylab, main=main, col=col, zlim=zlim, zlab=zlab, rowColors=rowColors, colColors=colColors, ...)
    
    #closing device, if it's not X11
    if (!is.null(dev)) dev.off()
}


mypretty <- function(x){
    cs <- pretty(x, 4, 4)
    lim <- range(x)
    csin <- cs[cs>= lim[1] & cs <= lim[2]]
    if (length(csin)==0) return(cs)
    csin
}

