#METHODS TO TYPESET HTML
htmlHeader <- function(title=NULL){
	header <- "<!DOCTYPE html> <html> "
	if (!is.null(title)) header <- c(header, paste0("<head><title>", title, "</title></head>"), paste0("<body><h1><center>", title, "</center></h1><br>"))
	else header <- c(header, "<body>")
	header
}

htmlFooter <- function(){
	"</body> </html>"
}

htmlSection <- function(name, contents){
	paste0(
	"<div>",
	"<h2>",name,"</h2><br>",
	contents,
	"</div>")
}

#convert a matrix of html elements to a html table
htmlMatToTable <- function(htmls, valign=NULL){
	if (length(htmls)<=0) return("")
	paste0(collapse="\n",
	c("<table>",
	apply(htmls, 1, htmlTableRow, valign=valign),
	"</table>")
	)
}

htmlTableRow <- function(contents, valign=NULL){
	if (!is.null(valign)){
		rowHtml <- paste0(collapse="", "<td valign=\"", valign, "\">", contents, "</td>")
	} else {
		rowHtml <- paste0(collapse="", "<td>", contents, "</td>")
	}
	paste0("<tr>", rowHtml, "</tr>")
}

#creates the html for an image where you can click
#and it sends you to the supporting data for that image
htmlImgLink <- function(imgpath, linkpath){
	htmlLink(paste0("<img src=\"", basename(imgpath), "\">"), linkpath)
}

#creates an html link to the destination specified by linkpath
#inside the <a> tag the content 'content' is displayed
htmlLink <- function(content, linkpath){
	paste0("<a href=\"", basename(linkpath), "\">",content, "</a>")
}

