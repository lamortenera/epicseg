#make a list of count matrices
nmats <- 4; nc <- 500; nr <- 5
clist <- lapply(1:nmats, function(r) matrix(rpois(nr*nc, lambda=500), nrow=nr))
for (i in seq_along(clist)) rownames(clist[[i]]) <- paste0("mark", 1:nr)
names(clist) <- paste0("dataset",1:nmats)


#make some matching regions
binsize <- 200
gr <- GRanges(seqnames="chr1", IRanges(start=1, width=binsize*nc))

#segment
segmentation <- segment(clist, gr, 5, verbose=F)
