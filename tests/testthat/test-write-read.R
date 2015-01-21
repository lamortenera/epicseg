context("A: write and read")

test_that("writing and reading counts",{
	counts <- matrix(nrow=10, rpois(1000, lambda=100))
	rownames(counts) <- 1:nrow(counts)
	
	#write and read in text format
	path <- tempfile()
	epicseg:::writeCounts(counts, path)
	wrote_and_read_counts <- epicseg:::readCounts(path)
	expect_identical(counts, wrote_and_read_counts)
	
	#write and read in Rdata archive
	path <- paste0(path, ".rda")
	epicseg:::writeCounts(counts, path)
	wrote_and_read_counts <- epicseg:::readCounts(path)
	expect_identical(counts, wrote_and_read_counts)
})
