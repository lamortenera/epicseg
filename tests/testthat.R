library(testthat)
library(epicseg)

filt <- NULL
#we don't run the CLI tests on Windows
if (.Platform$OS.type != "unix") filt <- "^[A-Za-z]"
test_check("epicseg", filter=filt, reporter="summary")
