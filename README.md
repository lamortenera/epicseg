epicseg
=======

[![Build Status](https://travis-ci.org/lamortenera/epicseg.svg?branch=master)](https://travis-ci.org/lamortenera/epicseg)


Chromatin segmentation in R
## Installation

`epicseg` depends on Bioconductor packages, CRAN packages, and another package from github.
For the installation, most of the work is done by the function `devtools::install_github`. Because lately this function cannot resolve Bioconductor dependencies anymore (see this issue: https://github.com/hadley/devtools/issues/700), we will need to install some Bioconductor packages manually. Below is an example on how to do that. 

When I last modified this README, the Bioconductor dependencies were `IRanges`, `GenomicRanges`, `bamsignals` and `edgeR`. This list might go out of sync as `epicseg` evolves, have a look at the `DESCRIPTION` file to have the complete list (or just wait for error messages during installation) and modify the lines below accordingly, if needed.

```R
source("http://bioconductor.org/biocLite.R")
biocLite("IRanges", "GenomicRanges", "bamsignals", "edgeR")
```

Install and load the `devtools` package to be able to directly install R packages hosted on github :
```R
install.packages("devtools")
library(devtools)
```

To install `epicseg` type:

```R
install_github("lamortenera/kfoots")
install_github("lamortenera/epicseg")
```

To use the command line interface to epicseg, you need to create an executable to be used with Rscript:

```R
library(epicseg)
epicseg:::getLauncher("epicseg.R")
```

This creates a launcher called "epicseg.R" in your current working directory.
To check out the epicseg command line interface you need to have RScript.

If you have it, just type at the terminal, in the same folder where the file "epicseg.R" was created:

```bash
Rscript epicseg.R
```

And follow the instructions
