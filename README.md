## EpiCSeg  
[![Build Status](https://travis-ci.org/lamortenera/epicseg.svg?branch=master)](https://travis-ci.org/lamortenera/epicseg)

Chromatin segmentation in R. For a full description and for citing us, see the following [article](http://www.genomebiology.com/2015/16/1/151) on Genome Biology:

> A Mammana, HR Chung (2015). Chromatin segmentation based on a probabilistic model for read counts explains a large portion of the epigenome. Genome Biology 2015, 16:151

For using EpiCSeg from the command line, you can find the full manual [HERE!](https://cdn.rawgit.com/lamortenera/epicseg/master/inst/manual.html)


### Installation

`epicseg` needs R 3.2 (or newer) and depends on Bioconductor packages, CRAN packages, and another package from github. 
For the installation, most of the work is done by the function `devtools::install_github`. Because lately this function cannot resolve Bioconductor dependencies anymore (see this issue: https://github.com/hadley/devtools/issues/700), we will need to install some Bioconductor packages manually.

The Bioconductor dependencies are `IRanges`, `GenomicRanges`, `bamsignals` and `edgeR`. At the interactive R terminal, type:

```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("IRanges", "GenomicRanges", "bamsignals", "edgeR"))
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

### Usage from the command line

To use EpiCSeg from the command line, you need to:

1. create a launcher to be used with Rscript. This is done
by typing `epicseg:::getLauncher("epicseg.R")` at the R interactive 
terminal, which will create the file `epicseg.R` in your working directory. 
You can move and rename this file the way you want. 
2. To use it, type `Rscript epicseg.R subprogram arguments`.
In UNIX you can also simply do `./epicseg.R subprogram arguments` provided that
you have execution permission on the file `.epicseg.R`. 
3. To see what the available subprograms are, simply type: 
`Rscript epicseg.R` 
4. To see which arguments each subprogram needs, you can type: 
`Rscript epicseg.R subprogram`


