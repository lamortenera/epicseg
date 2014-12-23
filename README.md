epicseg
=======

Chromatin segmentation in R
## Installation

Install and load the `devtools` package to be able to directly install R packages hosted on github :

```R
install.packages("devtools")
load("devtools")
```

To install `epicseg` type:

```R
devtools::install_github("lamortenera/bamsignals")
devtools::install_github("lamortenera/kfoots")
devtools::install_github("lamortenera/epicseg")
```

To use the command line interface to epicseg, you need to create an exectuable to be used with Rscript:

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
