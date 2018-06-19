# Gimpute
Gimpute: An efficient genetic data processing and imputation pipeline


# Getting started  
## Prerequisites
Gimpute runs on any 64-bit x86 Linux distribution and it requires the following tools:

* [R](https://www.r-project.org/) 
* [PLINK](https://www.cog-genomics.org/plink2) 
* [GCTA64](http://cnsgenomics.com/software/gcta/#Download) 
* [SHAPEIT](http://www.shapeit.fr/) 
* [IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html) 
* [GTOOL](http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html) 

## Installation 
Development version from Github:
```{r eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("transbioZI/Gimpute")
```
This function`install_github()` requires that you build from source, namely, `make` and compilers must be installed on the system.

For more information, please look at [Gimpute tutorial](https://github.com/transbioZI/Gimpute/blob/master/vignettes/GimputeTutorial.Rmd).
