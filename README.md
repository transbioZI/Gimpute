# Gimpute
Gimpute: An efficient genetic data imputation pipeline 

# Getting started  
## Prerequisites
Gimpute runs in any 64-bit x86 Linux distribution and it requires the following tools:

* [R](https://www.r-project.org/) 
* [PLINK](https://www.cog-genomics.org/plink2) 
* [GCTA64](http://cnsgenomics.com/software/gcta/download.html) 
* [SHAPEIT](http://www.shapeit.fr/) 
* [IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html) 
* [GTOOL](http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html) 

## Installation 
Development version from Github:
```{r eval=FALSE}
library("devtools")
install_github("Junfang/Gimpute")
```
This function`install_github()` requires that you build from source, namely, `make` and compilers must be installed on the system.
 
