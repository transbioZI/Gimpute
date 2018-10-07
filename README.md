# Gimpute
Gimpute: An efficient genetic data processing and imputation pipeline

# Getting started  

## Installation 
Development version from Github:

1.) Install Gimpute from the command line
```{r eval=FALSE}
git clone https://github.com/transbioZI/Gimpute
R CMD build Gimpute
R CMD INSTALL Gimpute_*.tar.gz
```
2.) Install Gimpute in R (Recommended)
```{r eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("transbioZI/Gimpute", build_vignettes=TRUE)
```

Gimpute runs on any 64-bit x86 Linux distribution. Additional dependencies are described in the tutorial.

## Tutorial
Please check [Gimpute tutorial](https://github.com/transbioZI/Gimpute/blob/master/vignettes/GimputeTutorial.Rmd).

The best view of the tutorial is in HTML format with a table of contents by executing the following R codes after you have downloaded the tutorial (GimputeTutorial.Rmd in vignettes directory).

```{r eval=FALSE}
install.packages("rmarkdown")
library("rmarkdown")
render("GimputeTutorial.Rmd")
```
## Reference

Chen, J., Lippold, D., Frank, J., Rayner, W., Meyer-Lindenberg, A., Schwarz, E., & Stegle, O. (2018). Gimpute: An efficient genetic data imputation pipeline. Bioinformatics.
