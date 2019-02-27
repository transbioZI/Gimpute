## Gimpute: An efficient genetic data processing and imputation pipeline

### Function 
In order to ensure the reliability and reproducibility of genome-wide association study (GWAS) data, and enable meta-analysis across cohorts from different genotyping arrays. we set up an efficient, automatic and comprehensive genotype data processing and imputation pipeline termed __Gimpute__. It consists of pre-processing (genetic variant information updating/matching/liftOver, quality control of genetic variants and study samples, the alignment of variants to the imputation references), pre-phasing and imputation, as well as post-imputation quality control. 


### Installation 

Install Gimpute in R:
```{r eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("transbioZI/Gimpute", build_vignettes=TRUE)
```

Gimpute runs on any 64-bit x86 Linux distribution. Additional dependencies are described in the [tutorial](https://github.com/transbioZI/Gimpute/blob/master/vignettes/GimputeTutorial.Rmd).

### Tutorial
The detailed instruction is explained in the [Gimpute tutorial](https://github.com/transbioZI/Gimpute/blob/master/vignettes/GimputeTutorial.Rmd), along with a complete [running example](https://github.com/transbioZI/Gimpute/blob/master/tests/runTests.R). 

The best view of the tutorial is in HTML format with a table of contents by executing the following R codes after you have downloaded the tutorial (GimputeTutorial.Rmd in vignettes directory).

```{r eval=FALSE}
install.packages("rmarkdown")
library("rmarkdown")
render("GimputeTutorial.Rmd")
```
### Citation

Chen, J., et al. (2018). Gimpute: an efficient genetic data imputation pipeline. Bioinformatics.