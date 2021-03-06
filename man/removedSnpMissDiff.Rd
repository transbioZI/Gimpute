% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypeQC.R
\name{removedSnpMissDiff}
\alias{removedSnpMissDiff}
\title{Remove SNPs with difference in SNP missingness between cases and controls.}
\usage{
removedSnpMissDiff(
  plink,
  inputPrefix,
  snpMissDifCutOff,
  outputPrefix,
  groupLabel
)
}
\arguments{
\item{plink}{an executable program in either the current working directory 
or somewhere in the command path.}

\item{inputPrefix}{the prefix of the input PLINK binary files.}

\item{snpMissDifCutOff}{the cutoff of the difference in missingness between 
cases and controls.}

\item{outputPrefix}{the prefix of the output PLINK binary files.}

\item{groupLabel}{a string value indicating the outcome label: "control",
or, "case" or "caseControl" for both existing groups. For more details, see 
\code{\link{getGroupLabel}}.}
}
\value{
The output PLINK binary files.
}
\description{
Remove SNPs with difference in SNP missingness between cases and controls. 
To test for differential call rates between cases and controls for each SNP
}
\details{
Only if both case-control groups exist in the input genotype data, 
differential SNPs are removed.
}
\examples{
 
## In the current working directory
bedFile <- system.file("extdata", "genoUpdatedData.bed", package="Gimpute")
bimFile <- system.file("extdata", "genoUpdatedData.bim", package="Gimpute") 
famFile <- system.file("extdata", "genoUpdatedData.fam", package="Gimpute")
system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
inputPrefix <- "genoUpdatedData" 
snpMissDifCutOff <- 0.02
outputPrefix <- "2_09_removedSnpMissDiff" 
groupLabel <- "control"
## Not run: Requires an executable program PLINK, e.g.
## plink <- "/home/tools/plink"
## removedSnpMissDiff(plink, inputPrefix, snpMissDifCutOff, 
##                    outputPrefix, groupLabel)
}
\seealso{
\code{\link{getGroupLabel}}.
}
\author{
Junfang Chen
}
