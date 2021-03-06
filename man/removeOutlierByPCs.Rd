% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypeQC.R
\name{removeOutlierByPCs}
\alias{removeOutlierByPCs}
\title{Remove population outliers}
\usage{
removeOutlierByPCs(
  plink,
  gcta,
  inputPrefix,
  nThread = 20,
  cutoff = NULL,
  cutoffSign,
  inputPC4subjFile,
  outputPC4outlierFile,
  outputPCplotFile,
  outputPrefix
)
}
\arguments{
\item{plink}{an executable program in either the current working directory 
or somewhere in the command path.}

\item{gcta}{an executable program in either the current working directory 
or somewhere in the command path.}

\item{inputPrefix}{the prefix of the input PLINK binary files.}

\item{nThread}{the number of threads used for computation. 
The default is 20.}

\item{cutoff}{the cutoff that distinguishes the outliers from ordinary
population using PCA. If it is null, then there are no outliers or 
outliers are not required to be removed. The default is NULL.}

\item{cutoffSign}{the cutoff sign: 'greater' or 'smaller' that determines 
if the outliers should be greater or smaller than the cutoff value.}

\item{inputPC4subjFile}{the pure text file that stores all the subject IDs 
and their corresponding eigenvalues of the first two principle components.}

\item{outputPC4outlierFile}{the pure text file that stores the outlier IDs 
and their corresponding eigenvalues of the first two principle components.}

\item{outputPCplotFile}{the plot file for visualizing the first two 
principle components of all subjects without population outliers.}

\item{outputPrefix}{the prefix of the output PLINK binary files.}
}
\value{
1.) The output PLINK binary files after outlier removal. 
2.) The output pure text file (if any) for storing removed outlier IDs 
and their corresponding PCs. 3.) The plot file (if any) for visualizing 
the first two principle components after outlier removal.
}
\description{
Remove population outliers by using principle component analysis.
}
\details{
This function is used for removing population outliers. 
If the outliers are necessary to be removed, then one uses the eigenvalues 
from the first principle component as a criterion to find out the outliers 
by assigning an appropriate cutoff.
}
\examples{
 
## In the current working directory
bedFile <- system.file("extdata", "QCdata.bed", package="Gimpute")
bimFile <- system.file("extdata", "QCdata.bim", package="Gimpute") 
famFile <- system.file("extdata", "QCdata.fam", package="Gimpute")
system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
inputPrefix <- "QCdata"  
cutoff <- NULL ## no outlier to be removed
cutoffSign <- "greater" ## not used if cutoff == NULL
inputPC4subjFile <- "2_13_eigenvalAfterQC.txt"
outputPC4outlierFile <- "2_13_eigenval4outliers.txt"
outputPCplotFile <- "2_13_removedOutliers.png"
outputPrefix <- "2_13_removedOutliers" 
## Not run: Requires an executable program PLINK and GCTA, e.g.
## plink <- "/home/tools/plink"
## gcta <- "/home/tools/gcta64"
## removeOutlierByPCs(plink, gcta, inputPrefix, nThread=20, 
##                    cutoff, cutoffSign, inputPC4subjFile, 
##                    outputPC4outlierFile, outputPCplotFile, outputPrefix)
}
\author{
Junfang Chen
}
