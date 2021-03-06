% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypeQC.R
\name{genoQC}
\alias{genoQC}
\title{Quality control for genotype data}
\usage{
genoQC(
  plink,
  inputPrefix,
  snpMissCutOffpre = 0.05,
  sampleMissCutOff = 0.02,
  Fhet = 0.2,
  cutoffSubject,
  cutoffSNP,
  snpMissCutOffpost = 0.02,
  snpMissDifCutOff = 0.02,
  femaleChrXmissCutoff = 0.05,
  pval4autoCtl = 1e-06,
  pval4femaleXctl = 1e-06,
  outputPrefix,
  keepInterFile = TRUE
)
}
\arguments{
\item{plink}{an executable program in either the current working directory 
or somewhere in the command path.}

\item{inputPrefix}{the prefix of the input PLINK binary files.}

\item{snpMissCutOffpre}{the cutoff of the missingness for removing SNPs 
before subject removal. The default is 0.05.}

\item{sampleMissCutOff}{the cutoff of the missingness for removing 
subjects/instances. The default is 0.02.}

\item{Fhet}{the cutoff of the autosomal heterozygosity deviation. 
The default is 0.2.}

\item{cutoffSubject}{the cutoff determines that families (subjects) with 
more than the predefined cutoff of Mendel errors by considering all SNPs 
will be removed. The default is 0.05.}

\item{cutoffSNP}{the cutoff indicates that SNPs with more than the 
predefined cutoff of Mendel error rate will be excluded 
(i.e. based on the number of trios/duos). The default is 0.1.}

\item{snpMissCutOffpost}{the cutoff of the missingness for removing SNPs 
after subject removal. The default is 0.02.}

\item{snpMissDifCutOff}{the cutoff of the difference in missingness between 
cases and controls. The default is 0.02.}

\item{femaleChrXmissCutoff}{the cutoff of the missingness in female 
chromosome X SNPs. The default is 0.05.}

\item{pval4autoCtl}{the p-value cutoff for controlling HWE test in either 
control or case subjects. Only autosomal SNPs are considered. 
The default is 0.000001}

\item{pval4femaleXctl}{the p-value cutoff for controlling HWE test in 
female control subjects. Only chromosome X SNPs are considered. 
The default is 0.000001}

\item{outputPrefix}{the prefix of the output PLINK binary files after QC.}

\item{keepInterFile}{a logical value indicating if the intermediate 
processed files should be kept or not. The default is TRUE.}
}
\value{
The output PLINK binary files after QC.
}
\description{
Perform quality control on the genotype data.
}
\details{
The original PLINK files are implicitly processed by the following 
default steps: 
1.) Set all heterozygous alleles of SNPs on male chrX as missing;
2.) SNP missingness < 0.05 (before sample removal);
3.) Subject missingness < 0.02;  
4.) Remove subjects with |Fhet| >= 0.2;
5.) Reset paternal and maternal codes;
6.) SNP missingness < 0.02 (after sample removal);
7.) Remove SNPs with difference >= 0.02 of SNP missingness 
between cases and controls;
8.) Remove subjects or SNPs with Mendel errors for family based data.
9.) Remove chrX SNPs with missingness >= 0.05 in females.
(Optional, if no chrX data);
10.) Remove autosomal SNPs with HWE p < 10-6 in controls;
11.) Remove chrX SNPs with HWE p < 10-6 in female controls. 
(Optional, if no chrX data).
}
\examples{
 
## In the current working directory
bedFile <- system.file("extdata", "genoUpdatedData.bed", package="Gimpute")
bimFile <- system.file("extdata", "genoUpdatedData.bim", package="Gimpute") 
famFile <- system.file("extdata", "genoUpdatedData.fam", package="Gimpute")
system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
inputPrefix <- "genoUpdatedData" 
outputPrefix <- "2_13_removedSnpHweFemaleX"  
## Not run: Requires an executable program PLINK, e.g.
## plink <- "/home/tools/plink"
## genoQC(plink, inputPrefix, 
##        snpMissCutOffpre=0.05, 
##        sampleMissCutOff=0.02, 
##        Fhet=0.2, cutoffSubject, cutoffSNP,
##        snpMissCutOffpost=0.02, 
##        snpMissDifCutOff=0.02,
##        femaleChrXmissCutoff=0.05, 
##        pval4autoCtl=0.000001, 
##        pval4femaleXctl=0.000001, 
##        outputPrefix, keepInterFile=TRUE)
}
\references{
Schizophrenia Working Group of the Psychiatric Genomics, C. 
(2014). Biological insights from 108 schizophrenia-associated genetic loci. 
Nature 511(7510): 421-427.
}
\author{
Junfang Chen
}
