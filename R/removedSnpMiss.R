
  


##########################################   
##########################################
#' Remove SNPs with missing values
#'
#' @description
#' Remove SNPs with missingness of greater than a certain threshold before/after removing subjects.

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param snpMissCutOff the cutoff of the missingness for removing SNPs.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.

#' @return  The output PLINK format files after removing SNPs with pre-defined missing values.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 
removedSnpMiss <- function(plink, snpMissCutOff, inputPrefix, outputPrefix){
 
	system( paste0(plink, " --bfile ", inputPrefix, " --geno ", snpMissCutOff, " --make-bed --out ", outputPrefix) )  
	 
}  
   
  