

##########################################   
##########################################
#' Set haploid heterozygous SNPs as missing 
#'
#' @description
#' Set all heterozygous alleles of chromosome X SNPs in male as missing.

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.
#' @return  The output PLINK format files.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
  
setHeteroHaploMissing <- function(plink, inputPrefix, outputPrefix){
 
	system( paste0(plink, " --bfile ", inputPrefix, " --set-hh-missing --make-bed --out ", outputPrefix) )  
 
}  
