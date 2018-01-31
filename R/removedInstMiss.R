


##########################################   
########################################## 
#' Remove subjects with missing values
#'
#' @description
#' Remove Subjects or instances with missingness of greater than a certain threshold.
#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param sampleMissCutOff the cutoff of the missingness for removing subjects/instances.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.
#' @return  The output PLINK format files after removing subjects with pre-defined missing values.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

 
removedInstMiss <- function(plink, sampleMissCutOff, inputPrefix, outputPrefix){
 
	system( paste0(plink, " --bfile ", inputPrefix, " --mind ", sampleMissCutOff, " --make-bed --out ", outputPrefix) ) 
	system( paste0("rm ", outputPrefix, ".irem") )

}   
