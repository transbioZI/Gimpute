
##########################################   
##########################################
#' Remove probes or SNPs
#'
#' @description
#' Remove probes or SNPs that may be duplicated, or with unexpected probe names, which should be defined in advance. 
#' For instance, one can define these probes in the configuration folder. 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param excludedProbeIdsFile a pure text file that stores the probe IDs, which should be removed. 
#' If it is null, then copy and paste the input PLINK files from the last step.

#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after removing unwanted probe IDs.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 


## step 6
## Remove the excluded probes of the chip and check that there are not two probes with the same name afterwards
## > remove AFFX, cnvi, or similar snps
 
removedExclProbe <- function(plink, inputPrefix, excludedProbeIdsFile, outputPrefix){
 
	if (!is.null(excludedProbeIdsFile)) {

		system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", excludedProbeIdsFile, " --make-bed --out ", outputPrefix) )
		## remove .txt
		system( paste0("rm ", excludedProbeIdsFile) )
	} else { 
		## copy/rename all snp info updated plink files
		system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
		system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
		system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )
	    }

}
  