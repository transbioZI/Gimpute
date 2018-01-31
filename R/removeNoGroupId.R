

##########################################   
########################################## removeNoGroupId.R
#' Remove samples without group information
#'
#' @description
#' Remove samples without group/outcome/phenotype information, which is coded as -9 in PLINK .fam file.

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param metaDataFile a pure text file that stores the meta information of the samples.  
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after removing samples without group information/IDs.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

removeNoGroupId <- function(plink, inputPrefix, outputPrefix){
 
	fam = read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
	## 1. check if any sample without phenotypes
	noGroupIds = fam[which(fam[,6] == -9), 1:2]
	print(dim(noGroupIds))
	## if any then remove and afterwards add phenotypes
	noGroupIdsfn = paste0(outputPrefix, ".txt")
	write.table(noGroupIds, file=noGroupIdsfn, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")   
	system( paste0(plink, " --bfile ", inputPrefix, " --remove ", noGroupIdsfn, " --make-bed --out ", outputPrefix) )  

	system( paste0("rm ",  noGroupIdsfn) )
}
 
