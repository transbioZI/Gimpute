
  
##########################################   
########################################## 
#' Remove samples with improper ancestry
#'
#' @description
#' Remove samples with improper ancestry or keep samples with your own purpose.

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param metaDataFile a pure text file that stores the meta information of the samples.  
#' @param ancestrySymbol an indicator that shows the symbol of targeted ancestry. Such as 'EA' stands for the European, 'AA' for African American. 
#' If it is null, then all samples are selected.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after ancestry check.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

## step 5
## wrong or improper ancestry instances define in the meta data >> 
## ancestrySymbol == 'EA' In this case, only keep EA 
## If you want to impute your data set using Multi-Population Reference Panels, then you don't have to exclude improper ancestry.
 
removedWrongAnceInst <- function(plink, inputPrefix, metaDataFile, ancestrySymbol, outputPrefix){

	if (is.null(ancestrySymbol)) {
		## copy/rename all snp info updated plink files
		system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
		system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
		system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )

	} else {
		metaData = read.table(metaDataFile, stringsAsFactors=FALSE, header=TRUE) 
		ids = metaData[which(metaData[,"ance"]==ancestrySymbol),"IID"]
		fam = read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
		famEA = fam[is.element(fam[,2], ids), 1:2]
		plinkFormatDat = famEA
		write.table(plinkFormatDat, file=paste0(outputPrefix,".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")  

		system( paste0(plink, " --bfile ", inputPrefix, " --keep ", paste0(outputPrefix,".txt"), " --make-bed --out ", outputPrefix) )
		system( paste0("rm ",  outputPrefix, ".txt") )

	}	
	
}
