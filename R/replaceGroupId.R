

##########################################   
##########################################  replaceGroupId.R
#' Replace group and geneder presentation into proper PLINK format
#'
#' @description
#' Replace group and gender presentation into proper and consistent PLINK format.

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param metaDataFile a pure text file that stores the meta information of the samples.  
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after presenting gender/group/outcome information into proper PLINK format. 
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

## 1. find the shared Ids between plink files and metadata file 
##    a. update the GroupId to the right format
## 	  b. Group should be 1 and 2. (1=unaff, 2=aff, 0=miss) instead of case 1 control 0
# metaData: data.frame with at least 2 columns, 1st column is sample IDs, 2nd is the group, case 1 control 0.
  

## 1. check the phenotypes, if there is any sample without phenotype, then remove it.
## for some datasets, there is no such case
# metaData: data.frame with at least 2 columns, 1st column is sample IDs, 2nd is the group, case 1 control 0.
 

replaceGroupId <- function(plink, inputPrefix, metaDataFile, outputPrefix){
 
	fam = read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
	metaData = read.table(metaDataFile, stringsAsFactors=FALSE, header=TRUE) 
	## to make sure only compare with the same set of plIndID --> IID
	interIDs = intersect(fam[,2], metaData[,"IID"]) 	# the shared IDs 
	metadataSub = metaData[is.element(metaData[,"IID"], interIDs),]

	## one must keep the order of .fam unchanged!!
	metadataSubsort = metadataSub[match(fam[,2], metadataSub[,"IID"]),]
	missIDsIndex = which(is.na(metadataSubsort[,1])==TRUE) ## IID is in the 1st col; 
	fam[,6] = metadataSubsort[,"group"] + 1
	fam[missIDsIndex, 6] = -9  ## keep the IDs of no missing group info as -9
 	newPheno = fam[,c(1,2,6)]  ##  a pheno file that contains 3 columns (one row per individual)
 	write.table(newPheno, file="pheno.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")  

	## Alternate phenotype files
	system( paste0(plink, " --bfile ", inputPrefix, " --pheno pheno.txt --make-bed --out ", outputPrefix) )
	system( "rm pheno.txt" )

}
