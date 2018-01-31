
  
##########################################   
##########################################  
#' Remove duplicated sample IDs
#'
#' @description
#' Remove duplicated sample IDs, which should be defined in advance. For instance, one can define these duplicated sample IDs in the configuration folder. 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param dupSampleIDFile a pure text file that stores the duplicated sample IDs, which should be removed. 
#' If it is null, then copy and paste the input PLINK files from the last step.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after removing duplicated sample IDs. 
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

#### clean/prepare the raw plink file  
# in the FAM file
# colnames(fam) = c("famid", "id", "paternalID", "maternalID", "sex", "phenotype") 
# check the followings
# 1. any .dup IDs (both famid and id), if so remove IDs with .dup
# 2. any duplicated IDs (only for id), if so remove them
# 3. correct sex format?  (1=male; 2=female; other=unknown) 
# 4. correct phenotype ? Assuming a disease phenotype (1=unaff, 2=aff, 0=miss) 

## remove .dup individuals 
 
removeDupID <- function(plink, dupSampleIDFile, inputPrefix, outputPrefix){

	if (!is.null(dupSampleIDFile)){

		famv0 = read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE)  
		dupIDs = read.table(file=dupSampleIDFile, stringsAsFactors=FALSE)  
		## write into plink format .txt for removal
		dupIDs4plink = famv0[is.element(famv0[,2], dupIDs[,1]), 1:2] 
		dupIDsfn = paste0(outputPrefix, ".txt")
		write.table(dupIDs4plink, file=dupIDsfn, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 

		system( paste0(plink, " --bfile ", inputPrefix, " --remove ", dupIDsfn, " --make-bed --out ", outputPrefix) )
		system( paste0("rm ", dupSampleIDFile) )
		system( paste0("rm ", dupIDsfn) )

		## remove the raw genotype plink files
		system( paste0("rm ", inputPrefix, ".*") )

	} else { 
		## copy/rename all snp info updated plink files
		system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
		system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
		system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )

		## remove the raw genotype plink files
		system( paste0("rm ", inputPrefix, ".*") )
	  
	    } 
}

