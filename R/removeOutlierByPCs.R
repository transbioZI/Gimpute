


 

  
 
##
######################################################
######################################################
#' Remove population outliers by principle components
#'
#' @description
#' This function decides if one would remove population outliers or not. If the outliers are necessary to be removed 
#' then one uses the eigenvalues from the first principle component as a criterion for finding out the outliers by defining a proper cutoff.
  

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param cutoff the cutoff that would distinguish the outliers from ordinary population. 
#' If it is null, then no outliers or not necessary to remove outliers.
#' If it has only a value, then one also has to define the cutoffSign as below.
#' If it has two values, then one doesn't have to define the cutoffSign.
#' @param cutoffSign the cutoff sign: 'greater' or 'smaller' that would determine if the outliers should be greater or smaller than the cutoff value.
#' if the cutoff score has two values, then no need to define the cutoffSign.
#' @param inputPC4subjFile the pure text file that stores all the subject IDs after QC and their corresponding eigenvalues of the first two principle components, if any.
#' @param outputPC4outlierFile the pure text file that stores the outlier IDs and their corresponding eigenvalues of the first two principle components, if any.
#' @param outputPCplotFile the plot file for visualizing the first two principle components of all subjects without population outliers, if any.

#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files..
#' @param outputSNPs the output pure text file that stores the removed monomorphic SNPs.

#' @return  The output PLINK format files after outlier removal. The output pure text file (if any) for storing removed outlier IDs and their corresponding PCs. 
#' The plot file (if any) for visualizing the first two principle components after outlier removal.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 
 
removeOutlierByPCs <- function(plink, gcta, inputPrefix, cutoff, cutoffSign, inputPC4subjFile, outputPC4outlierFile, outputPCplotFile, outputPrefix) {

	## if no outliers or no need to remove PC outliers. 
	if ( is.null(cutoff) ==TRUE ){ 
		## copy/rename all snp info updated plink files
		system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
		system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
		system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )
	} 
		
	else { 

		subjID_PCs =read.table(file=inputPC4subjFile, stringsAsFactors=FALSE)   
		if (length(cutoff) > 1) { ## if the outliers should be removed on both side of the cluster
				outliersPC1v1 = subjID_PCs[which(subjID_PCs[,3] <= cutoff[1]), ] ## detected by PC1
				outliersPC1v2 = subjID_PCs[which(subjID_PCs[,3] >= cutoff[2]), ] ## detected by PC1
				subjID_PCs4outlier = rbind(outliersPC1v1, outliersPC1v2)

			} else {
			  	  if (cutoffSign == "smaller"){ 
			  	  	subjID_PCs4outlier = subjID_PCs[which(subjID_PCs[,3] <= cutoff), ] ## detected by PC1
			  	  	} else if (cutoffSign == "greater"){
			  	  		subjID_PCs4outlier = subjID_PCs[which(subjID_PCs[,3] >= cutoff), ] ## detected by PC1
			  	  	 }
				}	 
	 	# str(subjID_PCs4outlier)
	 	subjID_PCs4outlierSorted = subjID_PCs4outlier[order(subjID_PCs4outlier[,3]), ] ## sorted by first PC.
	 	write.table(subjID_PCs4outlierSorted, file=outputPC4outlierFile, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
	 	# str(subjID_PCs4outlierSorted)
	 	subjID4outlierTmp  = subjID_PCs4outlierSorted[,1:2]
	 	str(subjID4outlierTmp)
	 	subjID4outlierTmpFile = "subjID4outlierTmp.txt"
		write.table(subjID4outlierTmp, file=subjID4outlierTmpFile, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
		system( paste0(plink, " --bfile ", inputPrefix, " --remove ", subjID4outlierTmpFile, " --make-bed --out ", outputPrefix) )
		system( paste0("rm ", subjID4outlierTmpFile) )

 
		## Plot first two PCs again
		outputPC4subjFiletmp = 'outputPC4subjFile.txt' ## PCs for the retained subjects 
		plotPCA4plink(gcta, inputPrefix=outputPrefix, outputPC4subjFiletmp, outputPCplotFile)
		system( paste0("rm ", outputPC4subjFiletmp))
 
	}	
}
 






##########################################   
##########################################
# #' Remove SNPs with missing values
# #'
# #' @description
# #' Remove SNPs with missingness of greater than a certain threshold before removing subjects.

# #' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
# #' @param snpMissCutOff the cutoff of the missingness for removing SNPs.
# #' @param inputPrefix the prefix of the input PLINK format files.
# #' @param outputPrefix the prefix of the output PLINK format files.

# #' @return  The output PLINK format files after removing SNPs with missing certain missing values.
# #' @export 

# #' @author Junfang Chen <junfang.chen3@gmail.com> 
# #' @examples 
 
# ## plink2 --bfile unpruned_data --make-king-table --king-cutoff 0.177 --king-table-filter 0.177 --make-bed --out pruned_data
# ## https://groups.google.com/forum/#!topic/plink2-users/F-b4XRF8CSc
  
# removedInstRelated <- function(plink, kinshipValue, inputPrefix, outputPrefix, outputPrefix.ID){ 

# 	# new step6. 6. Replace the paternal ID and maternal ID of instances  by the value zero  
# 	### special case for plink 
# 	# plink2 = '/data/noether/tools/plink/plink2'
# 	# system( paste0(plink2, " --bfile ", inputPrefix, " --make-king-table --king-cutoff 0.177  --king-table-filter 0.11 --make-bed --out ", outputPrefix) ) 
# 	## >> this doesn't work; 

# 	## copy/rename all snp info updated plink files
# 	system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
# 	system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
# 	system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )
# 	system( paste0("touch ", outputPrefix.ID, ".txt") )  ## empty

# }   
	 