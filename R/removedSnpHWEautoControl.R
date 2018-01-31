

##########################################   
##########################################
#' Hardy weinberg equilibrium test for autosomal SNPs in controls.
#'
#' @description
#' Remove autosomal SNPs in controls.
 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param pval the p-value cutoff for controlling HWE test in control subjects. Only autosomal SNPs are tested here. 
#' @param outputFile_pVal the output pure text file that stores autosomal SNPs and their sorted HWE p-values.
#' @param outputSNPfile the output pure text file that stores the removed SNPs.
#' @param outputPrefix the prefix of the output PLINK format files.

#' @return  The output PLINK format files.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 


##  

removedSnpHWEautoControl <- function(plink, inputPrefix, pval, outputFile_pVal, outputSNPfile, outputPrefix){ 
 
	outputPrefix.tmp = paste0(outputPrefix, "tmp")
	system(  paste0(plink, " --bfile ", inputPrefix, " --filter-controls --hardy --autosome --make-bed --out ", outputPrefix.tmp) )  
	###  read HWE p values 
	hweCheck <- read.table(file=paste0(outputPrefix.tmp, ".hwe"), header=T, stringsAsFactors=F) 
	hweControl = hweCheck[which(hweCheck$TEST == "UNAFF"), ] # ## for controls 
	snpHweValuesAutoCt = subset(hweControl, select=c(SNP, P))
	snpHweValuesAutoCt = snpHweValuesAutoCt[order(snpHweValuesAutoCt[,'P']),]
	# head(snpHweValuesAutoCt)
	write.table(snpHweValuesAutoCt, file=outputFile_pVal, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	removedSNPs_hweControl = hweControl[which(hweControl$P <= pval), "SNP"]
	write.table(removedSNPs_hweControl, file=outputSNPfile, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

	########## exclude SNPs 
	system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", outputSNPfile, " --make-bed --out ", outputPrefix) ) 
	system( paste0("rm ", outputPrefix.tmp, ".*") )

}

