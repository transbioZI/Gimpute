
#' Remove monomorphic SNPs 
#'
#' @description
#' Remove monomorphic SNPs after QC and alignment to the imputation reference.

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input plink files.
#' @param outputPrefix the prefix of the output plink files (after removing monomorphic SNPs).
#' @param outputSNPs the output pure text file that stores the removed monomorphic SNPs.

#' @return  The output plink files (after removing monomorphic SNPs) and a pure text file with removed monomorphic SNPs.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 

removedMonoSnp <- function(plink, inputPrefix, outputPrefix, outputSNPs){  

		## input  
		bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=F)
		monoSNPs = bim[which(bim[,5]==0),2]  
		write.table(monoSNPs, file=outputSNPs, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 

		## exclude those not well aligned snps   
		system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", outputSNPs, " --make-bed --out ", outputPrefix)) 

}

