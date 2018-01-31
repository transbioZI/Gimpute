
##########################################   
##########################################	
#' Remove subjects with great autosomal heterozygosity deviation
#'
#' @description
#' Remove instances with great autosomal heterozygosity deviation. i.e. with |Fhet| >= 0.2.
#' This analysis will automatically skip haploid markers (male X and Y chromosome markers).

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param Fhet the cutoff of the autosomal heterozygosity deviation.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.

#' @return  The output PLINK format files after removing subjects with great autosomal heterozygosity deviation.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 

removedInstFhet <- function(plink, Fhet, inputPrefix, outputPrefix){ 

	system( paste0(plink, " --bfile ", inputPrefix, " --het --out ", outputPrefix) )
	#  F inbreeding coefficient estimate
	autoHet = read.table(file=paste0(outputPrefix, ".het"), header=T)  
	fhet = autoHet[, "F"]
	qc_data_fhet = autoHet[which(abs(fhet) >= Fhet), c(1, 2)]  
	## the individual IDs to be removed  
	write.table(qc_data_fhet, file=paste0(outputPrefix, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	##############
	# To remove certain individuals 
	system( paste0(plink, " --bfile ", inputPrefix, " --remove ", paste0(outputPrefix, ".txt"), " --make-bed --out ", outputPrefix) )
	system( paste0("rm ", outputPrefix, ".het") )
	system( paste0("rm ", outputPrefix, ".txt") )
 
}   