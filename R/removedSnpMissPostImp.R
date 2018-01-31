##########################################################################
########################################################################## removedSnpMissPostImp.R
#' Remove SNPs after post imputation  
#'
#' @description
#' Remove SNPs which have a non missing value for less than a predefined number of instances.    
#' 
#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param missCutoff  the cutoff of the least number of instances for a SNP that is not missing.
#' @param snpWithManyMissSNPfile the output file of SNPs with pre-defined missing values.
#' @param outputPrefix  the prefix of the PLINK format files. 

#' @return  The PLINK format files after post imputation quality control and a pure text file contains SNPs with pre-defined missing values.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples   


removedSnpMissPostImp <- function(plink, inputPrefix, missCutoff, snpWithManyMissSNPfile, outputPrefix){ 

	## get the missing info 
	system(paste0(plink, " --bfile ", inputPrefix, " --missing --out ", inputPrefix)) 

	missSNPinfo = read.table(paste0(inputPrefix, ".lmiss"), stringsAsFactors=F, h=T)
	missSNPinfo[,6] <- missSNPinfo[,"N_GENO"] - missSNPinfo[,"N_MISS"] 
	snpWithManyMissSNPs <- missSNPinfo[which(missSNPinfo[,6] < missCutoff), "SNP"] 
	write.table(snpWithManyMissSNPs, file=snpWithManyMissSNPfile, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
	system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", snpWithManyMissSNPfile, " --make-bed --out ", outputPrefix) )
	system( "rm *.imiss *.lmiss *.log") 

}
# 

 
 