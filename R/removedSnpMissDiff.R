

##########################################   
##########################################  
#' Remove SNPs with difference in SNP missingness between cases and controls. 
#'
#' @description
#' Remove SNPs with difference in SNP missingness between cases and controls. 
#' To test for differential call rates between cases and controls for each SNP
#' This only works for the genotype data set with cases and controls. 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param snpMissDifCutOff the cutoff of the difference in missingness between cases and controls. 
#' @param outputPrefix the prefix of the output PLINK format files.
#' @param caseControl a logical value indicating whether the data set contains both case-control subjects.

#' @return The output PLINK format files.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 
# SNP missingness < 0.02 (after sample removal);
removedSnpMissDiff <- function(plink, inputPrefix, snpMissDifCutOff, outputPrefix, caseControl){

	if (caseControl==FALSE){  ## this is only for the control data set

			system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
			system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
			system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )

		} else if (caseControl==TRUE) {

  			outputPrefix.tmp = paste0(outputPrefix, "tmp")
			system (paste0(plink, " --bfile ", inputPrefix, " --test-missing --out ", outputPrefix.tmp) ) 
			# Writing case/control missingness test to [ *.missing ]  
			###### add R command to compute differential call rates between cases and controls for each SNP 
			ccmissing = read.table(file=paste0(outputPrefix.tmp, ".missing"), header=T, sep="")  
			whmiss = which(abs(ccmissing[, "F_MISS_A"] - ccmissing[, "F_MISS_U"]) >= snpMissDifCutOff) 
			SNPmissCC = ccmissing[whmiss, "SNP"]
			SNPdifCallrate = paste0(outputPrefix, ".txt")
			write.table(SNPmissCC, file=SNPdifCallrate, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

			########## exclude SNPs 
			system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", SNPdifCallrate, " --make-bed --out ", outputPrefix) )  
			system( paste0("rm ", outputPrefix.tmp, ".*") )
			system( paste0("rm ", SNPdifCallrate) )
		}
}   
 