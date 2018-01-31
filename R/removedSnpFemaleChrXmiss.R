

##########################################   
##########################################  removedSnpFemaleChrXmiss.R 
#' remove chromosome X SNPs in females
#'
#' @description
#' #.Remove chrX SNPs with a pre-defined cutoff for missingness in females.

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param femaleChrXmissCutoff the cutoff of the missingness in female chromosome X SNPs.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.

#' @return  The output PLINK format files.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 


removedSnpFemaleChrXmiss <- function(plink, femaleChrXmissCutoff, inputPrefix, outputPrefix){ 
 
  	### additional QC (female-chrX SNPs, missingness ok?)  
 	outputPrefix.tmp1 = paste0(outputPrefix, "tmp1")
	outputPrefix.tmp2 = paste0(outputPrefix, "tmp2")
	system( paste0(plink, " --bfile ", inputPrefix, " --filter-females --chr 23 --make-bed --out ", outputPrefix.tmp1) )
	system( paste0(plink, " --bfile ", inputPrefix, " --filter-females --chr 23 --geno ", femaleChrXmissCutoff, " --make-bed --out ", outputPrefix.tmp2) )

	 ## check if equal  
	femaleChrXorig = read.table(paste0(outputPrefix.tmp1, ".bim"), stringsAsFactors=FALSE) 
	femaleChrXMiss = read.table(paste0(outputPrefix.tmp2, ".bim"), stringsAsFactors=FALSE)  
	snps2removed = setdiff(femaleChrXorig[,2], femaleChrXMiss[,2]) 
	snps2removedfile = paste0(outputPrefix, ".txt")
	write.table(snps2removed, file=snps2removedfile, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 

	system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", snps2removedfile, " --make-bed --out ", outputPrefix) )
	system( paste0("rm ", outputPrefix.tmp1, ".*") )
	system( paste0("rm ", outputPrefix.tmp2, ".*") )
	system( paste0("rm ", snps2removedfile) ) 

}

