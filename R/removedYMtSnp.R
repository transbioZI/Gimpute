
 
  
##########################################   
##########################################  
#' Remove SNPs on the chromosome Y and mitochondria
#'
#' @description
#' Remove SNPs on the chromosome Y and mitochondria (in this case, chromosome==24 & chromosome==26 are removed). 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after splitting chromosome X into pseudoautosomal region and non-pseudoautosomal region.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 
 
removedYMtSnp <- function(plink, inputPrefix, outputPrefix){

 	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
 	snpChrY =  bim[which(bim[,1]== 24), 2] ##  
	snpChrMt = bim[which(bim[,1] == 26), 2] 
	snpChrYMt = c(snpChrY, snpChrMt) 
	write.table(snpChrYMt, file=paste0(outputPrefix, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", paste0(outputPrefix, ".txt"), " --make-bed --out ", outputPrefix) )  
	system( paste0("rm ", paste0(outputPrefix, ".txt")) )
}

  