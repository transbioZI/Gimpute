 
#' Remaming chromosomes
#'
#' @description
 
#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after renaming chromosomes.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 
##########################################   
########################################## replace chr name
## If the plink files contain names of the chromosomes X, Y, 
## XY (pseudo-autosomal region of X) and MT (mitochondrial) 
## replace these by the following numbers: 
## X by 23, Y by 24, XY by 25 and MT by 26.
## plink: --update-chr [filename] {chr col. number} {variant ID col.}
renamedSnpChr <- function(plink, inputPrefix, outputPrefix){

	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=F)
	whX = which(bim[,1]=="X")
	whY = which(bim[,1]=="Y")
	whXY = which(bim[,1]=="XY")
	whMT = which(bim[,1]=="MT")

 	tmpBIM = bim
 	tmpBIM[whX,1] <- 23
 	tmpBIM[whY,1] <- 24
 	tmpBIM[whXY,1] <- 25
 	tmpBIM[whMT,1] <- 26
 	
 	updateSNPchr = tmpBIM[c(whX, whY, whXY, whMT), c(2,1)] ## change to --> 1st SNP; 2nd chr names
  	write.table(updateSNPchr, file=paste0(outputPrefix, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
 	system( paste0(plink, " --bfile ", inputPrefix,  " --update-chr ", paste0(outputPrefix, ".txt"), " 2 1 --make-bed --out ", outputPrefix) )  
 	# system( paste0("rm ", outputPrefix, ".txt") )  
}

 