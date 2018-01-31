

##########################################   
########################################## 
#' Split chromosome X into pseudoautosomal region and non-pseudoautosomal region
#'
#' @description
#' Split chromosome X into pseudoautosomal region and non-pseudoautosomal region.

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after splitting chromosome X into pseudoautosomal region and non-pseudoautosomal region.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 

## step 10 

changedXyChr <- function(plink, inputPrefix, outputPrefix){
 	## split X chr into PAR(chr25) and non-PAR (chr23)
 	system( paste0(plink, " --bfile ", inputPrefix, "  --split-x hg19 --make-bed --out ", outputPrefix ))
}

