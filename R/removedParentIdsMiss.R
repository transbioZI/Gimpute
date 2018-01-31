
##########################################   
##########################################
#' Reset paternal and maternal codes  
#'
#' @description
#' Reset paternal and maternal codes of non-founders if parents not present. 
#' Replace the paternal ID and maternal ID of instances (childs) by the
#' value zero if the paternal ID and the maternal ID do not belong to any
#' instance (parent) with the same family ID as the child. 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPrefix the prefix of the output PLINK format files.

#' @return  The output PLINK format files.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

# By default, if parental IDs are provided for a sample, they are not treated as a founder even if neither parent is in the dataset. 
# With no modifiers, --make-founders clears both parental IDs whenever at least one parent is not in the dataset, 
# and the affected samples are now considered founders. 
# This is done by invoking PLINK command '--make-founders'.
 
removedParentIdsMiss <- function(plink, inputPrefix, outputPrefix){ 

	# Remove the parent IDs which do not belong to instances.
	system( paste0(plink, " --bfile ", inputPrefix, " --make-founders require-2-missing --make-bed --out ", outputPrefix) ) 
 
}   
 