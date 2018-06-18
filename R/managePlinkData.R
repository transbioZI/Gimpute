
 

#' Rename PLINK binary files  
#'
#' @description
#' Rename a set of PLINK binary files (.BED, .BIM and .FAM).   
#' 
#' @param inputPrefix the prefix of the input PLINK binary files. 
#' @param outputPrefix  the prefix of the output PLINK binary files. 

#' @return  Renamed PLINK binary files. 
#' @details The original input files can be retained using the action "copy"
#' or removed by using "move".

#' @export 
#' @author Junfang Chen 


renamePlinkBFile <- function(inputPrefix, outputPrefix, action){

	if (action == "copy"){
		system(paste0("scp ", inputPrefix, ".bed ", outputPrefix, ".bed"))
		system(paste0("scp ", inputPrefix, ".bim ", outputPrefix, ".bim"))
		system(paste0("scp ", inputPrefix, ".fam ", outputPrefix, ".fam"))
	} else if (action == "move"){
		system(paste0("mv ", inputPrefix, ".bed ", outputPrefix, ".bed"))
		system(paste0("mv ", inputPrefix, ".bim ", outputPrefix, ".bim"))
		system(paste0("mv ", inputPrefix, ".fam ", outputPrefix, ".fam"))
	}
}
 

