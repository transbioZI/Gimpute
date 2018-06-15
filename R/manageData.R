
 

#' Rename PLINK binary files  
#'
#' @description
#' Rename a set of PLINK format files. Namely, .BED, .BIM and .FAM   
#' 
#' @param inputPrefix the prefix of the input PLINK format files. 
#' @param outputPrefix  the prefix of the output PLINK format files. 

#' @return  Renamed PLINK binary files. 
#' @export 

#' @author Junfang Chen <junfang.chen@zi-mannheim.de> 
#' @details The original input files can be retained using the action "copy"
#' or removed by using "move".

renamePlinkBFile <- function(inputPrefix, outputPrefix, action){

	if (action == "copy"){
		system(paste0("scp ", inputPrefix, ".bed", outputPrefix, ".bed"))
		system(paste0("scp ", inputPrefix, ".bim", outputPrefix, ".bim"))
		system(paste0("scp ", inputPrefix, ".fam", outputPrefix, ".fam"))
	} else if (action == "move"){
		system(paste0("mv ", inputPrefix, ".bed", outputPrefix, ".bed"))
		system(paste0("mv ", inputPrefix, ".bim", outputPrefix, ".bim"))
		system(paste0("mv ", inputPrefix, ".fam", outputPrefix, ".fam"))
	}
}
 

