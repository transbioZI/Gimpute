
 

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
 


#' Get the outcome label of the genotype data
#'
#' @description
#' Get the group label from the PLINK FAM file.
#' 
#' @param inputPrefix the PLINK FAM file. 
 
#' @return  The group label of the genotype data: "control" or "case" or 
#' "caseControl" indicating both groups exist.
#' @details If the input FAM file also has missing outcomes, which are 
#' shown in the sixth column of FAM file as "0", then an error is given.

#' @export 
#' @author Junfang Chen 


getGroupLabel <- function(inputFAMfile){ 

	## check case control label: (1=unaff, 2=aff, 0=miss)
	fam <- read.table(file=inputFAMfile, stringsAsFactors=FALSE)
	groupIDs <- names(table(fam[,6]))

	if (length(groupIDs) == 1) {
		if (groupIDs == "1") {label <- "control"} 
		if (groupIDs == "2") {label <- "case"} 
		print(label)
	} else if (length(groupIDs) == 2) {
		label <- "caseControl"
		print(label)
	} else {
		print("ERROR: more than two groups!")
	}

	return(label)
}