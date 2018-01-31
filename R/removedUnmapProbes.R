

## step7 
##########################################   
##########################################  outputPrefix >  remove unmapped 2 chipRef
#' Remove unmapped probes or SNPs
#'
#' @description
#' Remove probes or SNPs that are not mapped to the chip annotation information. This chip annotation file is defined in the configuration folder. 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param chipAnnoFile a pure text file that stores the chip annotation information (Affymetrix, Illumination, PsychChip and so on). 
#' This file can be found http://www.well.ox.ac.uk/~wrayner/strand/.
#' @param outputSNPunmapFile a pure text file that stores the probe/SNP IDs, which are not mapped to the chip annotation file.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after removing unmapped probe IDs.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

## re-generate chip annotation file/prepare chipAnnoFile 
## generate chip annotation file in a correct format
 

removedUnmapProbes <- function(plink, inputPrefix, chipAnnoFile, outputPrefix, outputSNPunmapFile){

 	annoFile = "chipAnnoRefb37.txt"
 	if (chipType=="affymetrix") { 
 		prepareChipAnnoFile4affymetrix(inputFile=chipAnnoFile, outputFile=annoFile)
 	} else if (chipType=="illumina"){ 
 		prepareChipAnnoFile4Illumina(inputFile=chipAnnoFile, outputFile=annoFile)
 	  }
 
	## find the overlapping
	chipAnno = read.table(file=annoFile, header=TRUE, stringsAsFactors=FALSE)
	# system( paste0("rm ", annoFile))
	## check the overlapping SNPs
	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
	interSNPs = intersect(bim[,2], chipAnno[,"chipSnpID"])  
	unmapped = setdiff(bim[,2], interSNPs)
	write.table(unmapped, file=outputSNPunmapFile, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
	system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", outputSNPunmapFile, " --make-bed --out ", outputPrefix) )
 
}
 