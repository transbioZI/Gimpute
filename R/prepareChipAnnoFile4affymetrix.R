


##########################################   
##########################################   
#' prepare Affymetrix chip annotation file 
#'
#' @description
#' Prepare Affymetrix chip annotation file into the format of interest.

#' @param inputFile an input pure text file that shows the chip annotation file which can be downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.
#' @param outputFile an output pure text file that stores the chip annotation information in used-defined format (Affymetrix, Illumination, PsychChip and so on). 
 
#' @return  a pure text file that stores the prepared chip annotation information in used-defined format. 
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

prepareChipAnnoFile4affymetrix <- function(inputFile, outputFile){

	inputNew = paste0(inputFile, "New")
	system( paste0("sed 1d ", inputFile, " > ", inputNew) )
	chipAnnoRefraw = read.table(file=inputNew, stringsAsFactors=FALSE) 
	# colnames(chipAnnoRefraw) = c("chipSnpID", "chr", "pos", "strand") ## for Illumina
	colnames(chipAnnoRefraw) = c("chipSnpID", "rsID", "chr", "pos", "strand")

	## remove SNPs with strange strand 
	whUnknown = which(chipAnnoRefraw[,"strand"]=="---") 
	chipAnnoRefraw2 = chipAnnoRefraw[-whUnknown,]

	## only see 3 different cases (if 25--> XY)
	whX = which(chipAnnoRefraw2[,"chr"]=="X")
	whY = which(chipAnnoRefraw2[,"chr"]=="Y")
	whMT = which(chipAnnoRefraw2[,"chr"]=="MT")

	chipAnnoRefraw2[whX,"chr"] = 23 
	chipAnnoRefraw2[whY,"chr"] = 24
	chipAnnoRefraw2[whMT,"chr"] = 26
  
  	whAFF = grep("AFFX-SNP", chipAnnoRefraw2[,"chipSnpID"])
	# str(whAFF)
	chipAnnoRefraw3 = chipAnnoRefraw2[-whAFF,]
	write.table(chipAnnoRefraw3, file=outputFile, quote=F, row.names=F, col.names=TRUE, eol="\r\n", sep="\t")

}
 
 