





 




##########################################   
########################################## 
#' prepare Illumina chip annotation file 
#'
#' @description
#' Prepare Illumina chip annotation file into the format of interest.

#' @param inputFile an input pure text file that shows the chip annotation file which can be downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.
#' @param outputFile an output pure text file that stores the chip annotation information in used-defined format (Affymetrix, Illumination, PsychChip and so on). 
 
#' @return  a pure text file that stores the prepared chip annotation information in used-defined format. 
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 


## prepareChipAnnoFile4Illumina

# colnames(chipAnnoRefraw) = c("chipSnpID", "chr", "pos", "match2genom", "strand", "allele")
## the column name must be fixed
## chipSnpID is set in the way so that it's comparable to Affymetrix chip

## 1. don't have to remove "cnvi", since we will remove this type of SNPs in .raw plink files
## 2. don't have to remove 'unknown' SNPs, since they won't map to our raw plink files.
 
prepareChipAnnoFile4Illumina <- function(inputFile, outputFile){

	chipAnnoRefraw = read.table(file=inputFile, stringsAsFactors=FALSE)
	# print(colnames(chipAnnoRefraw))
	## c("chipSnpID", "chr", "pos", "match2genom", "strand", "allele")
	chipAnnoRefraw <- chipAnnoRefraw[,-c(4, 6)]
	colnames(chipAnnoRefraw) = c("chipSnpID", "chr", "pos", "strand")

	## only see 3 different cases (if 25--> XY)
	whX = which(chipAnnoRefraw[,"chr"]=="X")
	whY = which(chipAnnoRefraw[,"chr"]=="Y")
	whMT = which(chipAnnoRefraw[,"chr"]=="MT")

	chipAnnoRefraw[whX,"chr"] = 23 
	chipAnnoRefraw[whY,"chr"] = 24
	chipAnnoRefraw[whMT,"chr"] = 26 

	chipAnnoNew = chipAnnoRefraw
	write.table(chipAnnoNew, file=outputFile, quote=F, row.names=F, col.names=TRUE, eol="\r\n", sep="\t")

}

  