
##########################################   
########################################## 
#' Remove duplicated probes or SNPs
#'
#' @description
#' Remove duplicated probes or SNPs that have same rs-names but different SNP-A IDs found in chip annotation information (only for Affymetrix chip). 
#' This chip annotation file is defined in the configuration folder. 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param chipAnnoFile a pure text file that stores the chip annotation information (Affymetrix, Illumination, PsychChip and so on). 
#' This file can be found http://www.well.ox.ac.uk/~wrayner/strand/.
#' @param chipType a string name that defines the type of the chip annotation file: 'Illumina', 'affymetrix' or 'psychChip'.
#' @param outputSNPdupFile a pure text file that stores the duplicated probe/SNP IDs, which can be found with the help of the chip annotation file.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after removing duplicated probe IDs.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 
## remove double SNPs, dif SNP-A IDs but same rs-names 


removedDoubleProbes <- function(plink, inputPrefix, chipAnnoFile, chipType, outputSNPdupFile, outputPrefix){
 	
 	annoFile = "chipAnnoRefb37.txt"
 	if (chipType=="affymetrix") { 
 		prepareChipAnnoFile4affymetrix(inputFile=chipAnnoFile, outputFile=annoFile)
 	} else if (chipType=="illumina"){ 
 		prepareChipAnnoFile4Illumina(inputFile=chipAnnoFile, outputFile=annoFile)
 	  }
 
	## find the overlapping
	chipAnno = read.table(file=annoFile, header=TRUE, stringsAsFactors=FALSE) 

	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=F) 
	## cbind (combine) bim and chip annotation files  
	chipAnnoV1 = chipAnno[is.element(chipAnno[,"chipSnpID"], bim[,2]),]
	chipAnnoV1sort = chipAnnoV1[match(bim[,2], chipAnnoV1[,"chipSnpID"]),]
	comb = cbind(bim, chipAnnoV1sort)
	############# remove SNPs which have a duplicated rs-name or position (i.e. bp and chr) in this file
	## remove SNPs with duplicated position first
	chrNames = names(table(comb[,1]))
	dupPos = c()
	for (i in chrNames) { 
		print(i)
		subData = comb[which(comb[,1]==i), ]
  		subDup = subData[duplicated(subData[,"pos"]) | duplicated(subData[,"pos"], fromLast=TRUE), ] 
 		dupPos = rbind(dupPos, subDup)
	} 
	# print(dupPos)	
	snpWithdupPos = dupPos[,"chipSnpID"]
	## return all the duplicated rs-names (not only either one)

	if (chipType=="affymetrix") { 
 		snpdup = comb[duplicated(comb[,"rsID"]) | duplicated(comb[,"rsID"], fromLast=TRUE), "chipSnpID"] 
 	} else if (chipType=="illumina"){ 
 		snpdup = comb[duplicated(comb[,"V2"]) | duplicated(comb[,"V2"], fromLast=TRUE), "chipSnpID"]  
 	 }

 	allDupSNPs = c(snpWithdupPos, snpdup) 
 	allDupSNPs = unique(allDupSNPs)
	write.table(allDupSNPs, file=outputSNPdupFile, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
	cmd = paste0(plink, " --bfile ", inputPrefix, " --exclude ", outputSNPdupFile, " --make-bed --out ", outputPrefix)  
	system(cmd)

}
