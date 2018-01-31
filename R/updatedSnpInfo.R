
##########################################   
##########################################  updateSNPinfo
#' Update the probes or SNPs information
#'
#' @description
#' Remove SNP information including SNP name, genomic position, chromosome location and the strand information, with the help of the chip annotation file. 
#' This chip annotation file is defined in the configuration folder. 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param chipAnnoFile a pure text file that stores the chip annotation information (Affymetrix, Illumination, PsychChip and so on). 
#' This file can be found http://www.well.ox.ac.uk/~wrayner/strand/.
#' @param chipType a string name that defines the type of the chip annotation file: 'Illumina', 'affymetrix' or 'psychChip'.
#' @param outputPrefix the prefix of the output PLINK format files.
 
#' @return  The output PLINK format files after updating SNP information.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 
updatedSnpInfo <- function(plink, inputPrefix, chipAnnoFile, chipType, outputPrefix){
 	
 	annoFile = "chipAnnoRefb37.txt"
 	if (chipType=="affymetrix") { 
 		prepareChipAnnoFile4affymetrix(inputFile=chipAnnoFile, outputFile=annoFile)
 	} else if (chipType=="illumina"){ 
 		prepareChipAnnoFile4Illumina(inputFile=chipAnnoFile, outputFile=annoFile)
 	  }
 
	## find the overlapping
	chipAnno = read.table(file=annoFile, header=TRUE, stringsAsFactors=FALSE)
	system( paste0("rm ", annoFile)) ## not used anymore

	# ## cbind (combine) bim and chip annotation files 
	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
	interSNPs = intersect(bim[,2], chipAnno[,"chipSnpID"]) 
	bimV1 = bim[is.element(bim[,2], interSNPs),]
	chipAnnoV1 = chipAnno[is.element(chipAnno[,"chipSnpID"], interSNPs),]
	chipAnnoV1sort = chipAnnoV1[match(bimV1[,2], chipAnnoV1[,"chipSnpID"]),]
	comV2 = cbind(bimV1, chipAnnoV1sort) 
 
	## Update main info  
	if (chipType=="affymetrix") { 
 			updateSNP2rs = subset(comV2, select=c(V2, rsID))
			updateSNPchr = subset(comV2, select=c(rsID, chr))
			updateSNPpos = subset(comV2, select=c(rsID, pos)) 
			updateSNPbackward = comV2[which(comV2[,"strand"]=="-"), "rsID"]  ## note strand 
 	} else if (chipType=="illumina"){ 
	    	updateSNP2rs = subset(comV2, select=c(V2, V2)) ## no need to change but for consistency with Affymtrix
			updateSNPchr = subset(comV2, select=c(V2, chr))
			updateSNPpos = subset(comV2, select=c(V2, pos)) 
			updateSNPbackward = comV2[which(comV2[,"strand"]=="-"), "V2"]  ## note strand 
 	 }

 	inputPrefix.rs = "1_09.updatedSnp2rs" ## tobeRemoved
	inputPrefix.chr = "1_09.updatedSnpchr" ## tobeRemoved
	inputPrefix.pos = "1_09.updatedSnppos" ## tobeRemoved
	inputPrefix.strand = "1_09.updatedSnpstrand" ## tobeRemoved
	 
	write.table(updateSNP2rs, file=paste0(inputPrefix.rs, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	write.table(updateSNPchr, file=paste0(inputPrefix.chr, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	write.table(updateSNPpos, file=paste0(inputPrefix.pos, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	write.table(updateSNPbackward, file=paste0(inputPrefix.strand, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	
	## update rs, chr, and pos one by one 	## flip to the forward strand
	system( paste0(plink, " --bfile ", inputPrefix,    " --update-name ", paste0(inputPrefix.rs, ".txt"), " 2 1  --make-bed --out ", inputPrefix.rs) )  
	system( paste0(plink, " --bfile ", inputPrefix.rs,  " --update-chr ", paste0(inputPrefix.chr, ".txt"), " 2 1 --make-bed --out ", inputPrefix.chr) )  
	system( paste0(plink, " --bfile ", inputPrefix.chr, " --update-map ", paste0(inputPrefix.pos, ".txt"), " 2 1 --make-bed --out ", inputPrefix.pos) )   
	system( paste0(plink, " --bfile ", inputPrefix.pos, " --flip ", paste0(inputPrefix.strand, ".txt"), " --make-bed --out ", inputPrefix.strand) )  

	## copy/rename all snp info updated plink files
	system( paste0("cp ", inputPrefix.strand, ".bed ", outputPrefix, ".bed") )
	system( paste0("cp ", inputPrefix.strand, ".bim ", outputPrefix, ".bim") )
	system( paste0("cp ", inputPrefix.strand, ".fam ", outputPrefix, ".fam") )
 
 	## remove all tmp files (rs, chr, pos)
	system( paste0("rm ", inputPrefix.rs, ".*") )  
	system( paste0("rm ", inputPrefix.chr, ".*") )  
	system( paste0("rm ", inputPrefix.pos, ".*") )  
 	system( paste0("rm ", inputPrefix.strand, ".*"))

}

