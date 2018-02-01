
 
 
##########################################   
##########################################
#' Remove male subjects with haploid heterozygous SNPs
#'
#' @description
#' Determine the frequency of male subjects that have heterozygous SNPs on chromosome X and a reasonable cutoff to remove those affect males.


#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param hhSubjCutOff the cutoff for removing male subjects with haploid heterozygous SNPs on the chromosome X.
#' @param outputPrefix the prefix of the output PLINK format files.
#' @param outputFile_subjHetFreqAll the output pure text file that stores male subjects that have heterozygous SNPs with their frequency (if any), i.e. the number of .hh SNPs in this male. 
#' @param outputFile_subjHetFreqRetained the output pure text file that stores male subjects that have heterozygous SNPs with their frequency after 'improper' subject removal
#' (if any).
#' @param outputFile_SNPhhFreqAll the output pure text file that stores all heterozygous SNPs with their frequency (the number of males for this SNP), if any. 

#' @return  The output PLINK format files and two pure text files (if any) with heterozygous SNPs and their respective frequency. 
#' @export 


removedMaleHetX <- function(plink, inputPrefix, hhSubjCutOff, outputPrefix, outputFile_subjHetFreqAll, outputFile_subjHetFreqRetained, outputFile_SNPhhFreqAll){		

			## just to get .hh file and .fam file 
			system( paste0(plink, " --bfile ", inputPrefix, " --chr 23 --filter-males --make-bed --out male23nonPAR" ))

			if ( is.na(file.size("male23nonPAR.hh")) ){  
				## copy/rename all snp info updated plink files
				system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
				system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
				system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") ) 
				system(paste0('touch ', outputFile_subjHetFreqAll) )
				system(paste0('touch ', outputFile_subjHetFreqRetained) )
	  
		    } else { 	

					# .hh (heterozygous haploid and nonmale Y chromosome call list)
					# Produced automatically when the input data contains heterozygous calls 
					# where they should not be possible (haploid chromosomes, male X/Y),
					# or there are nonmissing calls for nonmales on the Y chromosome.

					# A text file with one line per error (sorted primarily by variant ID, secondarily by sample ID) with the following three fields:
					# Family ID  Within-family ID Variant ID

				hh = read.table("male23nonPAR.hh", stringsAsFactors=F)
				fam = read.table("male23nonPAR.fam", stringsAsFactors=F)
								
				hetInstFreq = table(hh[,2]) 
				# head(sort(hetInstFreq,decreasing=T))   str(unique(hh[,2])) 
				str(unique(hh[,2]))
				mostFakeInst = hetInstFreq[which(hetInstFreq>= hhSubjCutOff)] 
				# str(names(mostFakeInst))
				# str(unique(names(mostFakeInst)))
				mostFakeInst4plink = fam[is.element(fam[,2], names(mostFakeInst)), 1:2]

				write.table(mostFakeInst4plink, file="mostFakeInst4plink.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
				## remove these fake SNPs
				system( paste0(plink, " --bfile ", inputPrefix, " --remove mostFakeInst4plink.txt --make-bed --out ", outputPrefix) )
				system( paste0("rm mostFakeInst4plink.txt") )

				## generate hetSNPsFreq in .txt file 
				InstWithHetSNPs = as.data.frame(hetInstFreq, stringsAsFactors=FALSE)
				# str(InstWithHetSNPs)
				## sort -nr -k 2 
				InstWithHetSNPs = InstWithHetSNPs[order(InstWithHetSNPs[,2], decreasing=T),] 
				write.table(InstWithHetSNPs, file=outputFile_subjHetFreqAll, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 

				## remaining males with heterozygous SNPs 
				InstWithHetSNPsub = InstWithHetSNPs[!is.element(InstWithHetSNPs[,1], names(mostFakeInst)), ]
				# str(InstWithHetSNPsub)
				write.table(InstWithHetSNPsub, file=outputFile_subjHetFreqRetained, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
			# head(InstWithHetSNPsub)

	  		} 

			system( paste0("rm male23nonPAR.*") )
 
			## output a text file that stores all heterozygous SNPs with their frequency after the above steps; 
			system( paste0(plink, " --bfile ", outputPrefix, " --chr 23 --filter-males --make-bed --out male23nonPAR" ))

			if ( is.na(file.size("male23nonPAR.hh")) ){  
				system(paste0('touch ', outputFile_SNPhhFreqAll) ) ## with empty output	  
		    } else { 	 
				hh = read.table("male23nonPAR.hh", stringsAsFactors=F)  
				hetSNPsFreq = table(hh[,3])  
				hetSNPsWithInstNumber = as.data.frame(hetSNPsFreq, stringsAsFactors=FALSE)
				hetSNPsWithInstNumber = hetSNPsWithInstNumber[order(hetSNPsWithInstNumber[,2], decreasing=T),] 
				## all heterozygous SNPs 
				write.table(hetSNPsWithInstNumber, file=outputFile_SNPhhFreqAll, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")   
				}
	 		system( paste0("rm male23nonPAR.*") )	
}	
 
