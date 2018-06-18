#######################################################################
#
# Package Name: Gimpute
# Description:
#   Gimpute -- An efficient genetic data imputation pipeline
#
# Gimpute R package, An efficient genetic data imputation pipeline
# Copyright (C) 2018  Junfang Chen (junfang.chen@zi-mannheim.de)
# All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 

##########################################################################
#
   


#' Remove heterozygous SNPs in male chromosome X
#'
#' @description
#' Remove heterozygous SNPs in haploid male chromosome X only if chromosome X 
#' data exists.


#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param hhCutOff the cutoff for removing male haploid heterozygous SNPs 
#' on the chromosome X.
#' @param outputPrefix the prefix of the output PLINK binary files.
#' @param outputHetSNPfile the output pure text file that stores 
#' all heterozygous SNPs with their frequency (the number of males for the 
#' corresponding SNP), if any. Lines are sorted by descending number.
#' @param outputRetainSNPfile the output pure text file that stores 
#' retained heterozygous SNPs with their frequency (the number of males for 
#' the corresponding SNP), if any. Lines are sorted by descending number.

#' @return 1.) The output PLINK binary files. 2.) A pure text files (if any)
#' with two columns: SNPs and their corresponding frequency. 3.) After SNP 
#' removal, a pure text files (if any) with two columns: SNPs and their 
#' corresponding frequency.
#' @details  A haploid heterozygous is a male genotype that is heterozygous, 
#' which could be an error given the haploid nature of the male X chromosome.
#' In principle, one has to remove all heterozygous SNPs of chromosome X in males. 
#' However, too many SNPs might be removed in some data sets. 
#' Therefore a small percentage of such SNPs in the data set is allowed.

#' @export  
#' @author Junfang Chen 
##' @examples   
   
  
removedSnpHetX <- function(plink, inputPrefix, hhCutOff, outputPrefix, 
						   outputHetSNPfile, outputRetainSNPfile){

	bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
	chrDist <- table(bim[,1])
	chr23check <- is.element(names(chrDist), 23)
	if (chr23check == TRUE) { 
		## just to get .hh file and .fam file 
		system(paste0(plink, " --bfile ", inputPrefix, 
			   " --chr 23 --filter-males --make-bed --out male23nonPAR"))

		if (is.na(file.size("male23nonPAR.hh"))) { 
			## copy/rename all snp info updated plink files
			renamePlinkBFile(inputPrefix, outputPrefix, action="copy")  
	    } else { 	 

			## *.hh: A text file with one line per error (sorted primarily by 
			## variant ID, secondarily by sample ID) with the following three fields:
			# Family ID  Within-family ID Variant ID
			hh <- read.table("male23nonPAR.hh", stringsAsFactors=FALSE)
			fam <- read.table("male23nonPAR.fam", stringsAsFactors=FALSE)

			hetSNPsFreq <- table(hh[,3])
			# hetSNPFreqFreq <- table(hetSNPs) 

			cutoff4removeHetSNP <- nrow(fam) * hhCutOff
			mostFakeSNPs <- hetSNPsFreq[which(hetSNPsFreq >= cutoff4removeHetSNP)] 
			mostFakeSNPs <- names(mostFakeSNPs)
			str(mostFakeSNPs)
			write.table(mostFakeSNPs, file="mostFakeSNPs.txt", quote=FALSE, 
					    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 
			## remove these fake SNPs
			system(paste0(plink, " --bfile ", inputPrefix, 
				   " --exclude mostFakeSNPs.txt --make-bed --out ", outputPrefix))
			## remove unwanted files
			system(paste0("rm mostFakeSNPs.txt"))

			## generate hetSNPsFreq in .txt file 
			hetSNPinstNum <- as.data.frame(hetSNPsFreq, stringsAsFactors=FALSE)
			hetSNPinstNum <- hetSNPinstNum[order(hetSNPinstNum[,2], decreasing=TRUE),] 
			## all heterozygous SNPs 
			write.table(hetSNPinstNum, file=outputHetSNPfile, quote=FALSE, 
						row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 

			## remaining heterozygous SNPs 
			hetSNPinstNumSub <- hetSNPinstNum[!is.element(hetSNPinstNum[,1], mostFakeSNPs), ]
			write.table(hetSNPinstNumSub, file=outputRetainSNPfile, 
					    quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 
	    } 

		system(paste0("rm male23nonPAR.* ", inputPrefix,".*")) 

	} else { 
		## copy/rename plink files
		renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
		system(paste0("rm ", inputPrefix,".*"))
	}
}	

 
  
 
##########################################   
##########################################
#' Remove male subjects with haploid heterozygous SNPs 
#'
#' @description
#' Determine the frequency of male subjects that have heterozygous SNPs on 
#' chromosome X and a reasonable cutoff to remove those affect males, if 
#' chromosome X data exists.


#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param hhSubjCutOff the cutoff for removing male subjects with haploid 
#' heterozygous SNPs on the chromosome X.
#' @param outputPrefix the prefix of the output PLINK binary files.
#' @param outputSubjHetFile the output pure text file that stores male subjects 
#' that have heterozygous SNPs with their frequency (if any), i.e. the number of 
#' .hh SNPs in this male. Lines are sorted by descending number.
#' @param outputRetainSubjectFile the output pure text file that stores
#' male subjects that have heterozygous SNPs with their frequency after 
#' subject removal (if any). Lines are sorted by descending number.
#' @param outputHetSNPfile the output pure text file that stores all 
#' heterozygous SNPs with their frequency (the number of males for this SNP)
#' , if any. Lines are sorted by descending number.

#' @return 1.) The output PLINK binary files. 2.) A pure text file with 
#' two columns: heterozygous male subjects and their corresponding heterozygous
#' SNPs. 3.) After subject removal, a pure text file consisting of two columns:
#' heterozygous male subjects and their corresponding heterozygous SNPs. 
#' A pure text file with two columns: all heterozygous SNPs and their frequency.

#' @details  A haploid heterozygous is a male genotype that is heterozygous, 
#' which could be an error given the haploid nature of the male X chromosome.
#' In principle, one has to remove all males that have heterozygous SNPs on the 
#' chromosome X. However, too many males might be removed in some data sets. 
#' Therefore a small percentage of such males in the data set is allowed.
 
#' @export 
#' @author Junfang Chen 
 

removedMaleHetX <- function(plink, inputPrefix, hhSubjCutOff, outputPrefix, 
							outputSubjHetFile, outputRetainSubjectFile, 
							outputHetSNPfile){		

	bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
	chrDist <- table(bim[,1])
	chr23check <- is.element(names(chrDist), 23)
	if (chr23check == TRUE) { 
		## just to get .hh file and .fam file 
		system(paste0(plink, " --bfile ", inputPrefix, 
			   " --chr 23 --filter-males --make-bed --out male23nonPAR"))
		if ( is.na(file.size("male23nonPAR.hh")) ){  
			## copy/rename all snp info updated plink files
			renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
			system(paste0("touch ", outputSubjHetFile) )
			system(paste0("touch ", outputRetainSubjectFile) )
  
	    } else { 	

			## .hh: A text file with one line per error (sorted primarily by  
			## variant ID,secondarily by sample ID) with the following three fields:
			## Family ID  Within-family ID Variant ID
			hh <- read.table("male23nonPAR.hh", stringsAsFactors=FALSE)
			fam <- read.table("male23nonPAR.fam", stringsAsFactors=FALSE)
							
			hetInstFreq <- table(hh[,2])  
			str(unique(hh[,2]))
			mostFakeInst <- hetInstFreq[which(hetInstFreq >= hhSubjCutOff)]  
			mostFakeInst4plink <- fam[is.element(fam[,2], names(mostFakeInst)), 1:2]

			write.table(mostFakeInst4plink, file="mostFakeInst4plink.txt", quote=FALSE, 
				        row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 
			## remove these fake SNPs
			system(paste0(plink, " --bfile ", inputPrefix, 
				   " --remove mostFakeInst4plink.txt --make-bed --out ", outputPrefix))
			system(paste0("rm mostFakeInst4plink.txt"))

			## generate hetSNPsFreq in .txt file 
			InstHetSNP <- as.data.frame(hetInstFreq, stringsAsFactors=FALSE) 
			InstHetSNP <- InstHetSNP[order(InstHetSNP[,2], decreasing=T),] 
			write.table(InstHetSNP, file=outputSubjHetFile, quote=FALSE, 
					    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 

			## remaining males with heterozygous SNPs 
			InstHetSNPsub <- InstHetSNP[!is.element(InstHetSNP[,1], names(mostFakeInst)), ] 
			write.table(InstHetSNPsub, file=outputRetainSubjectFile, 
					    quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 
  		} 

		system(paste0("rm male23nonPAR.*")) 
		## output a text file that stores all heterozygous SNPs 
		## with their frequency after the above steps; 
		system(paste0(plink, " --bfile ", outputPrefix, 
			   " --chr 23 --filter-males --make-bed --out male23nonPAR"))

		if (is.na(file.size("male23nonPAR.hh"))) {  
			system(paste0("touch ", outputHetSNPfile)) ## with empty output	  
	    } else { 	 
			hh <- read.table("male23nonPAR.hh", stringsAsFactors=FALSE)  
			hetSNPsFreq <- table(hh[,3])  
			hetSNPinstNum <- as.data.frame(hetSNPsFreq, stringsAsFactors=FALSE)
			hetSNPinstNum <- hetSNPinstNum[order(hetSNPinstNum[,2], decreasing=TRUE),] 
			## all heterozygous SNPs 
			write.table(hetSNPinstNum, file=outputHetSNPfile, 
					    quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")   
		}
 		system(paste0("rm male23nonPAR.*"))	
 	} else { 
		## copy/rename plink files
		renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
	}
}	
 




##########################################   
##########################################
#' Set haploid heterozygous SNPs as missing 
#'
#' @description
#' Set all heterozygous alleles of chromosome X SNPs in male as missing.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.
#' @return The output PLINK binary files after setting haploid heterozygous 
#' SNPs as missing.

#' @export 

#' @author Junfang Chen 
###' @examples  
  
setHeteroHaploMissing <- function(plink, inputPrefix, outputPrefix){
	 
	system(paste0(plink, " --bfile ", inputPrefix, 
		   " --set-hh-missing --make-bed --out ", outputPrefix))  
 
}  




##########################################   
##########################################
#' Remove SNPs with missing values
#'
#' @description
#' Remove SNPs with missingness of greater than a certain threshold.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param snpMissCutOff the cutoff of the missingness for removing SNPs.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after removing SNPs with pre-defined 
#' missing values.

#' @export 
#' @author Junfang Chen 
##' @examples  
 
removedSnpMiss <- function(plink, snpMissCutOff, inputPrefix, outputPrefix){
 
	system(paste0(plink, " --bfile ", inputPrefix, " --geno ", snpMissCutOff, 
		   " --make-bed --out ", outputPrefix))  
	 
}  
   




##########################################   
########################################## 
#' Remove subjects with missing values
#'
#' @description
#' Remove Subjects or instances with missingness of greater than a certain 
#' threshold.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param sampleMissCutOff the cutoff of the missingness for removing 
#' subjects/instances.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after removing subjects with  
#' a pre-defined removal cutoff.

#' @export 
#' @author Junfang Chen 
##' @examples  

 
removedInstMiss <- function(plink, sampleMissCutOff, inputPrefix, outputPrefix){
 
	system(paste0(plink, " --bfile ", inputPrefix, " --mind ", sampleMissCutOff, 
		   " --make-bed --out ", outputPrefix)) 
	system(paste0("rm ", outputPrefix, ".irem"))

}   



##########################################   
##########################################	
#' Remove subjects with abnormal autosomal heterozygosity deviation
#'
#' @description
#' Remove subjects with great autosomal heterozygosity deviation. 

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param Fhet the cutoff of the autosomal heterozygosity deviation.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after removing subjects with great 
#' autosomal heterozygosity deviation.
#' @details If the cutoff of the autosomal heterozygosity deviation is set to 
#' be greater than 0.2, i.e. |Fhet| >= 0.2, then this analysis will 
#' automatically skip haploid markers (male X and Y chromosome markers).

#' @export 
#' @author Junfang Chen 
##' @examples  
 

removedInstFhet <- function(plink, Fhet, inputPrefix, outputPrefix){ 

	system(paste0(plink, " --bfile ", inputPrefix, " --het --out ", outputPrefix))
	##  F inbreeding coefficient estimate
	autoHet <- read.table(file=paste0(outputPrefix, ".het"), header=T)  
	fhet <- autoHet[, "F"]
	qc_data_fhet <- autoHet[which(abs(fhet) >= Fhet), c(1, 2)]  
	## the individual IDs to be removed  
	write.table(qc_data_fhet, file=paste0(outputPrefix, ".txt"), quote=FALSE, 
				row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
  
	## To remove certain individuals 
	system(paste0(plink, " --bfile ", inputPrefix, " --remove ", 
		   paste0(outputPrefix, ".txt"), " --make-bed --out ", outputPrefix))
	system(paste0("rm ", outputPrefix, ".het"))
	system(paste0("rm ", outputPrefix, ".txt"))
 
}   



##########################################   
##########################################
#' Reset paternal and maternal codes  
#'
#' @description
#' Reset paternal and maternal codes of non-founders if parents not present. 
#' Replace the paternal ID and maternal ID of subjects (childs) by the
#' value zero if the paternal ID and the maternal ID do not belong to any
#' subject (parent) with the same family ID as the child. 

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files.
#' @details By default, if parental IDs are provided for a sample, 
#' they are not treated as a founder even if neither parent is 
#' in the dataset.  With no modifiers, --make-founders clears 
#' both parental IDs whenever at least one parent is not in the dataset, 
#' and the affected samples are now considered founders. 



#' @export  
#' @author Junfang Chen 
##' @examples  

# 
 
removedParentIdsMiss <- function(plink, inputPrefix, outputPrefix){ 

	# Remove the parent IDs which do not belong to subjects
	system(paste0(plink, " --bfile ", inputPrefix, 
		   " --make-founders require-2-missing --make-bed --out ", outputPrefix)) 
 
}   
 

##########################################   
##########################################  
#' Remove SNPs with difference in SNP missingness between cases and controls. 
#'
#' @description
#' Remove SNPs with difference in SNP missingness between cases and controls. 
#' To test for differential call rates between cases and controls for each SNP
#' This only works for the genotype data with cases and controls. 

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param snpMissDifCutOff the cutoff of the difference in missingness between 
#' cases and controls. 
#' @param outputPrefix the prefix of the output PLINK binary files.
#' @param caseControl a logical value indicating whether the data contains 
#' both case-control subjects.

#' @return The output PLINK binary files.

#' @export 
#' @author Junfang Chen 
##' @examples  


removedSnpMissDiff <- function(plink, inputPrefix, snpMissDifCutOff, 
							   outputPrefix, caseControl){

	if (caseControl == FALSE){  
			## this is only for the control data set
			renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
	} else if (caseControl == TRUE) { 
		outputPrefix.tmp <- paste0(outputPrefix, "tmp")
		system (paste0(plink, " --bfile ", inputPrefix, 
				" --test-missing --out ", outputPrefix.tmp) ) 
		## Write case/control missingness test to [ *.missing ]  
		## compute differential call rates 
		## between cases and controls for each SNP 
		ccmissing <- read.table(file=paste0(outputPrefix.tmp, ".missing"), 
								header=T, sep="")  
		fmiss <- abs(ccmissing[, "F_MISS_A"] - ccmissing[, "F_MISS_U"])
		whmiss <- which(fmiss >= snpMissDifCutOff) 
		SNPmissCC <- ccmissing[whmiss, "SNP"]
		SNPdifCallrate <- paste0(outputPrefix, ".txt")
		write.table(SNPmissCC, file=SNPdifCallrate, quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

		## exclude SNPs 
		system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
			   SNPdifCallrate, " --make-bed --out ", outputPrefix))  
		system(paste0("rm ", outputPrefix.tmp, ".*"))
		system(paste0("rm ", SNPdifCallrate))
	}
}   
 


##########################################   
##########################################    
#' remove chromosome X SNPs in females
#'
#' @description
#' Remove SNPs on the chromosome X with a pre-defined cutoff for 
#' missingness in females.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param femaleChrXmissCutoff the cutoff of the missingness 
#' in female chromosome X SNPs.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files.

#' @export 
#' @author Junfang Chen 
##' @examples  


removedSnpFemaleChrXmiss <- function(plink, femaleChrXmissCutoff, 
	 								 inputPrefix, outputPrefix){ 

	bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
	chrDist <- table(bim[,1])
	chr23check <- is.element(names(chrDist), 23)
	if (chr23check == TRUE) { 
	  	## additional QC (female-chrX SNPs, missingness ok?)  
	 	outputPrefix.tmp1 <- paste0(outputPrefix, "tmp1")
		outputPrefix.tmp2 <- paste0(outputPrefix, "tmp2")
		system(paste0(plink, " --bfile ", inputPrefix, 
			   " --filter-females --chr 23 --make-bed --out ", 
			   outputPrefix.tmp1))
		system(paste0(plink, " --bfile ", inputPrefix, 
			   " --filter-females --chr 23 --geno ", femaleChrXmissCutoff, 
			   " --make-bed --out ", outputPrefix.tmp2) )

		 ## check if equal  
		femaleChrXorig <- read.table(paste0(outputPrefix.tmp1, ".bim"), 
									 stringsAsFactors=FALSE) 
		femaleChrXMiss <- read.table(paste0(outputPrefix.tmp2, ".bim"), 
									 stringsAsFactors=FALSE)  
		snps2removed <- setdiff(femaleChrXorig[,2], femaleChrXMiss[,2]) 
		snps2removedfile <- paste0(outputPrefix, ".txt")
		write.table(snps2removed, file=snps2removedfile, quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 

		system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
			   snps2removedfile, " --make-bed --out ", outputPrefix))
		system(paste0("rm ", outputPrefix.tmp1, ".* ", outputPrefix.tmp2, ".*")) 
		system(paste0("rm ", snps2removedfile)) 
	} else { 
		## copy/rename plink files
		renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
	}
}





##########################################   
##########################################
#' Hardy weinberg equilibrium test for autosomal SNPs in controls.
#'
#' @description
#' Remove autosomal SNPs deviating from Hardy weinberg equilibrium in controls.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param pval the p-value cutoff for controlling HWE test in control subjects. 
#' Only autosomal SNPs are considered. 
#' @param outputPvalFile the output pure text file that stores autosomal SNPs and 
#' their sorted HWE p-values.
#' @param outputSNPfile the output pure text file that stores the removed SNPs, 
#' one per line.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files.

#' @export 
#' @author Junfang Chen 
##' @examples  


##  

removedSnpHWEautoControl <- function(plink, inputPrefix, pval, outputPvalFile, 
									 outputSNPfile, outputPrefix){ 
 
	outputPrefix.tmp <- paste0(outputPrefix, "tmp")
	system(paste0(plink, " --bfile ", inputPrefix, 
		   " --filter-controls --hardy --autosome --make-bed --out ", 
		   outputPrefix.tmp))  
	## read HWE p values 
	hweCheck <- read.table(file=paste0(outputPrefix.tmp, ".hwe"), 
						   header=TRUE, stringsAsFactors=FALSE) 
	## for controls 
	hweControl <- hweCheck[which(hweCheck$TEST == "UNAFF"), ]  
	snpHweValAutoCt <- subset(hweControl, select=c(SNP, P))
	snpHweValAutoCt <- snpHweValAutoCt[order(snpHweValAutoCt[,"P"]),] 
	write.table(snpHweValAutoCt, file=outputPvalFile, quote=FALSE, 
				row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	removedSNPhweControl <- hweControl[which(hweControl$P <= pval), "SNP"]
	write.table(removedSNPhweControl, file=outputSNPfile, quote=FALSE, 
				row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

	## exclude SNPs 
	system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
		   outputSNPfile, " --make-bed --out ", outputPrefix)) 
	system(paste0("rm ", outputPrefix.tmp, ".*"))

}



##########################################   
##########################################  
#' Hardy weinberg equilibrium test for chromosome X SNPs in female controls. 
#'
#' @description
#' Hardy weinberg equilibrium test for SNPs on the chromosome X in 
#' female controls.  

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param pval the p-value cutoff for controlling HWE test in female control 
#' subjects. Only chromosome X SNPs are considered. 
#' @param outputPvalFile the output pure text file that stores chromosome X SNPs 
#' and their sorted HWE p-values.
#' @param outputSNPfile the output pure text file that stores the removed SNPs, 
#' one per line.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files.

#' @export 
#' @author Junfang Chen 
##' @examples  

##   
 
removedSnpFemaleChrXhweControl <- function(plink, inputPrefix, pval, 
										   outputPvalFile, outputSNPfile, 
										   outputPrefix){ 
	
	bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
	chrDist <- table(bim[,1])
	chr23check <- is.element(names(chrDist), 23)
	if (chr23check == TRUE) { 
	 
		outputPrefix.tmp <- paste0(outputPrefix, "tmp") 
		system(paste0(plink, " --bfile ", inputPrefix, 
			   " --filter-females --filter-controls --chr 23 --hardy ", 
			   " --make-bed --out ", outputPrefix.tmp) )
		## read p values
		hweCheck <- read.table(file=paste0(outputPrefix.tmp, ".hwe"), 
							   header=TRUE, stringsAsFactors=FALSE) 
		## for controls 
		hweControl <- hweCheck[which(hweCheck$TEST == "UNAFF"), ] # 
		snpHweValChrXCt <- subset(hweControl, select=c(SNP, P))
		snpHweValChrXCt <- snpHweValChrXCt[order(snpHweValChrXCt[,"P"]),]
		write.table(snpHweValChrXCt, file=outputPvalFile, quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		## remove bad SNPs
		removedSNPhweControl <- hweControl[which(hweControl$P <= pval), "SNP"]
		write.table(removedSNPhweControl, file=outputSNPfile, quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

		system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
				outputSNPfile, " --make-bed --out ", outputPrefix))
		system(paste0("rm ", outputPrefix.tmp, ".*")) 

	} else {
		## copy/rename plink files
		renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
	}
}
 

 




##########################################   
##########################################
#' Population outlier detection
#'
#' @description
#' Principle component analysis (PCA) on the genotype data is performed 
#' to detect population outliers, and the first two PCs 
#' are plotted for the visualization. 

#' @param gcta an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param nThread the number of threads used for computation. 
#' @param outputPC4subjFile the pure text file that stores all the subject IDs 
#' and their corresponding eigenvalues of the first two principle components.
#' @param outputPCplotFile the plot file for visualizing the first two 
#' principle components of all investigated subjects.

#' @return The output pure text file and plot file for storing first two 
#' principle components of study subjects.
#' @details Before population outlier detection, it's better to perform QC 
#' on the genotype data. 
#' Only autosomal genotypes are used for principle component analysis. 

#' @export 
#' @import lattice  

#' @author Junfang Chen 
##' @examples  
 

plotPCA4plink <- function(gcta, inputPrefix, nThread, 
						  outputPC4subjFile, outputPCplotFile){ 

	autosomefn <- paste0(inputPrefix, "Autosome")
	system(paste0(gcta, " --bfile ", inputPrefix, 
		   " --make-grm --autosome --out ", autosomefn, " --thread-num ", nThread))
	system(paste0(gcta, " --grm ", autosomefn, " --pca 20 --out ", 
		   autosomefn, " --thread-num ", nThread))

	eigen <- read.table(file=paste0(autosomefn,".eigenvec"), stringsAsFactors=FALSE)
	pcs <- eigen[,1:4] ## first two PCs in the 3rd and 4th column.
	write.table(pcs, outputPC4subjFile, quote=FALSE, row.names=FALSE, 
			    col.names=FALSE, eol="\r\n", sep=" ")
	pcWithGroup <- cbind(pcs, stringsAsFactors=FALSE)

	png(outputPCplotFile, width=8, height=6, units="in", res=800)
	print( xyplot(pcWithGroup[,4] ~ pcWithGroup[,3], data=pcWithGroup, 
	       auto.key=list(space="right"),  
	       jitter.x=TRUE, jitter.y=TRUE, xlab="PC1", ylab="PC2") )
	dev.off()
	## remove unwanted files
	system(paste0("rm ", autosomefn, ".*"))
}


 

 
##
######################################################
######################################################
#' Remove population outliers
#'
#' @description
#' Remove population outliers by using principle component analysis.
  

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param gcta an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param cutoff the cutoff that distinguishes the eigenvalues of the outliers  
#' from ordinary population. If it is null, then there are no outliers or 
#' outliers are not required to be removed.  
#' @param cutoffSign the cutoff sign: 'greater' or 'smaller' that determines 
#' if the outliers should be greater or smaller than the cutoff value.
#' @param inputPC4subjFile the pure text file that stores all the subject IDs 
#' and their corresponding eigenvalues of the first two principle components.
#' @param outputPC4outlierFile the pure text file that stores the outlier IDs 
#' and their corresponding eigenvalues of the first two principle components.
#' @param outputPCplotFile the plot file for visualizing the first two 
#' principle components of all subjects without population outliers.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return 1.) The output PLINK binary files after outlier removal. 
#' 2.) The output pure text file (if any) for storing removed outlier IDs 
#' and their corresponding PCs. 3.) The plot file (if any) for visualizing 
#' the first two principle components after outlier removal.

#' @details This function is used for removing population outliers. 
#' If the outliers are necessary to be removed, then one uses the eigenvalues 
#' from the first principle component as a criterion to find out the outliers 
#' by assigning an appropriate cutoff. 
  

#' @export 
#' @import lattice  
#' @author Junfang Chen 
  
 
removeOutlierByPCs <- function(plink, gcta, inputPrefix, cutoff, cutoffSign, 
							   inputPC4subjFile, outputPC4outlierFile, 
							   outputPCplotFile, outputPrefix){

	## if no outliers or no need to remove PC outliers. 
	if (is.null(cutoff) == TRUE) { 
		## copy/rename all snp info updated plink files
		renamePlinkBFile(inputPrefix, outputPrefix, action="copy")
	} else {  
		subjID_PCs <- read.table(file=inputPC4subjFile, stringsAsFactors=FALSE)   
		if (length(cutoff) > 1) { 
			## if the outliers should be removed on both side of the cluster
			## detected by PC1
			outliersPC1v1 <- subjID_PCs[which(subjID_PCs[,3] <= cutoff[1]), ] 
			outliersPC1v2 <- subjID_PCs[which(subjID_PCs[,3] >= cutoff[2]), ]  
			subjID_PCs4outlier <- rbind(outliersPC1v1, outliersPC1v2)
		} else {
	  	    if (cutoffSign == "smaller"){ 
	  	    	## detected by PC1
	  	  	    subjID_PCs4outlier <- subjID_PCs[which(subjID_PCs[,3] <= cutoff), ] 
	  	    } else if (cutoffSign == "greater"){
	  	    	## detected by PC1
	  	  	    subjID_PCs4outlier <- subjID_PCs[which(subjID_PCs[,3] >= cutoff), ]
	  	    }
		}	  
		## sorted by first PC.
	 	subjID_PCs4outlierSorted <- subjID_PCs4outlier[order(subjID_PCs4outlier[,3]), ] 
	 	write.table(subjID_PCs4outlierSorted, file=outputPC4outlierFile, quote=FALSE, 
	 				row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	 	subjID4outlierTmp  <- subjID_PCs4outlierSorted[,1:2]
	 	subjID4outlierTmpFile <- "subjID4outlierTmp.txt"
		write.table(subjID4outlierTmp, file=subjID4outlierTmpFile, quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		system(paste0(plink, " --bfile ", inputPrefix, " --remove ", 
			   subjID4outlierTmpFile, " --make-bed --out ", outputPrefix))
		system(paste0("rm ", subjID4outlierTmpFile)) 
		## Plot first two PCs again
		outputPC4subjFiletmp <- "outputPC4subjFile.txt" ## PCs for the retained subjects 
		plotPCA4plink(gcta, inputPrefix=outputPrefix, outputPC4subjFiletmp, outputPCplotFile)
		system(paste0("rm ", outputPC4subjFiletmp))
	}	
}
 




