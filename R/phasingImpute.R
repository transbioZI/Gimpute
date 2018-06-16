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
## removedMonoSnp.R
##########################################################################
# 

#' Exclude monomorphic SNPs 
#'
#' @description
#' Detect monomorphic SNPs from PLINK BIM file and exclude them if any.

#' @param plink an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK files.
#' @param outputPrefix the prefix of the output PLINK files 
#' after removing monomorphic SNPs.
#' @param outputSNPfile the output pure text file that stores 
#' the removed monomorphic SNPs, one per line.

#' @return  The output PLINK files after removing monomorphic SNPs 
#' and a pure text file with removed monomorphic SNPs.
#' @export 

#' @author Junfang Chen 


removedMonoSnp <- function(plink, inputPrefix, outputPrefix, outputSNPfile){  

	## input BIM file  
	bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)
	monoSNPs <- bim[which(bim[,5] == 0),2]  
	write.table(monoSNPs, file=outputSNPfile, quote=FALSE, row.names=FALSE, 
		 		col.names=FALSE, eol="\r\n", sep=" ") 
  
	system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", outputSNPfile, 
		   " --make-bed --out ", outputPrefix)) 
}


 

##########################################################################
## chrWiseSplit.R
########################################################################## 

#' Split genome-wide genotyping data into chromosome-wide PLINK files.
#'
#' @description
#' Split the whole genome genotyping data chromosome-wise; 
#' allow parallel computating for all chromosomes.

#' @param plink an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK files before splitting.
#' @param outputPrefix the prefix of the output PLINK files after splitting  
#' separately for each chromosome, appended with the chromosome codes.
#' @param chrX_PAR1suffix  if chromosome 25 is available and with PAR1, 
#' then generate the suffix with X_PAR1 for chrX_PAR1. 
#' @param chrX_PAR2suffix  if chromosome 25 is available and with PAR2, 
#' then generate the suffix with X_PAR2 for chrX_PAR2.
#' @param nCore the number of cores used for parallel computation. 
#' The default value is 25.  

#' @return  The output PLINK files for each chromosome and possibly 
#' also the logical value for the pseudo-autosomal region (PAR)
#' indicating if PAR exists in the input genotyping data or not.   

#' @details If chromosome 25 is also available, namely the pseudo-autosomal
#' region of chromosome X, then further split chr25 (PAR or Chr_XY)
#' into PAR1 and PAR2 according to the genomic coordination GRCh37
#' from https://en.wikipedia.org/wiki/Pseudoautosomal_region.
#' The locations of the PARs within GRCh37 are:  
#' PAR1	X	60001	2699520; 
#' PAR2	X	154931044	155260560.   

#' @export 
#' @import doParallel  

#' @author Junfang Chen


chrWiseSplit <- function(plink, inputPrefix, chrX_PAR1suffix, 
						 chrX_PAR2suffix, nCore=25){ 

    ## check which chromosomes are available to be splitted from the .bim file
	bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)
	chrs <- as.integer(names(table(bim[,1])))
 
	chrslist <- as.list(chrs)
	mclapply(chrslist, function(i){
		cmd <- paste0(plink, " --bfile ", inputPrefix, " --chr ", i, 
			          " --make-bed --out ", inputPrefix, i)  
		system(cmd)
	}, mc.cores=nCore)

	## if chromosome 25 is also available then re-arrange it
 	if (is.element(25, chrs)){  

 		print("PAR is available in chrX!") 
		bim25 <- read.table(paste0(inputPrefix, "25.bim"), stringsAsFactors=FALSE) 
		pos4PAR1 <- c(60001, 2699520) 
		## first check for PAR1 and afterwards for PAR2
		if ( length(which(bim25[,4] <= pos4PAR1[2]))!= 0 ){ 
	   
	   		print("PAR1 is available in chrX!")
	   		bimPos4par1 <- which(bim25[,4]<= pos4PAR1[2])
			rs4PAR1 <- bim25[bimPos4par1,2]
			write.table(rs4PAR1, file="rs4PAR1.txt", quote=FALSE, 
						row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
			system(paste0(plink, " --bfile ", inputPrefix, 
				   "25 --extract rs4PAR1.txt --make-bed --out ", 
				   inputPrefix, chrX_PAR1suffix))
			par1 <- TRUE  
			## check for PAR2, if any SNPs out of PAR1, PAR2 also available 
			if (length(bimPos4par1)<nrow(bim25)){
				print("PAR2 is available in chrX!") 
				## just to exclude rs4PAR1.txt
				system(paste0(plink, " --bfile ", inputPrefix, 
					   "25 --exclude rs4PAR1.txt --make-bed --out ", 
					   inputPrefix, chrX_PAR2suffix))
				par2 <- TRUE 
			} else { 
				print("PAR2 is NOT available in chrX!") 
			}

		} else { 
			print("PAR2 is available in chrX! But NOT PAR1, all chr25 on PAR2")
			par1 <- FALSE
			par2 <- TRUE 
			renamePlinkBFile(inputPrefix=paste0(inputPrefix, "25"), 
							 outputPrefix=paste0(inputPrefix, chrX_PAR2suffix), 
							 action="copy") 
		} 
	} else {  
		print("PAR is NOT available in chrX!") 
		par1 <- FALSE
		par2 <- FALSE
	} 

	return(par=list(par1, par2)) 
}
 






##########################################################################
## chunk4eachChr.R
########################################################################## 
#' Chunk each chromosome into multiple segments
#'
#' @description
#' Chunk each chromosome genotyping data into multiple segments 
#' by a predefined window size.

#' @param inputPrefix the prefix of the input PLINK files for each chromosome.
#' @param outputPrefix the prefix of the output pure text files that keep 
#' all the chunks for each chromosome separately.
#' @param chrs specifiy the chromosome codes for chunking.
#' @param windowSize  the window size of each chunk.

#' @return  The output pure text files include all the chunks 
#' for each chromosome separately. 
#' @export 

#' @author Junfang Chen 

 

chunk4eachChr <- function(inputPrefix, outputPrefix, chrs, windowSize){  

	for (i in chrs){ 

		bimfilename <- paste0(inputPrefix, i, ".bim")
		bimdata <- read.table(file=bimfilename, sep="\t", stringsAsFactors=FALSE)
		position <- bimdata[,4]
		posStart <- head(position,1)
		posEnd <- tail(position,1)
		chunkStart <- seq(posStart, posEnd, windowSize)
		chunkEnd <- chunkStart + windowSize -1
		chunkMatrix <- cbind(chunkStart, chunkEnd)

		## positions are only within a chunk
		if (nrow(chunkMatrix) == 1){
			chunks <- chunkMatrix
		} else {  
			##  it may happen that only a few SNPs from the last chunk; 
			## but if the last chunk is large, then specify -allow_large_regions 
			chunks <- head(chunkMatrix, -1) ## merge last-second to last  
			chunks[nrow(chunks), 2] <- posEnd
		}
		
		## check if any chunk with NO snps within it   
		SNPcountsPerChunk <- c() 
		for (j in seq_len(nrow(chunks))){
			chunkbottom <- chunks[j,1]
			chunkup <- chunks[j,2]
			## ## which fall within chunk 
			tmp <- which(position >= chunkbottom & position <= chunkup)  
			tmp <- length(tmp)
			SNPcountsPerChunk <- c(SNPcountsPerChunk, tmp) 
		} 

	 	wh0 <- which(SNPcountsPerChunk == 0)
	 	print(paste0("chr",i))
	 	print(wh0)
	 	chunkLength <- nrow(chunks) - length(wh0)
	 	## remove such chunks if wh0  

	 	if (length(wh0) != 0){ chunks <- chunks[-wh0,] }
	 	print(nrow(chunks) == chunkLength)
		chunkfilename <- paste0(outputPrefix, i, ".txt")
		write.table(chunks, file=chunkfilename, quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	}
}



	

##########################################################################
## prePhasingByShapeit.R
########################################################################## 
#' Prephasing genotypes using SHAPEIT
#'
#' @description
#' Perform prephasing for study genotypes by SHAPEIT for the autosomal and 
#' sex chromosome haplotypes using a reference panel (pre-set).

#' @param shapeit an executable program in either the 
#' current working directory or somewhere in the command path.
#' @param chrs specifiy the chromosome codes for phasing.
#' @param dataDIR the directory where genotype PLINK files are located.
#' @param prefix4plinkEachChr the prefix of PLINK files for each chromosome.
#' @param impRefDIR the directory where the imputation reference files 
#' are located.
#' @param phaseDIR the directory where resulting pre-phased files 
#' will be located.
#' @param nThread the number of threads used for computation.
#' @param effectiveSize this parameter controls the effective population size.
#' @param nCore the number of cores used for computation. This can be tuned 
#' along with nThread.
 

#' @return  The pre-phased haplotypes for given chromosomes.  
#' @details If ChrX is available then it is done differently by passing the flag 
#' --chrX to SHAPEIT.

#' @export 
#' @import doParallel  

#' @author Junfang Chen 
  

prePhasingByShapeit <- function(shapeit, chrs, dataDIR, 
								prefix4plinkEachChr, impRefDIR, phaseDIR, 
								nThread, effectiveSize, nCore){

	chrslist <- as.list(chrs)
	mclapply(chrslist, function(i){

		# GWAS data files in PLINK binary format
		GWASDATA_BED <- paste0(dataDIR, prefix4plinkEachChr, i, ".bed ") 
		GWASDATA_BIM <- paste0(dataDIR, prefix4plinkEachChr, i, ".bim ")
		GWASDATA_FAM <- paste0(dataDIR, prefix4plinkEachChr, i, ".fam ")
		# reference data files
		GENMAP_FILE <- paste0(impRefDIR, "genetic_map_chr", i, 
							  "_combined_b37.txt ")
		HAPS_FILE <- paste0(impRefDIR, "ALL_1000G_phase1integrated_v3_chr", i, 
							"_impute_macGT1.hap.gz ") 
		LEGEND_FILE <- paste0(impRefDIR, "ALL_1000G_phase1integrated_v3_chr", i, 
							  "_impute_macGT1.legend.gz ")
		SAMPLE_FILE <- paste0(impRefDIR, "ALL_1000G_phase1integrated_v3.sample ")

		# main output file
		OUTPUT_HAPS <- paste0(phaseDIR, "chr", i, ".haps ")     
		OUTPUT_SAMPLE <- paste0(phaseDIR, "chr", i, ".sample ")     
		OUTPUT_LOG <- paste0(phaseDIR, "chr", i, ".log ")    

		if (i != 23){  ## prePhasing for the autosome
			system(paste0(shapeit, 
			" --input-bed ", GWASDATA_BED, GWASDATA_BIM, GWASDATA_FAM, " \ ", 
			" --input-map ", GENMAP_FILE, " \ ",  
			"--input-ref ", HAPS_FILE, LEGEND_FILE, SAMPLE_FILE, " \ ", 
			"--thread ", nThread, " \ ", 
			"--effective-size ", effectiveSize, " \ ", 
			"--output-max ", OUTPUT_HAPS, OUTPUT_SAMPLE, " \ ", 
			"--output-log ", OUTPUT_LOG) )
		} else if (i == 23){
			system(paste0(shapeit, 
			" --input-bed ", GWASDATA_BED, GWASDATA_BIM, GWASDATA_FAM, " \ ", 
			" --input-map ", GENMAP_FILE, " \ ", 
			" --chrX \ ", 
			"--input-ref ", HAPS_FILE, LEGEND_FILE, SAMPLE_FILE, " \ ", 
			"--thread ", nThread, " \ ", 
			"--effective-size ", effectiveSize, " \ ", 
			"--output-max ", OUTPUT_HAPS, OUTPUT_SAMPLE, " \ ", 
			"--output-log ", OUTPUT_LOG) ) 
		}

	}, mc.cores=nCore)	 
} 

 



##########################################################################
## imputedByImpute2.R
########################################################################## 
#' Impute genotypes using IMPUTE2
#'
#' @description
#' Perform imputation by IMPUTE2 for the autosomal and sex chromosome 
#' prephased known haplotypes with a reference panel.

#' @param impute2 an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param chrs specifiy the chromosome codes for imputation.
#' @param prefixChunk  the prefix of the chunk files for each chromosome, 
#' along with the proper location directory.
#' @param phaseDIR the directory where prephased haplotypes are located.
#' @param impRefDIR the directory where the imputation reference files 
#' are located.
#' @param imputedDIR the directory where imputed files will be located.
#' @param prefix4plinkEachChr the prefix of IMPUTE2 files for each chunk.
#' @param nCore the number of cores used for computation.
#' @param effectiveSize this parameter controls the effective population size.
#' @param XPAR a logical value indicating whether --chrX flag should be 
#' passed for prephasing using SHAPEIT.
#' --chrX flag, specifically for chrX imputation'
#' @return  The imputed files for all chunks from given chromosomes.  
#' @export 
#' @import doParallel  

#' @author Junfang Chen 
#' @examples 
 

imputedByImpute2 <- function(impute2, chrs, prefixChunk, phaseDIR, 
						     impRefDIR, imputedDIR, prefix4plinkEachChr, 
						     nCore, effectiveSize){ 
 
	for (i in chrs){ 	

		chunkfn <- paste0(prefixChunk, i, ".txt")
		chunks <- read.table(chunkfn, sep=" ")
 
		chunklist <- as.list(seq_len(nrow(chunks)))
		mclapply(chunklist, function(j){

			chunkSTART <- chunks[j,1]
			chunkEND   <- chunks[j,2] 
			## Input: haplotypes from SHAPEIT phasing (method B)
			GWAS_HAPS_FILE <- paste0(phaseDIR, "chr", i, ".haps ") 
			GWAS_SAMP_FILE <- paste0(phaseDIR, "chr", i, ".sample ") 
			## reference data files
			## For other reference panels you want to modify the following setting  
			GENMAP_FILE <- paste0(impRefDIR, "genetic_map_chr", i, 
								  "_combined_b37.txt ")
			HAPS_FILE <- paste0(impRefDIR, "ALL_1000G_phase1integrated_v3_chr", i, 
							    "_impute_macGT1.hap.gz ") 
			LEGEND_FILE <- paste0(impRefDIR, "ALL_1000G_phase1integrated_v3_chr", i, 
							      "_impute_macGT1.legend.gz ")
 
			## main output file    
			OUTPUT_FILE <- paste0(imputedDIR, prefix4plinkEachChr, i, 
								  ".pos", chunkSTART, "-", chunkEND, ".impute2 ")   
			################## impute genotypes from GWAS haplotypes 
			autosomeCode = seq_len(22)
			if (is.element(i, autosomeCode)) { 
				## impute for the autosomes
				system(paste0(impute2, 
				" -iter 30  \ ", 
				" -burnin 10  \ ", 
				" -k_hap 500  \ ", 
				" -use_prephased_g  \ ",  
				" -m ", GENMAP_FILE, " \ ",  
				" -h ", HAPS_FILE, " \ ", 
				" -l ", LEGEND_FILE, " \ ", 
				" -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
				" -Ne ", effectiveSize, " \ ", 
				" -int ", chunkSTART, " ", chunkEND, " \ ", 
				" -buffer 1000  \ ",
				" -o ", OUTPUT_FILE, " \ ", 
				" -allow_large_regions \ ",
				" -seed 367946 \ " ))
			} else if (is.element(i, c("X_PAR1", "X_PAR2"))){  
				## impute for chrX PAR >> with an additional flag: --Xpar.
				system(paste0(impute2, 
				" -iter 30  \ ", 
				" -burnin 10  \ ", 
				" -k_hap 500  \ ", 
				" -use_prephased_g  \ ", 
				" -Xpar \ ",     #################
				" -m ", GENMAP_FILE, " \ ",  
				" -h ", HAPS_FILE, " \ ", 
				" -l ", LEGEND_FILE, " \ ", 
				" -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
				" -Ne ", effectiveSize, " \ ", 
				" -int ", chunkSTART, " ", chunkEND, " \ ", 
				" -buffer 1000  \ ",
				" -o ", OUTPUT_FILE, " \ ", 
				" -allow_large_regions \ ",
				" -seed 367946 \ " ))
			} else if (i == 23){ 
				## impute for chrX 
				## >> with an additional flag: --chrX, and sample_known_haps_g
				system(paste0(impute2, 
				" -iter 30  \ ", 
				" -burnin 10  \ ", 
				" -k_hap 500  \ ", 
				" -use_prephased_g  \ ", 
				" -chrX \ ",   #################     
				" -m ", GENMAP_FILE, " \ ",  
				" -h ", HAPS_FILE, " \ ", 
				" -l ", LEGEND_FILE, " \ ", 
				" -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
				" -sample_known_haps_g ", GWAS_SAMP_FILE, " \ ",   
				" -Ne ", effectiveSize, " \ ", 
				" -int ", chunkSTART, " ", chunkEND, " \ ", 
				" -buffer 1000  \ ",
				" -o ", OUTPUT_FILE, " \ ", 
				" -allow_large_regions \ ",
				" -seed 367946 \ " ))
			}

		}, mc.cores=nCore)  
 	} 
} 		
 
 




  
##########################################################################
## convertImpute2ByGtool.R 
##########################################################################  
#' Convert IMPUTE2 format files into PLINK format
#'
#' @description
#' Convert all chunks of IMPUTE2 format files into PLINK format using GTOOL.

#' @param gtool an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param chrs specifiy the chromosome codes for conversion. 
#' @param prefixChunk  the prefix of the chunk files for each chromosome, 
#' along with the location directory.
#' @param phaseDIR the directory where pre-phased files are located.
#' @param imputedDIR the directory where the imputated files are located.
#' @param prefix4plinkEachChr the prefix of the input IMPUTE2 files and 
#' also the output PLINK files for each chunk.
#' @param suffix4imputed the suffix of the IMPUTE2 format file that stores 
#' the imputed value.
#' @param postImputeDIR the directory where converted PLINK files will be located. 
#' @param nCore the number of cores used for computation.  
 
#' @return  The converted PLINK format files for each chunk from IMPUTE2 results.
#' @export 
#' @import doParallel  

#' @author Junfang Chen 
  

convertImpute2ByGtool <- function(gtool, chrs, prefixChunk, 
							      phaseDIR, imputedDIR, prefix4plinkEachChr, 
							      suffix4imputed, postImputeDIR, nCore){

	for (i in chrs){ 
		chunkfn <- paste0(prefixChunk, i, ".txt")
		chunks <- read.table(chunkfn, sep=" ") 
		chunklist <- as.list(seq_len(nrow(chunks)))

		mclapply(chunklist, function(j){ 
			chunkSTART <- chunks[j,1]
			chunkEND   <- chunks[j,2] 
			## INPUT data files
			SAM_FILE <- paste0(phaseDIR, "chr", i, ".sample")  
			GEN_FILE <- paste0(imputedDIR, prefix4plinkEachChr, i, 
							   ".pos", chunkSTART, "-", chunkEND, suffix4imputed) 
			## output PLINK files
			PED_FILE <- paste0(postImputeDIR, prefix4plinkEachChr, i, ".pos", 
							   chunkSTART, "-", chunkEND, ".ped") 
			MAP_FILE <- paste0(postImputeDIR, prefix4plinkEachChr, i, ".pos", 
				  			   chunkSTART, "-", chunkEND, ".map") 
			## converting chunk-wise
			system( paste0(gtool, " -G ",
							" --g ", GEN_FILE, " \ ", 
							" --s ", SAM_FILE, " \ ",  
							"--phenotype plink_pheno \ ", 
							"--chr ", i, " \ ", 
							"--ped ", PED_FILE, " \ ", 
							"--map ", MAP_FILE, " \ ", 
							"--snp ") )  
		}, mc.cores=nCore)
	} 
} 
  
	 
   
  


##########################################################################
## mergePlinkData.R
########################################################################## 
#' Merge chunk-wise PLINK files 
#'
#' @description
#' Merge all chunk-wise PLINK files into chromosome-wise PLINK files then 
#' assemble into a genome-wide PLINK file set. 

#' @param plink an executable program in either the current working 
#' directory or somewhere in the command path.
#' @param chrs specifiy the chromosome codes to be merged. 
#' @param prefix4plinkEachChr the prefix of the input chunk-wise PLINK files. 
#' @param prefix4mergedPlink  the prefix of the final output PLINK 
#' format files. 
#' @param nCore the number of cores used for computation.  

#' @return  The merged genome-wide PLINK format files.
#' @details Create a file containing a list chunk-wise PLINK PED and MAP 
#' file names. The prefix of these files must already indicate in which 
#' chromosome they belong to and files from the same chromosome will be 
#' combined. Then all chromosomal PLINK files are assembled together 
#' into one whole genome PLINK file set.
#' @export 
#' @import doParallel  

#' @author Junfang Chen 


 
 
mergePlinkData <- function(plink, chrs, prefix4plinkEachChr, 
						   prefix4mergedPlink, nCore){ 

	## firstly, only consider chromosomes from 1:23; as Xpar chrs 
	## are slightly different for processing.
	pureAutoChrs <- setdiff(chrs, c("X_PAR1", "X_PAR2")) 
	chrslist <- as.list(pureAutoChrs)   
	mclapply(chrslist, function(i){

		pedFile_chr <- system(paste0("ls ", prefix4plinkEachChr, i, ".*.ped"), 
							  intern=TRUE)
		mapFile_chr <- system(paste0("ls ", prefix4plinkEachChr, i, ".*.map"), 
							  intern=TRUE)	
		pedmap_chr <- paste0(pedFile_chr, " ", mapFile_chr)
		fA <- gsub(".ped", "", pedFile_chr[1])
		pedmap_tobeMerged <- pedmap_chr[-1]
		filesetname <- paste0("fileset_chr", i, ".txt")
		write.table(pedmap_tobeMerged, file=filesetname, quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	    system(paste0(plink, " --file ", fA, " --merge-list ", filesetname, 
	    	   " --make-bed --out gwasImputed_chr", i)) 

	}, mc.cores=nCore)

	## combine chrX_PAR and convert into chr25 
	if (is.element(c("X_PAR1"), chrs) | is.element(c("X_PAR2"), chrs) ){  

		pedFile_chr <- system(paste0("ls ", prefix4plinkEachChr, "X_PAR*.ped"), 
						      intern=TRUE)
		mapFile_chr <- system(paste0("ls ", prefix4plinkEachChr, "X_PAR*.map"), 
							  intern=TRUE)	
		pedmap_chr <- paste0(pedFile_chr, " ", mapFile_chr)
		fA <- gsub(".ped", "", pedFile_chr[1])
		pedmap_tobeMerged <- pedmap_chr[-1]
		filesetname <- paste0("fileset_chr25.txt")
		write.table(pedmap_tobeMerged, file=filesetname, quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		 
		arg <- paste0(plink, " --file ", fA, " --merge-list ", filesetname, 
					  " --allow-extra-chr --make-bed --out gwasImputed_oldchr25")
		system(arg) 
		## update chr code for XPAR --> 25
		bim <- read.table("gwasImputed_oldchr25.bim", stringsAsFactors=FALSE)
		updateSNPchr <- cbind(bim[,2], rep(25, length=nrow(bim))) 
		write.table(updateSNPchr, file="gwasImputed_newchr25.txt", quote=FALSE, 
					row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 


		system(paste0(plink, " --bfile gwasImputed_oldchr25 --allow-extra-chr ", 
			   " --update-chr gwasImputed_newchr25.txt 2 1",
			   " --make-bed --out gwasImputed_chr25"))  
		system("rm gwasImputed_oldchr25.* gwasImputed_newchr25.txt")
		# system( paste0("rm ", filesetname))
	}	 
	## combine all bed files
	bedFile_chr <- system(paste0("ls gwasImputed_chr*.bed"), intern=TRUE)
	bimFile_chr <- system(paste0("ls gwasImputed_chr*.bim"), intern=TRUE)	
	famFile_chr <- system(paste0("ls gwasImputed_chr*.fam"), intern=TRUE)	
	bfile_chr <- paste0(bedFile_chr, " ", bimFile_chr, " ", famFile_chr)
	fA <- paste0(gsub(".bed", "", bedFile_chr[1]))
	tobeMerged <- bfile_chr[-1]
	mergefilesetname <- paste0("mergeGwasImputed.txt")
	write.table(tobeMerged, file=mergefilesetname, quote=FALSE, 
				row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
 
	system(paste0(plink, " --bfile ", fA, " --merge-list ", mergefilesetname, 
		   " --make-bed --out ", prefix4mergedPlink))
 
} 

 

  
 

##########################################################################
## filterImputeData.R
########################################################################## 
#' Filter genetic variants    
#'
#' @description
#' Filter out genetic variants accoring to the imputation quality score.
#' 
#' @param plink an executable program in either the current working 
#' directory or somewhere in the command path.
#' @param suffix4impute2info the suffix of input IMPUTE2 generated files that 
#' store the imputation quality score for each variant from .impute2_info files.
#' @param outputInfoFile the output file of impute2 info scores consisting of 
#' two columns: all imputed SNPs and their info scores.  
#' @param infoScore the cutoff of filtering imputation quality score for 
#' each variant.   
#' @param badImputeSNPfile the output file of SNPs with bad info scores.  
#' @param inputPrefix the prefix of the input imputed PLINK format files. 
#' @param outputPrefix the prefix of the output filtered PLINK format files. 

#' @return A pure text file contains the info scores of all imputed SNPs with 
#' two columns: SNP names and the corresponding info scores. 
#' A pure text file with all excluded SNPs having bad info scores. 
#' The filtered PLINK format imputed files, 
#' @details Filter genetic variants accoring to the imputation quality score with 
#' the help of .impute2_info files generated by IMPUTE2. 
#' Often, we keep variants with imputation info score of greater than 0.6.    
#' Note that imputed SNPs with more than two alleles are not considered. 

#' @export 
#' @author Junfang Chen 
   

filterImputeData <- function(plink, suffix4impute2info, outputInfoFile, 
							 infoScore, badImputeSNPfile, inputPrefix, 
							 outputPrefix){ 

	## read each .impute2_info file, remove 1st line, add to another file and repeat   
	## get all impute2_info files for each chunk
	files <- system(paste0("ls *", suffix4impute2info), intern=TRUE) 
	for (i in seq_len(length(files))) { 
		## impute2infoAllvariants.txt is the temporal file
	 	system(paste0("sed 1d ", files[i], "  >> impute2infoAllvariants.txt")) 
	}  
	 
	## only keep SNPs and SNPs with two alleles  
	## impute2infoUpdateTmp.txt > temporal file
	arg1 = paste0("grep 'rs' impute2infoAllvariants.txt ")
	arg2 = paste0("awk '{if(length($4) == 1 && length($5) == 1) print}' ")
	arg3 = paste0("awk '{print $2, $7}' > impute2infoUpdateTmp.txt")
	system(paste0(arg1, " | ", arg2, " | ", arg3))
	# system("grep 'rs' impute2infoAllvariants.txt | 
	# awk '{if(length($4) == 1 && length($5) == 1) print}' | 
	# awk '{print $2, $7}' > impute2infoUpdateTmp.txt") 
	system(paste0("mv impute2infoUpdateTmp.txt ", outputInfoFile))
	## added colnames 
	impute2info <- read.table(file=outputInfoFile, stringsAsFactors=FALSE)  
	colnames(impute2info) <- c("rs_id", "info") 

	##  filtering   
	snpWithBadInfo <- impute2info[which(impute2info[, "info"] < infoScore), 1]  
	write.table(snpWithBadInfo, file=badImputeSNPfile, quote=FALSE, 
				row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	## extract filtered SNPs  
	system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
		   badImputeSNPfile, " --make-bed --out ", outputPrefix)) 
	system("rm impute2infoAllvariants.txt")
	
}


##########################################################################
## removedSnpMissPostImp.R
########################################################################## 
#' Remove SNPs after post imputation  
#'
#' @description
#' Remove SNPs which have a non missing value for less than a predefined 
#' number of instances.    
#' 
#' @param plink an executable program in either the current working 
#' directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param missCutoff  the cutoff of the least number of instances for 
#' a SNP that is not missing.
#' @param outputSNPfile the output file of SNPs with pre-defined 
#' missing values.
#' @param outputPrefix  the prefix of the PLINK format files. 

#' @return  The PLINK format files after post imputation quality control 
#' and a pure text file contains SNPs with pre-defined missing values.
#' @export 

#' @author Junfang Chen 


removedSnpMissPostImp <- function(plink, inputPrefix, missCutoff, 
								  outputSNPfile, outputPrefix){ 

	## get the missing info 
	system(paste0(plink, " --bfile ", inputPrefix, 
		   " --missing --out ", inputPrefix)) 

	missSNPinfo <- read.table(paste0(inputPrefix, ".lmiss"), 
							  stringsAsFactors=FALSE, header=TRUE)
	missSNPinfo[,6] <- missSNPinfo[,"N_GENO"] - missSNPinfo[,"N_MISS"] 
	snpWithManyMissSNPs <- missSNPinfo[which(missSNPinfo[,6] < missCutoff), "SNP"] 
	write.table(snpWithManyMissSNPs, file=outputSNPfile, quote=FALSE, 
				row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
		   outputSNPfile, " --make-bed --out ", outputPrefix))
	system( "rm *.imiss *.lmiss *.log") 

}
# 

 
 