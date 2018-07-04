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

#' Exclude monomorphic SNPs 
#'
#' @description
#' Detect monomorphic SNPs from PLINK BIM file and exclude them if any.

#' @param plink an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files 
#' after removing monomorphic SNPs.
#' @param outputSNPfile the output pure text file that stores 
#' the removed monomorphic SNPs, one per line, if any.

#' @return The output PLINK binary files after removing monomorphic SNPs 
#' and a pure text file with removed monomorphic SNPs.
#' @export 

#' @author Junfang Chen 
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "alignedData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "alignedData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "alignedData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, " ."))   
#' system(paste0("scp ", bimFile, " ."))   
#' system(paste0("scp ", famFile, " ."))   
#' inputPrefix <- "alignedData"
#' outputPrefix <- "removedMonoSnp" 
#' outputSNPfile <- "monoSNP.txt"  
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedMonoSnp(plink, inputPrefix, outputPrefix, outputSNPfile)

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

#' Split genome-wide genotyping data into chromosome-wide PLINK binary files.
#'
#' @description
#' Split the whole genome genotyping data chromosome-wise; 
#' allow parallel computating for all chromosomes.

#' @param plink an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files 
#' before splitting. 
#' @param chrXPAR1suffix  if chromosome 25 is available and with PAR1, 
#' then generate the suffix with X_PAR1 for chrX_PAR1. 
#' @param chrXPAR2suffix  if chromosome 25 is available and with PAR2, 
#' then generate the suffix with X_PAR2 for chrX_PAR2.
#' @param nCore the number of cores used for parallel computation. 
#' The default value is 25.  

#' @return The output PLINK binary files for each chromosome with the same 
#' prefix as the inputPrefix but appended with the chromosome codes, and  
#' possibly also the logical value for the pseudo-autosomal region (PAR)
#' indicating if PAR exists in the input genotyping data or not.   

#' @details If chromosome 25 is also available, namely the pseudo-autosomal
#' region of chromosome X, then further split chr25 (PAR or Chr_XY)
#' into PAR1 and PAR2 according to the genomic coordination GRCh37
#' from https://en.wikipedia.org/wiki/Pseudoautosomal_region.
#' The locations of the PARs within GRCh37 are:  
#' PAR1    X    60001    2699520; 
#' PAR2    X    154931044    155260560.   
#' @export 
#' @import doParallel  

#' @author Junfang Chen
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "alignedData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "alignedData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "alignedData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, " ."))   
#' system(paste0("scp ", bimFile, " ."))   
#' system(paste0("scp ", famFile, " ."))    
#' inputPrefix <- "alignedData"
#' chrXPAR1suffix <- "X_PAR1"
#' chrXPAR2suffix <- "X_PAR2"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## chrWiseSplit(plink, inputPrefix, chrXPAR1suffix, chrXPAR2suffix, nCore)
 

chrWiseSplit <- function(plink, inputPrefix, chrXPAR1suffix, 
                         chrXPAR2suffix, nCore=25){ 

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
        bim25 <- read.table(paste0(inputPrefix, "25.bim"), 
                            stringsAsFactors=FALSE) 
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
                   inputPrefix, chrXPAR1suffix))
            par1 <- TRUE  
            ## check for PAR2, if any SNPs out of PAR1, PAR2 also available 
            if (length(bimPos4par1)<nrow(bim25)){
                print("PAR2 is available in chrX!") 
                ## just to exclude rs4PAR1.txt
                system(paste0(plink, " --bfile ", inputPrefix, 
                       "25 --exclude rs4PAR1.txt --make-bed --out ", 
                       inputPrefix, chrXPAR2suffix))
                par2 <- TRUE 
            } else { 
                par2 <- FALSE
                print("PAR2 is NOT available in chrX!") 
            }

        } else { 
            print("PAR2 is available in chrX! But NOT PAR1, all chr25 on PAR2")
            par1 <- FALSE
            par2 <- TRUE 
            renamePlinkBFile(inputPrefix=paste0(inputPrefix, "25"), 
                             outputPrefix=paste0(inputPrefix, chrXPAR2suffix), 
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

#' @param inputPrefix the prefix of the input PLINK .bim file for 
#' each chromosome, without the chromosome codes.
#' @param outputPrefix the prefix of the output pure text files that keep 
#' all the chunks for each chromosome separately.
#' @param chrs specifiy the chromosome codes for chunking.
#' @param windowSize the window size of each chunk. 
#' The default value is 3000000.

#' @return The output pure text files include all the chunks 
#' for each chromosome separately. 

#' @export 
#' @author Junfang Chen 
#' @examples  
#' ## In the current working directory
#' bimFile <- system.file("extdata", "gwas_data_chr23.bim", package="Gimpute") 
#' system(paste0("scp ", bimFile, " ."))    
#' inputPrefix <- "gwas_data_chr"
#' outputPrefix <- "chunks_chr"
#' chrs <- 23
#' print(chrs)   
#' chunk4eachChr(inputPrefix, outputPrefix, chrs, windowSize=3000000)

chunk4eachChr <- function(inputPrefix, outputPrefix, chrs, windowSize=3000000){  

    for (i in chrs){ 

        bim <- paste0(inputPrefix, i, ".bim")
        bimdata <- read.table(file=bim, sep="\t", stringsAsFactors=FALSE)
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
            ## it may happen that only a few SNPs from the last chunk; 
            ## if the last chunk is large, then specify -allow_large_regions 
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
## .prePhasingByShapeit.R
########################################################################## 
#' Prephasing genotypes using SHAPEIT
#'
#' @description
#' Perform prephasing for study genotypes by SHAPEIT for the autosomal and 
#' sex chromosome haplotypes using a reference panel (pre-set).

#' @param shapeit an executable program in either the 
#' current working directory or somewhere in the command path.
#' @param chrs specifiy the chromosome codes for phasing.
#' @param dataDIR the directory where genotype PLINK binary files are located.
#' @param prefix4eachChr the prefix of PLINK binary files for 
#' each chromosome.
#' @param referencePanel a string indicating the type of imputation 
#' reference panels is used: c("1000Gphase1v3_macGT1", "1000Gphase3").
#' @param impRefDIR the directory where the imputation reference files 
#' are located.
#' @param phaseDIR the directory where resulting pre-phased files 
#' will be located.
#' @param nThread the number of threads used for computation.
#' The default value is 40. 
#' @param effectiveSize this parameter controls the effective population size. 
#' The default value is 20000.
#' @param nCore the number of cores used for computation. This can be tuned 
#' along with nThread. The default value is 1

#' @return The pre-phased haplotypes for given chromosomes.  
#' @details If ChrX is available then it is done differently by passing 
#' the flag --chrX to SHAPEIT.

#' @export 
 
#' @author Junfang Chen    
#' @seealso \code{\link{phaseImpute2}}.


.prePhasingByShapeit <- function(shapeit, chrs, dataDIR, prefix4eachChr,  
                                 referencePanel, impRefDIR, phaseDIR, 
                                 nThread=40, effectiveSize=20000, nCore=1){

    chrslist <- as.list(chrs)
    mclapply(chrslist, function(i){

        # GWAS data files in PLINK binary format
        GWASDATA_BED <- paste0(dataDIR, prefix4eachChr, i, ".bed ") 
        GWASDATA_BIM <- paste0(dataDIR, prefix4eachChr, i, ".bim ")
        GWASDATA_FAM <- paste0(dataDIR, prefix4eachChr, i, ".fam ")
        ## >>> ref panel  
        GENMAP_FILE <- paste0(impRefDIR, "genetic_map_chr", i, 
                              "_combined_b37.txt ")
        GENMAP.chrXnonPAR <- paste0(impRefDIR, "genetic_map_chr",
                              "X_nonPAR_combined_b37.txt ")
        GENMAP.chrXPAR1 <- paste0(impRefDIR, "genetic_map_chr", 
                              "X_PAR1_combined_b37.txt ")
        GENMAP.chrXPAR2 <- paste0(impRefDIR, "genetic_map_chr", 
                              "X_PAR2_combined_b37.txt ")

        if (referencePanel == "1000Gphase1v3_macGT1"){  
            ## autosome
            HAPS_FILE <- paste0(impRefDIR, 
                                "ALL_1000G_phase1integrated_v3_chr", i, 
                                "_impute_macGT1.hap.gz ") 
            LEGEND_FILE <- paste0(impRefDIR, 
                                  "ALL_1000G_phase1integrated_v3_chr", i, 
                                  "_impute_macGT1.legend.gz ") 
            SAMPLE_FILE <- paste0(impRefDIR, 
                                  "ALL_1000G_phase1integrated_v3.sample ") 
            ## chrX_nonPAR
            HAPS.chrXnonPAR <- paste0(impRefDIR, 
                                      "ALL_1000G_phase1integrated_v3_chr", 
                                      "X_nonPAR_impute_macGT1.hap.gz ") 
            LEGEND.chrXnonPAR <- paste0(impRefDIR, 
                                        "ALL_1000G_phase1integrated_v3_chr", 
                                        "X_nonPAR_impute_macGT1.legend.gz ") 
            ## .chrXPAR1  
            HAPS.chrXPAR1 <- paste0(impRefDIR, 
                                "ALL_1000G_phase1integrated_v3_chr",
                                "X_PAR1_impute_macGT1.hap.gz ") 
            LEGEND.chrXPAR1  <- paste0(impRefDIR, 
                                  "ALL_1000G_phase1integrated_v3_chr", 
                                  "X_PAR1_impute_macGT1.legend.gz ") 
            ## .chrXPAR2
            HAPS.chrXPAR2  <- paste0(impRefDIR, 
                                "ALL_1000G_phase1integrated_v3_chr", 
                                "X_PAR2_impute_macGT1.hap.gz ") 
            LEGEND.chrXPAR2  <- paste0(impRefDIR, 
                                  "ALL_1000G_phase1integrated_v3_chr", 
                                  "X_PAR2_impute_macGT1.legend.gz ") 

        } else if (referencePanel == "1000Gphase3"){
            HAPS_FILE <- paste0(impRefDIR, "1000GP_Phase3_chr", i, 
                                ".hap.gz ") 
            ## autosome
            LEGEND_FILE <- paste0(impRefDIR, "1000GP_Phase3_chr", i, 
                                  ".legend.gz ")
            SAMPLE_FILE <- paste0(impRefDIR, "1000GP_Phase3.sample ")

            ## chrX_nonPAR
            HAPS.chrXnonPAR <- paste0(impRefDIR, 
                                      "1000GP_Phase3_chrX_NONPAR.hap.gz ") 
            LEGEND.chrXnonPAR <- paste0(impRefDIR, 
                                        "1000GP_Phase3_chr", 
                                        "X_NONPAR.legend.gz ") 
            ## .chrXPAR1  
            HAPS.chrXPAR1 <- paste0(impRefDIR, 
                                    "1000GP_Phase3_chrX_PAR1.hap.gz ") 
            LEGEND.chrXPAR1 <- paste0(impRefDIR, 
                                      "1000GP_Phase3_chrX_PAR1.legend.gz ") 
            ## .chrXPAR2
            HAPS.chrXPAR2 <- paste0(impRefDIR, 
                                    "1000GP_Phase3_chrX_PAR2.hap.gz ") 
            LEGEND.chrXPAR2 <- paste0(impRefDIR, 
                                      "1000GP_Phase3_chrX_PAR2.legend.gz ") 
        } else { print("Wrong reference panel during phasing!!")}
        ## <<< ref panel 
        # main output file
        OUTPUT_HAPS <- paste0(phaseDIR, "chr", i, ".haps ")     
        OUTPUT_SAMPLE <- paste0(phaseDIR, "chr", i, ".sample ")     
        OUTPUT_LOG <- paste0(phaseDIR, "chr", i, ".log ")    

        autosomeCode = seq_len(22)
        if (is.element(i, autosomeCode)) { ## prePhasing for the autosome
            system(paste0(shapeit, 
            " --input-bed ", GWASDATA_BED, GWASDATA_BIM, GWASDATA_FAM, " \ ", 
            " --input-map ", GENMAP_FILE, " \ ",  
            "--input-ref ", HAPS_FILE, LEGEND_FILE, SAMPLE_FILE, " \ ", 
            "--thread ", nThread, " \ ", 
            "--effective-size ", effectiveSize, " \ ", 
            "--output-max ", OUTPUT_HAPS, OUTPUT_SAMPLE, " \ ", 
            "--output-log ", OUTPUT_LOG) )
        } else if (i == "X_PAR1"){
            system(paste0(shapeit, 
            " --input-bed ", GWASDATA_BED, GWASDATA_BIM, GWASDATA_FAM, " \ ", 
            " --input-map ", GENMAP.chrXPAR1, " \ ",  
            "--input-ref ", HAPS.chrXPAR1, LEGEND.chrXPAR1, SAMPLE_FILE, " \ ", 
            "--thread ", nThread, " \ ", 
            "--effective-size ", effectiveSize, " \ ", 
            "--output-max ", OUTPUT_HAPS, OUTPUT_SAMPLE, " \ ", 
            "--output-log ", OUTPUT_LOG) )       
        } else if (i == "X_PAR2"){
            system(paste0(shapeit, 
            " --input-bed ", GWASDATA_BED, GWASDATA_BIM, GWASDATA_FAM, " \ ", 
            " --input-map ", GENMAP.chrXPAR2, " \ ",  
            "--input-ref ", HAPS.chrXPAR2, LEGEND.chrXPAR2, SAMPLE_FILE, " \ ", 
            "--thread ", nThread, " \ ", 
            "--effective-size ", effectiveSize, " \ ", 
            "--output-max ", OUTPUT_HAPS, OUTPUT_SAMPLE, " \ ", 
            "--output-log ", OUTPUT_LOG) )       
        } else if (i == 23){
            system(paste0(shapeit, 
            " --input-bed ", GWASDATA_BED, GWASDATA_BIM, GWASDATA_FAM, " \ ", 
            " --input-map ", GENMAP.chrXnonPAR, " \ ", 
            " --chrX \ ",  ## special case
            "--input-ref ", HAPS.chrXnonPAR, LEGEND.chrXnonPAR, SAMPLE_FILE, " \ ", 
            "--thread ", nThread, " \ ", 
            "--effective-size ", effectiveSize, " \ ", 
            "--output-max ", OUTPUT_HAPS, OUTPUT_SAMPLE, " \ ", 
            "--output-log ", OUTPUT_LOG) ) 
        }

    }, mc.cores=nCore)     
} 





##########################################################################
## .imputedByImpute2.R
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
#' @param referencePanel a string indicating the type of imputation 
#' reference panels is used: c("1000Gphase1v3_macGT1", "1000Gphase3").
#' @param impRefDIR the directory where the imputation reference files 
#' are located.
#' @param imputedDIR the directory where imputed files will be located.
#' @param prefix4eachChr the prefix of IMPUTE2 files for each chunk.
#' @param nCore the number of cores used for computation.
#' @param effectiveSize this parameter controls the effective population size.
#' Commonly denoted as Ne. A universal -Ne value of 20000 is suggested.
#' @param XPAR a logical value indicating whether --chrX flag should be 
#' passed for prephasing using SHAPEIT.
#' --chrX flag, specifically for chrX imputation'
#' @return The imputed files for all chunks from given chromosomes.  
#' @export 
#' @import doParallel  

#' @author Junfang Chen 
#' @seealso \code{\link{phaseImpute2}}.



.imputedByImpute2 <- function(impute2, chrs, prefixChunk, phaseDIR, 
                              referencePanel, impRefDIR, imputedDIR, 
                              prefix4eachChr, nCore, effectiveSize=20000){ 

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

            ## >>> ref panel  
            GENMAP_FILE <- paste0(impRefDIR, "genetic_map_chr", i, 
                                  "_combined_b37.txt ")
            GENMAP.chrXnonPAR <- paste0(impRefDIR, "genetic_map_chr",
                                  "X_nonPAR_combined_b37.txt ")
            GENMAP.chrXPAR1 <- paste0(impRefDIR, "genetic_map_chr", 
                                  "X_PAR1_combined_b37.txt ")
            GENMAP.chrXPAR2 <- paste0(impRefDIR, "genetic_map_chr", 
                                  "X_PAR2_combined_b37.txt ")
 
            if (referencePanel == "1000Gphase1v3_macGT1"){  
                ## autosome
                HAPS_FILE <- paste0(impRefDIR, 
                                    "ALL_1000G_phase1integrated_v3_chr", i, 
                                    "_impute_macGT1.hap.gz ") 
                LEGEND_FILE <- paste0(impRefDIR, 
                                      "ALL_1000G_phase1integrated_v3_chr", i, 
                                      "_impute_macGT1.legend.gz ") 
                ## chrX_nonPAR
                HAPS.chrXnonPAR <- paste0(impRefDIR, 
                                          "ALL_1000G_phase1integrated_v3_chr", 
                                          "X_nonPAR_impute_macGT1.hap.gz ") 
                LEGEND.chrXnonPAR <- paste0(impRefDIR, 
                                            "ALL_1000G_phase1integrated_v3_chr", 
                                            "X_nonPAR_impute_macGT1.legend.gz ") 
                ## .chrXPAR1  
                HAPS.chrXPAR1 <- paste0(impRefDIR, 
                                    "ALL_1000G_phase1integrated_v3_chr",
                                    "X_PAR1_impute_macGT1.hap.gz ") 
                LEGEND.chrXPAR1  <- paste0(impRefDIR, 
                                      "ALL_1000G_phase1integrated_v3_chr", 
                                      "X_PAR1_impute_macGT1.legend.gz ") 
                ## .chrXPAR2
                HAPS.chrXPAR2  <- paste0(impRefDIR, 
                                    "ALL_1000G_phase1integrated_v3_chr", 
                                    "X_PAR2_impute_macGT1.hap.gz ") 
                LEGEND.chrXPAR2  <- paste0(impRefDIR, 
                                      "ALL_1000G_phase1integrated_v3_chr", 
                                      "X_PAR2_impute_macGT1.legend.gz ") 

            } else if (referencePanel == "1000Gphase3"){
                HAPS_FILE <- paste0(impRefDIR, "1000GP_Phase3_chr", i, 
                                    ".hap.gz ") 
                ## autosome
                LEGEND_FILE <- paste0(impRefDIR, "1000GP_Phase3_chr", i, 
                                      ".legend.gz ") 
                ## chrX_nonPAR
                HAPS.chrXnonPAR <- paste0(impRefDIR, 
                                          "1000GP_Phase3_chrX_NONPAR.hap.gz ") 
                LEGEND.chrXnonPAR <- paste0(impRefDIR, 
                                            "1000GP_Phase3_chr", 
                                            "X_NONPAR.legend.gz ") 
                ## .chrXPAR1  
                HAPS.chrXPAR1 <- paste0(impRefDIR, 
                                        "1000GP_Phase3_chrX_PAR1.hap.gz ") 
                LEGEND.chrXPAR1 <- paste0(impRefDIR, 
                                          "1000GP_Phase3_chrX_PAR1.legend.gz ") 
                ## .chrXPAR2
                HAPS.chrXPAR2 <- paste0(impRefDIR, 
                                        "1000GP_Phase3_chrX_PAR2.hap.gz ") 
                LEGEND.chrXPAR2 <- paste0(impRefDIR, 
                                          "1000GP_Phase3_chrX_PAR2.legend.gz ") 
            } else { print("Wrong reference panel during imputation!!")}
            ## <<< ref panel 

            ## main output file    
            OUTPUT_FILE <- paste0(imputedDIR, prefix4eachChr, i, 
                                  ".pos", chunkSTART, 
                                  "-", chunkEND, ".impute2 ")   
            ##  impute genotypes from GWAS haplotypes 
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
            } else if (i == "X_PAR1"){  
                ## impute for chrX PAR >> with an additional flag: --Xpar.
                system(paste0(impute2, 
                " -iter 30  \ ", 
                " -burnin 10  \ ", 
                " -k_hap 500  \ ", 
                " -use_prephased_g  \ ", 
                " -Xpar \ ",     ########## special
                " -m ", GENMAP.chrXPAR1, " \ ",  
                " -h ", HAPS.chrXPAR1, " \ ", 
                " -l ", LEGEND.chrXPAR1, " \ ", 
                " -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
                " -Ne ", effectiveSize, " \ ", 
                " -int ", chunkSTART, " ", chunkEND, " \ ", 
                " -buffer 1000  \ ",
                " -o ", OUTPUT_FILE, " \ ", 
                " -allow_large_regions \ ",
                " -seed 367946 \ " ))
            } else if (i == "X_PAR2"){  
                ## impute for chrX PAR >> with an additional flag: --Xpar.
                system(paste0(impute2, 
                " -iter 30  \ ", 
                " -burnin 10  \ ", 
                " -k_hap 500  \ ", 
                " -use_prephased_g  \ ", 
                " -Xpar \ ",     ########## special
                " -m ", GENMAP.chrXPAR2, " \ ",  
                " -h ", HAPS.chrXPAR1, " \ ", 
                " -l ", LEGEND.chrXPAR1, " \ ", 
                " -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
                " -Ne ", effectiveSize, " \ ", 
                " -int ", chunkSTART, " ", chunkEND, " \ ", 
                " -buffer 1000  \ ",
                " -o ", OUTPUT_FILE, " \ ", 
                " -allow_large_regions \ ",
                " -seed 367946 \ " ))
            } else if (i == 23){ 
                ## impute for chrX 
                ## >>  additional flag: --chrX, and sample_known_haps_g
                system(paste0(impute2, 
                " -iter 30  \ ", 
                " -burnin 10  \ ", 
                " -k_hap 500  \ ", 
                " -use_prephased_g  \ ", 
                " -chrX \ ",   ########### special
                " -m ", GENMAP.chrXnonPAR, " \ ",  
                " -h ", HAPS.chrXnonPAR, " \ ", 
                " -l ", LEGEND.chrXnonPAR, " \ ", 
                " -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
                " -sample_known_haps_g ", GWAS_SAMP_FILE, " \ ",   
                " -Ne ", effectiveSize, " \ ", 
                " -int ", chunkSTART, " ", chunkEND, " \ ", 
                " -buffer 1000  \ ",
                " -o ", OUTPUT_FILE, " \ ", 
                " -allow_large_regions \ ",
                " -seed 367946 \ " ))
            } else { print(" wrong chromosome code during phasing!!") }
        }, mc.cores=nCore)  
    } 
}







##########################################################################
## .convertImpute2ByGtool.R 
##########################################################################  
#' Convert IMPUTE2 format files into PLINK format
#'
#' @description
#' Convert all chunks of IMPUTE2 format files into binary PLINK format using 
#' GTOOL.

#' @param gtool an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param chrs specifiy the chromosome codes for conversion. 
#' @param prefixChunk  the prefix of the chunk files for each chromosome, 
#' along with the location directory.
#' @param phaseDIR the directory where pre-phased files are located.
#' @param imputedDIR the directory where the imputated files are located.
#' @param prefix4eachChr the prefix of the input IMPUTE2 files and 
#' also the output PLINK binary files for each chunk.
#' @param suffix4imputed the suffix of the IMPUTE2 format file that stores 
#' the imputed value.
#' @param postImputeDIR the directory where converted PLINK binary files 
#' will be located.  
#' @param threshold threshold for merging genotypes from GEN probability. 
#' Default 0.9. 
#' @param nCore the number of cores used for computation.  
 
#' @return The converted binary PLINK format files for each chunk from IMPUTE2 
#' results.
##' @export 
#' @import doParallel  

#' @author Junfang Chen 

 


.convertImpute2ByGtool <- function(gtool, chrs, prefixChunk, 
                                   phaseDIR, imputedDIR, prefix4eachChr, 
                                   suffix4imputed, postImputeDIR, 
                                   threshold, nCore){

    for (i in chrs){ 
        chunkfn <- paste0(prefixChunk, i, ".txt")
        chunks <- read.table(chunkfn, sep=" ") 
        chunklist <- as.list(seq_len(nrow(chunks)))

        mclapply(chunklist, function(j){ 
            chunkSTART <- chunks[j,1]
            chunkEND   <- chunks[j,2] 
            ## INPUT data files
            SAM_FILE <- paste0(phaseDIR, "chr", i, ".sample")  
            GEN_FILE <- paste0(imputedDIR, prefix4eachChr, i, 
                               ".pos", chunkSTART, 
                               "-", chunkEND, suffix4imputed) 
            ## output PLINK binary files
            PED_FILE <- paste0(postImputeDIR, prefix4eachChr, i, ".pos", 
                               chunkSTART, "-", chunkEND, ".ped") 
            MAP_FILE <- paste0(postImputeDIR, prefix4eachChr, i, ".pos", 
                                 chunkSTART, "-", chunkEND, ".map") 
            ## converting chunk-wise
            system( paste0(gtool, " -G ",
                            " --g ", GEN_FILE, " \ ", 
                            " --s ", SAM_FILE, " \ ",  
                            "--phenotype plink_pheno \ ", 
                            "--chr ", i, " \ ", 
                            "--ped ", PED_FILE, " \ ", 
                            "--map ", MAP_FILE, " \ ", 
                            "--threshold ", threshold, " \ ", 
                            "--snp ") )  
        }, mc.cores=nCore)
    } 
} 


 


##########################################################################
## .mergePlinkData.R
########################################################################## 
#' Merge chunk-wise PLINK files 
#'
#' @description
#' Merge all chunk-wise PLINK binary files into chromosome-wise PLINK 
#' binary files then assemble into a genome-wide PLINK binary file set. 

#' @param plink an executable program in either the current working 
#' directory or somewhere in the command path.
#' @param chrs specifiy the chromosome codes to be merged. 
#' @param prefix4eachChr the prefix of the input chunk-wise PLINK 
#' files. 
#' @param prefix4mergedPlink  the prefix of the final output PLINK 
#' binary files. 
#' @param nCore the number of cores used for computation.  

#' @return The merged genome-wide PLINK binary files.
#' @details Create a file containing a list chunk-wise PLINK PED and MAP 
#' file names. The prefix of these files must already indicate in which 
#' chromosome they belong to and files from the same chromosome will be 
#' combined. Then all chromosomal PLINK files are assembled together 
#' into one whole genome PLINK binary file set.

##' @export 
#' @import doParallel  
#' @author Junfang Chen 


.mergePlinkData <- function(plink, chrs, prefix4eachChr, 
                           prefix4mergedPlink, nCore){ 

    ## firstly, only consider chromosomes from 1:23; as Xpar chrs 
    ## are slightly different for processing.
    pureAutoChrs <- setdiff(chrs, c("X_PAR1", "X_PAR2")) 
    chrslist <- as.list(pureAutoChrs)   
    mclapply(chrslist, function(i){

        pedFile_chr <- system(paste0("ls ", prefix4eachChr, i, ".*.ped"), 
                              intern=TRUE)
        mapFile_chr <- system(paste0("ls ", prefix4eachChr, i, ".*.map"), 
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

        pedFile_chr <- system(paste0("ls ", prefix4eachChr, "X_PAR*.ped"), 
                              intern=TRUE)
        mapFile_chr <- system(paste0("ls ", prefix4eachChr, "X_PAR*.map"), 
                              intern=TRUE)    
        pedmap_chr <- paste0(pedFile_chr, " ", mapFile_chr)
        fA <- gsub(".ped", "", pedFile_chr[1])
        pedmap_tobeMerged <- pedmap_chr[-1]
        filesetname <- paste0("fileset_chr25.txt")
        write.table(pedmap_tobeMerged, file=filesetname, quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
         
        arg <- paste0(plink, " --file ", fA, " --merge-list ", filesetname, 
                      " --allow-extra-chr --make-bed --out gwasImputedOld25")
        system(arg) 
        ## update chr code for XPAR --> 25
        bim <- read.table("gwasImputedOld25.bim", stringsAsFactors=FALSE)
        updateSNPchr <- cbind(bim[,2], rep(25, length=nrow(bim))) 
        write.table(updateSNPchr, file="gwasImputed_newchr25.txt", quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 


        system(paste0(plink, " --bfile gwasImputedOld25 --allow-extra-chr ", 
               " --update-chr gwasImputed_newchr25.txt 2 1",
               " --make-bed --out gwasImputed_chr25"))  
        # system("rm gwasImputedOld25.* gwasImputed_newchr25.txt")
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
## .filterImputeData.R
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
#' each variant. The default value is 0.6. 
#' @param badImputeSNPfile the output file of SNPs with bad info scores.  
#' @param inputPrefix the prefix of the input imputed PLINK binary files. 
#' @param outputPrefix the prefix of the output filtered PLINK binary files. 

#' @return A pure text file contains the info scores of all imputed SNPs with 
#' two columns: SNP names and the corresponding info scores. 
#' A pure text file with all excluded SNPs having bad info scores. 
#' The filtered PLINK binary imputed files, 
#' @details Filter genetic variants accoring to the imputation quality score 
#' with the help of .impute2_info files generated by IMPUTE2. 
#' Often, we keep variants with imputation info score of greater than 0.6.    
#' Note that imputed SNPs with more than two alleles are not considered. 

##' @export 
#' @author Junfang Chen 


.filterImputeData <- function(plink, suffix4impute2info, outputInfoFile, 
                             infoScore=0.6, badImputeSNPfile, inputPrefix, 
                             outputPrefix){ 

    ## read each .impute2_info file, remove 1st line, add to another file 
    ## and repeat get all impute2_info files for each chunk
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
    system(paste0("mv impute2infoUpdateTmp.txt ", outputInfoFile))
    ## added colnames 
    impute2info <- read.table(file=outputInfoFile, stringsAsFactors=FALSE)  
    colnames(impute2info) <- c("rs_id", "info") 

    ## filtering   
    snpWithBadInfo <- impute2info[which(impute2info[, "info"] < infoScore), 1]  
    write.table(snpWithBadInfo, file=badImputeSNPfile, quote=FALSE, 
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
    ## extract filtered SNPs  
    system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
           badImputeSNPfile, " --make-bed --out ", outputPrefix)) 
    # system("rm impute2infoAllvariants.txt")

}




#' Filter genetic variants    
#'
#' @description
#' Filter out genetic variants accoring to the info score.
#' 
#' @param plink an executable program in either the current working 
#' directory or somewhere in the command path. 
#' @param outputInfoFile the output file of info scores consisting of 
#' two columns: SNP names and their info scores.  
#' @param infoScore the cutoff of filtering imputation quality score for 
#' each variant. The default value is 0.6. 
#' @param badImputeSNPfile the output file of SNPs with bad info scores.  
#' @param inputPrefix the prefix of the input imputed PLINK binary files. 
#' @param outputPrefix the prefix of the output filtered PLINK binary files. 

#' @return A pure text file contains the info scores of all imputed SNPs with 
#' two columns: SNP names and the corresponding info scores. 
#' A pure text file with all excluded SNPs having bad info scores. 
#' The filtered PLINK binary imputed files, 
#' @details Filter genetic variants accoring to the imputation quality score 
#' with the help of .impute2_info files generated by IMPUTE2. 
#' Often, we keep variants with imputation info score of greater than 0.6.    
#' Note that imputed SNPs with more than two alleles are not considered. 

#' @export 
#' @author Junfang Chen 


.filterImputeData2 <- function(plink, outputInfoFile, infoScore=0.6, 
                               badImputeSNPfile, inputPrefix, 
                               outputPrefix){ 
 
    ## with colnames: rsid, info
    impute2info <- read.table(file=outputInfoFile, 
                              stringsAsFactors=FALSE, header=TRUE)
    impute2info[,2] <- round(impute2info[,2], 3)
 
    ## filtering   
    snpWithBadInfo <- impute2info[which(impute2info[, "info"] < infoScore), 1]  
    write.table(snpWithBadInfo, file=badImputeSNPfile, quote=FALSE, 
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
    ## extract filtered SNPs  
    system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
           badImputeSNPfile, " --make-bed --out ", outputPrefix))  

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
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param missCutoff  the cutoff of the least number of instances for 
#' a SNP that is not missing. The default is 20.
#' @param outputSNPfile the output file of SNPs with pre-defined 
#' missing values.
#' @param outputPrefix  the prefix of the PLINK binary files. 

#' @return The PLINK binary files after post imputation quality control 
#' and a pure text file contains SNPs with pre-defined missing values.

#' @export 
#' @author Junfang Chen 
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "alignedData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "alignedData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "alignedData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, " ."))   
#' system(paste0("scp ", bimFile, " ."))   
#' system(paste0("scp ", famFile, " ."))    
#' inputPrefix <- "alignedData" 
#' outputPrefix <- "removedSnpMissPostImp" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedSnpMissPostImp(plink, inputPrefix, missCutoff=20, outputPrefix)


removedSnpMissPostImp <- function(plink, inputPrefix, missCutoff, 
                                  outputSNPfile, outputPrefix){ 

    ## get the missing info 
    system(paste0(plink, " --bfile ", inputPrefix, 
           " --missing --out ", inputPrefix)) 

    missSNPinfo <- read.table(paste0(inputPrefix, ".lmiss"), 
                              stringsAsFactors=FALSE, header=TRUE)
    missSNPinfo[,6] <- missSNPinfo[,"N_GENO"] - missSNPinfo[,"N_MISS"] 
    manyMissSNPs <- missSNPinfo[which(missSNPinfo[,6] < missCutoff), "SNP"] 
    write.table(manyMissSNPs, file=outputSNPfile, quote=FALSE, 
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
    system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
           outputSNPfile, " --make-bed --out ", outputPrefix))
    system( "rm *.imiss *.lmiss *.log") 

}


##########################################################################
## phaseImpute.R
########################################################################## 
#' Phasing and imputation 
#'
#' @description
#' Perform phasing, imputation and conversion from IMPUTE2 format into PLINK 
#' binary files. 

#' @param inputPrefix the prefix of the input PLINK binary files for 
#' the imputation.
#' @param outputPrefix the prefix of the output PLINK binary files  
#' after imputation and filtering out bad imputed variants.
#' @param prefix4final the prefix of the output PLINK binary files  
#' after imputation.
#' @param plink an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param shapeit an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param impute2 an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param gtool an executable program in either the current 
#' working directory or somewhere in the command path. 
#' @param windowSize  the window size of each chunk. 
#' The default value is 3000000. 
#' @param effectiveSize this parameter controls the effective population size. 
#' Commonly denoted as Ne. A universal -Ne value of 20000 is suggested.
#' @param nCore4phase the number of cores used for phasing. This can be tuned 
#' along with nThread. The default value is 1 
#' @param nThread the number of threads used for computation.
#' The default value is 40. 
#' @param nCore4impute the number of cores used for imputation. 
#' The default value is 40. 
#' @param threshold threshold for merging genotypes from GEN probability. 
#' Default 0.9. 
#' @param nCore4gtool the number of cores used for computation. 
#' The default value is 40. 
#' @param infoScore the cutoff of filtering imputation quality score for 
#' each variant. The default value is 0.6. 
#' @param outputInfoFile the output file of impute2 info scores consisting of 
#' two columns: all imputed SNPs and their info scores.   
#' @param referencePanel a string indicating the type of imputation 
#' reference panels is used: c("1000Gphase1v3_macGT1", "1000Gphase3").
#' @param impRefDIR the directory where the imputation reference files 
#' are located.  
#' @param tmpImputeDir the name of the temporary directory used for 
#' storing phasing and imputation results.
#' @param keepTmpDir a logical value indicating if the directory 
#' 'tmpImputeDir' should be kept or not. The default is TRUE.

#' @return 1.) The filtered imputed PLINK binary files; 
#' 2.) The final PLINK binary files including bad imputed variants;
#' 3.) A pure text file contains the info scores of all imputed SNPs with 
#' two columns: SNP names and the corresponding info scores. 
#' @details The whole imputation process mainly consists of the following 
#' steps: 
#' 1.) Phasing the input PLINK data using an existing imputation reference;  
#' 2.) Imputing the input PLINK data using phased results and an existing 
#' reference data; 
#' 3.) Converting IMPUTE2 format data into PLINK format.
#' 4.) Combining all imputed data into whole-genome PLINK binary files.
#' 5.) Filtering out imputed variants with bad imputation quality. 
#' Parallel computing in R is supported.

#' @export 
#' @import doParallel 
#' @author Junfang Chen 

#' @references  
#' \enumerate{
#'   \item Howie, B., et al. (2012). Fast and accurate genotype imputation 
#'         in genome-wide association studies through pre-phasing. Nat Genet 
#'         44(8): 955-959.
#'   \item Howie, B. N., et al. (2009). A flexible and accurate genotype 
#'         imputation method for the next generation of genome-wide association 
#'         studies. PLoS Genet 5(6): e1000529.
#' }

#' @examples 
#' ## In the current working directory
#' bedFile <- system.file("extdata", "alignedData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "alignedData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "alignedData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, " ."))   
#' system(paste0("scp ", bimFile, " ."))   
#' system(paste0("scp ", famFile, " ."))   
#' inputPrefix <- "alignedData"  
#' outputPrefix <- "gwasImputedFiltered"
#' prefix4final <- "gwasImputed"   
#' outputInfoFile <- "infoScore.txt"
#' tmpImputeDir <- "tmpImpute"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## phaseImpute2(inputPrefix, outputPrefix, prefix4final,
#' ##             plink, shapeit, impute2, gtool, 
#' ##             windowSize=3000000, effectiveSize=20000, 
#' ##             nCore4phase=1, nThread=40, 
#' ##             nCore4impute=40, nCore4gtool=40, 
#' ##             infoScore=0.6, outputInfoFile, referencePanel, 
#' ##             impRefDIR, tmpImputeDir, keepTmpDir=TRUE)


phaseImpute2 <- function(inputPrefix, outputPrefix, prefix4final,
                        plink, shapeit, impute2, gtool, 
                        windowSize=3000000, effectiveSize=20000, 
                        nCore4phase=1, nThread=40, 
                        nCore4impute=40, threshold=0.9, 
                        nCore4gtool=40, infoScore=0.6, outputInfoFile,
                        referencePanel, impRefDIR, 
                        tmpImputeDir, keepTmpDir=TRUE){

    ## One must create directories for storing tmp imputation output files 
    ## The name of these directories must be fixed for the sake of 
    ## the subsequent steps.
    system(paste0("mkdir ", tmpImputeDir))
    setwd(tmpImputeDir) ## 
    ## sub-directories  
    system("mkdir 1-dataFiles")
    system("mkdir 2-chunkFile") 
    system("mkdir 3-phaseResults")
    system("mkdir 4-imputeResults")
    system("mkdir 5-postImpute")
    system("mkdir 6-finalResults")  
    # define directories
    dataDIR <- "1-dataFiles/"  
    chunkDIR <- "2-chunkFile/"
    phaseDIR <- "3-phaseResults/"  
    imputedDIR <- "4-imputeResults/"  
    postImputeDIR <- "5-postImpute/" 
    finalImputeDIR <- "6-finalResults/"  
    setwd("..")  
    ## step 2.1 
    ## copy plink files without monomorphic SNPs; prepare for the imputation.
    prefix4eachChr <- "gwas_data_chr"  
    system(paste0("scp ", inputPrefix, ".* ./", tmpImputeDir, "/", dataDIR))
    setwd(paste0("./", tmpImputeDir, "/", dataDIR))  
    renamePlinkBFile(inputPrefix, outputPrefix=prefix4eachChr, action="move")
    bimCurrent <- read.table(file=paste0(prefix4eachChr, ".bim"), 
                             stringsAsFactors=FALSE)  
    currentChr <- names(table(bimCurrent[,1]))
    print(currentChr)  
    chrXPAR1suffix <- "X_PAR1"
    chrXPAR2suffix <- "X_PAR2"
    ## nCore is chosen as the number of chromosomes available 
    PAR <- chrWiseSplit(plink, inputPrefix=prefix4eachChr, chrXPAR1suffix, 
                        chrXPAR2suffix, nCore=length(currentChr))
    print(PAR)  
    if (PAR[[1]]) {par1 <- "X_PAR1"} else {par1 <- NULL}
    if (PAR[[2]]) {par2 <- "X_PAR2"} else {par2 <- NULL}
    ## step 2.2
    chunkPrefix <- "chunks_chr" 
    chrs <- c(currentChr, par1, par2)    
    chunk4eachChr(inputPrefix=prefix4eachChr, 
                  outputPrefix=chunkPrefix, chrs, windowSize) 

    setwd("..") 
    system(paste0("mv ", dataDIR, chunkPrefix, "*.txt  ", chunkDIR)) 
    ## step 2.3     
    .prePhasingByShapeit(shapeit, chrs, dataDIR, 
                         prefix4eachChr, referencePanel, 
                         impRefDIR, phaseDIR, nThread, 
                         effectiveSize, nCore=nCore4phase)
    ## step 2.4   
    prefixChunk <- paste0(chunkDIR, chunkPrefix)        
    .imputedByImpute2(impute2, chrs, prefixChunk, phaseDIR, referencePanel, 
                      impRefDIR, imputedDIR, prefix4eachChr, 
                      nCore4impute, effectiveSize)
    ## step 2.5   
    ## extract only SNPs (without INDELs)
    #######################################################
    setwd(imputedDIR)  
    ## extract only SNPs starting with "rs";  .
    ls <- system("ls gwas*.impute2", intern=TRUE)
    snpPrefix <- "rs" 
    biglists <- as.list(ls)
    mclapply(biglists, function(i){ 
        arg1 <- paste0(" awk '{if(length($4) == 1 && length($5) == 1) print}'")
        arg2 <- paste0(i, "noINDEL.impute2")   
        system(paste0("grep '", snpPrefix, "' ", i, " | ", arg1, " > ", arg2))
    }, mc.cores=40) ## by default  
    setwd("..") 
    suffix4imputed <- ".impute2noINDEL.impute2"   
    .convertImpute2ByGtool(gtool, chrs, prefixChunk, phaseDIR, imputedDIR, 
                          prefix4eachChr, suffix4imputed, 
                          postImputeDIR, threshold, nCore4gtool)

    ## step 2.6  
    ####################################################### 
    ## Modify missing genotype format.
    setwd(postImputeDIR)   
    ## replace 'N' in the .ped files into 0 > missing values.
    chrslist <- as.list(chrs) 
    fn <- mclapply(chrslist, function(i){
        system(paste0("sed -i 's/N/0/g' ", prefix4eachChr, i, ".*ped "))
    }, mc.cores=length(chrslist)) 
    prefixMerge <- "gwasMerged" 
    .mergePlinkData(plink, chrs, prefix4eachChr, 
                    prefixMerge, nCore=length(chrslist))
    ## fam IDs may be changed: a.) if IDs have 'N'; 
    ## b.) IID, FID may be switched.
    ## >> update this as below 
    ## the original PLINK files before imputation
    setwd("..")
    system(paste0("scp ", dataDIR, prefix4eachChr, ".fam ", postImputeDIR)) 
    setwd(postImputeDIR) 
    ## update FAM IDs in the imputed PLINK files
    famOrig <- read.table(paste0(prefix4eachChr, ".fam"), 
                          stringsAsFactors=FALSE) 
    famImpute <- read.table(paste0(prefixMerge,".fam"), 
                            stringsAsFactors=FALSE) 
    ## changes ID codes for individuals specified in recoded.txt, 
    ## which should be in the format of 4 cols per row: 
    ## old FID, old IID, new FID, new IID, e.g.
    recodMat <- cbind(famImpute[,c("V1", "V2")], famOrig[,c("V1", "V2")]) 
    write.table(recodMat, file="recoded.txt", quote=FALSE, 
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")  

    system(paste0(plink, " --bfile ", prefixMerge, 
           " --update-ids recoded.txt --make-bed --out ", prefix4final)) 
    #######################################################
    setwd("..")
    system(paste0("mv ", postImputeDIR, prefix4final, "* ", finalImputeDIR)) 
    system(paste0("mv ", imputedDIR, "*.impute2_info ", finalImputeDIR)) 
    ## step 2.7    
    setwd(finalImputeDIR)
    suffix4impute2info <- ".impute2_info"
    badImputeSNPfile <- "badImputeSNPs.txt" 
    .filterImputeData(plink, suffix4impute2info, 
                     outputInfoFile, infoScore, badImputeSNPfile, 
                     inputPrefix=prefix4final, outputPrefix)
    setwd("..")
    setwd("..")
    if (keepTmpDir == FALSE){
        system(paste0("rm -r ", tmpImputeDir))
    } else if (keepTmpDir == TRUE){
        print("Keep the temporary imputation folder.")
    }

}











