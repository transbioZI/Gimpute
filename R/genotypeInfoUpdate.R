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
##########################################   
##########################################  
#' Remove duplicated sample IDs
#'
#' @description
#' Remove duplicated sample IDs from PLINK binary files. The duplicated IDs
#' are defined in a plain text file in advance. 

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param dupSampleIDFile a pure text file that stores the duplicated sample IDs, 
#' each ID per line. If it is null, then duplicate the input PLINK files from 
#' the last step as the output files.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after removing duplicated sample IDs. 
#' @export 

#' @author Junfang Chen 
#' @examples 

removeDupID <- function(plink, dupSampleIDFile, inputPrefix, outputPrefix){

    if (!is.null(dupSampleIDFile)){ 

        famv0 <- read.table(file=paste0(inputPrefix, ".fam"), 
                            stringsAsFactors=FALSE)  
        dupIDs <- read.table(file=dupSampleIDFile, stringsAsFactors=FALSE)  
        ## write into plink format .txt for removal 
        dupIDs4plink <- famv0[is.element(famv0[,2], dupIDs[,1]), c("V1", "V2")] 
        dupIDsfn <- paste0(outputPrefix, ".txt")
        write.table(dupIDs4plink, file=dupIDsfn, quote=FALSE, row.names=FALSE, 
                    col.names=FALSE, eol="\r\n", sep=" ") 
        system(paste0(plink, " --bfile ", inputPrefix, 
               " --remove ", dupIDsfn, " --make-bed --out ", outputPrefix))
        system(paste0("rm ", dupSampleIDFile, " ", dupIDsfn)) 
        # system(paste0("rm ", inputPrefix, ".*") )

    } else { 
        ## copy/rename all snp info updated plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy")
        ## remove the raw genotype plink files
        # system(paste0("rm ", inputPrefix, ".*"))
    } 
}



##########################################   
##########################################   
#' Update group and geneder information
#'
#' @description
#' Replace group and gender information in the PLINK binary files by using 
#' the information from the metadata file.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param metaDataFile a pure text file that stores the meta information of 
#' the samples.  
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after updating the gender and grouping 
#' information. 

#' @details Find the shared sample IDs between PLINK input files and metadata file. 
#' Use the information from the metadata file as the reference and update the group 
#' information/the outcome in the PLINK file. Group label should be 1 and 2. 
#' (1=unaff, 2=aff, 0=miss); missing phenotype will be indicated as -9.

#' @export  
#' @author Junfang Chen 



updateGroupIdAndSex <- function(plink, inputPrefix, metaDataFile, outputPrefix){
         
    fam <- read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
    metaData <- read.table(metaDataFile, stringsAsFactors=FALSE, header=TRUE) 
    ## to make sure only compare with the same set of plIndID --> IID
    interIDs <- intersect(fam[,2], metaData[,"IID"])     # the shared IDs 
    metadataSub <- metaData[is.element(metaData[,"IID"], interIDs),]

    ## one must keep the order of .fam unchanged!
    metadataSubsort <- metadataSub[match(fam[,2], metadataSub[,"IID"]),]
    ## IID is in the 1st col; 
    missIDsIndex <- which(is.na(metadataSubsort[,1]) == TRUE) 
    fam[,6] <- metadataSubsort[,"group"] + 1 # update to the correct format.
    ## keep the IDs of no missing group info as -9
    fam[missIDsIndex, 6] <- -9  
    ## a pheno file that contains 3 columns (one row per individual)
    newPheno <- fam[,c(1,2,6)]  
    write.table(newPheno, file="pheno.txt", quote=FALSE, row.names=FALSE, 
                col.names=FALSE, eol="\r\n", sep=" ")  

    ## Alternate phenotype files
    outputPrefixVgroup <- paste0(outputPrefix, "_chGroup")
    system(paste0(plink, " --bfile ", inputPrefix, 
           " --pheno pheno.txt --make-bed --out ", outputPrefixVgroup) )
    system("rm pheno.txt" ) 

    ## update the sex information according to the metaData file 
    ## keep the IDs of no missing sex info as 0
    fam[,5] <- metadataSubsort[,"sex"]  
    fam[missIDsIndex, 5] <- 0  
    ##  a pheno file that contains 3 columns (one row per individual)
    updateSex <- fam[,c(1,2,5)]  
    write.table(updateSex, file="updateSex.txt", quote=FALSE, 
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")  
 
    ## --update-sex expects a file with FIDs and IIDs in the first two columns, 
    ## the 3rd column is the sex information.
    system(paste0(plink, " --bfile ", outputPrefixVgroup, 
           " --update-sex updateSex.txt --make-bed --out ", outputPrefix) )
    system("rm updateSex.txt" )
    system(paste0("rm ", outputPrefixVgroup, "*"))
}




##########################################   
########################################## 
#' Remove samples without group information
#'
#' @description
#' Remove samples without group/outcome/phenotype information, which is coded 
#' as -9 in the PLINK .FAM file.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param metaDataFile a pure text file that stores the meta information of 
#' the samples.  
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after removing samples without group 
#' information.
#' @export 

#' @author Junfang Chen 


removeNoGroupId <- function(plink, inputPrefix, outputPrefix){
 
    fam <- read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
    ## 1. check if any sample without phenotypes
    noGroupIds <- fam[which(fam[,6] == -9), c("V1", "V2")] 
    ## if any then remove and afterwards add phenotypes
    noGroupIdsfn <- paste0(outputPrefix, ".txt")
    write.table(noGroupIds, file=noGroupIdsfn, quote=FALSE, 
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")   
    system(paste0(plink, " --bfile ", inputPrefix, " --remove ", 
           noGroupIdsfn, " --make-bed --out ", outputPrefix) )  
    system(paste0("rm ",  noGroupIdsfn) )
}


##########################################   
########################################## 
#' Remove samples with incorrect ancestry
#'
#' @description
#' Remove samples with the incorrect ancestry or keep samples at your own chioce.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param metaDataFile a pure text file that stores the meta information of 
#' the samples.  
#' @param ancestrySymbol an indicator that shows the symbol of genetic ancestry. 
#' If it is null, then all samples are selected. 
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after checking the ancestry information. 
#' @details ancestrySymbol, such as 'EA' stands for the European, 
#' 'AA' for African American. See the metaDataFile for more details. 
#' @export  
#' @author Junfang Chen  


removedWrongAnceInst <- function(plink, inputPrefix, metaDataFile, 
                                  ancestrySymbol, outputPrefix){

    if (is.null(ancestrySymbol)) { 
        ## copy/rename all snp info updated plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
    } else {
        metaData <- read.table(metaDataFile, stringsAsFactors=FALSE, header=TRUE) 
        ids <- metaData[which(metaData[,"ance"] == ancestrySymbol),"IID"]
        fam <- read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
        famEA <- fam[is.element(fam[,2], ids), c("V1", "V2")]
        plinkFormatDat <- famEA
        write.table(plinkFormatDat, file=paste0(outputPrefix,".txt"), quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")  

        system(paste0(plink, " --bfile ", inputPrefix, " --keep ", 
               paste0(outputPrefix,".txt"), " --make-bed --out ", outputPrefix) )
        system(paste0("rm ",  outputPrefix, ".txt") )
    }    

}



##########################################   
##########################################
#' Remove improper SNPs
#'
#' @description
#' Remove SNPs that may be duplicated, or with unexpected SNP names.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param excludedProbeIdsFile a pure text file that stores the SNP IDs, 
#' one per line, which need to be removed. If it is null, then no SNPs are removed.

#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after removing unwanted SNP IDs.
#' @details excludedProbeIdsFile should be defined in a plain text file 
#' in advance. Improper SNPs such as AFFX, cnvi etc.. with unexpect format
#' must be excluded.

#' @export 
#' @author Junfang Chen 
 
 
removedExclProbe <- function(plink, inputPrefix, excludedProbeIdsFile, outputPrefix){

    if (!is.null(excludedProbeIdsFile)) {
        system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
               excludedProbeIdsFile, " --make-bed --out ", outputPrefix) )
        ## remove .txt
        system(paste0("rm ", excludedProbeIdsFile) )
    } else { 
        ## copy/rename all snp info updated plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy")   
    }
}


##########################################   
##########################################   
#' Remove SNPs not in the chip annotation file

#' @description
#' Check if all SNPs are included in the chip annotation file. 
#' If some of the input SNPs are not included, then remove them.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param chipAnnoFile a pure text file that stores the chip annotation
#' information. 
#' @param outputSNPfile a pure text file that stores the SNP IDs, 
#' one per line, which are not mapped to the chip annotation file.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output text file contains the removed SNP IDs, one per line.
#' The PLINK binary files after removing unmapped SNP IDs.
#' @details If the chip annotation file is not available for your study, 
#' you can download it from http://www.well.ox.ac.uk/~wrayner/strand/.

#' @export 
#' @author Junfang Chen 


removedUnmapProbes <- function(plink, inputPrefix, chipAnnoFile, 
                               outputPrefix, outputSNPfile){

    if (!is.null(chipAnnoFile)){ 

        annoFile <- "chipAnnoRefb37.txt"
        if (chipType == "affymetrix") { 
            prepareChipAnnoFile4affymetrix(inputFile=chipAnnoFile, outputFile=annoFile)
        } else if (chipType == "illumina"){ 
            prepareChipAnnoFile4Illumina(inputFile=chipAnnoFile, outputFile=annoFile)
        } else if (chipType == "PsychChip"){ 
            prepareChipAnnoFile4PsychChip(inputFile=chipAnnoFile, outputFile=annoFile)
        }

        ## find the overlapping
        chipAnno <- read.table(file=annoFile, header=TRUE, stringsAsFactors=FALSE)
        ## check the overlapping SNPs
        bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
        interSNPs <- intersect(bim[,2], chipAnno[,"chipSnpID"])  
        unmapped <- setdiff(bim[,2], interSNPs) 
        write.table(unmapped, file=outputSNPfile, quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
        system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", outputSNPfile, 
               " --make-bed --out ", outputPrefix) )
    } else { 
        ## copy/rename plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy")  
    }
}



##########################################   
########################################## 
#' Remove duplicated SNPs
#'
#' @description
#' Remove duplicated SNPs that have same rs-names or duplicated genomic position.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param chipAnnoFile an input chip annotation file. 
#' If the chip annotation file is not available for your study, it can be 
#' downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.
#' @param chipType a string name that defines the type of the chip annotation 
#' file: 'illumina', 'affymetrix' or 'psychChip'.
#' @param outputSNPdupFile a pure text file that stores the duplicated SNP IDs, 
#' which are detected by the use of the chip annotation file.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after removing duplicated SNP IDs.
#' @details Duplicated SNPs have two levels of meaning: 1.) SNPs have same 
#' rs-names but different versions of SNP ID ound in chip annotation file. 
#' e.g. SNP-A IDs for Affymetrix chip. 2.) SNPs with duplicated genomic 
#' position: the combination of base pair position and chromosomal location. 

#' @export  
#' @author Junfang Chen 
##' @examples 


removedDoubleProbes <- function(plink, inputPrefix, chipAnnoFile, 
                                chipType, outputSNPdupFile, outputPrefix){

     if (!is.null(chipAnnoFile)){ 

        annoFile <- "chipAnnoRefb37.txt"
        if (chipType == "affymetrix") { 
            prepareChipAnnoFile4affymetrix(inputFile=chipAnnoFile, outputFile=annoFile)
        } else if (chipType == "illumina"){ 
            prepareChipAnnoFile4Illumina(inputFile=chipAnnoFile, outputFile=annoFile)
        } else if (chipType == "PsychChip"){ 
            prepareChipAnnoFile4PsychChip(inputFile=chipAnnoFile, outputFile=annoFile)
        }

        ## find the overlapping
        chipAnno <- read.table(file=annoFile, header=TRUE, stringsAsFactors=FALSE)  
        bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
        ## cbind (combine) bim and chip annotation files  
        sharedSNP = intersect(chipAnno[,"chipSnpID"], bim[,2])
        # chipAnnoV1 <- chipAnno[is.element(chipAnno[,"chipSnpID"], sharedSNP,]
        chipAnnoV1sort <- chipAnno[match(sharedSNP, chipAnno[,"chipSnpID"]),]
        bimV1 <- bim[match(sharedSNP, bim[,2]), ]
        comb <- cbind(bimV1, chipAnnoV1sort)   

        ## remove SNPs with duplicated position first
        chrNames <- names(table(comb[,1]))
        dupPos <- c()
        for (i in chrNames) { 
            print(i)
            subData <- comb[which(comb[,1] == i), ]  
            ## Remove the 1st duplicated ID (or 2nd if there are three replicates)
            subDup <- subData[duplicated(subData[,"pos"], fromLast=TRUE), ] 
            # print(dim(subDup))
             dupPos <- rbind(dupPos, subDup)
        }  
        snpWithdupPos <- dupPos[,"chipSnpID"]  
        ## return all the duplicated rs-names (not only either one)

        if (chipType == "affymetrix") { 
            whDup <- duplicated(comb[,"rsID"]) | duplicated(comb[,"rsID"], fromLast=TRUE)
            snpdup <- comb[whDup, "chipSnpID"] 
        } else if (chipType == "illumina"){ 
            whDup <- duplicated(comb[,"V2"]) | duplicated(comb[,"V2"], fromLast=TRUE)
            snpdup <- comb[whDup, "chipSnpID"]  
        } else if (chipType == "PsychChip"){ 
            whDup <- duplicated(comb[,"V2"]) | duplicated(comb[,"V2"], fromLast=TRUE)
            snpdup <- comb[whDup, "chipSnpID"]  
        }

        allDupSNPs <- c(snpWithdupPos, snpdup) 
        allDupSNPs <- unique(allDupSNPs)
        write.table(allDupSNPs, file=outputSNPdupFile, quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 
        system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
               outputSNPdupFile, " --make-bed --out ", outputPrefix))

    } else { 
        ## copy/rename plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy")  
    }
}



##########################################   
##########################################  updateSNPinfo
#' Update the SNP information
#'
#' @description
#' Update SNP information including SNP name, base-pair position, 
#' chromosomal location and the strand information.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param chipAnnoFile a pure text file that stores the chip annotation 
#' information. 
#' If the chip annotation file is not available for your study, it can be 
#' downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.

#' @param chipType a string name defines the type of the chip annotation file: 
#' 'illumina', 'affymetrix' or 'psychChip'.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after updating SNP information.
#' @details The SNP information in the chip annotation file is used as
#' the reference. 

#' @export 
#' @author Junfang Chen 

updatedSnpInfo <- function(plink, inputPrefix, 
                           chipAnnoFile, chipType, outputPrefix){

    if (!is.null(chipAnnoFile)){ 
        annoFile <- "chipAnnoRefb37.txt"
        if (chipType == "affymetrix") { 
            prepareChipAnnoFile4affymetrix(inputFile=chipAnnoFile, 
                                           outputFile=annoFile)
        } else if (chipType == "illumina"){ 
            prepareChipAnnoFile4Illumina(inputFile=chipAnnoFile, 
                                         outputFile=annoFile)
        } else if (chipType == "PsychChip"){ 
            prepareChipAnnoFile4PsychChip(inputFile=chipAnnoFile, 
                                          outputFile=annoFile)
        }

        ## find the overlapping
        chipAnno <- read.table(file=annoFile, 
                               header=TRUE, stringsAsFactors=FALSE)
        system(paste0("rm ", annoFile)) ## 

        # ## cbind (combine) bim and chip annotation files 
        bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
        interSNPs <- intersect(bim[,2], chipAnno[,"chipSnpID"]) 
        bimV1 <- bim[is.element(bim[,2], interSNPs),]
        chipAnnoV1 <- chipAnno[is.element(chipAnno[,"chipSnpID"], interSNPs),]
        annoV1sort <- chipAnnoV1[match(bimV1[,2], chipAnnoV1[,"chipSnpID"]),]
        comV2 <- cbind(bimV1, annoV1sort) 

        ## Update main info  
        if (chipType == "affymetrix") {      
            updateSNP2rs <- subset(comV2, select=c(V2, rsID))
            updateSNPchr <- subset(comV2, select=c(rsID, chr))
            updateSNPpos <- subset(comV2, select=c(rsID, pos)) 
            ## strand 
            updateSNPbackward <- comV2[which(comV2[,"strand"] == "-"), "rsID"] 
        } else if (chipType == "illumina"){ 
            updateSNP2rs <- subset(comV2, select=c(V2, V2))  
            updateSNPchr <- subset(comV2, select=c(V2, chr))
            updateSNPpos <- subset(comV2, select=c(V2, pos)) 
            ## strand
            updateSNPbackward <- comV2[which(comV2[,"strand"] == "-"), "V2"]  
        } else if (chipType == "PsychChip"){ 
            updateSNP2rs <- subset(comV2, select=c(V2, V2))  
            updateSNPchr <- subset(comV2, select=c(V2, chr))
            updateSNPpos <- subset(comV2, select=c(V2, pos)) 
            ## strand 
            updateSNPbackward <- comV2[which(comV2[,"strand"] == "-"), "V2"]  
        }

        inputPrefix.rs <- "1_09.updatedSnp2rs" ## tobeRemoved
        inputPrefix.chr <- "1_09.updatedSnpchr" ## tobeRemoved
        inputPrefix.pos <- "1_09.updatedSnppos" ## tobeRemoved
        inputPrefix.strand <- "1_09.updatedSnpstrand" ## tobeRemoved
        write.table(updateSNP2rs, file=paste0(inputPrefix.rs, ".txt"), 
                    quote=FALSE, row.names=FALSE, col.names=FALSE, 
                    eol="\r\n", sep=" ") 
        write.table(updateSNPchr, file=paste0(inputPrefix.chr, ".txt"), 
                    quote=FALSE, row.names=FALSE, col.names=FALSE, 
                    eol="\r\n", sep=" ") 
        write.table(updateSNPpos, file=paste0(inputPrefix.pos, ".txt"), 
                    quote=FALSE, row.names=FALSE, col.names=FALSE, 
                    eol="\r\n", sep=" ") 
        write.table(updateSNPbackward, file=paste0(inputPrefix.strand, ".txt"), 
                    quote=FALSE, row.names=FALSE, col.names=FALSE, 
                    eol="\r\n", sep=" ") 
        ## update rs, chr, and pos one by one     ## flip to the forward strand
        system(paste0(plink, " --bfile ", inputPrefix, " --update-name ", 
               inputPrefix.rs, ".txt 2 1  --make-bed --out ", inputPrefix.rs))  
        system(paste0(plink, " --bfile ", inputPrefix.rs,  " --update-chr ", 
               inputPrefix.chr, ".txt 2 1 --make-bed --out ", inputPrefix.chr))  
        system(paste0(plink, " --bfile ", inputPrefix.chr, " --update-map ", 
               inputPrefix.pos, ".txt 2 1 --make-bed --out ", inputPrefix.pos))   
        system(paste0(plink, " --bfile ", inputPrefix.pos, " --flip ", 
               inputPrefix.strand, ".txt --make-bed --out ", inputPrefix.strand))  
        ## copy/rename all snp info updated plink files
        renamePlinkBFile(inputPrefix.strand, outputPrefix, action="copy")   
         ## remove all tmp files (rs, chr, pos)
        system(paste0("rm ", inputPrefix.rs, ".* ", inputPrefix.chr, ".*"))   
        system(paste0("rm ", inputPrefix.pos, ".* ", inputPrefix.strand, ".*")) 

    } else { 
        ## copy/rename plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
    }
}


##########################################   
########################################## 
#' Split chromosome X into pseudoautosomal region and non-pseudoautosomal region.
#'
#' @description
#' Split chromosome X into pseudoautosomal region and non-pseudoautosomal region, 
#' if chromosome X data is available.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after splitting chromosome X into 
#' pseudoautosomal region and non-pseudoautosomal region.
#' @details Genomic coordinate system is on genome build hg19.

#' @export 

#' @author Junfang Chen 
##' @examples 


splitXchr <- function(plink, inputPrefix, outputPrefix){

    bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
    chrDist <- table(bim[,1])
    chr23check <- is.element(names(chrDist), 23)
    if (chr23check == TRUE) {
        ## split X chr into PAR(chr25) and non-PAR (chr23)
        system(paste0(plink, " --bfile ", inputPrefix, 
               " --split-x hg19 --make-bed --out ", outputPrefix))
    } else { 
        ## copy/rename plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
    }
}


##########################################   
##########################################  
#' Remove SNPs on the chromosome Y and mitochondrial DNA
#'
#' @description
#' Remove SNPs on the chromosome Y and mitochondrial DNA. 
#  If no such SNPs, the output files remain the same as the input.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after removing SNPs 
#' on the chromosome Y and mitochondrial DNA.
#' @details Note that if chromosome Y and mitochondrial DNA are available,
#' they must be coded as 24 and 26, respectively. 

#' @export 
#' @author Junfang Chen 
# #' @examples 

removedYMtSnp <- function(plink, inputPrefix, outputPrefix){

    bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
    snpChrY <-  bim[which(bim[,1]== 24), 2]   
    snpChrMt <- bim[which(bim[,1] == 26), 2] 
    snpChrYMt <- c(snpChrY, snpChrMt) 
    write.table(snpChrYMt, file=paste0(outputPrefix, ".txt"), quote=FALSE,
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 
    system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
           paste0(outputPrefix, ".txt"), " --make-bed --out ", outputPrefix))  
    system(paste0("rm ", paste0(outputPrefix, ".txt")))
}

##########################################   
########################################## 
##########################################   
########################################## 



##########################################   
##########################################   
#' prepare Affymetrix chip annotation file 
#'
#' @description
#' Prepare Affymetrix chip annotation file into the format of interest.


#' @param inputFile an input pure text file that contains the chip annotation
#' information. 
#' @param outputFile an output pure text file that stores the chip annotation 
#' information in a pre-defined format. 
 
#' @return a pure text file that stores the prepared chip annotation 
#' information in a user-defined format. 
#' @details If the chip annotation file is not available for your study, it can be 
#' downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.

#' @export 

#' @author Junfang Chen 
#' @examples 

prepareChipAnnoFile4affymetrix <- function(inputFile, outputFile){

    inputNew <- paste0(inputFile, "New")
    system(paste0("sed 1d ", inputFile, " > ", inputNew) )
    chipAnnoRefraw <- read.table(file=inputNew, stringsAsFactors=FALSE)  
    colnames(chipAnnoRefraw) <- c("chipSnpID", "rsID", "chr", "pos", "strand")

    ## remove SNPs with strange strand 
    whUnknown <- which(chipAnnoRefraw[,"strand"] == "---") 
    chipAnnoRefraw2 <- chipAnnoRefraw[-whUnknown,]

    ## only see 3 different cases (if 25--> XY)
    whX <- which(chipAnnoRefraw2[,"chr"] == "X")
    whY <- which(chipAnnoRefraw2[,"chr"] == "Y")
    whMT <- which(chipAnnoRefraw2[,"chr"] == "MT")

    chipAnnoRefraw2[whX,"chr"] <- 23 
    chipAnnoRefraw2[whY,"chr"] <- 24
    chipAnnoRefraw2[whMT,"chr"] <- 26
  
    whAFF <- grep("AFFX-SNP", chipAnnoRefraw2[,"chipSnpID"]) 
    chipAnnoRefraw3 <- chipAnnoRefraw2[-whAFF,]
    write.table(chipAnnoRefraw3, file=outputFile, quote=FALSE, 
                row.names=FALSE, col.names=TRUE, eol="\r\n", sep="\t")

}


##########################################   
########################################## 
#' prepare Illumina chip annotation file 
#'
#' @description
#' Prepare Illumina chip annotation file into the format of interest.


#' @param inputFile an input pure text file that contains the chip annotation
#' information. 
#' @param outputFile an output pure text file that stores the chip annotation 
#' information in a pre-defined format. 
 
#' @return a pure text file that stores the prepared chip annotation 
#' information in a user-defined format. 
#' @details If the chip annotation file is not available for your study, it can be 
#' downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.


#' @export  
#' @author Junfang Chen 
# #' @examples 


prepareChipAnnoFile4Illumina <- function(inputFile, outputFile){

    chipAnnoRefraw <- read.table(file=inputFile, stringsAsFactors=FALSE)
    # print(colnames(chipAnnoRefraw))
    ## c("chipSnpID", "chr", "pos", "match2genom", "strand", "allele")
    chipAnnoRefraw <- chipAnnoRefraw[,-c(4, 6)]
    colnames(chipAnnoRefraw) <- c("chipSnpID", "chr", "pos", "strand")

    ## only see 3 different cases (if 25--> XY)
    whX <- which(chipAnnoRefraw[,"chr"] == "X")
    whY <- which(chipAnnoRefraw[,"chr"] == "Y")
    whMT <- which(chipAnnoRefraw[,"chr"] == "MT")

    chipAnnoRefraw[whX,"chr"] <- 23 
    chipAnnoRefraw[whY,"chr"] <- 24
    chipAnnoRefraw[whMT,"chr"] <- 26 

    chipAnnoNew <- chipAnnoRefraw
    write.table(chipAnnoNew, file=outputFile, quote=FALSE, 
                row.names=FALSE, col.names=TRUE, eol="\r\n", sep="\t")

}



##########################################   
########################################## 
#' prepare PsychChip chip annotation file 
#'
#' @description
#' Prepare PsychChip chip annotation file into the format of interest.

#' @param inputFile an input pure text file that contains the chip annotation
#' information. 
#' @param outputFile an output pure text file that stores the chip annotation 
#' information in a pre-defined format. 
 
#' @return a pure text file that stores the prepared chip annotation 
#' information in a user-defined format. 
#' @details If the chip annotation file is not available for your study, it can be 
#' downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.

#' @export 
#' @author Junfang Chen 
##' @examples 

prepareChipAnnoFile4PsychChip <- function(inputFile, outputFile){

    chipAnnoRefraw <- read.table(file=inputFile, stringsAsFactors=FALSE) 
    colnames(chipAnnoRefraw) <- c("chipSnpID", "chr", "pos", 
                                  "match2genom", "strand", "allele") 
    chipAnnoRefraw <- chipAnnoRefraw[,-c(4, 6)]
    colnames(chipAnnoRefraw) <- c("chipSnpID", "chr", "pos", "strand")
 
    ## only see 3 different cases (if 25--> XY)
    whX <- which(chipAnnoRefraw[,"chr"] == "X")
    whY <- which(chipAnnoRefraw[,"chr"] == "Y")
    whMT <- which(chipAnnoRefraw[,"chr"] == "MT") 
    chipAnnoRefraw[whX,"chr"] <- 23 
    chipAnnoRefraw[whY,"chr"] <- 24
    chipAnnoRefraw[whMT,"chr"] <- 26  

    write.table(chipAnnoRefraw, file=outputFile, quote=FALSE, 
                row.names=FALSE, col.names=TRUE, eol="\r\n", sep="\t")

}



##########################################   
##########################################  



#' Update genotype information
#'
#' @description
#' Update genotype data information of the original PLINK binary files involving  
#' subject metadata information remapping and SNP information rearrangement and 
#' conversion according to the annotation file.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param dupSampleIDFile a pure text file that stores the duplicated sample IDs, 
#' each ID per line. If it is null, then duplicate the input PLINK files from 
#' the last step as the output files. 
#' @param ancestrySymbol an indicator that shows the symbol of genetic ancestry. 
#' If it is null, then all samples are selected. 
#' @param excludedProbeIdsFile a pure text file that stores the SNP IDs, 
#' one per line, which need to be removed. If it is null, no SNPs are removed.
#' @param chipAnnoFile a pure text file that stores the chip annotation
#' information. 
#' @param chipType a string name that defines the type of the chip annotation 
#' file: 'illumina', 'affymetrix' or 'psychChip'.

#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after genotype information remapping. 
#' @details The original PLINK files are implicitly processed by the following steps: 
#' 1.) remove duplicated subjects;
#' 2.) update group ID and sex information;
#' 3.) remove not labelled subjects;
#' 4.) remove subjects with wrong ancestry;
#' 5.) remove incorrectly annotated SNPs;
#' 6.) remove SNPs that are not in the annotation file;
#' 7.) remove duplicated SNPs;
#' 8.) update SNP genomic position and strand information;
#' 9.) splitchromosome X into pseudoautosomal region (PAR) and non-PAR;
#' 10.) remove SNPs on the chromosome Y and mitochondrial DNA.
#' The metadata information file and the chip annotation file are used as the 
#' reference for the update.
#' If the chip annotation file is not available for your study, it can be 
#' downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.

#' @export  
#' @author Junfang Chen 
##' @examples 


updateGenoInfo <- function(plink, inputPrefix, metaDataFile, dupSampleIDFile,
                           ancestrySymbol, excludedProbeIdsFile, chipAnnoFile,
                           chipType, outputPrefix){

    ## step 2
    outputPrefix2 <- "1_02_removedExclInst" 
    removeDupID(plink, dupSampleIDFile, inputPrefix, outputPrefix=outputPrefix2)
    # step 3 replace group IDs 
    metaDataFile <- "1_01_metaData.txt" 
    outputPrefix3 <- "1_03_replacedGroupAndSex"
    updateGroupIdAndSex(plink, inputPrefix=outputPrefix2, 
                        metaDataFile, outputPrefix=outputPrefix3)
    # step 4 remove instances without group IDs 
    outputPrefix4 <- "1_04_removedNoGroupId"
    removeNoGroupId(plink, inputPrefix=outputPrefix3, outputPrefix=outputPrefix4)

    ## step5 remove instances with improper ancestry  
    outputPrefix5 <- "1_05_removedWrongAnceInst"
    removedWrongAnceInst(plink, inputPrefix=outputPrefix4, metaDataFile,  
                         ancestrySymbol, outputPrefix=outputPrefix5)
    ## step 6   
    outputPrefix6 <- "1_06_removedExclProbe"  
    removedExclProbe(plink, inputPrefix=outputPrefix5, 
                     excludedProbeIdsFile, outputPrefix=outputPrefix6) 
    ## step 7 (Optional, if chip annotation file is not given) 
    outputPrefix7 <- "1_07_removedUnmapProbes"   
    outputSNPfile7 <- "1_07_probesUnmapped2ChipRef.txt"
    removedUnmapProbes(plink, inputPrefix=outputPrefix6, chipAnnoFile,
                       outputPrefix=outputPrefix7, outputSNPfile=outputSNPfile7)
    ## step 8 (Optional, if chip annotation file is not given) 
    outputSNPdupFile8 <- "1_08_probesDouble.txt"
    outputPrefix8 <- "1_08_removedDoubleProbes"   
    removedDoubleProbes(plink, inputPrefix=outputPrefix7, chipAnnoFile, 
                        chipType, outputSNPdupFile=outputSNPdupFile8, 
                        outputPrefix=outputPrefix8) 
    ## step 9 (Optional, if chip annotation file is not given) 
    outputPrefix9 <- "1_09_updatedSnpInfo"   
    updatedSnpInfo(plink, inputPrefix=outputPrefix8, 
                   chipAnnoFile, chipType, outputPrefix=outputPrefix9)
    ## step 10      
    outputPrefix10 <- "1_10_splitXchr"
    splitXchr(plink, inputPrefix=outputPrefix9, outputPrefix=outputPrefix10)
    ## step 11  
    outputPrefix11 <- "1_11_removedYMtSnp"
    removedYMtSnp(plink, inputPrefix=outputPrefix10, outputPrefix=outputPrefix11)
    ## remove intermediate files
    system(paste0("rm ", outputPrefix2, ".*"))
    system(paste0("rm ", outputPrefix3, ".*"))
    system(paste0("rm ", outputPrefix4, ".*"))
    system(paste0("rm ", outputPrefix5, ".*"))
    system(paste0("rm ", outputPrefix6, ".*"))
    system(paste0("rm ", outputPrefix7, ".*"))
    system(paste0("rm ", outputPrefix8, ".*"))
    system(paste0("rm ", outputPrefix9, ".*"))
    system(paste0("rm ", outputPrefix10, ".*"))

    if (file.exists(outputSNPfile7)) {  
        system(paste0("rm ", outputSNPfile7))
    }

    if (file.exists(outputSNPdupFile8)) {  
        system(paste0("rm ", outputSNPdupFile8))
    }    
 
}
