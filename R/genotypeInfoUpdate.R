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
#' Remove samples in PLINK files
#'
#' @description
#' Remove sample IDs that are useless such as duplicated or related IDs from 
#' PLINK binary files. These IDs are defined in a plain text file. 

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param removedSampIDFile a pure text file that stores the useless 
#' sample IDs, each ID per line. If it is null, then duplicate the input 
#' PLINK files as the output files.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after removing certain sample IDs. 
#' @export 

#' @author Junfang Chen 
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix
#' removedSampIDFile <- system.file("extdata", "excludedSampIDs.txt", 
#'                                  package="Gimpute")
#' outputPrefix <- "1_02_removedExclInst" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removeSampID(plink, removedSampIDFile, inputPrefix, outputPrefix)

removeSampID <- function(plink, removedSampIDFile, inputPrefix, outputPrefix){

    if (!is.null(removedSampIDFile)){ 

        famv0 <- read.table(file=paste0(inputPrefix, ".fam"), 
                            stringsAsFactors=FALSE)  
        rmvIDs <- read.table(file=removedSampIDFile, stringsAsFactors=FALSE)  
        ## write into plink format .txt for removal 
        rmvIDmat <- famv0[is.element(famv0[,2], rmvIDs[,1]), c("V1", "V2")] 
        rmvIDfn <- paste0(outputPrefix, ".txt")
        write.table(rmvIDmat, file=rmvIDfn, quote=FALSE, row.names=FALSE, 
                    col.names=FALSE, eol="\r\n", sep=" ") 
        system(paste0(plink, " --bfile ", inputPrefix, 
               " --remove ", rmvIDfn, " --make-bed --out ", outputPrefix))
    } else { 
        ## copy/rename all snp info updated plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy")
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
#' the samples. This file must contain at least the following content 
#' (column names are in parentheses):
#' family ID in the PLINK files (FID), individual ID in the PLINK files (IID), 
#' ID in the description files (descID), self identified ancestry 
#' (ance; e.g. AFR: African, AMR: Ad Mixed American, EAS: East Asian, 
#' EUR: European, SAS: South Asian), sex (sex; 1 = male, 2 = female), 
#' age (age), group (group; 0 = control/unaffected, 1 = case/affected). 
#' All unknown and missing values are represented by the value NA. 
#' Lines with a missing value for FID or IID are not contained.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after updating the gender and grouping 
#' information. 

#' @details Find the shared sample IDs between PLINK input files and 
#' metadata file. Use the information from the metadata file as the reference 
#' and update the group information/the outcome in the PLINK file. 
#' Group label should be 1 and 2. (1=unaff, 2=aff, 0=miss); 
#' missing phenotype will be indicated as -9.

#' @export  
#' @author Junfang Chen 

#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' metaDataFile <- system.file("extdata", "1_01_metaData.txt",
#'                             package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix
#' outputPrefix <- "1_03_replacedGroupAndSex" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## updateGroupIdAndSex(plink, inputPrefix, metaDataFile, outputPrefix)


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
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after removing samples without group 
#' information.
#' @export 

#' @author Junfang Chen 

#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix
#' outputPrefix <- "1_04_removedNoGroupId" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removeNoGroupId(plink, inputPrefix, outputPrefix)



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
#' Remove samples with the incorrect ancestry or keep samples 
#' at your own chioce.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param metaDataFile a pure text file that stores the meta information of 
#' the samples. This file must contain at least the following content 
#' (column names are in parentheses):
#' family ID in the PLINK files (FID), individual ID in the PLINK files (IID), 
#' ID in the description files (descID), self identified ancestry 
#' (ance; e.g. AFR: African, AMR: Ad Mixed American, EAS: East Asian, 
#' EUR: European, SAS: South Asian), sex (sex; 1 = male, 2 = female), 
#' age (age), group (group; 0 = control/unaffected, 1 = case/affected). 
#' All unknown and missing values are represented by the value NA. 
#' Lines with a missing value for FID or IID are not contained.
#' @param ancestrySymbol an indicator that shows the symbol of genetic ancestry. 
#' If it is null, then all samples are selected. 
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after checking 
#' the ancestry information. 
#' @details ancestrySymbol, such as 'EUR' stands for the European, 
#' 'EAS' for East Asian. See the metaDataFile for more details. 
#' @export  
#' @author Junfang Chen  

#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' metaDataFile <- system.file("extdata", "1_01_metaData.txt", 
#'                             package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix
#' ancestrySymbol <- "EAS"
#' outputPrefix <- "1_05_removedWrongAnceInst" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedWrongAnceInst(plink, inputPrefix, metaDataFile,  
#' ##                      ancestrySymbol, outputPrefix)


removedWrongAnceInst <- function(plink, inputPrefix, metaDataFile, 
                                  ancestrySymbol, outputPrefix){

    if (is.null(ancestrySymbol)) { 
        ## copy/rename all snp info updated plink files
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 
    } else {
        metaData <- read.table(metaDataFile, 
                               stringsAsFactors=FALSE, header=TRUE) 
        ids <- metaData[which(metaData[,"ance"] == ancestrySymbol),"IID"]
        fam <- read.table(file=paste0(inputPrefix, ".fam"), 
                          stringsAsFactors=FALSE) 
        famEA <- fam[is.element(fam[,2], ids), c("V1", "V2")]
        plinkFormatDat <- famEA
        write.table(plinkFormatDat, file=paste0(outputPrefix,".txt"), 
                    quote=FALSE, row.names=FALSE, 
                    col.names=FALSE, eol="\r\n", sep=" ")  

        system(paste0(plink, " --bfile ", inputPrefix, " --keep ", 
               outputPrefix, ".txt --make-bed --out ", outputPrefix))
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
#' one per line, which need to be removed. If it is null, 
#' then no SNPs are removed.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after removing unwanted SNP IDs.
#' @details excludedProbeIdsFile should be defined in a plain text file 
#' in advance. Improper SNPs such as AFFX, cnvi etc.. with unexpect format
#' must be excluded.

#' @export 
#' @author Junfang Chen 
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' excludedProbeIdsFile <- system.file("extdata", "excludedProbeIDs.txt",
#'                                     package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix
#' outputPrefix <- "1_06_removedExclProbe" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedExclProbe(plink, inputPrefix, excludedProbeIdsFile, outputPrefix)

removedExclProbe <- function(plink, inputPrefix, 
                             excludedProbeIdsFile, outputPrefix){

    if (!is.null(excludedProbeIdsFile)) {
        system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
               excludedProbeIdsFile, " --make-bed --out ", outputPrefix) ) 
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
#' @param chipType a string name defines the type of the chip annotation file: 
#' 'SNPIDstudy', and 'rsIDstudy'. The detail is described in 
#' \code{\link{prepareAnnoFile4affy}}.
#' @param outputSNPfile a pure text file that stores the SNP IDs, 
#' one per line, which are not mapped to the chip annotation file.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output text file contains the removed SNP IDs, one per line.
#' The PLINK binary files after removing unmapped SNP IDs.
#' @details If the chip annotation file is not available for your study, 
#' you can download it from http://www.well.ox.ac.uk/~wrayner/strand/.

#' @export 
#' @author Junfang Chen 
#' @seealso \code{\link{prepareAnnoFile4affy}}
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' chipAnnoFile <- system.file("extdata", "chipAnno.txt", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix
#' outputSNPfile <- "1_07_removedUnmapProbes"   
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedUnmapProbes(plink, inputPrefix, chipAnnoFile, chipType, 
#' ##                    outputPrefix, outputSNPfile)

removedUnmapProbes <- function(plink, inputPrefix, chipAnnoFile, chipType, 
                               outputPrefix, outputSNPfile){

    if (!is.null(chipAnnoFile)){  

        chipAnno <- read.table(file=chipAnnoFile, header=TRUE, 
                               stringsAsFactors=FALSE)
        ## check the overlapping SNPs
        bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
        if (chipType == "SNPIDstudy"){
            interSNPs <- intersect(bim[,2], chipAnno[,"SNPIDstudy"])  
        } else if (chipType == "rsIDstudy"){
            interSNPs <- intersect(bim[,2], chipAnno[,"rsIDstudy"])   
        } else {
            print("Wrong chipType or column names in the annotation file!")
        }  
        unmapped <- setdiff(bim[,2], interSNPs) 
        write.table(unmapped, file=outputSNPfile, quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
        system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
               outputSNPfile, " --make-bed --out ", outputPrefix))
    } else {  
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy")  
    }
}



##########################################   
########################################## 
#' Remove duplicated SNPs
#'
#' @description
#' Remove duplicated SNPs that have same rs-names or duplicated genomic  
#' position.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param chipAnnoFile an input chip annotation file. 
#' If the chip annotation file is not available for your study, it can be 
#' downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.
#' @param chipType a string name defines the type of the chip annotation file: 
#' 'SNPIDstudy', and 'rsIDstudy'. The detail is described in 
#' \code{\link{prepareAnnoFile4affy}}.
#' @param outputSNPdupFile a pure text file that stores the duplicated SNP IDs, 
#' which are detected by the use of the chip annotation file.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after removing duplicated SNP IDs.
#' @details Duplicated SNPs have two levels of meaning: 1.) SNPs have same 
#' rs-names but different versions of SNP ID in chip annotation file. 
#' e.g. SNP-A IDs for Affymetrix chip. 2.) SNPs with duplicated genomic 
#' position: the combination of base pair position and chromosomal location. 

#' @export  
#' @author Junfang Chen 
#' @seealso \code{\link{prepareAnnoFile4affy}}
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' chipAnnoFile <- system.file("extdata", "chipAnno.txt", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))
#' inputPrefix <- "controlData"     
#' chipType <- "rsIDstudy"
#' outputSNPdupFile <- "snpDup.txt"
#' outputPrefix <- "removedDoubleProbes"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedDoubleProbes(plink, inputPrefix, chipAnnoFile,
#' ##                     chipType, outputSNPdupFile, outputPrefix)



removedDoubleProbes <- function(plink, inputPrefix, chipAnnoFile, 
                                chipType, outputSNPdupFile, outputPrefix){

     if (!is.null(chipAnnoFile)){ 

        chipAnno <- read.table(file=chipAnnoFile, header=TRUE, 
                               stringsAsFactors=FALSE)  
        bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)
        if (chipType == "SNPIDstudy"){
            interSNPs <- intersect(bim[,2], chipAnno[,"SNPIDstudy"])  
            chipAnnoV1 <- chipAnno[match(interSNPs, chipAnno[,"SNPIDstudy"]),]
        } else if (chipType == "rsIDstudy"){
            interSNPs <- intersect(bim[,2], chipAnno[,"rsIDstudy"])    
            chipAnnoV1 <- chipAnno[match(interSNPs, chipAnno[,"rsIDstudy"]),]
        } else {
            print("Wrong chipType or column names in the annotation file!")
        }  

        bimV1 <- bim[match(interSNPs, bim[,2]), ]
        comb <- cbind(bimV1[,seq_len(4)], chipAnnoV1)    
        ## remove SNPs with duplicated pos 
        chrNames <- names(table(comb[,"chr"]))
        dupPos <- c()
        for (i in chrNames) { 
            print(i)
            chrData <- comb[which(comb[,"chr"] == i), ]  
            ## Remove 1st duplicated ID (or 2nd if 3 replicates)
            subDup <- chrData[duplicated(chrData[,"pos"], fromLast=TRUE), ]  
            dupPos <- rbind(dupPos, subDup)
        }  
        ## return all the duplicated rs-names (not only either one) 
        if (chipType == "SNPIDstudy") { 
            snpWithdupPos <- dupPos[,"SNPIDstudy"]
            whDup <- duplicated(comb[,"rs"]) | 
                     duplicated(comb[,"rs"], fromLast=TRUE)
            snpdup <- comb[whDup, "SNPIDstudy"] 
        } else if (chipType == "rsIDstudy"){ 
            snpWithdupPos <- dupPos[,"rsIDstudy"]
            whDup <- duplicated(comb[,"V2"]) | 
                     duplicated(comb[,"V2"], fromLast=TRUE)
            snpdup <- comb[whDup, "rsIDstudy"]  
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
#' 'SNPIDstudy', and 'rsIDstudy'. The detail is described in 
#' \code{\link{prepareAnnoFile4affy}}.

#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after updating SNP information.
#' @details The SNP information in the chip annotation file is used as
#' the reference. 

#' @export 
#' @author Junfang Chen  
#' @seealso \code{\link{prepareAnnoFile4affy}}
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' chipAnnoFile <- system.file("extdata", "chipAnno.txt", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))
#' inputPrefix <- "controlData"     
#' chipType <- "rsIDstudy"
#' outputPrefix <- "updatedSnpInfo"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## updatedSnpInfo(plink, inputPrefix, chipAnnoFile, 
#' ##                chipType, outputPrefix, outputSNPfile)


updatedSnpInfo <- function(plink, inputPrefix, 
                           chipAnnoFile, chipType, outputPrefix){

    if (!is.null(chipAnnoFile)){   

        chipAnno <- read.table(file=chipAnnoFile, 
                               header=TRUE, stringsAsFactors=FALSE)  
        bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)

        if (chipType == "SNPIDstudy"){
            interSNPs <- intersect(bim[,2], chipAnno[,"SNPIDstudy"])  
            chipAnnoV1 <- chipAnno[match(interSNPs, chipAnno[,"SNPIDstudy"]),]
        } else if (chipType == "rsIDstudy"){
            interSNPs <- intersect(bim[,2], chipAnno[,"rsIDstudy"])    
            chipAnnoV1 <- chipAnno[match(interSNPs, chipAnno[,"rsIDstudy"]),]
        } else {
            print("Wrong chipType or column names in the annotation file!")
        }  

        bimV1 <- bim[match(interSNPs, bim[,2]), ]
        comV2 <- cbind(bimV1[,], chipAnnoV1)     
        ## Update geno info  
        if (chipType == "SNPIDstudy") {      
            updateSNP2rs <- subset(comV2, select=c(V2, rs))
            updateSNPchr <- subset(comV2, select=c(rs, chr))
            updateSNPpos <- subset(comV2, select=c(rs, pos)) 
            ## strand 
            updateSNPbackward <- comV2[which(comV2[,"strand"] == "-"), "rs"] 
        } else if (chipType == "rsIDstudy"){ 
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
        ## update rs, chr, and pos one by one  ## flip to the forward strand
        ## --update-name [filename] {new ID col. number} {old ID col.} {skip}
        system(paste0(plink, " --bfile ", inputPrefix, " --update-name ", 
               inputPrefix.rs, ".txt 2 1  --make-bed --out ", inputPrefix.rs))  
        system(paste0(plink, " --bfile ", inputPrefix.rs,  " --update-chr ", 
               inputPrefix.chr, ".txt 2 1 --make-bed --out ", inputPrefix.chr))  
        system(paste0(plink, " --bfile ", inputPrefix.chr, " --update-map ", 
               inputPrefix.pos, ".txt 2 1 --make-bed --out ", inputPrefix.pos))   
        system(paste0(plink, " --bfile ", inputPrefix.pos, " --flip ", 
               inputPrefix.strand, ".txt --make-bed --out ", 
               inputPrefix.strand))  
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
#' Split chromosome X into pseudoautosomal region and 
#' non-pseudoautosomal region.
#' @description
#' Split chromosome X into pseudoautosomal region and 
#' non-pseudoautosomal region, if chromosome X data is available.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after splitting chromosome X into 
#' pseudoautosomal region and non-pseudoautosomal region.
#' @details Genomic coordinate system is on genome build hg19.

#' @export  
#' @author Junfang Chen 
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix 
#' outputPrefix <- "1_10_splitXchr" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removeSampID(plink, inputPrefix, outputPrefix)

splitXchr <- function(plink, inputPrefix, outputPrefix){

    bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)    
    chrCodes <- names(table(bim[,1]))
    chr23check <- length(grep(23, chrCodes))
    if (chr23check == 1) { 
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
#' @examples 
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix 
#' outputPrefix <- "1_11_removedYMtSnp" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedYMtSnp(plink, inputPrefix, outputPrefix)

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
#' prepare Affymetrix chip annotation file 
#'
#' @description
#' Prepare Affymetrix chip annotation file into the format of interest.


#' @param inputFile an input pure text file that contains the chip annotation
#' information. 
#' @param outputFile an output pure text file that stores the chip annotation 
#' information in a user-defined format. 
#' @param chipType a string name defines the type of the chip annotation file: 
#' 'SNPIDstudy', and 'rsIDstudy'.

#' @return a pure text file that stores the prepared chip annotation 
#' information in a user-defined format. 
#' @details If the chip annotation file is not available for your study, 
#' it can be downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.
#' The chip annotation file is organized into two different types: 
#' \enumerate{
#'   \item If the snp name of your study genotype data starts with "SNP_", 
#'         then the chip type "SNPIDstudy" is used; Usually, Affymetrix chip 
#'         data belongs to this category. The prepared output annotation file  
#'         must at least consist of the following column names: 
#'         SNPIDstudy, rs, chr, pos, strand.
#'   \item If the snp name of your study genotype data starts with "rs", then
#'         the chip type "rsIDstudy" is used; The prepared output annotation 
#'         file must at least consist of the following column names:  
#'         SNPIDstudy, rs, chr, pos, strand. Illumina chip is often specified 
#'         in this format.
#' }
#' The column "strand" must only have two kinds of values "-" and "+". 
#' Variants with all other values should be excluded.

#' @export 

#' @author Junfang Chen 
##' @examples 

prepareAnnoFile4affy <- function(inputFile, outputFile, chipType){

    inputNew <- paste0(inputFile, "NewFile")
    system(paste0("sed 1d ", inputFile, " > ", inputNew)) ## remove 1st line
    chipAnnoRaw <- read.table(file=inputNew, stringsAsFactors=FALSE)  

    if (chipType == "SNPIDstudy"){
        colnames(chipAnnoRaw) <- c("SNPIDstudy", "rs", "chr", "pos", "strand")
    } else if (chipType == "rsIDstudy"){
        chipAnnoRaw <- chipAnnoRaw[,-1] ## remove SNP_A ID
        colnames(chipAnnoRaw) <- c("rsIDstudy", "chr", "pos", "strand")
    } else {
        print("Wrong chipType or column names in the annotation file!")
    }  

    ## remove SNPs with strange strand 
    whUnknown <- which(chipAnnoRaw[,"strand"] == "---") 
    chipAnnoRawV2 <- chipAnnoRaw[-whUnknown,]

    ## only see 3 different cases (if 25--> XY)
    whX <- which(chipAnnoRawV2[,"chr"] == "X")
    whY <- which(chipAnnoRawV2[,"chr"] == "Y")
    whMT <- which(chipAnnoRawV2[,"chr"] == "MT")

    chipAnnoRawV2[whX,"chr"] <- 23 
    chipAnnoRawV2[whY,"chr"] <- 24
    chipAnnoRawV2[whMT,"chr"] <- 26

    write.table(chipAnnoRawV2, file=outputFile, quote=FALSE, 
                row.names=FALSE, col.names=TRUE, eol="\r\n", sep="\t")
}



##########################################   
##########################################  



#' Update genotype information
#'
#' @description
#' Update genotype information of the original PLINK binary files involving  
#' subject metadata information remapping and SNP information rearrangement and 
#' conversion according to the annotation file.

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param metaDataFile a pure text file that stores the meta information of 
#' the samples. This file must contain at least the following content 
#' (column names are in parentheses):
#' family ID in the PLINK files (FID), individual ID in the PLINK files (IID), 
#' ID in the description files (descID), self identified ancestry 
#' (ance; e.g. AFR: African, AMR: Ad Mixed American, EAS: East Asian, 
#' EUR: European, SAS: South Asian), sex (sex; 1 = male, 2 = female), 
#' age (age), group (group; 0 = control/unaffected, 1 = case/affected). 
#' All unknown and missing values are represented by the value NA. 
#' Lines with a missing value for FID or IID are not contained.

#' @param removedSampIDFile a pure text file that stores the useless sample IDs, 
#' each ID per line. If it is null, then duplicate the input PLINK files from 
#' the last step as the output files. 
#' @param ancestrySymbol an indicator that shows the symbol of genetic ancestry. 
#' If it is null, then all samples are selected. 
#' @param excludedProbeIdsFile a pure text file that stores the SNP IDs, 
#' one per line, which need to be removed. If it is null, no SNPs are removed.
#' @param chipAnnoFile a pure text file that stores the chip annotation
#' information. 
#' @param chipType a string name defines the type of the chip annotation file: 
#' 'SNPIDstudy', and 'rsIDstudy'. The detail is described in 
#' \code{\link{prepareAnnoFile4affy}}.

#' @param keepInterFile a logical value indicating if the intermediate 
#' processed files should be kept or not. The default is TRUE.
#' @param outputPrefix the prefix of the output PLINK binary files.
 
#' @return The output PLINK binary files after genotype information remapping. 
#' @details The original PLINK files are implicitly processed by the 
#' following steps: 
#' 1.) remove duplicated subjects;
#' 2.) update group ID and sex information;
#' 3.) remove not labelled subjects;
#' 4.) remove subjects with wrong ancestry;
#' 5.) remove incorrectly annotated SNPs;
#' 6.) remove SNPs that are not in the annotation file;
#' 7.) remove duplicated SNPs;
#' 8.) update SNP genomic position and strand information;
#' 9.) split chromosome X into pseudoautosomal region (PAR) and non-PAR;
#' 10.) remove SNPs on the chromosome Y and mitochondrial DNA.
#' The metadata information file and the chip annotation file are used as the 
#' reference for the update.
#' If the chip annotation file is not available for your study, it can be 
#' downloaded from http://www.well.ox.ac.uk/~wrayner/strand/.

#' @export  
#' @author Junfang Chen   
#' @references Purcell, Shaun, et al. PLINK: a tool set for whole-genome 
#' association and population-based linkage analyses. The American Journal of 
#' Human Genetics 81.3 (2007): 559-575.

#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" 
#' metaDataFile <- system.file("extdata", "1_01_metaData.txt",  
#'                             package="Gimpute")
#' excludedProbeIdsFile <- system.file("extdata", "excludedProbeIDs.txt", 
#'                                     package="Gimpute")
#' removedSampIDFile <- system.file("extdata", "excludedSampIDs.txt", 
#'                                  package="Gimpute")
#' chipAnnoFile <- system.file("extdata", "chipAnno.txt", package="Gimpute")
#' ancestrySymbol <- "EUR"
#' outputPrefix <- "1_11_removedYMtSnp" 
#' metaDataFile <- "1_01_metaData.txt"
#' chipType <- "rsIDstudy"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"  
#' ## updateGenoInfo(plink, inputPrefix, metaDataFile, removedSampIDFile,
#' ##                ancestrySymbol, excludedProbeIdsFile, chipAnnoFile,
#' ##                chipType, outputPrefix, keepInterFile=TRUE)
 


updateGenoInfo <- function(plink, inputPrefix, metaDataFile, removedSampIDFile,
                           ancestrySymbol, excludedProbeIdsFile, chipAnnoFile,
                           chipType, outputPrefix, keepInterFile=TRUE){

    ## step 2
    outputPrefix2 <- "1_02_removedExclInst" 
    removeSampID(plink, removedSampIDFile, inputPrefix, 
                 outputPrefix=outputPrefix2)
    # step 3 replace group IDs 
    metaDataFile <- "1_01_metaData.txt" 
    outputPrefix3 <- "1_03_replacedGroupAndSex"
    updateGroupIdAndSex(plink, inputPrefix=outputPrefix2, 
                        metaDataFile, outputPrefix=outputPrefix3)
    # step 4 remove instances without group IDs 
    outputPrefix4 <- "1_04_removedNoGroupId"
    removeNoGroupId(plink, inputPrefix=outputPrefix3, 
                    outputPrefix=outputPrefix4)

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
    removedUnmapProbes(plink, inputPrefix=outputPrefix6, chipAnnoFile, chipType,
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
    removedYMtSnp(plink, inputPrefix=outputPrefix10, 
                  outputPrefix=outputPrefix11)

    if (keepInterFile==FALSE){ 
        # remove intermediate files
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
    ## remove unwanted files
    system(paste0("rm  *.log "))     
}
