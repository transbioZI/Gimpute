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
#' In principle, one has to remove all heterozygous SNPs of chromosome X 
#' in males. 
#' However, too many SNPs might be removed in some data sets. 
#' Therefore a small percentage of such SNPs in the data set is allowed.

#' @export  
#' @author Junfang Chen  
#' @examples 
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" 
#' hhCutOff <- 0.005 ##  can be tuned
#' outputPrefix <- "2_01_removedSnpHetX" 
#' outputHetSNPfile <- "2_01_snpHHfreqAll.txt"
#' outputRetainSNPfile <- "2_01_snpHHfreqRetained.txt"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedSnpHetX(plink, inputPrefix, hhCutOff, outputPrefix, 
#' ##                outputHetSNPfile, outputRetainSNPfile)


removedSnpHetX <- function(plink, inputPrefix, hhCutOff, outputPrefix, 
                           outputHetSNPfile, outputRetainSNPfile){

    bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
    chrCodes <- names(table(bim[,1]))
    chr23check <- length(grep(23, chrCodes))

    if (chr23check == 1) {  
        ## just to get .hh file and .fam file 
        system(paste0(plink, " --bfile ", inputPrefix, 
               " --chr 23 --filter-males --make-bed --out male23nonPAR"))

        if (is.na(file.size("male23nonPAR.hh"))) { 
            ## copy/rename all snp info updated plink files
            renamePlinkBFile(inputPrefix, outputPrefix, action="copy")  
        } else {      

            ## *.hh: A text file with one line per error (sorted primarily by 
            ## variant ID, 2nd by sample ID) with the following 3 fields:
            # Family ID  Within-family ID Variant ID
            hh <- read.table("male23nonPAR.hh", stringsAsFactors=FALSE)
            fam <- read.table("male23nonPAR.fam", stringsAsFactors=FALSE)

            hetSNPsFreq <- table(hh[,3])
            # hetSNPFreqFreq <- table(hetSNPs) 

            cut4removeHetSNP <- nrow(fam) * hhCutOff
            mostFakeSNPs <- hetSNPsFreq[which(hetSNPsFreq >= cut4removeHetSNP)] 
            mostFakeSNPs <- names(mostFakeSNPs)
            str(mostFakeSNPs)
            write.table(mostFakeSNPs, file="mostFakeSNPs.txt", quote=FALSE, 
                        row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 
            ## remove these fake SNPs
            system(paste0(plink, " --bfile ", inputPrefix, 
                   " --exclude mostFakeSNPs.txt --make-bed --out ", 
                   outputPrefix))
            ## remove unwanted files
            system(paste0("rm mostFakeSNPs.txt"))

            ## generate hetSNPsFreq in .txt file 
            hetSNPinstNum <- as.data.frame(hetSNPsFreq, stringsAsFactors=FALSE)
            hetSNPinstNum <- hetSNPinstNum[order(hetSNPinstNum[,2], 
                                                 decreasing=TRUE),] 
            ## all heterozygous SNPs 
            write.table(hetSNPinstNum, file=outputHetSNPfile, quote=FALSE, 
                        row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 

            ## remaining heterozygous SNPs 
            hetSNPinstNumSub <- hetSNPinstNum[!is.element(hetSNPinstNum[,1], 
                                                          mostFakeSNPs), ]
            write.table(hetSNPinstNumSub, file=outputRetainSNPfile, 
                        quote=FALSE, row.names=FALSE, 
                        col.names=FALSE, eol="\r\n", sep=" ") 
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
#' that have heterozygous SNPs with their frequency (if any), i.e. the number 
#' of .hh SNPs in this male. Lines are sorted by descending number.
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
#' @examples 
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" 
#' hhSubjCutOff <- 15 ##  can be tuned
#' outputPrefix <- "2_02_removedInstHetX" 
#' outputSubjHetFile <- "2_02_instHetXfreqAll.txt" 
#' outputRetainSubjectFile <- "2_02_instHetXfreqRetained.txt"  
#' outputHetSNPfile <- "2_02_snpHHfreqAll.txt"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedMaleHetX(plink, inputPrefix, hhSubjCutOff,
#' ##                 outputPrefix, outputSubjHetFile, 
#' ##                 outputRetainSubjectFile, outputHetSNPfile)


removedMaleHetX <- function(plink, inputPrefix, hhSubjCutOff, outputPrefix, 
                            outputSubjHetFile, outputRetainSubjectFile, 
                            outputHetSNPfile){        

    bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
    chrCodes <- names(table(bim[,1]))
    chr23check <- length(grep(23, chrCodes))
    if (chr23check == 1) { 
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
            ## variant ID,secondarily by sample ID) with the following 3 fields:
            ## Family ID  Within-family ID Variant ID
            hh <- read.table("male23nonPAR.hh", stringsAsFactors=FALSE)
            fam <- read.table("male23nonPAR.fam", stringsAsFactors=FALSE)
                            
            hetInstFreq <- table(hh[,2])  
            str(unique(hh[,2]))
            mostFakeInst <- hetInstFreq[which(hetInstFreq >= hhSubjCutOff)]
            whFakeID <- is.element(fam[,2], names(mostFakeInst))  
            mostFakeInstID <- fam[whFakeID, c("V1", "V2")]

            write.table(mostFakeInstID, file="mostFakeInst4plink.txt", 
                        quote=FALSE, row.names=FALSE, 
                        col.names=FALSE, eol="\r\n", sep=" ") 
            ## remove these fake SNPs
            system(paste0(plink, " --bfile ", inputPrefix, 
                   " --remove mostFakeInst4plink.txt --make-bed --out ", 
                   outputPrefix))
            system(paste0("rm mostFakeInst4plink.txt"))

            ## generate hetSNPsFreq in .txt file 
            InstHetSNP <- as.data.frame(hetInstFreq, stringsAsFactors=FALSE) 
            InstHetSNP <- InstHetSNP[order(InstHetSNP[,2], decreasing=TRUE),] 
            write.table(InstHetSNP, file=outputSubjHetFile, quote=FALSE, 
                        row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 

            ## remaining males with heterozygous SNPs 
            InstHetSNPsub <- InstHetSNP[!is.element(InstHetSNP[,1], 
                                                    names(mostFakeInst)), ] 
            write.table(InstHetSNPsub, file=outputRetainSubjectFile, 
                        quote=FALSE, row.names=FALSE, 
                        col.names=FALSE, eol="\r\n", sep=" ") 
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
            hetSNPinstNum <- hetSNPinstNum[order(hetSNPinstNum[,2], 
                                                 decreasing=TRUE),] 
            ## all heterozygous SNPs 
            write.table(hetSNPinstNum, file=outputHetSNPfile, 
                        quote=FALSE, row.names=FALSE, 
                        col.names=FALSE, eol="\r\n", sep=" ")   
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
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" 
#' outputPrefix <- "2_03_setHeteroHaploMissing" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## setHeteroHaploMissing(plink, inputPrefix, outputPrefix)

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
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' snpMissCutOff <- 0.05
#' inputPrefix <- "controlData" 
#' outputPrefix <- "2_04_removedSnpMissPre" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedSnpMiss(plink, snpMissCutOff, inputPrefix, outputPrefix)
 
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
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' sampleMissCutOff <- 0.02
#' inputPrefix <- "controlData" 
#' outputPrefix <- "2_05_removedInstMiss" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedInstMiss(plink, sampleMissCutOff, inputPrefix, outputPrefix)
 

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
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' Fhet <- 0.2
#' inputPrefix <- "controlData" 
#' outputPrefix <- "2_06_removedInstFhet" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedInstFhet(plink, Fhet, inputPrefix, outputPrefix)
 

removedInstFhet <- function(plink, Fhet, inputPrefix, outputPrefix){ 

    system(paste0(plink, " --bfile ", inputPrefix, 
           " --het --out ", outputPrefix))
    ##  F inbreeding coefficient estimate
    autoHet <- read.table(file=paste0(outputPrefix, ".het"), header=TRUE)  
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
#' @examples
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" 
#' outputPrefix <- "2_07_removedParentIdsMiss" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedParentIdsMiss(plink, inputPrefix, outputPrefix)

 
removedParentIdsMiss <- function(plink, inputPrefix, outputPrefix){ 

    # Remove the parent IDs which do not belong to subjects
    system(paste0(plink, " --bfile ", inputPrefix, 
           " --make-founders require-2-missing --make-bed --out ", 
           outputPrefix)) 
 
}

##########################################   
##########################################  
#' Remove SNPs with difference in SNP missingness between cases and controls. 
#'
#' @description
#' Remove SNPs with difference in SNP missingness between cases and controls. 
#' To test for differential call rates between cases and controls for each SNP
 
#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param snpMissDifCutOff the cutoff of the difference in missingness between 
#' cases and controls. 
#' @param outputPrefix the prefix of the output PLINK binary files.
#' @param groupLabel a string value indicating the outcome label: "control", or, 
#' "case" or "caseControl" for both existing groups.
 
#' @return The output PLINK binary files.
#' @details Only if both case-control groups exist in the input genotype data, 
#' differential SNPs are removed. 

#' @export 
#' @author Junfang Chen 
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" 
#' snpMissDifCutOff <- 0.02
#' outputPrefix <- "2_09_removedSnpMissDiff" 
#' groupLabel <- "control"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedSnpMissDiff(plink, inputPrefix, snpMissDifCutOff, 
#' ##                    outputPrefix, groupLabel)

removedSnpMissDiff <- function(plink, inputPrefix, snpMissDifCutOff, 
                               outputPrefix, groupLabel){


    if (groupLabel != "caseControl"){   
        ## this is only for the control data set
        renamePlinkBFile(inputPrefix, outputPrefix, action="copy") 

    } else if (groupLabel == "caseControl") { 
        outputPrefix.tmp <- paste0(outputPrefix, "tmp")
        system (paste0(plink, " --bfile ", inputPrefix, 
                " --test-missing --out ", outputPrefix.tmp) ) 
        ## Write case/control missingness test to [ *.missing ]  
        ## compute differential call rates 
        ## between cases and controls for each SNP 
        ccmissing <- read.table(file=paste0(outputPrefix.tmp, ".missing"), 
                                header=TRUE, sep="")  
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
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' femaleChrXmissCutoff <- 0.05
#' inputPrefix <- "controlData"  
#' outputPrefix <- "2_10_removedSnpFemaleChrXmiss" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedSnpFemaleChrXmiss(plink, femaleChrXmissCutoff, 
#' ##                          inputPrefix, outputPrefix)


removedSnpFemaleChrXmiss <- function(plink, femaleChrXmissCutoff, 
                                      inputPrefix, outputPrefix){ 

    bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
    chrCodes <- names(table(bim[,1]))
    chr23check <- length(grep(23, chrCodes))

    if (chr23check == 1) { 
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
#' Hardy Weinberg Equilibrium test for autosomal SNPs
#'
#' @description
#' Remove autosomal SNPs deviating from Hardy Weinberg Equilibrium (HWE).
 
#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param pval the p-value cutoff for controlling HWE test in either control or 
#' case subjects. Only autosomal SNPs are considered. The default value is
#' 0.000001.
#' @param outputPvalFile the output pure text file that stores autosomal SNPs  
#' and their sorted HWE p-values.
#' @param outputSNPfile the output pure text file that stores the removed SNPs, 
#' one per line.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after HWE test on the autosome.
#' @details 

#' @export 
#' @author Junfang Chen 
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' groupLabel <- "control"
#' inputPrefix <- "controlData" ## Specify the input PLINK file prefix 
#' outputPvalFile <- "2_11_snpHwePvalAuto.txt"
#' outputSNPfile <- "2_11_snpRemovedHweAuto.txt" 
#' outputPrefix <- "2_11_removedSnpHweAuto" 
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedSnpHWEauto(groupLabel, plink, inputPrefix, pval=0.000001,
#' ##                   outputPvalFile, outputSNPfile, outputPrefix)

removedSnpHWEauto <- function(groupLabel, plink, inputPrefix, 
                              pval=0.000001, outputPvalFile, 
                              outputSNPfile, outputPrefix){ 

    if (groupLabel == "control"){ 
        groupStatus <- "filter-controls"
        affection <- "UNAFF" 
    } else if (groupLabel == "case"){
        groupStatus <- "filter-cases"
        affection <- "AFF" 
    } else {
        "ERROR: during HWE test on autosome! Wrong label warning!"
    }

    outputPrefix.tmp <- paste0(outputPrefix, "tmp")
    system(paste0(plink, " --bfile ", inputPrefix, " --",
           groupStatus, " --hardy --autosome --make-bed --out ", 
           outputPrefix.tmp))  

    ## read HWE p values 
    hweCheck <- read.table(file=paste0(outputPrefix.tmp, ".hwe"), 
                           header=TRUE, stringsAsFactors=FALSE) 

    hwe <- hweCheck[which(hweCheck$TEST == affection), ]  
    snpHweValAuto <- subset(hwe, select=c(SNP, P))
    snpHweValAuto <- snpHweValAuto[order(snpHweValAuto[,"P"]),] 
    write.table(snpHweValAuto, file=outputPvalFile, quote=FALSE, 
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
    removedSNPhwe <- hwe[which(hwe$P <= pval), "SNP"]
    write.table(removedSNPhwe, file=outputSNPfile, quote=FALSE, 
                row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

    ## exclude SNPs 
    system(paste0(plink, " --bfile ", inputPrefix, " --exclude ", 
           outputSNPfile, " --make-bed --out ", outputPrefix)) 
    system(paste0("rm ", outputPrefix.tmp, ".*"))

}

##########################################   
##########################################  
#' Hardy Weinberg Equilibrium test for chromosome X SNPs in female controls. 
#'
#' @description
#' Hardy Weinberg Equilibrium test for SNPs on the chromosome X in 
#' female controls.  

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param pval the p-value cutoff for controlling HWE test in female control 
#' subjects. Only chromosome X SNPs are considered. 
#' The default value is 0.000001.
#' @param outputPvalFile the output pure text file that stores chromosome X 
#' SNPs and their sorted HWE p-values.
#' @param outputSNPfile the output pure text file that stores the removed SNPs, 
#' one per line.
#' @param outputPrefix the prefix of the output PLINK binary files.

#' @return The output PLINK binary files after HWE test on chromosomal X 
#' in female controls.

#' @export 
#' @author Junfang Chen 
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))   
#' inputPrefix <- "controlData"  
#' outputPvalFile <- "2_12_snpHwePvalfemaleXct.txt" 
#' outputSNPfile <- "2_12_snpRemovedHweFemaleXct.txt" 
#' outputPrefix <- "2_12_removedSnpHweFemaleX"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## removedSnpFemaleChrXhweControl(plink, inputPrefix, pval=0.000001,
#' ##                                outputPvalFile, outputSNPfile, 
#' ##                                outputPrefix)

removedSnpFemaleChrXhweControl <- function(plink, inputPrefix, pval=0.000001, 
                                           outputPvalFile, outputSNPfile, 
                                           outputPrefix){ 

    bim <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
    chrCodes <- names(table(bim[,1]))
    chr23check <- length(grep(23, chrCodes))

    if (chr23check == 1) { 
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
#' The default is 20.
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
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" 
#' outputPC4subjFile <- "2_13_eigenvalAfterQC.txt"
#' outputPCplotFile <- "2_13_eigenvalAfterQC.png" ## png format
#' ## Not run: Requires an executable program GCTA, e.g.
#' ## gcta <- "/home/tools/gcta64"
#' ## plotPCA4plink(gcta, inputPrefix, nThread=20, 
#' ##               outputPC4subjFile, outputPrefix)

plotPCA4plink <- function(gcta, inputPrefix, nThread=20, 
                          outputPC4subjFile, outputPCplotFile){ 

    autosomefn <- paste0(inputPrefix, "Autosome")
    system(paste0(gcta, " --bfile ", inputPrefix, 
           " --make-grm --autosome --out ", 
           autosomefn, " --thread-num ", nThread))
    system(paste0(gcta, " --grm ", autosomefn, " --pca 20 --out ", 
           autosomefn, " --thread-num ", nThread))

    eigen <- read.table(file=paste0(autosomefn,".eigenvec"), 
                        stringsAsFactors=FALSE)
    pcs <- eigen[,seq_len(4)] ## first two PCs in the 3rd and 4th column.
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
#' @param nThread the number of threads used for computation. 
#' The default is 20.
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
#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData"  
#' cutoff <- NULL ## no outlier to be removed
#' cutoffSign <- "greater" ## not used if cutoff == NULL
#' inputPC4subjFile <- "2_13_eigenvalAfterQC.txt"
#' outputPC4outlierFile <- "2_13_eigenval4outliers.txt"
#' outputPCplotFile <- "2_13_removedOutliers.png"
#' outputPrefix <- "2_13_removedOutliers" 
#' ## Not run: Requires an executable program PLINK and GCTA, e.g.
#' ## plink <- "/home/tools/plink"
#' ## gcta <- "/home/tools/gcta64"
#' ## removeOutlierByPCs(plink, gcta, inputPrefix, nThread=20, 
#' ##                    cutoff, cutoffSign, inputPC4subjFile, 
#' ##                    outputPC4outlierFile, outputPCplotFile, outputPrefix)

removeOutlierByPCs <- function(plink, gcta, inputPrefix, nThread=20, cutoff, 
                               cutoffSign, inputPC4subjFile, 
                               outputPC4outlierFile, 
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
            outlierID <- rbind(outliersPC1v1, outliersPC1v2)
        } else { 
            if (cutoffSign == "smaller"){ 
                ## detected by PC1
                outlierID <- subjID_PCs[which(subjID_PCs[,3] <= cutoff), ] 
            } else if (cutoffSign == "greater"){
                ## detected by PC1
                outlierID <- subjID_PCs[which(subjID_PCs[,3] >= cutoff), ]
            }
        }      
        ## sorted by first PC.
        outlierIDSorted <- outlierID[order(outlierID[,3]), ] 
        write.table(outlierIDSorted, file=outputPC4outlierFile, quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
        subjID4outlierTmp  <- outlierIDSorted[,c("V1", "V2")]
        subjID4outlierTmpFile <- "subjID4outlierTmp.txt"
        write.table(subjID4outlierTmp, file=subjID4outlierTmpFile, quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
        system(paste0(plink, " --bfile ", inputPrefix, " --remove ", 
               subjID4outlierTmpFile, " --make-bed --out ", outputPrefix))
        system(paste0("rm ", subjID4outlierTmpFile)) 
        ## Plot first two PCs again; PCs for the retained subjects 
        outputPC4subjFiletmp <- "outputPC4subjFile.txt" 
        plotPCA4plink(gcta, inputPrefix=outputPrefix, nThread, 
                      outputPC4subjFiletmp, outputPCplotFile)
        system(paste0("rm ", outputPC4subjFiletmp))
    } 
}



######################################################
######################################################
#' Quality control for genotype data
#'
#' @description
#' Perform quality control on the genotype data. 

#' @param plink an executable program in either the current working directory 
#' or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param snpMissCutOffpre the cutoff of the missingness for removing SNPs 
#' before subject removal. The default is 0.05.
#' @param sampleMissCutOff the cutoff of the missingness for removing 
#' subjects/instances. The default is 0.02.
#' @param Fhet the cutoff of the autosomal heterozygosity deviation. 
#' The default is 0.2.
#' @param snpMissCutOffpost the cutoff of the missingness for removing SNPs 
#' after subject removal. The default is 0.02. 
#' @param snpMissDifCutOff the cutoff of the difference in missingness between 
#' cases and controls. The default is 0.02.
#' @param femaleChrXmissCutoff the cutoff of the missingness in female 
#' chromosome X SNPs. The default is 0.05.
#' @param pval4autoCtl the p-value cutoff for controlling HWE test in either 
#' control or case subjects. Only autosomal SNPs are considered. 
#' The default is 0.000001
#' @param pval4femaleXctl the p-value cutoff for controlling HWE test in 
#' female control subjects. Only chromosome X SNPs are considered. 
#' The default is 0.000001
#' @param outputPrefix the prefix of the output PLINK binary files after QC.
#' @param keepInterFile a logical value indicating if the intermediate 
#' processed files should be kept or not. The default is TRUE.

#' @return The output PLINK binary files after QC.
#' @details The original PLINK files are implicitly processed by the following 
#' default steps: 
#' 1.) Set all heterozygous alleles of SNPs on male chrX as missing;
#' 2.) SNP missingness < 0.05 (before sample removal);
#' 3.) Subject missingness < 0.02;  
#' 4.) Remove subjects with |Fhet| >= 0.2;
#' 5.) Reset paternal and maternal codes;
#' 6.) SNP missingness < 0.02 (after sample removal);
#' 7.) Remove SNPs with difference >= 0.02 of SNP missingness 
#' between cases and controls;
#' 8.) Remove chrX SNPs with missingness >= 0.05 in females.
#' (Optional, if no chrX data);
#' 9.) Remove autosomal SNPs with HWE p < 10-6 in controls;
#' 10.) Remove chrX SNPs with HWE p < 10-6 in female controls. 
#' (Optional, if no chrX data). 

#' @export  
#' @author Junfang Chen 
#' @references Schizophrenia Working Group of the Psychiatric Genomics, C. 
#' (2014). Biological insights from 108 schizophrenia-associated genetic loci. 
#' Nature 511(7510): 421-427. 

#' @examples  
#' ## In the current working directory
#' bedFile <- system.file("extdata", "controlData.bed", package="Gimpute")
#' bimFile <- system.file("extdata", "controlData.bim", package="Gimpute") 
#' famFile <- system.file("extdata", "controlData.fam", package="Gimpute")
#' system(paste0("scp ", bedFile, bimFile, famFile, " ."))  
#' inputPrefix <- "controlData" 
#' outputPrefix <- "2_12_removedSnpHweFemaleX"  
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## genoQC(plink, inputPrefix, 
#' ##        snpMissCutOffpre=0.05, 
#' ##        sampleMissCutOff=0.02, 
#' ##        Fhet=0.2, 
#' ##        snpMissCutOffpost=0.02, 
#' ##        snpMissDifCutOff=0.02,
#' ##        femaleChrXmissCutoff=0.05, 
#' ##        pval4autoCtl=0.000001, 
#' ##        pval4femaleXctl=0.000001, 
#' ##        outputPrefix, keepInterFile=TRUE)


genoQC <- function(plink, inputPrefix, snpMissCutOffpre=0.05, 
                   sampleMissCutOff=0.02, Fhet=0.2, 
                   snpMissCutOffpost=0.02, snpMissDifCutOff=0.02,
                   femaleChrXmissCutoff=0.05, pval4autoCtl=0.000001, 
                   pval4femaleXctl=0.000001, 
                   outputPrefix, keepInterFile=TRUE) {

    ## check case control status
    groupLabel <- getGroupLabel(inputFAMfile=paste0(inputPrefix, ".fam"))
    ## step 3 
    # 3. Set all heterozygous alleles of SNPs of the chr 23 for males
    outputPrefix3 <- "2_03_setHeteroHaploMissing" 
    setHeteroHaploMissing(plink, inputPrefix, outputPrefix=outputPrefix3) 
    ## step 4  SNP missingness < 0.05 (before sample removal);  
    outputPrefix4 <- "2_04_removedSnpMissPre" 
    removedSnpMiss(plink, snpMissCutOff=snpMissCutOffpre, 
                   inputPrefix=outputPrefix3, outputPrefix=outputPrefix4)        
    ## step 5 
    # subject missingness < 0.02;  
    outputPrefix5 <- "2_05_removedInstMiss" 
    removedInstMiss(plink, sampleMissCutOff,  
                    inputPrefix=outputPrefix4, outputPrefix=outputPrefix5)
    ## step 6  
    outputPrefix6 <- "2_06_removedInstFhet" 
    removedInstFhet(plink, Fhet, 
                    inputPrefix=outputPrefix5, outputPrefix=outputPrefix6)
    ## step 7
    outputPrefix7 <- "2_07_removedParentIdsMiss" 
    removedParentIdsMiss(plink, inputPrefix=outputPrefix6, 
                         outputPrefix=outputPrefix7)
    ## step 8  
    outputPrefix8 <- "2_08_removedSnpMissPost" 
    removedSnpMiss(plink, snpMissCutOff=snpMissCutOffpost, 
                   inputPrefix=outputPrefix7, outputPrefix=outputPrefix8)
    ## step 9  
    outputPrefix9 <- "2_09_removedSnpMissDiff" 
    removedSnpMissDiff(plink, inputPrefix=outputPrefix8, 
                       snpMissDifCutOff, outputPrefix=outputPrefix9, groupLabel) 
    ## step 10
    outputPrefix10 <- "2_10_removedSnpFemaleChrXmiss" 
    removedSnpFemaleChrXmiss(plink, femaleChrXmissCutoff, 
                             inputPrefix=outputPrefix9, 
                             outputPrefix=outputPrefix10) 
    ## step 11
    outputPvalFile <- "2_11_snpHwePvalAuto.txt" 
    outputSNPfile <-  "2_11_snpRemovedHweAuto.txt" 
    outputPrefix11 <- "2_11_removedSnpHweAuto" 
    if (groupLabel == "control" | groupLabel == "caseControl"){
        removedSnpHWEauto(groupLabel="control", plink, 
                          inputPrefix=outputPrefix10, 
                          pval=pval4autoCtl, outputPvalFile, 
                          outputSNPfile, outputPrefix=outputPrefix11)
    ## HWE test is only performed on control data 
    } else { print("ERROR: HWE test on autosome!") }
    ## step 12 
    outputPvalFile <- "2_12_snpHwePvalfemaleXct.txt" 
    outputSNPfile <- "2_12_snpRemovedHweFemaleXct.txt"  
    removedSnpFemaleChrXhweControl(plink, inputPrefix=outputPrefix11, 
                                   pval=pval4femaleXctl, outputPvalFile,
                                   outputSNPfile, outputPrefix=outputPrefix)
    if (keepInterFile==FALSE){ 
        ## remove intermediate files 
        system(paste0("rm ", outputPrefix3, ".*"))
        system(paste0("rm ", outputPrefix4, ".*"))
        system(paste0("rm ", outputPrefix5, ".*"))
        system(paste0("rm ", outputPrefix6, ".*"))
        system(paste0("rm ", outputPrefix7, ".*"))
        system(paste0("rm ", outputPrefix8, ".*"))
        system(paste0("rm ", outputPrefix9, ".*"))
        system(paste0("rm ", outputPrefix10, ".*"))
        system(paste0("rm ", outputPrefix11, ".*")) 

        if (file.exists(outputPvalFile)) {  
            system(paste0("rm ", outputPvalFile))
        }
        if (file.exists(outputSNPfile)) {  
            system(paste0("rm ", outputSNPfile))
        }    
    }
    ## remove unwanted files
    system(paste0("rm  *.log "))   
} 