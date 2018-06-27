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
##########################################################################
########################################################################## 
## imputedByGenipe.R
#' Impute genotypes using Genipe
#'
#' @description
#' Perform imputation by Genipe for the autosomal and sex chromosome prephased 
#' known haplotypes with a reference panel. 
#' Note that pre-phasing using SHAPEIT is done without the refernce haplotypes.

#' @param chrs  specifiy the chromosomes for imputation. There are 
#' four different options ("autosomes", "1, or 2, or 3...or 22", "23", "25"). 
#' 1.) 'autosomes': will impute chromosome 1 to 22 together; 2.) 'chrs' 
#' belongs to one of (1, 2,...22) then the imputation is done just for 
#' one autosomal chromosome.
#' 3.) '23' will do the imputation for the non-pseudoautosomal region of 
#' chromosome 23. 
#' 4.) '25' imputes for the pseudoautosomal regions of chromosome 23.
#' @param impRefDir the directory where the imputation reference files 
#' are located.
#' @param inputPrefix the prefix of the input PLINK binary files.
#' @param shapeit an executable SHAPEIT binary program in either the current 
#' working directory or somewhere in the command path.
#' @param impute2 an executable IMPUTE2 binary program in either the current 
#' working directory or somewhere in the command path.
#' @param fastaFile the human reference files for the initial strand check.
#' @param segmentSize the length of a single segment for imputation.
#' @param thread4impute2 the number of threads for the imputation.  
#' @param thread4shapeit the number of threads for phasing.

#' @return The imputed files using genipe.
#' @export 

#' @author Junfang Chen 
#' @references  Lemieux Perreault, L. P., et al. (2016). genipe: an automated 
#' genome-wide imputation pipeline with automatic reporting and statistical 
#' tools. Bioinformatics, 32(23), 3661-3663.

###' @examples 


imputedByGenipe <- function(chrs, impRefDir, inputPrefix, 
                            shapeit, impute2, plink, fastaFile, 
                            segmentSize, thread4impute2, thread4shapeit){ 
    
    # ## IMPUTE2 autosomal reference   
    HAPS_FILE <- paste0(impRefDir, "1000GP_Phase3_chr{chrom}.hap.gz")   
    LEGEND_FILE <- paste0(impRefDir, "1000GP_Phase3_chr{chrom}.legend.gz") 
    GENMAP_FILE <- paste0(impRefDir, "genetic_map_chr{chrom}_combined_b37.txt")
    sampleFile <- paste0(impRefDir, "1000GP_Phase3.sample")

    ## chrX non-PAR
    HAPSfileChrXnonPAR <- paste0(impRefDir, 
                                 "1000GP_Phase3_chrX_NONPAR.hap.gz")   
    LEGENDXnonPAR <- paste0(impRefDir, "1000GP_Phase3_chrX_NONPAR.legend.gz") 
    GENMAPfileXnonPAR <- paste0(impRefDir, 
                                "genetic_map_chrX_nonPAR_combined_b37.txt") 
    
    ## chrX X_PAR1
    HAPS_FILEchrX_PAR1 <- paste0(impRefDir, "1000GP_Phase3_chrX_PAR1.hap.gz")   
    LEGENDfileX_PAR1 <- paste0(impRefDir, "1000GP_Phase3_chrX_PAR1.legend.gz") 
    GENMAPX_PAR1 <- paste0(impRefDir, "genetic_map_chrX_PAR1_combined_b37.txt") 

    ## chrX X_PAR2
    HAPS_FILEchrX_PAR2 <- paste0(impRefDir, "1000GP_Phase3_chrX_PAR2.hap.gz")   
    LEGENDfileX_PAR2 <- paste0(impRefDir, "1000GP_Phase3_chrX_PAR2.legend.gz") 
    GENMAPX_PAR2 <- paste0(impRefDir, "genetic_map_chrX_PAR2_combined_b37.txt") 

    autosomeCode = seq_len(22) 

    if (chrs == "autosomes"){ 
        system( paste0("genipe-launcher ", 
        " --chrom autosomes \ ", 
        " --bfile ", inputPrefix, " \ ",
        " --shapeit-bin ", shapeit, " \ ", 
        " --impute2-bin ", impute2, " \ ", 
        " --plink-bin ", plink, " \ ", 
        " --reference ", fastaFile, " \ ", 
        " --hap-template  ", HAPS_FILE, " \ ", 
        " --legend-template ", LEGEND_FILE, " \ ", 
        " --map-template ", GENMAP_FILE, " \ ", 
        " --sample-file ", sampleFile, " \ ", 
        " --segment-length ", segmentSize, " \ ", 
        " --thread ", thread4impute2, " \ ", 
        " --shapeit-thread ", thread4shapeit, " \ ", 
        " --report-title Tutorial \ ",  
        " --report-number 'TestReport' \ " ) )  
    } else if (is.element(chrs, autosomeCode)) { ## one chromosome at a time
        system( paste0("genipe-launcher ", 
        " --chrom ", chrs, " \ ", 
        " --bfile ", inputPrefix, " \ ",
        " --shapeit-bin ", shapeit, " \ ", 
        " --impute2-bin ", impute2, " \ ", 
        " --plink-bin ", plink, " \ ", 
        " --reference ", fastaFile, " \ ",  
        " --hap-template  ", HAPS_FILE, " \ ", 
        " --legend-template ", LEGEND_FILE, " \ ", 
        " --map-template ", GENMAP_FILE, " \ ", 
        " --sample-file ", sampleFile, " \ ", 
        " --segment-length ", segmentSize, " \ ", 
        " --thread ", thread4impute2, " \ ", 
        " --shapeit-thread ", thread4shapeit, " \ ", 
        " --report-title Tutorial \ ",  
        " --report-number 'TestReport' \ "  
        ) )   
    } else if (chrs == 23) { 
        system( paste0("genipe-launcher ", 
        " --chrom 23 \ ", 
        " --bfile ", inputPrefix, " \ ",
        " --shapeit-bin ", shapeit, " \ ", 
        " --impute2-bin ", impute2, " \ ", 
        " --plink-bin ", plink, " \ ", 
        " --reference ", fastaFile, " \ ", 
        " --hap-nonPAR  ", HAPSfileChrXnonPAR, " \ ",  ## chrX non-PAR
        " --legend-nonPAR ", LEGENDXnonPAR, " \ ", ## chrX non-PAR 
        " --map-nonPAR ", GENMAPfileXnonPAR, " \ ", ## chrX non-PAR
        " --hap-template  ", HAPS_FILE, " \ ", 
        " --legend-template ", LEGEND_FILE, " \ ", 
        " --map-template ", GENMAP_FILE, " \ ", 
        " --sample-file ", sampleFile, " \ ", 
        " --segment-length ", segmentSize, " \ ", 
        " --thread ", thread4impute2, " \ ", 
        " --shapeit-thread ", thread4shapeit, " \ ", 
        " --report-title Tutorial \ ",  
        " --report-number 'TestReport' \ "  
        ) )   
    } else if (chrs == 25){
        system( paste0("genipe-launcher ", 
        " --chrom 25 \ ", 
        " --bfile ", inputPrefix, " \ ",
        " --shapeit-bin ", shapeit, " \ ", 
        " --impute2-bin ", impute2, " \ ", 
        " --plink-bin ", plink, " \ ", 
        " --reference ", fastaFile, " \ ", 
        " --hap-PAR1  ", HAPS_FILEchrX_PAR1, " \ ",  ## chrX PAR1
        " --legend-PAR1 ", LEGENDfileX_PAR1, " \ ", ## chrX PAR1 
        " --map-PAR1 ", GENMAPX_PAR1, " \ ", ## chrX PAR1
        " --hap-PAR2  ", HAPS_FILEchrX_PAR2, " \ ",  ## chrX PAR2
        " --legend-PAR2 ", LEGENDfileX_PAR2, " \ ", ## chrX PAR2 
        " --map-PAR2 ", GENMAPX_PAR2, " \ ", ## chrX PAR2
        " --hap-template  ", HAPS_FILE, " \ ", 
        " --legend-template ", LEGEND_FILE, " \ ", 
        " --map-template ", GENMAP_FILE, " \ ", 
        " --sample-file ", sampleFile, " \ ", 
        " --segment-length ", segmentSize, " \ ", 
        " --thread ", thread4impute2, " \ ", 
        " --shapeit-thread ", thread4shapeit, " \ ", 
        " --report-title Tutorial \ ",  
        " --report-number 'TestReport' \ "  
        ) )  
    }
}



##########################################################################
########################################################################## 
## mergeByGenipe.R
#' Merge imputed files using Genipe

#' @description
#' Concatenate IMPUTE2 output files and retrieve some statistics. 
#' This is automatically called by the main genipe pipeline to merge IMPUTE2 
#' files generated for all the genomic segments. 
#' For details, see Genipe IMPUTE2 merger options.

#' @param chr  specifiy the chromosome segment to be merged, on which the 
#' imputation was made.
#' @param inputImpute2 the output from IMPUTE2.
#' @param probability the probability threshold for no calls. [<0.9]
#' @param completionRate the completion rate threshold for site exclusion. 
#' [<0.98]
#' @param info the measure of the observed statistical information associated 
#' with the allele frequency estimate threshold for site exclusion. (<0.00)
#' @param outputPrefix the prefix for the imputed output files. 

#' @return The merged imputed files using genipe.
#' @export 

#' @author Junfang Chen 
#' @references  Lemieux Perreault, L. P., et al. (2016). genipe: an automated 
#' genome-wide imputation pipeline with automatic reporting and statistical 
#' tools. Bioinformatics, 32(23), 3661-3663.

##' @examples 

mergeByGenipe <- function(inputImpute2, chr, probability, 
                          completionRate, info, outputPrefix){ 
    system( paste0("impute2-merger ", 
    " --impute2 ", inputImpute2, " \ ",
    " --chr ", chr, "\ ",   
    " --probability  ", probability, " \ ", 
    " --completion ", completionRate, " \ ", 
    " --info  ", info, " \ ",  
    " --prefix ", outputPrefix ) )    
}






##########################################################################
########################################################################## 
## extractByGenipe.R


#' Extract imputed markers using Genipe

#' @description
#' Extract imputed markers located in a specific genomic region using Genipe. 
#' Note that, 1.) 'bed' PLINK binary format is specifically used for the 
#' output format.
#' 2.) The markers of the whole chromosome are extracted together. 
#' For the filtering of maf and info will be done during post imputation.

#' @param inputImpute2 the output from IMPUTE2.
#' @param inputMAP the output PLINK MAP file from Genipe, which will be used 
#' for generating markers in a text file (only one column without column name).
# #' @param index  only perform the indexation.
#' @param outputPrefix the prefix of the output files. [impute2_extractor]
#' @param format the output format. Can specify either "impute2" for 
#' probabilities (same as impute2 format, i.e. 3 values per sample), "dosage" 
#' for dosage values (one value between 0 and 2 by sample), 
#'  "calls" for hard calls, or "bed" for Plink binary format (with hard calls). 
#' [impute2]
# #' @param long write the output file in the long format 
#' (one line per sample per marker). This option is only compatible with 
#' the "calls" and "dosage" format (option "– format").
#' @param prob the probability threshold used when creating a file 
#' in the dosage or call format. [0.9]
##' @param extractMarkerFile file containing marker names to extract.
##' @param genomic CHR:START-END The range to extract (e.g. 22 1000000 1500000). 
#' Can be use in combination with "--rate", "--maf" and "--info".
##' @param maf extract markers with a minor allele frequency equal or higher 
#' than the specified threshold. Can be use in combination with "–rate", 
#' "–info" and "–genomic".
##' @param rate extract markers with a completion rate equal or higher to 
#' the specified threshold. Can be use in combination with "--maf", "--info" 
#' and "--genomic".
##' @param info Extract markers with an information equal or higher to 
#' the specified threshold. Can be use in combination with "--maf", "--rate" 
#' and "--genomic".

#' @return The extracted imputed files using genipe.
#' @export 

#' @author Junfang Chen 
#' @references  Lemieux Perreault, L. P., et al. (2016). genipe: an automated 
#' genome-wide imputation pipeline with automatic reporting and statistical 
#' tools. Bioinformatics, 32(23), 3661-3663.

##' @examples  
 

extractByGenipe <- function(inputImpute2, inputMAP, 
                            outputPrefix, format, prob){ 

    ## extract all the markers from the input map file 
    tmpMarkerFile <- "markers.txt"
    system( paste0("awk '{print $2}' ", inputMAP, " > ", tmpMarkerFile) )
    system( paste0("impute2-extractor ", 
        " --impute2 ", inputImpute2, " \ ",  
        " --out ", outputPrefix, " \ ",  
        " --format  ", format, " \ ",  
        " --extract  ", tmpMarkerFile, " \ "
    ) )
    system(paste0("rm ", tmpMarkerFile))    
}
