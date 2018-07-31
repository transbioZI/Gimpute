

##########################################################################
## computeInfoByQctool
########################################################################## 
#' Calculate the info score by QCTOOL  
#'
#' @description
#' Calculate the info score by QCTOOL for a set of SNPs from .GEN files.
#' 
#' @param qctool an executable program in either the current working 
#' directory or somewhere in the command path.
#' @param inputSuffix the suffix of input .GEN files within current directory. 
#' @param outputInfoFile the output info scores file consisting of 
#' two columns: all SNPs from .GEN files and their info scores.    

#' @return A pure text file contains the info scores of all SNPs from .GEN  
#' files with two columns: SNP names and the corresponding info scores.  

#' @details These .GEN files may come from the output of impute4 results.
#' The intermediate generated files are retained in the directory.

#' @export 
#' @author Junfang Chen 
#' @examples 
#' ## In the current working directory
#' inputSuffix <- "noINDEL.gen"
#' outputInfoFile <- "infoScore.txt"
#' ## computeInfoByQctool(qctool, inputSuffix, outputInfoFile)



computeInfoByQctool <- function(qctool, inputSuffix, outputInfoFile){

    files <- system(paste0("ls *", inputSuffix), intern=TRUE) 
    for (i in seq_len(length(files))) {  
        out <- paste0(files[i], ".txt")
        system(paste0(qctool, " -g ", files[i], " -snp-stats -osnp ", out))
        out2 <- paste0("tmpInfo", i, ".txt")
        arg1 <- paste0("grep -o '^[^#]*' ", out) ## reomve comments
        arg2 <- paste0("awk '{print $2, $17}' > ", out2)  
        system(paste0(arg1, " | ", arg2))
    }  

    ## note, no other file names starting with tmpInfo*
    a <- "'NR==1 {header=$_} FNR==1 && NR!=1 { $_ ~ $header getline; } {print}'"
    system(paste0("awk ", a, " tmpInfo* > ", outputInfoFile))
}



 
##########################################################################
## .imputedByImpute4.R
########################################################################## 
#' Impute genotypes using IMPUTE4
#'
#' @description
#' Perform imputation by IMPUTE4 for the autosomal prephased known haplotypes 
#' with a reference panel.

#' @param impute4 an executable program in either the current 
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
#' @return The imputed files for all chunks from given chromosomes, except   
#' sex chromosomes. 

#' @export 
#' @import doParallel  

#' @author Junfang Chen 
#' @seealso \code{\link{phaseImpute4}}.


.imputedByImpute4 <- function(impute4, chrs, prefixChunk, phaseDIR, 
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
            GENMAP.chrXPAR1 <- paste0(impRefDIR, "genetic_map_chr", 
                                  "X_PAR1_combined_b37.txt ")
            GENMAP.chrXPAR2 <- paste0(impRefDIR, "genetic_map_chr", 
                                  "X_PAR2_combined_b37.txt ") 
            if (referencePanel == "1000Gphase1v3_macGT1"){ 
                ## autosome
                impPrefix <- "ALL_1000G_phase1integrated_v3_chr"

                HAPS_FILE <- paste0(impRefDIR, impPrefix, i, 
                                    "_impute_macGT1.hap.gz ") 
                LEGEND_FILE <- paste0(impRefDIR, impPrefix, i, 
                                      "_impute_macGT1.legend.gz ")  
                ## .chrXPAR1  
                HAPS.chrXPAR1 <- paste0(impRefDIR, impPrefix,
                                    "X_PAR1_impute_macGT1.hap.gz ") 
                LEGEND.chrXPAR1  <- paste0(impRefDIR, impPrefix, 
                                      "X_PAR1_impute_macGT1.legend.gz ") 
                ## .chrXPAR2
                HAPS.chrXPAR2  <- paste0(impRefDIR, impPrefix, 
                                    "X_PAR2_impute_macGT1.hap.gz ") 
                LEGEND.chrXPAR2  <- paste0(impRefDIR, impPrefix, 
                                      "X_PAR2_impute_macGT1.legend.gz ") 

            } else if (referencePanel == "1000Gphase3"){
                HAPS_FILE <- paste0(impRefDIR, "1000GP_Phase3_chr", i, 
                                    ".hap.gz ") 
                ## autosome
                LEGEND_FILE <- paste0(impRefDIR, "1000GP_Phase3_chr", i, 
                                      ".legend.gz ")
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
            ## Output file GEN format   
            OUTPUT_FILE <- paste0(imputedDIR, prefix4eachChr, i, 
                                  ".pos", chunkSTART, 
                                  "-", chunkEND)   ## .gen output  
            ################## impute genotypes from GWAS haplotypes 
            autosomeCode = seq_len(22)
            if (is.element(i, autosomeCode)) { 
                ## impute for the autosomes
                system(paste0(impute4,  
                " -no_maf_align \ ",   
                " -m ", GENMAP_FILE, " \ ",  
                " -h ", HAPS_FILE, " \ ", 
                " -l ", LEGEND_FILE, " \ ", 
                " -g ", GWAS_HAPS_FILE, " \ ", 
                " -Ne ", effectiveSize, " \ ", 
                " -int ", chunkSTART, " ", chunkEND, " \ ", 
                " -buffer 1000  \ ",
                " -o ", OUTPUT_FILE, " \ " ))
            } else if (i == "X_PAR1") {  
                ## impute for chrX PAR >> with an additional flag: --Xpar.
                system(paste0(impute4,   
                " -no_maf_align \ ",   
                " -m ", GENMAP_FILE, " \ ",  
                " -h ", HAPS.chrXPAR1, " \ ", 
                " -l ", LEGEND.chrXPAR1, " \ ", 
                " -g ", GWAS_HAPS_FILE, " \ ", 
                " -Ne ", effectiveSize, " \ ", 
                " -int ", chunkSTART, " ", chunkEND, " \ ", 
                " -buffer 1000  \ ",
                " -o ", OUTPUT_FILE, " \ " ))
            } else if (i == "X_PAR2") {  
                ## impute for chrX PAR >> with an additional flag: --Xpar.
                system(paste0(impute4,   
                " -no_maf_align \ ",   
                " -m ", GENMAP_FILE, " \ ",  
                " -h ", HAPS.chrXPAR2, " \ ", 
                " -l ", LEGEND.chrXPAR2, " \ ", 
                " -g ", GWAS_HAPS_FILE, " \ ", 
                " -Ne ", effectiveSize, " \ ", 
                " -int ", chunkSTART, " ", chunkEND, " \ ", 
                " -buffer 1000  \ ",
                " -o ", OUTPUT_FILE, " \ " ))
            } else if (i == 23){ 
                ## impute for chrX 
                print("chrX option >> Not available for now!")
            }
        }, mc.cores=nCore)  
    } 
}



 