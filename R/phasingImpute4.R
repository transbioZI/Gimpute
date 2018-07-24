

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





#' Phasing and imputation 
#'
#' @description
#' Perform phasing, imputation and conversion from IMPUTE2 format into PLINK 
#' binary files. 

#' @param inputPrefix the prefix of the input PLINK binary files for 
#' the imputation. 
#' @param outputPrefix the prefix of the output PLINK binary files  
#' after imputation.
#' @param plink an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param shapeit an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param impute4 an executable program in either the current 
#' working directory or somewhere in the command path.
#' @param qctool an executable program in either the current working 
#' directory or somewhere in the command path.
#' @param gtool an executable program in either the current 
#' working directory or somewhere in the command path. 
#' @param windowSize  the window size of each chunk. 
#' The default value is 3000000. 
#' @param effectiveSize this parameter controls the effective population size. 
#' Commonly denoted as Ne. A universal -Ne value of 20000 is suggested.
#' @param nCore the number of cores used for splitting chromosome by PLINK, 
#' phasing, imputation, genotype format modification, genotype conversion, and 
#' merging genotypes. The default value is 40. 
#' @param threshold threshold for merging genotypes from GEN probability. 
#' Default 0.9. 
#' @param outputInfoFile the output file of info scores consisting of 
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
#' @details Currently, chromosome X is not supported for the impute4.
#' The whole imputation process mainly consists of the following 
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
#'   \item Bycroft, C., et al. Genome-wide genetic data on~ 500,000 UK Biobank 
#'         participants. BioRxiv (2017): 166298.
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
#' outputPrefix <- "gwasImputed"   
#' outputInfoFile <- "infoScore.txt"
#' tmpImputeDir <- "tmpImpute"
#' ## Not run: Requires an executable program PLINK, e.g.
#' ## plink <- "/home/tools/plink"
#' ## phaseImpute4(inputPrefix, outputPrefix, prefix4final,
#' ##             plink, shapeit, impute4, qctool, gtool,
#' ##             windowSize=3000000, effectiveSize=20000, 
#' ##             nCore=40, threshold=0.9, outputInfoFile,
#' ##             referencePanel, impRefDIR, tmpImputeDir, keepTmpDir=TRUE)



phaseImpute4 <- function(inputPrefix, outputPrefix,
                         plink, shapeit, impute4, qctool, gtool, 
                         windowSize=3000000, effectiveSize=20000, 
                         nCore=40, threshold=0.9, outputInfoFile,
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
                        chrXPAR2suffix, nCore)
    print(PAR)  
    if (PAR[[1]]) {par1 <- "X_PAR1"} else {par1 <- NULL}
    if (PAR[[2]]) {par2 <- "X_PAR2"} else {par2 <- NULL}
    ## step 2.2
    chunkPrefix <- "chunks_chr" 
    chrs <- c(currentChr, par1, par2)    
    if (is.element(23, chrs) == TRUE) {
        chrs <- setdiff(chrs, 23) ## chrX is not available for impute4.
    } 
    chunk4eachChr(inputPrefix=prefix4eachChr, 
                  outputPrefix=chunkPrefix, chrs, windowSize) 

    setwd("..") 
    system(paste0("mv ", dataDIR, chunkPrefix, "*.txt  ", chunkDIR)) 
    ## step 2.3     
    .prePhasingByShapeit(shapeit, chrs, dataDIR, 
                         prefix4eachChr, referencePanel, 
                         impRefDIR, phaseDIR, nThread=nCore, 
                         effectiveSize, nCore=1) 

    ## step 2.4  ############################# 
    prefixChunk <- paste0(chunkDIR, chunkPrefix)        
    .imputedByImpute4(impute4, chrs, prefixChunk, phaseDIR, 
                      referencePanel, impRefDIR, 
                      imputedDIR, prefix4eachChr, 
                      nCore, effectiveSize)
    ## step 2.5  #############################  
    ## extract only SNPs (without INDELs)
    #######################################################
    setwd(imputedDIR)  
    ## extract only SNPs starting with "rs";  .
    
    ls <- system("ls gwas*.gen", intern=TRUE)
    snpPrefix <- "rs" 
    biglists <- as.list(ls)
    mclapply(biglists, function(i){ 
        arg1 <- paste0(" awk '{if(length($4) == 1 && length($5) == 1) print}'")
        arg2 <- paste0(i, "NoINDEL.gen")   
        system(paste0("grep '", snpPrefix, "' ", i, " | ", arg1, " > ", arg2))
    }, mc.cores=nCore) ## by default 

    inputSuffix <- ".genNoINDEL.gen" 
    ## compute info score for each chunk, then combine all info scores
    computeInfoByQctool(qctool, inputSuffix, outputInfoFile)
    
    setwd("..") 
    suffix4imputed <- ".genNoINDEL.gen"  
    .convertImpute2ByGtool(gtool, chrs, prefixChunk, phaseDIR, imputedDIR, 
                           prefix4eachChr, suffix4imputed, 
                           postImputeDIR, threshold, nCore)
    ## step 2.6  
    ####################################################### 
    ## Modify missing genotype format.
    setwd(postImputeDIR)   
    ## replace 'N' in the .ped files into 0 > missing values.
    chrslist <- as.list(chrs) 
    fn <- mclapply(chrslist, function(i){
        system(paste0("sed -i 's/N/0/g' ", prefix4eachChr, i, ".*ped "))
    }, mc.cores=nCore) ## by default
    prefixMerge <- "gwasMerged" 
    .mergePlinkData(plink, chrs, prefix4eachChr, 
                    prefixMerge, nCore)
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
           " --update-ids recoded.txt --make-bed --out ", outputPrefix)) 
    #######################################################
    setwd("..")
    setwd("..")
    system(paste0("mv ", tmpImputeDir, "/", postImputeDIR, outputPrefix, ".* .")) 
    system(paste0("mv ", tmpImputeDir, "/", imputedDIR, outputInfoFile,  " .")) 

    if (keepTmpDir == FALSE){
        system(paste0("rm -r ", tmpImputeDir))
    } else if (keepTmpDir == TRUE){
        print("Keep the intermediate imputation folder.")
    }

}

