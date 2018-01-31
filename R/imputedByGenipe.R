##########################################################################
########################################################################## imputedByGenipe.R
#' Impute genotypes using Genipe
#'
#' @description
#' Perform imputation by Genipe for the autosomal and sex chromosome prephased known haplotypes with a reference panel. 
#' Note that pre-phasing using SHAPEIT is done without the refernce haplotypes.

#' @param chrs  specifiy the chromosomes for imputation. There are four different options ("autosomes", "1, or 2, or 3...or 22", "23", "25"). 
#' 1.) 'autosomes': will impute chromosome 1 to 22 together; 2.) 'chrs' belongs to one of (1, 2,...22) then the imputation is done just for one autosomal chromosome.
#' 3.) '23' will do the imputation for the non-pseudoautosomal region of chromosome 23. 4.) '25' imputes for the pseudoautosomal regions of chromosome 23.
#' @param impRefDIR4genipe the directory where the imputation reference files are located.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param shapeit an executable SHAPEIT binary program in either the current working directory or somewhere in the command path.
#' @param impute2 an executable IMPUTE2 binary program in either the current working directory or somewhere in the command path.
#' @param fastaFile the human reference files for the initial strand check.
#' @param segmentSize the length of a single segment for imputation.
#' @param thread4impute2 the number of threads for the imputation.  
#' @param thread4shapeit the number of threads for phasing.
  
#' @return  The imputed files using genipe.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
  

imputedByGenipe <- function(chrs, impRefDIR4genipe, inputPrefix, shapeit, impute2, plink, fastaFile, segmentSize, thread4impute2, thread4shapeit){ 
 
       
        # ## IMPUTE2 autosomal reference
        # HAPS_FILE = paste0(impRefDIR4genipe, '1000GP_Phase3_chr', i, '.hap.gz')  
        # LEGEND_FILE = paste0(impRefDIR4genipe, '1000GP_Phase3_chr', i, '.legend.gz')  
        # GENMAP_FILE = paste0(impRefDIR4genipe, 'genetic_map_chr', i, '_combined_b37.txt')   
        HAPS_FILE = paste0(impRefDIR4genipe, '1000GP_Phase3_chr{chrom}.hap.gz')   
        LEGEND_FILE = paste0(impRefDIR4genipe, '1000GP_Phase3_chr{chrom}.legend.gz') 
        GENMAP_FILE = paste0(impRefDIR4genipe, 'genetic_map_chr{chrom}_combined_b37.txt')
        sampleFile = paste0(impRefDIR4genipe, '1000GP_Phase3.sample')

        ## chrX non-PAR
        HAPS_FILEchrXnonPAR = paste0(impRefDIR4genipe, '1000GP_Phase3_chrX_NONPAR.hap.gz')   
        LEGEND_FILEchrXnonPAR = paste0(impRefDIR4genipe, '1000GP_Phase3_chrX_NONPAR.legend.gz') 
        GENMAP_FILEchrXnonPAR = paste0(impRefDIR4genipe, 'genetic_map_chrX_nonPAR_combined_b37.txt') 
        
        ## chrX X_PAR1
        HAPS_FILEchrX_PAR1 = paste0(impRefDIR4genipe, '1000GP_Phase3_chrX_PAR1.hap.gz')   
        LEGEND_FILEchrX_PAR1= paste0(impRefDIR4genipe, '1000GP_Phase3_chrX_PAR1.legend.gz') 
        GENMAP_FILEchrX_PAR1 = paste0(impRefDIR4genipe, 'genetic_map_chrX_PAR1_combined_b37.txt') 

        ## chrX X_PAR2
        HAPS_FILEchrX_PAR2 = paste0(impRefDIR4genipe, '1000GP_Phase3_chrX_PAR2.hap.gz')   
        LEGEND_FILEchrX_PAR2 = paste0(impRefDIR4genipe, '1000GP_Phase3_chrX_PAR2.legend.gz') 
        GENMAP_FILEchrX_PAR2 = paste0(impRefDIR4genipe, 'genetic_map_chrX_PAR2_combined_b37.txt') 


        if (chrs=='autosomes'){ 
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
        } else if (is.element(chrs, 1:22)) { ## one chromosome at a time
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
            " --hap-nonPAR  ", HAPS_FILEchrXnonPAR, " \ ",  ## chrX non-PAR
            " --legend-nonPAR ", LEGEND_FILEchrXnonPAR, " \ ", ## chrX non-PAR 
            " --map-nonPAR ", GENMAP_FILEchrXnonPAR, " \ ", ## chrX non-PAR
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
        } else if (chrs==25){
            system( paste0("genipe-launcher ", 
            " --chrom 25 \ ", 
            " --bfile ", inputPrefix, " \ ",
            " --shapeit-bin ", shapeit, " \ ", 
            " --impute2-bin ", impute2, " \ ", 
            " --plink-bin ", plink, " \ ", 
            " --reference ", fastaFile, " \ ", 
            " --hap-PAR1  ", HAPS_FILEchrX_PAR1, " \ ",  ## chrX PAR1
            " --legend-PAR1 ", LEGEND_FILEchrX_PAR1, " \ ", ## chrX PAR1 
            " --map-PAR1 ", GENMAP_FILEchrX_PAR1, " \ ", ## chrX PAR1
            " --hap-PAR2  ", HAPS_FILEchrX_PAR2, " \ ",  ## chrX PAR2
            " --legend-PAR2 ", LEGEND_FILEchrX_PAR2, " \ ", ## chrX PAR2 
            " --map-PAR2 ", GENMAP_FILEchrX_PAR2, " \ ", ## chrX PAR2
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
 
