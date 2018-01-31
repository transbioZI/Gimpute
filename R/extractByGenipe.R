
 


##########################################################################
########################################################################## extractByGenipe.R
#' Extract imputed markers using Genipe

#' @description
#' Extract imputed markers located in a specific genomic region using Genipe. Note that, 1.) 'bed' PLINK binary format is specifically used for the output format.
#' 2.) The markers of the whole chromosome are extracted together. For the filtering of maf and info will be done during post imputation.

#' @param inputImpute2 the output from IMPUTE2.
#' @param inputMAP the output PLINK MAP file from Genipe, which will be used for generating markers in a text file (only one column without column name).
# #' @param index  only perform the indexation.
#' @param outputPrefix the prefix of the output files. [impute2_extractor]
#' @param format the output format. Can specify either ‘impute2’ for probabilities 
#'  (same as impute2 format, i.e. 3 values per sample), ‘dosage’ for dosage values (one value between 0 and 2 by sample), 
#'  ‘calls’ for hard calls, or ‘bed’ for Plink binary format (with hard calls). [impute2]
# #' @param long write the output file in the long format (one line per sample per marker). This option is only compatible with the ‘calls’ and ‘dosage’ format (option ‘– format’).
#' @param prob the probability threshold used when creating a file in the dosage or call format. [0.9]
##' @param extractMarkerFile file containing marker names to extract.
##' @param genomic CHR:START-END The range to extract (e.g. 22 1000000 1500000). Can be use in combination with ‘--rate’, ‘--maf’ and ‘--info’.
##' @param maf extract markers with a minor allele frequency equal or higher than the specified threshold. Can be use in combination with ‘–rate’, ‘–info’ and ‘–genomic’.
##' @param rate extract markers with a completion rate equal or higher to the specified threshold. Can be use in combination with ‘--maf’, ‘--info’ and ‘--genomic’.
##' @param info Extract markers with an information equal or higher to the specified threshold. Can be use in combination with ‘--maf’, ‘--rate’ and ‘--genomic’.
  
#' @return  The extracted imputed files using genipe.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples  


# impute2-extractor  --impute2 imputedChr2.impute2   --out imputedChr2   --format  bed   --extract  imputedChr2.marker 

extractByGenipe <- function(inputImpute2, inputMAP, outputPrefix, format, prob){ 
        
        ## extract all the markers from the input map file 
        tmpMarkerFile = 'markers.txt'
        system( paste0("awk '{print $2}' ", inputMAP, " > ", tmpMarkerFile) )

        system( paste0("impute2-extractor ", 
        " --impute2 ", inputImpute2, " \ ", 
        # " --index ", " \ ",  
        " --out ", outputPrefix, " \ ",  
        " --format  ", format, " \ ",  
        " --extract  ", tmpMarkerFile, " \ "
        ) )
        system(paste0('rm ', tmpMarkerFile))    
}
 

# extractByGenipeGenomicRegion <- function(inputImpute2, outputPrefix, format, prob, genomicRange, maf, rate, info){             

#             system( paste0("impute2-extractor ", 
#             " --impute2 ", inputImpute2, " \ ",
#             " --index ", " \ ",  
#             " --out ", outputPrefix, " \ ",  
#             " --format  ", format, " \ ",  
#             " --genomic  ", genomicRange, " \ ",
#             " --maf  ", maf, " \ ",
#             " --rate  ", rate, " \ ",
#             " --info  ", info, " \ "
#             ) )    
# }



 
 

 
