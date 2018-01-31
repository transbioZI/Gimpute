

##########################################################################
########################################################################## mergeByGenipe.R
#' Merge imputed files using Genipe

#' @description
#' Concatenate IMPUTE2 output files and retrieve some statistics. 
#' This is automatically called by the main genipe pipeline to merge IMPUTE2 files generated for all the genomic segments. 
#' For details, see Genipe IMPUTE2 merger options.

#' @param chr  specifiy the chromosome segment to be merged, on which the imputation was made.
#' @param inputImpute2 the output from IMPUTE2.
#' @param probability the probability threshold for no calls. [<0.9]
#' @param completionRate the completion rate threshold for site exclusion. [<0.98]
#' @param info the measure of the observed statistical information associated with the allele frequency estimate threshold for site exclusion. (<0.00)
#' @param outputPrefix the prefix for the imputed output files. 
 
#' @return  The merged imputed files using genipe.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
  

mergeByGenipe <- function(inputImpute2, chr, probability, completionRate, info, outputPrefix){ 
 
        
            system( paste0("impute2-merger ", 
            " --impute2 ", inputImpute2, " \ ",
            " --chr ", chr, "\ ",   
            " --probability  ", probability, " \ ", 
            " --completion ", completionRate, " \ ", 
            " --info  ", info, " \ ",  
            " --prefix ", outputPrefix ) )    
}
 

