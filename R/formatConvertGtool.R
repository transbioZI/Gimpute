

  
##########################################################################
########################################################################## formatConvertGtool.R  
#' Convert IMPUTE2 format files into PLINK format
#'
#' @description
#' Convert all chunks of IMPUTE2 format files into PLINK format using GTOOL.

#' @param gtool an executable GTOOL program in either the current working directory or somewhere in the command path.
#' @param chrs  specifiy the chromosomes for conversion 
#' @param prefixChunk  the prefix of the chunk files for each chromosome, along with the proper location directory.
#' @param phaseDIR the directory where pre-phased files are located.
#' @param imputedDIR the directory where the imputated files are located.
#' @param prefix4plinkEachChr the prefix of the input IMPUTE2 files and also the output PLINK files for each chunk.
#' @param suffix4imputed the suffix of the IMPUTE2 format file that stores the imputed value.
#' @param postImputeDIR the directory where converted PLINK files will be located. 
#' @param nCore the number of cores used for computation.  
 

#' @return  The converted PLINK format files for each chunk.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
  

formatConvertGtool <- function(gtool, chrs, prefixChunk, phaseDIR, imputedDIR, prefix4plinkEachChr, suffix4imputed, postImputeDIR, nCore){

		
		for (i in chrs){ 

			chunkfn = paste0(prefixChunk, i, ".txt")
			chunks = read.table(chunkfn, sep=" ") 
			chunklist = as.list(1:nrow(chunks))

			mclapply(chunklist, function(j){ 
				chunkSTART = chunks[j,1]
				chunkEND   = chunks[j,2] 
  				# INPUT data files
				SAM_FILE = paste0(phaseDIR, 'chr', i, '.sample')  
				GEN_FILE = paste0(imputedDIR, prefix4plinkEachChr, i, '.pos', chunkSTART, '-', chunkEND, suffix4imputed) 
				# output PLINK files
				PED_FILE = paste0(postImputeDIR, prefix4plinkEachChr, i, '.pos', chunkSTART, '-', chunkEND, '.ped') 
				MAP_FILE = paste0(postImputeDIR, prefix4plinkEachChr, i, '.pos', chunkSTART, '-', chunkEND, '.map') 
				
				## converting chunk-wise
				system( paste0(gtool, " -G ",
								" --g ", GEN_FILE, " \ ", 
								" --s ", SAM_FILE, " \ ",  
								"--phenotype plink_pheno \ ", 
								"--chr ", i, " \ ", 
								"--ped ", PED_FILE, " \ ", 
								"--map ", MAP_FILE, " \ ", 
								"--snp ") )  
			}, mc.cores=nCore)
		} 
} 
  
	 
   