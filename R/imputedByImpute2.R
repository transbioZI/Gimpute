

##########################################################################
########################################################################## imputedByImpute2.R
#' Impute genotypes using IMPUTE2
#'
#' @description
#' Perform imputation by IMPUTE2 for the autosomal and sex chromosome prephased known haplotypes with a reference panel.

#' @param impute2 an executable IMPUTE2 program in either the current working directory or somewhere in the command path.
#' @param chrs  specifiy the chromosomes for imputation.
#' @param prefixChunk  the prefix of the chunk files for each chromosome, along with the proper location directory.
#' @param phaseDIR the directory where prephased haplotypes are located.
#' @param impRefDIR the directory where the imputation reference files are located.
#' @param imputedDIR the directory where imputed files will be located.
#' @param prefix4plinkEachChr the prefix of IMPUTE2 files for each chunk.
#' @param nCore the number of cores used for computation.
#' @param effectiveSize this parameter controls the effective population size.
#' @param XPAR a logical value indicating whether --chrX flag should be passed for prephasing using SHAPEIT.
#' --chrX flag, specifically for chrX imputation'
#' @return  The imputed files for all chunks from given chromosomes.  
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 

imputedByImpute2 <- function(impute2, chrs, prefixChunk, phaseDIR, impRefDIR, imputedDIR, prefix4plinkEachChr, nCore, effectiveSize){ 
 
 
	for (i in chrs){ 	

		chunkfn = paste0(prefixChunk, i, ".txt")
		chunks = read.table(chunkfn, sep=" ")
		# dim(chunks)  
		chunklist = as.list(1:nrow(chunks))
		mclapply(chunklist, function(j){

			chunkSTART = chunks[j,1]
			chunkEND   = chunks[j,2] 
			# Input: haplotypes from SHAPEIT phasing (method B)
			GWAS_HAPS_FILE = paste0(phaseDIR, 'chr', i, '.haps ') 
			GWAS_SAMP_FILE = paste0(phaseDIR, 'chr', i, '.sample ') 
			# reference data files
			## For other reference panels you want to modify the following setting  
			GENMAP_FILE = paste0(impRefDIR, 'genetic_map_chr', i, '_combined_b37.txt ')
			HAPS_FILE = paste0(impRefDIR, 'ALL_1000G_phase1integrated_v3_chr', i, '_impute_macGT1.hap.gz ') 
			LEGEND_FILE = paste0(impRefDIR, 'ALL_1000G_phase1integrated_v3_chr', i, '_impute_macGT1.legend.gz ')
 
			# main output file    
			OUTPUT_FILE = paste0(imputedDIR, prefix4plinkEachChr, i, '.pos', chunkSTART, '-', chunkEND, '.impute2 ')   

			################## impute genotypes from GWAS haplotypes 
			if ( is.element(i, c(1:22)) ){ 
				## impute for the autosomes
				system(paste0(impute2, 
				" -iter 30  \ ", 
				" -burnin 10  \ ", 
				" -k_hap 500  \ ", 
				" -use_prephased_g  \ ",  
				" -m ", GENMAP_FILE, " \ ",  
				" -h ", HAPS_FILE, " \ ", 
				" -l ", LEGEND_FILE, " \ ", 
				" -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
				" -Ne ", effectiveSize, " \ ", 
				" -int ", chunkSTART, " ", chunkEND, " \ ", 
				" -buffer 1000  \ ",
				" -o ", OUTPUT_FILE, " \ ", 
				" -allow_large_regions \ ",
				" -seed 367946 \ " ))
			} else if ( is.element(i, c("X_PAR1", "X_PAR2")) ){  
				## impute for chrX PAR >> with an additional flag: --Xpar.
				system(paste0(impute2, 
				" -iter 30  \ ", 
				" -burnin 10  \ ", 
				" -k_hap 500  \ ", 
				" -use_prephased_g  \ ", 
				" -Xpar \ ",     #################
				" -m ", GENMAP_FILE, " \ ",  
				" -h ", HAPS_FILE, " \ ", 
				" -l ", LEGEND_FILE, " \ ", 
				" -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
				" -Ne ", effectiveSize, " \ ", 
				" -int ", chunkSTART, " ", chunkEND, " \ ", 
				" -buffer 1000  \ ",
				" -o ", OUTPUT_FILE, " \ ", 
				" -allow_large_regions \ ",
				" -seed 367946 \ " ))
			} else if ( i==23 ){ 
				## impute for chrX  >> with an additional flag: --chrX, and sample_known_haps_g
				system(paste0(impute2, 
				" -iter 30  \ ", 
				" -burnin 10  \ ", 
				" -k_hap 500  \ ", 
				" -use_prephased_g  \ ", 
				" -chrX \ ",   #################     
				" -m ", GENMAP_FILE, " \ ",  
				" -h ", HAPS_FILE, " \ ", 
				" -l ", LEGEND_FILE, " \ ", 
				" -known_haps_g ", GWAS_HAPS_FILE, " \ ", 
				" -sample_known_haps_g ", GWAS_SAMP_FILE, " \ ",  #################
				" -Ne ", effectiveSize, " \ ", 
				" -int ", chunkSTART, " ", chunkEND, " \ ", 
				" -buffer 1000  \ ",
				" -o ", OUTPUT_FILE, " \ ", 
				" -allow_large_regions \ ",
				" -seed 367946 \ " ))
			}

		}, mc.cores=nCore)  
 	} 
} 		
 
 