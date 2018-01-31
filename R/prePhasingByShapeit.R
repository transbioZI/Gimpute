

	

##########################################################################
########################################################################## prePhasingByShapeit.R
#' Prephasing genotypes using SHAPEIT
#'
#' @description
#' Perform prephasing for study genotypes by SHAPEIT for the autosomal and sex chromosome haplotypes using a reference panel (pre-set).
#' If ChrX is available then it is done differently by passing the flag --chrX to SHAPEIT.

#' @param shapeit an executable SHAPEIT program in either the current working directory or somewhere in the command path.
#' @param chrs  specifiy the chromosomes for phasing.
#' @param dataDIR the directory where genotype PLINK files are located.
#' @param prefix4plinkEachChr the prefix of PLINK files for each chromosome.
#' @param impRefDIR the directory where the imputation reference files are located.
#' @param phaseDIR the directory where resulting pre-phased files will be located.
#' @param nThread the number of threads used for computation.
#' @param effectiveSize this parameter controls the effective population size.
#' @param nCore the number of cores used for computation. This can be tuned along with nThread.
 

#' @return  The pre-phased haplotypes for given chromosomes.  
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 

prePhasingByShapeit <- function(shapeit, chrs, dataDIR, prefix4plinkEachChr, impRefDIR, phaseDIR, nThread, effectiveSize, nCore){

	chrslist = as.list(chrs)
	mclapply(chrslist, function(i){
 
			# GWAS data files in PLINK binary format
			GWASDATA_BED = paste0(dataDIR, prefix4plinkEachChr, i, '.bed ') 
			GWASDATA_BIM = paste0(dataDIR, prefix4plinkEachChr, i, '.bim ')
			GWASDATA_FAM = paste0(dataDIR, prefix4plinkEachChr, i, '.fam ')
			# reference data files
			GENMAP_FILE = paste0(impRefDIR, 'genetic_map_chr', i, '_combined_b37.txt ')
			HAPS_FILE = paste0(impRefDIR, 'ALL_1000G_phase1integrated_v3_chr', i, '_impute_macGT1.hap.gz ') 
			LEGEND_FILE = paste0(impRefDIR, 'ALL_1000G_phase1integrated_v3_chr', i, '_impute_macGT1.legend.gz ')
			SAMPLE_FILE = paste0(impRefDIR, 'ALL_1000G_phase1integrated_v3.sample ')

			# main output file
			OUTPUT_HAPS = paste0(phaseDIR, 'chr', i, '.haps ')     
			OUTPUT_SAMPLE = paste0(phaseDIR, 'chr', i, '.sample ')     
			OUTPUT_LOG = paste0(phaseDIR, 'chr', i, '.log ')    

			if (i!=23){  ## prePhasing for the autosome
				system(paste0(shapeit, 
				" --input-bed ", GWASDATA_BED, GWASDATA_BIM, GWASDATA_FAM, " \ ", 
				" --input-map ", GENMAP_FILE, " \ ",  
				"--input-ref ", HAPS_FILE, LEGEND_FILE, SAMPLE_FILE, " \ ", 
				"--thread ", nThread, " \ ", 
				"--effective-size ", effectiveSize, " \ ", 
				"--output-max ", OUTPUT_HAPS, OUTPUT_SAMPLE, " \ ", 
				"--output-log ", OUTPUT_LOG) )
			} else if (i==23){
				system(paste0(shapeit, 
				" --input-bed ", GWASDATA_BED, GWASDATA_BIM, GWASDATA_FAM, " \ ", 
				" --input-map ", GENMAP_FILE, " \ ", 
				" --chrX \ ", 
				"--input-ref ", HAPS_FILE, LEGEND_FILE, SAMPLE_FILE, " \ ", 
				"--thread ", nThread, " \ ", 
				"--effective-size ", effectiveSize, " \ ", 
				"--output-max ", OUTPUT_HAPS, OUTPUT_SAMPLE, " \ ", 
				"--output-log ", OUTPUT_LOG) ) 
			}

	}, mc.cores=nCore)	 
} 

 

