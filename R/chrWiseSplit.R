
 

##########################################################################
########################################################################## chrWiseSplit.R


#' Split genome-wide genotyping data into separate files by chromosome
#'
#' @description
#' Split the whole genome genotyping data chromosome-wise; allow parallel computating for all chromosomes.
#' if chromosome 25 is also available, further split chr25 (PAR or Chr_XY) into PAR1 and PAR2 according to the genomic coordination GRCh37 

#' from https://en.wikipedia.org/wiki/Pseudoautosomal_region.
#' The locations of the PARs within GRCh37 are:  
#' PAR1	X	60001	2699520 
#' PAR2	X	154931044	155260560 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input plink files before splitting.
#' @param outputPrefix the prefix of the output plink files after splitting and the chromosome number will be appended separately.
#' @param chrX_PAR1suffix  if chromosome 25 is available and with PAR1, then generate the suffix with X_PAR1 for chrX_PAR1 
#' @param chrX_PAR2suffix  if chromosome 25 is available and with PAR2, then generate the suffix with X_PAR2 for chrX_PAR2 

#' @return  The output plink files for each chromosome and possibly also the suffix of chrX_PAR.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 



chrWiseSplit <- function(plink, inputPrefix, chrX_PAR1suffix, chrX_PAR2suffix){ 


   
    ## check which chromosomes are available to be splitted from the .bim file
	bim = read.table(paste0(inputPrefix, '.bim'), stringsAsFactors=F)
	chrs = as.integer(names(table(bim[,1])))
 
	chrslist = as.list(chrs)
	mclapply(chrslist, function(i){
		cmd = paste0(plink, " --bfile ", inputPrefix, " --chr ", i, " --make-bed --out ", inputPrefix, i)  
		system(cmd)
	}, mc.cores=length(chrs))

	## if chromosome 25 is also available then re-arrange it
 	if (is.element(25, chrs)){  

 		print('PAR is available in chrX!') 
		bim25 = read.table(paste0(inputPrefix, "25.bim"), stringsAsFactors=F) 
		pos4PAR1= c(60001, 2699520) 
		## first check for PAR1 and afterwards for PAR2
		if ( length(which(bim25[,4]<=pos4PAR1[2]))!=0 ){ 
	   
	   		print('PAR1 is available in chrX!')
	   		bimPos4par1 = which(bim25[,4]<=pos4PAR1[2])
			rs4PAR1 = bim25[bimPos4par1,2]
			write.table(rs4PAR1, file="rs4PAR1.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
			system( paste0(plink, " --bfile ", inputPrefix, "25 --extract rs4PAR1.txt --make-bed --out ", inputPrefix, chrX_PAR1suffix) )
			par1 = TRUE  
			## check for PAR2, if any SNPs out of PAR1, then PAR2 also available 
			if (length(bimPos4par1)<nrow(bim25)){
				print('PAR2 is available in chrX!') 
				## just to exclude rs4PAR1.txt
				system( paste0(plink, " --bfile ", inputPrefix, "25 --exclude rs4PAR1.txt --make-bed --out ", inputPrefix, chrX_PAR2suffix) )
				par2 = TRUE 
			} else {print('PAR2 is NOT available in chrX!') }

		} else { 
			print('PAR2 is available in chrX! But NOT PAR1, all chr25 on PAR2')
			par1 = FALSE
			par2 = TRUE
			system( paste0("cp ", inputPrefix, "25.bed ", inputPrefix, chrX_PAR2suffix, ".bed") )
			system( paste0("cp ", inputPrefix, "25.bim ", inputPrefix, chrX_PAR2suffix, ".bim") )
			system( paste0("cp ", inputPrefix, "25.fam ", inputPrefix, chrX_PAR2suffix, ".fam") )
		} 
	} else {  
		print('PAR is NOT available in chrX!') 
		par1 = FALSE
		par2 = FALSE
	} 

	return(par=list(par1, par2)) 
}
 

