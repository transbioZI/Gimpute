
#' Remove heterozygous SNPs in haploid male chromosome X
#'
#' @description
#' In principle, one has to remove all heterozygous SNPs of chromosome X in males. However, too many SNPs might be removed in some data sets. 
#' We want to accept a small percentage of such SNPs in the data set so that we do not lose too much SNPs. 

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param hhCutOff the cutoff for removing male haploid heterozygous SNPs on the chromosome X.
#' @param outputPrefix the prefix of the output PLINK format files.
#' @param outputFile_SNPhhFreqAll the output pure text file that stores all heterozygous SNPs with their frequency (the number of males for this SNP), if any. 
#' @param outputFile_SNPhhFreqRetained the output pure text file that stores retained heterozygous SNPs with their frequency (the number of males for this SNP), if any. 

#' @return  The output PLINK format files and two pure text files (if any) with heterozygous SNPs and their respective frequency. 
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 
   
# Output file with the number of instances with heterozygous alleles for
# each SNP of the chromosome 23 before SNP removal (each line contains
# a SNP name and the respective number, lines are sorted descending by
# number):
# outputSNPs: 
# Output file with the number of instances with heterozygous alleles for
# each SNP of the chromosome 23 after SNP removal:
# # Determine for each SNP of the chromosome 23 from the genotype data
# the number of male instances which have the value one as the minor
# allele count for that SNP and remove all SNPs which number is higher
# than 0.5 % of the number of male instances.

 
#  A haploid heterozygous is a male genotype that is heterozygous, 
#  which is an error given the haploid nature of the male X chromosome.
removedSnpHetX <- function(plink, inputPrefix, hhCutOff, outputPrefix, outputFile_SNPhhFreqAll, outputFile_SNPhhFreqRetained){

		

			## just to get .hh file and .fam file 
			system( paste0(plink, " --bfile ", inputPrefix, " --chr 23 --filter-males --make-bed --out male23nonPAR" ))

			if ( is.na(file.size("male23nonPAR.hh")) ){ 
				## copy/rename all snp info updated plink files
				system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
				system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
				system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") ) 
	  
		    } else { 	

				# .hh (heterozygous haploid and nonmale Y chromosome call list)
				# Produced automatically when the input data contains heterozygous calls 
				# where they should not be possible (haploid chromosomes, male X/Y),
				# or there are nonmissing calls for nonmales on the Y chromosome.

				# A text file with one line per error (sorted primarily by variant ID, secondarily by sample ID) with the following three fields:
				# Family ID  Within-family ID Variant ID
				hh = read.table("male23nonPAR.hh", stringsAsFactors=F)
				fam = read.table("male23nonPAR.fam", stringsAsFactors=F)

				hetSNPsFreq = table(hh[,3])
				# hetSNPFreqFreq = table(hetSNPs)
				# head(sort(hetSNPsFreq, decreasing=TRUE),11)

				cutoff4removeHetSNP = nrow(fam)*hhCutOff
				mostFakeSNPs = hetSNPsFreq[which(hetSNPsFreq>= cutoff4removeHetSNP)] 
				mostFakeSNPs = names(mostFakeSNPs)
				str(mostFakeSNPs)
				write.table(mostFakeSNPs, file="mostFakeSNPs.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
				## remove these fake SNPs
				system( paste0(plink, " --bfile ", inputPrefix, " --exclude mostFakeSNPs.txt --make-bed --out ", outputPrefix) )
				## remove unwanted files
				system( paste0("rm mostFakeSNPs.txt") )

				## generate hetSNPsFreq in .txt file 
				hetSNPsWithInstNumber = as.data.frame(hetSNPsFreq, stringsAsFactors=FALSE)
				hetSNPsWithInstNumber = hetSNPsWithInstNumber[order(hetSNPsWithInstNumber[,2], decreasing=T),] 
				## all heterozygous SNPs 
				write.table(hetSNPsWithInstNumber, file=outputFile_SNPhhFreqAll, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 

				## remaining heterozygous SNPs 
				hetSNPsWithInstNumberSub = hetSNPsWithInstNumber[!is.element(hetSNPsWithInstNumber[,1], mostFakeSNPs), ]
				write.table(hetSNPsWithInstNumberSub, file=outputFile_SNPhhFreqRetained, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
		    } 

			system( paste0("rm male23nonPAR.*") )
			system( paste0("rm ", inputPrefix,".*") )
}	
