
  


##########################################################################
########################################################################## mergeImputeData.R
#' Merge chunk-wise PLINK files 
#'
#' @description
#' Merge all chunk-wise PLINK files into chromosome-wise PLINK files then assemble into one PLINK file set. 
#' Create a file containing a list chunk ped and map file names.
#' At last, combine all chrs (combine the first 23 chrs; then Xpar).

#' @param plink an executable PLINK program in either the current working directory or somewhere in the command path.
#' @param chrs  specifiy the chromosomes to be merged. 
#' @param prefix4plinkEachChr the prefix of the input chunk-wise IMPUTE2 files. 
#' @param prefix4imputedPlink  the prefix of the final imputed PLINK format files. 
#' @param nCore the number of cores used for computation.  

#' @return  The merged PLINK format files.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

 
 
mergeImputeData <- function(plink, chrs, prefix4plinkEachChr, prefix4imputedPlink, nCore){ 

 
	## firstly, only consider chromosomes from 1:23; as Xpar chrs are slightly different for processing.
	pureAutoChrs = setdiff(chrs, c("X_PAR1", "X_PAR2")) 
	chrslist = as.list(pureAutoChrs)   
	mclapply(chrslist, function(i){
		
		pedFile_chr = system(paste0("ls ", prefix4plinkEachChr, i, ".*.ped"), intern=TRUE)
		mapFile_chr = system(paste0("ls ", prefix4plinkEachChr, i, ".*.map"), intern=TRUE)	
		pedmap_chr = paste0(pedFile_chr, " ", mapFile_chr)
		fA = gsub(".ped", "", pedFile_chr[1])
		pedmap_tobeMerged = pedmap_chr[-1]
		filesetname = paste0("fileset_chr", i, ".txt")
		write.table(pedmap_tobeMerged, file=filesetname, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	    system( paste0(plink, " --file ", fA, " --merge-list ", filesetname, " --make-bed --out gwasImputed_chr", i) ) 
	    # system( paste0('rm ', filesetname))
	}, mc.cores=nCore)


	############################################################################## combine chrX_PAR and convert into chr25
	##  

	if (is.element(c("X_PAR1"), chrs) | is.element(c("X_PAR2"), chrs) ){  ## tobeImproved     the condition has length > 1 and only the first element will be used

		pedFile_chr = system(paste0("ls ", prefix4plinkEachChr, "X_PAR*.ped"), intern=TRUE)
		mapFile_chr = system(paste0("ls ", prefix4plinkEachChr, "X_PAR*.map"), intern=TRUE)	
		pedmap_chr = paste0(pedFile_chr, " ", mapFile_chr)
		fA = gsub(".ped", "", pedFile_chr[1])
		pedmap_tobeMerged = pedmap_chr[-1]
		filesetname = paste0("fileset_chr25.txt")
		write.table(pedmap_tobeMerged, file=filesetname, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		 
		arg <- paste0(plink, " --file ", fA, " --merge-list ", filesetname, " --allow-extra-chr --make-bed --out gwasImputed_oldchr25")
		system(arg) 
		## update chr code for XPAR --> 25
		bim = read.table("gwasImputed_oldchr25.bim", stringsAsFactors=F)
		updateSNPchr = cbind(bim[,2], rep(25, length=nrow(bim))) 
		write.table(updateSNPchr, file="gwasImputed_newchr25.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
		system( paste0(plink, " --bfile gwasImputed_oldchr25 --allow-extra-chr --update-chr gwasImputed_newchr25.txt 2 1 --make-bed --out gwasImputed_chr25") )  
		system( 'rm gwasImputed_oldchr25.* gwasImputed_newchr25.txt ')
		# system( paste0('rm ', filesetname))
	}	 
	

	############################################################################## combine all bed files
	bedFile_chr = system(paste0("ls gwasImputed_chr*.bed"), intern=TRUE)
	bimFile_chr = system(paste0("ls gwasImputed_chr*.bim"), intern=TRUE)	
	famFile_chr = system(paste0("ls gwasImputed_chr*.fam"), intern=TRUE)	
	bfile_chr = paste0(bedFile_chr, " ", bimFile_chr, " ", famFile_chr)
	fA = paste0(gsub(".bed", "", bedFile_chr[1]))
	tobeMerged = bfile_chr[-1]
	mergefilesetname = paste0("mergeGwasImputed.txt")
	write.table(tobeMerged, file=mergefilesetname, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	arg <- paste0(plink, " --bfile ", fA, " --merge-list ", mergefilesetname, " --make-bed --out ", prefix4imputedPlink)
	system(arg)

	# system( paste0('rm ', mergefilesetname)) 
} 

 
