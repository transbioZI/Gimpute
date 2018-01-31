
 

##########################################   
##########################################
#' Population outlier detection via PCA 
#'
#' @description
#' Perform principle component analysis (PCA) on quality controlled genotype data set in order to detect population outliers, 
#' and plot the PCs for visualization. Note that, only autosomal genotypes are used for computing principle components. 
  

#' @param gcta an executable GCTA program in either the current working directory or somewhere in the command path.
#' @param inputPrefix the prefix of the input PLINK format files.
#' @param outputPC4subjFile the pure text file that stores all the subject IDs after QC and their corresponding eigenvalues of the first two principle components.
#' @param outputPCplotFile the plot file for visualizing the first two principle components of all subjects.

#' @return  The output pure text file and plot file for storing first two principle components of study subjects.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 

## http://cnsgenomics.com/software/gcta/#Download

plotPCA4plink <- function(gcta, inputPrefix, outputPC4subjFile, outputPCplotFile){ 

	autosomefn = paste0(inputPrefix, "Autosome")
	system( paste0(gcta, " --bfile ", inputPrefix, " --make-grm --autosome --out ", autosomefn, " --thread-num 30") )
	system( paste0(gcta, " --grm ", autosomefn, " --pca 20 --out ", autosomefn, " --thread-num 30") )

	eigen = read.table(file=paste0(autosomefn,".eigenvec"), stringsAsFactors=F)
	pcs = eigen[,1:4] ## first two PCs
	write.table(pcs, outputPC4subjFile, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
	pcWithGroup = cbind(pcs, stringsAsFactors=F)

	png(outputPCplotFile, width=8, height=6, units="in", res=800)
	print( xyplot(pcWithGroup[,4] ~ pcWithGroup[,3], data=pcWithGroup, 
	       auto.key=list(space="right"),  
	       jitter.x=TRUE, jitter.y=TRUE, xlab="PC1", ylab="PC2") )
	dev.off()
	## remove unwanted files
	system( paste0("rm ", autosomefn, ".*") )
}


 