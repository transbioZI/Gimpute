



##########################################################################
########################################################################## chunk4eachChr.R
#' Chunk each chromosome into multiple segments
#'
#' @description
#' Chunk each chromosome into multiple segments by a predefined window size.

#' @param inputPrefix the prefix of the input plink files for each chromosome.
#' @param outputPrefix the prefix of the output pure text files that keep all the chunks for each chromosome.
#' @param chrs  specifiy the chromosomes for chunking.
#' @param windowSize  the window size of each chunk.

#' @return  The output pure text files that keep all the chunks for each chromosome.
#' @export 

#' @author Junfang Chen <junfang.chen3@gmail.com> 
#' @examples 
 

chunk4eachChr <- function(inputPrefix, outputPrefix, chrs, windowSize){  

		for (i in chrs){ 

			bimfilename = paste0(inputPrefix, i, ".bim")
			bimdata = read.table(file=bimfilename, sep="\t", stringsAsFactors=F)
			position = bimdata[,4]
			posStart = head(position,1)
			posEnd = tail(position,1)
			chunkStart = seq(posStart, posEnd, windowSize)
			chunkEnd = chunkStart + windowSize -1
			chunkMatrix = cbind(chunkStart, chunkEnd)

			## positions are only within a chunk
			if (nrow(chunkMatrix) == 1){
				chunks = chunkMatrix
			} else {  
				##  it may happen that only a few SNPs from the last chunk; but if the last chunk is large, then specify -allow_large_regions 
				chunks = head(chunkMatrix, -1) ## merge last-second to last  
				chunks[nrow(chunks), 2] <- posEnd
			}
			

			## check if any chunk with NO snps within it   
			SNPcountsPerChunk <- c() 
			for (j in 1:nrow(chunks)){
				chunkbottom = chunks[j,1]
				chunkup = chunks[j,2]
				tmp = which(position >= chunkbottom & position <= chunkup)  ## which fall within chunk
				tmp = length(tmp)
				SNPcountsPerChunk <- c(SNPcountsPerChunk, tmp) 
			} 


		 	wh0 = which(SNPcountsPerChunk==0)
		 	print(paste0("chr",i))
		 	print(wh0)
		 	chunkLength = nrow(chunks) - length(wh0)
		 	## remove such chunks if wh0  

		 	if (length(wh0)!=0){ chunks = chunks[-wh0,]  }
		 	print(nrow(chunks) == chunkLength)

			chunkfilename = paste0(outputPrefix, i, ".txt")
			write.table(chunks, file=chunkfilename, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		}
}


