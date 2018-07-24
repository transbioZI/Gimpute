

#' Post imputation quality control
#'
#' @description
#' Perform quality control and data management after imputation. 

#' @param inputPrefix the prefix of the final imputed PLINK binary files.
#' @param out2 the prefix of the final imputed PLINK binary files with the index. 
#' @param out3 the prefix of well imputed PLINK binary files with the index. 
#' @param out4 the prefix of well imputed PLINK binary files after removing any 
#' SNPs with the same positions (if any), the index is also appended.
#' @param out5 the prefix of well imputed PLINK binary files with the index 
#' after adding previously identified monomorphic SNPs if any.
#' @param out6 the prefix of final well imputed PLINK binary files with the index. 
#' @param infoScore the cutoff of filtering imputation quality score for 
#' each variant. The default value is 0.6. 
#' @param outputInfoFile the output file of impute2 info scores consisting of 
#' two columns: all imputed SNPs and their info scores.   
#' @param inputPrefix4aligned2impRef the prefix of the output PLINK binary 
#' files after removing SNPs whose alleles are not in the imputation reference,
#' taking their genomic positions into account.
#' @param missCutoff  the cutoff of the least number of instances for 
#' a SNP that is not missing. The default is 20.
#' @param outRemovedSNPfile the output file of SNPs with pre-defined 
#' missing values that are removed.
#' @param outRetainSNPfile the output file of SNPs that are retained. 


postImpQC <- function(inputPrefix, out1, out2, out3, out4, out5, out6,
                      outputInfoFile, infoScore=0.6, 
                      inputPrefix4aligned2impRef, missCutoff=20, 
                      outRemovedSNPfile, outRetainSNPfile, referencePanel){

    ## step 2 
    ## Final imputed results, including bad imputed genotypes. 
    renamePlinkBFile(inputPrefix=inputPrefix, outputPrefix=out2, action="move")
    system(paste0("rm ", inputPrefix, "*")) ## remove redundant 
    ## step 3
    ## Filtered imputed data set; Remove imputed SNPs with (info < 0.6), 
    ## only retain "Good" SNPs.    
    .filterImputeData2(plink, outputInfoFile, infoScore, 
                       inputPrefix=out2, outputPrefix=out3)  

    ## step 4 
    ## Remove previously identified monomorphic SNPs in the imputed dataset. 
    ## Note that snps with same genomic position but can have different snp name.
    ## if no monomorphic SNPs:
    if (file.size(paste0(outputMonoSNPfile)) == 0 ){
        renamePlinkBFile(inputPrefix=out3, outputPrefix=out4, action="copy")   
    } else { 
        ## extract PLINK files contain only monomorphic SNPs from 
        ## the original aligned (lifted and QC-ed) data set.
        system(paste0(plink, " --bfile ", inputPrefix4aligned2impRef, 
               " --extract ", outputMonoSNPfile, " --make-bed --out ", 
               inputPrefix4aligned2impRef, "Tmp")) 
        bim1 <- read.table(paste0(inputPrefix4aligned2impRef, "Tmp.bim"), 
                           stringsAsFactors=F)
        system(paste0("awk '{print $1, $2, $4}' ", 
               out3, ".bim > tmpFilterImp.txt"))
        bim2 <- read.table("tmpFilterImp.txt", stringsAsFactors=F) 
        colnames(bim1) <- c("chr", "rsID", "gd", "pos", "a0", "a1") 
        colnames(bim2) <- c("chr", "rsID", "pos")
        outputFile <- "tmp.txt"
        .snpSharedPos(inputFile1=bim1, inputFile2=bim2, outputFile, nCore=25) 
        system(paste0(plink, " --bfile ", out3, 
               " --exclude tmp.txt --make-bed --out ", out4))
        system("rm tmpFilterImp.txt tmp.txt")
    }
     
    ## step 5
    ## Add previously identified monomorphic SNPs in the imputed dataset. 
     ## if no monomorphic SNPs:
    if ( file.size(paste0(outputMonoSNPfile))==0 ){ 
            renamePlinkBFile(inputPrefix=out3, outputPrefix=out5, action="copy")  
    } else { 
        ## merge both datasets
        system(paste0(plink, " --bfile ", out4, " --bmerge ", 
               inputPrefix4aligned2impRef, ".bed ", 
               inputPrefix4aligned2impRef, ".bim ", 
               inputPrefix4aligned2impRef, ".fam ", "--make-bed --out ", out5))
    }  

    ## step 6
    ## Remove SNPs which have a non missing value for less then 20 instances. 
    removedSnpMissPostImp(plink, inputPrefix=out5, missCutoff, 
                          outRemovedSNPfile, outRetainSNPfile, outputPrefix=out6)
    ## special case: 
    ## 1.) if the input ref panel is from phase3. 
    ## (2) Impute4 modified imputed output rsID if the original SNPs are not given.
    ## e.g. rs456706:11001104:C:T << rsID  rs456706
    if (referencePanel == "1000Gphase3"){
        outRename <- "4_6_renameSNP"
        renamePlinkBFile(inputPrefix=out6, outputPrefix=outRename, action="move") 
        snp2renameV0 <- read.table(outRetainSNPfile, stringsAsFactors=FALSE)
        str(snp2renameV0)
        whDiff <- grep(":", snp2renameV0[,1])
        str(whDiff)
        if (length(whDiff) > 0){         
            snp2renameV1 <- snp2renameV0[whDiff,]
            tmp = strsplit(snp2renameV1, ":", fixed = TRUE)
            snpNew <- unlist(lapply(tmp, function(s){s[[1]]})) 
            snpNewOld <- cbind(snpNew, snp2renameV0[whDiff,1])
            str(snpNewOld)
            ## rename the plink files. 
            snpRenameFile <- "4_6_renameSNP.txt"
            write.table(snpNewOld, file=snpRenameFile, quote=FALSE, 
                        row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
            system(paste0(plink, " --bfile ", outRename, " --update-name ", 
                   snpRenameFile, " 1 2 --make-bed --out ", outRename, "tmp"))
            renamePlinkBFile(inputPrefix=paste0(outRename, "tmp"), 
                             outputPrefix=out6, action="move") 
            system(paste0("rm ", outRename, "*")) 
        }     
    }

}
