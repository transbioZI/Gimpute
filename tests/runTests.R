
# File   : runTests.R
# Author : Junfang Chen
# Version0: 28 Jun 2016
# VersionX: 24 Jul 2018
   


library(Gimpute)

## Start an R session in a directory where you'd like to generate the data.
system("mkdir 1-genoUpdate")
system("mkdir 2-genoQC")
system("mkdir 3-checkAlign")
system("mkdir 4-imputation")
system("mkdir 5-reductAndExpand")
system("mkdir 6-finalResults")

## Define the directory where you place the imputation reference files 
## Reference panel 1
referencePanel <- "1000Gphase1v3_macGT1" ## indicator
impRefDIR1kGp1v3 <- "/data/noether/dataRawReadOnly/reference/1000GP_phase1v3/"
impRefDIR <- impRefDIR1kGp1v3

 # Reference panel 2
referencePanel <- "1000Gphase3" ## indicator 
impRefDIR1kGp3 <- "/data/noether/dataRawReadOnly/reference/1000GP_Phase3/"
impRefDIR <- paste0(impRefDIR1kGp3, "1000GP_Phase3/")


## Genotyping chip annotation file 
chipAnnoFile <- system.file("extdata", "coriellAffyChip.txt", 
                                 package="Gimpute")
## Self-defined configuration files
removedSampIDFile <- system.file("extdata", "excludedSampIDs.txt", 
                                 package="Gimpute")
excludedProbeIdsFile <- system.file("extdata", "excludedProbeIDs.txt", 
                                    package="Gimpute")

## Define required tools
plink <- "/home/junfang.chen/Gimpute/tools/plink"
gcta <- "/home/junfang.chen/Gimpute/tools/gcta64" 
shapeit <- "/home/junfang.chen/Gimpute/tools/shapeit"
impute2 <- "/home/junfang.chen/Gimpute/tools/impute2"
gtool <- "/home/junfang.chen/Gimpute/tools/gtool"
## Gimpute has the following dependencies: 
impute4 <- "/home/junfang.chen/Gimpute/tools/impute4.1_r291.2"
qctool <- "/home/junfang.chen/Gimpute/tools/qctool"

library(lattice)
library(doParallel)
  

############################################################

## Run the following code, only if you have the above tools 
## and the imputation reference files, configuration files.

############################################################


  


  
############################################################
## code chunk number 1: SNP information update  
############################################################
runTimeList <- list()
t4genoUpdateTmp <- proc.time()   
## step 0
## Load PLINK binary files and additional files from Gimpute.
setwd("./1-genoUpdate/") 
bedFile <- system.file("extdata", "dataChr21.bed", package="Gimpute")
bimFile <- system.file("extdata", "dataChr21.bim", package="Gimpute") 
famFile <- system.file("extdata", "dataChr21.fam", package="Gimpute")
system(paste0("scp ", bedFile, " ."))   
system(paste0("scp ", bimFile, " ."))   
system(paste0("scp ", famFile, " ."))   

metadataFile <- system.file("extdata", "1_01_metaData.txt", package="Gimpute")
removedSampIDFile <- system.file("extdata", "excludedSampIDs.txt", 
                                 package="Gimpute")
excludedProbeIdsFile <- system.file("extdata", "excludedProbeIDs.txt", 
                                    package="Gimpute")
## Genotyping chip annotation file 
chipAnnoFile <- system.file("extdata", "coriellAffyChip.txt", 
                                 package="Gimpute")

## step 1  
system(paste0("scp ", metadataFile, " ."))  
 
 
############################################################ 
## pipeline function
inputPrefix <- "dataChr21"
ancestrySymbol <- "EUR"
outputPrefix <- "1_11_removedYMtSnp" 
metaDataFile <- "1_01_metaData.txt"
chipType <- "rsIDstudy"
updateGenoInfo(plink, inputPrefix, metaDataFile, removedSampIDFile,
               ancestrySymbol, excludedProbeIdsFile, chipAnnoFile,
               chipType, outputPrefix, keepInterFile=TRUE)
setwd("..")   

## runtime
t4genoUpdate <- proc.time() - t4genoUpdateTmp
print(t4genoUpdate)
runTimeList$t4genoUpdate <- t4genoUpdate

t4genoQCTmp <- proc.time() 

############################################################
### code chunk number 2: Quality Control
############################################################
## step 0
## copy the last output plink files from 1-genoUpdate
inputPrefix4QC <- "1_11_removedYMtSnp"
system(paste0("scp ./1-genoUpdate/", inputPrefix4QC, ".*", " ./2-genoQC/"))
setwd("./2-genoQC/")

## step 1
inputPrefix <- inputPrefix4QC
hhCutOff <- 0.005 ##  can be tuned
outputPrefix <- "2_01_removedSnpHetX" 
outputHetSNPfile <- "2_01_snpHHfreqAll.txt"
outputRetainSNPfile <- "2_01_snpHHfreqRetained.txt"
removedSnpHetX(plink, inputPrefix, hhCutOff, outputPrefix, 
               outputHetSNPfile, outputRetainSNPfile)


## step 2  2_02_removedHetXInst
inputPrefix <- "2_01_removedSnpHetX"
hhSubjCutOff <- 15
outputPrefix <- "2_02_removedInstHetX"
outputSubjHetFile <- "2_02_instHetXfreqAll.txt" 
outputRetainSubjectFile <- "2_02_instHetXfreqRetained.txt"  
outputHetSNPfile <- "2_02_snpHHfreqAll.txt"
removedMaleHetX(plink, inputPrefix, hhSubjCutOff, 
                outputPrefix, outputSubjHetFile, 
                outputRetainSubjectFile, outputHetSNPfile) 

inputPrefix <- "2_02_removedInstHetX"
outputPrefix <- "2_12_removedSnpHweFemaleX" 
genoQC(plink, inputPrefix, 
       snpMissCutOffpre=0.05, 
       sampleMissCutOff=0.02, 
       Fhet=0.2, 
       snpMissCutOffpost=0.02, 
       snpMissDifCutOff=0.02,
       femaleChrXmissCutoff=0.05, 
       pval4autoCtl=0.000001, 
       pval4femaleXctl=0.000001, outputPrefix, keepInterFile=TRUE)
 

## step 13 
inputPrefix <- "2_12_removedSnpHweFemaleX" 
outputPC4subjFile <- "2_13_eigenvalAfterQC.txt"
outputPCplotFile <- "2_13_eigenvalAfterQC.png"
nThread <- 20
plotPCA4plink(gcta, inputPrefix, nThread, outputPC4subjFile, outputPCplotFile)


############################################################

## Run the following code, only after detecting the potential 
## PCA outliers. Otherwise, set cutoff as "NULL"
############################################################


## remove outliers 
cutoff <-  c(-0.4, 0.2)
cutoffSign <- "greater" ## not used if cutoff == NULL, and with two values 
inputPC4subjFile <- "2_13_eigenvalAfterQC.txt"
outputPC4outlierFile <- "2_13_eigenval4outliers.txt"
outputPCplotFile <- "2_13_removedOutliers.png"
outputPrefix <- "2_13_removedOutliers" 
removeOutlierByPCs(plink, gcta, inputPrefix, nThread=20, 
                   cutoff, cutoffSign, inputPC4subjFile, 
                   outputPC4outlierFile, outputPCplotFile, outputPrefix) 

setwd("..")

## runtime
t4genoQC <- proc.time() - t4genoQCTmp
print(t4genoQC)
runTimeList$t4genoQC <- t4genoQC

t4checkAlignTmp <- proc.time()   

############################################################
### code chunk number 3: check the alignment
############################################################   
system("cp ./2-genoQC/2_13_removedOutliers.* ./3-checkAlign/ ")
setwd("./3-checkAlign/")
renamePlinkBFile(inputPrefix="2_13_removedOutliers", 
                 outputPrefix="3_1_QCdata", action="move")


inputFile <- paste0(impRefDIR,"*.legend.gz")  
## To recode the intermediate disk usage >> keep in the same directory.
bimReferenceFile <- "bimImputeRef.txt"

## If not available  >>> 
.prepareLegend2bim(inputFile, referencePanel, 
                   outputFile=bimReferenceFile, ncore=25) 
## If not available  <<< 

inputPrefix <- "3_1_QCdata" 
out2.snp <- "3_2_snpSameNameDiffPos"
out2 <- "3_2_removedSnpSameNameDiffPos"
out3 <- "3_3_removedSnpMissPos"
out3.snp <- "3_3_snpMissPos"
out4 <- "3_4_removedSnpDiffAlleles"
out4.snp <- "3_4_snpDiffAlleles"
out4.snpRetained <- "3_4_snpImpRefAlleles"
checkAlign2ref(plink, inputPrefix, referencePanel, bimReferenceFile, 
               out2, out2.snp, out3, out3.snp, 
               out4, out4.snp, out4.snpRetained, nCore=25)
setwd("..") 

## runtime
t4checkAlign <- proc.time() - t4checkAlignTmp
print(t4checkAlign)
runTimeList$t4checkAlign <- t4checkAlign

t4imputeTmp <- proc.time()   

############################################################
### code chunk number 4: Imputation
############################################################
## step 1 
## Remove monomorphic SNPs from lifted/QC-ed data  
inputPrefix4aligned2impRef <- "3_4_removedSnpDiffAlleles" 
outputPrefix <- "4_1_removedMonoSnp"
outputMonoSNPfile <- "4_1_snpMonoRemoved.txt" # will be used in step 4,5.

## copy plink files from last step; 
system(paste0("cp ./3-checkAlign/", 
       inputPrefix4aligned2impRef, ".* ./4-imputation/"))
## remove Monomorphic SNPs
setwd("4-imputation")
removedMonoSnp(plink, inputPrefix=inputPrefix4aligned2impRef, 
               outputPrefix, outputSNPfile=outputMonoSNPfile) 
# step 2 
#########################################################################
######################################################################### 
# imputation main pipeline
 
inputPrefix <- "4_1_removedMonoSnp"  
outputPrefix <- "gwasImputedFiltered"
prefix4final <- "gwasImputed"   
outputInfoFile <- "infoScore.txt"
tmpImputeDir <- paste0("tmp", referencePanel)
phaseImpute2(inputPrefix, outputPrefix, prefix4final,
            plink, shapeit, impute2, gtool, 
            windowSize=3000000, effectiveSize=20000, 
            nCore=40, threshold=0.9, infoScore=0.6, outputInfoFile, 
            referencePanel, impRefDIR, tmpImputeDir, keepTmpDir=TRUE)


# effectiveSize=20000
# nCore4phase=1
# nThread=40 
# nCore4impute=40
# threshold=0.9
# nCore4gtool=40
# infoScore=0.6

# ## alternatively
tmpImputeDir <- paste0("imp4tmp", referencePanel)
phaseImpute4(inputPrefix, outputPrefix, prefix4final,
            plink, shapeit, impute4, qctool, gtool, 
            windowSize=3000000, effectiveSize=20000, 
            nCore=40, threshold=0.9, infoScore=0.6, outputInfoFile, 
            referencePanel, impRefDIR, tmpImputeDir, keepTmpDir=TRUE)
  
##################################################### After imputation

## runtime
t4impute <- proc.time() - t4imputeTmp
print(t4impute)
runTimeList$t4impute <- t4impute

t4postImputeTmp <- proc.time()   

## step 2 
## Final imputed results, including bad imputed genotypes.
imputedDatasetfn <- "4_2_imputedDataset"
system(paste0("scp ./", tmpImputeDir, "/6-finalResults/gwasImputed.* ."))
renamePlinkBFile(inputPrefix="gwasImputed", 
                 outputPrefix="4_2_imputedDataset", action="move")
## step 3
## Filtered imputed data set; Remove imputed SNPs with (info < 0.6), 
## only retain "Good" SNPs. 
wellImputedfn <- "4_3_wellImputeData" 
snpImputedInfoScoreFile <- "4_3_snpImputedInfoScore.txt"
system(paste0("scp ./", tmpImputeDir, 
       "/6-finalResults/gwasImputedFiltered.* . "))
renamePlinkBFile(inputPrefix="gwasImputedFiltered", 
                 outputPrefix=wellImputedfn, action="move") 
system(paste0("scp ./", tmpImputeDir, "/6-finalResults/", 
       outputInfoFile, " ", snpImputedInfoScoreFile))

## step 4 
## Remove previous identified monomorphic SNPs in the imputed dataset. 
## Note that snps with same genomic position but can have different snp name.
removedMonoSnpAfter <- "4_4_removedMonoSnpAfter"
## if no monomorphic SNPs:
if (file.size(paste0(outputMonoSNPfile)) == 0 ){
    renamePlinkBFile(inputPrefix=wellImputedfn, 
                     outputPrefix=removedMonoSnpAfter, action="copy")   
} else { 
    ## extract PLINK files contain only monomorphic SNPs from 
    ## the original aligned (lifted and QC-ed) data set.
    system(paste0(plink, " --bfile ", inputPrefix4aligned2impRef, 
           " --extract ", outputMonoSNPfile, " --make-bed --out ", 
           inputPrefix4aligned2impRef, "Tmp")) 
    bim1 <- read.table(paste0(inputPrefix4aligned2impRef, "Tmp.bim"), 
                       stringsAsFactors=F)
    system(paste0("awk '{print $1, $2, $4}' ", 
           wellImputedfn, ".bim > tmpFilterImp.txt"))
    bim2 <- read.table("tmpFilterImp.txt", stringsAsFactors=F) 
    colnames(bim1) <- c("chr", "rsID", "gd", "pos", "a0", "a1") 
    colnames(bim2) <- c("chr", "rsID", "pos")
    outputFile <- "tmp.txt"
    .snpSharedPos(inputFile1=bim1, inputFile2=bim2, outputFile, nCore=25) 
    system(paste0(plink, " --bfile ", wellImputedfn, 
           " --exclude tmp.txt --make-bed --out ", removedMonoSnpAfter))
    system("rm tmpFilterImp.txt tmp.txt")
}
 
## step 5
## Add previous identified monomorphic SNPs in the imputed dataset.
addedMonoSnpAfter <- "4_5_addedMonoSnpAfter" 
 ## if no monomorphic SNPs:
if ( file.size(paste0(outputMonoSNPfile))==0 ){ 
        renamePlinkBFile(inputPrefix=wellImputedfn, 
                         outputPrefix=addedMonoSnpAfter, action="copy")  
} else { 
    ## merge both datasets
    system(paste0(plink, " --bfile ", removedMonoSnpAfter, " --bmerge ", 
           inputPrefix4aligned2impRef, ".bed ", 
           inputPrefix4aligned2impRef, ".bim ", 
           inputPrefix4aligned2impRef, ".fam ", 
           "--make-bed --out ", addedMonoSnpAfter))
}  

## step 6
## Remove SNPs which have a non missing value for less then 20 instances. 

inputPrefix <- addedMonoSnpAfter  
missCutoff <- 20
outputPrefix <- "4_6_removedSnpMissPostImp"
outRemovedSNPfile <- "4_6_snpRemovedMissPostImp.txt"
outRetainSNPfile <- "4_6_snpRetainMissPostImp.txt"
removedSnpMissPostImp(plink, inputPrefix, missCutoff, 
                      outRemovedSNPfile, outRetainSNPfile, outputPrefix)
 

## special case: 
## 1.) if the input ref panel is from phase3. 
## (2) Impute4 modified imputed output rsID if the original SNPs are not given.
## e.g. rs456706:11001104:C:T << rsID  rs456706

if (referencePanel == "1000Gphase3"){

    outRename <- "4_6_renameSNP"
    renamePlinkBFile(inputPrefix=outputPrefix, 
                     outputPrefix=outRename, action="move") 
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
        snpRenameFile <- "4_6_snpRename.txt"
        write.table(snpNewOld, file=snpRenameFile, quote=FALSE, 
                    row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
        system(paste0(plink, " --bfile ", outRename, " --update-name ", 
               snpRenameFile, " 1 2 --make-bed --out ", outRename, "tmp"))
        renamePlinkBFile(inputPrefix=paste0(outRename, "tmp"), 
                         outputPrefix=outputPrefix, action="move")  
    }     
}
   
setwd("..") 


############################################################
### code chunk number 5: Data subset and expansion 
############################################################
inputPrefix <- "4_6_removedSnpMissPostImp"
inputOriginalQCed <- "3_1_QCdata"
reducedToSpecificfn <- "5_1_reducedToSpecific"
extSpecificDiffAllelefn <- "5_2_extSpecificDiffAllele"
extSpecificMissPosfn <- "5_3_extSpecificMissPos"
extSpecificDiffPosfn <- "5_4_extSpecificDiffPos" 
dir5 <- "./5-reductAndExpand/"
## imputed dataset
system(paste0("scp ./4-imputation/", inputPrefix, ".* ", dir5)) 
## original but QC-ed dataset
system(paste0("scp ./3-checkAlign/", inputOriginalQCed, ".* ", dir5))

system(paste0("scp ./3-checkAlign/3_4_snpImpRefAlleles.txt ", dir5))
system(paste0("scp ./3-checkAlign/3_4_snpDiffAlleles.txt ", dir5))
system(paste0("scp ./3-checkAlign/3_3_snpMissPos.txt ", dir5))
system(paste0("scp ./3-checkAlign/3_2_snpSameNameDiffPos.txt ", dir5))

setwd("5-reductAndExpand/") 

## 1. Reduce the imputed dataset to the SNPs before imputation. 

if (referencePanel == "1000Gphase3"){ 
    ## Before doing this, make sure there are no variant with 3+ alleles present.
    ## Remove if any.
    bimOriQCed <- read.table(paste0(inputOriginalQCed, ".bim"), 
                             stringsAsFactors=FALSE)
    bimPostImp <- read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)
    ## check duplicated.
    dupPostImp <- bimPostImp[duplicated(bimPostImp[,2]),]
    str(dupPostImp)
    allele3 <- intersect(dupPostImp[,2], bimOriQCed[,2])
    str(allele3)
    ## remove these duplicated and with 3+alleles 
    write.table(allele3, file="snpAllele3.txt", quote=FALSE, 
                        row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

    outAllele2 <- paste0(inputPrefix, "Allele2")
    system(paste0(plink, " --bfile ", inputPrefix, 
           " --exclude snpAllele3.txt --make-bed --out ",  outAllele2)) 
     
} else { outAllele2 <- inputPrefix }

system(paste0(plink, " --bfile ", outAllele2, 
       " --extract 3_4_snpImpRefAlleles.txt --make-bed --out ", 
       reducedToSpecificfn)) 
 
## 2. Add the SNPs with different alleles with their values 
## from the dataset before removing SNPs. 
if ( file.size(paste0("3_4_snpDiffAlleles.txt")) == 0 ){ 
     renamePlinkBFile(inputPrefix=reducedToSpecificfn, 
                      outputPrefix=extSpecificDiffAllelefn, action="copy")   
} else { 
    system(paste0(plink, " --bfile ", inputOriginalQCed, 
           " --extract 3_4_snpDiffAlleles.txt --make-bed --out tmp")) 
    system(paste0(plink, " --bfile ", reducedToSpecificfn, 
           " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", 
           extSpecificDiffAllelefn)) 
    system("rm tmp.*")
} 



## 3. Add the SNPs with missing positions with their values 
## from the dataset before removing SNPs. 
if ( file.size(paste0("3_3_snpMissPos.txt")) == 0 ){  
    renamePlinkBFile(inputPrefix=extSpecificDiffAllelefn, 
                     outputPrefix=extSpecificMissPosfn, action="copy")  
} else {
    system(paste0(plink, " --bfile ", inputOriginalQCed, 
           " --extract 3_3_snpMissPos.txt --make-bed --out tmp")) 
    system(paste0(plink, " --bfile ", extSpecificDiffAllelefn, 
           " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", 
           extSpecificMissPosfn)) 
    system("rm tmp.*")
}
      
## 4. Add the SNPs with different positions by their values 
## from the dataset before removing SNPs. 
if ( file.size(paste0("3_2_snpSameNameDiffPos.txt")) == 0 ){  
    renamePlinkBFile(inputPrefix=extSpecificMissPosfn, 
                     outputPrefix=extSpecificDiffPosfn, action="copy") 
} else {
    system(paste0(plink, " --bfile ", inputOriginalQCed, 
           " --extract 3_2_snpSameNameDiffPos.txt --make-bed --out tmp")) 
    system(paste0(plink, " --bfile ", extSpecificMissPosfn, 
           " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", 
           extSpecificDiffPosfn)) 
    system("rm tmp.*")
}
 

# system(paste0("rm ", inputPrefix, "* ", inputOriginalQCed, "*"))  
# system(paste0("rm  *.txt *.log")) 
setwd("..")

############################################################
### code chunk number 6: Final result
############################################################

dir6 <- "./6-finalResults/"
system(paste0("scp ./1-genoUpdate/1_01_metaData.txt ", dir6))
system(paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.* ", dir6))  
system(paste0("scp ./5-reductAndExpand/5_4_extSpecificDiffPos.* ", dir6))
setwd(dir6)
renamePlinkBFile(inputPrefix="4_6_removedSnpMissPostImp", 
                 outputPrefix="imputedSnpsDataset", action="move")
renamePlinkBFile(inputPrefix="5_4_extSpecificDiffPos", 
                 outputPrefix="specificSnpsDataset", action="move")



## runtime
t4postImpute <- proc.time() - t4postImputeTmp
print(t4postImpute)
runTimeList$t4postImpute <- t4postImpute 

totallist <- lapply(runTimeList, function(time){time[[3]]})
print(totallist)  

runTimeList$totalMin <- sum(unlist(totallist))/60
runTimeList$totalSec <- sum(unlist(totallist))

print(runTimeList)

############################################################
### code chunk number 6: Extending pipeline
############################################################
 
############################################################

## Run the following code, only after you have installed genipe
## and the specified imputation reference files.

############################################################

# ## additional configuration files for Genipe
# mainRefGenipe <- "/data/noether/dataProcessResults/10_Common/"
# impRefDir <- paste0(mainRefGenipe, "imputeReference/1000G_Phase3_2014/")
# fastaFile <- paste0(mainRefGenipe, "imputeReference/hg19/hg19.fasta")

# ## Impute genotypes using Genipe
# chrs <- 21
# inputPrefix <- "dataChr21"
# thread4impute2 <- 20 ## tune by yourself
# thread4shapeit <- 30
# segmentSize <- 3000000
# # imputedByGenipe(chrs, impRefDir, inputPrefix, shapeit, impute2, 
# #                plink, fastaFile, segmentSize, thread4impute2, thread4shapeit) 

# ## merge chunked genomic imputed results
# ## example
# chr <- 21 
# inputImpute2 <- "chr2.33000001_36000000.impute2"
# probability <- 0.9
# completionRate <- 0.98
# # info <- 0.6
# outputPrefix <- paste0("imputedChr", chr)
# # mergeByGenipe(inputImpute2, chr, probability, 
# #              completionRate, info, outputPrefix)
 
# ## extract imputed markers using Genipe 
# chr <- 3
# inputImpute2 <- paste0("chr", chr,".imputed.impute2")
# inputMAP <- paste0("chr", chr,".imputed.map")
# format <- "bed"
# prob <- 0.9
# outputPrefix <- paste0("imputedChr", chr)  
# # extractByGenipe(inputImpute2, inputMAP, outputPrefix, format, prob)
   
#  