
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

#  # Alternative Reference panel 2
# referencePanel <- "1000Gphase3" ## indicator 
# impRefDIR1kGp3 <- "/data/noether/dataRawReadOnly/reference/1000GP_Phase3/"
# impRefDIR <- paste0(impRefDIR1kGp3, "1000GP_Phase3/")

## Define required tools
toolDIR <- "/home/junfang.chen/Gimpute/tools/"
plink <- paste0(toolDIR, "plink")
gcta <- paste0(toolDIR, "gcta64") 
shapeit <- paste0(toolDIR, "shapeit") 
gtool <- paste0(toolDIR, "gtool")
## Gimpute has the following dependencies:  
qctool <- paste0(toolDIR, "qctool")

imputeTool <- "impute2"
if (imputeTool == "impute2"){
    impute <- paste0(toolDIR, "impute2")
} else if (imputeTool == "impute4"){
    impute <- paste0(toolDIR, "impute4.1_r291.2")
} else {
    print("Wrong imputeTool or no imputation tool is provided!")
}


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
nCores <- detectCores() 
plotPCA4plink(gcta, inputPrefix, nThread=nCores, 
              outputPC4subjFile, outputPCplotFile)


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
removeOutlierByPCs(plink, gcta, inputPrefix, nThread=nCores, 
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
                   outputFile=bimReferenceFile, ncore=nCores) 
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
               out4, out4.snp, out4.snpRetained, nCore=nCores)
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
outputPrefix <- "4_2_imputedDataset"   
outputInfoFile <- "4_2_snpImputedInfoScore.txt"
tmpImputeDir <- paste0("tmp", referencePanel)
phaseImpute(inputPrefix, outputPrefix,
             plink, shapeit, imputeTool, impute, qctool, gtool, 
             windowSize=3000000, effectiveSize=20000, 
             nCore=nCores, threshold=0.9, outputInfoFile, 
             referencePanel, impRefDIR, tmpImputeDir, keepTmpDir=TRUE)

 
# # ## alternatively
# tmpImputeDir <- paste0("imp4tmp", referencePanel)
# phaseImpute4(inputPrefix, outputPrefix,
#              plink, shapeit, impute4, qctool, gtool, 
#              windowSize=3000000, effectiveSize=20000, 
#              nCore=nCores, threshold=0.9, outputInfoFile, 
#              referencePanel, impRefDIR, tmpImputeDir, keepTmpDir=TRUE)
  
##################################################### After imputation
## runtime
t4impute <- proc.time() - t4imputeTmp
print(t4impute)
runTimeList$t4impute <- t4impute
t4postImputeTmp <- proc.time()   


inputPrefix4aligned2impRef <- "3_4_removedSnpDiffAlleles" 
inputPrefix <- "4_2_imputedDataset" 
out1 <- "4_3_wellImputeData"
out2 <- "4_4_removedMonoSnpAfter"
out3 <- "4_5_addedMonoSnpAfter"
out4 <- "4_6_removedSnpMissPostImp" 
outRemovedSNPfile <- "4_6_snpRemovedMissPostImp.txt"
outRetainSNPfile <- "4_6_snpRetainMissPostImp.txt"

postImpQC(inputPrefix, out1, out2, out3, out4,
          outputInfoFile, infoScore=0.6, inputPrefix4aligned2impRef, 
          missCutoff=20, outRemovedSNPfile, 
          outRetainSNPfile, referencePanel)
   
setwd("..") 


############################################################
### code chunk number 5: Data subset and expansion 
############################################################

inputPrefix <- "4_6_removedSnpMissPostImp"
inputQCprefix <- "3_1_QCdata"
snpRefAlleleFile <- "3_4_snpImpRefAlleles.txt"
snpDiffAlleleFile <- "3_4_snpDiffAlleles.txt"
snpMissPosFile <- "3_3_snpMissPos.txt"
snpSameNameDifPosFile <- "3_2_snpSameNameDiffPos.txt" 

dir5 <- "./5-reductAndExpand/"
## imputed dataset
system(paste0("scp ./4-imputation/", inputPrefix, ".* ", dir5)) 
system(paste0("scp ./3-checkAlign/", inputQCprefix, ".* ", dir5))
system(paste0("scp ./3-checkAlign/", snpRefAlleleFile, " ", dir5))
system(paste0("scp ./3-checkAlign/", snpDiffAlleleFile, " ", dir5))
system(paste0("scp ./3-checkAlign/", snpMissPosFile, " ", dir5))
system(paste0("scp ./3-checkAlign/", snpSameNameDifPosFile, " ", dir5))

setwd(dir5)
reducedToSpecificfn <- "5_1_reducedToSpecific"
specificDiffAllelefn <- "5_2_specificDiffAllele"
specificMissPosfn <- "5_3_specificMissPos"
specificDiffPosfn <- "5_4_specificDiffPos" 

reductExpand(referencePanel, inputPrefix, inputQCprefix, 
             snpRefAlleleFile, snpDiffAlleleFile, 
             snpMissPosFile, snpSameNameDifPosFile, 
             reducedToSpecificfn, specificDiffAllelefn, 
             specificMissPosfn, specificDiffPosfn)

setwd("..")

############################################################
### code chunk number 6: Final result
############################################################

dir6 <- "./6-finalResults/"
system(paste0("scp ./1-genoUpdate/1_01_metaData.txt ", dir6))
system(paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.* ", dir6))  
system(paste0("scp ./5-reductAndExpand/5_4_specificDiffPos.* ", dir6))
setwd(dir6)
renamePlinkBFile(inputPrefix="4_6_removedSnpMissPostImp", 
                 outputPrefix="imputedSnpsDataset", action="move")
renamePlinkBFile(inputPrefix="5_4_specificDiffPos", 
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