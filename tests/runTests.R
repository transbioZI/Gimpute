
# File   : runTests.R
# Author : Junfang Chen
# Version0: 28 Jun 2016
# VersionX: 14 Jun 2018
  



## Define the directory where you place the required/downloaded tools.
## required tools
plink <- "/home/junfang.chen/Gimpute/tools/plink"
gcta <- "/home/junfang.chen/Gimpute/tools/gcta64" 
shapeit <- "/home/junfang.chen/Gimpute/tools/shapeit"
impute2 <- "/home/junfang.chen/Gimpute/tools/impute2"
gtool <- "/home/junfang.chen/Gimpute/tools/gtool"

## configuration files/. 
## configuration directory: where you put all your files in this folder for the configuration. 
chipAnnoFile= NULL 
dupSampleIDFile <- NULL ## if it's not NULL, then the file stores the duplicated sample IDs should be placed in this directory.
excludedProbeIdsFile <- NULL ## if it's not NULL, then the file stores the probes which have to be excluded should be placed in this directory.

## The directory where you specify your imputation reference files:
impRefDIR <- "/data/noether/datatmp-nobackup/tb2refDatabase/imputeRef/1000Gphase1/"
 
## global parameters/variables
ancestrySymbol <- NULL
caseControl <- FALSE 

  

############################################################
## code chunk number 1: SNP information update  
############################################################

library(Gimpute)

## go to the directory where you"d like to generate the data
# make the following directories and sub-directories  
system("mkdir 1-conversion")
system("mkdir 2-QC")
system("mkdir 3-lifting")
system("mkdir 4-imputation")
system("mkdir 5-reductAndExpand")
system("mkdir 6-finalResults")
 
## step 1
## copy plink files and meta information file
bedFile <- system.file("extdata", "study.bed", package="Gimpute")
bimFile <- system.file("extdata", "study.bim", package="Gimpute") 
famFile <- system.file("extdata", "study.fam", package="Gimpute")
metadataFile <- system.file("extdata", "1_01_metaData.txt", package="Gimpute")

system( paste0("scp ", bedFile, " ./1-conversion/") ) 
system( paste0("scp ", bimFile, " ./1-conversion/") ) 
system( paste0("scp ", famFile, " ./1-conversion/") ) 
system( paste0("scp ", metadataFile, " ./1-conversion/") )  
setwd("./1-conversion/") 
  
 
## step 2   
## dupSampleIDFile; defined as a global variable
inputPrefix <- "study"
outputPrefix <- "1_02_removedExclInst" 
removeDupID(plink, dupSampleIDFile, inputPrefix, outputPrefix)
 

# step 3 replace group IDs 
metaDataFile <- "1_01_metaData.txt"
inputPrefix <- "1_02_removedExclInst"
outputPrefix <- "1_03_replacedGroupAndSex"
replaceGroupIdAndSex(plink, inputPrefix, metaDataFile, outputPrefix)

# step 4 remove instances without group IDs
metaDataFile <- "1_01_metaData.txt"
inputPrefix <- "1_03_replacedGroupAndSex"
outputPrefix <- "1_04_removedNoGroupId"
removeNoGroupId(plink, inputPrefix, outputPrefix)

## step5 remove instances with improper ancestry 
metaDataFile <- "1_01_metaData.txt" 
inputPrefix <- "1_04_removedNoGroupId"
outputPrefix <- "1_05_removedWrongAnceInst"
removedWrongAnceInst(plink, inputPrefix, metaDataFile, 
					 ancestrySymbol, outputPrefix)

## step 6 
inputPrefix <- "1_05_removedWrongAnceInst"
excludedProbeIdsFile <- excludedProbeIdsFile    
outputPrefix <- "1_06_removedExclProbe"  
removedExclProbe(plink, inputPrefix, excludedProbeIdsFile, outputPrefix) 

 
## step 7 (Optional, if chip annotation file is not given)
inputPrefix <- "1_06_removedExclProbe"
chipAnnoFile <- chipAnnoFile  
chipType <- "illumina"
outputPrefix <- "1_07_removedUnmapProbes"   
outputSNPunmapFile <- "1_07_probesUnmapped2ChipRef.txt"
removedUnmapProbes(plink, inputPrefix, chipAnnoFile, 
				   outputPrefix, outputSNPunmapFile)
 

## step 8 (Optional, if chip annotation file is not given)
inputPrefix <- "1_07_removedUnmapProbes" 
chipAnnoFile <- chipAnnoFile  
chipType <- "illumina"
outputSNPdupFile <- "1_08_probesDouble.txt"
outputPrefix <- "1_08_removedDoubleProbes"   
removedDoubleProbes(plink, inputPrefix, chipAnnoFile, 
				    chipType, outputSNPdupFile, outputPrefix)
 
## step 9 (Optional, if chip annotation file is not given)
inputPrefix <- "1_08_removedDoubleProbes" 
chipAnnoFile <- chipAnnoFile  
chipType <- "illumina"

outputPrefix <- "1_09_updatedSnpInfo"   
updatedSnpInfo(plink, inputPrefix,  chipAnnoFile, chipType, outputPrefix)


## step 10 	
inputPrefix <- "1_09_updatedSnpInfo"
outputPrefix <- "1_10_changedXyChr"
changedXyChr(plink, inputPrefix, outputPrefix)


## step 11 
inputPrefix <- "1_10_changedXyChr"
outputPrefix <- "1_11_removedYMtSnp"
removedYMtSnp(plink, inputPrefix, outputPrefix)


## remove unwanted files
system( paste0("rm  *.log *.hh") ) 
## change dir to the main directory
setwd("..")   





############################################################
### code chunk number 2: Quality Control
############################################################


## step 0
## copy the last output plink files from 1-conversion
inputPrefix4QC <- "1_11_removedYMtSnp"
system( paste0("scp ./1-conversion/", inputPrefix4QC, ".*", " ./2-QC/") )
setwd("./2-QC/")


## step 1 
inputPrefix <- inputPrefix4QC
hhCutOff <- 0.005 ##  can be tuned
outputPrefix <- "2_01_removedSnpHetX" 
outputFile_SNPhhFreqAll <- "2_01_snpHHfreqAll.txt"
outputFile_SNPhhFreqRetained <- "2_01_snpHHfreqRetained.txt"
removedSnpHetX(plink, inputPrefix, hhCutOff, outputPrefix, 
		 	   outputFile_SNPhhFreqAll, outputFile_SNPhhFreqRetained)


## step 2  2_02_removedHetXInst

inputPrefix <- "2_01_removedSnpHetX"
hhSubjCutOff <- 15
outputPrefix <- "2_02_removedInstHetX"
outputFile_subjHetFreqAll <- "2_02_instHetXfreqAll.txt" 
outputFile_subjHetFreqRetained <- "2_02_instHetXfreqRetained.txt"  
outputFile_SNPhhFreqAll <- "2_02_snpHHfreqAll.txt"
removedMaleHetX(plink, inputPrefix, hhSubjCutOff, 
				outputPrefix, outputFile_subjHetFreqAll, 
				outputFile_subjHetFreqRetained, outputFile_SNPhhFreqAll) 

## step 3 
# 3. Set all heterozygous alleles of SNPs of the chromosome 23 for males
inputPrefix <- "2_02_removedInstHetX" 
outputPrefix <- "2_03_setHeteroHaploMissing" 
setHeteroHaploMissing(plink, inputPrefix, outputPrefix)
 

## step 4   SNP missingness < 0.05 (before sample removal);  
inputPrefix <- "2_03_setHeteroHaploMissing" 
snpMissCutOff <- 0.05 #
outputPrefix <- "2_04_removedSnpMissPre" 
removedSnpMiss(plink, snpMissCutOff, inputPrefix, outputPrefix)

## step 5 
# subject missingness < 0.02; 
inputPrefix <- "2_04_removedSnpMissPre" 
sampleMissCutOff <- 0.02
outputPrefix <- "2_05_removedInstMiss" 
removedInstMiss(plink, sampleMissCutOff, inputPrefix, outputPrefix)
 
## step 6 
inputPrefix <- "2_05_removedInstMiss" 
Fhet <- 0.2 #
outputPrefix <- "2_06_removedInstFhet" 
removedInstFhet(plink, Fhet, inputPrefix, outputPrefix)
  	

## step 7
inputPrefix <- "2_06_removedInstFhet"  
outputPrefix <- "2_07_removedParentIdsMiss" 
removedParentIdsMiss(plink, inputPrefix, outputPrefix)


## step 8 
snpMissCutOff <- 0.02
inputPrefix <- "2_07_removedParentIdsMiss"  
outputPrefix <- "2_08_removedSnpMissPost" 
removedSnpMiss(plink, snpMissCutOff, inputPrefix, outputPrefix)


## step 9  ## caseControl
 ##  Remove SNPs with difference >= 0.02 of SNP missingness between cases and controls.
inputPrefix <- "2_08_removedSnpMissPost"  
outputPrefix <- "2_09_removedSnpMissDiff" 
snpMissDifCutOff <- 0.02  
removedSnpMissDiff(plink, inputPrefix, 
				   snpMissDifCutOff, outputPrefix, caseControl) 

 
## step 10
femaleChrXmissCutoff <- 0.05
inputPrefix <- "2_09_removedSnpMissDiff"  
outputPrefix <- "2_10_removedSnpFemaleChrXmiss" 
removedSnpFemaleChrXmiss(plink, femaleChrXmissCutoff, inputPrefix, outputPrefix)
  

## step 11

pval <- 0.000001
inputPrefix <- "2_10_removedSnpFemaleChrXmiss"  
outputFile_pVal <- "2_11_snpHwePvalAutoCt.txt" 
outputSNPfile <-  "2_11_snpRemovedHweAutoCt.txt" 
outputPrefix <- "2_11_removedSnpHweAutoCt" 

removedSnpHWEautoControl(plink, inputPrefix, pval, 
						 outputFile_pVal, outputSNPfile, outputPrefix)


## step 12 
pval <- 0.000001
inputPrefix <- "2_11_removedSnpHweAutoCt"   
outputFile_pVal <- "2_12_snpHwePvalfemaleXct.txt" 
outputSNPfile <- "2_12_snpRemovedHweFemaleXct.txt" 
outputPrefix <- "2_12_removedSnpHweFemaleXct" 
removedSnpFemaleChrXhweControl(plink, inputPrefix, pval, 
							   outputFile_pVal, outputSNPfile, outputPrefix)
 
 
 
 
################################################
## get the ethnic group info
setwd("..")
metaDataFile <- "1_01_metaData.txt"
system( paste0("scp ./1-conversion/", metaDataFile, " ./2-QC/") )
setwd("./2-QC/")
################################################ 

## step 13 

inputPrefix <- "2_12_removedSnpHweFemaleXct" ## the output from step 12
outputPC4subjFile <- "2_13_eigenvalAfterQC.txt"
outputPCplotFile <- "2_13_eigenvalAfterQC.png"
plotPCA4plink(gcta, inputPrefix, outputPC4subjFile, outputPCplotFile)

 
## remove outliers
inputPrefix <- "2_12_removedSnpHweFemaleXct"
cutoff <-  NULL 
cutoffSign <- "greater" ## not used if cutoff == NULL


inputPC4subjFile <- "2_13_eigenvalAfterQC.txt"
outputPC4outlierFile <- "2_13_eigenval4outliers.txt"
outputPCplotFile <- "2_13_removedOutliers.png"
outputPrefix <- "2_13_removedOutliers" 

removeOutlierByPCs(plink, gcta, inputPrefix, 
				   cutoff, cutoffSign, inputPC4subjFile, 
				   outputPC4outlierFile, outputPCplotFile, outputPrefix)
  
 
######################## 
## remove unwanted files
## remove meta file and metaDataFile+AA.txt
system( paste0("rm  ", metaDataFile) ) 
system( paste0("rm  *.log *.hh") )
## change dir 
setwd("..") 

 


############################################################
### code chunk number 3: align/lift study genotypes to the imputation reference
############################################################

 
## step 1 copy the QC-PLINK files to next section   
system("cp ./2-QC/2_13_removedOutliers.* ./3-lifting/ ")
setwd("./3-lifting/")
renamePlinkBFile(inputPrefix="2_13_removedOutliers", 
				 outputPrefix="3_1_liftedDataset", action="move")
 
inputFile <- paste0(impRefDIR,"*.legend.gz")  
bimReferenceFile <- paste0(impRefDIR, "bimImputeRef.txt")
## take less than 1 minute
prepareLegend2bim(inputFile, outputFile=bimReferenceFile, ncore=25) 
  
## 2. Remove SNPs for which the name has a different position (i.e. combination of
# base pair position and chromosome) in the imputation reference files. 
inputPrefix <- "3_1_liftedDataset" 
out2.snp <- "3_2_snpSameNameDiffPos"
out2 <- "3_2_removedSnpSameNameDiffPos"
  
# 3. Remove SNPs which bp and chr position are not contained 
## in the imputation reference files.
out3 <- "3_3_removedSnpMissPos"
out3.snp <- "3_3_snpMissPos"
  
## 4. Remove SNPs which have an allele which is not in the imputation reference 
## files for that SNP.
out4 <- "3_4_removedSnpDiffAlleles"
out4.snp <- "3_4_snpDiffAlleles"
out4.snpRetained <- "3_4_snpImpRefAlleles"
 
checkAlign2ref(plink, inputPrefix, bimReferenceFile, out2, out2.snp, 
			   out3, out3.snp, out4, out4.snp, out4.snpRetained, nCore=25)


system( paste0("rm  *.log *.hh") ) 
setwd("..")





############################################################
### code chunk number 4: Imputation
############################################################

## step 1 
## Remove monomorphic SNPs from lifted/QC-ed data 
## will also be used in step 4 and 5;
inputPrefix4aligned2impRef <- "3_4_removedSnpDiffAlleles" 
outputPrefix <- "4_1_removedMonoSnp"
outputMonoSNPfile <- "4_1_snpMonoRemoved.txt" # will be used in step 4 and 5.

## copy plink files from last step; 
system(paste0("cp ./3-lifting/", inputPrefix4aligned2impRef, ".* ./4-imputation/"))
## remove Monomorphic SNPs
setwd("4-imputation")
removedMonoSnp(plink, inputPrefix=inputPrefix4aligned2impRef, 
			   outputPrefix, outputSNPfile=outputMonoSNPfile)
 
# ## remove unwanted plink files << 
## will also be used in step 4 and 5; after that you can remove.
# system(paste0("rm ", inputPrefix, "*"))


## step 2 
##########################################################################
########################################################################## 
## imputation main pipeline

## sub-steps
# step 2.1 chrWiseSplit.R
# step 2.2 chunk4eachChr.R
# step 2.3 prePhasingByShapeit.R
# step 2.4 imputedByImpute2.R 
# step 2.5 convertImpute2ByGtool 
# step 2.6 mergePlinkData.R 
# step 2.7 filterImputeData.R 
 
## One must create directories for storing temporal imputation output files 
## The name of these directories must be fixed for the sake of the subsequent steps.
tmp4imputeDIR <- "tmpImpute"
system( paste0("mkdir ", tmp4imputeDIR))
setwd( tmp4imputeDIR ) ## 
## sub-directories  
system( "mkdir 1-dataFiles")
system( "mkdir 2-chunkFile") 
system( "mkdir 3-phaseResults")
system( "mkdir 4-imputeResults")
system( "mkdir 5-postImpute")
system( "mkdir 6-finalResults") 
setwd("..")  


########################################
######################################## chrWiseSplit.R
## step 2.1

# copy plink files without monomorphic SNPs; prepare for the imputation.
outputPrefix <- "4_1_removedMonoSnp"
outputPrefixTmp  <- "gwas_data_chr"  
system( paste0("scp 4_1_removedMonoSnp.*  ./", tmp4imputeDIR, "/1-dataFiles/"))
 
  
setwd(paste0("./", tmp4imputeDIR, "/1-dataFiles/"))  

renamePlinkBFile(inputPrefix="4_1_removedMonoSnp", 
				 outputPrefix="gwas_data_chr", action="move")

bimCurrent <- read.table(file=paste0(outputPrefixTmp, ".bim"), 
						 stringsAsFactors=FALSE)  
currentChr <- names(table(bimCurrent[,1]))
print(currentChr) 
 
inputPrefix  <- outputPrefixTmp 
outputPrefix  <- outputPrefixTmp ## the chromosome number will be appended
chrX_PAR1suffix <- "X_PAR1"
chrX_PAR2suffix <- "X_PAR2"
PAR <- chrWiseSplit(plink, inputPrefix, chrX_PAR1suffix, chrX_PAR2suffix, nCore=25)
 
PAR 
 
if (PAR[[1]]) {par1 <- "X_PAR1"} else {par1 <- NULL}
if (PAR[[2]]) {par2 <- "X_PAR2"} else {par2 <- NULL}

  

########################################
######################################## chunk4eachChr.R
## step 2.2

inputPrefix <- "gwas_data_chr"
outputPrefix <- "chunks_chr" 
chrs <- c(currentChr, par1, par2)   
windowSize <- 3000000 
chunk4eachChr(inputPrefix, outputPrefix, chrs, windowSize) 

setwd("..") 
system( paste0("mv ./1-dataFiles/", outputPrefix, "*.txt  ./2-chunkFile/"))



## step 2.3  

# define directories
dataDIR <- "./1-dataFiles/" 
phaseDIR <- "./3-phaseResults/"
 
prefix4plinkEachChr <- "gwas_data_chr"  
# chrs <- c(1, 11, 23, par1, par2)   
nThread <- 40
effectiveSize <- 20000 
nCore <- 1 
prePhasingByShapeit(shapeit, chrs, dataDIR, prefix4plinkEachChr, 
					impRefDIR, phaseDIR, nThread, effectiveSize, nCore)

 

# step 2.4   
# define directories
prefixChunk <- "./2-chunkFile/chunks_chr"  ## 
phaseDIR <- "./3-phaseResults/" 
imputedDIR <- "./4-imputeResults/"

prefix4plinkEachChr <- "gwas_data_chr"
# chrs <- c(1,23, par1, par2)     
nCore <- 30 ## try to tune it for your own data size
effectiveSize <- 20000  

imputedByImpute2(impute2, chrs, prefixChunk, phaseDIR, impRefDIR, 
				 imputedDIR, prefix4plinkEachChr, nCore, effectiveSize)


  
 
# step 2.5   
 ## extract only SNPs (without INDELs)
#######################################################
setwd("./4-imputeResults")  
## extract only SNPs starting with "rs";  .
ls <- system("ls gwas*.impute2", intern=T)
snpPrefix <- "rs" 

biglists <- as.list(ls)
mclapply(biglists, function(i){ 
	arg1 <- paste0(" awk '{if(length($4) == 1 && length($5) == 1) print}'")
	arg2 <- paste0(i, "noINDEL.impute2")   
	system(paste0("grep '", snpPrefix, "' ", i, " | ", arg1, " > ", arg2))
}, mc.cores=38) 

setwd("..")  
####################################################### <<<
 
 
# chrs <- c(1,23, par1, par2)
   
prefixChunk <- "./2-chunkFile/chunks_chr"  ## 
phaseDIR <- "./3-phaseResults/" 
imputedDIR <- "./4-imputeResults/"
prefix4plinkEachChr <- "gwas_data_chr"
suffix4imputed <- ".impute2noINDEL.impute2"
# suffix4imputed <- ".impute2"
postImputeDIR <- "./5-postImpute/" 

nCore <- 30   
convertImpute2ByGtool(gtool, chrs, prefixChunk, phaseDIR, imputedDIR, 
				      prefix4plinkEachChr, suffix4imputed, postImputeDIR, nCore)

 


# step 2.6  
####################################################### Modify missing genotype format.
setwd("./5-postImpute/")  
# replace 'N' in the .ped files into 0 > missing values.
chrslist <- as.list(chrs)
prefix4plinkEachChr <- "gwas_data_chr" ## just for parallel computing
fn <- mclapply(chrslist, function(i){
	system(paste0("sed -i 's/N/0/g' ", prefix4plinkEachChr, i, ".*ped "))
}, mc.cores=nCore)
## Note: check if also any "N" in *.fam files. If so, change back after merging. 
####################################################### <<< 

prefix4plinkEachChr <- "gwas_data_chr"
prefix4mergedPlink <- "gwasImputed"
par1 <- "X_PAR1"
par2 <- "X_PAR2"
# chrs <- c(1,23, par1, par2)   
mergePlinkData(plink, chrs, prefix4plinkEachChr, prefix4mergedPlink, nCore)
 
 

#######################################################
setwd("..")
system("mv ./5-postImpute/gwasImputed* ./6-finalResults/") 
system("mv ./4-imputeResults/*.impute2_info ./6-finalResults/") 

# step 2.7    
setwd("./6-finalResults")
suffix4impute2info <- ".impute2_info"
outputInfoFile <- "impute2infoUpdated.txt"
infoScore <- 0.6
badImputeSNPfile <- "badImputeSNPs.txt"
inputPrefix <- prefix4mergedPlink
outputPrefix <- "gwasImputedFiltered" 
filterImputeData(plink, suffix4impute2info, 
				 outputInfoFile, infoScore, badImputeSNPfile, 
				 inputPrefix, outputPrefix)

setwd("..")
setwd("..")

###################################################### ###### After imputation

## step 2 
## Final imputed results >> 
  
## output file name change to: 
imputedDatasetfn <- "4_2_imputedDataset"
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/gwasImputed.* .") )
renamePlinkBFile(inputPrefix="gwasImputed", 
				 outputPrefix="4_2_imputedDataset", action="move")

## step 3
## Filtered imputed data set; Remove imputed SNPs with (info < 0.6), 
## only retain "Good" SNPs.
prefix4filteredPlink <- "gwasImputedFiltered"
filteredImputedDatasetfn <- "4_3_removedSnpInfoPostImp" 
snpWithBadInfoFile <- "4_3_snpRemovedInfoPostImp.txt"
snpImputedInfoScoreFile <- "4_3_snpImputedInfoScore.txt"
 
system(paste0("scp ./", tmp4imputeDIR, "/6-finalResults/gwasImputedFiltered.* . "))
renamePlinkBFile(inputPrefix="gwasImputedFiltered", 
				 outputPrefix="4_3_removedSnpInfoPostImp", action="move")  

system(paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", 
	   outputInfoFile, " ", snpImputedInfoScoreFile))
system(paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", 
	   badImputeSNPfile, " ", snpWithBadInfoFile))
 
## step 4 
## Remove previous identified monomorphic SNPs in the imputed dataset. 
## Note that snps with same genomic position but can have different snp name.
filteredImputedDatasetfn <- "4_3_removedSnpInfoPostImp" 
removedMonoSnpAfter <- "4_4_removedMonoSnpAfter"

## if no monomorphic SNPs:
if ( file.size(paste0(outputMonoSNPfile)) == 0 ){
		renamePlinkBFile(inputPrefix=filteredImputedDatasetfn, 
				 	 	 outputPrefix=removedMonoSnpAfter, action="copy")   

} else { 
	## extract PLINK files contain only monomorphic SNPs from 
	## the original aligned (lifted and QC-ed) data set.
	system(paste0(plink, " --bfile ", inputPrefix4aligned2impRef, 
		   " --extract ", outputMonoSNPfile, " --make-bed --out ", 
		   inputPrefix4aligned2impRef, "Tmp") ) 
	bim1 <- read.table(paste0(inputPrefix4aligned2impRef, "Tmp.bim"), 
					   stringsAsFactors=F)
	system(paste0("awk '{print $1, $2, $4}' ", 
		   filteredImputedDatasetfn, ".bim > tmpFilterImp.txt"))
	bim2 <- read.table("tmpFilterImp.txt", stringsAsFactors=F) 
	colnames(bim1) <- c("chr", "rsID", "gd", "pos", "a0", "a1") 
	colnames(bim2) <- c("chr", "rsID", "pos")
	outputFile <- "tmp.txt"
	snpSharedPos(inputFile1=bim1, inputFile2=bim2, outputFile, nCore=25) 
	system(paste0(plink, " --bfile ", filteredImputedDatasetfn, 
		   " --exclude tmp.txt --make-bed --out ", removedMonoSnpAfter))
	system("rm tmpFilterImp.txt tmp.txt")
} 

 
## step 5
## Add previous identified monomorphic SNPs in the imputed dataset.
addedMonoSnpAfter <- "4_5_addedMonoSnpAfter" 
 ## if no monomorphic SNPs:
if ( file.size(paste0(outputMonoSNPfile))==0 ){ 
		renamePlinkBFile(inputPrefix=filteredImputedDatasetfn, 
				         outputPrefix=addedMonoSnpAfter, action="copy")  
} else { 
  
	## merge both datasets
	system(paste0(plink, " --bfile ", removedMonoSnpAfter, " --bmerge ", 
		   inputPrefix4aligned2impRef, ".bed ", 
		   inputPrefix4aligned2impRef, ".bim ", 
		   inputPrefix4aligned2impRef, ".fam ", 
		   "--make-bed --out ", addedMonoSnpAfter))
	## remove tmp files
	# system( paste0("rm tmp.txt") )
	# system( paste0("rm ", inputPrefix4aligned2impRef, "*") )
}  


## step 6
## Remove SNPs which have a non missing value for less then 20 instances. 

inputPrefix <- addedMonoSnpAfter  
missCutoff <- 20
outputPrefix <- "4_6_removedSnpMissPostImp"
outputSNPfile <- "4_6_snpRemovedMissPostImp.txt"
removedSnpMissPostImp(plink, inputPrefix, missCutoff, 
					  outputSNPfile, outputPrefix)

   
setwd("..")
  


############################################################
### code chunk number 5: Data subset and expansion 
############################################################


  
inputPrefix <- "4_6_removedSnpMissPostImp"
inputOriginalQCed <- "3_1_liftedDataset"

reducedToSpecificfn <- "5_1_reducedToSpecific"
extSpecificDiffAllelefn <- "5_2_extSpecificDiffAllele"
extSpecificMissPosfn <- "5_3_extSpecificMissPos"
extSpecificDiffPosfn <- "5_4_extSpecificDiffPos" 
 

## go to the directory/ dataset name
## imputed dataset
system(paste0("scp ./4-imputation/", inputPrefix, ".* ./5-reductAndExpand/")) 

## original but QC-ed dataset
system(paste0("scp ./3-lifting/", inputOriginalQCed, ".* ./5-reductAndExpand/ "))

system(paste0("scp ./3-lifting/3_4_snpImpRefAlleles.txt ./5-reductAndExpand/"))
system(paste0("scp ./3-lifting/3_4_snpDiffAlleles.txt ./5-reductAndExpand/"))
system(paste0("scp ./3-lifting/3_3_snpMissPos.txt ./5-reductAndExpand/"))
system(paste0("scp ./3-lifting/3_2_snpSameNameDiffPos.txt ./5-reductAndExpand/"))


setwd("5-reductAndExpand/")
# 1. Reduce the imputed dataset to the SNPs before imputation. 
system(paste0(plink, " --bfile ", inputPrefix, 
	   " --extract 3_4_snpImpRefAlleles.txt --make-bed --out ", 
	   reducedToSpecificfn)) 

# 2. Add the SNPs with different alleles with their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_4_snpDiffAlleles.txt")) == 0 ){ 
 	renamePlinkBFile(inputPrefix=reducedToSpecificfn, 
				     outputPrefix=extSpecificDiffAllelefn, action="copy")   
} else { 
	system(paste0(plink, " --bfile ", inputOriginalQCed, 
		   " --extract 3_4_snpDiffAlleles.txt --make-bed --out tmp")) 
	system(paste0(plink, " --bfile ", reducedToSpecificfn, 
		   " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", 
		   extSpecificDiffAllelefn) ) 
	system("rm tmp.*")
  } 
  	

# 3. Add the SNPs with missing positions with their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_3_snpMissPos.txt")) == 0 ){  
	renamePlinkBFile(inputPrefix=extSpecificDiffAllelefn, 
					 outputPrefix=extSpecificMissPosfn, action="copy")  
} else {
	system(paste0(plink, " --bfile ", inputOriginalQCed, 
		   " --extract 3_3_snpMissPos.txt --make-bed --out tmp")) 
	system(paste0(plink, " --bfile ", extSpecificDiffAllelefn, 
		   " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", 
		   extSpecificMissPosfn) ) 
	system("rm tmp.*")
  }
  	
# 4. Add the SNPs with different positions by their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_2_snpSameNameDiffPos.txt")) == 0 ){  
	renamePlinkBFile(inputPrefix=extSpecificMissPosfn, 
				     outputPrefix=extSpecificDiffPosfn, action="copy") 
} else {
	system(paste0(plink, " --bfile ", inputOriginalQCed, 
		   " --extract 3_2_snpSameNameDiffPos.txt --make-bed --out tmp")) 
	system(paste0(plink, " --bfile ", extSpecificMissPosfn, 
		   " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", 
		   extSpecificDiffPosfn) ) 
	system("rm tmp.*")
  }
 

system(paste0("rm ", inputPrefix, "* ", inputOriginalQCed, "*"))  
system(paste0("rm  *.txt *.log *.hh")) 
setwd("..")
 

 


############################################################
### code chunk number 6: Final result
############################################################

## imputed dataset
system(paste0("scp ./1-conversion/1_01_metaData.txt ./6-finalResults/metaData.txt " ))
 
## go to the directory/ dataset name 
system(paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.* ./6-finalResults/"))  
system(paste0("scp ./5-reductAndExpand/5_4_extSpecificDiffPos.* ./6-finalResults/"))

setwd("./6-finalResults/")
renamePlinkBFile(inputPrefix="4_6_removedSnpMissPostImp.", 
				 outputPrefix="imputedSnpsDataset", action="move")

renamePlinkBFile(inputPrefix="5_4_extSpecificDiffPos.", 
				 outputPrefix="specificSnpsDataset", action="move")
 


############################################################
### code chunk number 6: Extending pipeline
############################################################

## additional configuration files for Genipe
impRefDir <- "/data/noether/dataProcessResults/10_Common/imputeReference/1000G_Phase3_2014/"
fastaFile <- "/data/noether/dataProcessResults/10_Common/imputeReference/hg19/hg19.fasta"

 
## Impute genotypes using Genipe
# chrs <- "autosomes"
# chrs <- 23
chrs <- 22
inputPrefix <- "/data/noether/datatmp-nobackup/3_1_BrainCloud/v4/4-imputation/tmp4genipe/3_4_removedSnpDiffAlleles"
thread4impute2 <- 20 ## tune by yourself
thread4shapeit <- 30
segmentSize <- 3000000
imputedByGenipe(chrs, impRefDir, inputPrefix, shapeit, impute2, 
			    plink, fastaFile, segmentSize, thread4impute2, thread4shapeit) 


## merge chunked genomic imputed results
## example
chr <- 2 
inputImpute2 <- "chr2.33000001_36000000.impute2"
probability <- 0.9
completionRate <- 0.98
# info <- 0.6
outputPrefix <- paste0("imputedChr", chr)
mergeByGenipe(inputImpute2, chr, probability, completionRate, info, outputPrefix)
 
 
## extract imputed markers using Genipe 
chr <- 3
inputImpute2 <- paste0("chr", chr,".imputed.impute2")
inputMAP <- paste0("chr", chr,".imputed.map")
format <- "bed"
prob <- 0.9
outputPrefix <- paste0("imputedChr", chr)  
extractByGenipe(inputImpute2, inputMAP, outputPrefix, format, prob)
   
 