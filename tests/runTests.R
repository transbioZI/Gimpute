
# File   : runTests.R
# Author : Junfang Chen
# Version0: 28 Jun 2016
# VersionX: 31 Jan 2018
 


############################################################
### code chunk number 0: create all necessary directories   
############################################################
# make the following directories and sub-directories if they are not created yet.
# system('mkdir 0-rawData')
# system('mkdir 1-conversion')
# system('mkdir 2-QC')
# system('mkdir 3-lifting')
# system('mkdir 4-imputation')
# system('mkdir 5-reductAndExpand')
# system('mkdir 6-finalResults')

# ## go to ./0-rawData
# setwd('0-rawData')
# system('mkdir plinkFiles') ## Original plink files 
# system('mkdir sampleInfo') ## Meta data information
# ## Now: you want to put the necessary PLINK files and meta data in the above two folders.
# setwd('..')



## Define the directory where you place the required/downloaded tools and additional scripts (perl scripts).
GimputeDIR = "/home/junfang.chen/Gimpute/"

## required tools
plink = paste0(GimputeDIR,"tools/plink")
gcta = paste0(GimputeDIR,"tools/gcta64") 
shapeit = paste0(GimputeDIR,"tools/shapeit") 
impute2 = paste0(GimputeDIR,"tools/impute2")
gtool = paste0(GimputeDIR,"tools/gtool")

## required perl scripts // chmod u+x * for the permission
perl4snpsWithSamePos = paste0(GimputeDIR,"extdata/inst/snpsWithSamePos.pl")
alignScript4compareName = paste0(GimputeDIR,"extdata/inst/dataSetImpRefCompareName.pl ")
alignScript = paste0(GimputeDIR,"extdata/inst/dataSetImpRefComparePos.pl ")

## required global libraries/tools/configuration files/.
 

## required libraries  
library(Gimpute)  
library(doParallel) ## parallel computing
library(lattice)    ## PCA plot

## configuration directory: where you put all your files in this folder for the configuration.
chipAnnoFile= "/home/junfang.chen/Gimpute/config/Illumina/Human1M-Duov3_B-b37.Illmn.strand"
dupSampleIDFile = NULL ## if it's not NULL, then the file stores the duplicated sample IDs should be placed in this directory.
excludedProbeIdsFile = NULL ## if it's not NULL, then the file stores the probes which have to be excluded should be placed in this directory.

## The directory where you specify your imputation reference files:
impRefDIR = "/data/noether/datatmp-nobackup/tb2refDatabase/imputeRef/1000Gphase1/"
 
## global parameters/variables
ancestrySymbol = NULL
caseControl = FALSE 



############################################################
### code chunk number 1: SNP information update  
############################################################

## pre-defined parameters
# rawPlinkFiles: the prefix of input plink file name



########################### go to main directory
 
## step 1
## copy plink files and meta information file
rawPlinkFiles = 'snp121'
system( paste0("scp ./0-rawData/plinkFiles/", rawPlinkFiles, ".*  ./1-conversion/") ) 
system( paste0("scp ./0-rawData/sampleInfo/1_01_metaData.txt ./1-conversion/") ) 
setwd("./1-conversion/") 
  
 
## step 2   
## dupSampleIDFile; defined as global variable
inputPrefix = rawPlinkFiles
outputPrefix = "1_02_removedExclInst" 
removeDupID(plink, dupSampleIDFile, inputPrefix, outputPrefix)
 

# step 3 replace group IDs 
metaDataFile = "1_01_metaData.txt"
inputPrefix = "1_02_removedExclInst"
outputPrefix = "1_03_replacedGroupAndSex"
replaceGroupId(plink, inputPrefix, metaDataFile, outputPrefix)

 
# step 4 remove instances without group IDs
metaDataFile = "1_01_metaData.txt"
inputPrefix = "1_03_replacedGroupAndSex"
outputPrefix = "1_04_removedNoGroupId"
removeNoGroupId(plink, inputPrefix, outputPrefix)

## step5 remove instances with improper ancestry 
metaDataFile = "1_01_metaData.txt" 
inputPrefix = "1_04_removedNoGroupId"
outputPrefix = "1_05_removedWrongAnceInst"
removedWrongAnceInst(plink, inputPrefix, metaDataFile, ancestrySymbol, outputPrefix)

## step 6 
inputPrefix = "1_05_removedWrongAnceInst"
excludedProbeIdsFile = excludedProbeIdsFile  ## excludedProbeIds is defined in svn
outputPrefix = "1_06_removedExclProbe" ## 
removedExclProbe(plink, inputPrefix, excludedProbeIdsFile, outputPrefix) 

 
## step 7 
inputPrefix = "1_06_removedExclProbe"
chipAnnoFile = chipAnnoFile ## defined in SVN
chipType = 'illumina'
outputPrefix = "1_07_removedUnmapProbes"   
outputSNPunmapFile = "1_07_probesUnmapped2ChipRef.txt"
removedUnmapProbes(plink, inputPrefix, chipAnnoFile, outputPrefix, outputSNPunmapFile)
 

## step 8 
inputPrefix = "1_07_removedUnmapProbes" 
chipAnnoFile = chipAnnoFile ## defined in SVN
chipType = 'illumina'
outputSNPdupFile = "1_08_probesDouble.txt"
outputPrefix = "1_08_removedDoubleProbes"   
removedDoubleProbes(plink, inputPrefix, chipAnnoFile, chipType, outputSNPdupFile, outputPrefix)


## step 9
inputPrefix = "1_08_removedDoubleProbes" 
chipAnnoFile = chipAnnoFile ## defined in SVN
chipType = 'illumina'

outputPrefix = "1_09_updatedSnpInfo"   
updatedSnpInfo(plink, inputPrefix,  chipAnnoFile, chipType, outputPrefix)


## step 10 
inputPrefix = "1_09_updatedSnpInfo"
outputPrefix = "1_10_changedXyChr"
changedXyChr(plink, inputPrefix, outputPrefix)


## step 11 
inputPrefix = "1_10_changedXyChr"
outputPrefix = "1_11_removedYMtSnp"
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
inputPrefix4QC = "1_11_removedYMtSnp"
system( paste0("scp ./1-conversion/", inputPrefix4QC, ".*", " ./2-QC/") )
setwd("./2-QC/")

## step 1 
inputPrefix = inputPrefix4QC
hhCutOff = 0.005 ##  can be tuned
outputPrefix = "2_01_removedSnpHetX" 
outputFile_SNPhhFreqAll = "2_01_snpHHfreqAll.txt"
outputFile_SNPhhFreqRetained = "2_01_snpHHfreqRetained.txt"
removedSnpHetX(plink, inputPrefix, hhCutOff, outputPrefix, outputFile_SNPhhFreqAll, outputFile_SNPhhFreqRetained)


## step 2  2_02_removedHetXInst

inputPrefix = "2_01_removedSnpHetX"
hhSubjCutOff = 15
outputPrefix = "2_02_removedInstHetX"
outputFile_subjHetFreqAll = "2_02_instHetXfreqAll.txt" 
outputFile_subjHetFreqRetained = "2_02_instHetXfreqRetained.txt"  
outputFile_SNPhhFreqAll = "2_02_snpHHfreqAll.txt"
removedMaleHetX(plink, inputPrefix, hhSubjCutOff, outputPrefix, outputFile_subjHetFreqAll, outputFile_subjHetFreqRetained, outputFile_SNPhhFreqAll) 

## step 3 
# 3. Set all heterozygous alleles of SNPs of the chromosome 23 for males
inputPrefix = "2_02_removedInstHetX" 
outputPrefix = "2_03_setHeteroHaploMissing" 
setHeteroHaploMissing(plink, inputPrefix, outputPrefix)
 

## step 4   SNP missingness < 0.05 (before sample removal);  
inputPrefix = "2_03_setHeteroHaploMissing" 
snpMissCutOff = 0.05 #
outputPrefix = "2_04_removedSnpMissPre" 
removedSnpMiss(plink, snpMissCutOff, inputPrefix, outputPrefix)

## step 5 
# subject missingness < 0.02; 
inputPrefix = "2_04_removedSnpMissPre" 
sampleMissCutOff = 0.02
outputPrefix = "2_05_removedInstMiss" 
removedInstMiss(plink, sampleMissCutOff, inputPrefix, outputPrefix)
 
## step 6 
inputPrefix = "2_05_removedInstMiss" 
Fhet = 0.2 #
outputPrefix = "2_06_removedInstFhet" 
removedInstFhet(plink, Fhet, inputPrefix, outputPrefix)
  	

## step 7
inputPrefix = "2_06_removedInstFhet"  
outputPrefix = "2_07_removedParentIdsMiss" 
removedParentIdsMiss(plink, inputPrefix, outputPrefix)


## step 8 
snpMissCutOff = 0.02
inputPrefix = "2_07_removedParentIdsMiss"  
outputPrefix = "2_08_removedSnpMissPost" 
removedSnpMiss(plink, snpMissCutOff, inputPrefix, outputPrefix)


## step 9  ## caseControl
 ##  Remove SNPs with difference >= 0.02 of SNP missingness between cases and controls.
inputPrefix = "2_08_removedSnpMissPost"  
outputPrefix = "2_09_removedSnpMissDiff" 
snpMissDifCutOff = 0.02  
removedSnpMissDiff(plink, inputPrefix, snpMissDifCutOff, outputPrefix, caseControl) 

 
## step 10
femaleChrXmissCutoff = 0.05
inputPrefix = "2_09_removedSnpMissDiff"  
outputPrefix = "2_10_removedSnpFemaleChrXmiss" 
removedSnpFemaleChrXmiss(plink, femaleChrXmissCutoff, inputPrefix, outputPrefix)
  

## step 11

pval = 0.000001
inputPrefix = "2_10_removedSnpFemaleChrXmiss"  
outputFile_pVal = "2_11_snpHwePvalAutoCt.txt" 
outputSNPfile =  "2_11_snpRemovedHweAutoCt.txt" 
outputPrefix = "2_11_removedSnpHweAutoCt" 

removedSnpHWEautoControl(plink, inputPrefix, pval, outputFile_pVal, outputSNPfile, outputPrefix)


## step 12 
pval = 0.000001
inputPrefix = "2_11_removedSnpHweAutoCt"   
outputFile_pVal = "2_12_snpHwePvalfemaleXct.txt" 
outputSNPfile = "2_12_snpRemovedHweFemaleXct.txt" 
outputPrefix = "2_12_removedSnpHweFemaleXct" 
removedSnpFemaleChrXhweControl(plink, inputPrefix, pval, outputFile_pVal, outputSNPfile, outputPrefix)
 
 

## step 13 
pval = 0.000001
inputPrefix = "2_11_removedSnpHweAutoCt"   ## the output from the step before last step. 
outputFile_pVal = "2_13_snpHwePvalfemaleXall.txt" 
outputSNPfile = "2_13_snpRemovedHweFemaleXall.txt" 
outputPrefix = "2_13_removedSnpHweFemaleXall" 
removedSnpFemaleChrXhweAll(plink, inputPrefix, pval, outputFile_pVal, outputSNPfile, outputPrefix)
  

 
################################################
## in order to get the ethnic group info
setwd('..')
metaDataFile = "1_01_metaData.txt"
system( paste0("scp ./1-conversion/", metaDataFile, " ./2-QC/") )
setwd("./2-QC/")
################################################ 

## step 14  

inputPrefix = "2_13_removedSnpHweFemaleXall" ## the output from step 12
outputPC4subjFile = "2_14_eigenvalAfterQC.txt"
outputPCplotFile = "2_14_eigenvalAfterQC.png"
plotPCA4plink(gcta, inputPrefix, outputPC4subjFile, outputPCplotFile)

 
## remove outliers
inputPrefix = "2_13_removedSnpHweFemaleXall"
cutoff =  NULL 
cutoffSign = 'greater' ## not used if cutoff==NULL


inputPC4subjFile = "2_14_eigenvalAfterQC.txt"
outputPC4outlierFile = "2_14_eigenval4outliers.txt"
outputPCplotFile = "2_14_removedOutliers.png"
outputPrefix = "2_14_removedOutliers" 

removeOutlierByPCs(plink, gcta, inputPrefix, cutoff, cutoffSign, inputPC4subjFile, outputPC4outlierFile, outputPCplotFile, outputPrefix)
 
 

 
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

 
## step 1  
system( paste0("cp ./2-QC/2_12_removedSnpHweFemaleXct.bed ./3-lifting/3_1_liftedDataset.bed") )
system( paste0("cp ./2-QC/2_12_removedSnpHweFemaleXct.fam ./3-lifting/3_1_liftedDataset.fam") )
system( paste0("cp ./2-QC/2_12_removedSnpHweFemaleXct.bim ./3-lifting/3_1_liftedDataset.bim") )
setwd("./3-lifting/")

 
 
# 2. Remove SNPs for which the name (rs-ID) has a different position (i.e. combination of
# base pair position and chromosome) in the imputation reference files. 

inputPrefix = "3_1_liftedDataset" 
out2.snp = "3_2_snpDiffNamePos"
out2 = "3_2_removedSnpDiffNamePos"
 
inputBIMfn = paste0(inputPrefix, ".bim ") ## bim file from step1 
system( paste0(alignScript4compareName, inputBIMfn, impRefDIR, " pos > ", out2.snp, ".txt"))
system( paste0(plink, " --bfile ",  inputPrefix, " --exclude ", out2.snp, ".txt --make-bed --out ", out2) )  


# 3. Remove SNPs which bp and chr position are not contained in the imputation reference files.

out3 = "3_3_removedSnpMissPos"
out3.snp = "3_3_snpMissPos"

bim2fn = paste0(out2, ".bim ")  
system( paste0(alignScript, bim2fn, impRefDIR, " pos > ", out3.snp, ".txt"))
system( paste0(plink, " --bfile ",  out2, " --exclude ", paste0(out3.snp, ".txt"), " --make-bed --out ", out3) ) 



## 4. Remove SNPs which have an allele which is not in the imputation reference files for that SNP.
out4 = "3_4_removedSnpDiffAlleles"
out4.snp = "3_4_snpDiffAlleles"
out4.snpRetained = "3_4_snpImpRefAlleles"
 
bim3fn = paste0(out3, ".bim ")  
system( paste0(alignScript, bim3fn, impRefDIR, " alleles > ", out4.snp, ".txt"))
system( paste0(plink, " --bfile ",  out3, " --exclude ", paste0(out4.snp, ".txt"), " --make-bed --out ", out4) )  

##  Remove SNPs which have an allele which is not in the imputation reference files for that SNP.
snpDifAllele = read.table(paste0(out4.snp, ".txt"), stringsAsFactors=F)
snpDifAllele = snpDifAllele[,1] 
inputBIM = read.table(paste0(out3, ".bim"), stringsAsFactors=F)
snpRetained = setdiff(inputBIM[,2], snpDifAllele)
write.table(snpRetained, file=paste0(out4.snpRetained, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")

 



system( paste0("rm  *.log *.hh") ) 

setwd("..")





############################################################
### code chunk number 4: Imputation
############################################################

## step 1 
## Remove monomorphic SNPs from lifted/QC-ed data 

inputPrefix4aligned2impRef = "3_4_removedSnpDiffAlleles" ## will also be used in step 4 and 5;
outputPrefix = "4_1_removedMonoSnp"
outputMonoSNPfile = "4_1_snpMonoRemoved.txt" # will be used in step 4 and 5.

## copy plink files from last step; 
system( paste0("cp ./3-lifting/", inputPrefix4aligned2impRef, ".bed  ./4-imputation/") )
system( paste0("cp ./3-lifting/", inputPrefix4aligned2impRef, ".fam  ./4-imputation/") )
system( paste0("cp ./3-lifting/", inputPrefix4aligned2impRef, ".bim  ./4-imputation/") )

## remove Monomorphic SNPs
setwd('4-imputation')
removedMonoSnp(plink, inputPrefix4aligned2impRef, outputPrefix, outputMonoSNPfile)

# ## remove unwanted plink files << ## will also be used in step 4 and 5; after that you can remove.
# system(paste0('rm ', inputPrefix, '*'))


## step 2 
##########################################################################
########################################################################## imputation main pipeline

## sub-steps
# step 2.1 chrWiseSplit.R
# step 2.2 chunk4eachChr.R
# step 2.3 prePhasingByShapeit.R
# step 2.4 imputedByImpute2.R 
# step 2.5 formatConvertGtool.R 
# step 2.6 mergeImputeData.R 
# step 2.7 filterImputeData.R 
 
## One must create directories for storing temporal imputation output files 
## The name of these directories must be fixed for the sake of the subsequent steps.
tmp4imputeDIR = 'tmpImpute'
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
outputPrefix = "4_1_removedMonoSnp"
outputPrefixTmp  = 'gwas_data_chr'  

system( paste0("scp ", outputPrefix, ".bed  ./", tmp4imputeDIR, "/1-dataFiles/", outputPrefixTmp, ".bed"))
system( paste0("scp ", outputPrefix, ".fam  ./", tmp4imputeDIR, "/1-dataFiles/", outputPrefixTmp, ".fam"))
system( paste0("scp ", outputPrefix, ".bim  ./", tmp4imputeDIR, "/1-dataFiles/", outputPrefixTmp, ".bim"))
  

setwd(paste0('./', tmp4imputeDIR, '/1-dataFiles/'))
inputPrefix  = outputPrefixTmp 
outputPrefix  = outputPrefixTmp ## Same as the input files; but the chromosome number will be appended
chrX_PAR1suffix = 'X_PAR1'
chrX_PAR2suffix = 'X_PAR2'
PAR = chrWiseSplit(plink, inputPrefix, chrX_PAR1suffix, chrX_PAR2suffix)
 
PAR 
 
if (PAR[[1]]) {par1 = 'X_PAR1'} else {par1 = NULL}
if (PAR[[2]]) {par2 = 'X_PAR2'} else {par2 = NULL}

  

########################################
######################################## chunk4eachChr.R
## step 2.2

inputPrefix = 'gwas_data_chr'
outputPrefix = 'chunks_chr' 
chrs = c(1:23, par1, par2)   
windowSize = 3000000 
chunk4eachChr(inputPrefix, outputPrefix, chrs, windowSize) 

setwd('..') 
system( paste0("mv ./1-dataFiles/", outputPrefix, "*.txt  ./2-chunkFile/"))



## step 2.3  

# define directories
dataDIR = "./1-dataFiles/"
impRefDIR = "/data/noether/datatmp-nobackup/tb2refDatabase/imputeRef/1000Gphase1/"
phaseDIR = "./3-phaseResults/"
 
prefix4plinkEachChr = 'gwas_data_chr'  
# chrs = c(1, 11, 23, par1, par2)   
nThread = 40
effectiveSize = 20000 
nCore = 1 
prePhasingByShapeit(shapeit, chrs, dataDIR, prefix4plinkEachChr, impRefDIR, phaseDIR, nThread, effectiveSize, nCore)

 

# step 2.4   
# define directories
prefixChunk = "./2-chunkFile/chunks_chr"  ## 
phaseDIR = "./3-phaseResults/"
impRefDIR = "/data/noether/datatmp-nobackup/tb2refDatabase/imputeRef/1000Gphase1/"
imputedDIR = "./4-imputeResults/"

prefix4plinkEachChr = 'gwas_data_chr'
# chrs = c(1,23, par1, par2)     
nCore = 30 ## try to tune it for your own data size
effectiveSize = 20000  

imputedByImpute2(impute2, chrs, prefixChunk, phaseDIR, impRefDIR, imputedDIR, prefix4plinkEachChr, nCore, effectiveSize)


  
 
# step 2.5   
####################################################### extract only SNPs (without INDELs)
setwd("./4-imputeResults")  
## extract only SNPs starting with 'rs';  .
ls = system("ls gwas*.impute2", intern=T)
variantPrefix = 'rs' 

biglists = as.list(ls)
mclapply(biglists, function(i){
	arg1 = paste0(i, "noINDEL.impute2")
	arg2 = paste0("grep '", variantPrefix, "' ", i, " | awk '{if(length($4)==1 && length($5)==1) print}' > ", arg1)
	# print(arg2)
	system(arg2)
}, mc.cores= 38) 
setwd("..")  
####################################################### <<<


# chrs = c(1,23, par1, par2)
   
prefixChunk = "./2-chunkFile/chunks_chr"  ## 
phaseDIR = "./3-phaseResults/" 
imputedDIR = "./4-imputeResults/"
prefix4plinkEachChr = 'gwas_data_chr'
suffix4imputed = ".impute2noINDEL.impute2"
# suffix4imputed = ".impute2"
postImputeDIR = './5-postImpute/' 

nCore = 30   
formatConvertGtool(gtool, chrs, prefixChunk, phaseDIR, imputedDIR, prefix4plinkEachChr, suffix4imputed, postImputeDIR, nCore)

 


# step 2.6  
####################################################### Modify missing genotype format.
setwd("./5-postImpute/")  
# replace 'N' in the .ped files into 0 > missing values.
chrslist = as.list(chrs)
prefix4plinkEachChr = 'gwas_data_chr' ## just for parallel computing
fn = mclapply(chrslist, function(i){
	system(paste0("sed -i 's/N/0/g' ", prefix4plinkEachChr, i, ".*ped "))
}, mc.cores=nCore)
## Note: check if also any "N" in *.fam files. If so, change back after merging. 
####################################################### <<< 

prefix4plinkEachChr = 'gwas_data_chr'
prefix4imputedPlink = 'gwasImputed'
par1 = 'X_PAR1'
par2 = 'X_PAR2'
# chrs = c(1,23, par1, par2)   
mergeImputeData(plink, chrs, prefix4plinkEachChr, prefix4imputedPlink, nCore)
 
 

#######################################################
setwd("..")
system("mv ./5-postImpute/gwasImputed* ./6-finalResults/") 
system("mv ./4-imputeResults/*.impute2_info ./6-finalResults/") 

# step 2.7    
setwd("./6-finalResults")
suffix4impute2info = ".impute2_info"
impute2infoFile = "impute2infoUpdated.txt"
infoScore = 0.6
badImputeSNPfile = "badImputeSNPs.txt"
prefix4imputedPlink = "gwasImputed"
prefix4imputedFilterPlink = "gwasImputedFiltered" 
filterImputeData(plink, suffix4impute2info, impute2infoFile, infoScore, badImputeSNPfile, prefix4imputedPlink, prefix4imputedFilterPlink)

setwd("..")
setwd("..")

########################################################################## After imputation

## step 2 
## Final imputed results >> 
 
prefix4imputedPlink = "gwasImputed"
## output file name change to: 
imputedDatasetfn = "4_2_imputedDataset"
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedPlink, ".bed ", imputedDatasetfn, ".bed") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedPlink, ".bim ", imputedDatasetfn, ".bim") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedPlink, ".fam ", imputedDatasetfn, ".fam") )


## step 3
## Filtered imputed data set; Remove imputed SNPs with (info < 0.6), only retain 'Good' SNPs.
prefix4imputedFilterPlink = "gwasImputedFiltered"
filteredImputedDatasetfn = "4_3_removedSnpInfoPostImp" 
snpWithBadInfoFile = "4_3_snpRemovedInfoPostImp.txt"
snpImputedInfoScoreFile = "4_3_snpImputedInfoScore.txt"
 
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedFilterPlink, ".bed ", filteredImputedDatasetfn, ".bed") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedFilterPlink, ".bim ", filteredImputedDatasetfn, ".bim") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedFilterPlink, ".fam ", filteredImputedDatasetfn, ".fam") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", impute2infoFile, " ", snpImputedInfoScoreFile) )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", badImputeSNPfile, " ", snpWithBadInfoFile) )


## step 4 
## Remove previous identified monomorphic SNPs in the imputed dataset.
filteredImputedDatasetfn = "4_3_removedSnpInfoPostImp" 
removedMonoSnpAfter = "4_4_removedMonoSnpAfter"
# perl4snpsWithSamePos = "/home/junfang.chen/groupSVN/trunk/datasets/bin/dataProcess/4-imputation/snpsWithSamePos.pl"

## if no monomorphic SNPs:
if ( file.size(paste0(outputMonoSNPfile))==0 ){ 
		system(paste0("scp ", filteredImputedDatasetfn, ".bed ", removedMonoSnpAfter, ".bed") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".bim ", removedMonoSnpAfter, ".bim") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".fam ", removedMonoSnpAfter, ".fam") )  
} else { 
	## extract PLINK files contain only monomorphic SNPs from the original aligned (lifted and QC-ed) data set.
	system( paste0(plink, " --bfile ", inputPrefix4aligned2impRef, " --extract ", outputMonoSNPfile, " --make-bed --out ", inputPrefix4aligned2impRef, "Tmp") ) 
	system( paste0(perl4snpsWithSamePos, " ", filteredImputedDatasetfn, ".bim ", inputPrefix4aligned2impRef, "Tmp", ".bim > tmp.txt" ))
	system( paste0(plink, " --bfile ", filteredImputedDatasetfn, " --exclude tmp.txt --make-bed --out ", removedMonoSnpAfter) )
} 


## step 5
## Add previous identified monomorphic SNPs in the imputed dataset.
addedMonoSnpAfter = "4_5_addedMonoSnpAfter" 

 ## if no monomorphic SNPs:
if ( file.size(paste0(outputMonoSNPfile))==0 ){ 
		system(paste0("scp ", filteredImputedDatasetfn, ".bed ", addedMonoSnpAfter, ".bed") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".bim ", addedMonoSnpAfter, ".bim") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".fam ", addedMonoSnpAfter, ".fam") )  
} else { 
  
	## merge both datasets
	system( paste0(plink, " --bfile ", removedMonoSnpAfter, " --bmerge ", 
		inputPrefix4aligned2impRef, ".bed ", inputPrefix4aligned2impRef, ".bim ", inputPrefix4aligned2impRef, ".fam --make-bed --out ", addedMonoSnpAfter) )
	## remove tmp files
	# system( paste0("rm tmp.txt") )
	# system( paste0("rm ", inputPrefix4aligned2impRef, "*") )
}  


## step 6
## Remove SNPs which have a non missing value for less then 20 instances. 

inputPrefix = addedMonoSnpAfter  
missCutoff = 20
outputPrefix = "4_6_removedSnpMissPostImp"
snpWithManyMissSNPfile = "4_6_snpRemovedMissPostImp.txt"
removedSnpMissPostImp(plink, inputPrefix, missCutoff, snpWithManyMissSNPfile, outputPrefix)

   
setwd("..")
  


############################################################
### code chunk number 5: Data subset and expansion 
############################################################


  
inputPrefix = "4_6_removedSnpMissPostImp"
inputOriginalQCed = "3_1_liftedDataset"

reducedToSpecificfn = "5_1_reducedToSpecific"
extSpecificDiffAllelefn = "5_2_extSpecificDiffAllele"
extSpecificMissPosfn = "5_3_extSpecificMissPos"
extSpecificDiffPosfn = "5_4_extSpecificDiffPos" 
 

## go to the directory/ dataset name
## imputed dataset
system( paste0("scp ./4-imputation/", inputPrefix, ".bed ./5-reductAndExpand/ " ) )
system( paste0("scp ./4-imputation/", inputPrefix, ".fam ./5-reductAndExpand/ " ) )
system( paste0("scp ./4-imputation/", inputPrefix, ".bim ./5-reductAndExpand/ " ) )

## original but QC-ed dataset
system( paste0("scp ./3-lifting/", inputOriginalQCed, ".bed ./5-reductAndExpand/ " ) )
system( paste0("scp ./3-lifting/", inputOriginalQCed, ".fam ./5-reductAndExpand/ " ) )
system( paste0("scp ./3-lifting/", inputOriginalQCed, ".bim ./5-reductAndExpand/ " ) )
## 
system( paste0("scp ./3-lifting/3_4_snpImpRefAlleles.txt ./5-reductAndExpand/ ") )
system( paste0("scp ./3-lifting/3_4_snpDiffAlleles.txt ./5-reductAndExpand/ ") )
system( paste0("scp ./3-lifting/3_3_snpMissPos.txt ./5-reductAndExpand/ ") )
system( paste0("scp ./3-lifting/3_2_snpDiffNamePos.txt ./5-reductAndExpand/ ") )


setwd("5-reductAndExpand/")
# 1. Reduce the imputed dataset to the SNPs before imputation. 
system(paste0(plink, " --bfile ", inputPrefix, " --extract 3_4_snpImpRefAlleles.txt --make-bed --out ", reducedToSpecificfn) ) 

# 2. Add the SNPs with different alleles with their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_4_snpDiffAlleles.txt"))==0 ){  
	system(paste0("scp ", reducedToSpecificfn, ".bed ", extSpecificDiffAllelefn, ".bed") ) 
	system(paste0("scp ", reducedToSpecificfn, ".bim ", extSpecificDiffAllelefn, ".bim") ) 
	system(paste0("scp ", reducedToSpecificfn, ".fam ", extSpecificDiffAllelefn, ".fam") ) 
} else { 
	system(paste0(plink, " --bfile ", inputOriginalQCed, " --extract 3_4_snpDiffAlleles.txt --make-bed --out tmp") ) 
	system(paste0(plink, " --bfile ", reducedToSpecificfn, " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", extSpecificDiffAllelefn) ) 
	system("rm tmp.*")
  } 
  	

# 3. Add the SNPs with missing positions with their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_3_snpMissPos.txt"))==0 ){  
	system(paste0("scp ", extSpecificDiffAllelefn, ".bed ", extSpecificMissPosfn, ".bed") ) 
	system(paste0("scp ", extSpecificDiffAllelefn, ".bim ", extSpecificMissPosfn, ".bim") ) 
	system(paste0("scp ", extSpecificDiffAllelefn, ".fam ", extSpecificMissPosfn, ".fam") ) 
} else {
	system(paste0(plink, " --bfile ", inputOriginalQCed, " --extract 3_3_snpMissPos.txt --make-bed --out tmp") ) 
	system(paste0(plink, " --bfile ", extSpecificDiffAllelefn, " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", extSpecificMissPosfn) ) 
	system("rm tmp.*")
  }
  	
# 4. Add the SNPs with different positions by their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_2_snpDiffNamePos.txt"))==0 ){  
	system(paste0("scp ", extSpecificMissPosfn, ".bed ", extSpecificDiffPosfn, ".bed") ) 
	system(paste0("scp ", extSpecificMissPosfn, ".bim ", extSpecificDiffPosfn, ".bim") ) 
	system(paste0("scp ", extSpecificMissPosfn, ".fam ", extSpecificDiffPosfn, ".fam") ) 
} else {
	system(paste0(plink, " --bfile ", inputOriginalQCed, " --extract 3_2_snpDiffNamePos.txt --make-bed --out tmp") ) 
	system(paste0(plink, " --bfile ", extSpecificMissPosfn, " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", extSpecificDiffPosfn) ) 
	system("rm tmp.*")
  }
 

system( paste0("rm ", inputPrefix, "*") ) 
system( paste0("rm ", inputOriginalQCed, "*") ) 
 
 
system( paste0("rm  *.txt *.log *.hh") ) 
setwd("..")
 

 


############################################################
### code chunk number 6: Final result
############################################################

## imputed dataset
system( paste0("scp ./1-conversion/1_01_metaData.txt ./6-finalResults/metaData.txt " ) )
 
## go to the directory/ dataset name 
system( paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.bed ./6-finalResults/imputedSnpsDataset.bed") )
system( paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.bim ./6-finalResults/imputedSnpsDataset.bim") )
system( paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.fam ./6-finalResults/imputedSnpsDataset.fam") )

system( paste0("scp ./5-reductAndExpand/5_4_extSpecificDiffPos.bed ./6-finalResults/specificSnpsDataset.bed") )
system( paste0("scp ./5-reductAndExpand/5_4_extSpecificDiffPos.bim ./6-finalResults/specificSnpsDataset.bim") )
system( paste0("scp ./5-reductAndExpand/5_4_extSpecificDiffPos.fam ./6-finalResults/specificSnpsDataset.fam") )




############################################################
### code chunk number 6: Extending pipeline
############################################################

## additional configuration files for Genipe
impRefDIR4genipe = "/data/noether/dataProcessResults/10_Common/imputeReference/1000G_Phase3_2014/"
fastaFile = "/data/noether/dataProcessResults/10_Common/imputeReference/hg19/hg19.fasta"



## Impute genotypes using Genipe
# chrs ="autosomes"
# chrs =23
chrs =22
inputPrefix = "/data/noether/datatmp-nobackup/3_1_BrainCloud/v4/4-imputation/tmp4genipe/3_4_removedSnpDiffAlleles"
thread4impute2 = 20 ## tune by yourself
thread4shapeit = 30
segmentSize = 3000000
imputedByGenipe(chrs, impRefDIR4genipe, inputPrefix, shapeit, impute2, plink, fastaFile, segmentSize, thread4impute2, thread4shapeit) 


## merge chunked genomic imputed results
## example
chr = 2 
inputImpute2 = 'chr2.33000001_36000000.impute2'
probability = 0.9
completionRate = 0.98
# info = 0.6
outputPrefix = paste0('imputedChr', chr)
mergeByGenipe(inputImpute2, chr, probability, completionRate, info, outputPrefix)
 
 
## extract imputed markers using Genipe 
chr = 3
inputImpute2 = paste0('chr', chr,'.imputed.impute2')
inputMAP = paste0('chr', chr,'.imputed.map')
format = 'bed'
prob = 0.9
outputPrefix = paste0('imputedChr', chr)  
extractByGenipe(inputImpute2, inputMAP, outputPrefix, format, prob)
   
 