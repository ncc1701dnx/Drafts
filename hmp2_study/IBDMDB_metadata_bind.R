getwd()
# Load the needed packages
# mzR and xcms and edgeR are indispensable
# Write date 04 July ~ 05 July
library(stats) # In adjustRtime
library(Spectra)
library(mzR)
library(ggplot2)
library(RaMS)
#library(readMzXmlData)#use to quickly read mzXML files if not use xcms
library(xcms) # LC-MS data pre-processing 
#library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(MsFeatures)
library(rlang) # dependency for dplyr
library(dplyr)
#library(usethis) # dependency for devtool ignore after xcmsViewer installed
#library(devtools) # use for install xcmsViewer
#library(omicsViewer) # mandatory for xcmsViewer
#library(xcmsViewer) # MS data annotation
#library(parallel) # required for annoMS1 function from xcmsViewer
#library(BiocParallel) # required for annoMS1 function from xcmsViewer
library(MAIT) # dependency for xcmsViewer
library(MSnbase) # In most case already active in baseR, but just for sure 
library(stringr) # use for subset of stings

IBDMDB_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/hmp2_metadata.csv", header = T)
IBDMDB_c18neg_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/C18n_RawFileInventory.csv", header = T)
IBDMDB_c18pos_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/C8p_RawFileInventory.csv", header = T)
IBDMDB_HILneg_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/HILn_RawFileInventory.csv", header = T)
IBDMDB_HILpos_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/HILp_RawFileInventory.csv", header = T)
#let's check which marker is for patients' cross-sample name
length(unique(IBDMDB_metadata$External.ID))
# in total 2895 uniques
length(unique(IBDMDB_metadata$Participant.ID))
# in total 131 uniques
table(IBDMDB_metadata$data_type)
# which lable can be used as across sample match lable in mbx datasets
length(which(IBDMDB_c18neg_metadata$Sample.ID %in% IBDMDB_c18pos_metadata$Sample.ID))
length(which(IBDMDB_c18neg_metadata$Sample.ID %in% IBDMDB_HILpos_metadata$Sample.ID))
# the Sample.ID seems okay, almost all of them are identical
# 使用Sample.ID作为MBX数据的matching lable
length(unique(IBDMDB_c18neg_metadata$Sample.ID))
strsplit()
length(unique(str_sub(IBDMDB_metadata$External.ID, start = -5)))
which(IBDMDB_metadata$data_type == "metabolomics")
MBX_IBDMDB_metadata <- IBDMDB_metadata[which(IBDMDB_metadata$data_type == "metabolomics"),]
name_metadata <- unlist(sapply(MBX_IBDMDB_metadata$External.ID, function(x){str_sub(x, start = -5)}))
name_MBX <- unlist(sapply(IBDMDB_c18neg_metadata$Sample.ID, function(x){
  strsplit(x, "-")[[1]][2]
}))
table(MBX_IBDMDB_metadata$Age.when.started.smoking)
length(which(name_metadata %in% name_MBX))
class(name_MBX)
class(name_metadata)
colnames(IBDMDB_metadata)
table(MBX_IBDMDB_metadata$Tube.A..Metabolomics)
length(which(MBX_IBDMDB_metadata$Tube.A..Metabolomics %in% IBDMDB_c18pos_metadata$Sample.ID))
# the Tube.A..Metabolomics is the match for metadata and inventory

# 原metadata table实在是太差了，MBX数据里面没有任何的病人用药信息
# 用唯一matchable的ID：Participant.ID来计算总和的病人用药情况，然后插入MBX的数据里
# Firstly is about AZA:
table(IBDMDB_metadata$Azathioprine..Imuran..Azasan.)
#check who never used
aza_contrl <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Azathioprine..Imuran..Azasan. == "Never taken")]
which(MBX_IBDMDB_metadata$Participant.ID %in% aza_contrl)
MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.[which(MBX_IBDMDB_metadata$Participant.ID %in% aza_contrl)] = "Never taken"
#check who used before
aza_prior <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Azathioprine..Imuran..Azasan. == "Taken prior to baseline")]
which(MBX_IBDMDB_metadata$Participant.ID %in% aza_prior)
MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.[which(MBX_IBDMDB_metadata$Participant.ID %in% aza_prior)] = "Taken prior to baseline"
#check who are using
aza_treat <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Azathioprine..Imuran..Azasan. == "Current")]
which(MBX_IBDMDB_metadata$Participant.ID %in% aza_treat)
MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.[which(MBX_IBDMDB_metadata$Participant.ID %in% aza_treat)] = "Current"
table(MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.)
# Then 6MP
table(MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.)
table(IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.)
#check who never used
mp_contrl <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Mercaptopurine..Purinethol..6MP. == "Never taken")]
which(MBX_IBDMDB_metadata$Participant.ID %in% mp_contrl)
MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.[which(MBX_IBDMDB_metadata$Participant.ID %in% mp_contrl)] = "Never taken"
#check who used before
mp_prior <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Mercaptopurine..Purinethol..6MP. == "Taken prior to baseline")]
which(MBX_IBDMDB_metadata$Participant.ID %in% mp_prior)
MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.[which(MBX_IBDMDB_metadata$Participant.ID %in% mp_prior)] = "Taken prior to baseline"
#check who are using
mp_treat <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Mercaptopurine..Purinethol..6MP. == "Current")]
which(MBX_IBDMDB_metadata$Participant.ID %in% mp_treat)
MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.[which(MBX_IBDMDB_metadata$Participant.ID %in% mp_treat)] = "Current"
table(MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.)
# Now save this changed MBX metadata table
write.table(MBX_IBDMDB_metadata, "R_Works/outputs/MBX_IBDMDB_metadata.csv", sep = ",", quote = F,col.names = T, row.names = F)

# Change the inventory table from .raw to .mzXML 
IBDMDB_c18neg_metadata$Raw.file.name <- gsub(".raw", ".mzXML", IBDMDB_c18neg_metadata$Raw.file.name)
IBDMDB_c18pos_metadata$Raw.file.name <- gsub(".raw", ".mzXML", IBDMDB_c18pos_metadata$Raw.file.name)
IBDMDB_HILneg_metadata$Raw.file.name <- gsub(".raw", ".mzXML", IBDMDB_HILneg_metadata$Raw.file.name)
IBDMDB_HILpos_metadata$Raw.file.name <- gsub(".raw", ".mzXML", IBDMDB_HILpos_metadata$Raw.file.name)
# Well we wanna select the UC patients' file first since the mzXML files are too big and conversion is slow
IBDMDB_c18neg_metadata$Raw.file.name <- gsub(".mzXML", ".raw",IBDMDB_c18neg_metadata$Raw.file.name)
IBDMDB_c18pos_metadata$Raw.file.name <- gsub(".mzXML", ".raw",IBDMDB_c18pos_metadata$Raw.file.name)
IBDMDB_HILneg_metadata$Raw.file.name <- gsub(".mzXML", ".raw",IBDMDB_HILneg_metadata$Raw.file.name)
IBDMDB_HILpos_metadata$Raw.file.name <- gsub(".mzXML", ".raw",IBDMDB_HILpos_metadata$Raw.file.name)

# Write a for loop to bind inventory data and the changed MBX metadata togather
bind_metadata_MBXinventory <- function(Inventory, MBX_Metadata){
azaList <- c()
mpList <- c()
diaList <- c()
for (i in 1:nrow(Inventory)) {
  #i = 1 #test
  if(Inventory$Sample.ID[i] %in% MBX_Metadata$Tube.A..Metabolomics){
    tt <- which(MBX_Metadata$Tube.A..Metabolomics == Inventory$Sample.ID[i])
    azaList <- c(azaList, MBX_Metadata$Azathioprine..Imuran..Azasan.[tt])
    mpList <- c(mpList, MBX_Metadata$Mercaptopurine..Purinethol..6MP.[tt])
    diaList <- c(diaList, MBX_Metadata$diagnosis[tt])
  } else {
    print(paste0(Inventory$Sample.ID[i], " not found in metadata"))
    azaList <- c(azaList, NA)
    mpList <- c(mpList, NA)
    diaList <- c(diaList, NA)
  }
}
Inventory$Azathioprine..Imuran..Azasan. =  azaList
Inventory$Mercaptopurine..Purinethol..6MP. = mpList
Inventory$diagnosis = diaList
Inventory
}
bind_metadata_MBXinventory(IBDMDB_c18pos_metadata,MBX_IBDMDB_metadata)
IBDMDB_c18pos_metadata <- bind_metadata_MBXinventory(IBDMDB_c18pos_metadata,MBX_IBDMDB_metadata)
IBDMDB_HILneg_metadata <- bind_metadata_MBXinventory(IBDMDB_HILneg_metadata,MBX_IBDMDB_metadata)
IBDMDB_HILpos_metadata <- bind_metadata_MBXinventory(IBDMDB_HILpos_metadata, MBX_IBDMDB_metadata)
table(MBX_IBDMDB_metadata$diagnosis)
rm("azaList","mpList","diaList","tt")
nrow(MBX_IBDMDB_metadata)
which(IBDMDB_c18neg_metadata$diagnosis == "UC" & IBDMDB_c18neg_metadata$Azathioprine..Imuran..Azasan. == "Current")
which(IBDMDB_c18neg_metadata$diagnosis == "UC" & IBDMDB_c18neg_metadata$Mercaptopurine..Purinethol..6MP. == "Current")

head(IBDMDB_c18neg_metadata$Raw.file.name)

UC_c18neg_metadata <- subset(IBDMDB_c18neg_metadata, diagnosis == "UC")
UC_c18pos_metadata <- subset(IBDMDB_c18pos_metadata, diagnosis == "UC")
UC_HILneg_metadata <- subset(IBDMDB_HILneg_metadata, diagnosis == "UC")
UC_HILpos_metadata <- subset(IBDMDB_HILpos_metadata, diagnosis == "UC")
filename_c18neg <- UC_c18neg_metadata$Raw.file.name
filename_c18pos <- UC_c18pos_metadata$Raw.file.name
filename_HILneg <- UC_HILneg_metadata$Raw.file.name
filename_HILpos <- UC_HILpos_metadata$Raw.file.name
filename_c18neg <- sapply(filename_c18neg, function(x){paste0("/data/C18neg/rawfiles/",x)})
filename_c18pos <- sapply(filename_c18pos, function(x){paste0("/data/C18pos/rawfiles/",x)})
filename_HILneg <- sapply(filename_HILneg, function(x){paste0("/data/HILneg/rawfiles/",x)})
filename_HILpos <- sapply(filename_HILpos, function(x){paste0("/data/HILpos/rawfiles/",x)})

writeLines(filename_c18neg, "R_Works/outputs/filelist_c18neg.txt")
writeLines(filename_c18pos, "R_Works/outputs/filelist_c18pos.txt")
writeLines(filename_HILneg, "R_Works/outputs/filelist_HILneg.txt")
writeLines(filename_HILpos, "R_Works/outputs/filelist_HILpos.txt")
length(which(UC_c18neg_metadata$Azathioprine..Imuran..Azasan. == "Current"))
length(which(UC_c18neg_metadata$Mercaptopurine..Purinethol..6MP. == "Current"))


