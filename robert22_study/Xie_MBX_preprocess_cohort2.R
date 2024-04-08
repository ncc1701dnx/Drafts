#!/usr/bin/env Rscript

# Load the needed packages
# mzR and xcms and edgeR are indispensable
library(mzR)
library(ggplot2)
library(RaMS)
library(readMzXmlData)#use to quickly read mzXML files if not use xcms
library(xcms) # LC-MS data pre-processing 
#library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
#library(pheatmap)
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


# We define the /path/to/the/metadata.csv firstly. Remember must specify to .csv
args <- commandArgs(trailingOnly = TRUE)
metadata_path <- args[1]

# If you don't use the Xie_MBX_ionazation_mode_separate.R firstly, please enable the following codes
# Active them by change if(F) to if(T)
if(F){
  
  positive_path <- args[2]
  negative_path <- args[3]
  
}

# or use this following codes if you used the Xie_MBX_ionazation_mode_separate.R
# Since we will group their findings togather
if(T){
  
  positive_files_path <- "dir_Ionization/positive_files.txt"
  negative_files_path <- "dir_Ionization/negative_files.txt"
  
  positive_files <- readLines(positive_files_path)
  negative_files <- readLines(negative_files_path)
  
}

# Read the metadata we generated use Xie_MBX_metadata_merge.R
# !!!MUST!!!DO!!! if you caompare the cohort 1 and cohort 2 todather

# or you can use cohort 1 and 2 separatly
#metadata <- read.csv("metadata/metadata_cohort2.csv",stringsAsFactors = FALSE, header = T)
metadata <- read.csv(args[1], stringsAsFactors = FALSE, header = T)

# For cohort 2, we need to specify those patients that have UC
metadata <- metadata[grep("UC",metadata_cohort2$Diagnosis),]

# Separate the treatment group and the control group and form correct metadata
#treatments <- metadata_UC[which(metadata_UC$Current_IM_Type_1AZA_26MP_3MTX%in%c("1","2")),]
#controls <- metadata_UC[which(metadata_UC$IM_Exposure=="0"),]
treatments <- metadata[which(metadata$Current_IM_Type_1AZA_26MP_3MTX%in%c("1","2")),]
controls <- metadata[which(metadata$IM_Exposure=="0"),]
treatments$current_IM <- rep("AZA/6MP", dim(treatments)[1])
controls$current_IM <- rep("control", dim(controls)[1])
treatments <- subset(treatments, !('Metabolomics_FileName' %in% "not applicable"))
controls <- subset(controls, !('Metabolomics_FileName' %in% "not applicable"))
metadata <- merge(treatments, controls, 
                  by = c('Metabolomics_FileName', 'current_IM'), all = T)

# Define a function firstly to process files with different ion polarity
# Define it firs because we will write function for cohort2

# Read the refference metadata as it's the xcms example said
data(hmdb_msdial)
pd
# Define a function to process files
?file.path
process_files <- function(files, pos_or_neg) {
  # Merge the files and metadata
  ff <- list.files(ifelse(pos_or_neg = "pos", positive_files, negative_files))
  pd <- data.frame(file = file.path("mzXML", ff), stringsAsFactors = F)
  ## In cohort 1 subset file name use strsplit("_") and sub(".mzXML")
  pd$lable <- basename(pd$file)
  pd$group <- merge(pd, metadata, by.x = "file_path", by.y = "Metabolomics_FileName")
  
  # Run pruned xcms set
  res <- runPrunedXcmsSet(
    files = pd$file_path, # path to the files
    pheno = pd, # phenotype data
    tmpdir = "Rds/temp", # folder to store temporary files, needs to be an empty folder or not existing then a new folder will be created
    postfilter = c(5, 3000), # at least 5 peaks has intensity higher than 3000
    noisefilter = c(ms1.noise = 1000, # ms1 spectrum, intensity lower than this will be considered as noise
                    ms1.maxPeaks = Inf,  # ms1 spectrum, maximum number of peaks to retain
                    ms1.maxIdenticalInt = 20, # noise peaks are usually same intensity, maximum number identical peaks to be considered as non-noise
                    ms2.noise = 50, # same but for ms2
                    ms2.maxPeaks = 100, 
                    ms2.maxIdenticalInt = 6,
                    BPPARAM=bpparam()), # parallel processing
    keepMS1 = T, # should the MS1 be kept, if yes, we can generate chromatogram give an intensity range. But it will make the data much larger. 
    mode = polarity, # ionization mode, "pos" or "neg"
    ref = hmdb_msdial, # feature annotation database()
    RTAdjustParam = PeakGroupsParam(minFraction = 0.5, span = 0.5), #  retention time adjustment parameter, will be passed to "adjustRtimePeakGroups"
    ppmtol = 20, # mass tolerance in PPM
    mclapplyParam = list(fun_parallel = mclapply, mc.cores = 1) # parallel, set mc.core to 4 when Linux, mc.core to 1 when Windows
  )
  
  return(res)
}

# Now ;et's process the positively and negativelu charged files
positive_data <- process_files(positive_files, "positive")
negative_data <- process_files(negative_files, "negative")

# Save the data to rds to run them in the diffrential separation
saveRDS(positive_data, "positive_data.rds")
saveRDS(negative_data, "negative_data.rds")


# End of cohort 2 differential analysis
# Perform statistical analysis then
# I will write another script for analysis, before that I haven't decide the write.out format