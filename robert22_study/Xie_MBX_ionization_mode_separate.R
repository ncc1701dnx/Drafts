#!/usr/bin/env Rscript

# Use it as Rscript Xie_MBX_ionazation_mode_separate.R "/your/real/path/to/.mzXML/files/"
# if use .mzML file, change the parameter in the loop

# Load the needed packages
# mzR and xcms are indispensable
library(mzR)
library(ggplot2)
library(RaMS)
#library(readMzXmlData)#use to quickly read mzXML files if not use xcms
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
library(MAIT) #dependency for xcmsViewer

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Argument is the !real-path! to the data file
dir_path <- args[1]
#dir_path <- "experiment_data/cohort2_metabolomics_check/"

#Firstly Read all file in the filepath 
files <- list.files(dir_path, pattern = "\\.mzXML$", full.names = TRUE)## Change to .mzML if needed

# Now we need firstly check the ionazation mode
# Despite almost all of them should be positive (·•᷄ࡇ•)
# Initialize lists to store files by ionization mode
positive_files <- list()
negative_files <- list()

# Loop over each file
  for (file in files) {
    # Open those .mzXML file use xcms package
    mzxml <- openMSfile(file)
    
    # Get header information for the first spectrum to get the ionization mode 
    # The ionization mode are for all metabolites in one run the same
    header <- header(mzxml, 1)
    
    # Check the polarity and add the file to the given list
    if (header$polarity == 1) {
      positive_files <- c(positive_files, file)
    } else if (header$polarity == -1) {
      negative_files <- c(negative_files, file)
    }
    
    # Close the mzXML file
    # ! Highly necessary, or the system will crash if you have opened too many spectra!!!!!!
    close(mzxml)
  }
  
  # Create the dir_Ionization directory if it doesn't exist
  # If you wanna run for cohort 1 and cohort 2 in one run please tell me
  if (!dir.exists("dir_Ionization")) {
    dir.create("dir_Ionization")
  }
  
  # Write the lists of files as tables in the aformentioned directory
  write.table(positive_files, file = "dir_Ionization/positive_files.txt", row.names = FALSE)
  write.table(negative_files, file = "dir_Ionization/negative_files.txt", row.names = FALSE)
  
  #End of ionization mode determination