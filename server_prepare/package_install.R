### Package install script. If need to install any package, type here.
###
### Packages for MBX analysis with XCMS
getwd()
setwd("/home/rstudio/R_Works/")
install.packages("ggplot2")
install.packages("stringr")
install.packages("XML")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("Rcpp")
install.packages("KernSmooth")
library(KernSmooth)
BiocManager::install("affyio")
BiocManager::install("affy")
install.packages("RCurl")
install.packages("ncdf4")
BiocManager::install("GenomeInfoDb")
BiocManager::install("Spectra")
BiocManager::install("mzR", force = T)
install.packages("ragg")
install.packages("RaMS", dependencies = T)
BiocManager::install("xcms")
install.packages("RColorBrewer")
install.packages("pander")
install.packages("magrittr")
install.packages("pheatmap")
BiocManager::install("SummarizedExperiment",  force = T)
BiocManager::install("MsFeatures", force = T)
install.packages("rlang")
install.packages("dplyr")
BiocManager::install("MAIT", force = T)
BiocManager::install("MSnbase", force = T) 
install.packages("rcartocolor")

library.path <- .libPaths()
library("xcms", lib.loc = library.path)[1]
.libPaths()[1]
.libPaths(new = "/home/rstudio/R/x86_64-pc-linux-gnu-library/4.3/")

BiocManager::install("MsBackendMgf", force = T)### For write XCMS Spectra object out as .mgf files
library('MsBackendMgf')