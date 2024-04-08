#!/usr/bin/env Rscript

#######
##########
################################ Package loading step ##############################################################################################
# Load the needed packages
# mzR and xcms and edgeR are indispensable
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
library('MsBackendMgf') # Write Spectra object out as .mgf
################################# End of package loading ############################################################################################
#############
#######

#######
###############
############################### Read file and reorganize file based on metadata #####################################################################
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Argument is the !real-path! to the data file
dir_path <- args[1]

## for cohort1 ### I have moved the files to new folder to prevent the permission problem ##################################
dir_path <- "/nfs/data/metabolomics_chenrui2/cohort1/massive.ucsd.edu/MSV000082094/updates/2019-01-10_rhmills_712042dc/raw/"

#Firstly Read all file in the filepath 
files <- list.files(dir_path, pattern = "\\.mzXML$", 
                    recursive = T, full.names = TRUE)## Change to .mzML if needed

## for cohort2 ### I have already put all files regardless of their plates in one folder ###################################
dir_path <- "/nfs/data/metabolomics_chenrui2/cohort2/mzXML_files/"

#Firstly Read all file in the filepath 
files <- list.files(dir_path, pattern = "\\.mzXML$", 
                    recursive = T, full.names = TRUE)## Change to .mzML if needed

# Now we need firstly check the ionization mode
# Despite almost all of them are positive (·•᷄ࡇ•)
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
# Write it out in case we need.
write.table(positive_files, file = "dir_Ionization/positive_files.txt", row.names = FALSE)
write.table(negative_files, file = "dir_Ionization/negative_files.txt", row.names = FALSE)

# End of ionization mode determination

# Start of data preprocessing

# Start of peak detection and de-noising

# Read the 
# We define the /path/to/the/metadata.csv firstly. Remember must specify to .csv
# cohort1
metadata <- read.csv("Metadata/metadata_cohort1.csv",stringsAsFactors = FALSE, header = T)
table(metadata$current_IM)
# cohort2
metadata <- read.csv("Metadata/metadata_cohort2.csv",stringsAsFactors = FALSE, header = T)
table(metadata$Current_IM_Type_1AZA_26MP_3MTX)
table(metadata$Diagnosis)
#metadata <- read.csv(args[2], stringsAsFactors = FALSE, header = T)

###### Cohort 2
if (T) {  #cohort 2
  # For cohort 2, we need to specify those patients that have UC
  #metadata <- metadata[grep("UC",metadata$Diagnosis),] ## 如果你只需要研究UC病人
  metadata <- metadata[metadata$Diagnosis %in% c("UC", "CD"), ]#如果需要研究所有IBD
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
  
  # Now we firstly read those .mzXML files into the the xcms readable file
  #files <- list.files(path = "experiment_data/cohort2_metabolomics_check/", recursive = T, full.names = T)
  # Run a for loop to read all files exist in metadata and mark their respectively group
  sample_group <- c()
  sample_names <- c()
  ff <- c()
  for (file in files) {
    if (basename(file) %in% metadata$Metabolomics_FileName){
      ifelse(metadata$current_IM[which(metadata$Metabolomics_FileName == basename(file))]=="AZA/6MP",
             sample_group <- c(sample_group, "AZA/6MP"), sample_group <- c(sample_group,"Control"))
      sample_names <- c(sample_names, sub(basename(file), pattern = ".mzXML", replacement = "", fixed = T))
      ff <- c(ff, file)
    } else {
      print(paste0(basename(file),"Do not exist in metadata"))
    }
  }
  files <- ff
}


###### Cohort 1
if (T) {  #Cohort 1
  # Separate the treatment group and the control group and form correct metadata
  treatments <- metadata[which(metadata$current_IM%in%c("AZA","6MP")),]
  controls <- metadata[which(metadata$IM_type=="0"),]
  treatments$current_IM <- rep("AZA/6MP", dim(treatments)[1])
  controls$current_IM <- rep("control", dim(controls)[1])
  metadata <- merge(treatments, controls, 
                    by = c('vial_name', 'current_IM'), all = T)
  
  # Now we firstly read those .mzXML files into the the xcms readable file
  # In terms of Cohort 1, we have to subset the name in the filelist to be identical to metadata
  # Run a for loop to read all files exist in metadata and mark their respectively group
  sample_group <- c()
  sample_names <- c()
  ff <- c()
  for (file in files) {
    if (strsplit(basename(file), "_")[[1]][1] %in% metadata$vial_name){
      ifelse(metadata$current_IM[which(metadata$vial_name==strsplit(basename(file), "_")[[1]][1])]=="AZA/6MP",
             sample_group <- c(sample_group, "AZA/6MP"), sample_group <- c(sample_group,"Control"))
      sample_names <- c(sample_names, sub(basename(file), pattern = ".mzXML", replacement = "", fixed = T))
      ff <- c(ff, file)
    } else {
      print(paste0(basename(file),"Do not exist in metadata"))
    }
  }
  files <- ff
}

############################ End of file reorganize, we now have grouped the files #################################################################
##############
#######

#######
##############
########################### Parameter try and settings lead to QC ################################################################################## 
# Set pd for phenodata(metadata) of xcms analysis 
length(which(sample_group == "AZA/6MP"))
dim(controls)
dim(treatments)
pd <- data.frame(sample_name = sample_names,
                 sample_group = sample_group,
                 stringsAsFactors = F)
# Bind the list with path/to/each/file and respectively phenodata
rawdata <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd),mode = "onDisk")
head(rtime(rawdata))

# Plot the raw XIC figure
group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("Control", "AZA/6MP")
str(rawdata)
bpis <- chromatogram(rawdata, aggregationFun = "max")
plot(bpis, col = group_colors[rawdata$sample_group])
# Check the difference between two groups 
tic(rawdata)
tc <- split(tic(rawdata), f = fromFile(rawdata))
boxplot(tc, col = group_colors[rawdata$sample_group],
        ylab = "intensity", main = "Total ion current",ylim = c(0, 6e+3), 
        names = sub(basename(files), pattern = ".mzXML", replacement = "", fixed = T),
        las = 2)#,ylim = c(5000, 1e+5))
legend("topleft", legend = c("Control","AZA/6MP") , 
       col = group_colors , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = F, inset = c(0.0, 0.00))
## Bin the BPC
bpis_bin <- MSnbase::bin(bpis, binSize = 2)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- rawdata$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = rawdata$sample_group)
rownames(ann) <- rawdata$sample_name

## Perform the cluster analysis
dim(cormat)
is.na(cormat) %>% table()
na.omit(cormat, invert = F)
pheatmap(cormat, annotation = ann,
         annotation_color = list(group = group_colors), na_col = "grey90")


# Now we filter noises and fetch all detected non-noise peaks
# First try which parameters should be use
# We take the standard hexakis (1H,1H,3H-tetrafluoropropoxy) phosphazene ions (Synquest Laboratories, m/z 922.0098)
# ppm set as 60 as the tolerance for MS1 is 0.03 Da
# Set to F to ignore all code inside this if loop
if (T) {
  # Define the rt and mz range
  rtr <- c(0,840)
  mzr <- c(167.029669,167.049713)
  # Extract and check the chromatogram
  chr_raw <- chromatogram(rawdata, mz = mzr, rt = rtr)
  plot(chr_raw, col = group_colors[chr_raw$sample_group])
}
# It seems the peakwidth could be 4,7,10,12,, we set 3 to 20
# Noise for ms 1 is set as the orginal paper suggested 
# Prefilter is set as the whole range of XIC, can set more accurate if more information available
cwp <- CentWaveParam(peakwidth = c(0.275, 2.75), noise = 1000, ppm = 10)
#cwp <- CentWaveParam(peakwidth = c(0.275, 5), noise = 1000, ppm = 10)
?CentWaveParam
xdata <- findChromPeaks(rawdata, param = cwp)
chr_xdata <- chromatogram(xdata, mz = mzr, rt = rtr)
plot(chr_xdata, col = group_colors[chr_xdata$sample_group])
head(chromPeaks(xdata))
chromPeakData(xdata)
# We can see that their are many NA peaks in those data but we don't fill those gap here

# Now we provide the summary data to our file, to check if we have good value data
summary_fun <- function(x)
  c(peak_count = nrow(x), rt = quantile(x[, "rtmax"] - x[, "rtmin"]))
TB <- lapply(split.data.frame(
  chromPeaks(xdata), f = chromPeaks(xdata)[, "sample"]),
  FUN = summary_fun)
TB <- do.call(rbind, TB)
rownames(TB) <- basename(fileNames(xdata))

pandoc.table(
  TB,
  caption = paste0("Summary statistics on identified chromatographic",
                   " peaks. Shown are number of identified peaks per",
                   " sample and widths/duration of chromatographic ",
                   "peaks."))
# Seems it's better if we can use the fill and refineChromPeaks function,
# but we do not do it now since it can interfer to our results

# We then provide some plots that can be more straightforward
# Fistly provide the chromatographic peaks graph for all files
# If too many file, then shouldn't plot for each of them, set default to F
if (F){
  for (plotNUM in 1:(length(files))) {
    plotChromPeaks(xdata, file = plotNUM)
  }
}
# For almost all of them, there are no many peaks in the middle of the chromatograph

# Then we do a correlation heatmap of peak existens in all files 
plotChromPeakImage(xdata)
# Yes most of those peaks are in the last minute of run
# Identical to the XIC graph we plot at first time

# Extract a list of per-sample peak intensities (in log2 scale)
# Check the difference between groups again
ints <- split(log2(chromPeaks(xdata)[, "into"]),
              f = chromPeaks(xdata)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
        ylab = expression(log[2]~intensity), main = "Peak intensities",
        names = sub(basename(files), pattern = ".mzXML", replacement = "", fixed = T), las = 2)
legend("topleft", legend = c("Control","AZA/6MP") , 
       col = group_colors , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = F, inset = c(0.0, 0.00))
grid(nx = NA, ny = NULL)
# We can see after the peak detection and de-noise, the outliers before are now the main peaks(note the log2)

# End of Peak detection and de-noise

# Start of peak alignment within group

# Accoruding to the data we have, we have average 60-70 peaks in 800 seconds,so a bin of 2-3 is fairly ok
# However, we have 20 peaks in last one minute. So may need more accurate. We set binSize =1 as a start point
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.30)) ## based on the paper
head(adjustedRtime(xdata))
# before adjust:
#F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006 
#0.533    0.662    0.769    0.874    0.977    1.079 
# after adjust BinSize 1.0:
#F1.S0001  F1.S0002  F1.S0003  F1.S0004  F1.S0005  F1.S0006 
#0.5290000 0.6583341 0.7656112 0.8708831 0.9741499 1.0764141 
# after adjust BinSize 0.4:
#F1.S0001  F1.S0002  F1.S0003  F1.S0004  F1.S0005  F1.S0006 
#0.5290000 0.6583341 0.7656112 0.8708831 0.9741499 1.0764141 

## Get the base peak chromatograms after adjustment
if (F){
  bpis_adj <- chromatogram(xdata, aggregationFun = "max", include = "none", mz = mzr, rt = rtr)
  plot(bpis_adj, col = group_colors[bpis_adj$sample_group])
  ## Plot also the difference of adjusted to raw retention time.
  plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group])
  # And cmpare the XIC before and after the alignment
  plot(bpis_adj, col = group_colors[bpis_adj$sample_group])
  plot(bpis, col = group_colors[rawdata$sample_group])
  # NO better than before, more accuracy needed(binSize to 0.4?)
  # NO big difference, set = 1.0 is ok
  # Set plot frame back to one figure a time
  par(mfrow = c(1,1))
}
if (F){
  ## Try plot the chromatograme in set 
  mzr <- c(330, 360)
  chr_mzr <- chromatogram(xdata,mz = mzr)
  ## Define the parameters for the peak density method
  pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                          minFraction = 0.4, bw = 10)
  sample_colors <- group_colors[chr_mzr$sample_group]
  plotChromPeakDensity(chr_mzr, col = sample_colors, param = pdp,
                       peakBg = sample_colors[chromPeaks(chr_mzr)[, "sample"]],
                       peakCol = sample_colors[chromPeaks(chr_mzr)[, "sample"]],
                       peakPch = 16)
}
# bw = 10 is okay

# Let's start the performing of correspondence
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.4, bw = 10)
xdata <- groupChromPeaks(xdata, param = pdp)
res <- quantify(xdata, value = "into")
#colData(res)
# Feature annotations inspection with rowData and/or featureDif
rowData(res)
featureDefinitions(xdata)
# Then this is the table we need to identifie which compound are not enlisted in control but in AZA/6MP
# Try to write out those file
# We can fill the peaks in different group respectively since we have (6) spectra with very little peaks
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())
head(featureValues(xdata))
featureDefinitions(xdata)
## Missing values before filling in peaks
apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))
## Missing values after filling in peaks
apply(featureValues(xdata), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))
# Write filled spectra into res
assays(res)$raw_filled <- featureValues(xdata, filled = TRUE)
# End of peak alignment within group

# Start of data analysis

# Feature abundance inspection with assay
assayNames(res)
head(assay(res))
# General per-feature summary with featureSummary
head(featureSummary(xdata, group = xdata$sample_group))#rsd for relative standard deviation
# This can also been write out 

# PCA analysis
## Extract the features and log2 transform them
ft_ints <- log2(assay(res))
## Perform the PCA omitting all features with an NA in any of the
## samples. Also, the intensities are mean centered.
pc <- prcomp(t(na.omit(ft_ints)), center = TRUE)
## Plot the PCA
print(pc$rotation)
cols <- group_colors[xdata$sample_group]
pcSummary <- summary(pc)
par(mfrow = c(1,1))
plot(pc$x[, 1], pc$x[,2], pch = 21, main = "",
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = "darkgrey", bg = cols, cex = 2)
grid()

text(pc$x[, 1], pc$x[,2], labels = xdata$sample_name, col = "darkgrey",
     pos = 3, cex = 1)

# Check the chromatography here
floor(runif(4, min=0, max=800))
feature_chroms <- featureChromatograms(xdata, features = floor(runif(4, min=0, max=800)))
feature_chroms
# And plot the extracted ion chromatograms.
sample_colors <- group_colors[xdata$sample_group]
plot(feature_chroms, col = sample_colors,
     peakBg = sample_colors[chromPeaks(feature_chroms)[, "sample"]])
######################################### End of QC step #######################################################################################
############
#######


#######
#############
######################################## Grep of MS2 Level information #########################################################################
# Check for ms 1 data and ms 2 data existency
table(msLevel(xdata))
# Set precusor m/z value as well as the intensity value for each ms 2 spectra
xdata %>%
  filterMsLevel(2L) %>%
  precursorMz() %>%
  head()
xdata %>%
  filterMsLevel(2L) %>%
  precursorIntensity() %>%
  head()

# Check ms 2 level spectra
xdata_spectra <- chromPeakSpectra(
  xdata, msLevel = 2L, return.type = "Spectra")
table(xdata_spectra)

# Get MS1 level precursor feature name
xdata_spectra$peak_id
length(unique(xdata_spectra$peak_id))
# show 11961 peak ids in cohort1
xdata_spectra$
## check MS2 data qualities
length(xdata_spectra)
head(xdata_spectra)
class(xdata_spectra)
plotSpectra(xdata_spectra[100])
## 写出去，用sirius进行分析
## 
#writeMgfData(xdata_spectra, file = "/nfs/data/metabolomics_chenrui2/cohort1/MS2_mgfs/Cohort1_MS2.mgf") ## 看来spectra对象不支持直接写出.mgf
##
## 尝试将Spectra转化为MSNbase吧
#xdata_spectra_MSNbase <- as(xdata_spectra, "MSnExp") ## no working
#showMethods("writeMgfData")# 检查writeMgf支持什么
#class(xdata_spectra)
## 试点儿别的
#xdata_spectra_MSNbase <- as(xdata_spectra, "Spectrum") ## no working too
export(xdata_spectra, MsBackendMgf(), file = "/nfs/data/metabolomics_chenrui2/cohort1/MS2_mgfs/Cohort1_MS2.mgf") ## This can really do!
### 不过峰太多了，我们需要尝试仅将需要的化合物的峰提取出来
### 具体请查找work_0811_grep_RobertsMS2.R



# Grab one feature at rtime 787.835 to analysis
str(xdata_spectra)
chromPeaks(xdata, mz = 644.0109, ppm = 0) 
ex_id <- rownames(chromPeaks(xdata, mz = 644.0109, ppm = 0))
ex_spectra <- xdata_spectra[xdata_spectra$peak_id == 'CP488']
ex_spectrum <- combineSpectra(ex_spectra, FUN = combinePeaks, ppm = 0,
                              peaks = "intersect", minProp = 0.8,
                              intensityFun = median, mzFun = median,
                              f = ex_spectra$peak_id)
# 这不是正确方法！
if(F){
ms2_specs <- featureSpectra(
  xdata,
  msLevel = 2,
  expandRt = 5,
  expandMz = 2,
  ppm = 10,
  skipFilled = FALSE,
  return.type = c("MSpectra", "list")
)
}
##################################### End of MS2 Level Information extract #############################################################################
############
######

# Check process history
processHistory(xdata)

#提取某一步的数据
ph <- processHistory(xdata, type = "Retention time correction")

#提取数据的参数
processParam(ph[[1]])

#提取某个文件的数据
subs <- filterFile(xdata, file = c(2, 4))

#提取数据并留取保留时间
subs <- filterFile(xdata, keepAdjustedRtime = TRUE)

#按保留时间提取数据
subs <- filterRt(xdata, rt = c(3000, 3500))
range(rtime(subs))

#提取某个文件的所有数据
one_file <- filterFile(xdata, file = 3)
one_file_2 <- xdata[fromFile(xdata) == 3]

# Check the peak file 
head(chromPeaks(xdata))

#导出数据
result <- cbind(as.data.frame(featureDefinitions(xdata)),featureValues(xdata, value = "into"))
View(result)
write.table(result,file="outputs/cohrot1_result.csv",sep="\t",quote=F)

#PCA
values <- groupval(xdata)

data <- t(values)

pca.result <- pca(data)

library(pcaMethods)# Should be a aprt from MAIT

pca.result <- pca(data)

loadings <- pca.result@loadings

scores <- pca.result@scores

plotPcs(pca.result, type = "scores",col=as.integer(sampclass(xdata)) + 1)

#MDS
library(MASS)

for (r in 1:ncol(data)){
  
  data[,r] <- data[,r] / max(data[,r])

data.dist <- dist(data)
}

mds <- isoMDS(data.dist)

plot(mds$points, type = "n")

text(mds$points, labels = rownames(data),col=as.integer(sampclass(xdata))+ 1 )


# Save all variables and environments
save.image(file = "RDs/cohort2_try0601.RData")
