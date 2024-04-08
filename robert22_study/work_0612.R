
# Cohort 1
if(F){
dir_path <- "/nfs/data/metabolomics_chenrui/cohort1/massive.ucsd.edu/MSV000082094/updates/2019-01-10_rhmills_712042dc/raw/"
#Firstly Read all file in the filepath 
files <- list.files(dir_path, pattern = "\\.mzXML$", 
                    recursive = T, full.names = TRUE)## Change to .mzML if needed
# We define the /path/to/the/metadata.csv firstly. Remember must specify to .csv
metadata <- read.csv("Metadata/metadata_cohort2.csv",stringsAsFactors = FALSE, header = T)
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
# Set pd for phenodata(metadata) of xcms analysis 
length(which(sample_group == "AZA/6MP"))
length(which(metadata$current_IM%in%c("AZA", "6MP")))
pd <- data.frame(sample_name = sample_names,
                 sample_group = sample_group,
                 stringsAsFactors = F)
# Bind the list with path/to/each/file and respectively phenodata
rawdata <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd),mode = "onDisk")

}


# 用函数查看这个产物需要的是否合适
# examine # 6-TUA
plot_formula(rawdata = rawdata, formula = "C5H4N4O2S1", compound_name = "6TUA")
# 这里发现

cwp <- CentWaveParam(peakwidth = c(3, 30), noise = 1000, ppm = 60)
xdata <- findChromPeaks(rawdata, param = cwp)
rtr <- c(160,240)
mzr <- c(202.01,202.07)
chr_xdata <- chromatogram(xdata, mz = mzr, rt = rtr)
plot(chr_xdata, col = group_colors[chr_xdata$sample_group])
head(chromPeaks(xdata))
chromPeakData(xdata)

# Then we do a correlation heatmap of peak existens in all files 
plotChromPeakImage(xdata)
# Adjust runtime
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.4, bw = 10)
xdata <- groupChromPeaks(xdata, param = pdp)

# fill missing peaks xdata
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())
# plot again
plot(chr_xdata, col = group_colors[chr_xdata$sample_group])
# reset rtr to fetch all possible features
rtr <- c(0,840)
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
xdata_spectra
# Put xdata into a data.frame named results
result <- cbind(as.data.frame(featureDefinitions(xdata)),featureValues(xdata, value = "into"))

plot_volcano(result = result, controls = controls, treatments = treatments)
View(results)
metadata$Metagenomics_sequencing_name
