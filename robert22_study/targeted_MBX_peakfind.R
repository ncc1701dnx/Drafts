getwd()

## 先对于cohort 1的数据进行分析
dir_path <- "/nfs/data/metabolomics_chenrui/cohort1/massive.ucsd.edu/MSV000082094/updates/2019-01-10_rhmills_712042dc/raw/"
files <- list.files(dir_path, pattern = "\\.mzXML$", 
                    recursive = T, full.names = TRUE)## Change to .mzML if needed

## 再对于cohort 2的数据进行分析
dir_path <- "/nfs/data/cohort2_mbx_peaks/"
files <- list.files(dir_path, pattern = "\\.mzXML$", 
                    recursive = T, full.names = TRUE)## Change to .mzML if needed

# read metadata table of cohort 1
metadata_cohort1 <- read.csv("Metadata/metadata_cohort1.csv",stringsAsFactors = FALSE, header = T)
# read metadata table of cohort 2
metadata_cohort2 <- read.csv("Metadata/metadata_cohort2.csv",stringsAsFactors = FALSE, header = T)
###### Cohort 1 data preprocess
if (T) {  #Cohort 1
  # Separate the treatment group and the control group and form correct metadata
  treatments_cohort1 <- metadata_cohort1[which(metadata_cohort1$current_IM%in%c("AZA","6MP")),]
  controls_cohort1 <- metadata_cohort1[which(metadata_cohort1$IM_type=="0"),]
  treatments_cohort1$current_IM <- rep("AZA/6MP", dim(treatments_cohort1)[1])
  controls_cohort1$current_IM <- rep("control", dim(controls_cohort1)[1])
  metadata_cohort1 <- merge(treatments_cohort1, controls_cohort1, 
                    by = c('vial_name', 'current_IM'), all = T)
  
  # Now we firstly read those .mzXML files into the the xcms readable file
  # In terms of Cohort 1, we have to subset the name in the filelist to be identical to metadata
  # Run a for loop to read all files exist in metadata and mark their respectively group
  sample_group <- c()
  sample_names <- c()
  ff <- c()
  for (file in files) {
    if (strsplit(basename(file), "_")[[1]][1] %in% metadata_cohort1$vial_name){
      ifelse(metadata_cohort1$current_IM[which(metadata_cohort1$vial_name==strsplit(basename(file), "_")[[1]][1])]=="AZA/6MP",
             sample_group <- c(sample_group, "AZA/6MP"), sample_group <- c(sample_group,"Control"))
      sample_names <- c(sample_names, sub(basename(file), pattern = ".mzXML", replacement = "", fixed = T))
      ff <- c(ff, file)
    } else {
      print(paste0(basename(file),"Do not exist in metadata"))
    }
  }
  files <- ff
}

###### Cohort 2
if (T) {  #cohort 2
  # For cohort 2, we need to specify those patients that have UC
  metadata_cohort2 <- metadata_cohort2[grep("UC",metadata_cohort2$Diagnosis),]
  
  # Separate the treatment group and the control group and form correct metadata
  #treatments <- metadata_UC[which(metadata_UC$Current_IM_Type_1AZA_26MP_3MTX%in%c("1","2")),]
  #controls <- metadata_UC[which(metadata_UC$IM_Exposure=="0"),]
  treatments_cohort2 <- metadata_cohort2[which(metadata_cohort2$Current_IM_Type_1AZA_26MP_3MTX%in%c("1","2")),]
  controls_cohort2 <- metadata_cohort2[which(metadata_cohort2$IM_Exposure=="0"),]
  treatments_cohort2$current_IM <- rep("AZA/6MP", dim(treatments_cohort2)[1])
  controls_cohort2$current_IM <- rep("control", dim(controls_cohort2)[1])
  treatments_cohort2 <- subset(treatments_cohort2, !('Metabolomics_FileName' %in% "not applicable"))
  controls_cohort2 <- subset(controls_cohort2, !('Metabolomics_FileName' %in% "not applicable"))
  metadata_cohort2 <- merge(treatments_cohort2, controls_cohort2, 
                    by = c('Metabolomics_FileName', 'current_IM'), all = T)
  
  # Now we firstly read those .mzXML files into the the xcms readable file
  #files <- list.files(path = "experiment_data/cohort2_metabolomics_check/", recursive = T, full.names = T)
  # Run a for loop to read all files exist in metadata and mark their respectively group
  sample_group <- c()
  sample_names <- c()
  ff <- c()
  for (file in files) {
    if (basename(file) %in% metadata_cohort2$Metabolomics_FileName){
      ifelse(metadata_cohort2$current_IM[which(metadata_cohort2$Metabolomics_FileName == basename(file))]=="AZA/6MP",
             sample_group <- c(sample_group, "AZA/6MP"), sample_group <- c(sample_group,"Control"))
      sample_names <- c(sample_names, sub(basename(file), pattern = ".mzXML", replacement = "", fixed = T))
      ff <- c(ff, file)
    } else {
      print(paste0(basename(file),"Do not exist in metadata"))
    }
  }
  files <- ff
}

# 检查一下处理好的数据，并用pd把这两个处理的数据bind起来
length(which(sample_group == "AZA/6MP"))
length(which(sample_group == "Control"))
pd <- data.frame(sample_name = sample_names,
                 sample_group = sample_group,
                 stringsAsFactors = F)
# Bind the list with path/to/each/file and respectively phenodata
rawdata <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd),mode = "onDisk")

# 用写好的函数（MBX_functions.R里面的）来处理每一个化合物
# 6MeMP
plot_formula(rawdata = rawdata,formula = "C6H6N4S1", compound_name = "6-MeMP")
# 6TUA
plot_formula(rawdata = rawdata, formula = "C5H4N4O2S1", compound_name = "6TUA")
# 6-MeTIMP
plot_formula(rawdata = rawdata, formula = "C11H15N4O7P1S1", compound_name = "6-MeTIMP")
# 6-TIMP
plot_formula(rawdata = rawdata, formula = "C10H13N4O7P1S1", compound_name = "6-TIMP")
# AZA
plot_formula(rawdata = rawdata, formula = "C9H7N7O2S1", compound_name = "AZA")
# 6-MP
plot_formula(rawdata = rawdata, formula = "C5H4N4S1", compound_name = "6-MP")
# 6-TGMP
plot_formula(rawdata = rawdata, formula = "C10H14N5O7P1S1", compound_name = "6-TGMP")
# 6-TGDP
plot_formula(rawdata = rawdata, formula = "C10H15N5O10P2S1", compound_name = "6-TGDP")
# 6-TGTP
plot_formula(rawdata = rawdata, formula = "C10H16N5O13P3S1", compound_name = "6-TGTP")

# 然后根据图片设定rawdata的chrome peak finding 的参数
# SET ppm = 60 for cohort 1 but ppm = 100 for cohort2
cwp <- CentWaveParam(peakwidth = c(5,30), noise = 1000, ppm = 60)
# for most compounds set lower limit as 3
xdata <- findChromPeaks(rawdata, param = cwp)
# 对于xdata进行后续处理
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.4, bw = 10)
xdata <- groupChromPeaks(xdata, param = pdp)
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())
# 对于处理好的数据，可以跑一下图片看看我们选出来的峰现在有没有阴影（Chrome peak被找到）
# for 6-MeMP M+H
rtr <- c(0,840)
mzr <- c(380.0050,380.0809)
chr_xdata <- chromatogram(xdata, mz = mzr, rt = rtr)
 # group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
 # names(group_colors) <- c("Control", "AZA/6MP")
plot(chr_xdata, col = group_colors[chr_xdata$sample_group])
# 有阴影就可以通过filterMz函数来找是否有需要的peak或者feature
# subs_MeMP <- filterRt(xdata, rt = c(170, 190)) # 因为Rt不够准确，不建议使用
subs_MeMP <- filterMz(xdata, mz = mzr)
# 检查subset出来的数据
table(msLevel(subs_MeMP))
table(msLevel(xdata))
# Set precusor m/z value as well as the intensity value for each ms 2 spectra
# 这一步将峰信息导入到MS2 level spectra里面，如果有MS 2 色谱的话可以直接提取信息
subs_MeMP %>%
  filterMsLevel(2L) %>%
  precursorMz() %>%
  head()
subs_MeMP %>%
  filterMsLevel(2L) %>%
  precursorIntensity() %>%
  head()

# Check ms 2 level spectra
MeMP_spectra <- chromPeakSpectra(
 subs_MeMP, msLevel = 2L, return.type = "Spectra")
# 如果没有MS 2 level spectra或者太少，会报warning
# 查看MS 2 spectra信息，这里发现虽然有MS 2 spectra，但是Rt太早了，不是我们想要的
MeMP_spectra
# 查看能否找到一个feature（也就是被xcms计算为一个在很多sample中都出现的峰）
featureSpectra(subs_MeMP)
featureDefinitions(subs_MeMP)
# 如果feature里面没有，在chrome peak里面找，这些都是出现比较少的峰
head(chromPeaks(subs_MeMP))
chromPeaks(subs_MeMP)
# 在这里我们找到了这个峰，注意对比M/Z和Rtime
# 没有MS 2 spectra

# 接下来我们需要找到这个chrome peak （CP44411）是哪一个sample的peak
# 然后做成表方便输出
which(featureDefinitions(subs_MeMP)$mzmed )
MBX_targets <- data.frame("compound" = "6-MeMP","formula" =  "C6H6N4S1", "chrome peak" = "CP23072", "feature" = NA, "mz" = 184.0517718, "rt" = 235.38201904,
                          "MS2" = NA, "adducts" = "NH4" ,"files" = "FIT270E3CAL_BC5_01_52126.mzXML")
  MBX_targets[2,] <- c("6TUA", "C5H4N4O2S1", NA, NA, NA, NA, NA, NA,NA)
MBX_targets[, "adducts"] <- c("H", "NH4")
  MBX_targets[3,] <- c("6-MeTIMP", "C11H15N4O7P1S1", NA, NA, NA, NA, NA, NA,NA)
MBX_targets[4,] <- c("6-TIMP", "C10H13N4O7P1S1", NA, NA, NA, NA, NA, NA,NA)
MBX_targets[5,] <- c("AZA", "C9H7N7O2S1", NA,NA,NA,NA,NA,NA,NA)
MBX_targets[6,] <- c("6-MP", "C5H4N4S1", NA,NA,NA,NA,NA,NA,NA)
MBX_targets[7,] <- c("6-TGMP", "C10H14N5O7P1S1", NA,NA,NA,NA,NA,NA,NA)
MBX_targets[8,] <- c("6-TGDP", "C10H15N5O10P2S1",NA,NA,NA,NA,NA,NA,NA)
MBX_targets[9,] <- c("6-TGTP", "C10H16N5O13P3S1",NA,NA,NA,NA,NA,NA,NA)

# 查询某个feature在哪些Sample中出现
featureSummary(subs_MeMP, perSampleCounts = T)

peaksCount(subs_MeMP)
chromPeakData(subs_MeMP)
fileNames(subs_MeMP)
chromPeaks(subs_MeMP)

# 写一个函数自动给出输入的CP的来源文件名
# 需要这个函数的原因是没有grep cp的来源的函数
# 如果有Feature，用featureSummary(subs_MeMP, perSampleCounts = T)
# 查询某个feature在哪些Sample中出现
# featureSummary(subs_MeMP, perSampleCounts = T)
cp_file_name <- function(xdata, cp_name){
  sample_names <- c()
  for (i in 1:length(cp_name)){
    # Get the peak matrix
    peak_data <- chromPeaks(xdata)
    # Find the row number of the peak in the peak matrix
    # Important since we do not use the subset data but rather the whole xdata as input
    peak_row <- which(rownames(peak_data) == cp_name[i])
    # Access the peak row
    sample_index <- peak_data[peak_row, "sample"]
    # Find the directory of upload file
    sample_dir <- fileNames(subs_MeMP)[sample_index]
    # Subset to the file name. Should be the last after strsplit
    sample_name <- strsplit(sample_dir, "/")[[1]][length(strsplit(sample_dir, "/")[[1]])]
    sample_names <- c(sample_names, sample_name)
  }
  sample_names
}
cp_file_name(xdata, 'CP23072')
tt2 <- cp_file_name(xdata, c("CP27056", "CP27180", "CP45817", "CP45840"))
MBX_targets[,"files"] <- c(tt1, NA, "FIT100E3CAL_GC9_01_12077.mzXML FIT100E3CAL_GC9_01_12077.mzXML", "FIT154E3CAL_GA4_01_12046.mzXML FIT154E3CAL_GA4_01_12046.mzXML",NA,NA)

MBX_targets_cohort1 <- MBX_targets
MBX_targets_cohort2 <- MBX_targets
write.csv(MBX_targets_cohort1, "outputs/targeted_MBX_cohort1_results.csv", quote = F)
write.csv(MBX_targets_cohort2, "outputs/targeted_MBX_cohort2_results.csv", quote = F)
