str(result)
result$AZA.6MP
which(result$Control == 0)
which(result$AZA.6MP == 0)
result$Control[which(result$AZA.6MP == 0)]

### 6 mercaptopurine
if (T) {
  # Define the rt and mz range
  rtr <- c(0,840)
  mzr <- c(167.029669,167.049713)
  # Extract and check the chromatogram
  chr_raw <- chromatogram(rawdata, mz = mzr, rt = rtr)
  plot(chr_raw, col = group_colors[chr_raw$sample_group])
}
if (T) {
  # Define the rt and mz range
  rtr <- c(150,200)
  mzr <- c(167.029669,167.049713)
  # Extract and check the chromatogram
  chr_raw <- chromatogram(rawdata, mz = mzr, rt = rtr)
  plot(chr_raw, col = group_colors[chr_raw$sample_group])
}
if (T) {
  # Define the rt and mz range
  rtr <- c(480,560)
  mzr <- c(167.029669,167.049713)
  # Extract and check the chromatogram
  chr_raw <- chromatogram(rawdata, mz = mzr, rt = rtr)
  plot(chr_raw, col = group_colors[chr_raw$sample_group])
  title(main = "\t\t\t\t\tPlot", cex.main = 2, col.main = "coral")
}

# 找 6-MeMP
cwp <- CentWaveParam(peakwidth = c(3, 20), noise = 1000, ppm = 60)
xdata <- findChromPeaks(rawdata, param = cwp)
chr_xdata <- chromatogram(xdata, mz = mzr, rt = rtr)
plot(chr_xdata, col = group_colors[chr_xdata$sample_group])

subs_6memp <- filterRt(xdata, rt = rtr)
subs_6memp <- filterMz(subs_6memp, mz = mzr)
plotChromPeakImage(subs_6memp)
ints <- split(log2(chromPeaks(subs_6memp)[, "into"]),
              f = chromPeaks(subs_6memp)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[subs_6memp$sample_group],
        ylab = expression(log[2]~intensity), main = "Peak intensities",
        names = sub(basename(files), pattern = ".mzXML", replacement = "", fixed = T), las = 2)
legend("topleft", legend = c("Control","AZA/6MP") , 
       col = group_colors , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = F, inset = c(0.0, 0.00))
grid(nx = NA, ny = NULL)
pdp_6memp <- PeakDensityParam(sampleGroups = subs_6memp$sample_group,
                        minFraction = 0.4, bw = 20)
subs_6memp <- groupChromPeaks(subs_6memp, param = pdp_6memp)
head(featureValues(subs_6memp))
head(featureSummary(subs_6memp, group = subs_6memp$sample_group))
plot(chromatogram(subs_6memp, mz = mzr, rt = rtr), col = group_colors[subs_6memp$sample_group])
subs_feature_6memp <- featureSummary(subs_6memp, group = subs_6memp$sample_group)
colnames(subs_feature)
which(subs_feature[,11]==0)
featureNames(subs_6memp)
# 峰太少了，无法group成一个feature

# 找6-TGMP
if (T) {
  # Define the rt and mz range
  rtr <- c(0,840)
  mzr <- c(380.004427,380.080435)
  # Extract and check the chromatogram
  chr_raw <- chromatogram(rawdata, mz = mzr, rt = rtr)
  plot(chr_raw, col = group_colors[chr_raw$sample_group])
}

# 找6-TGDP
if (T) {
  # Define the rt and mz range
  rtr <- c(0,840)
  mzr <- c(459.962761,460.054763)
  # Extract and check the chromatogram
  chr_raw <- chromatogram(rawdata, mz = mzr, rt = rtr)
  plot(chr_raw, col = group_colors[chr_raw$sample_group])
}

# Find 6-TIMP
### 6 mercaptopuring
if (T) {
  # Define the rt and mz range
  rtr <- c(0,840)
  mzr <- c(365.009630,365.053434)
  # Extract and check the chromatogram
  chr_raw <- chromatogram(rawdata, mz = mzr, rt = rtr)
  plot(chr_raw, col = group_colors[chr_raw$sample_group])
}
