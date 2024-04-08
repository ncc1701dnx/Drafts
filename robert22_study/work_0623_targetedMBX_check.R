str(xdata)
featureDefinitions(xdata)
head(adjustedRtime(xdata))
xdata %>%
  filterMsLevel(2L) %>%
  precursorMz() %>%
  head()
xdata %>%
  filterMsLevel(2L) %>%
  precursorIntensity() %>%
  head()
xdata_spectra <- chromPeakSpectra(
  xdata, msLevel = 2L, return.type = "Spectra")
xdata_spectra$peak_id
featureDefinitions(xdata_spectra)
chromPeaks(xdata_spectra, mz = c(397.0301,397.1093), ppm = 60)
table(mslevel(xdata_spectra))
metadata_cohort1$vial_name
strsplit(MBX_targets_cohort1$files, "_")[[1]][1]
metadata_cohort1$mayo_endoscopic_score[which(metadata_cohort1$vial_name == strsplit(MBX_targets_cohort1$files, "_")[[1]][1])]
table(metadata_cohort1$historic_extent_category)
getwd()
tt <- sapply(strsplit(as.character(MBX_targets_cohort1$files), "_"), `[`, 1)
sapply(tt, function(x){
  metadata_cohort1$mayo_endoscopic_score[which(metadata_cohort1$vial_name == x)]
})
getwd()
table(metadata_cohort2$Current_IM)
which(metadata_cohort2$Current_IM == 1)
metadata_cohort2$Sample_name[which(metadata_cohort2$Current_IM == 1)]
metadata_cohort2$UCEIS[which(metadata_cohort2$Sample_name == "FIT270E3CAL")]
metadata_cohort2$UCEIS[which(metadata_cohort2$Current_IM == 1)]
metadata_cohort2$UCEIS[which(metadata_cohort2$Current_IM == 0)]
metadata_cohort1$UCEIS_endoscopic_score[which(metadata_cohort1$current_IM == "not applicable")]
metadata_cohort1$UCEIS_endoscopic_score[which(metadata_cohort1$current_IM%in%c("AZA","6MP"))]
metadata_cohort1$vial_name[which(metadata_cohort1$current_IM%in%c("AZA","6MP"))]
which(metadata_cohort2$Sample_name %in% c("FIT002CCAL", "FIT288E3CAL", "FIT242E3CAL", "FIT233E3CAL", "FIT211E3CAL") )
metadata_cohort2$[which(metadata_cohort2$Current_IM == 1)]

# Your data
Treatments <- c("Missing", "Missing", "0", "3", "0", "4", "0", "0", "2", "7", "0", "4", "2", "2", "0", "0", "0")
Controls <- c("0", "4", "6", "1", "0", "0", "3", "0", "0", "1", "6", "2", "2", "0", "0", "1", "1", "0", "3", "3", "0", "1", "0", "1", "4", "5", "4", "5", "0", "5", "1", "4", "1", "Missing", "5", "4", "2", "1", "5", "6", "2", "6", "1", "5", "0", "5")
length(Treatments)
length(Controls)
Treatments <- c(3, 2, 0, 1, 0, 1, 0, 0, 0, 0, 0)
Controls <- c(2, 3, 2 ,2, 2, 2 ,2 ,3, 2, 3, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1)


Treatments <- c(0,0,0,0)
Controls <- c("Missing", "Missing", "3", "4", "2", "7", "0", "4", "2", "2", "0", "0", "0")
# wilcox t-test
wlcx.test_result <- wilcox.test(Treatments, Controls)
# Convert to numeric, change "Missing" to NA
Treatments <- as.numeric(ifelse(Treatments == "Missing", NA, Treatments))
Controls <- as.numeric(ifelse(Controls == "Missing", NA, Controls))

# t-test
print(t.test_result)

df <- data.frame(
  value = c(Treatments, Controls),
  group = factor(c(rep('Have Positive compound', length(Treatments)), rep('No positive compound', length(Controls))))
)
# Boxplot with p-value annotation
ggplot(df, aes(x = group, y = value)) +
  stat_boxplot(geom = 'errorbar') +
  geom_boxplot() +
  scale_fill_manual(values = c("darkgreen", "darkcyan"))+
  geom_point(position = position_jitter(w = 0.1, h = 0), alpha = 0.7, size = 2.5) + 
  geom_text(
    aes(x = 1.5, y = max(value, na.rm = TRUE),
        label = paste0("p-value = ", round(wlcx.test_result$p.value, digits = 4))),
    vjust = -1  # adjust this value to move the label up or down
  ) +
  theme_bw() + 
  theme(text=element_text(size=13),
                     axis.title.x=element_text(size=12),
                     panel.spacing = unit(2.5, "lines"),
                     axis.ticks.x=element_blank(),
                     axis.text.x = element_text(size=9, angle=30, vjust = 1, hjust=1),
                     plot.title = element_text(size=15,face="bold",hjust = 0.5),
                     legend.position="right",
                     strip.text = element_text(face = "bold"),
                     panel.margin = unit(0, "lines"))  + 
  #scale_fill_brewer(palette="Dark2") +
  labs(x = "Group", y = "Value", title = "Boxplot of AZA/6MP patients for cohort1")
