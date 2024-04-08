# This Script is for pre-inspection of IBDMDB MBX data
MBX_results <- read.csv("R_Works/inputs/iHMP_metabolomics.csv", header = T)

treatments <- subset(MBX_IBDMDB_metadata, Mercaptopurine..Purinethol..6MP. == "Current" | Azathioprine..Imuran..Azasan. == "Current")
controls <- subset(MBX_IBDMDB_metadata,Mercaptopurine..Purinethol..6MP. == "Never taken" & Azathioprine..Imuran..Azasan. == "Never taken")

controls_ID <- controls$External.ID
treatments_ID <- treatments$External.ID

treatments_ID %in% colnames(MBX_results)

controls_ID %in% colnames(MBX_results)
print(colnames(MBX_results))

controls_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], controls_ID)]
treatments_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], treatments_ID)]

table(treatments_MBX_results$Method)
MBX_neg_treatments <- subset(treatments_MBX_results, Method == "C18-neg" | Method == "HILIC-neg")
MBX_pos_treatments <- subset(treatments_MBX_results, Method == "C8-pos" | Method == "HILIC-pos")

# Perform t-test for each metabolite peak number
ttest_results <- lapply(names(MBX_results)[-1], function(metabolite) {
  t.test(MBX_results[MBX_results$ID %in% control_group_ids, metabolite], 
         data_df[data_df$ID %in% treatment_group_ids, metabolite])
})

####重要！ 分割表格前一定记得加上行名方便日后对齐
row.names(MBX_results) = MBX_results$Compound

## Set negative and positive to analyze different adducts
controls_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], controls_ID)]
treatments_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], treatments_ID)]
MBX_neg_treatments <- subset(treatments_MBX_results, Method == "C18-neg" | Method == "HILIC-neg")
MBX_pos_treatments <- subset(treatments_MBX_results, Method == "C18-pos" | Method == "HILIC-pos")
MBX_neg_controls <- subset(controls_MBX_results, Method == "C18-neg" | Method == "HILIC-neg")
MBX_pos_controls <- subset(controls_MBX_results, Method == "C18-pos" | Method == "HILIC-pos")

names(MBX_results)[-c(1:7)]
control_group <- MBX_results[, controls_ID]

dim(control_group)
treatment_group <- MBX_results[, treatments_ID]
treatment_group[rownames(MBX_results)[1],]
dim(treatment_group)
log(c(NA,NA,NA,NA))
ttest_all_results_withna <- lapply(rownames(MBX_results), function(metabolite) {
  control_mm <- control_group[metabolite,]
  treatment_mm <- treatment_group[metabolite,]
  control_mm <- log(as.numeric(control_mm))
  treatment_mm <- log(as.numeric(treatment_mm))
  if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {NULL}
  else {t.test(control_mm, treatment_mm, na.action=na.omit)}
})


length(p_values_withna)
ttest_all_results_withna[[1]]$estimate[2]-ttest_all_results_withna[[1]]$estimate[1]
p_values_withna <- unlist(sapply(ttest_all_results_withna, function(result){
  result$p.value
}))
mean_dif_whithna <- unlist(sapply(ttest_all_results_withna, function(result){
  result$estimate[2]/result$estimate[1]
}))

ttest_and_diff_withna <- data.frame(diff = mean_dif_whithna, p_val = p_values_withna)
ttest_and_diff_withna$significance <- ifelse(log2(ttest_and_diff_withna$diff) > 1 &
                                               -log10(ttest_and_diff_withna$p_val) > 2, "Significantly Abundant", "No Significance")
table(ttest_and_diff_withna$signif)
length(which(abs(log2(ttest_and_diff_withna$diff)) > 1 & -log10(ttest_and_diff_withna$p_val) > 2))
## only 156 metabolites are significantly abundant in treatment group
## But 3199 in total significantly abundant
## Means 3043 significantly abundant in control group

# rownames(ttest_and_diff) <- rownames(MBX_results)
ggplot(ttest_and_diff_withna, aes(x = log2(diff), y = -log10(p_val), col = significance)) +
  geom_point() +
  scale_colour_brewer(palette = "Set2")+
  xlab("Difference in means (Log2 Fold Change)") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot") +
  theme(text=element_text(size=20),
        legend.position="right",
        strip.text = element_text(face = "bold"))
  #theme_bw()

## 因为很多metabolite是某个group的专属，所以现在设置NA为 0.5*Min_Value
all(is.na(control_group)[1,])
dim(control_group[1,])
ttest_all_results <- lapply(rownames(MBX_results), function(metabolite) {
  control_mm <- control_group[metabolite,]
  treatment_mm <- treatment_group[metabolite,]
  control_mm <- as.numeric(control_mm)
  treatment_mm <- as.numeric(treatment_mm)
  if(all(is.na(control_mm)) & (length(which(!is.na(treatment_mm))) >= 4)){
   control_mm = 0.5*min(treatment_mm, na.rm = T)
   treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)}
  if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {NULL}
  else {t.test(control_mm, treatment_mm, na.action=na.omit)}
})


length(which(ttest_and_diff$signif == "Significantly_Abundant"))
length(p_values)
ttest_all_results[[1]]$estimate[2]-ttest_all_results[[1]]$estimate[1]
p_values <- unlist(sapply(ttest_all_results, function(result){
  result$p.value
}))
mean_dif <- unlist(sapply(ttest_all_results, function(result){
  result$estimate[2]/result$estimate[1]
}))

ttest_and_diff <- data.frame(diff = mean_dif, p_val = p_values)
ttest_and_diff$significant <- ifelse(log2(ttest_and_diff$diff) > 1 & 
                                       -log10(ttest_and_diff$p_val) > 2, "Significantly Abundant", "No Significance")
table(ttest_and_diff$signif)
length(which(abs(log2(ttest_and_diff$diff)) > 1.5 & -log10(ttest_and_diff$p_val) > 2))
ggplot(ttest_and_diff, aes(x = log2(diff), y = -log10(p_val), col = signif)) +
  geom_point() +
  scale_colour_brewer(palette = "Set2")+
  xlab("Difference in means (Log2 Fold Change)") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot") +
  theme(text=element_text(size=20),
        legend.position="right",
        strip.text = element_text(face = "bold"))
  #theme_bw()

