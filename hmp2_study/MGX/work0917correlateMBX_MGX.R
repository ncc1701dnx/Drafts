### 1：首先进行分组：
### Treatment and Control
controls_ID <- meta_data$External.ID[which(meta_data$aza6mp == 0)]
treatments_ID <- meta_data$External.ID[which(meta_data$aza6mp == 1)]
length(controls_ID) ## 772
length(treatments_ID) ## 78
length(unique(meta_data$Participant.ID[which(meta_data$aza6mp == 0)])) ##90
length(unique(meta_data$Participant.ID[which(meta_data$aza6mp == 1)])) ##13

### IBD and non-IBD
controls_ID <- meta_data$External.ID[which(meta_data$diagnosis == "nonIBD")]
treatments_ID <- meta_data$External.ID[which(meta_data$diagnosis == "CD" | meta_data$diagnosis == "UC")]

### 2：然后处理成function所需input：
controls_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% controls_ID))]
treatments_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% treatments_ID))]
# 现在将MBX结果分割成不同的质谱处理模式
## Set negative and positive to analyze different adducts
MBX_C18neg_treatments <- subset(treatments_MBX_results, Method == "C18-neg" )
MBX_HILneg_treatments <- subset(treatments_MBX_results, Method == "HILIC-neg")
MBX_C18pos_treatments <- subset(treatments_MBX_results, Method == "C8-pos")
MBX_HILpos_treatments <- subset(treatments_MBX_results, Method == "HILIC-pos")
table(controls_MBX_results$Method)
MBX_C18neg_controls <- subset(controls_MBX_results, Method == "C18-neg")
MBX_HILneg_controls <- subset(controls_MBX_results, Method == "HILIC-neg")
MBX_C18pos_controls <- subset(controls_MBX_results, Method == "C8-pos")
MBX_HILpos_controls <- subset(controls_MBX_results, Method == "HILIC-pos")
## 创建一个data list，包含所有需要的MBX result 分表格
dataset_list <- list(
  "MBX_C18neg_treatments" = MBX_C18neg_treatments,
  "MBX_HILneg_treatments" = MBX_HILneg_treatments,
  "MBX_C18pos_treatments" = MBX_C18pos_treatments,
  "MBX_HILpos_treatments" = MBX_HILpos_treatments,
  "MBX_C18neg_controls" = MBX_C18neg_controls,
  "MBX_HILneg_controls" = MBX_HILneg_controls,
  "MBX_C18pos_controls" = MBX_C18pos_controls,
  "MBX_HILpos_controls" = MBX_HILpos_controls
)

### 接着就可以做区分compound的画图分析工作了
## 6TIMP
IBDMDB_6TIMP <- IBDMDB_MBX_check(formula = "C10H13N4O7P1S1", popular = F)
IBDMDB_6TIMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TIMP)
IBDMDB_MBX_plot(IBDMDB_6TIMP, "6TIMP")
IBDMDBresults_boxplot(IBDMDB_6TIMP, "6TIMP")
IBDMDB_compound <- IBDMDB_6TIMP
compound_name <- "6TIMP"
rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_result <- dplyr::inner_join(rr_result, MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
i = 1
rr_result <- rr_result[i, -(1:9)]
length(which(colnames(rr_result) %in% GMPS_Treatment$sampleID))
colnames(rr_result)
dim(rr_result)

#small_GMPS_MGX <- filtered_data_GMPS[, -c(1:16)]
colnames(GENE_Treatment)
small_GMPS_MGX <- GMPS_Treatment[, -8] # -"Participant.ID"

rr_long <- rr_result %>%
  gather(key = "sampleID", value = "MBX_abundance")

merged_data_GMPS_6TGTP <- small_GMPS_MGX %>%
  left_join(rr_long, by = "sampleID")

merged_data_GMPS_6TGTP <- merged_data_GMPS_6TGTP %>%
  filter(!is.na(rpkm_value) & !is.na(MBX_abundance))

unique(merged_data_GMPS_6TGTP$sampleID)
dim(merged_data_GMPS_6TGTP)

cor_result <- cor(merged_data_GMPS_6TGTP$rpkm_value, merged_data_GMPS_6TGTP$MBX_abundance, method = "pearson")
print(cor_result)

plot_pearson <- ggplot(merged_data_GMPS_6TGTP, aes(x = rpkm_value, y = MBX_abundance)) +
  geom_point(aes(color = aza6mp), alpha = 0.5) +  # 将点的颜色映射到aza6mp列
  geom_smooth(method = "lm", se = FALSE, color = "#C4961A") +  # 添加线性回归线
  theme_minimal() +
  labs(title = "Scatter plot with Regression Line",
       x = "RPKM Value",
       y = "MBX Abundance",
       caption = paste("Pearson correlation: ", round(cor_result, 2))) +
  theme(text=element_text(size=30), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold")) +
  theme_bw()

print(plot_pearson)

# 6TUA
IBDMDB_6TUA <- IBDMDB_MBX_check(formula = "C5H4N4O2S1", popular = F)
IBDMDB_6TUA$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TUA)
IBDMDB_MBX_plot(IBDMDB_6TUA, "6TUA")
IBDMDBresults_boxplot(IBDMDB_6TUA, "6TUA")

IBDMDB_compound <- IBDMDB_6TUA
compound_name <- "6TUA"
rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_result <- dplyr::inner_join(rr_result, MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
i = 1
rr_result <- rr_result[i, -(1:9)]
length(which(colnames(rr_result) %in% GENE_Treatment$sampleID))
colnames(rr_result)
dim(rr_result)

#small_GMPS_MGX <- filtered_data_GMPS[, -c(1:16)]
colnames(GENE_Treatment)
small_GENE_MGX <- GENE_Treatment[, -c(1,2,11)] # -"Participant.ID", "gc_id", "msp_id"

rr_long <- rr_result %>%
  gather(key = "sampleID", value = "MBX_abundance")

merged_data <- small_GENE_MGX %>%
  left_join(rr_long, by = "sampleID")

merged_data <- merged_data %>%
  filter(!is.na(rpkm_value) & !is.na(MBX_abundance))

unique(merged_data$sampleID)
dim(merged_data)

cor_result <- cor(merged_data$rpkm_value, merged_data$MBX_abundance, method = "pearson")
print(cor_result)

plot_pearson <- ggplot(merged_data, aes(x = rpkm_value, y = MBX_abundance)) +
  #geom_point(aes(color = aza6mp), alpha = 0.5) +  # 将点的颜色映射到aza6mp列
  geom_point(color = "#20793c", alpha = 10, size = 3) +  # 将点的颜色映射到aza6mp列
  geom_smooth(method = "lm", se = FALSE, color = "#C4961A") +  # 添加线性回归线
  theme_minimal() +
  labs(title = "Scatter plot with Regression Line",
       x = "RPKM Value",
       y = "MBX Abundance",
       caption = paste("Pearson correlation: ", round(cor_result, 2)), size = 30) +
  theme(axis.title.x = element_blank())+  # 隐藏X轴的标签
  theme(axis.title.y = element_blank())+
  #theme(axis.text.y = element_text(size = 15, face = "bold")) +
  theme(text=element_text(size=30), #调整字号
        #legend.position="right",
        strip.text = element_text(face = "bold"))

print(plot_pearson)


# 6-TGTP
IBDMDB_6_TGTP <- IBDMDB_MBX_check(formula = "C10H16N5O13P3S1", popular = F)
IBDMDB_6_TGTP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6_TGTP)
IBDMDB_MBX_plot(IBDMDB_6_TGTP, "6TGTP")
IBDMDBresults_boxplot(IBDMDB_6_TGTP, "6TGTP")



IBDMDB_compound <- IBDMDB_6TUA
compound_name <- "6TUA"
rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_result <- dplyr::inner_join(rr_result, MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
i = 1
rr_result <- rr_result[i, -(1:9)]
length(which(colnames(rr_result) %in% GENE_Treatment$sampleID))
colnames(rr_result)
dim(rr_result)

#small_GMPS_MGX <- filtered_data_GMPS[, -c(1:16)]
colnames(GENE_Treatment)
small_GENE_MGX <- GENE_Treatment[, -c(1,2,11)] # -"Participant.ID", "gc_id", "msp_id"

rr_long <- rr_result %>%
  gather(key = "sampleID", value = "MBX_abundance")

merged_data <- small_GENE_MGX %>%
  left_join(rr_long, by = "sampleID")

merged_data <- merged_data %>%
  filter(!is.na(rpkm_value) & !is.na(MBX_abundance))

unique(merged_data$sampleID)
dim(merged_data)

cor_result <- cor(merged_data$rpkm_value, merged_data$MBX_abundance, method = "pearson")
print(cor_result)

plot_pearson <- ggplot(merged_data, aes(x = rpkm_value, y = MBX_abundance)) +
  #geom_point(aes(color = aza6mp), alpha = 0.5) +  # 将点的颜色映射到aza6mp列
  geom_point(color = "#20793c", alpha = 10, size = 3) +  # 将点的颜色映射到aza6mp列
  geom_smooth(method = "lm", se = FALSE, color = "#C4961A") +  # 添加线性回归线
  theme_minimal() +
  labs(title = "Scatter plot with Regression Line",
       x = "RPKM Value",
       y = "MBX Abundance",
       caption = paste("Pearson correlation: ", round(cor_result, 2)), size = 30) +
  theme(axis.title.x = element_blank())+  # 隐藏X轴的标签
  theme(axis.title.y = element_blank())+
  #theme(axis.text.y = element_text(size = 15, face = "bold")) +
  theme(text=element_text(size=30), #调整字号
        #legend.position="right",
        strip.text = element_text(face = "bold"))

print(plot_pearson)
