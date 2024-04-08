rm(selected_spectra)
rm(selected_feature)
selected_spectra <- xdata_spectra[which(spectraNames(xdata_spectra) == "F07.S0200")]
selected_spectra <- xdata_spectra[which(spectraNames(xdata_spectra) == "F10.S1421")]
selected_spectra <- xdata_spectra[which(spectraNames(xdata_spectra) == "F04.S2013")]
mz_value <- c(round(precursorMz(selected_spectra)[[1]],digits = 0)-10,round(precursorMz(selected_spectra)[[1]],digits = 0)+10)
mz_value <- c(140,160)
rt_value <- round(rtime(selected_spectra)[[1]], digits = 1)
selected_feature <- xdata
selected_feature <- filterMz(xdata, mz = mz_value)
selected_feature <- filterRt(selected_feature, rt = rt_value)
#selected_feature <- findChromPeaks(selected_feature, param = cwp)
#selected_feature <- groupChromPeaks(selected_feature, param = pdp)
??filterChromPeaks()
abundance_matrix <- featureValues(selected_feature)
dim(abundance_matrix)

# cohort1的vial_name里面的仅仅是检测标签名字而不是abundance matrix的文件名
# cohort2 不需要改变
cols <- sapply(colnames(abundance_matrix), function(x){
  strsplit(x, split = "_")[[1]][1]
})
colnames(abundance_matrix) <- cols

# 从 abundance_matrix 中提取两组的数据
#cohort1
treatments_data <- abundance_matrix[, treatments$vial_name[which(treatments$vial_name %in% colnames(abundance_matrix))]]
controls_data <- abundance_matrix[, controls$vial_name[which(controls$vial_name %in% colnames(abundance_matrix))]]
#cohort2
treatments_data <- abundance_matrix[, treatments$Metabolomics_FileName[which(treatments$Metabolomics_FileName %in% colnames(abundance_matrix))]]
controls_data <- abundance_matrix[, controls$Metabolomics_FileName[which(controls$Metabolomics_FileName %in% colnames(abundance_matrix))]]

abundance_matrix[1,]
colnames(treatments_data)

# 将NA值转变为0以方便做test

#abundance_matrix[is.na(abundance_matrix)] <- 0
abundance_matrix <- abundance_matrix %>%
  as.data.frame() %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  as.matrix()

# 检查每个特征在两组之间的差异
p_values <- apply(abundance_matrix, 1, function(row_data) {
  print("one run")
  wilcox.test(row_data[colnames(treatments_data)], row_data[colnames(controls_data)], na.rm = T)$p.value
})
wilcox.test(abundance_matrix[1, colnames(treatments_data)], abundance_matrix[1, colnames(controls_data)], na.rm = T)
# 查看显著差异的特征
significant_features <- which(p_values < 0.05)

# 修改格式并选择缺失值较少的行作为样本来画图
significant_abundance_matrix <- abundance_matrix[significant_features,]
significant_abundance_matrix <- significant_abundance_matrix[rowSums(significant_abundance_matrix != 0) > 10, ]

nrow(significant_abundance_matrix)
# 将选择的行修改为长数据
# 预处理数据
# 合并 `controls` 和 `treatments` 到一个数据框中
all_samples <- data.frame(vial_name = c(controls$vial_name, treatments$vial_name),
                          group = c(rep("Control", length(controls$vial_name)), 
                                    rep("Treatment", length(treatments$vial_name))))

# 将数据转换为长格式 #3, 11
long_data <- significant_abundance_matrix[3, ] %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(all_samples, by = c("sample_id" = "vial_name")) %>%
  gather(key = "Feature", value = "Abundance", -sample_id, -group)

# 绘图
p <- ggplot(long_data, aes(x = group, y = Abundance)) + 
  geom_boxplot(lwd = 1.5, outlier.shape = NA) +
  geom_jitter(size = 3, aes(color = group), width = 0.1, height = 0) +
  labs(title = "Boxplot and Scatterplot for Feature 3",
       y = "Abundance") +
  theme_minimal()

# 计算Wilcox测试的p值
wilcox_test <- wilcox.test(Abundance ~ group, data = long_data)
wilcox_test$p.value
# 在图上添加p值
y_max <- max(significant_abundance_matrix[3,])
p <- p + geom_segment(aes(x = 1, xend = 2, y = y_max + 0.05 * y_max, yend = y_max + 0.05 * y_max), color = "black") +
  geom_text(aes(x = 1.5, y = y_max + 0.1 * y_max, label = paste("p =", wilcox_test$p.value)), color = "black", size = 8) +
  theme(text=element_text(size=20),
        legend.position="right",
        strip.text = element_text(face = "bold"))

plot(p)


## 写成一个函数
# cohort1
plot_MS2_cohort1 <- function(matrix){
  for (i in 1:nrow(matrix)) {
    long_data <- significant_abundance_matrix[i, ] %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("sample_id") %>%
      left_join(all_samples, by = c("sample_id" = "vial_name")) %>%
      gather(key = "Feature", value = "Abundance", -sample_id, -group)
    
    # 绘图
    p <- ggplot(long_data, aes(x = group, y = Abundance)) + 
      geom_boxplot(lwd = 1.5, outlier.shape = NA) +
      geom_jitter(size = 3, aes(color = group), width = 0.1, height = 0) +
      labs(title = paste0("Boxplot and Scatterplot for Feature ", i),
           y = "Abundance") +
      theme_minimal()
    
    # 计算Wilcox测试的p值
    wilcox_test <- wilcox.test(Abundance ~ group, data = long_data)
    wilcox_test$p.value
    # 在图上添加p值
    y_max <- max(significant_abundance_matrix[i,])
    p <- p + geom_segment(aes(x = 1, xend = 2, y = y_max + 0.05 * y_max, yend = y_max + 0.05 * y_max), color = "black") +
      geom_text(aes(x = 1.5, y = y_max + 0.1 * y_max, label = paste("p =", wilcox_test$p.value)), color = "black", size = 8) +
      theme(text=element_text(size=20),
            legend.position="right",
            strip.text = element_text(face = "bold"))
    
    plot(p)
  }
}

# cohort2
plot_MS2_cohort2 <- function(matrix){
  for (i in 1:nrow(matrix)) {
    long_data <- significant_abundance_matrix[i, ] %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("sample_id") %>%
      left_join(all_samples, by = c("sample_id" = "Metabolomics_FileName")) %>%
      gather(key = "Feature", value = "Abundance", -sample_id, -group)
    
    # 绘图
    p <- ggplot(long_data, aes(x = group, y = Abundance)) + 
      geom_boxplot(lwd = 1.5, outlier.shape = NA) +
      geom_jitter(size = 3, aes(color = group), width = 0.1, height = 0) +
      labs(title = paste0("Boxplot and Scatterplot for Feature ", i),
           y = "Abundance") +
      theme_minimal()
    
    # 计算Wilcox测试的p值
    wilcox_test <- wilcox.test(Abundance ~ group, data = long_data)
    wilcox_test$p.value
    # 在图上添加p值
    y_max <- max(significant_abundance_matrix[i,])
    p <- p + geom_segment(aes(x = 1, xend = 2, y = y_max + 0.05 * y_max, yend = y_max + 0.05 * y_max), color = "black") +
      geom_text(aes(x = 1.5, y = y_max + 0.1 * y_max, label = paste("p =", wilcox_test$p.value)), color = "black", size = 8) +
      theme(text=element_text(size=20),
            legend.position="right",
            strip.text = element_text(face = "bold"))
    
    plot(p)
  }
}

plot_MS2_cohort2(matrix = significant_abundance_matrix)
plot_MS2_cohort1(matrix = significant_abundance_matrix)
D