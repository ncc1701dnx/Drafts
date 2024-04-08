
# 首先是一个函数，用来从xdata中计算目标化合物的全部m/z区间并合并成表
cal_mz_range_filter_feature <- function(formula, ion_charge_mode = "pos", popular = F, ppm =10){
  ads <- cal_direct_range(formula = formula,  ion_charge_mode = ion_charge_mode, popular = popular, ppm = ppm)
  PossibleFeatures_list <- list()
  for (i in 1:nrow(ads)) {
    print(i)
    isError <- FALSE
    ad <- ads[i,]
    mz_value <- c(ad$LowerBound, ad$UpperBound)
    selected_feature <- filterMz(xdata, mz = mz_value)
    tryCatch({
      # featureValues code here. May potentially throw an error when no feature found in this range
      abundance_matrix <- featureValues(selected_feature)  
      # Additional processing if the function call is successful.
      PossibleFeatures_list[[i]] <- abundance_matrix
    }, error = function(e) {
      print(paste0("error", " in ",i ))
      isError <- TRUE 
    })
    if (isError) {
      next
    }
  }
  PossibleFeatures <- do.call(rbind, PossibleFeatures_list)
  PossibleFeatures
}

# 使用例
RobertFeatures5ASA <- cal_mz_range_filter_feature(formula = "C7H7N1O3",  ion_charge_mode = "pos", popular = F, ppm = 1000)

### 接下来的几步需要手动执行，并修改参数
# 将NA值转变为0以方便做test
#
RobertFeatures <- RobertFeatures5ASA

RobertFeatures <- RobertFeatures %>%
  as.data.frame() %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  as.matrix()

# cohort1的vial_name里面的仅仅是检测标签名字而不是abundance matrix的文件名
# cohort2 不需要改变
cols <- sapply(colnames(RobertFeatures), function(x){
  strsplit(x, split = "_")[[1]][1]
})
colnames(RobertFeatures) <- cols

# 因为不是全部患者都是用药用户，所以我们subset一下用药患者
RobertFeatures <- RobertFeatures[, treatments$vial_name[which(treatments$vial_name %in% colnames(RobertFeatures))]]

remissions <- treatments[which(treatments$histologic_remission == 1), ]
noremissions <- treatments[which(treatments$histologic_remission == 0), ]

remission_data <- RobertFeatures[, remissions$vial_name[which(remissions$vial_name %in% colnames(RobertFeatures))]]
noremission_data <- RobertFeatures[, noremissions$vial_name[which(noremissions$vial_name %in% colnames(RobertFeatures))]]

## 合并所有数据到一个大dataframe
all_samples <- data.frame(vial_name = c(remissions$vial_name, noremissions$vial_name),
                          group = c(rep("remissions", length(remissions$vial_name)), 
                                    rep("noremissions", length(noremissions$vial_name))))

## 接下来的函数就是将上述结果合并并画散点图计算p值
# cohort1
plot_MS1_cohort1 <- function(matrix){
  for (i in 1:nrow(matrix)) {
    #i =32
    print(i)
    long_data <- matrix[i, ] %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("sample_id") %>%
      left_join(all_samples, by = c("sample_id" = "vial_name")) %>%
      gather(key = "Feature", value = "Abundance", -sample_id, -group)
    
    # 绘图
    p <- ggplot(long_data, aes(x = group, y = Abundance)) + 
      #geom_boxplot(lwd = 1.5, outlier.shape = NA) +
      geom_jitter(size = 4, aes(color = group), width = 0.1, height = 0) +
      labs(title = paste0("Boxplot and Scatterplot for Feature ", i),
           y = "Abundance") +
      theme_minimal()
    
    # 计算Wilcox测试的p值
    wilcox_test <- wilcox.test(Abundance ~ group, data = long_data)
    wilcox_test$p.value
    # 在图上添加p值
    y_max <- max(matrix[i,])
    p <- p + geom_segment(aes(x = 1, xend = 2, y = y_max + 0.05 * y_max, yend = y_max + 0.05 * y_max), color = "black") +
      geom_text(aes(x = 1.5, y = y_max + 0.1 * y_max, label = paste("p =", wilcox_test$p.value)), color = "black", size = 8) +
      theme(text=element_text(size=25, face = "bold"),
            legend.position="right",
            strip.text = element_text(face = "bold"))
    
    plot(p)
  }
}



# 使用例：
plot_MS1_cohort1(matrix = RobertFeatures)

# 6-MP
RobertFeatures6MP <- cal_mz_range_filter_feature(formula = "C5H4N4S1",  ion_charge_mode = "pos", popular = F, ppm = 1000)
dim(RobertFeatures6MP)
RobertFeatures <- RobertFeatures6MP
plot_MS1_cohort1(matrix = RobertFeatures)
