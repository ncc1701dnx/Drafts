
### 将需要的sampleID，全部导入成一个ID（还是做了区分因为有些病人ID和sampleID并不是全都需要的：比如6MP只需要UC病人）
all_ID <- c(controls_ID, treatments_ID)
all_MBX_results <- MBX_results[, c(1:7, which(colnames(MBX_results) %in% all_ID))]
dim(all_MBX_results)

### MBX results 用ion charge mode分类，因为adducts需要计算不同的ion mode
MBX_C18neg_samples <- subset(MBX_results, Method == "C18-neg" )
MBX_HILneg_samples <- subset(MBX_results, Method == "HILIC-neg")
MBX_C18pos_samples <- subset(MBX_results, Method == "C8-pos")
MBX_HILpos_samples <- subset(MBX_results, Method == "HILIC-pos")
class(MBX_C18neg_samples)

### 合并成一个表
dataset_list <- list(
  "MBX_C18neg" = MBX_C18neg_samples,
  "MBX_HILneg" = MBX_HILneg_samples,
  "MBX_C18pos" = MBX_C18pos_samples,
  "MBX_HILpos" = MBX_HILpos_samples
)
names(dataset_list["MBX_C18neg"])

IBDMDB_binary_check <- function(formula, popular = T, ppm = 10){
  ## 这个函数将对于给出化合物的所有
  
  ################# 首先是用 ion charge为negative来跑两次 ####################
  
  ion_charge_mode = "neg"
  print("跑negative ion mode")
  
  ### 先跑C18的数据
  print("跑C18 neg的资料")
  
  subset <- MBX_C18neg_samples 
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode, popular = popular, ppm = ppm)
  
  for (x in 1:nrow(com_range)) {
    print(paste0("x=",x))
    #x =1
    tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
    if(length(tt) == 0){next}
    
    results <- subset[tt, ]
  }
  
  
  ### 然后跑HIL的数据
  print("跑HIL neg的数据")
  subset <- MBX_HILneg_samples
  for (x in 1:nrow(com_range)) {
    print(paste0("x=",x))
    tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
    if(length(tt) == 0){next}
    
    results <- rbind(results, subset[tt, ])
  }
  
  ################# 然后是用 ion charge为positive来跑两次 ####################
  ion_charge_mode = "pos"
  print("跑positive ion mode")
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode, popular = popular, ppm = ppm)
  
  ### 先跑C18的数据
  print("跑C18 pos的资料")
  subset <- MBX_C18pos_samples ## 记得在这里和下面都要写上两个【【】】，因为是提取列表
  
  for (x in 1:nrow(com_range)) {
    x = 1
    print(paste0("x=",x))
    tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
    if(length(tt) == 0){next}
    
    results <- rbind(results, subset[tt, ])
  }
  
  ### 然后跑HIL的数据
  subset_list <- subset_list_HIL
  print("跑HIL pos的数据")
  subset <- MBX_HILpos_samples
  
  for (x in 1:nrow(com_range)) {
    #x = 4
    print(paste0("x=",x))
    tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
    if(length(tt) == 0){next}
    
    results <- rbind(results, subset[tt, ])
  }
  
  results
}

# 画散点图专用函数
plot_and_save <- function(data, row_num) {
  p <- ggplot(data, aes(x = variable, y = value)) +
    #geom_boxplot(outlier.shape = NA, lwd = 1.5) + 
    #geom_jitter(position=position_jitter(0.2), size=4) +
    geom_jitter(size=2, width = 0.1, height = 0) +
    xlab("Samples") +
    ylab("Log e Abundance") +
    ggtitle(paste0("Observation", row_num)) +
    theme(text=element_text(size=22),
          legend.position="right",
          strip.text = element_text(face = "bold"))+
    theme_minimal()
  
  # 你可以使用ggsave将每张图保存到文件中
  # ggsave(filename = paste0("Observation_", row_num, ".png"), plot = p)
  print(p)
}

# 主循环
for (i in 1:nrow(binaryTable5ASA)) {
  #i = 44
  current_row <- binaryTable5ASA[i, -c(1:7)]
  current_row <- log(current_row)
  long_data <- as.data.frame(t(current_row))
  colnames(long_data) <- c("value")
  long_data$variable <- rownames(long_data)
  
  plot_and_save(long_data, i)
}

# 将循环写成函数
IBDMDB_plot_binary <- function(binaryTable){
  for (i in 1:nrow(binaryTable)) {
    #binaryTable = binaryTable6MP
    #i = 10
    print(i)
    current_row <- binaryTable[i, -c(1:7)]
    current_row <- log(current_row)
    long_data <- as.data.frame(t(current_row))
    colnames(long_data) <- c("value")
    long_data$variable <- rownames(long_data)
    
    plot_and_save(long_data, i)
  }
}


### 实战阶段

## 先用5ASA确认下函数有用
controls_ID <- meta_data$External.ID[which(meta_data$any5asa== 0)]
treatments_ID <- meta_data$External.ID[which(meta_data$any5asa == 1)]
length(controls_ID) ## 676
length(treatments_ID) ## 174
all_ID <- c(controls_ID, treatments_ID)
length(all_ID)

binaryTable5ASA <- IBDMDB_binary_check("C7H7N1O3", popular = F)
dim(binaryTable5ASA)
which(binaryTable5ASA$Metabolite == "154.0502_3.83")#第44号
IBDMDB_plot_binary(binaryTable5ASA)

binaryTableNA5ASA <- IBDMDB_binary_check("C9H9N1O4", popular = F)
dim(binaryTableNA5ASA)
which(binaryTableNA5ASA$Metabolite == "196.0609_2.81")# 34
which(binaryTableNA5ASA$Metabolite == "196.0602_5.19")# 36
which(MBX_results$Metabolite == "196.0579_7.11")
IBDMDB_plot_binary(binaryTableNA5ASA)

### 使用6MP进行调查
controls_ID <- meta_data$External.ID[which(meta_data$aza6mp == 0)]
treatments_ID <- meta_data$External.ID[which(meta_data$aza6mp == 1)]
length(controls_ID) ## 772
length(treatments_ID) ## 78
all_ID <- c(controls_ID, treatments_ID)
length(all_ID)

binaryTable6MP <- IBDMDB_binary_check("C5H4N4S1", popular = F, ppm = 50)
dim(binaryTable6MP)
IBDMDB_plot_binary(binaryTable6MP)
binaryTable6MP[81,]

binaryTable6TIMP <- IBDMDB_binary_check("C10H13N4O7P1S1", popular = F)
dim(binaryTable6TIMP)
IBDMDB_plot_binary(binaryTable6TIMP)

binaryTable6TUA <- IBDMDB_binary_check("C5H4N4O2S1", popular = F, ppm = 50)
dim(binaryTable6TUA)
IBDMDB_plot_binary(binaryTable6TUA)
