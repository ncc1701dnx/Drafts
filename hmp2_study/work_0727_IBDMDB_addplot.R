### 6TGMP
IBDMDB_compound <- IBDMDB_6_TGTP
compound_name <- "6TGTP"
  rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
  rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  which(rr_result$RT == 10.23)
  i =1
  as.numeric(rr_control[i, -(1:9)])
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 100000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 100000)+9)]/12)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50000)+9)]/3)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 10000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 10000)+9)]/5)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 5000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 5000)+9)]/600)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 3000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 3000)+9)]/400)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 1000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 1000)+9)]/300)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 300)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 300)+9)]/10)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50)+9)]/4)
  
### 5ASA 
  IBDMDB_compound <- IBDMDB_5_ASA
  compound_name <- "5ASA"
  rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
  rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  i =1
  as.numeric(rr_control[i, -(1:9)])
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 100000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 100000)+9)]/600)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50000)+9)]/300)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 10000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 10000)+9)]/500)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 5000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 5000)+9)]/600)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 3000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 3000)+9)]/400)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 1000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 1000)+9)]/300)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 300)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 300)+9)]/10)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50)+9)]/4)
### N acetyl 5-ASA 
  IBDMDB_compound <- IBDMDB_NAcetyl_5_ASA
  compound_name <- "N acetyl 5-ASA"
  rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
  rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  which(rr_result$RT == 6.19)
  i =3
  as.numeric(rr_control[i, -(1:9)])
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 100000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 100000)+9)]/600)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 50000)+9)]/300)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 10000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 10000)+9)]/100)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 5000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 5000)+9)]/600)
  rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 3000)+9)] <- trunc(rr_control[i,(which(as.numeric(rr_control[i,-(1:9)]) > 3000)+9)]/300)
  
### 按照师兄的说法，我们用我们的workflow画一下X154.0502_3.83这个峰
X154 <- MBX_results[which(MBX_results$Metabolite == "154.0502_3.83"), -c(1:7)]

rr_treatment <- X154[,which(colnames(X154) %in% treatments_ID)]
rr_control <- X154[,which(colnames(X154) %in% controls_ID)]
rr_treatment <- unlist(rr_treatment)
rr_control <- unlist(rr_control)

  ###### Start of Plot Preparation ####
  for(i in 1:nrow(rr_result)){
    # 这里是后面画图的时候图片的名字
    i =5
    print(i)
    plot_title <- IBDMDB_compound[which(IBDMDB_compound$significance == "Significantly Abundant"),][1:4]
    pt <- paste(as.character(plot_title[i,]), collapse = " ")
    #pt <- "X154.0502_3.83"
    # 数据转换，log transform， Na eliminate，一如之前的函数
    control_mm <- rr_control[i,-(1:9)]
    #control_mm <- as.numeric(rr_control)
    control_mm_2 <- control_mm
    treatment_mm <- rr_treatment[i,-(1:9)]
    #treatment_mm <- as.numeric(rr_treatment)
    if(any(is.na(control_mm)) == T & all(is.na(control_mm)) == F ){
      control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
    }
    if(all(is.na(control_mm)) & all(is.na(treatment_mm)) == F){
      control_mm[which(is.na(control_mm))] = 0.5*min(treatment_mm, na.rm = T)
      treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)
    }
    if(any(is.na(treatment_mm)) == T){
      treatment_mm[which(is.na(treatment_mm))] = 0.5*min(control_mm_2, na.rm = T) 
    }
    control_mm <- log(control_mm)
    treatment_mm <- log(treatment_mm)
    if_paper_treatment <- c()
    if_paper_treatment <- ifelse(colnames(control_mm) %in% paper_treatment_5ASA, "Yes", "No")
    if_paper_treatment <- c(if_paper_treatment, ifelse(colnames(treatment_mm) %in% paper_treatment_5ASA, "Yes", "No"))
    # 重新计算p value
    #t_result <- wilcox.test(control_mm, treatment_mm)
    control_mm <- as.numeric(control_mm)
    treatment_mm <- as.numeric(treatment_mm)
    t_result <- t.test(control_mm, treatment_mm)
    p_value <- t_result$p.value
    
    # 建立ggplot2主表格
    combined_data <- data.frame(
      abundance = c(control_mm, treatment_mm),
      group = c(rep("Control", length(control_mm)),rep("Treatment", length(treatment_mm)))
    )
    combined_data$if_paper_treat <- if_paper_treatment
    
    #设立因子级别
    combined_data$group <- factor(combined_data$group, levels = c("Control", "Treatment"))
    
    # ggplot2 画画图 
    # 根据p值确定标签的显示方式
    label_p <- ifelse(p_value < 0.001, "< 0.001", formatC(p_value, format = "e", digits = 2))
    
    # 手动添加色彩
    #my_colors <- brewer.pal(11, "Spectral")
    #my_colors <- my_colors[c( 9, 10, 11, 7)] # 手动定义颜色
    #my_colors <- c("#C4961A", "steelblue")
    # IBDMDB_newuser_5ASA
    #my_colors <- c("cyan4", "royalblue4")
    my_colors <- c("#00798c", "#66a182")
    
    colnames(combined_data)
    # 使用ggplot2绘图
    p <- ggplot(combined_data, aes(x=group, y=abundance, color=if_paper_treat)) + 
      geom_boxplot(outlier.shape = NA, lwd = 1.5) + 
      geom_jitter(size=4, width = 0.1, height = 0) +
      scale_color_manual(values = my_colors) +
      xlab("Groups") +
      ylab("Log e Abundance") +
      theme_minimal()
    
    # 计算两个组的数据的最大值，以确定标注线的位置
    y_max <- max(combined_data$abundance)
    
    # 在图上增加一个标注线和文本
    p <- p + 
      geom_segment(aes(x = 1, xend = 2, y = y_max + 0.05 * y_max, yend = y_max + 0.05 * y_max), color = "black") +
      geom_text(aes(x = 1.5, y = y_max + 0.1 * y_max, label = paste("p =", p_value)), color = "black", size = 8) +
      ggtitle(paste0(compound_name, " ", pt))+
      theme(text=element_text(size=22),
            legend.position="right",
            strip.text = element_text(face = "bold"))
    plot(p)
  }

  #### 5 asa 作者图作图
  # 加载必要的库
  library(ggplot2)
  
  # 创建数据框
  data_5asa <- data.frame(
    binary154 = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    Participant.ID = c('C3004', 'C3004', 'C3031', 'C3031', 'H4014', 'H4014', 'H4015', 'H4015', 'H4035', 'H4035', 'H4040', 'H4040', 'M2008', 'M2008', 'M2028', 'M2028', 'M2071', 'M2071', 'P6009', 'P6009', 'P6010', 'P6010', 'P6012', 'P6012', 'P6016', 'P6016'),
    visit_num = c(13, 20, 11, 14, 13, 19, 8, 13, 23, 29, 18, 29, 4, 9, 13, 19, 12, 16, 4, 25, 8, 12, 4, 7, 6, 11),
    SampleID = c('CSM5MCXL', 'CSM67UDN', 'CSM79HRG', 'CSM7KOO9', 'HSM6XRRV', 'HSM6XRVO', 'HSM5MD73', 'HSM6XRS8', 'HSMA33OZ', 'HSMA33MI', 'HSMA33OJ', 'HSMA33M8', 'MSM5LLDI', 'MSM5LLDS', 'MSM6J2IG', 'MSM6J2Q3', 'MSMA26AZ', 'MSMB4LZ4', 'PSM6XBRK', 'PSM7J1CU', 'PSM6XBSK', 'PSM6XBUG', 'PSM6XBSE', 'PSM6XBVM', 'PSM7J19B', 'PSM7J19J'),
    X154.0502_3.83 = c(283552, 507071413, 172164, 129818420, 290011, 1530601385, 168868, 3121179635, 2352694, 82747705, 390204, 1468731677, 1773492, 643195341, 557418, 109827300, 1036302, 323254933, 169904, 456525486, 254310, 879384985, 1059704, 12783515, 268885, 1171158465)
  )
  
  data_5asa$binary154 <- ifelse(data_5asa$binary154 == 1, "Post", "Pre")
  data_5asa$binary154 <- as.factor(data_5asa$binary154)
  # 使用ggplot2绘制图形
  p <- ggplot(data_5asa, aes(x = factor(binary154), y = X154.0502_3.83, color = factor(binary154))) +
    geom_boxplot(lwd = 1.5) +
   geom_point(position = position_dodge(width = 0.75), size = 5) + # 不使用jitter
    xlab("Groups") +
    ylab("Log e Abundance") +
    labs(color = "Groups")+
    theme_minimal()+
    ggtitle("5-ASA")+
    theme(text=element_text(size=22),
          legend.position="right",
          strip.text = element_text(face = "bold"))
  print(p)
  
  
  ### 6TGTP
  IBDMDB_compound <- IBDMDB_6_TGTP
  compound_name <- "6TGTP"
  rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
  rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  rr_result <- dplyr::inner_join(rr_result, MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  i = 1
  rr_result <- rr_result[i, -(1:9)]
  length(which(colnames(rr_result) %in% filtered_data_GMPS$sampleID))
  colnames(rr_result)
  
small_GMPS_MGX <- filtered_data_GMPS[, -c(1:16)]
small_GMPS_MGX <- small_GMPS_MGX[, -8] # -"Participant.ID"

rr_long <- rr_result %>%
  gather(key = "sampleID", value = "MBX_abundance")
  
merged_data_GMPS_6TGTP <- small_GMPS_MGX %>%
  left_join(rr_long, by = "sampleID")

merged_data_GMPS_6TGTP <- merged_data_GMPS_6TGTP %>%
  filter(!is.na(rpkm_value) & !is.na(MBX_abundance))

cor_result <- cor(merged_data_GMPS_6TGTP$rpkm_value, merged_data_GMPS_6TGTP$MBX_abundance, method = "pearson")
print(cor_result)

plot_pearson <- ggplot(merged_data_GMPS_6TGTP, aes(x = rpkm_value, y = MBX_abundance)) +
  geom_point(aes(color = aza6mp), alpha = 0.5) +  # 将点的颜色映射到aza6mp列，根据你的数据情境你可以选择是否这样做
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # 添加线性回归线
  theme_minimal() +
  labs(title = "Scatter plot with Regression Line",
       x = "RPKM Value",
       y = "MBX Abundance",
       caption = paste("Pearson correlation: ", round(cor_result, 2)))

print(plot_pearson)


 ## 6TIMP
IBDMDB_compound <- IBDMDB_6TIMP
compound_name <- "6TIMP"
rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_result <- dplyr::inner_join(rr_result, MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
i = 2
rr_result <- rr_result[i, -(1:9)]
length(which(colnames(rr_result) %in% GMPS_Treatment$sampleID))
colnames(rr_result)
dim(rr_result)

#small_GMPS_MGX <- filtered_data_GMPS[, -c(1:16)]
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
  geom_point(aes(color = aza6mp), alpha = 0.5) +  # 将点的颜色映射到aza6mp列，根据你的数据情境你可以选择是否这样做
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # 添加线性回归线
  theme_minimal() +
  labs(title = "Scatter plot with Regression Line",
       x = "RPKM Value",
       y = "MBX Abundance",
       caption = paste("Pearson correlation: ", round(cor_result, 2)))

print(plot_pearson)
