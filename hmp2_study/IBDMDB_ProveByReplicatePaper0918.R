
library(ggplot2)
#library(readMzXmlData)#use to quickly read mzXML files if not use xcms
#library(faahKO)
library(RColorBrewer)
library(rlang) # dependency for dplyr
library(dplyr)
library(tidyr)
#library(usethis) # dependency for devtool ignore after xcmsViewer installed
#library(devtools) # use for install xcmsViewer
#library(omicsViewer) # mandatory for xcmsViewer
#library(xcmsViewer) # MS data annotation
#library(parallel) # required for annoMS1 function from xcmsViewer
#library(BiocParallel) # required for annoMS1 function from xcmsViewer
 # In most case already active in baseR, but just for sure 
#library(stringr) # use for subset of stings
GMPS_Treatment$group <- "Treatment"
GMPS_Control$group <- "Control"

combined_data <- rbind(GMPS_Treatment, GMPS_Control)
test_result <- t.test(rpkm_value ~ group, data = combined_data)
p <- ggplot(combined_data, aes(x = group, y = rpkm_value)) +
  geom_boxplot(aes(fill = group)) +
  labs(title = "Comparison of RPKM values between Treatment and Control", x = "Group", y = "RPKM value") +
  #theme_minimal() +
  geom_text(aes(x = 1.5, y = max(rpkm_value), label = paste0("p-value = ", round(test_result$p.value, 5))),
            vjust = 2) +
  theme(text=element_text(size=30), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold")) +
  theme_bw()

print(p)



### 5ASA
IBDMDB_compound <- IBDMDB_5_ASA
compound_name <- "5ASA"
rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_result <- dplyr::inner_join(rr_result, MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
i = 16
paper_result_5ASA <- rr_result[i, -c(1:9)]
paper_treatment_5ASA <- rr_result[i, -c(1:9)]
paper_treatment_5ASA <- paper_treatment_5ASA[which(paper_treatment_5ASA > 10000000)]
paper_treatment_5ASA <- colnames(paper_treatment_N5ASA)


### N5ASA

i = 5
IBDMDB_compound <- IBDMDB_NAcetyl_5_ASA
compound_name <- "N acetyl 5-ASA"
rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
rr_result <- dplyr::inner_join(rr_result, MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
paper_treatment_N5ASA <- rr_result[i, -c(1:9)]
paper_treatment_N5ASA <- paper_treatment_N5ASA[which(paper_treatment_N5ASA > 10000000)]
paper_treatment_N5ASA <- colnames(paper_treatment_N5ASA)
###### Plot group-wise ####
for(i in 1:nrow(rr_result)){
  # 这里是后面画图的时候图片的名字
  i = 13
  print(i)
  plot_title <- IBDMDB_compound[which(IBDMDB_compound$significance == "Significantly Abundant"),][1:4]
  pt <- paste(as.character(plot_title[i,]), collapse = " ")
  pt <- paste(pt, i, collapse = " ")
  #pt <- "X154.0502_3.83"
  # 数据转换，log transform， Na eliminate，一如之前的函数
  control_mm <- rr_control[i,-(1:9)]
  #control_mm <- as.numeric(rr_control)
  treatment_mm <- rr_treatment[i,-(1:9)]
  names_control_mm <- colnames(control_mm)
  names_treatment_mm <- colnames(treatment_mm)
  control_mm <- as.numeric(control_mm)
  treatment_mm <- as.numeric(treatment_mm)
  control_mm_2 <- control_mm
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
  if_paper_treatment_control <- c()
  if_paper_treatment_treatment <- c()
  rm(prop_c)
  rm(prop_t)
  if_paper_treatment_control = ifelse(names_control_mm %in% paper_treatment_N5ASA, "Yes", "No")
  prop_c = round(length(which(if_paper_treatment_control == "Yes"))/length(if_paper_treatment_control), digits = 2)
  if_paper_treatment_treatment = ifelse(names_treatment_mm %in% paper_treatment_N5ASA, "Yes", "No")
  prop_t = round(length(which(if_paper_treatment_treatment == "Yes"))/length(if_paper_treatment_treatment), digits = 2)
  if_paper_treatment = c(if_paper_treatment_control, if_paper_treatment_treatment)
  # 重新计算p value
  #t_result <- wilcox.test(control_mm, treatment_mm)
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
  my_colors <- c("#C5A06B", "#5E9CA8")
  
  colnames(combined_data)
  # 使用ggplot2绘图
  p <- ggplot(combined_data, aes(x=group, y=abundance)) + 
    geom_boxplot(outlier.shape = NA, lwd = 1.5, color="black") +  # Set a uniform color for boxplots
    geom_jitter(aes(color=if_paper_treat), size=2, width = 0.1, height = 0) +  # Color points based on if_paper_treat
    scale_color_manual(values = my_colors) +
    xlab("Groups") +
    ylab("Log e Abundance") +
    theme_minimal()
  
  # 计算两个组的数据的最大值，以确定标注线的位置
  y_max <- max(combined_data$abundance)
  
  # 在图上增加一个标注线和文本
  p <- p + 
    geom_segment(aes(x = 1, xend = 2, y = y_max + 0.05 * y_max, yend = y_max + 0.05 * y_max), color = "black") +
    geom_text(aes(x = 1.5, y = y_max + 0.1 * y_max, label = paste(prop_c, " ", "p =", p_value, " ", prop_t)), color = "black", size = 8) +
    ggtitle(paste0(compound_name, " ", pt))+
    theme(text=element_text(size=22),
          legend.position="right",
          strip.text = element_text(face = "bold"))
  plot(p)
}


###### Plot non-group with correlation ####
for(i in 1:nrow(rr_result)){
  # 这里是后面画图的时候图片的名字
  i =3
  print(i)
  plot_title <- IBDMDB_compound[which(IBDMDB_compound$significance == "Significantly Abundant"),][1:4]
  pt <- paste(as.character(plot_title[i,]), collapse = " ")
  pt <- paste(pt, i, collapse = " ")
  #pt <- "X154.0502_3.83"
  # 数据转换，log transform， Na eliminate，一如之前的函数
  control_mm <- rr_control[i,-(1:9)]
  #control_mm <- as.numeric(rr_control)
  treatment_mm <- rr_treatment[i,-(1:9)]
  names_control_mm <- colnames(control_mm)
  names_treatment_mm <- colnames(treatment_mm)
  control_mm <- as.numeric(control_mm)
  treatment_mm <- as.numeric(treatment_mm)
  control_mm_2 <- control_mm
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
  #control_mm <- log(control_mm)
  #treatment_mm <- log(treatment_mm)
  all_mm_names <- c(names_control_mm, names_treatment_mm)
  all_mm <- c(control_mm, treatment_mm)
  all_mm <- as.matrix(all_mm)
  all_mm <- t(all_mm)
  colnames(all_mm) <- all_mm_names
  all_mm <- as.data.frame(all_mm)
  
  all_mm_long <- all_mm %>% gather(key = "sample", value = "value_all_mm")
  paper_result_5ASA_long <- paper_result_5ASA %>% gather(key = "sample", value = "value_paper_treatment")
  
  combined_data <- left_join(all_mm_long, paper_result_5ASA_long, by = "sample")
  
  # 重新计算p value
  # if.t.test
  t_result <- t.test(treatment_mm, control_mm)
  p_value <- t_result$p.value
  
  # 计算Pearson相关性
  cor_value <- cor(combined_data$value_all_mm, combined_data$value_paper_treatment, method = "pearson")
  # 画散点图并拟合
  plot_cor <- combined_data %>% 
    ggplot(aes(x = value_all_mm, y = value_paper_treatment)) +
    geom_point() +
    ggtitle(paste0("Pearson:", " ", round(cor_value, 2), " ", i)) +
    #geom_smooth(method = "lm", se = FALSE, color = "blue") +
    theme(text=element_text(size=24), #调整字号
          legend.position="right",
          strip.text = element_text(face = "bold")) +
    theme_minimal() 
 
  # 计算Pearson相关性
  #cor_value <- cor(combined_data$value_all_mm, combined_data$value_paper_treatment, method = "pearson")
  
  # 将cor_value添加到图上
  #plot_cor <- plot_cor + geom_text(aes(x = 3, y = 0,label = paste("Pearson:", round(cor_value, 2))), color = "black", size = 6)
  print(plot_cor)
}

###### Plot non-group with correlation ####
for(i in 1:nrow(rr_result)){
  # 这里是后面画图的时候图片的名字
  #i =1
  print(i)
  plot_title <- IBDMDB_compound[which(IBDMDB_compound$significance == "Significantly Abundant"),][1:4]
  pt <- paste(as.character(plot_title[i,]), collapse = " ")
  pt <- paste(pt, i, collapse = " ")
  #pt <- "X154.0502_3.83"
  # 数据转换，log transform， Na eliminate，一如之前的函数
  control_mm <- rr_control[i,-(1:9)]
  #control_mm <- as.numeric(rr_control)
  treatment_mm <- rr_treatment[i,-(1:9)]
  names_control_mm <- colnames(control_mm)
  names_treatment_mm <- colnames(treatment_mm)
  control_mm <- as.numeric(control_mm)
  treatment_mm <- as.numeric(treatment_mm)
  control_mm_2 <- control_mm
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
  #control_mm <- log(control_mm)
  #treatment_mm <- log(treatment_mm)
  all_mm_names <- c(names_control_mm, names_treatment_mm)
  all_mm <- c(control_mm, treatment_mm)
  all_mm <- as.matrix(all_mm)
  all_mm <- t(all_mm)
  colnames(all_mm) <- all_mm_names
  all_mm <- as.data.frame(all_mm)
  
  all_mm_long <- all_mm %>% gather(key = "sample", value = "value_all_mm")
  paper_result_5ASA_long <- paper_result_5ASA %>% gather(key = "sample", value = "value_paper_treatment")
  
  combined_data <- left_join(all_mm_long, paper_result_5ASA_long, by = "sample")
  
  # 计算Pearson相关性
  cor_value <- cor(combined_data$value_all_mm, combined_data$value_paper_treatment, method = "pearson")
  # 画散点图并拟合
  plot_cor <- combined_data %>% 
    ggplot(aes(x = value_all_mm, y = value_paper_treatment)) +
    geom_point() +
    ggtitle(paste0("Pearson:", " ", round(cor_value, 2), " ", i)) +
    #geom_smooth(method = "lm", se = FALSE, color = "blue") +
    theme(text=element_text(size=24), #调整字号
          legend.position="right",
          strip.text = element_text(face = "bold")) +
    theme_minimal() 
  
  # 计算Pearson相关性
  #cor_value <- cor(combined_data$value_all_mm, combined_data$value_paper_treatment, method = "pearson")
  
  # 将cor_value添加到图上
  #plot_cor <- plot_cor + geom_text(aes(x = 3, y = 0,label = paste("Pearson:", round(cor_value, 2))), color = "black", size = 6)
  print(plot_cor)
}


###### Give the table of possible 5-ASA results ####
PeakNum <- c()
PeakMZ <- c()
PeakRT <- c()
PeakIonChargeMode <- c()
PeakAdduct <- c()
pval <- c()
corval <- c()
FoldChange <- c()
for(i in 1:nrow(rr_result)){
  # 这里是后面画图的时候图片的名字
  #i =16
  print(i)
  PeakNum <- c(PeakNum, i)
  plot_title <- IBDMDB_compound[which(IBDMDB_compound$significance == "Significantly Abundant"),][1:4]
  plot_title <- plot_title[i, ]
  #pt <- "X154.0502_3.83"
  # 数据转换，log transform， Na eliminate，一如之前的函数
  control_mm <- rr_control[i,-(1:9)]
  #control_mm <- as.numeric(rr_control)
  treatment_mm <- rr_treatment[i,-(1:9)]
  names_control_mm <- colnames(control_mm)
  names_treatment_mm <- colnames(treatment_mm)
  control_mm <- as.numeric(control_mm)
  treatment_mm <- as.numeric(treatment_mm)
  control_mm_2 <- control_mm
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
  #control_mm <- log(control_mm)
  #treatment_mm <- log(treatment_mm)
  all_mm_names <- c(names_control_mm, names_treatment_mm)
  all_mm <- c(control_mm, treatment_mm)
  all_mm <- as.matrix(all_mm)
  all_mm <- t(all_mm)
  colnames(all_mm) <- all_mm_names
  all_mm <- as.data.frame(all_mm)
  
  all_mm_long <- all_mm %>% gather(key = "sample", value = "value_all_mm")
  paper_result_5ASA_long <- paper_result_5ASA %>% gather(key = "sample", value = "value_paper_treatment")
  
  combined_data <- left_join(all_mm_long, paper_result_5ASA_long, by = "sample")
  
  # 输出M/Z等到数个string中去
  PeakMZ <- c(PeakMZ, plot_title$m_z)
  PeakRT <- c(PeakRT, plot_title$RTime)
  PeakIonChargeMode <- c(PeakIonChargeMode, plot_title$ion_charge)
  PeakAdduct <- c(PeakAdduct, plot_title$adduct)
  
  # 重新计算p value & Fold Change
  # if.t.test
  t_result <- t.test(treatment_mm, control_mm)
  p_value <- t_result$p.value
  pval <- c(pval, p_value)
  fold_change <- round(median(treatment_mm)/median(control_mm), digits = 1)
  FoldChange <- c(FoldChange, fold_change)
  # 计算Pearson相关性
  cor_value <- cor(combined_data$value_all_mm, combined_data$value_paper_treatment, method = "pearson")
  corval <- c(corval, cor_value)
}

table_5ASA <- data.frame("Number" = PeakNum, "M/Z" = PeakMZ, "RTime" = PeakRT, "Ion Charge Mode" = PeakIonChargeMode, 
                         "Adduct" = PeakAdduct, "P Value" = pval, "Fold Change" = FoldChange, "Correlation" = corval)

table_5ASA <- table_5ASA %>% 
  arrange(P.Value)
#table_5ASA$Number <- c(1:nrow(table_5ASA))
table_5ASA$Number[1]
print(table_5ASA)  # 显示前几行查看排序结果
table_5ASA[, -1]
write.csv(table_5ASA, file = "outputs/table5ASA.csv")
