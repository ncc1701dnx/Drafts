### 最开始，先把需要的包装上，以防R崩溃

# Load the needed packages
# mzR and xcms and edgeR are indispensable
# Write date 04 July ~ 05 July
library(stats) # In adjustRtime
library(Spectra)
library(mzR)
library(ggplot2)
library(RaMS)
#library(readMzXmlData)#use to quickly read mzXML files if not use xcms
library(xcms) # LC-MS data pre-processing 
#library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(MsFeatures)
library(rlang) # dependency for dplyr
library(dplyr)
#library(usethis) # dependency for devtool ignore after xcmsViewer installed
#library(devtools) # use for install xcmsViewer
#library(omicsViewer) # mandatory for xcmsViewer
#library(xcmsViewer) # MS data annotation
#library(parallel) # required for annoMS1 function from xcmsViewer
#library(BiocParallel) # required for annoMS1 function from xcmsViewer
library(MAIT) # dependency for xcmsViewer
library(MSnbase) # In most case already active in baseR, but just for sure 
library(stringr) # use for subset of stings

###########################################################################################################
#### 读取MBX结果文件
###########################################################################################################

# 读取MBX结果文件
MBX_results <- read.csv("inputs/iHMP_metabolomics.csv", header = T)
MBX_results <- MBX_results %>%
  mutate(Metabolite = if_else(Metabolite == "", 
                              paste0(m.z, "_", RT), 
                              Metabolite))
names(MBX_results) # in total 553

## 将MBX结果分成两个组别
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

##########################################################################################################
#### 现在是将函数load进来
##########################################################################################################


## 记得跑MBX_funtions.R中的函数，并且将以下两个list跑进去
adducts_pos <- list(
  "M + H" = 1.0078,
  "M + NH4" = 18.0344,
  "M + Na" = 22.9898,
  "M + K" = 38.9637,
  "M + H2O + H" = 19.0184,
  "M - H2O + H" = -17.0027,
  "M - 2H2O + H" = -35.0133,
  "M - H2O + NH4" = 0.0238,
  "M + H2O + Na" = 41.0003,
  "M + CH3OH + H" = 33.0340,
  "M + CH3OH + Na" = 55.0160,
  "M + CH3OH + Na + H2O" = 73.0265,
  "M + CH3OH + K" = 70.9899,
  "M + CH3CN + H" = 42.0344,
  "M + CH3CN + Na" = 64.0163,
  "M + H2O + CH3OH + H" = 51.0446,
  "M + Ni" = 57.9353,
  "M + Mo" = 97.9054,
  "M + Fe" = 55.9349,
  "M + Cu" = 62.9296,
  "M + CO2" = 43.9898,
  "M + CO" = 27.9949,
  "M + SO2" = 63.9619,
  "M + 2CH3CN + H" = 83.0609,
  "M + 2CH3CN + Na" = 105.0429,
  "M + 2CH3CN + Cu" = 144.9827,
  "M + 2CH3CN + Ni" = 139.9884,
  "M + 3CH3CN + Ni" = 181.0150,
  "M + 3CH3HCOO + Fe" = 235.9983,
  "M + 3CH3HCOO + Ni" = 237.9987,
  "M + K + CH3OH" = 70.9899
)

adducts_neg <- list(
  "M - H" = -1.0078,
  "M - 2H" = -2.0156,
  "M + H2O - H" = 17.0027,
  "M + CH3OH - H" = 31.0184,
  "M - 3H" = -3.0234,
  "M + CH3CN - H" = 40.0187,
  "M + Cl" = 34.9689,
  "M + Br" = 78.9183,
  "M + HCOO" = 44.9977,
  "M + CH3COO" = 59.0133,
  "M + HCOOH + HCOO" = 91.0031,
  "M + CH3COOH + HCOO" = 105.0188,
  "M + HCOOH + CH3COO" = 105.0188,
  "M + HSO4" = 96.9596,
  "M + H2PO4" = 96.9691,
  "M + CF3COO" = 112.9850
)

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

######################################################################################################################################
## 现在写四个函数，分别是1：自动做表，2：自动加significance，3：自动绘volcano plot，4：自动画significant的化合物boxplot
## IBDMDB_MBX_check这个函数的用处是对于给定的化合物，自动生成表格并且做t test
######################################################################################################################################
## 1：现在输入这个自动做表的函数：
IBDMDB_MBX_check <- function(formula, popular = T){
  ## 这个函数看起来很长， 但是实际上就是一个循环不同输入跑四次，别怕
  
  ################# 首先是用 ion charge为negative来跑两次 ####################
  
  ion_charge_mode = "neg"
  print("跑negative ion mode")
  subset_list_C18 <- dataset_list[c("MBX_C18neg_treatments", "MBX_C18neg_controls")]
  subset_list_HIL <- dataset_list[c("MBX_HILneg_treatments", "MBX_HILneg_controls")]
  
  ### 第一次跑之前需要设置这些为零
  ion_charge <- c()
  adduct <- c()
  m_z <- c()
  RT <- c()
  compound <- c()
  t_treatment <- c()
  t_control <- c()
  mean_treatment <- c()
  mean_control <- c()
  p_vals <- c()
  fold_changes <- c()
  
  ### 先跑C18的数据
  subset_list <- subset_list_C18
  print("跑C18的资料")
  
  subset_name <- names(subset_list)[1]
  subset <- subset_list[[subset_name]] ## 记得在这里和下面都要写上两个【【】】，因为是提取列表
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode, popular = popular)
  subset_treatment <- subset_list[[names(subset_list)[1]]] ## 如上 [[]]
  subset_controls <- subset_list[[names(subset_list)[2]]]  ## 如上 [[]]
  
  for (x in 1:nrow(com_range)) {
    print(paste0("x=",x))
    #x =1
    tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
    if(length(tt) == 0){next}
    
    for(i in 1:length(tt)){
      print(paste0("i=",i))
      #i = 3
      tt_t <- length(which(!is.na(subset_treatment[tt[i], -(1:7)])))
      tt_c <- length(which(!is.na(subset_controls[tt[i], -(1:7)])))
      compound <- c(compound, formula)
      adduct <- c(adduct, names(adducts_neg[x]))
      ion_charge <- c(ion_charge, ifelse(length(subset$Method[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$Method[tt[i]]))
      m_z <- c(m_z, ifelse(length(subset$m.z[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$m.z[tt[i]]))
      RT <- c(RT, ifelse(length(subset$RT[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$RT[tt[i]]))
      ## 以下是计算t test以及给出p value和fold change的两个if循环
      control_mm <- as.numeric(subset_controls[tt[i],-(1:7)])
      control_mm_2 <- control_mm
      treatment_mm <- as.numeric(subset_treatment[tt[i],-(1:7)])
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
      if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {
        t_test_result <- NA
        p_val <- NA
        fold_change <- NA}
      else {
        # t-test.if.all.patients
        t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
        p_val <- t_test_result$p.value
        fold_change <- t_test_result$estimate[2]/t_test_result$estimate[1]
        # wilcox.if.new.users
        #t_test_result <- wilcox.test(control_mm, treatment_mm)
        #p_val <- t_test_result$p.value
        #fold_change <- mean(exp(treatment_mm))/mean(exp(control_mm))
      }
      p_vals <- c(p_vals, p_val)
      fold_changes <- c(fold_changes, fold_change)
      ## t test计算完成
      mean_treatment <- c(mean_treatment, mean(as.numeric(subset_treatment[tt[i],-(1:7)]), na.rm = T))
      mean_control <- c(mean_control, mean(as.numeric(subset_controls[tt[i],-(1:7)]), na.rm = T))
      t_treatment <- c(t_treatment, tt_t)
      t_control <- c(t_control, tt_c)
    }
  }
  
  
  ### 然后跑HIL的数据
  subset_list <- subset_list_HIL
  print("跑HIL的数据")
  subset_name <- names(subset_list)[1]
  subset <- subset_list[[subset_name]] ## 记得在这里和下面都要写上两个【【】】，因为是提取列表
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode, popular = popular)
  subset_treatment <- subset_list[[names(subset_list)[1]]] ## 如上 [[]]
  subset_controls <- subset_list[[names(subset_list)[2]]]  ## 如上 [[]]
  
  for (x in 1:nrow(com_range)) {
    print(paste0("x=",x))
    tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
    if(length(tt) == 0){next}
    
    for(i in 1:length(tt)){
      print(paste0("i=",i))
      tt_t <- length(which(!is.na(subset_treatment[tt[i], -(1:7)])))
      tt_c <- length(which(!is.na(subset_controls[tt[i], -(1:7)])))
      compound <- c(compound, formula)
      adduct <- c(adduct, names(adducts_neg[x]))
      ion_charge <- c(ion_charge, ifelse(length(subset$Method[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$Method[tt[i]]))
      m_z <- c(m_z, ifelse(length(subset$m.z[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$m.z[tt[i]]))
      RT <- c(RT, ifelse(length(subset$RT[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$RT[tt[i]]))
      control_mm <- as.numeric(subset_controls[tt[i],-(1:7)])
      control_mm_2 <- control_mm
      treatment_mm <- as.numeric(subset_treatment[tt[i],-(1:7)])
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
      if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {
        t_test_result <- NA
        p_val <- NA
        fold_change <- NA}
      else {
        # t-test.if.all.patients
        t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
        p_val <- t_test_result$p.value
        fold_change <- t_test_result$estimate[2]/t_test_result$estimate[1]
        # wilcox.if.new.users
        #t_test_result <- wilcox.test(control_mm, treatment_mm)
        #p_val <- t_test_result$p.value
        #fold_change <- mean(exp(treatment_mm))/mean(exp(control_mm))
      }
      p_vals <- c(p_vals, p_val)
      fold_changes <- c(fold_changes, fold_change)
      mean_treatment <- c(mean_treatment, mean(as.numeric(subset_treatment[tt[i],-(1:7)]), na.rm = T))
      mean_control <- c(mean_control, mean(as.numeric(subset_controls[tt[i],-(1:7)]), na.rm = T))
      t_treatment <- c(t_treatment, tt_t)
      t_control <- c(t_control, tt_c)
    }
  }
  
  ################# 然后是用 ion charge为positive来跑两次 ####################
  ion_charge_mode = "pos"
  print("跑positive ion mode")
  subset_list_C18 <- dataset_list[c("MBX_C18pos_treatments", "MBX_C18pos_controls")]
  subset_list_HIL <- dataset_list[c("MBX_HILpos_treatments", "MBX_HILpos_controls")]
  
  ### 先跑C18的数据
  subset_list <- subset_list_C18
  print("跑C18的资料")
  subset_name <- names(subset_list)[1]
  subset <- subset_list[[subset_name]] ## 记得在这里和下面都要写上两个【【】】，因为是提取列表
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode, popular = popular)
  subset_treatment <- subset_list[[names(subset_list)[1]]] ## 如上 [[]]
  subset_controls <- subset_list[[names(subset_list)[2]]]  ## 如上 [[]]
  
  for (x in 1:nrow(com_range)) {
    x = 1
    print(paste0("x=",x))
    tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
    if(length(tt) == 0){next}
    
    for(i in 1:length(tt)){
      i =1
      print(paste0("i=",i))
      tt_t <- length(which(!is.na(subset_treatment[tt[i], -(1:7)])))
      tt_c <- length(which(!is.na(subset_controls[tt[i], -(1:7)])))
      compound <- c(compound, formula)
      adduct <- c(adduct, names(adducts_pos[x]))
      ion_charge <- c(ion_charge, ifelse(length(subset$Method[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$Method[tt[i]]))
      m_z <- c(m_z, ifelse(length(subset$m.z[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$m.z[tt[i]]))
      RT <- c(RT, ifelse(length(subset$RT[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$RT[tt[i]]))
      control_mm <- as.numeric(subset_controls[tt[i],-(1:7)])
      control_mm_2 <- control_mm
      treatment_mm <- as.numeric(subset_treatment[tt[i],-(1:7)])
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
      if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {
        t_test_result <- NA
        p_val <- NA
        fold_change <- NA}
      else {
        # t-test.if.all.patients
        t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
        p_val <- t_test_result$p.value
        fold_change <- t_test_result$estimate[2]/t_test_result$estimate[1]
        # wilcox.if.new.users
        #t_test_result <- wilcox.test(control_mm, treatment_mm)
        #p_val <- t_test_result$p.value
        #fold_change <- mean(exp(treatment_mm))/mean(exp(control_mm))
      }
      p_vals <- c(p_vals, p_val)
      fold_changes <- c(fold_changes, fold_change)
      mean_treatment <- c(mean_treatment, mean(as.numeric(subset_treatment[tt[i],-(1:7)]), na.rm = T))
      mean_control <- c(mean_control, mean(as.numeric(subset_controls[tt[i],-(1:7)]), na.rm = T))
      t_treatment <- c(t_treatment, tt_t)
      t_control <- c(t_control, tt_c)
    }
  }
  
  ### 然后跑HIL的数据
  subset_list <- subset_list_HIL
  print("跑HIL的数据")
  subset_name <- names(subset_list)[1]
  subset <- subset_list[[subset_name]] ## 记得在这里和下面都要写上两个【【】】，因为是提取列表
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode, popular = popular)
  subset_treatment <- subset_list[[names(subset_list)[1]]] ## 如上 [[]]
  subset_controls <- subset_list[[names(subset_list)[2]]]  ## 如上 [[]]
  
  for (x in 1:nrow(com_range)) {
    #x = 4
    print(paste0("x=",x))
    tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
    if(length(tt) == 0){next}
    
    for(i in 1:length(tt)){
      #i = 0
      print(paste0("i=",i))
      tt_t <- length(which(!is.na(subset_treatment[tt[i], -(1:7)])))
      tt_c <- length(which(!is.na(subset_controls[tt[i], -(1:7)])))
      compound <- c(compound, formula)
      adduct <- c(adduct, names(adducts_pos[x]))
      ion_charge <- c(ion_charge, ifelse(length(subset$Method[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$Method[tt[i]]))
      m_z <- c(m_z, ifelse(length(subset$m.z[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$m.z[tt[i]]))
      RT <- c(RT, ifelse(length(subset$RT[tt[i]]) == 0 | is.na(subset$Method[tt[i]]), NA, subset$RT[tt[i]]))
      control_mm <- as.numeric(subset_controls[tt[i],-(1:7)])
      control_mm_2 <- control_mm
      treatment_mm <- as.numeric(subset_treatment[tt[i],-(1:7)])
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
      if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {
        t_test_result <- NA
        p_val <- NA
        fold_change <- NA}
      else {
        # t-test.if.all.patients
        t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
        p_val <- t_test_result$p.value
        fold_change <- t_test_result$estimate[2]/t_test_result$estimate[1]
        # wilcox.if.new.users
        #t_test_result <- wilcox.test(control_mm, treatment_mm)
        #p_val <- t_test_result$p.value
        #fold_change <- mean(exp(treatment_mm))/mean(exp(control_mm))
      }
      p_vals <- c(p_vals, p_val)
      fold_changes <- c(fold_changes, fold_change)
      mean_treatment <- c(mean_treatment, mean(as.numeric(subset_treatment[tt[i],-(1:7)]), na.rm = T))
      mean_control <- c(mean_control, mean(as.numeric(subset_controls[tt[i],-(1:7)]), na.rm = T))
      t_treatment <- c(t_treatment, tt_t)
      t_control <- c(t_control, tt_c)
    }
  }
  
  compoundwise_IBDMDB <- data.frame(ion_charge = ion_charge, adduct = adduct, m_z = m_z, RTime = RT, formula = compound, 
                                    Num_treatment_21 = t_treatment, Num_control_100 = t_control, 
                                    mean_value_treatment = mean_treatment, mean_value_controls = mean_control, 
                                    P_Values = p_vals, Fold_changes = fold_changes)
  compoundwise_IBDMDB
}
###################################### 第一个函数完成 ###################################################################################

###################################### 写第二个和第三个函数 #############################################################################
## 为了画图好看，加上一列significance的变量
IBDMDB_MBX_add_signif <- function(IBDMDB_MBX_check_result){
  IBDMDB_MBX_check_result$signif <- ifelse(log2(IBDMDB_MBX_check_result$Fold_changes) > 1 & -log10(IBDMDB_MBX_check_result$P_Values) > 2, 
                                           "Significantly Abundant", "No Significance")
}

## 开始写自动化作图的函数
IBDMDB_MBX_plot <- function(IBDMDB_MBX_check_result, compound){
  ggplot(na.omit(IBDMDB_MBX_check_result), aes(x = log2(Fold_changes), y = -log10(P_Values), col = significance)) +
    geom_point(size = 3) +
    scale_colour_brewer(palette = "Set2")+
    xlab("Difference in means (Log2 Fold Change)") +
    ylab("-Log10 P-value") +
    ggtitle(paste0(compound, " Volcano Plot")) +
    theme(text=element_text(size=22),
          legend.position="right",
          strip.text = element_text(face = "bold"))
  # +theme_bw()
}

##############################再写一个函数，能够根据给出的结果表格找到原始的MBX_results表格的数据画box plot ##############################
#IBDMDB_6MP$m_z[which(IBDMDB_6MP$significance == "Significantly Abundant")]
IBDMDBresults_boxplot <- function(IBDMDB_compound, compound_name){
  rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
  rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  
  ###### Start of Plot Preparation ####
  for(i in 1:nrow(rr_result)){
    # 这里是后面画图的时候图片的名字
    #i = 1
    print(i)
    plot_title <- IBDMDB_compound[which(IBDMDB_compound$significance == "Significantly Abundant"),][1:4]
    pt <- paste(as.character(plot_title[i,]), collapse = " ")
    
    # 数据转换，log transform， Na eliminate，一如之前的函数
    control_mm <- as.numeric(rr_control[i,-(1:9)])
    control_mm_2 <- control_mm
    treatment_mm <- as.numeric(rr_treatment[i,-(1:9)])
    if(any(is.na(control_mm)) == T & all(is.na(control_mm)) == F ){
      control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
    }
    if(all(is.na(control_mm)) & all(is.na(treatment_mm)) == F){
      control_mm = 0.5*min(treatment_mm, na.rm = T)
      treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)
    }
    if(any(is.na(treatment_mm)) == T){
      treatment_mm[which(is.na(treatment_mm))] = 0.5*min(control_mm_2, na.rm = T) 
    }
    control_mm <- log(control_mm)
    treatment_mm <- log(treatment_mm)
    # 重新计算p value
    # if.t.test
    t_result <- t.test(treatment_mm, control_mm)
    p_value <- t_result$p.value
    # if.wilcox.test
    #t_result <- wilcox.test(control_mm, treatment_mm)
    #p_value <- t_result$p.value
    
    # 建立ggplot2主表格
    combined_data <- data.frame(
      abundance = c(treatment_mm, control_mm),
      group = c(rep("Treatment", length(treatment_mm)), rep("Control", length(control_mm)))
    )
    
    # ggplot2 画画图 
    # 根据p值确定标签的显示方式
    label_p <- ifelse(p_value < 0.001, "< 0.001", formatC(p_value, format = "e", digits = 2))
    
    # 手动添加色彩
    #my_colors <- brewer.pal(11, "Spectral")
    #my_colors <- my_colors[c( 9, 10, 11, 7)] # 手动定义颜色
    # IBDMDB 6MP
    my_colors <- c("#00798c", "#66a182")
    # IBDMDB 5ASA
    #my_colors <- c("#C4961A", "steelblue")
    
    # 使用ggplot2绘图
    p <- ggplot(combined_data, aes(x=group, y=abundance, color=group)) + 
      geom_boxplot(outlier.shape = NA, lwd = 1.5) + 
      #geom_jitter(position=position_jitter(0.2), size=4) +
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
}


############ 所有的函数写作完成 #########################################################################################################

#########################################################################################################
### 准备函数需要的数据
#### 与 work_0727_IBDMDB_recheck.R 一起使用
#########################################################################################################
IBDMDB_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/hmp2_metadata.csv", header = T)
### metadata的转换，如何变成IDs，请使用IBDMDB_metadata_reorganize.R和work_0727_IBDMDB_recheck.R

MBX_results <- read.csv("inputs/iHMP_metabolomics.csv", header = T)

MBX_results <- MBX_results %>%
  mutate(Metabolite = if_else(Metabolite == "", 
                              paste0(m.z, "_", RT), 
                              Metabolite))

controls_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% controls_ID))]
treatments_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% treatments_ID))]

##################### 现在去work_0727_IBDMDB_recheck.R用这四个函数计算全部感兴趣的化合物的状态 ##########################################


