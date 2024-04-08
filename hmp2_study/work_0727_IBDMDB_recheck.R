####################### 5ASA
table(meta_data$oral5asa_loose)
table(meta_data$any5asa)
controls_ID <- meta_data$External.ID[which(meta_data$any5asa== 0)]
treatments_ID <- meta_data$External.ID[which(meta_data$any5asa == 1)]
length(controls_ID) ## 676
length(treatments_ID) ## 174
all_ID <- c(controls_ID, treatments_ID)
length(all_ID)

####################### AZA/6MP
table(meta_data$aza6mp)
controls_ID <- meta_data$External.ID[which(meta_data$aza6mp == 0)]
treatments_ID <- meta_data$External.ID[which(meta_data$aza6mp == 1)]
length(controls_ID) ## 772
length(treatments_ID) ## 78
all_ID <- c(controls_ID, treatments_ID)
length(unique(meta_data$Participant.ID[which(meta_data$aza6mp == 0)])) ##90
length(unique(meta_data$Participant.ID[which(meta_data$aza6mp == 1)])) ##13


###########################
table(meta_data$diagnosis)
controls_ID <- meta_data$External.ID[which(meta_data$diagnosis == "nonIBD")]
treatments_ID <- meta_data$External.ID[which(meta_data$diagnosis == "CD" | meta_data$diagnosis == "UC")]


controls_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% controls_ID))]
treatments_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% treatments_ID))]

dim(MBX_C18pos_treatments)
which(MBX_HILpos_controls$Metabolite == "154.0502_3.83")
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

head(MBX_C18pos_controls)

##################### 现在开始用这四个函数计算全部感兴趣的化合物的状态 ##################################################################


################# Firstly the AZA/6MP compounds
## 使用6-MeMP
IBDMDB_6MeMP <- IBDMDB_MBX_check(formula = "C6H6N4S1", popular = F)
IBDMDB_6MeMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6MeMP)
IBDMDB_MBX_plot(IBDMDB_6MeMP, "6MeMP")
IBDMDBresults_boxplot(IBDMDB_6MeMP, "6MeMP")
# 6TUA
IBDMDB_6TUA <- IBDMDB_MBX_check(formula = "C5H4N4O2S1", popular = F)
IBDMDB_6TUA$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TUA)
IBDMDB_MBX_plot(IBDMDB_6TUA, "6TUA")
IBDMDBresults_boxplot(IBDMDB_6TUA, "6TUA")
# 6-MeTIMP
IBDMDB_6MeTIMP <- IBDMDB_MBX_check(formula = "C11H15N4O7P1S1", popular = F)
IBDMDB_6MeTIMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6MeTIMP)
IBDMDB_MBX_plot(IBDMDB_6MeTIMP, "6MeTIMP")
IBDMDBresults_boxplot(IBDMDB_6MeTIMP, "6MeTIMP")

# 6-TIMP
IBDMDB_6TIMP <- IBDMDB_MBX_check(formula = "C10H13N4O7P1S1", popular = F)
IBDMDB_6TIMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TIMP)
IBDMDB_MBX_plot(IBDMDB_6TIMP, "6TIMP")
IBDMDBresults_boxplot(IBDMDB_6TIMP, "6TIMP")
# AZA
IBDMDB_AZA <- IBDMDB_MBX_check(formula = "C9H7N7O2S1", popular = F)
IBDMDB_AZA$significance <- IBDMDB_MBX_add_signif(IBDMDB_AZA)
IBDMDB_MBX_plot(IBDMDB_AZA, "AZA")
IBDMDBresults_boxplot(IBDMDB_AZA, "AZA")
# 6-MP
IBDMDB_6MP <- IBDMDB_MBX_check(formula = "C5H4N4S1", popular = F)
IBDMDB_6MP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6MP)
IBDMDB_MBX_plot(IBDMDB_6MP, "6MP")
IBDMDBresults_boxplot(IBDMDB_6MP,"6MP")
# 6-TGMP
IBDMDB_6TGMP <- IBDMDB_MBX_check(formula = "C10H14N5O7P1S1", popular = F)
IBDMDB_6TGMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TGMP)
IBDMDB_MBX_plot(IBDMDB_6TGMP, "6TGMP")
IBDMDBresults_boxplot(IBDMDB_6TGMP, "6TGMP")
# 6-TGDP
IBDMDB_6TGDP <- IBDMDB_MBX_check(formula = "C10H15N5O10P2S1", popular = F)
IBDMDB_6TGDP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TGDP)
IBDMDB_MBX_plot(IBDMDB_6TGDP, "6TGDP")
IBDMDBresults_boxplot(IBDMDB_6TGDP, "6TGDP")
# 6-TGTP
IBDMDB_6_TGTP <- IBDMDB_MBX_check(formula = "C10H16N5O13P3S1", popular = F)
IBDMDB_6_TGTP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6_TGTP)
IBDMDB_MBX_plot(IBDMDB_6_TGTP, "6TGTP")
IBDMDBresults_boxplot(IBDMDB_6_TGTP, "6TGTP")



###################### Then 5-ASA Compounds ########################
### Be careful! Use IBDMDB_5ASA_check.R to re construct the datasets firstly

# 5-ASA
IBDMDB_5_ASA <- IBDMDB_MBX_check(formula = "C7H7N1O3", popular = T)
IBDMDB_5_ASA$significance <- IBDMDB_MBX_add_signif(IBDMDB_5_ASA)
IBDMDB_MBX_plot(IBDMDB_5_ASA, "5-ASA")
IBDMDBresults_boxplot(IBDMDB_5_ASA, "5ASA")
# N acetyl 5-ASA
IBDMDB_NAcetyl_5_ASA <- IBDMDB_MBX_check(formula = "C9H9N1O4")
IBDMDB_NAcetyl_5_ASA$significance <- IBDMDB_MBX_add_signif(IBDMDB_NAcetyl_5_ASA)
IBDMDB_MBX_plot(IBDMDB_NAcetyl_5_ASA, "N acetyl 5-ASA")
IBDMDBresults_boxplot(IBDMDB_NAcetyl_5_ASA, "N acetyl 5ASA")
# N butyryl 5-ASA

IBDMDB_NAcetyl_5_ASA <- IBDMDB_MBX_check(formula = "C11H13N1O4")

# N-Propionyl Mesalazine
IBDMDB_NAcetyl_5_ASA <- IBDMDB_MBX_check(formula = "C10H11N1O4")
IBDMDB_MBX_plot(IBDMDB_NAcetyl_5_ASA, "N Propionyl 5-ASA")


# nicotinuric acid
IBDMBX_NICA <- IBDMDB_MBX_check(formula = "C8H8N2O3")
IBDMBX_NICA$significance <- IBDMDB_MBX_add_signif(IBDMBX_NICA)
IBDMDB_MBX_plot(IBDMBX_NICA, "nicotinuric acid")
IBDMDBresults_boxplot(IBDMBX_NICA, "NICA")

#################### Then Melanie's compounds ###################

### 因为不确定是比较哪一个比较多，因此在这里也修改significance的判定
## 为了画图好看，加上一列significance的变量
IBDMDB_MBX_add_signif <- function(IBDMDB_MBX_check_result){
  IBDMDB_MBX_check_result$signif <- ifelse(abs(log2(IBDMDB_MBX_check_result$Fold_changes)) > 1 & -log10(IBDMDB_MBX_check_result$P_Values) > 2, 
                                           "Significantly Abundant", "No Significance")
}
## 因为melannie的这些图需要对比的不是treatment和control的区别，而是IBD和UC，UC和CD，IBD和CD的区别
## 需要修改下画图的组名和x轴的text，在这块修改画图函数
IBDMDB_compound <- IBDMDB_6MP
### Box plot函数
IBDMDBresults_boxplot <- function(IBDMDB_compound, compound_name){
  rr_result <- IBDMDB_compound[,c(1:5)][which(IBDMDB_compound$significance == "Significantly Abundant"),]
  rr_treatment <- dplyr::inner_join(rr_result, treatments_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  rr_control <- dplyr::inner_join(rr_result, controls_MBX_results, by = c("ion_charge" = "Method", "m_z" = "m.z", "RTime" = "RT"))
  
  ###### Start of Plot Preparation ####
  for(i in 1:nrow(rr_result)){
    # 这里是后面画图的时候图片的名字
    print(i)
    plot_title <- IBDMDB_compound[which(IBDMDB_compound$significance == "Significantly Abundant"),][1:4]
    pt <- paste(as.character(plot_title[i,]), collapse = " ")
    
    # 数据转换，log transform， Na eliminate，一如之前的函数
    control_mm <- as.numeric(rr_control[i,-(1:9)])
    treatment_mm <- as.numeric(rr_treatment[i,-(1:9)])
    if(any(is.na(treatment_mm)) == T & all(is.na(treatment_mm)) == F){
      treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T) 
    }
    if(any(is.na(control_mm)) == T & all(is.na(control_mm)) == F ){
      control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
    }
    if(all(is.na(control_mm)) & (length(which(!is.na(treatment_mm))) >= 4)){
      control_mm = 0.5*min(treatment_mm, na.rm = T)
      treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)}
    if(all(is.na(treatment_mm)) & (length(which(!is.na(control_mm))) >= 4)){
      treatment_mm = 0.5*min(control_mm, na.rm = T)
    }
    control_mm <- log(control_mm)
    treatment_mm <- log(treatment_mm)
    # 重新计算p value
    t_result <- t.test(treatment_mm, control_mm)
    p_value <- t_result$p.value
    
    # 建立ggplot2主表格
    combined_data <- data.frame(
      abundance = c(treatment_mm, control_mm),
      group = c(rep("CD", length(treatment_mm)), rep("UC", length(control_mm)))
    )
    
    # ggplot2 画画图 
    # 根据p值确定标签的显示方式
    label_p <- ifelse(p_value < 0.001, "< 0.001", formatC(p_value, format = "e", digits = 2))
    
    # 使用ggplot2绘图
    p <- ggplot(combined_data, aes(x=group, y=abundance, color=group)) + 
      geom_boxplot(outlier.shape = NA, lwd = 1.5) + 
      geom_jitter(size=4, width = 0.2, height = 0) +
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

### Firstly the comparison between IBD and Non-IBD
table(meta_data$diagnosis)
controls_ID <- meta_data$External.ID[which(meta_data$diagnosis == "nonIBD")]
treatments_ID <- meta_data$External.ID[which(meta_data$diagnosis == "CD" | meta_data$diagnosis == "UC")]

controls_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% controls_ID))]
treatments_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% treatments_ID))]
## Dihydroaeruginoic acid

IBDMDB_dhia <- IBDMDB_MBX_check(formula = "C10H9N1O3S1", popular = F)
IBDMDB_dhia$significance <- IBDMDB_MBX_add_signif(IBDMDB_dhia)
IBDMDB_MBX_plot(IBDMDB_dhia, "Dihydroaeruginoic acid")
IBDMDBresults_boxplot(IBDMDB_dhia, "Dihydroaeruginoic acid")

## salicylic acid

IBDMDB_SA <- IBDMDB_MBX_check(formula = "C7H6O3", popular = F)
IBDMDB_SA$significance <- IBDMDB_MBX_add_signif(IBDMDB_SA)
IBDMDB_MBX_plot(IBDMDB_SA, "salicylic acid")
IBDMDBresults_boxplot(IBDMDB_SA, "salicylic acid")

## cysteine

IBDMDB_CYS <- IBDMDB_MBX_check(formula = "C3H7N1O2S1", popular = F)
IBDMDB_CYS$significance <- IBDMDB_MBX_add_signif(IBDMDB_CYS)
IBDMDB_MBX_plot(IBDMDB_CYS, "cysteine")
IBDMDBresults_boxplot(IBDMDB_CYS, "cysteine")


### Then the comparison between IBD and UC
table(meta_data$diagnosis)
controls_ID <- meta_data$External.ID[which(meta_data$diagnosis == "nonIBD")]
treatments_ID <- meta_data$External.ID[which(meta_data$diagnosis == "UC")]

controls_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% controls_ID))]
treatments_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% treatments_ID))]


### Then the comparison between IBD and CD
table(meta_data$diagnosis)
controls_ID <- meta_data$External.ID[which(meta_data$diagnosis == "nonIBD")]
treatments_ID <- meta_data$External.ID[which(meta_data$diagnosis == "CD")]

controls_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% controls_ID))]
treatments_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% treatments_ID))]

### Then the comparison between UC and CD
table(meta_data$diagnosis)
controls_ID <- meta_data$External.ID[which(meta_data$diagnosis == "UC")]
treatments_ID <- meta_data$External.ID[which(meta_data$diagnosis == "CD")]

controls_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% controls_ID))]
treatments_MBX_results <- MBX_results[ ,c(1:7, which(colnames(MBX_results) %in% treatments_ID))]

####  
## 检查下有哪些病人是新用户
#################################################
## 先转换下格式，马上再转换回去:)
meta_data$aza6mp <- as.numeric(meta_data$aza6mp)
meta_data$oral5asa_loose <- as.numeric(meta_data$oral5asa_loose)


### 先比较AZA和6-MP的
## loose means all patients have changed the drug usage will be considered
result_aza6mp_loose <- meta_data %>%
  # Group by Participant.ID
  group_by(Participant.ID) %>%
  
  # Calculate the number of changes in aza6mp for each Participant.ID
  mutate(change_in_aza6mp = sum(abs(diff(aza6mp)))) %>%
  
  # Ungroup to operate on the full data again
  ungroup() %>%
  
  # Filter rows for participants with changes in aza6mp
  filter(change_in_aza6mp > 0) %>%
  
  # Select required columns
  select(Participant.ID, visit_num, aza6mp, External.ID) %>%
  
  # Remove duplicate rows
  distinct()

## strict means only the new users of aza6mp will be considered

result_aza6mp_strict <- meta_data %>%
  group_by(Participant.ID) %>%
  mutate(prev_aza6mp = lag(aza6mp)) %>%
  filter(!is.na(prev_aza6mp) & prev_aza6mp == 0 & aza6mp == 1) %>%
  ungroup() %>%
  select(Participant.ID, visit_num, External.ID, aza6mp)

View(result_aza6mp_strict)
meta_data$aza6mp[which(meta_data$Participant.ID == result_aza6mp_strict$Participant.ID[4])]


### great! 现在比较New User的区别
IDs <- as.character(result_aza6mp_strict$Participant.ID)

treatments_ID <- meta_data %>% 
  filter(Participant.ID %in% IDs & aza6mp == 1)
treatments_ID <- treatments_ID$External.ID
controls_ID <- meta_data %>% 
  filter(Participant.ID %in% IDs & aza6mp == 0)
controls_ID <- controls_ID$External.ID

meta_data$any5asa
meta_data$aza6mp
### 下面是5-ASA的new user 比较

## loose means all patients have changed the drug usage will be considered
result_5asa_loose <- meta_data %>%
  # Group by Participant.ID
  group_by(Participant.ID) %>%
  
  # Calculate the number of changes in any5asa for each Participant.ID
  mutate(change_in_any5asa = sum(abs(diff(any5asa)))) %>%
  
  # Ungroup to operate on the full data again
  ungroup() %>%
  
  # Filter rows for participants with changes in any5asa
  filter(change_in_any5asa > 0) %>%
  
  # Select required columns
  select(Participant.ID, visit_num, any5asa, External.ID) %>%
  
  # Remove duplicate rows
  distinct()
## strict means only the new users of 5asa will be considered
View(result_5asa_loose)

result_5asa_strict <- meta_data %>%
  group_by(Participant.ID) %>%
  mutate(prev_any5asa = lag(any5asa)) %>%
  filter(!is.na(prev_any5asa) & prev_any5asa == 0 & any5asa == 1) %>%
  ungroup() %>%
  select(Participant.ID, visit_num, External.ID, any5asa)
View(result_5asa_strict)
result_5asa_strict$Participant.ID %in% result_5asa_loose$Participant.ID
### great! 现在比较New User的区别
IDs <- as.character(result_5asa_strict$Participant.ID)
length(which(meta_data$any5asa == 1))
length(treatments_ID)
treatments_ID <- meta_data %>% 
  filter(Participant.ID %in% IDs & any5asa == 1)
treatments_ID <- treatments_ID$External.ID
controls_ID <- meta_data %>% 
  filter(Participant.ID %in% IDs & any5asa == 0)
controls_ID <- controls_ID$External.ID

### MGX
MGXresult_aza6mp_strict <- meta_MGXdata %>%
  group_by(Participant.ID) %>%
  mutate(prev_aza6mp = lag(aza6mp)) %>%
  filter(!is.na(prev_aza6mp) & prev_aza6mp == 0 & aza6mp == 1) %>%
  ungroup() %>%
  select(Participant.ID, visit_num, External.ID, aza6mp)

View(MGXresult_aza6mp_strict)
meta_MGXdata$aza6mp[which(meta_MGXdata$Participant.ID == MGXresult_aza6mp_strict$Participant.ID[4])]


### great! 现在比较New User的区别
IDs <- as.character(MGXresult_aza6mp_strict$Participant.ID)

treatments_ID <- meta_MGXdata %>% 
  filter(Participant.ID %in% IDs & aza6mp == 1)
treatments_ID <- treatments_ID$External.ID
controls_ID <- meta_MGXdata %>% 
  filter(Participant.ID %in% IDs & aza6mp == 0)
controls_ID <- controls_ID$External.ID

class(treatments_ID)
writeLines(treatments_ID, "outputs/IBDMDB_MGX_aza6mp_new_user_treatmentID.txt")
writeLines(controls_ID, "outputs/IBDMDB_MGX_aza6mp_new_user_controlID.txt")
