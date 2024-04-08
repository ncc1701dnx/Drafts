which(MBX_HILpos_controls$Metabolite == "196.0609_2.81")#"196.0609_2.81", "154.0502_3.83"
which(MBX_HILpos_treatments$Metabolite == "196.0609_2.81")
control_mm <- MBX_HILpos_controls[which(MBX_HILpos_controls$Metabolite == "196.0609_2.81"), -c(1:7)]
treatment_mm <- MBX_HILpos_treatments[which(MBX_HILpos_treatments$Metabolite == "196.0609_2.81"), -c(1:7)]
control_mm <- as.numeric(control_mm)
treatment_mm <- as.numeric(treatment_mm)
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
t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
p_val <- t_test_result$p.value
fold_change <- round(median(treatment_mm)/median(control_mm), digits = 1)


####################### 5ASA #############################################################################
table(meta_data$oral5asa_loose)
table(meta_data$any5asa)
controls_ID <- meta_data$External.ID[which(meta_data$any5asa== 0)]
treatments_ID <- meta_data$External.ID[which(meta_data$any5asa == 1)]
length(controls_ID) ## 676
length(treatments_ID) ## 174
all_ID <- c(controls_ID, treatments_ID)
length(all_ID)

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

################## 5ASA new users #######################################################################
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

############################ 6MP #######################################################
controls_ID <- meta_data$External.ID[which(meta_data$aza6mp == 0)]
treatments_ID <- meta_data$External.ID[which(meta_data$aza6mp == 1)]
length(controls_ID) ## 772
length(treatments_ID) ## 78
#all_ID <- c(controls_ID, treatments_ID)
length(unique(meta_data$Participant.ID[which(meta_data$aza6mp == 0)])) ##90
length(unique(meta_data$Participant.ID[which(meta_data$aza6mp == 1)])) ##13

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


