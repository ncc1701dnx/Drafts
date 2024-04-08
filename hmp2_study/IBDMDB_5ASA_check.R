## work 7th July 2023
## This is for check 5ASA in IBDMDB dataset
## Derived from IBDMDB_functions_0709.R
## 请遵照上述script的要求进行包loading和数据处理
IBDMDB_metadat
names(IBDMDB_metadata)
# mesalazine is called mesalamine in the USA
table(IBDMDB_metadata$Apriso..mesalamine.) # a brand of mesalazine
# no current, omit this.
table(IBDMDB_metadata$Asacol..mesalamine.) # a brand of delayed release mesalazine
# current 26, never taken 468. 
table(IBDMDB_metadata$Lialda..mesalamine.) # a brand of delayed release mesalazine 
# current 36, never taken 466.
table(IBDMDB_metadata$Rowasa.enemas..mesalamine.enemas.) # a brand of mesalazine used from rectum
# current 8, never taken 491.
table(IBDMDB_metadata$Canasa.suppositories..mesalamine.suppositories.) # a brand of delayed release mesalazine used from rectum
# current 12, never taken 481.
## reconstruct the MBX_metadata table to make it have information for 5-ASA usage
ASA_MBX_IBDMDB_metadata <- IBDMDB_metadata[which(IBDMDB_metadata$data_type == "metabolomics"),]

########### 以下是重构IBDMDB的metadata的table的环节 ##################################################################################
# 原metadata table实在是太差了，MBX数据里面没有任何的病人用药信息
# 用唯一matchable的ID：Participant.ID来计算总和的病人用药情况，然后插入MBX的数据里
########################################################################################################################################
# 如同AZA和6MP。不过我们这里将所有5ASA病人合并到一起
# check who are using 5ASA
length(which(IBDMDB_metadata$Asacol..mesalamine. == "Current" | 
               IBDMDB_metadata$Lialda..mesalamine. == "Current" |
                IBDMDB_metadata$Rowasa.enemas..mesalamine.enemas. == "Current" |
                 IBDMDB_metadata$Canasa.suppositories..mesalamine.suppositories. == "Current"))
mesalazine_treat <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Asacol..mesalamine. == "Current" | 
                                                           IBDMDB_metadata$Lialda..mesalamine. == "Current" |
                                                           IBDMDB_metadata$Rowasa.enemas..mesalamine.enemas. == "Current" |
                                                           IBDMDB_metadata$Canasa.suppositories..mesalamine.suppositories. == "Current")]
which(ASA_MBX_IBDMDB_metadata$Participant.ID %in% mesalazine_treat)

# check who never used 5ASA
length(which(IBDMDB_metadata$Apriso..mesalamine. == "Never taken" &
               IBDMDB_metadata$Asacol..mesalamine. == "Never taken" &
               IBDMDB_metadata$Lialda..mesalamine. == "Never taken" &
               IBDMDB_metadata$Rowasa.enemas..mesalamine.enemas. == "Never taken" &
               IBDMDB_metadata$Canasa.suppositories..mesalamine.suppositories. == "Never taken"))
mesalazine_control <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Apriso..mesalamine. == "Never taken" &
                                                             IBDMDB_metadata$Asacol..mesalamine. == "Never taken" &
                                                             IBDMDB_metadata$Lialda..mesalamine. == "Never taken" &
                                                             IBDMDB_metadata$Rowasa.enemas..mesalamine.enemas. == "Never taken" &
                                                            IBDMDB_metadata$Canasa.suppositories..mesalamine.suppositories. == "Never taken")]
which(ASA_MBX_IBDMDB_metadata$Participant.ID %in% mesalazine_control)

### 在这里我们发现有些participantID在两个组别中都出现了，我们将这些conflict的ID去除以免干扰
length(unique(mesalazine_treat[which(mesalazine_treat %in% mesalazine_control)]))
mesalazine_conflick <- mesalazine_treat[which(mesalazine_treat %in% mesalazine_control)]
mesalazine_treat <- mesalazine_treat[which(!mesalazine_treat %in% mesalazine_conflick)]
mesalazine_control <- mesalazine_control[which(!mesalazine_control %in% mesalazine_conflick)]
## 这里我们创造一个新的ASA metadata的列，以免麻烦
length(unique(mesalazine_control))
length(unique(mesalazine_treat))
dim(ASA_MBX_IBDMDB_metadata)
ASA_MBX_IBDMDB_metadata$mesalamine_total <- rep(NA, nrow(ASA_MBX_IBDMDB_metadata))
## 现在将两个组别加到新建的列里
ASA_MBX_IBDMDB_metadata$mesalamine_total[which(ASA_MBX_IBDMDB_metadata$Participant.ID %in% mesalazine_treat)] = "Current"
ASA_MBX_IBDMDB_metadata$mesalamine_total[which(ASA_MBX_IBDMDB_metadata$Participant.ID %in% mesalazine_control)] = "Never taken"
## 不同于AZA/6MP，5ASA针对于整个IBD patients都有作用，因此不作细分

# 师兄建议，可以先检查下数据是否满足正态分布
# 先安装并加载需要的包
# if (!require("tidyverse")) install.packages("tidyverse")
# library(tidyverse)
# length(unique(MBX_results_long$PatientID))
# class(MBX_results_long$Abundance)
# rm(MBX_results_long)
# # 转置数据框，并且将数据转换为long format以方便做箱型图
# MBX_results_long <- MBX_results[, -(1:7)] %>%
#   t() %>%   # 转置数据
#   as.data.frame() %>%  # 将结果转化为数据框
#   rownames_to_column("PatientID") %>%  # 将行名变为一列
#   pivot_longer(-PatientID, names_to = "CompoundID", values_to = "Abundance")  # 转换为长格式
# 
# # 现在，我们可以使用ggplot2来创建一个箱线图
# ggplot(MBX_results_long, aes(x = PatientID, y = Abundance)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
#   labs(x = "Patient ID", y = "Compound Abundance", title = "Compound Abundance per Patient")
tt <- ASA_control_group[5,]
tt[1,] <- rep(10, 330)

## 在这里先把
table(ASA_MBX_IBDMDB_metadata$mesalamine_total)
ASA_MBX_IBDMDB_metadata$External.ID[which(ASA_MBX_IBDMDB_metadata$mesalamine_total == "Never taken")]
ASA_controls_ID <- ASA_MBX_IBDMDB_metadata$External.ID[which(ASA_MBX_IBDMDB_metadata$mesalamine_total == "Never taken")]
ASA_treatment_ID <- ASA_MBX_IBDMDB_metadata$External.ID[which(ASA_MBX_IBDMDB_metadata$mesalamine_total == "Current")]
ASA_control_group <- MBX_results[, ASA_controls_ID]
ASA_treatment_group <- MBX_results[, ASA_treatment_ID]
ASA_treatment_group[rownames(MBX_results)[1],]
length(ASA_controls_ID)
length(ASA_treatment_ID)
which(is.na(ASA_control_group))
### 这里新写一个函数，将NA value的赋值和log transform 全部先做好

if(any(is.na(as.numeric(ASA_treatment_group[4,])))){print("yes")}
rownames(MBX_results)[1]
ASA_treatment_group["HILp_QI21194",]
for(i in 1: nrow(MBX_results)){
  metabolite <- rownames(MBX_results)[i]
  print(i)
  treatment_mm <- ASA_treatment_group[metabolite,]
  
  treatment_mm <- as.numeric(treatment_mm)
  
  if(any(is.na(treatment_mm)) == T & all(is.na(treatment_mm)) == F){
    treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T) 
  }
  
  ASA_treatment_group[metabolite,] <- treatment_mm
}

# ASA_treatment_group <- lapply(rownames(MBX_results), function(metabolite){
#   
#   treatment_mm <- ASA_treatment_group[metabolite,]
#   
#   treatment_mm <- as.numeric(treatment_mm)
#   
#   if(any(is.na(treatment_mm)) & all(!is.na(treatment_mm))){
#     treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T) 
#   }
#   
#   ASA_treatment_group[metabolite,] <- treatment_mm
#   ASA_treatment_group
# })
for(i in 1: nrow(MBX_results)){
  metabolite <- rownames(MBX_results)[i]
  print(i)
  control_mm <- ASA_control_group[metabolite,]
  
  control_mm <- as.numeric(control_mm)
  
  if(any(is.na(control_mm)) == T & all(is.na(control_mm)) == F ){
    control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
  }
  
  ASA_control_group[metabolite,] <- control_mm
}

# ASA_control_group <- lapply(rownames(MBX_results), function(metabolite){
#   control_mm <- ASA_control_group[metabolite,]
#  
#   control_mm <- as.numeric(control_mm)
#  
#   if(any(is.na(control_mm)) & all(!is.na(control_mm))){
#     control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
#   }
#   
#   ASA_control_group[metabolite,] <- control_mm
#   ASA_control_group
# })

for(i in 1: nrow(MBX_results)){
  metabolite <- rownames(MBX_results)[i]
  print(i)
  control_mm <- ASA_control_group[metabolite,]
  control_mm <- as.numeric(control_mm)
  treatment_mm <- ASA_treatment_group[metabolite,]
  treatment_mm <- as.numeric(treatment_mm)
  
  if(all(is.na(control_mm)) & (length(which(!is.na(treatment_mm))) >= 4)){
    control_mm = 0.5*min(treatment_mm, na.rm = T)
    treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)}
  ASA_control_group[metabolite,] <- control_mm
}


for(i in 1: nrow(MBX_results)){
  metabolite <- rownames(MBX_results)[i]
  print(i)
  control_mm <- ASA_control_group[metabolite,]
  control_mm <- as.numeric(control_mm)
  treatment_mm <- ASA_treatment_group[metabolite,]
  treatment_mm <- as.numeric(treatment_mm)
  
  if(all(is.na(treatment_mm)) & (length(which(!is.na(control_mm))) >= 4)){
    treatment_mm = 0.5*min(control_mm, na.rm = T)
  }
  ASA_treatment_group[metabolite,] <- treatment_mm
}
# ASA_control_group <- lapply(rownames(MBX_results), function(metabolite){
#   control_mm <- ASA_control_group[metabolite,]
#   control_mm <- as.numeric(control_mm)
#   if(all(is.na(control_mm)) & (length(which(!is.na(treatment_mm))) >= 4)){
#     control_mm = 0.5*min(treatment_mm, na.rm = T)
#     treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)}
#   ASA_control_group[metabolite,] <- control_mm
#   ASA_control_group
# })

getwd()
write.csv(ASA_control_group, file = "outputs/ASA_control_group_results.csv", row.names = TRUE)
write.csv(ASA_treatment_group, file = "outputs/ASA_treatment_group_results.csv", row.names = TRUE)
###这里开始做t test和画图
ttest_all_results_ASA <- lapply(rownames(MBX_results), function(metabolite) {
  control_mm <- ASA_control_group[metabolite,]
  treatment_mm <- ASA_treatment_group[metabolite,]
  # control_mm <- as.numeric(control_mm)
  # treatment_mm <- as.numeric(treatment_mm)
  control_mm <- log(as.numeric(control_mm))
  treatment_mm <- log(as.numeric(treatment_mm))
  if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {NULL}
  else {t.test(control_mm, treatment_mm, na.action=na.omit)}
})

## 
p_values_ASA <- unlist(sapply(ttest_all_results_ASA, function(result){
  result$p.value
}))
p_values_ASA <- p.adjust(p_values_ASA, method = "BH")
mean_dif_ASA <- unlist(sapply(ttest_all_results_ASA, function(result){
  result$estimate[2]/result$estimate[1]
}))
ttest_and_diff_ASA <- data.frame(diff = mean_dif_ASA, p_val = p_values_ASA)
ttest_and_diff_ASA$significance <- ifelse(log2(ttest_and_diff_ASA$diff) > 1 &
                                               -log10(ttest_and_diff_ASA$p_val) > 3, "Significantly Abundant", "No Significance")

table(ttest_and_diff_ASA$significance)
# 3633 significantly abundant
# 78231 not

ggplot(ttest_and_diff_ASA, aes(x = log2(diff), y = -log10(p_val), col = significance)) +
  geom_point() +
  scale_colour_brewer(palette = "Dark2")+
  xlab("Difference in means (Log2 Fold Change)") +
  ylab("-Log10 P-value") +
  ggtitle("5ASA Volcano Plot") +
  theme(text=element_text(size=20),
        legend.position="right",
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        strip.text = element_text(face = "bold"))
#theme_bw()
print(-log10(p_values_ASA))

table(IBDMDB_metadata$data_type)
### 这里用5-ASA的数据替换原本AZA的数据以跑原来的函数
controls_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], ASA_controls_ID)]
treatments_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], ASA_treatment_ID)]
