### Work at 09.July.2023
### The main goal is to write functions to give a list to those suspisious compounds

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
getwd()

### 然后，读取需要的数据，这里主要是两个：IBDMDB的metadata以及所有的metabolites list

##### 本次一共有三种文件，一种是整个IBDMDB数据的metadata，包括MBX，MGX，MPX等等即IBDMDB_metadata
###### 另一种是4个处理方式：C18pos, C18neg, HILpos, HILneg的质谱文件数据即Inventory文件
####### 还有一种是最后MBX文件的质谱结果数据即MBX_results文件

########## 先对metadata table进行处理
# 这是metadata
IBDMDB_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/hmp2_metadata.csv", header = T)
names(IBDMDB_metadata)
# 变成只有MBX的metadata降低表格长度
MBX_IBDMDB_metadata <- IBDMDB_metadata[which(IBDMDB_metadata$data_type == "metabolomics"),]

########### 以下是重构IBDMDB的metadata的table的环节 ##################################################################################
# 原metadata table实在是太差了，MBX数据里面没有任何的病人用药信息
# 用唯一matchable的ID：Participant.ID来计算总和的病人用药情况，然后插入MBX的数据里
########################################################################################################################################
# Firstly is about AZA:
table(IBDMDB_metadata$Azathioprine..Imuran..Azasan.)
#check who never used
aza_contrl <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Azathioprine..Imuran..Azasan. == "Never taken")]
which(MBX_IBDMDB_metadata$Participant.ID %in% aza_contrl)
MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.[which(MBX_IBDMDB_metadata$Participant.ID %in% aza_contrl)] = "Never taken"
#check who used before
aza_prior <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Azathioprine..Imuran..Azasan. == "Taken prior to baseline")]
which(MBX_IBDMDB_metadata$Participant.ID %in% aza_prior)
MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.[which(MBX_IBDMDB_metadata$Participant.ID %in% aza_prior)] = "Taken prior to baseline"
#check who are using
aza_treat <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Azathioprine..Imuran..Azasan. == "Current")]
which(MBX_IBDMDB_metadata$Participant.ID %in% aza_treat)
MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.[which(MBX_IBDMDB_metadata$Participant.ID %in% aza_treat)] = "Current"
table(MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.)
length(which(aza_treat %in% aza_contrl))
# Then 6MP
table(MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.)
table(IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.)
#check who never used
mp_contrl <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Mercaptopurine..Purinethol..6MP. == "Never taken")]
which(MBX_IBDMDB_metadata$Participant.ID %in% mp_contrl)
MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.[which(MBX_IBDMDB_metadata$Participant.ID %in% mp_contrl)] = "Never taken"
#check who used before
mp_prior <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Mercaptopurine..Purinethol..6MP. == "Taken prior to baseline")]
which(MBX_IBDMDB_metadata$Participant.ID %in% mp_prior)
MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.[which(MBX_IBDMDB_metadata$Participant.ID %in% mp_prior)] = "Taken prior to baseline"
#check who are using
mp_treat <- IBDMDB_metadata$Participant.ID[which(IBDMDB_metadata$Mercaptopurine..Purinethol..6MP. == "Current")]
which(MBX_IBDMDB_metadata$Participant.ID %in% mp_treat)
MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.[which(MBX_IBDMDB_metadata$Participant.ID %in% mp_treat)] = "Current"
table(MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.)
length(which(mp_treat %in% mp_contrl))
length(mp_contrl)
length(aza_contrl)
# Now save this changed MBX metadata table
write.table(MBX_IBDMDB_metadata, "R_Works/outputs/MBX_IBDMDB_metadata.csv", sep = ",", quote = F,col.names = T, row.names = F)
########################################################################################################################################
tr1_mp <- treatments$Participant.ID
tr2_mp <- treatments$Participant.ID
length(tr1_mp)
length(unique(controls$Participant.ID))
unique(treatments$Participant.ID)
######################### 以下是分割缩小MBX metadata table的环节 ######################################################################
###### 进行的动作为 1：分割MBX数据由全部病人的数据变为仅有UC患者的数据
######## 2：分割出controls和treatments的ID以方便做t test
######################################################################################################################################
# 细分MBX的metadata，仅需要UC的patients
table(MBX_IBDMDB_metadata$diagnosis)
MBX_IBDMDB_metadata <- MBX_IBDMDB_metadata[which(MBX_IBDMDB_metadata$diagnosis == "UC"),]
nrow(MBX_IBDMDB_metadata)
# 细分MBX的metadata，分割出control和treatments
table(MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.)
table(MBX_IBDMDB_metadata$Azathioprine..Imuran..Azasan.)
treatments <- subset(MBX_IBDMDB_metadata,Mercaptopurine..Purinethol..6MP. == "Current" | Azathioprine..Imuran..Azasan. == "Current")
controls <- subset(MBX_IBDMDB_metadata,Mercaptopurine..Purinethol..6MP. == "Never taken" & Azathioprine..Imuran..Azasan. == "Never taken")
table(MBX_IBDMDB_metadata$Mercaptopurine..Purinethol..6MP.)
controls_ID <- controls$External.ID
treatments_ID <- treatments$External.ID
#########################################################################################################################################
length(treatments_ID)
length(controls_ID)

# 读取IBDMDB 质谱文件数据的MBX文件名
################### 如果我们不自己通过XCMS管道做MBX research的话可以不用读这些文件 ########################################
#IBDMDB_c18neg_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/C18n_RawFileInventory.csv", header = T)          #
#IBDMDB_c18pos_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/C8p_RawFileInventory.csv", header = T)           #
#IBDMDB_HILneg_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/HILn_RawFileInventory.csv", header = T)          #
#IBDMDB_HILpos_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/HILp_RawFileInventory.csv", header = T)          #
###########################################################################################################################


###### 接下来是处理MBX_results文件 ######################################################################################################
# 读取MBX的文件的results文件
## 然后将MBX results分割成treatment和control两个表格，可以进行t test#####################################################################
MBX_results <- read.csv("inputs/iHMP_metabolomics.csv", header = T)

####重要！ 分割表格前一定记得加上行名方便日后对齐
row.names(MBX_results) = MBX_results$Compound
length(row.names(MBX_results))
length(unique(row.names(MBX_results)))
# unique 81867, and total length 81867

#### 检查下controls和treatments能不能和results data对上
treatments_ID %in% colnames(MBX_results)
controls_ID %in% colnames(MBX_results)
# 全部可以对上
length(treatments_ID)
length(controls_ID)
# 将MBX文件名文件和metadata文件联系起来
## 需要在列名里面加上【1：7】而不仅仅是controls以及treatments，因为这样带元素信息（m/z，RT等）方便到时候做表
controls_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], controls_ID)]
treatments_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], treatments_ID)]
#########################################################################################################################################
## 如果不是自己从XCMS管道跑MBX数据，则不用跑下列函数
if(F){
bind_metadata_MBXinventory <- function(Inventory, MBX_Metadata){
  azaList <- c()
  mpList <- c()
  diaList <- c()
  for (i in 1:nrow(Inventory)) {
    #i = 1 #test
    if(Inventory$Sample.ID[i] %in% MBX_Metadata$Tube.A..Metabolomics){
      tt <- which(MBX_Metadata$Tube.A..Metabolomics == Inventory$Sample.ID[i])
      azaList <- c(azaList, MBX_Metadata$Azathioprine..Imuran..Azasan.[tt])
      mpList <- c(mpList, MBX_Metadata$Mercaptopurine..Purinethol..6MP.[tt])
      diaList <- c(diaList, MBX_Metadata$diagnosis[tt])
    } else {
      print(paste0(Inventory$Sample.ID[i], " not found in metadata"))
      azaList <- c(azaList, NA)
      mpList <- c(mpList, NA)
      diaList <- c(diaList, NA)
    }
  }
  Inventory$Azathioprine..Imuran..Azasan. =  azaList
  Inventory$Mercaptopurine..Purinethol..6MP. = mpList
  Inventory$diagnosis = diaList
  Inventory
}
# run这个函数
IBDMDB_c18pos_metadata <- bind_metadata_MBXinventory(IBDMDB_c18pos_metadata,MBX_IBDMDB_metadata)
IBDMDB_c18neg_metadata <- bind_metadata_MBXinventory(IBDMDB_c18neg_metadata,MBX_IBDMDB_metadata)
IBDMDB_HILneg_metadata <- bind_metadata_MBXinventory(IBDMDB_HILneg_metadata,MBX_IBDMDB_metadata)
IBDMDB_HILpos_metadata <- bind_metadata_MBXinventory(IBDMDB_HILpos_metadata, MBX_IBDMDB_metadata)
}

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
## 现在写三个函数，分别是1：自动做表，2：自动加significance，3：自动绘volcano plot
## IBDMDB_MBX_check这个函数的用处是对于给定的化合物，自动生成表格并且做t test
######################################################################################################################################
## 现在输入这个自动做表的函数：
IBDMDB_MBX_check <- function(formula){
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
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode)
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
      treatment_mm <- as.numeric(subset_treatment[tt[i],-(1:7)])
      if(any(is.na(treatment_mm)) == T & all(is.na(treatment_mm)) == F){
        treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T) 
      }
      if(any(is.na(control_mm)) == T & all(is.na(control_mm)) == F ){
        control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
      }
      if(all(is.na(control_mm)) & (length(which(!is.na(treatment_mm))) >= 4)){
        control_mm = 0.5*min(treatment_mm, na.rm = T)
        treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)
        }
      if(all(is.na(treatment_mm)) & (length(which(!is.na(control_mm))) >= 4)){
        treatment_mm = 0.5*min(control_mm, na.rm = T)
      }
      control_mm <- log(control_mm)
      treatment_mm <- log(treatment_mm)
      if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {
        t_test_result <- NA
        p_val <- NA
        fold_change <- NA}
      else {
        t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
        p_val <- t_test_result$p.value
        fold_change <- t_test_result$estimate[2]/t_test_result$estimate[1]
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
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode)
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
      treatment_mm <- as.numeric(subset_treatment[tt[i],-(1:7)])
      if(any(is.na(treatment_mm)) == T & all(is.na(treatment_mm)) == F){
        treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T) 
      }
      if(any(is.na(control_mm)) == T & all(is.na(control_mm)) == F ){
        control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
      }
      if(all(is.na(control_mm)) & (length(which(!is.na(treatment_mm))) >= 4)){
        control_mm = 0.5*min(treatment_mm, na.rm = T)
        treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)
        }
      if(all(is.na(treatment_mm)) & (length(which(!is.na(control_mm))) >= 4)){
        treatment_mm = 0.5*min(control_mm, na.rm = T)
      }
      control_mm <- log(control_mm)
      treatment_mm <- log(treatment_mm)
      if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {
        t_test_result <- NA
        p_val <- NA
        fold_change <- NA}
      else {
        t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
        p_val <- t_test_result$p.value
        fold_change <- t_test_result$estimate[2]/t_test_result$estimate[1]
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
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode)
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
      treatment_mm <- as.numeric(subset_treatment[tt[i],-(1:7)])
      if(any(is.na(treatment_mm)) == T & all(is.na(treatment_mm)) == F){
        treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T) 
      }
      if(any(is.na(control_mm)) == T & all(is.na(control_mm)) == F ){
        control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
      }
      if(all(is.na(control_mm)) & (length(which(!is.na(treatment_mm))) >= 4)){
        control_mm = 0.5*min(treatment_mm, na.rm = T)
        treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T)
        }
      if(all(is.na(treatment_mm)) & (length(which(!is.na(control_mm))) >= 4)){
        treatment_mm = 0.5*min(control_mm, na.rm = T)
      }
      control_mm <- log(control_mm)
      treatment_mm <- log(treatment_mm)
      if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {
        t_test_result <- NA
        p_val <- NA
        fold_change <- NA}
      else {
        t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
        p_val <- t_test_result$p.value
        fold_change <- t_test_result$estimate[2]/t_test_result$estimate[1]
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
  com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode)
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
      treatment_mm <- as.numeric(subset_treatment[tt[i],-(1:7)])
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
      if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {
        t_test_result <- NA
        p_val <- NA
        fold_change <- NA}
      else {
        t_test_result <- t.test(control_mm, treatment_mm, na.action=na.omit)
        p_val <- t_test_result$p.value
        fold_change <- t_test_result$estimate[2]/t_test_result$estimate[1]
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
    theme(text=element_text(size=20),
          legend.position="right",
          strip.text = element_text(face = "bold"))
    # +theme_bw()
}

###再写一个函数，能够根据给出的结果表格找到原始的MBX_results表格的数据画box plot
IBDMDB_6MP$m_z[which(IBDMDB_6MP$significance == "Significantly Abundant")]

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
  group = c(rep("Treatment", length(treatment_mm)), rep("Control", length(control_mm)))
)

# ggplot2 画画图 
# 根据p值确定标签的显示方式
label_p <- ifelse(p_value < 0.001, "< 0.001", formatC(p_value, format = "e", digits = 2))

# 使用ggplot2绘图
p <- ggplot(combined_data, aes(x=group, y=abundance, color=group)) + 
  geom_boxplot(outlier.shape = NA, lwd = 1.5) + 
  geom_jitter(position=position_jitter(0.2), size=4) +
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
  theme(text=element_text(size=20),
        legend.position="right",
        strip.text = element_text(face = "bold"))
plot(p)
}
}


############ 所有的函数写作完成 #########################################################################################################

##################### 现在开始用这三个函数计算全部感兴趣的化合物的状态 ##################################################################


################# Firstly the AZA/6MP compounds
## 使用6-MeMP
IBDMDB_6MeMP <- IBDMDB_MBX_check(formula = "C6H6N4S1")
IBDMDB_6MeMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6MeMP)
IBDMDB_MBX_plot(IBDMDB_6MeMP, "6MeMP")
IBDMDBresults_boxplot(IBDMDB_6MeMP, "6MeMP")
# 6TUA
IBDMDB_6TUA <- IBDMDB_MBX_check(formula = "C5H4N4O2S1")
IBDMDB_6TUA$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TUA)
IBDMDB_MBX_plot(IBDMDB_6TUA, "6TUA")
IBDMDBresults_boxplot(IBDMDB_6TUA, "6TUA")
# 6-MeTIMP
IBDMDB_6MeTIMP <- IBDMDB_MBX_check(formula = "C11H15N4O7P1S1")
IBDMDB_6MeTIMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6MeTIMP)
IBDMDB_MBX_plot(IBDMDB_6MeTIMP, "6MeTIMP")
IBDMDBresults_boxplot(IBDMDB_6MeTIMP, "6MeTIMP")
# 6-TIMP
IBDMDB_6TIMP <- IBDMDB_MBX_check(formula = "C10H13N4O7P1S1")
IBDMDB_6TIMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TIMP)
IBDMDB_MBX_plot(IBDMDB_6TIMP, "6TIMP")
IBDMDBresults_boxplot(IBDMDB_6TIMP, "6TIMP")
# AZA
IBDMDB_AZA <- IBDMDB_MBX_check(formula = "C9H7N7O2S1")
IBDMDB_AZA$significance <- IBDMDB_MBX_add_signif(IBDMDB_AZA)
IBDMDB_MBX_plot(IBDMDB_AZA, "AZA")
IBDMDBresults_boxplot(IBDMDB_AZA, "AZA")
# 6-MP
IBDMDB_6MP <- IBDMDB_MBX_check(formula = "C5H4N4S1")
IBDMDB_6MP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6MP)
IBDMDB_MBX_plot(IBDMDB_6MP, "6MP")
IBDMDBresults_boxplot(IBDMDB_6MP,"6MP")
# 6-TGMP
IBDMDB_6TGMP <- IBDMDB_MBX_check(formula = "C10H14N5O7P1S1")
IBDMDB_6TGMP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TGMP)
IBDMDB_MBX_plot(IBDMDB_6TGMP, "6TGMP")
IBDMDBresults_boxplot(IBDMDB_6TGMP, "6TGMP")
# 6-TGDP
IBDMDB_6TGDP <- IBDMDB_MBX_check(formula = "C10H15N5O10P2S1")
IBDMDB_6TGDP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6TGDP)
IBDMDB_MBX_plot(IBDMDB_6TGDP, "6TGDP")
IBDMDBresults_boxplot(IBDMDB_6TGDP, "6TGDP")
# 6-TGTP
IBDMDB_6_TGTP <- IBDMDB_MBX_check(formula = "C10H16N5O13P3S1")
IBDMDB_6_TGTP$significance <- IBDMDB_MBX_add_signif(IBDMDB_6_TGTP)
IBDMDB_MBX_plot(IBDMDB_6_TGTP, "6TGTP")
IBDMDBresults_boxplot(IBDMDB_6_TGTP, "6TGTP")



###################### Then 5-ASA Compounds ########################
### Be careful! Use IBDMDB_5ASA_check.R to re construct the datasets firstly

# 5-ASA
IBDMDB_5_ASA <- IBDMDB_MBX_check(formula = "C7H7N1O3")
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
  