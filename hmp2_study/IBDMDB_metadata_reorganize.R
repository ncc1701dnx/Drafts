## work 25.07.2023, a windy and rainy Tuesday

##这个主要的作用就是根据5ASA文章作者的说法，重做整个meta data table

########## Followed by the instruction of 5-ASA authors, we re-organize the metadata table of IBDMDB  #####
###########################################################################################################

### Load packaged needed

library(dplyr)
library(tidyr)
#install.packages("table1")
library(purrr)
library(table1) 
library(tidyr)
library(magrittr)
#install.packages("sm")
library(sm)
library(caret) 
library(grid)

### read the metadata_as_we already done
IBDMDB_metadata <- read.csv("/nfs/data/IBDMDB_MBX_data/metadata/hmp2_metadata.csv", header = T)
names(IBDMDB_metadata)
#select columns of interest
IBDMDB_metadata2 <- IBDMDB_metadata %>% select(1:10,"diagnosis","site_name","consent_age","Age.at.diagnosis", "Tube.A..Metabolomics",
                                 93,94,99,100,188,297,298,318,320,322,354,373,375,377,379,381,383,431,478,334,337,342,348,350,352,358,361,366,369)
# divide up data into different chunks
# 根据作者描述只有serology and methylome as well as host_genome have the drug usage data, so really important of chunk the data first 
###########################
meta_WGS <- IBDMDB_metadata2 %>% filter(data_type=="metagenomics")  # 师兄那里可能有用
meta_serology <- IBDMDB_metadata2 %>% filter(data_type=="serology")  #has drug data
meta_hostgene <- IBDMDB_metadata2 %>% filter(data_type=="host_genome") #has drug data
meta_metabolome <- IBDMDB_metadata2 %>% filter(data_type=="metabolomics")
meta_metatranscript <- IBDMDB_metadata2 %>% filter(data_type=="metatranscriptomics")  
meta_methylome <- IBDMDB_metadata2 %>% filter(data_type=="methylome") #has drug data  
## 检查drug usage 信息的完整性
table(meta_methylome$Azathioprine..Imuran..Azasan.)
nrow(meta_methylome)
table(meta_serology$Azathioprine..Imuran..Azasan.)
nrow(meta_serology)
table(meta_hostgene$Azathioprine..Imuran..Azasan.)
nrow(meta_hostgene)
table(meta_metabolome$Azathioprine..Imuran..Azasan.)
## 最后，将用药信息汇总
meta_drugdata <- IBDMDB_metadata2 %>% filter(data_type=="methylome"|data_type=="host_genome"|data_type=="serology")

#oral 5asa current or not (strict: includes only "current")  
meta_drugdata <- meta_drugdata %>% 
  mutate(oral5asa_strict = ifelse(Asacol..mesalamine. == "Current" | Pentasa..mesalamine. == "Current" | Lialda..mesalamine. == "Current" | Apriso..mesalamine.== "Current",1,
                             ifelse(Asacol..mesalamine. == "Never taken" | Pentasa..mesalamine. == "Never taken" | Lialda..mesalamine. == "Never taken" | Apriso..mesalamine.== "Never taken"|
                                      Asacol..mesalamine. == "Taken since last visit" | Pentasa..mesalamine. == "Taken since last visit" | Lialda..mesalamine. == "Taken since last visit" | Apriso..mesalamine.== "Taken since last visit"|
                                      Asacol..mesalamine. == "Taken prior to baseline" | Pentasa..mesalamine. == "Taken prior to baseline" | Lialda..mesalamine. == "Taken prior to baseline" | Apriso..mesalamine.== "Taken prior to baseline",
                                    0,"")))

#oral 5asa current or not (loose: includes "taken since last visit")
meta_drugdata <- meta_drugdata %>% 
  mutate(oral5asa_loose = ifelse(Asacol..mesalamine. == "Current" | Pentasa..mesalamine. == "Current" | Lialda..mesalamine. == "Current" | Apriso..mesalamine.== "Current"|
                                Asacol..mesalamine. == "Taken since last visit" | Pentasa..mesalamine. == "Taken since last visit" | Lialda..mesalamine. == "Taken since last visit" | Apriso..mesalamine.== "Taken since last visit",1,
                              ifelse(Asacol..mesalamine. == "Never taken" | Pentasa..mesalamine. == "Never taken" | Lialda..mesalamine. == "Never taken" | Apriso..mesalamine.== "Never taken"|
                                       Asacol..mesalamine. == "Taken prior to baseline" | Pentasa..mesalamine. == "Taken prior to baseline" | Lialda..mesalamine. == "Taken prior to baseline" | Apriso..mesalamine.== "Taken prior to baseline",
                                     0,"")))

#PR5asa (loose) 
meta_drugdata <- meta_drugdata %>% 
  mutate(pr5asa = ifelse(Canasa.suppositories..mesalamine.suppositories. == "Current" | Rowasa.enemas..mesalamine.enemas. == "Current"|
                           Canasa.suppositories..mesalamine.suppositories. == "Taken since last visit" | Rowasa.enemas..mesalamine.enemas. == "Taken since last visit",1,
                         ifelse(Canasa.suppositories..mesalamine.suppositories. == "Never taken" | Rowasa.enemas..mesalamine.enemas. == "Never taken"|
                                  Canasa.suppositories..mesalamine.suppositories. == "Taken prior to baseline" | Rowasa.enemas..mesalamine.enemas. == "Taken prior to baseline",
                                0,"")))

#Bond5asa (loose) 
meta_drugdata <- meta_drugdata %>% 
  mutate(bond5asa  = ifelse(Dipentum..olsalazine. == "Current" | Colozal..balasalizide. == "Current" | Sulfasalizine..Azulfidine. == "Current" |
                              Dipentum..olsalazine. == "Taken since last visit" | Colozal..balasalizide. == "Taken since last visit" | Sulfasalizine..Azulfidine. == "Taken since last visit" ,1,
                            ifelse(Dipentum..olsalazine. == "Never taken" | Colozal..balasalizide. == "Never taken" | Sulfasalizine..Azulfidine. == "Never taken" |
                                     Dipentum..olsalazine. == "Taken prior to baseline" | Colozal..balasalizide. == "Taken prior to baseline" | Lialda..mesalamine. == "Taken prior to baseline",
                                   0,"")))

#Any 5asa use (loose)
meta_drugdata <- meta_drugdata %>% mutate(any5asa = ifelse(oral5asa_loose == 1 | pr5asa  == 1 | bond5asa  == 1,1,ifelse(oral5asa_loose == 0 | pr5asa  == 0 | bond5asa  == 0,0,"")))

#Other drugs 
##############
#biologics
meta_drugdata <- meta_drugdata %>% 
  mutate(biologics = ifelse(Remicade..Infliximab. == "Current" | Humira..Adalimumab. == "Current" | Cimzia..Certlizumab. == "Current" | Tysabri..Natalizumab.== "Current"|
                              Remicade..Infliximab. == "Taken since last visit" | Humira..Adalimumab. == "Taken since last visit" | Cimzia..Certlizumab. == "Taken since last visit" | Tysabri..Natalizumab.== "Taken since last visit",1,
                            ifelse(Remicade..Infliximab. == "Never taken" | Humira..Adalimumab. == "Never taken" | Cimzia..Certlizumab. == "Never taken" | Tysabri..Natalizumab.== "Never taken"|
                                     Remicade..Infliximab. == "Taken prior to baseline" | Humira..Adalimumab. == "Taken prior to baseline" | Cimzia..Certlizumab. == "Taken prior to baseline" | Tysabri..Natalizumab.== "Taken prior to baseline",0,"")))

#steroids 
meta_drugdata <- meta_drugdata %>% 
  mutate(steroids = ifelse(Prednisone == "Current" | Entocort..Budesonide. == "Current" | Solumedrol..Medrol. == "Current"|
                             Prednisone == "Taken since last visit" | Entocort..Budesonide. == "Taken since last visit" | Solumedrol..Medrol. == "Taken since last visit",1,
                           ifelse(Prednisone == "Never taken" | Entocort..Budesonide. == "Never taken" | Solumedrol..Medrol. == "Never taken"|
                                    Prednisone == "Taken prior to baseline" | Entocort..Budesonide. == "Taken prior to baseline" | Solumedrol..Medrol. == "Taken prior to baseline",0,"")))

#immunomodulators
meta_drugdata <- meta_drugdata %>% 
  mutate(immunomod = ifelse(Azathioprine..Imuran..Azasan. == "Current" | Methotrexate == "Current" | Mercaptopurine..Purinethol..6MP. == "Current" |
                              Azathioprine..Imuran..Azasan. == "Taken since last visit" | Methotrexate == "Taken since last visit" | Mercaptopurine..Purinethol..6MP. == "Taken since last visit",1,
                            ifelse(Azathioprine..Imuran..Azasan. == "Never taken" | Methotrexate == "Never taken" | Mercaptopurine..Purinethol..6MP. == "Never taken"|
                                     Azathioprine..Imuran..Azasan. == "Taken prior to baseline" | Methotrexate == "Taken prior to baseline" | Mercaptopurine..Purinethol..6MP. == "Taken prior to baseline",0,"")))
#AZA and 6MP specific 
meta_drugdata <- meta_drugdata%>% 
  mutate(aza6mp = ifelse(Azathioprine..Imuran..Azasan. == "Current" | Mercaptopurine..Purinethol..6MP. == "Current" | 
                           Azathioprine..Imuran..Azasan. == "Taken since last visit" | Mercaptopurine..Purinethol..6MP. == "Taken since last visit", 1, 
                         ifelse(Azathioprine..Imuran..Azasan. == "Never taken" | Mercaptopurine..Purinethol..6MP. == "Never taken" | 
                                  Azathioprine..Imuran..Azasan. == "Taken prior to baseline" | Mercaptopurine..Purinethol..6MP. == "Taken prior to baseline",0,"")))

table(meta_drugdata$aza6mp)
table(meta_drugdata$any5asa)
length(unique(meta_drugdata$Participant.ID))
####################################################################################################################################################################################
################################  Data Clean and bind ##############################################################################################################################

###### 1: clean 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111

#######################################
# clean the dataset for AZA6MP analysis
#######################################
## 将participantID和visit num联合起来是很好的，这样就可以判断一个病人的信息是否重复
## week.num也是一样的，不过5ASA的文章作者用的是visit num 那我也就照做了，免得重复结果不对

nrow(meta_drugdata)
meta_drugdata$id <- paste(meta_drugdata$Participant.ID,meta_drugdata$visit_num,sep=".")
meta_drugdata$aza6mp[which(meta_drugdata$id == "M2025.31")]
meta_drugdata$aza6mp[which(meta_drugdata$id == "C3002.31")]
colnames(meta_drugdata[,"id"])
unique(meta_drugdata$id)
names(meta_drugdata)
meta_drugdata2 <- meta_drugdata %>% distinct(id, .keep_all = TRUE) %>% arrange(Participant.ID, visit_num) %>% select(id,45:53,21,22,2,3,9) #gets rid of duplicate rows/ If you set all drug info
#meta_drugdata2 <- meta_drugdata %>% distinct(id, .keep_all = TRUE) %>% arrange(Participant.ID, visit_num) %>% select(id,45,46,47,21,22,2,3,9) #gets rid of duplicate rows/ If you set only aza6mp info
## 我们这里使用44到52，因为这些是我们上一步加上的从5asa到aza6mp的用药信息
## 2是"External.ID"，3是"Participant.ID"，不要弄错了，后面关联表格会用到
## 9是"visit_num"，用于处理同一个病人的不同服药数据

###########################################################################
########### MBX MBX MBX MBX MBX MBX MBX MBX MBX MBX MBX MBX ###############
###########################################################################

#########################################
# also clean MBX data for AZA6MP analysis
#########################################

meta_metabolome$id <- paste(meta_metabolome$Participant.ID,meta_metabolome$visit_num,sep=".")
names(meta_metabolome)
unique(meta_metabolome$id)
## ...the unique id in meta_metabolome is more than those have drug usage information in metadata table
## hard...
## same as before, subset the data
meta_metabolome2 <- meta_metabolome[,c(2,3,9,11,12,13:19,33,34,45)]

###### 2: bind 22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222

##############################################
# Bind the MBX table with the drug usage table
##############################################
## 使用dplyr进行binding，这样即使id不匹配也会合并
## 有问题，这样进来的id是混乱的，时间线是不匹配的
## 用arrange函数，使用新创建的ParticipantID和visit num进行排序
colnames(meta_drugdata2)
colnames(meta_data_pre)
meta_data_pre <- full_join(meta_metabolome2, meta_drugdata2,by="id") %>% 
  mutate(Participant.ID=coalesce(Participant.ID.x,Participant.ID.y), visit_num = coalesce(visit_num.x,visit_num.y)) %>%
  arrange(Participant.ID, visit_num) %>% # 排序
  select(-Participant.ID.x,-Participant.ID.y,-visit_num.x,-visit_num.y) %>% # 排序完成后可以删掉新造的两个列
 rename(age_diagnostic=Age.at.diagnosis,hospital_last_2weeks=X5..In.the.past.2.weeks..have.you.been.hospitalized.,
         bowel_surgery=X6..Have.you.ever.had.bowel.surgery.,CRP_mg_L=CRP..mg.L.,ESR_mm_hr=ESR..mm.hr.)
################################################################################
unique(meta_data_pre$id)
names(meta_data_pre)
## 将我们需要的服药数据确认转换为数字型
meta_data_pre <- meta_data_pre %>% mutate_at(14:22,as.numeric) # all drug info 
#meta_data_pre <- meta_data_pre %>% mutate_at(14:15,as.numeric) # aza6mp only
#carry forward within each participant so that non-measured data matches with stool data 
#this will be done first forward, then for those missing prior to baseline, will carry backwards
################################################################################

## if all drugs
meta_data <- meta_data_pre %>% 
  group_by(Participant.ID) %>% 
  fill(c(diagnosis:bowel_surgery,oral5asa_strict:aza6mp),.direction="down")%>%
  fill(c(diagnosis:bowel_surgery,oral5asa_strict:aza6mp),.direction="up") 
## if aza6mp only
if(F){
meta_data <- meta_data_pre %>% 
  group_by(Participant.ID) %>% 
  fill(c(diagnosis:bowel_surgery,immunomod:aza6mp),.direction="down")%>%
  fill(c(diagnosis:bowel_surgery,immunomod:aza6mp),.direction="up") 
}
## one thing haven't be changed: the External ID
names(meta_data)
meta_data <- meta_data %>% 
  mutate(External.ID =coalesce(External.ID.x, External.ID.y)) %>%
  select(-External.ID.x, -External.ID.y)


#for non IBD patients with missing drug data, set to 0

# 必须要牢记的是我们这里的meta_data是MBX的metadata
meta_data[meta_data$diagnosis=="nonIBD",] %<>% mutate_at(vars(oral5asa_strict:aza6mp), ~replace(., is.na(.), 0)) # all drug info
meta_data[meta_data$diagnosis=="nonIBD",] %<>% mutate_at(vars(immunomod:aza6mp), ~replace(., is.na(.), 0)) # aza6mp only

cols <- c("oral5asa_strict","oral5asa_loose","aza6mp","diagnosis","biologics",
          "site_name","Antibiotics","bowel_surgery",
          "Participant.ID","steroids","pr5asa","bond5asa")
#meta_data[cols] <- lapply(meta_data[cols], factor) #converts these to a factor variable)

#meta_data$diagnosis <- relevel(meta_data$diagnosis, ref = "nonIBD") # Set nonIBD as reference

## Write it out, the server is not stable...Ahhhhhh!
write.csv(meta_data,"outputs/reorganized_MBX_metatable.csv")

meta_data_IBD <- meta_data %>% filter(diagnosis=="CD"|diagnosis=="UC")

###############################################
### 下面的是将meta data和MBX结果数据联系起来
###############################################

### 检查下有哪些病人存在MBX数据
names(meta_data) # in total 
MBX_results <- read.csv("inputs/iHMP_metabolomics.csv", header = T)
MBX_results <- MBX_results %>%
  mutate(Metabolite = if_else(Metabolite == "", 
                              paste0(m.z, "_", RT), 
                              Metabolite))
names(MBX_results) # in total 553
length(which(meta_data$External.ID %in% names(MBX_results))) # 546 match in this table

length(unique(meta_data$Participant.ID)) ## 106, almost same, as we have 130 patients in total
length(unique(names(MBX_results)[-(1:7)]))
length(unique(meta_data$External.ID))

###########################################################################
########### MGX MGX MGX MGX MGX MGX MGX MGX MGX MGX MGX MGX ###############
###########################################################################

### 第一步的前述步骤请去MBX的分界线的前面的部分找

#########################################
# and clean MGX data for AZA6MP analysis
#########################################

meta_WGS$id <- paste(meta_WGS$Participant.ID,meta_WGS$visit_num,sep=".")
names(meta_WGS)
unique(meta_WGS$id)
## ...the unique id in meta_metabolome is more than those have drug usage information in metadata table
## hard...
## same as before, subset the data
meta_WGS2 <- meta_WGS[,c(2,3,9,11,12,13:19,33,34,45)]
colnames(meta_WGS[,c(2,3,9,11,12,13:19,33,34,45)])
###### 2: bind 22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222

##############################################
# Bind the MGX table with the drug usage table
##############################################
## 使用dplyr进行binding，这样即使id不匹配也会合并
## 有问题，这样进来的id是混乱的，时间线是不匹配的
## 用arrange函数，使用新创建的ParticipantID和visit num进行排序
meta_MGXdata_pre <- full_join(meta_WGS2, meta_drugdata2,by="id") %>% 
  mutate(Participant.ID=coalesce(Participant.ID.x,Participant.ID.y), visit_num = coalesce(visit_num.x,visit_num.y)) %>%
  arrange(Participant.ID, visit_num) %>%
  select(-Participant.ID.x,-Participant.ID.y,-visit_num.x,-visit_num.y) %>% # 排序完成后可以删掉新造的两个列
  rename(age_diagnostic=Age.at.diagnosis,hospital_last_2weeks=X5..In.the.past.2.weeks..have.you.been.hospitalized.,
         bowel_surgery=X6..Have.you.ever.had.bowel.surgery.,CRP_mg_L=CRP..mg.L.,ESR_mm_hr=ESR..mm.hr.)

################################################################################
unique(meta_MGXdata_pre$id)
names(meta_MGXdata_pre)
## 将我们需要的服药数据确认转换为数字型
meta_MGXdata_pre <- meta_MGXdata_pre %>% mutate_at(14:22,as.numeric) 
#carry forward within each participant so that non-measured data matches with stool data 
#this will be done first forward, then for those missing prior to baseline, will carry backwards
################################################################################

meta_MGXdata <- meta_MGXdata_pre %>% 
  group_by(Participant.ID) %>% 
  fill(c(diagnosis:bowel_surgery,oral5asa_strict:aza6mp),.direction="down")%>%
  fill(c(diagnosis:bowel_surgery,oral5asa_strict:aza6mp),.direction="up") 
## one thing haven't be changed: the External ID
names(meta_MGXdata)
meta_MGXdata <- meta_MGXdata %>% 
  mutate(External.ID =coalesce(External.ID.x, External.ID.y)) %>%
  select(-External.ID.x, -External.ID.y)

#for non IBD patients with missing drug data, set to 0


meta_MGXdata[meta_MGXdata$diagnosis=="nonIBD",] %<>% mutate_at(vars(oral5asa_strict:aza6mp), ~replace(., is.na(.), 0))

cols <- c("oral5asa_strict","oral5asa_loose","aza6mp","diagnosis","biologics",
          "site_name","Antibiotics","bowel_surgery",
          "Participant.ID","steroids","pr5asa","bond5asa")
#meta_MGXdata[cols] <- lapply(meta_MGXdata[cols], factor) #converts these to a factor variable)

#meta_MGXdata$diagnosis <- relevel(meta_MGXdata$diagnosis, ref = "nonIBD") # Set nonIBD as reference

## Write it out, the server is not stable...Ahhhhhh!
write.csv(meta_MGXdata,"outputs/reorganized_MGX_metatable.csv")


###############################################################################################################################################################
####################### 开始做分析 ############################################################################################################################

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
length(IDs)
length(which(meta_data$any5asa == 1))
length(treatments_ID)
treatments_ID <- meta_data %>% 
  filter(Participant.ID %in% IDs & any5asa == 1)
treatments_ID <- treatments_ID$External.ID
controls_ID <- meta_data %>% 
  filter(Participant.ID %in% IDs & any5asa == 0)
controls_ID <- controls_ID$External.ID
dim(GMPS_rpkm)
