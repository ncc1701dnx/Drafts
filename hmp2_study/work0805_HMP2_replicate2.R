MBX_results$Metabolite[100:105]
MBX_results$m.z[100:105]
MBX_results$RT[100:105]
names(MBX_results)[1:10]

### 将MBX结果的表格按照HMP2发现5ASA作者的思路继续整理
MBX_results <- MBX_results %>%
  mutate(Metabolite = if_else(Metabolite == "", 
                              paste0(m.z, "_", RT), 
                              Metabolite))
which(c("carboxyibuprofen","metronidazole","salicylate", "sulfapyridine" ,"154.0502_3.83") %in% MBX_results$Metabolite)
which(MBX_results$Metabolite == "154.0502_3.83")

MBX_results$Metabolite

MBX_5asa_result <- MBX_results %>% select(-"Pooled.QC.sample.CV",-"m.z",-"HMDB...Representative.ID.",-"Method",-"RT",-"Compound") %>%
  filter(Metabolite == "154.0502_3.83")

View(MBX_5asa_rotated)
MBX_5asa_rotated = setNames(data.frame(t(MBX_5asa_result[,-1])), MBX_5asa_result[,1])
MBX_5asa_rotated <- cbind(Metabolite=row.names(MBX_5asa_rotated),MBX_5asa_rotated)
row.names(MBX_5asa_rotated) <- NULL
class(MBX_5asa_rotated$Metabolite)
MBX_5asa_rotated <- MBX_5asa_rotated %>% mutate_all(~replace(., is.na(.), 0))
MBX_5asa_rotated <- MBX_5asa_rotated %>% rename(External.ID = Metabolite)
MBX_5asa_rotated <- MBX_5asa_rotated %>% rename(X154.0502_3.83 = "154.0502_3.83") ### 纯数字名字导致很多问题
meta_data_IBD <- meta_data %>% filter(diagnosis=="CD"|diagnosis=="UC")
names(meta_data_IBD)
names(MBX_5asa_rotated)
MBX_5asa_bind <- inner_join(meta_data_IBD,MBX_5asa_rotated,by= "External.ID") %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0))
names(MBX_5asa_bind)

table(MBX_5asa_bind$binary154)


## loose means all patients have changed the drug usage will be considered
binary154_check <- MBX_5asa_bind %>%
  # Group by Participant.ID
  group_by(Participant.ID) %>%
  
  mutate(change_in_binary5asa = sum(abs(diff(binary154)))) %>%
  
  # Ungroup to operate on the full data again
  ungroup() %>%
  
  # Filter rows for participants with changes in any5asa
  filter(change_in_binary5asa > 0) %>%
  
  # Select required columns
  select(Participant.ID, visit_num, binary154, External.ID) %>%
  
  # Remove duplicate rows
  distinct()
View(binary154_check)
unique(binary154_check$Participant.ID)
## [1] "C3003" "C3004" "C3015" "C3016" "C3017" "C3031" "C3032" "E5004" "H4007" "H4014" "H4015" "H4035" "H4040" "H4042" "M2008" "M2014" "M2021" "M2028" "M2068" "M2071" "M2085" "P6005" "P6009" "P6010"
## [25] "P6012" "P6016" "P6024" "P6033"
unique(data_5asa$SampleID)
unique(data_5asa$Participant.ID)
## [1] "C3004" "C3031" "H4014" "H4015" "H4035" "H4040" "M2008" "M2028" "M2071" "P6009" "P6010" "P6012" "P6016"
unique(result_5asa_loose$Participant.ID)

all(unique(result_5asa_loose$Participant.ID) %in% unique(binary154_check$Participant.ID))
all(unique(data_5asa$Participant.ID) %in% unique(binary154_check$Participant.ID))

which(MBX_results$Metabolite == "196.0609_2.81")

###############################################################################################################################
MBX_154196 <- MBX_results %>% select(-"Pooled.QC.sample.CV",-"m.z",-"HMDB...Representative.ID.",-"Method",-"RT",-"Compound") %>%
  filter(Metabolite == "154.0502_3.83" | Metabolite == "196.0609_2.81")

MBX_154196_rotated = setNames(data.frame(t(MBX_154196[,-1])), MBX_154196[,1])
MBX_154196_rotated <- cbind(Metabolite=row.names(MBX_154196_rotated),MBX_154196_rotated)
row.names(MBX_154196_rotated) <- NULL

MBX_154196_rotated <- MBX_154196_rotated %>% mutate_all(~replace(., is.na(.), 0))
MBX_154196_rotated <- MBX_154196_rotated %>% rename(External.ID = Metabolite)
MBX_154196_rotated <- MBX_154196_rotated %>% rename(X154.0502_3.83 = "154.0502_3.83") %>% rename(X196.0609_2.81 = "196.0609_2.81") ### 纯数字名字导致很多问题
View(MBX_154196_rotated)

MBX_154196_bind <- inner_join(meta_data_IBD,MBX_154196_rotated,by= "External.ID") %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% 
  mutate(binary196 = if_else(X196.0609_2.81 > 10000000,1,0))
View(MBX_154196_bind)

binary154196_check <- MBX_154196_bind %>%
  # Group by Participant.ID
  group_by(Participant.ID) %>%
  
  mutate(change_in_binary5asa = sum(abs(diff(binary154)))) %>%
  
  mutate(change_in_binary196 = sum(abs(diff(binary196)))) %>%
  
  # Ungroup to operate on the full data again
  ungroup() %>%
  
  # Filter rows for participants with changes in any5asa
  filter(change_in_binary5asa > 0 & change_in_binary196 > 0) %>%
  
  # Select required columns
  select(Participant.ID, visit_num, binary154, binary196, External.ID, any5asa, oral5asa_loose, oral5asa_strict, diagnosis, X154.0502_3.83, X196.0609_2.81) %>%
  
  # Remove duplicate rows
  distinct()

View(binary154196_check)
##########################################################################################################################################################
length(unique(binary154196_check$Participant.ID))
length(unique(binary154_check$Participant.ID))
names()
data_5asa$SampleID
names(meta_data)
which(MBX_results$Metabolite == "181.0616_6.15")
cal_direct_range(formula = "C8H8N2O3", ion_charge_mode = "pos", ppm = 30, popular = T)


MBX_154196181 <- MBX_results %>% select(-"Pooled.QC.sample.CV",-"m.z",-"HMDB...Representative.ID.",-"Method",-"RT",-"Compound") %>%
  filter(Metabolite == "154.0502_3.83" | Metabolite == "196.0609_2.81" | Metabolite == "nicotinuric acid")

MBX_154196181_rotated = setNames(data.frame(t(MBX_154196181[,-1])), MBX_154196181[,1])
MBX_154196181_rotated <- cbind(Metabolite=row.names(MBX_154196181_rotated),MBX_154196181_rotated)
row.names(MBX_154196181_rotated) <- NULL

MBX_154196181_rotated <- MBX_154196181_rotated %>% mutate_all(~replace(., is.na(.), 0))
MBX_154196181_rotated <- MBX_154196181_rotated %>% rename(External.ID = Metabolite)
MBX_154196181_rotated <- MBX_154196181_rotated %>% rename(X154.0502_3.83 = "154.0502_3.83") %>% 
  rename(X196.0609_2.81 = "196.0609_2.81") %>% 
  rename(nicotinuric.acid = "nicotinuric acid")### 纯数字名字导致很多问题，此外没加.的话，空格很多函数无法识别
View(MBX_154196181_bind)
MBX_154196181_bind <- inner_join(meta_data_IBD,MBX_154196181_rotated,by= "External.ID") %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% 
  mutate(binary196 = if_else(X196.0609_2.81 > 10000000,1,0)) %>%
  mutate(binary181 = if_else(nicotinuric.acid > 10000000,1,0))

binary154196181_check <- MBX_154196181_bind %>%
  # Group by Participant.ID
  group_by(Participant.ID) %>%
  
  mutate(change_in_binary5asa = sum(abs(diff(binary154)))) %>%
  
  mutate(change_in_binary196 = sum(abs(diff(binary196)))) %>%
  
  mutate(change_in_binary181 = sum(abs(diff(binary181)))) %>%
  
  # Ungroup to operate on the full data again
  ungroup() %>%
  
  # Filter rows for participants with changes in any5asa
  filter(change_in_binary5asa > 0 & change_in_binary196 > 0) %>% #& change_in_binary181 > 0) %>%
  
  # Select required columns
  select(Participant.ID, visit_num, binary154, binary196, binary181, External.ID, any5asa, oral5asa_loose, oral5asa_strict, diagnosis, X154.0502_3.83, X196.0609_2.81, nicotinuric.acid) %>%
  
  # Remove duplicate rows
  distinct()

which(MBX_results$Metabolite == "nicotinuric acid")
which(MBX_results$Metabolite == "nicotinate")
View(binary154196181_check)
names(binary154196181_check)
binary154196181_check$binary154[1:10]
binary154196181_check$binary196
binary154196181_check$binary181
binary154196181_check$Participant.ID[1:10]

binary154196181_only <- binary154196181_check %>%
  group_by(Participant.ID) %>%
  mutate(prev_binary154 = lag(binary154, 1),
         prev_binary196 = lag(binary196, 1),
         prev_binary181 = lag(binary181, 1)) %>%
  filter((prev_binary154 == 0 & prev_binary196 == 0 & prev_binary181 == 0 &
            binary154 == 1 & binary196 == 1 & binary181 == 1) |
           (binary154 == 0 & binary196 == 0 & binary181 == 0 &
              lead(binary154) == 1 & lead(binary196) == 1 & lead(binary181) == 1)) %>%
  ungroup()

View(binary154196181_only)
unique(binary154196181_only$Participant.ID)
unique(data_5asa$Participant.ID)
length(unique(binary154196_check$Participant.ID))
length(unique(binary154196181_only$Participant.ID))
which(MBX_results$m.z == 154.0496)
MBX_results$RT[which(MBX_results$m.z == 154.0496)]
which(MBX_results$Metabolite == "154.0496_6.55")

## the "5-ASA"

MBX_1540496 <- MBX_results %>% select(-"Pooled.QC.sample.CV",-"m.z",-"HMDB...Representative.ID.",-"Method",-"RT",-"Compound") %>%
  filter(Metabolite == "154.0496_6.55")

MBX_1540496_rotated = setNames(data.frame(t(MBX_1540496[,-1])), MBX_1540496[,1])
MBX_1540496_rotated <- cbind(Metabolite=row.names(MBX_1540496_rotated),MBX_1540496_rotated)
row.names(MBX_1540496_rotated) <- NULL

MBX_1540496_rotated <- MBX_1540496_rotated %>% mutate_all(~replace(., is.na(.), 0))
MBX_1540496_rotated <- MBX_1540496_rotated %>% rename(External.ID = Metabolite)
MBX_1540496_rotated <- MBX_1540496_rotated %>% rename(X154.0496_6.55 = "154.0496_6.55")### 纯数字名字导致很多问题，此外没加.的话，空格很多函数无法识别
View(MBX_1540496_bind)
MBX_1540496_bind <- inner_join(meta_data_IBD,MBX_1540496_rotated,by= "External.ID") %>% mutate(binary154 = if_else(X154.0496_6.55 == 0,0,1))

binary1540496_check <- MBX_1540496_bind %>%
  # Group by Participant.ID
  group_by(Participant.ID) %>%
  
  mutate(change_in_binary154 = sum(abs(diff(binary154)))) %>%
  
  # Ungroup to operate on the full data again
  ungroup() %>%
  
  # Filter rows for participants with changes in any5asa
  filter(change_in_binary154 > 0) %>% #& change_in_binary181 > 0) %>%
  
  # Select required columns
  select(Participant.ID, visit_num, binary154, External.ID, any5asa, oral5asa_loose, oral5asa_strict, diagnosis,X154.0496_6.55) %>%
  
  # Remove duplicate rows
  distinct()

View(binary1540496_check)
unique(binary154_check$Participant.ID)
unique(binary1540496_check$Participant.ID)
unique(data_5asa$Participant.ID)
all(unique(data_5asa$Participant.ID) %in% unique(binary1540496_check$Participant.ID))
all(unique(data_5asa$Participant.ID) %in% unique(binary154_check$Participant.ID))


##### Use of the workflow from the HMP2 paper to research our data
###作用是使用6-MP作为标定，使用一百万丰度作为阈值，尝试找到符合我们的定义的化合物
which(MBX_results$Metabolite == "149.9994_4.86")

MBX_6mp <- MBX_results %>% select(-"Pooled.QC.sample.CV",-"m.z",-"HMDB...Representative.ID.",-"Method",-"RT",-"Compound") %>%
  filter(Metabolite == "149.9994_4.86")

MBX_6mp_rotated = setNames(data.frame(t(MBX_6mp[,-1])), MBX_6mp[,1])
MBX_6mp_rotated <- cbind(Metabolite=row.names(MBX_6mp_rotated),MBX_6mp_rotated)
row.names(MBX_6mp_rotated) <- NULL

MBX_6mp_rotated <- MBX_6mp_rotated %>% mutate_all(~replace(., is.na(.), 0))
MBX_6mp_rotated <- MBX_6mp_rotated %>% rename(External.ID = Metabolite)
MBX_6mp_rotated <- MBX_6mp_rotated %>% rename(X149.9994_4.86 = "149.9994_4.86")### 纯数字名字导致很多问题，此外没加.的话，空格很多函数无法识别

MBX_6mp_bind <- inner_join(meta_data_IBD,MBX_6mp_rotated,by= "External.ID") %>% mutate(binary6mp = if_else(X149.9994_4.86 == 0,0,1))
View(MBX_6mp_bind)

binary6mp_check <- MBX_6mp_bind %>%
  # Group by Participant.ID
  group_by(Participant.ID) %>%
  
  mutate(change_in_binary6mp = sum(abs(diff(binary6mp)))) %>%
  
  # Ungroup to operate on the full data again
  ungroup() %>%
  
  # Filter rows for participants with changes in any5asa
  filter(change_in_binary6mp > 0) %>% #& change_in_binary181 > 0) %>%
  
  # Select required columns
  select(Participant.ID, visit_num, binary6mp, External.ID, aza6mp, diagnosis,X149.9994_4.86) %>%
  
  # Remove duplicate rows
  distinct()

View(binary6mp_check)
unique(binary154_check$Participant.ID)
unique(binary1540496_check$Participant.ID)
unique(data_5asa$Participant.ID)
all(unique(data_5asa$Participant.ID) %in% unique(binary1540496_check$Participant.ID))
all(unique(data_5asa$Participant.ID) %in% unique(binary154_check$Participant.ID))