cal_direct_range("C7H7N1O3", ion_charge_mode = "pos", ppm = 10, popular = T)

names(MBX_results)
length(MBX_results$m.z)
cal_mz_range_filter_MBX <- function(formula, popular = F, ppm =10){
  ads <- cal_direct_range(formula = formula,  ion_charge_mode = "pos", popular = popular, ppm = ppm)
  filters <- logical(length(MBX_results$m.z))
  for (i in 1:nrow(ads)) {
    ad <- ads[i,]
    filters <- filters | (MBX_results$m.z >= ad$LowerBound & MBX_results$m.z <= ad$UpperBound)
  }
  filtered_MBX <- MBX_results[filters,]
  ads <- cal_direct_range(formula = formula,  ion_charge_mode = "neg", popular = popular, ppm = ppm)
  filters <- logical(length(MBX_results$m.z))
  for (i in 1:nrow(ads)) {
    ad <- ads[i,]
    filters <- filters | (MBX_results$m.z >= ad$LowerBound & MBX_results$m.z <= ad$UpperBound)
  }
  filtered_MBX <- rbind(filtered_MBX, MBX_results[filters,])
  filtered_MBX
}

MBX_results_5ASA <- cal_mz_range_filter_MBX("C7H7N1O3", ppm = 10, popular = T)
MBX_results_5ASA$Metabolite
MBX_results_6MP <- cal_mz_range_filter_MBX("C5H4N4S1", ppm = 10, popular = T)
class(MBX_results_6MP$Metabolite)
Metabolites_6MP <- MBX_results_6MP$Metabolite

which(MBX_results$Metabolite %in% Metabolites_6MP)

MBX_6mp <- MBX_results %>% select(-"Pooled.QC.sample.CV",-"m.z",-"HMDB...Representative.ID.",-"Method",-"RT",-"Compound") %>%
  filter(Metabolite %in% Metabolites_6MP)

MBX_6mp <- MBX_6mp %>%
  mutate(Metabolite = paste0("X", Metabolite))

names(MBX_6mp_rotated)[-1]
MBX_6mp_rotated = setNames(data.frame(t(MBX_6mp[,-1])), MBX_6mp[,1])
MBX_6mp_rotated <- cbind(Metabolite=row.names(MBX_6mp_rotated),MBX_6mp_rotated)
row.names(MBX_6mp_rotated) <- NULL

MBX_6mp_rotated <- MBX_6mp_rotated %>% mutate_all(~replace(., is.na(.), 0))
MBX_6mp_rotated <- MBX_6mp_rotated %>% rename(External.ID = Metabolite)


MBX_6mp_bind <- inner_join(meta_data_IBD,MBX_6mp_rotated,by= "External.ID")
MBX_6mp_transformed <- MBX_6mp_bind %>%
  mutate(across(starts_with("X"), 
                ~if_else(.x > 10000, 1, 0), 
                .names = "binary{.col}")
  ) %>% mutate(binary169.0177_3.2 = if_else(X169.0177_3.2 > 100000,1,0)) ## X169.0177_3.2这个和别的不一样，它的丰度更大
View(MBX_6mp_transformed)

binary6mp_check <- MBX_6mp_transformed %>%
  # Group by Participant.ID
  group_by(Participant.ID) %>%
  
  mutate(across(starts_with("binary"), 
                ~sum(abs(diff(.))), 
                .names = "diff_sum_{.col}")) %>%
  
  # Ungroup to operate on the full data again
  ungroup() %>%
  
  # Filter rows for participants with changes in any5asa
  filter(if_any(starts_with("diff_sum_"), ~ . > 0)) %>% #& change_in_binary181 > 0) %>%
  
  # Select required columns
  #select(Participant.ID, visit_num, binary6mp, External.ID, aza6mp, diagnosis,X149.9994_4.86) %>%
  select(-starts_with("X")) %>%
  select(-c(2,3,4,6,7:11))

  # Remove duplicate rows
  distinct()
View(binary6mp_check)
