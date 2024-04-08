library(stringr)
metadata_cohort1 <- read.csv("Metadata/metadata_cohort1.csv",stringsAsFactors = FALSE, header = T)
View(metadata_cohort1)
which(colnames(metadata_cohort1) == "histologic_remission")
metadata_cohort1$id
metadata_cohort1$id <- paste0("11549.", metadata_cohort1$id)
metadata_cohort1$age_diagnosis
metadata_cohort1$age
metadata_cohort1$histologic_remission
metadata_cohort1$histologic_remission <- ifelse(metadata_cohort1$histologic_remission == 1, "remission", 
                                                ifelse(metadata_cohort1$histologic_remission == 0, NA, NA))
metadata_cohort1$biologic_exposure
metadata_cohort1$height
metadata_cohort1$current_IM
metadataC1 <- data.frame(sample_id = metadata_cohort1$id, patient_id = paste0("P", seq_len(nrow(metadata_cohort1))), 
                         disease_group = metadata_cohort1$age, 
                         disease_state = metadata_cohort1$histologic_remission,
                         use_antibiotic = rep("Missing", nrow(metadata_cohort1)), gender = metadata_cohort1$height, 
                         age = metadata_cohort1$age_diagnosis, 
                         IM_type = metadata_cohort1$current_IM, cohort_id = rep("cohort1", nrow(metadata_cohort1)))

metadata_cohort2 <- read.csv("Metadata/metadata_cohort2.csv",stringsAsFactors = FALSE, header = T)
View(metadata_cohort2)
nrow(metadata_cohort2)
colnames(metadata_cohort2)
metadata_cohort2$tube_id
metadata_cohort2$tube_id <- paste0("12675.", metadata_cohort2$tube_id)
metadata_cohort2$Diagnosis
metadata_cohort2 <- metadata_cohort2 %>%
  mutate(Diagnosis = ifelse(Diagnosis == "Healthy_control", "Healthy", Diagnosis))
metadata_cohort2$Diagnosis[205] <- NA
metadata_cohort2$Antibiotics
metadata_cohort2$sex
metadata_cohort2 <- metadata_cohort2 %>%
  mutate(sex = case_when(
    sex == "male" ~ "M",
    sex == "female" ~ "F",
    TRUE ~ NA_character_
  ))
metadata_cohort2$host_age2
metadata_cohort2$host_age2 <- as.numeric(metadata_cohort2$host_age)
metadata_cohort2$IM_Exposure_Type_1AZA_26MP_3MTX
replacements <- c(
  "1" = "AZA",
  "2" = "6MP",
  "3" = "MTX"
)
replacements <- c("0" = "not applicable",
                  "Missing" = "not applicable")
metadata_cohort2 <- metadata_cohort2 %>%
  mutate(IM_Exposure_Type_1AZA_26MP_3MTX = str_replace_all(IM_Exposure_Type_1AZA_26MP_3MTX, replacements))


metadataC2 <- data.frame(sample_id = metadata_cohort2$tube_id, patient_id = paste0("P", 40+seq_len(nrow(metadata_cohort2))),
                         disease_group = metadata_cohort2$Diagnosis, 
                         disease_state = rep(NA, nrow(metadata_cohort2)), 
                         use_antibiotic = metadata_cohort2$Antibiotics, gender = metadata_cohort2$sex,
                         age = metadata_cohort2$host_age2, 
                         IM_type = metadata_cohort2$IM_Exposure_Type_1AZA_26MP_3MTX,
                         cohort_id = rep("cohort2", nrow(metadata_cohort2)))

colnames(metadataC2)
metadata_robert <-rbind(metadataC1, metadataC2)
metadata_robert$use_antibiotic
metadata_robert$use_antibiotic <- ifelse(metadata_robert$use_antibiotic %in% c("Not applicable", "Missing"), 0, 1)
metadata_robert
write.csv(metadata_robert, "outputs/metadata_robert.csv", row.names = FALSE)
metadataC2$sample_id

metadataC2[which(metadataC2$sample_id == "12675.107220"),]
missingMETA <- c("12675.107220", "12675.FIT010E3CAL", "12675.FIT019E3CAL.SWAB", "12675.FIT037CCAL", "12675.FIT191E3CAL.A",
                 "12675.FIT191E3CAL.B", "12675.FIT238E3CAL.SWAB", "12675.FIT256E3CAL", "12675.FIT257E3CAL", "12675.FIT271E3CAL",
                 "12675.FIT273E3CAL", "12675.FIT276E3CAL", "12675.FIT278E3CAL", "12675.FIT294E3CAL")
which(missingMETA %in% metadataC2$sample_id)
length(missingMETA)
missingMETA[which(missingMETA %in% metadataC2$sample_id)]
metadata_robert$sample_id == "12675.107220"
