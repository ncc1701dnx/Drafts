
#cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "neg", ppm = 60)
# 6MeMP
plot_formula(rawdata = rawdata,formula = "C6H6N4S1", compound_name = "6-MeMP")
cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,2]
cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "neg")
cal_mass(formula = "C6H6N4S1")
cal_adduct(compound_mass = 166.031317, ion_charge_mode = "pos", popular = T)
cal_direct_range("C6H6N4S1", ion_charge_mode = "pos", popular = T, ppm = 30)

# 6TUA
plot_formula(rawdata = rawdata, formula = "C5H4N4O2S1", compound_name = "6TUA")
cal_mass(formula = "C5H4N4O2S1")
cal_adduct(compound_mass = 184.0054962, ion_charge_mode = "pos", popular = T)
cal_direct_range("C5H4N4O2S1", ion_charge_mode = "pos", popular = T, ppm = 30)
# 6-MeTIMP
plot_formula(rawdata = rawdata, formula = "C11H15N4O7P1S1", compound_name = "6-MeTIMP")
# 6-TIMP
plot_formula(rawdata = rawdata, formula = "C10H13N4O7P1S1", compound_name = "6-TIMP")
cal_mass(formula = "C10H13N4O7P1S1")
cal_adduct(compound_mass = 364.0242557, ion_charge_mode = "pos", popular = T)
cal_direct_range("C10H13N4O7P1S1", ion_charge_mode = "pos", popular = T, ppm = 30)
# AZA
plot_formula(rawdata = rawdata, formula = "C9H7N7O2S1", compound_name = "AZA")
cal_mass(formula = "C9H7N7O2S1")
cal_adduct(compound_mass = 277.0381932, ion_charge_mode = "pos", popular = T)
cal_direct_range("C9H7N7O2S1", ion_charge_mode = "pos", popular = T, ppm = 30)
# 6-MP
plot_formula(rawdata = rawdata, formula = "C5H4N4S1", compound_name = "6-MP")
cal_mass(formula = "C5H4N4S1")
cal_adduct(compound_mass = 152.015667, ion_charge_mode = "pos", popular = T)
cal_direct_range("C5H4N4S1", ion_charge_mode = "pos", popular = T, ppm = 30)
# 6-TGMP
plot_formula(rawdata = rawdata, formula = "C10H14N5O7P1S1", compound_name = "6-TGMP")
# 6-TGDP
plot_formula(rawdata = rawdata, formula = "C10H15N5O10P2S1", compound_name = "6-TGDP")
# 6-TGTP
plot_formula(rawdata = rawdata, formula = "C10H16N5O13P3S1", compound_name = "6-TGTP")
cal_mass(formula = "C10H16N5O13P3S1")
cal_adduct(compound_mass = 538.9678153, ion_charge_mode = "pos", popular = T)
cal_direct_range("C10H16N5O13P3S1", ion_charge_mode = "pos", popular = T, ppm = 30)

MBX_pos_treatments[which(MBX_pos_treatments$m.z > cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,2] &
      MBX_pos_treatments$m.z < cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,3]), ]
length(which(is.na(MBX_pos_controls[which(MBX_pos_treatments$m.z > cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,2] &
                         MBX_pos_treatments$m.z < cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,3]), ])))

tt <- MBX_pos_controls[which(MBX_pos_treatments$m.z > cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,2] &
                               MBX_pos_treatments$m.z < cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,3]), ]
t1 <- MBX_pos_treatments[which(MBX_pos_treatments$m.z > cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,2] &
                                 MBX_pos_treatments$m.z < cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,3]), ]
0.5*min(t1[,-(1:7)], na.rm = T)
min(c(10,2,3,4,5,6,NA,1,4,NA), na.rm = T)
tt[1,] <- 0.5*min(t1[,-(1:7)], na.rm = T)
class(t1[,-(1:7)])

tt <- which(MBX_pos_treatments$m.z > cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,2] &
        MBX_pos_treatments$m.z < cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "pos")[1,3])
t_treatment <- length(which(is.na(MBX_pos_treatments[tt[1], -(1:7)])))
t_control <- length(which(is.na(MBX_pos_controls[tt[3], -(1:7)])))
dim(MBX_pos_controls)

t_treatment <- c()
t_control <- c()
for(i in 1:length(tt)){
  tt_t <- length(which(is.na(MBX_pos_treatments[tt[i], -(1:7)])))
  t_treatment <- c(t_treatment,tt_t)
  tt_c <- length(which(is.na(MBX_pos_controls[tt[i], -(1:7)])))
  t_control <- c(t_control,tt_c)
}

dim(MBX_pos_treatments)
dim(MBX_pos_controls)


#### 开始写函数
MBX_C18neg_treatments <- subset(treatments_MBX_results, Method == "C18-neg" )
MBX_HILneg_treatments <- subset(treatments_MBX_results, Method == "HILIC-neg")
MBX_C18pos_treatments <- subset(treatments_MBX_results, Method == "C8-pos")
MBX_HILpos_treatments <- subset(treatments_MBX_results, Method == "HILIC-pos")
table(controls_MBX_results$Method)
MBX_C18neg_controls <- subset(controls_MBX_results, Method == "C18-neg")
MBX_HILneg_controls <- subset(controls_MBX_results, Method == "HILIC-neg")
MBX_C18pos_controls <- subset(controls_MBX_results, Method == "C8-pos")
MBX_HILpos_controls <- subset(controls_MBX_results, Method == "HILIC-pos")

ion_charge <- c()
adduct <- c()
m_z <- c()
RT <- c()
t_treatment <- c()
t_control <- c()
mean_treatment <- c()
ion_charge_mode <- "neg"
com_range <- cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = ion_charge_mode)
for (x in 1:nrow(com_range)) {
  tt_C18neg <- which(MBX_C18neg_treatments$m.z > com_range[x,2] &
                       MBX_C18neg_treatments$m.z < com_range[x,3])
  for(i in 1:length(tt_C18neg)){
    tt_t_C18neg <- length(which(is.na(MBX_C18neg_treatments[tt_C18neg[i], -(1:7)])))
    t_treatment <- c(t_treatment,tt_t_C18neg)
    tt_c_C18neg <- length(which(is.na(MBX_C18neg_controls[tt_C18neg[i], -(1:7)])))
    t_control <- c(t_control,tt_c_C18neg)
    adduct <- c(adduct, names(adducts_pos[x]))
    ion_charge <- c(ion_charge, MBX_C18neg_treatments$Method[tt_C18neg[i]])
    m_z <- c(m_z, MBX_C18neg_treatments$m.z[tt_C18neg[i]])
    RT <- c(RT, MBX_C18neg_treatments$RT[tt_C18neg[i]])
    mean_treatment <- c(mean_treatment, mean(as.numeric(MBX_C18neg_treatments[tt_C18neg[1],-(1:7)]), na.rm = T))
  }
}

length(ion_charge)
### n
# create a list of all your datasets
ion_charge_mode <- "neg"
formula <- "C6H6N4S1"
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

# select the appropriate datasets based on ion_charge_mode
if (ion_charge_mode == "neg") {
  subset_list_C18 <- dataset_list[c("MBX_C18neg_treatments", "MBX_C18neg_controls")]
  subset_list_HIL <- dataset_list[c("MBX_HILneg_treatments", "MBX_HILneg_controls")]
} else if (ion_charge_mode == "pos") {
  subset_list_C18 <- dataset_list[c("MBX_C18pos_treatments", "MBX_C18pos_controls")]
  subset_list_HIL <- dataset_list[c("MBX_HILpos_treatments", "MBX_HILpos_controls")]
  sub_list <- list("subset_list_C18" = subset_list_C18,
                   "subset_list_HIL" = subset_list_HIL)
}

# # define your vectors
# ion_charge <- c()
# adduct <- c()
# m_z <- c()
# RT <- c()
# compound <- c()
# t_treatment <- c()
# t_control <- c()
# mean_treatment <- c()
# 
# length(subset_list_C18)
# 
# for(n in 1:2){
#   n =1
#   if (n ==1) {
#     subset_list <- subset_list_C18
#     names(subset_list) <- names(subset_list_C18)
#   } else {
#     subset_list <- subset_list_HIL
#     names(subset_list) <- subset_list_HIL
#   }
#   for (subset_name in names(subset_list)) {
#     subset_name = names(subset_list)[1]
#     subset <- subset_list[[subset_name]]
#     com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode)
#     
#     for (x in 1:nrow(com_range)) {
#       x =1
#       tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
#       
#       for(i in 1:length(tt)){
#         i = 1
#         tt_t <- length(which(is.na(subset[tt[i], -(1:7)])))
#         tt_c <- length(which(is.na(subset[tt[i], -(1:7)])))
#         compound <- c(compound, formula)
#         adduct <- c(adduct, names(adducts_pos[x]))
#         ion_charge <- c(ion_charge, ifelse(length(subset$Method[tt[i]]) == 0, NA, subset$Method[tt[i]]))
#         m_z <- c(m_z, ifelse(length(subset$m.z[tt[i]]) == 0, NA, subset$m.z[tt[i]]))
#         RT <- c(RT, ifelse(length(subset$RT[tt[i]]) == 0, NA, subset$RT[tt[i]]))
#         mean_treatment <- c(mean_treatment, mean(as.numeric(subset[tt[1],-(1:7)]), na.rm = T))
#         
#         if(grepl("treatments", subset_name)){
#           t_treatment <- c(t_treatment, tt_t)
#         } else if(grepl("controls", subset_name)){
#           t_control <- c(t_control, tt_c)
#         }
#       }
#     }
#   }
# }

IBDMDB_MBX_check <- function(formula){
  ion_charge_mode = "neg"
  # select the appropriate datasets based on ion_charge_mode
  # if (ion_charge_mode == "neg") {
  #   subset_list_C18 <- dataset_list[c("MBX_C18neg_treatments", "MBX_C18neg_controls")]
  #   subset_list_HIL <- dataset_list[c("MBX_HILneg_treatments", "MBX_HILneg_controls")]
  # } else if (ion_charge_mode == "pos") {
  #   subset_list_C18 <- dataset_list[c("MBX_C18pos_treatments", "MBX_C18pos_controls")]
  #   subset_list_HIL <- dataset_list[c("MBX_HILpos_treatments", "MBX_HILpos_controls")]
  #   sub_list <- list("subset_list_C18" = subset_list_C18,
  #                    "subset_list_HIL" = subset_list_HIL)
  # }
  subset_list_C18 <- dataset_list[c("MBX_C18neg_treatments", "MBX_C18neg_controls")]
  subset_list_HIL <- dataset_list[c("MBX_HILneg_treatments", "MBX_HILneg_controls")]
  sub_list <- list("subset_list_C18" = subset_list_C18,
                   "subset_list_HIL" = subset_list_HIL)
  # define your vectors
  ion_charge <- c()
  adduct <- c()
  m_z <- c()
  RT <- c()
  compound <- c()
  t_treatment <- c()
  t_control <- c()
  mean_treatment <- c()
  
  
  # loop over the subsets
  
  # if (ion_charge_mode == "neg") {
  #   subset_list <- dataset_list[c("MBX_C18neg_treatments", "MBX_HILneg_treatments", "MBX_C18neg_controls", "MBX_HILneg_controls")]
  # } else {
  #   subset_list <- dataset_list[c("MBX_C18pos_treatments", "MBX_HILpos_treatments", "MBX_C18pos_controls", "MBX_HILpos_controls")]
  # }
  subset_list <- dataset_list[c("MBX_C18neg_treatments", "MBX_HILneg_treatments", "MBX_C18neg_controls", "MBX_HILneg_controls")]
  for (subset_name in names(subset_list)) {
    subset <- subset_list[[subset_name]]
    com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode)
    
    for (x in 1:nrow(com_range)) {
      tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
      
      for(i in 1:length(tt)){
        tt_t <- length(which(is.na(subset[tt[i], -(1:7)])))
        tt_c <- length(which(is.na(subset[tt[i], -(1:7)])))
        compound <- c(compound, formula)
        adduct <- c(adduct, names(adducts_pos[x]))
        ion_charge <- c(ion_charge, ifelse(length(subset$Method[tt[i]]) == 0, NA, subset$Method[tt[i]]))
        m_z <- c(m_z, ifelse(length(subset$m.z[tt[i]]) == 0, NA, subset$m.z[tt[i]]))
        RT <- c(RT, ifelse(length(subset$RT[tt[i]]) == 0, NA, subset$RT[tt[i]]))
        mean_treatment <- c(mean_treatment, mean(as.numeric(subset[tt[1],-(1:7)]), na.rm = T))
        
        if(grepl("treatments", subset_name)){
          t_treatment <- c(t_treatment, tt_t)
        } else if(grepl("controls", subset_name)){
          t_control <- c(t_control, tt_c)
        }
      }
    }
  }
  
  ion_charge_mode = "pos"
  # select the appropriate datasets based on ion_charge_mode
  # if (ion_charge_mode == "neg") {
  #   subset_list_C18 <- dataset_list[c("MBX_C18neg_treatments", "MBX_C18neg_controls")]
  #   subset_list_HIL <- dataset_list[c("MBX_HILneg_treatments", "MBX_HILneg_controls")]
  # } else if (ion_charge_mode == "pos") {
  #   subset_list_C18 <- dataset_list[c("MBX_C18pos_treatments", "MBX_C18pos_controls")]
  #   subset_list_HIL <- dataset_list[c("MBX_HILpos_treatments", "MBX_HILpos_controls")]
  #   sub_list <- list("subset_list_C18" = subset_list_C18,
  #                    "subset_list_HIL" = subset_list_HIL)
  # }
  subset_list_C18 <- dataset_list[c("MBX_C18pos_treatments", "MBX_C18pos_controls")]
  subset_list_HIL <- dataset_list[c("MBX_HILpos_treatments", "MBX_HILpos_controls")]
  sub_list <- list("subset_list_C18" = subset_list_C18,
                   "subset_list_HIL" = subset_list_HIL)
  # loop over the subsets
  # 
  # if (ion_charge_mode == "neg") {
  #   subset_list <- dataset_list[c("MBX_C18neg_treatments", "MBX_HILneg_treatments", "MBX_C18neg_controls", "MBX_HILneg_controls")]
  # } else {
  #   subset_list <- dataset_list[c("MBX_C18pos_treatments", "MBX_HILpos_treatments", "MBX_C18pos_controls", "MBX_HILpos_controls")]
  # }
  subset_list <- dataset_list[c("MBX_C18pos_treatments", "MBX_HILpos_treatments", "MBX_C18pos_controls", "MBX_HILpos_controls")]
  for (subset_name in names(subset_list)) {
    subset <- subset_list[[subset_name]]
    com_range <- cal_direct_range(formula = formula, ion_charge_mode = ion_charge_mode)
    
    for (x in 1:nrow(com_range)) {
      tt <- which(subset$m.z > com_range[x,2] & subset$m.z < com_range[x,3])
      
      for(i in 1:length(tt)){
        tt_t <- length(which(is.na(subset[tt[i], -(1:7)])))
        tt_c <- length(which(is.na(subset[tt[i], -(1:7)])))
        compound <- c(compound, formula)
        adduct <- c(adduct, names(adducts_pos[x]))
        ion_charge <- c(ion_charge, ifelse(length(subset$Method[tt[i]]) == 0, NA, subset$Method[tt[i]]))
        m_z <- c(m_z, ifelse(length(subset$m.z[tt[i]]) == 0, NA, subset$m.z[tt[i]]))
        RT <- c(RT, ifelse(length(subset$RT[tt[i]]) == 0, NA, subset$RT[tt[i]]))
        mean_treatment <- c(mean_treatment, mean(as.numeric(subset[tt[1],-(1:7)]), na.rm = T))
        
        if(grepl("treatments", subset_name)){
          t_treatment <- c(t_treatment, tt_t)
        } else if(grepl("controls", subset_name)){
          t_control <- c(t_control, tt_c)
        }
      }
    }
  }
  compoundwise_IBDMDB <- data.frame(ion_charge = ion_charge, adduct = adduct, m_z = m_z, RTime = RT, formula = compound, 
                                    na_in_treatment = t_treatment, na_in_control = t_control, 
                                    mean_value_treatment = mean_treatment)
  compoundwise_IBDMDB
}




# 6MeMP
IBDMDB_6MeMP <- IBDMDB_MBX_check(formula = "C6H6N4S1")
