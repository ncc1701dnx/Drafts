
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



