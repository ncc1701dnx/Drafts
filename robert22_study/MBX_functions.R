# Function “cal_mass”
  # 对以后的维护者（很可能就是我，嗯）
  # 本函数的主题构思即为：一个formula一般是原子开头，然后是原子个数，因此可以对第i个元素直接判断, i如果是数字，那么i-1和i+1或许是数字或许是元素，
  # 那么以数字为判断标准，从数字开始前面一定有一个元素，把元素到下一个元素之间的判断出来就是数字，然后一个数字和一个元素组成的对就是目前的计算格
 
cal_mass <- function(formula) {
  # 设置R输出的小数点，不然默认是仅有七位有效数字
  options(digits = 10)
  # 输入的formula必须保证1：全大写，Na和Cl也必须首字母大写。2：每一个元素后面必须带数字
  # 设置monoisotopic mass表，如有需要后续按顺序加入
  elements <- c('C', 'H', 'O', 
                'N', 'P', 'S', 
                'Na', 'Mg', 'F', 
                'K', 'Ca', 'Cl', 
                'I', 'Br')
  masses <- c(12.0000000, 1.0078250, 15.9949146, 
              14.0030740, 30.9737615, 31.9720710, 
              22.9897693, 23.9850417, 18.9984032, 
              38.9637069, 39.9625909, 34.9688527, 
              126.9044727, 79.9040014)
  names(masses) <- elements
  
  total_mass <- 0
  current_element <- ""
  current_count <- ""
  
  for (i in 1:nchar(formula)) {
    char <- substr(formula, i, i) # 用[i]对于全string无用 
    # 用正则表达式regex \\d for digit, [A-Z] [a-z] for letters
    # 必须使用regex否则无法直接对输入string进行操作
    # 不能用letter作为判定，而使用数字digit的原因是，对于Na，Cl这种用letter无法判定
    # 一个例子： NaCl，或者OPNa，因此无法使用letter作为判定依据，因此可以不输入数字也就无从谈起
    if (grepl("\\d", char)) {
      current_count <- paste0(current_count, char)
      if(i == nchar(formula)){
        total_mass <- total_mass + masses[current_element] * as.numeric(current_count)
      }
    } else {
      # 对于数字进行判断，非数字跳过上一个循环，即为element
      if (current_element != "" && current_count != "") {
        # 对上个循环的element和digit进行计算
        total_mass <- total_mass + masses[current_element] * as.numeric(current_count)
        current_count <- ""
      }
      if(grepl("[A-Z]", char)){
        # 仅仅使用[A-Z]而不是用[A-Za-z]的原因是如果可以输入小写，后续需要goupper()，Na会变成NA，
        # 而如果加入判断语句和element库进行matching，Na会被判断成N和a (╥﹏╥) 
        current_element <- char
      }else{
        # 假设是小写，那默认是Na或者Cl一类的元素的第二个字母，直接粘贴
        current_element <- paste0(current_element, char)
      }
    }
  }
  print(paste0(formula, "'s exact monoisotopic mass is ", total_mass))
  return(round(total_mass, 7))
}

# Test the function
#mstt <- cal_mass()
#   print(cal_mass("Na1Cl1"))
# print(cal_mass("C10H14N2Na2O8"))
# class(cal_mass("C10H15N5O10P2S1"))

# Function cal_adduct
  # 本函数的作用是通过给予的monoisotopoc mass（也就是cal_mass函数算的）计算加合后的质量
  # 一般来说算几个通用的就行了，popular默认是T
  # 对于2M的我们就不计算了，毕竟用的少，而且如果要算还要专门为他们写循环，挺麻烦的



# 函数会每一步进行print
cal_adduct <- function(compound_mass, ion_charge_mode, popular=TRUE) {
  # 先把网上搜到的做成表
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
  print(paste("Ion charge mode:", ion_charge_mode))
  print(paste("Popular:", popular))
  
  if (ion_charge_mode == "pos") {
    print("Using positive polarity adducts.")
    polarity_adducts <- adducts_pos
  } else if (ion_charge_mode == "neg") {
    print("Using negative polarity adducts.")
    polarity_adducts <- adducts_neg
  } else {
    stop("Invalid ion charge mode. Please input 'pos' or 'neg'.")
  }
  
  if (popular) {
    print("Using only the four most popular adducts.")
    polarity_adducts <- polarity_adducts[1:4]
  }
  
  print("Adducts being used:")
  print(polarity_adducts)
  
  adduct_masses <- compound_mass + unlist(polarity_adducts)
  print(adduct_masses)
  
  result <- data.frame(
    Adduct = names(adduct_masses),
    Mass = adduct_masses
  )
  
  return(result)
}

# Test the function
#adductmass_tt <- cal_adduct(mstt, "pos", TRUE)

# Function cal_adduct_range
  # 计算加上ppm以后的range

cal_adduct_range <- function(adduct_masses, ppm) {
  # 先计算给定ppm造成的误差的值
  ppm_difference <- (ppm / 1e6) * adduct_masses$Mass
  
  # 加上、减去误差
  lower_bound <- adduct_masses$Mass - ppm_difference
  upper_bound <- adduct_masses$Mass + ppm_difference
  
  mass_range <- data.frame(
    Adduct = adduct_masses$Adduct,
    LowerBound = lower_bound,
    UpperBound = upper_bound
  )
  
  return(mass_range)
}

# Test the function
#ads_tt <- cal_adduct_range(adduct_masses = adductmass_tt, ppm = 100)
# adductmass <- cal_adduct(cal_mass("C10H15N5O10P2S1"), "pos", TRUE)
# ads <- cal_adduct_range(adductmass, 200)

# 写一个可以计算速通mass range的函数
cal_direct_range <- function(formula, ion_charge_mode, ppm=60, popular = T){
  conmass <- cal_mass(formula = formula)
  addmass <- cal_adduct(compound_mass = conmass, ion_charge_mode = ion_charge_mode, popular = popular)
  addrange <- cal_adduct_range(adduct_masses = addmass, ppm = ppm)
  addrange
}
# Test the function
#cal_direct_range(formula = "C6H6N4S1", ion_charge_mode = "neg", ppm = 60)


# 最后，需要模块化自动化画图
# Function plot_fml_adducts
plot_fml_adducts <- function(rawdata, adductmasses, compound_name = NA){
  group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
  names(group_colors) <- c("Control", "AZA/6MP")
  # adductmasses需要是之前的cal_adduct_range函数算出来的data.frame
  ads <- adductmasses
  for (i in 1:nrow(ads)) {
    ad <- ads[i,]
    # Define the rt and mz range
    rtr <- c(0,840)
    mzr <- c(ad$LowerBound,ad$UpperBound)
    # Extract and check the chromatogram
    chr_raw <- chromatogram(rawdata, mz = mzr, rt = rtr)
    plot(chr_raw, col = group_colors[chr_raw$sample_group])
    title(main = paste0("\t\t\t\t\t\t\t\t\t",compound_name," ",ad$Adduct), cex.main = 1, col.main = "coral")
  }
}

# Test the function
#plot_fml_adducts(rawdata = rawdata, adductmasses = ads_tt)


### All set, now , a large function to do altogather
plot_formula <- function(rawdata ,formula, ion_mode = 'pos', ppm = 100, popular = TRUE, compound_name = NA){
  ms <- cal_mass(formula = formula)
  adductmass <- cal_adduct(compound_mass = ms, ion_charge_mode = ion_mode, popular = popular)
  ads <- cal_adduct_range(adduct_masses = adductmass, ppm = ppm)
  plot_fml_adducts(rawdata = rawdata, adductmasses = ads, compound_name = compound_name)
}

# Test this function
# examine # 6-MeMP
#plot_formula(formula = "C6H8N4S1", compound_name = "6MeMP")

# 将所有的图片转存到本地
# 这个函数的思想是，R的图片和结果是存在一个临时文档里面的
# 这个文档可以用tempdir()呼出
# plots.dir.path <- list.files(tempdir(), pattern=".png", full.names = TRUE)
# file.copy(from=plots.dir.path, to="Images/Cohort1_MBX_Targeted/")


# Function plot_volcano
# 这个函数给选定的mzr弄出来的xdata导出的result数据画volcano plot图片
# result <- cbind(as.data.frame(featureDefinitions(xdata)),featureValues(xdata, value = "into"))
# 详情请看work_0601_cohort1_mbx.R
plot_volcano <- function(result, controls, treatments){
  for(i in 1: length(colnames(result)[12:length(colnames(result))])){
    nm <- colnames(result)[12:length(colnames(result))][i]
    nm <- strsplit(nm, "_")[[1]][1]
    colnames(result)[11+i] <- nm
  }
  part_result <- result[,12:length(colnames(result))]
  control_samples <- controls$vial_name
  treatment_samples <- treatments$vial_name
  results <- data.frame()
  control_data <- c()
  treatment_data <- c()
  for (i in 1:nrow(part_result)) {
    # Subset the data for the control and treatment groups
    for (x in 1:length(part_result[i,])) {
      pct <- part_result[i,]
      ifelse(colnames(pct[x]) %in% control_samples, co <- as.numeric(pct[x]), tr <- as.numeric(pct[x]))
      control_data <- c(control_data, co)
      treatment_data <- c(treatment_data, tr)
    }
    # Performing the t-test
    ttest_result <- t.test(control_data, treatment_data, na.rm = TRUE)
    # Calculate the log2ed fold change
    log2_fold_change <- log2(mean(treatment_data, na.rm = TRUE) / mean(control_data, na.rm = TRUE))
    # Calculate non-loged fold change
    fold_change <- mean(treatment_data, na.rm = TRUE) / mean(control_data, na.rm = TRUE)
    # Store the results
    results[i, "Feature"] <- rownames(part_result)[i]
    results[i, "p.value"] <- ttest_result$p.value
    results[i, "t"] <- ttest_result$statistic
    results[i, "df"] <- ttest_result$parameter
    results[i, "confidence.low"] <- ttest_result$conf.int[1]
    results[i, "confidence.high"] <- ttest_result$conf.int[2]
    results[i, "fold_change"] <- fold_change
    results[i, "log2_fold_change"] <- log2_fold_change
  }
  results$log10_pvalue <- -log10(results$p.value)
  results$significant <- results$p.value < 0.05
  # Create a new column for significant features (p<0.01) to label
  results$label <- ifelse(results$p.value < 0.01, as.character(results$Feature), '')
  # Create a new column for the color group
  results$color_group <- ifelse(results$p.value < 0.01, "p<0.01",
                                ifelse(results$p.value < 0.05, "p<0.05", "ns"))
  # Volcano plot with color and labels
  ggplot(results, aes(x = fold_change, y = log10_pvalue, color = color_group, label = label)) +
    geom_point() +
    geom_text(aes(label = label),check_overlap = TRUE, hjust = 0,vjust = 0) +  # Add labels(Not work)
    xlab("Fold Change") +
    ylab("-log10(p-value)") +
    ggtitle("Volcano plot") +
    scale_color_manual(values = c("p<0.01" = "#990000", "p<0.05" = "#009999", "ns" = "#333333"))
}


# 写一个函数自动给出输入的CP的来源文件名
# 需要这个函数的原因是没有grep cp的来源的函数
# 如果有Feature，用featureSummary(subs_MeMP, perSampleCounts = T)
# 查询某个feature在哪些Sample中出现
# featureSummary(subs_MeMP, perSampleCounts = T)
cp_file_name <- function(xdata, cp_name){
  sample_names <- c()
  for (i in 1:length(cp_name)){
    # Get the peak matrix
    peak_data <- chromPeaks(xdata)
    # Find the row number of the peak in the peak matrix
    # Important since we do not use the subset data but rather the whole xdata as input
    peak_row <- which(rownames(peak_data) == cp_name[i])
    # Access the peak row
    sample_index <- peak_data[peak_row, "sample"]
    # Find the directory of upload file
    sample_dir <- fileNames(subs_MeMP)[sample_index]
    # Subset to the file name. Should be the last after strsplit
    sample_name <- strsplit(sample_dir, "/")[[1]][length(strsplit(sample_dir, "/")[[1]])]
    sample_names <- c(sample_names, sample_name)
  }
  sample_names
}

# Test the function
# tt1 <- cp_file_name(xdata, c("CP44411", "CP24887", "CP02509", "CP27190"))


### 两个用于根据m/z值subset经过XCMS处理的Spectra文件中的MS2数据的函数
## This function calculates the m/z range and give the filtered spectra object
cal_mz_range_filter_spectra <- function(formula, ion_charge_mode, popular = F, ppm =10){
  ads <- cal_direct_range(formula = formula,  ion_charge_mode = ion_charge_mode, popular = popular, ppm = ppm)
  filters <- logical(length(pepmasses))
  for (i in 1:nrow(ads)) {
    ad <- ads[i,]
    filters <- filters | (pepmasses >= ad$LowerBound & pepmasses <= ad$UpperBound)
  }
  ads[1,]$LowerBound
  filtered_spectra <- xdata_spectra[filters]
  filtered_spectra
}
#filtered_spectra_6MP <- cal_mz_range_filter_spectra(formula = "C5H4N4S1",  ion_charge_mode = "pos", popular = F, ppm = 10)


## This function uses above mentioned function and calculates for both pos and neg charge and write out
save_filter_spectra <- function(formula, compound_name, popular = F, ppm =10){
  fileout_dir <- "/nfs/data/metabolomics_chenrui2/Sirius/inputs/cohort1/" ## cohort1
  #fileout_dir <- "/nfs/data/metabolomics_chenrui2/Sirius/inputs/cohort2/" ## cohort2
  fileout_dir_pos <- paste0(fileout_dir, compound_name, "/pos.mgf")
  #fileout_dir_neg <- paste0(fileout_dir, compound_name, "/neg.mgf")
  ## calculate the mz range and subset the spectra object first
  filtered_spectra_pos <- cal_mz_range_filter_spectra(formula = formula,  ion_charge_mode = "pos", popular = popular, ppm = ppm)
  ##filtered_spectra_neg <- cal_mz_range_filter_spectra(formula = formula,  ion_charge_mode = "neg", popular = popular, ppm = ppm) ## Robert's data only has positive mode
  ## write the subseted spectra objects out
  export(filtered_spectra_pos, MsBackendMgf(), file = fileout_dir_pos)
  #export(filtered_spectra_neg, MsBackendMgf(), file = fileout_dir_neg)
}

