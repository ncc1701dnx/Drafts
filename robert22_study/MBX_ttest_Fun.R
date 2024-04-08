# t test making
treatments$vial_name
controls$vial_name
colnames(result)[12:length(colnames(result))]
strsplit(colnames(result)[12:length(colnames(result))], "_")[[1]][1]

for(i in 1: length(colnames(result)[12:length(colnames(result))])){
  nm <- colnames(result)[12:length(colnames(result))][i]
  nm <- strsplit(nm, "_")[[1]][1]
  colnames(result)[11+i] <- nm
}

truenames <- c("FIT009E3CAL_GB7_01_12062.mzXML",  "FIT018E3CAL_GC12_01_12080.mzXML", "FIT026E3CAL_GC8_01_12076.mzXML" , "FIT028E3CAL_GD4_01_12085.mzXML" , "FIT042E3CAL_GC5_01_12073.mzXML" , "FIT056E3CAL_GC11_01_12079.mzXML",
               "FIT061E3CAL_GA9_01_12051.mzXML",  "FIT062E3CAL_GD3_01_12084.mzXML",  "FIT069E3CAL_GA1_01_12043.mzXML" , "FIT072E3CAL_GD2_01_12083.mzXML" , "FIT075E3CAL_GC3_01_12071.mzXML" , "FIT096E3CAL_GA6_01_12048.mzXML" ,
               "FIT100E3CAL_GC9_01_12077.mzXML" , "FIT104E3CAL_GB1_01_12056.mzXML" , "FIT105E3CAL_GA10_01_12052.mzXML", "FIT109E3CAL_GA8_01_12050.mzXML" , "FIT117E3CAL_GA11_01_12053.mzXML" ,"FIT136E3CAL_GA12_01_12054.mzXML",
               "FIT139E3CAL_GB5_01_12060.mzXML",  "FIT143E3CAL_GB12_01_12067.mzXML" ,"FIT153E3CAL_GB9_01_12064.mzXML" , "FIT154E3CAL_GA4_01_12046.mzXML" , "FIT156E3CAL_GA7_01_12049.mzXML" , "FIT159E3CAL_GD1_01_12082.mzXML" ,
               "FIT167E3CAL_GB4_01_12059.mzXML",  "FIT168E3CAL_GB3_01_12058.mzXML" , "FIT178E3CAL_GA3_01_12045.mzXML" , "FIT180E3CAL_GA5_01_12047.mzXML" , "FIT182E3CAL_GB11_01_12066.mzXML" ,"FIT183E3CAL_GC2_01_12070.mzXML" ,
               "FIT184E3CAL_GC4_01_12072.mzXML" , "FIT185E3CAL_GB8_01_12063.mzXML" , "FIT186E3CAL_GC6_01_12074.mzXML" , "FIT200E3CAL_GC1_01_12069.mzXML")
length(truenames)
for(i in 1: length(colnames(result)[12:length(colnames(result))])){
  colnames(result)[11+i] <- truenames[i]
}

result
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
  # Perform the t-test
  ttest_result <- t.test(control_data, treatment_data, na.rm = TRUE)
  # Calculate the fold change
  log2_fold_change <- mean(log2(treatment_data), na.rm = TRUE) - mean(log2(control_data), na.rm = TRUE)
  fold_change <- mean(treatment_data, na.rm = TRUE) / mean(control_data, na.rm = TRUE)
  # Store the results
  results[i, "Feature"] <- rownames(part_result)[i]
  results[i, "p.value"] <- ttest_result$p.value
  results[i, "t"] <- ttest_result$statistic
  results[i, "df"] <- ttest_result$parameter
  results[i, "confidence.low"] <- ttest_result$conf.int[1]
  results[i, "confidence.high"] <- ttest_result$conf.int[2]
  results[i, "fold_change"] <- fold_change
  results[i, "log2_fold_change" <- log2_fold_change]
}
part_result[1,control_samples[1]]
print(results)
results$log10_pvalue <- -log10(results$p.value)
results$significant <- results$p.value < 0.05
ggplot(results, aes(x = fold_change, y = log10_pvalue)) + 
  geom_point() + 
  xlab("Fold Change") + 
  ylab("-log10(p-value)") + 
  ggtitle("Volcano plot") 
ggplot(results, aes(x = log2_fold_change, y = log10_pvalue, color = significant)) + 
  geom_point() + 
  xlab("Log2 Fold Change") + 
  ylab("-log10(p-value)") + 
  ggtitle("Volcano plot") +
  scale_color_manual(values = c("black", "red"))
# Create a new column for the color group
results$color_group <- ifelse(results$p.value < 0.01, "p<0.01",
                              ifelse(results$p.value < 0.05, "p<0.05", "ns"))

# Create a new column for significant features you want to label
results$label <- ifelse(results$p.value < 0.01, as.character(results$feature), '')

# Volcano plot with color and labels
ggplot(results, aes(x = log2_fold_change, y = log10_pvalue, color = color_group, label = label)) +
  geom_point() +
  geom_text(check_overlap = TRUE, vjust = 1.5) +  # Add labels
  xlab("Log2 Fold Change") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano plot") +
  scale_color_manual(values = c("p<0.01" = "red", "p<0.05" = "blue", "ns" = "black"))

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
    # Perform the t-test
    ttest_result <- t.test(control_data, treatment_data, na.rm = TRUE)
    # Calculate the log2ed fold change
    log2_fold_change <- mean(log2(treatment_data), na.rm = TRUE) - mean(log2(control_data), na.rm = TRUE)
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
  results$label <- ifelse(results$p.value < 0.01, as.character(results$feature), '')
  # Create a new column for the color group
  results$color_group <- ifelse(results$p.value < 0.01, "p<0.01",
                                ifelse(results$p.value < 0.05, "p<0.05", "ns"))
  # Volcano plot with color and labels
  ggplot(results, aes(x = log2_fold_change, y = log10_pvalue, color = color_group, label = label)) +
    geom_point() +
    geom_text(check_overlap = TRUE, vjust = 1.5) +  # Add labels
    xlab("Log2 Fold Change") +
    ylab("-log10(p-value)") +
    ggtitle("Volcano plot") +
    scale_color_manual(values = c("p<0.01" = "cyan", "p<0.05" = "coral", "ns" = "black"))
}
plot_volcano(result = result, controls = controls,treatments = treatments)
