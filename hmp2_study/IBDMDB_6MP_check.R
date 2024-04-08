## work 14th July 2023
## This is for check 6MP and AZA in IBDMDB dataset
## Derived from IBDMDB_functions_0709.R
## 请遵照上述script的要求进行包loading和数据处理

#### 检查下controls和treatments能不能和results data对上
treatments_ID %in% colnames(MBX_results)
controls_ID %in% colnames(MBX_results)
# 全部可以对上

## 检查下treatment和control冲突的名字
which(treatments_ID %in% controls_ID)
## 没有冲突

# 将MBX文件名文件和metadata文件联系起来
## 需要在列名里面加上【1：7】而不仅仅是controls以及treatments，因为这样带元素信息（m/z，RT等）方便到时候做表
controls_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], controls_ID)]
treatments_MBX_results <- MBX_results[ ,c(colnames(MBX_results)[1:7], treatments_ID)]


### 这里新写一个函数，将NA value的赋值和others全部先做好

## 修改自IBDMDB_5ASA_check.R
## 首先是去掉可以在Treatment的不是全部都是NA的行中的NA
## 使用自己的最小值赋值
if(any(is.na(as.numeric(treatments_MBX_results[4,])))){print("yes")}
rownames(MBX_results)[1]
treatments_MBX_results["HILp_QI21194",]
for(i in 1: nrow(MBX_results)){
  metabolite <- rownames(MBX_results)[i]
  print(i)
  treatment_mm <- treatments_MBX_results[metabolite,]
  
  treatment_mm <- as.numeric(treatment_mm)
  treatment_mm <- log(treatment_mm)
  
  if(any(is.na(treatment_mm)) == T & all(is.na(treatment_mm)) == F){
    treatment_mm[which(is.na(treatment_mm))] = 0.5*min(treatment_mm, na.rm = T) 
  }
  
  treatments_MBX_results[metabolite,] <- treatment_mm
}

## 对于control同理
for(i in 1: nrow(MBX_results)){
  metabolite <- rownames(MBX_results)[i]
  print(i)
  control_mm <- controls_MBX_results[metabolite,]
  
  control_mm <- as.numeric(control_mm)
  control_mm <- log(control_mm)
  
  if(any(is.na(control_mm)) == T & all(is.na(control_mm)) == F ){
    control_mm[which(is.na(control_mm))] = 0.5*min(control_mm, na.rm = T)
  }
  
  controls_MBX_results[metabolite,] <- control_mm
}

## 然后是对于全部都是NA的那种行，我们这里用另一组的有效值来赋值
## 首先做controls的赋值，从treatments中取值
for(i in 1: nrow(MBX_results)){
  metabolite <- rownames(MBX_results)[i]
  print(i)
  control_mm <- controls_MBX_results[metabolite,]
  control_mm <- as.numeric(control_mm)
  control_mm <- log(control_mm)
  treatment_mm <- treatments_MBX_results[metabolite,]
  treatment_mm <- as.numeric(treatment_mm)
  treatment_mm <- log(treatment_mm)
  
  if(all(is.na(control_mm)) & (length(which(!is.na(treatment_mm))) >= 4)){
    control_mm = 0.5*min(treatment_mm, na.rm = T)
    }
  controls_MBX_results[metabolite,] <- control_mm
}

## 然后做treatments的赋值，从controls中取值
for(i in 1: nrow(MBX_results)){
  metabolite <- rownames(MBX_results)[i]
  print(i)
  control_mm <- controls_MBX_results[metabolite,]
  control_mm <- as.numeric(control_mm)
  control_mm <- log(control_mm)
  treatment_mm <- treatments_MBX_results[metabolite,]
  treatment_mm <- as.numeric(treatment_mm)
  treatment_mm <- log(treatment_mm)
  
  if(all(is.na(treatment_mm)) & (length(which(!is.na(control_mm))) >= 4)){
    treatment_mm = 0.5*min(control_mm, na.rm = T)
  }
  treatments_MBX_results[metabolite,] <- treatment_mm
}

## 太不容易了，写出去保存免得服务器又崩掉了
write.csv(controls_MBX_results, file = "outputs/MP_AZA_control_group_results.csv", row.names = TRUE)
write.csv(treatments_MBX_results, file = "outputs/MP_AZA_treatment_group_results.csv", row.names = TRUE)

###这里开始做t test和画图
ttest_all_results_MP_AZA <- lapply(rownames(MBX_results), function(metabolite) {
  #control_mm <- log(as.numeric(controls_MBX_results[metabolite,]))
  #treatment_mm <- log(as.numeric(treatments_MBX_results[metabolite,]))
  
  control_mm <- as.numeric(controls_MBX_results[metabolite,])
  treatment_mm <- as.numeric(treatments_MBX_results[metabolite,])
  
  control_mm <- control_mm[!is.infinite(control_mm)]
  treatment_mm <- treatment_mm[!is.infinite(treatment_mm)]
  
  control_mm <- control_mm[!is.na(control_mm)]
  treatment_mm <- treatment_mm[!is.na(treatment_mm)]
  if(length(which(!is.na(control_mm))) < 4 | length(which(!is.na(treatment_mm))) < 4) {NULL}
  else {t.test(control_mm, treatment_mm, na.action=na.omit)}
})

## 
p_values_MP_AZA <- unlist(sapply(ttest_all_results_MP_AZA, function(result){
  result$p.value
}))
p_values_MP_AZA <- p.adjust(p_values_MP_AZA, method = "BH")
mean_dif_MP_AZA <- unlist(sapply(ttest_all_results_MP_AZA, function(result){
  result$estimate[2]/result$estimate[1]
}))
ttest_and_diff_MP_AZA <- data.frame(diff = mean_dif_MP_AZA, p_val = p_values_MP_AZA)
ttest_and_diff_MP_AZA$significance <- ifelse(log2(ttest_and_diff_MP_AZA$diff) > 1 &
                                            -log10(ttest_and_diff_MP_AZA$p_val) > 2, "Significantly Abundant", "No Significance")

table(ttest_and_diff_MP_AZA$significance)

ggplot(ttest_and_diff_MP_AZA, aes(x = log2(diff), y = -log10(p_val), col = significance)) +
  geom_point() +
  scale_colour_brewer(palette = "Dark2")+
  xlab("Difference in means (Log2 Fold Change)") +
  ylab("-Log10 P-value") +
  ggtitle("AZA/6MP Volcano Plot") +
  theme(text=element_text(size=20),
        legend.position="right",
        strip.text = element_text(face = "bold"))
