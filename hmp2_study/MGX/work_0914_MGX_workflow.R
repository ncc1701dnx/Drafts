library(dplyr)
library(tidyr)
library(ggplot2)
library(rcartocolor)
library(forcats) ## Force Lable in pie/bar plot in order


############################################################## IBDMDB HMP2 #################################################################################################
## read files in R
GMPS_FG <- read.table("MGX_out/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.FG_ann", sep = "\t", header = T, check.names = F)
GMPS_blastp <- read.table("MGX_out/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out", sep = "\t", check.names = F)
GMPS_msp <- read.table("MGX_out/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GMPS_rpkm <- read.table("MGX_out/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)


## add the column name of blastp output
colnames(GMPS_blastp) <- c("qseqid", "sseqid", "pident", "length", "qlength", "slength", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

## eliminate the mspid NA rows
GMPS_msp <- GMPS_msp %>%
  filter(!is.na(msp_id))

## calculate the numbers of words in gtdbtk 
strsplit(GMPS_msp$gtdbtk_ann, split = " ")
max(sapply(GMPS_msp$gtdbtk_ann, function(x){
  cc <- strsplit(x, split = " ")
  length(cc[[1]])
}))
# 2 maximum, so genus and species most likely

## write the code to extract the genus data from msp
GMPS_msp$genus <- sapply(strsplit(GMPS_msp$gtdbtk_ann, " "), '[', 1)

#### 从这个开始合并三个文件：msp,blastp,rpkm
GMPS_rpkm_long <- GMPS_rpkm %>%
  pivot_longer(cols = -sample_id, names_to = "gc_id", values_to = "rpkm_value")

# 1. 数据整合
combined_data_GMPS <- GMPS_blastp %>%
  #left_join(GMPS_msp, by = c("sseqid" = "gc_id")) %>%
  left_join(GMPS_rpkm_long, by = c("sseqid" = "gc_id"))

combined_data_GMPS <- GMPS_msp %>%
  #left_join(GMPS_msp, by = c("sseqid" = "gc_id")) %>%
  left_join(GMPS_rpkm_long, by = c("gc_id" = "gc_id"))

dim(combined_data_GMPS)
# 2. 过滤
colnames(combined_data_GMPS)
filtered_data_GMPS <- combined_data_GMPS %>%
  filter(!is.na(genus))

colnames(filtered_data_GMPS)
extracted_sampleID <- sapply(filtered_data_GMPS$sample_id, function(x) strsplit(x, split = "_")[[1]][1])

filtered_data_GMPS <- filtered_data_GMPS %>%
  mutate(sampleID = extracted_sampleID)

# subset meta_MGXdata，因为两个矩阵都太大了，因此不要全部合并
selected_meta_MGXdata <- meta_MGXdata %>%
  select(External.ID, aza6mp)

# join two matrix togather
filtered_data_GMPS <- filtered_data_GMPS %>%
  left_join(selected_meta_MGXdata, by = c("sampleID" = "External.ID"))

# delete NA values
filtered_data_GMPS <- filtered_data_GMPS %>%
  filter(!is.na(aza6mp))
dim(filtered_data_GMPS)
# delete unused columns 
colnames(filtered_data_GMPS)
#filtered_data_GMPS <- filtered_data_GMPS[, -c(1:16)]
dim(filtered_data_GMPS)

GMPS_Treatment <- filtered_data_GMPS %>%
  filter(aza6mp ==1)

dim(GMPS_Treatment)

GMPS_Control <- filtered_data_GMPS %>%
  filter(aza6mp ==0)

dim(GMPS_Control)

# 3. 统计
sum(GMPS_Treatment$rpkm_value)/nrow(GMPS_Treatment)
sum(GMPS_Control$rpkm_value)/nrow(GMPS_Control)

pie_plot_matrix <- GMPS_Treatment
pie_plot_matrix <- GMPS_Control
# 3.1. 计算每一个genus的rpkm值的总和
genus_rpkm_sum <- pie_plot_matrix %>%
  group_by(genus) %>%
  summarise(total_rpkm = sum(rpkm_value))
genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(GMPS_Treatment)
genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(GMPS_Control)
# 3.2. 计算总的rpkm值以及每个genus的rpkm值在总rpkm中的占比
#total_rpkm <- sum(genus_rpkm_sum$total_rpkm)
#genus_rpkm_sum <- genus_rpkm_sum %>%
#  mutate(rpkm_ratio = total_rpkm / total_rpkm)

# 3.3. 选取rpkm值占比最高的五个genus
top_7_genus <- genus_rpkm_sum %>%
  arrange(desc(total_rpkm)) %>%
  head(7)

# 3.3.1. 对每一个genus检查是否在top_7_genus$genus中
genus_rpkm_sum$genus <- ifelse(genus_rpkm_sum$genus %in% top_7_genus$genus, 
                               genus_rpkm_sum$genus, 
                               "others")
# 3.3.2. 按genus分组并计算total_rpkm的总和
genus_rpkm_sum <- genus_rpkm_sum %>%
  group_by(genus) %>%
  summarise(total_rpkm = sum(total_rpkm))


genus_rpkm_sum[7,1] <- "Muribaculaceae"
top_7_genus$genus

# 4.1. 自定义颜色
#carto_pal(n=10, "Earth") # referrence color
custom_colors <- c(
  "Alistipes" = "#A16928", 
  "Akkermansia" = "brown4",
  "Bifidobacterium" = "#C5A06B", 
  "Enterocloster" = "#66a182",
  "Escherichia" = "#DADEB5",
  "Phocaeicola" = "#B5C8B8",
  "Prevotella" = "#20793c",
  "Roseburia" = "#5E9CA8",
  "Muribaculaceae" = "#2887A1",
  "others" = "azure4"
)

# 4.2. Assign colors to missing Genus
# 从 custom_colors 创建一个数据框，包含所有可能的 genus
all_genus <- data.frame(genus = names(custom_colors))

# 合并数据框
genus_rpkm_sum <- all_genus %>%
  full_join(genus_rpkm_sum, by = "genus") %>%
  replace_na(list(total_rpkm = 0))

# Use forcats package to force order labels
genus_rpkm_sum$genus <- factor(genus_rpkm_sum$genus, levels = names(custom_colors))

# ggplot2 pie chart
if(F){
pie_without_unspecified <- ggplot(genus_rpkm_sum, aes(x = "", y = total_rpkm, fill = genus)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = custom_colors) +
  #scale_fill_carto_d(name = "Genus", palette = "Earth") + #使用好看的颜色 #d to discret
  # scale_fill_carto_d(name = "Genus", palette = "Pastel") +
  labs(x = "", y = "",title = "GMPS Genus in AZA/6MP Treatments") + 
  theme_bw() +
  theme(axis.ticks = element_blank()) + #去掉上下突出的横线
  theme(axis.text.x = element_blank()) + #去掉周围的文字
  #theme(legend.title = element_blank()) + # 去掉legend的标签
  #scale_fill_discrete(breaks = genus_counts_without_unspecified$genus, labels = genusPersent, col = ) + # force show customed lables（加上百分比）
  theme(panel.grid = element_blank()) + #去掉坐标轴 
  theme(panel.border = element_blank()) + #去掉外围方框
  theme(text=element_text(size=24), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold"))

print(pie_without_unspecified)
}


### Stacked Bar Plot

ggplot(genus_rpkm_sum, aes(x = "", y = total_rpkm, fill = genus)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = custom_colors, drop=FALSE) + 
  #coord_polar(theta = "y") +  # 转换为堆叠的圆环图(if percentage)，移除这行来得到常规的堆叠条形图
  labs(
    title = "Averaged Abundance Contributed by \
    Each Taxa (GMPS Treatments)",
    x = NULL,
    y = NULL,
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())+  # 隐藏X轴的标签
  theme(text=element_text(size=24), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold"))


############################################################### Robert22 ###################################################################################################
##read files
GMPS_FG <- read.table("MGX_out/Robert/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.FG_ann", sep = "\t", header = T, check.names = F)
GMPS_blastp <- read.table("MGX_out/Robert/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out", sep = "\t", check.names = F)
GMPS_msp <- read.table("MGX_out/Robert/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GMPS_rpkm <- read.table("MGX_out/Robert/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)


GMPS_msp <- read.table("MGX_out/Robert/XDH/Xdh_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
dim(GMPS_msp)
length(which(is.na(GMPS_msp$msp_id)))
