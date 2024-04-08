library(dplyr)
library(tidyr)
library(ggplot2)
library(rcartocolor)
library(forcats) ## Force Lable in pie/bar plot in order

GENE_msp[which(GENE_msp$genus == "SFEL01"),]


GENE_msp <- read.table("MGX_out/HMP2/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/HMP2/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)
# 0. Read and prepare data

## GMPS
GENE_msp <- read.table("MGX_out/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)

## HPRT
GENE_msp <- read.table("MGX_out/HPRT/uniprotkb_HPRT_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/HPRT/uniprotkb_HPRT_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)

# IMPDH
GENE_msp <- read.table("MGX_out/IMPDH/uniprotkb_IMPDH_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/IMPDH/uniprotkb_IMPDH_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)

# TPMT
GENE_msp <- read.table("MGX_out/TPMT/TPMT_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/TPMT/TPMT_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)

# XDH
GENE_msp <- read.table("MGX_out/XDH/Xdh_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/XDH/Xdh_AZA6MP_related.faa.IBD_HMP2.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)

## eliminate the mspid NA rows
GENE_msp <- GENE_msp %>%
  filter(!is.na(msp_id))

## write the code to extract the genus data from msp
GENE_msp$genus <- sapply(strsplit(GENE_msp$gtdbtk_ann, " "), '[', 1)

#### 此步骤开始合并2个文件：msp,rpkm，不再单独分析
GENE_rpkm_long <- GENE_rpkm %>%
  pivot_longer(cols = -sample_id, names_to = "gc_id", values_to = "rpkm_value")

# 1. 数据整合
combined_data_GENE <- GENE_msp %>%
  #left_join(GENE_msp, by = c("sseqid" = "gc_id")) %>%
  left_join(GENE_rpkm_long, by = c("gc_id" = "gc_id"))

dim(combined_data_GENE)

# 2. 过滤
colnames(combined_data_GENE)

# 2.1. filter NA values
filtered_data_GENE <- combined_data_GENE %>%
  filter(!is.na(genus))

dim(filtered_data_GENE)

# 2.2. Change Sample Names to match with MGX metadata data
extracted_sampleID <- sapply(filtered_data_GENE$sample_id, function(x) strsplit(x, split = "_")[[1]][1])

filtered_data_GENE <- filtered_data_GENE %>%
  mutate(sampleID = extracted_sampleID)

# 2.3 bind with MGX metadata

## Treatment or Control
if(T){

# subset meta_MGXdata，因为两个矩阵都太大了，因此不要全部合并
selected_meta_MGXdata <- meta_MGXdata %>%
  select(External.ID, aza6mp)

# join two matrix togather
filtered_data_GENE <- filtered_data_GENE %>%
  left_join(selected_meta_MGXdata, by = c("sampleID" = "External.ID"))

# delete NA values
filtered_data_GENE <- filtered_data_GENE %>%
  filter(!is.na(aza6mp))
dim(filtered_data_GENE)

# delete unused columns 
colnames(filtered_data_GENE)
#filtered_data_GENE <- filtered_data_GENE[, -c(1:16)]
#dim(filtered_data_GENE)

GENE_Treatment <- filtered_data_GENE %>%
  filter(aza6mp ==1)

dim(GENE_Treatment)

GENE_Control <- filtered_data_GENE %>%
  filter(aza6mp ==0)

dim(GENE_Control)
}

## IBD or NonIBD
if(F){
  
  selected_meta_MGXdata <- meta_MGXdata %>%
    select(External.ID, diagnosis)
  
  # join two matrix togather
  filtered_data_GENE <- filtered_data_GENE %>%
    left_join(selected_meta_MGXdata, by = c("sampleID" = "External.ID"))
  
  # delete NA values
  filtered_data_GENE <- filtered_data_GENE %>%
    filter(!is.na(diagnosis))
  dim(filtered_data_GENE)
  
  # delete unused columns 
  colnames(filtered_data_GENE)
  
  
  GENE_Treatment <- filtered_data_GENE %>%
    filter(diagnosis == "CD" | diagnosis == "UC")
  
  GENE_Control <- filtered_data_GENE %>%
    filter(diagnosis == "nonIBD")
  
  dim(GENE_Treatment)
  
  dim(GENE_Control)
}

# 3. 画图前的修改和准备 
# Assign value for plot
pie_plot_matrix <- GENE_Treatment
pie_plot_matrix <- GENE_Control

# 3.1. 计算每一个genus的rpkm值的总和
genus_rpkm_sum <- pie_plot_matrix %>%
  group_by(genus) %>%
  summarise(total_rpkm = sum(rpkm_value))
genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(GENE_Treatment)
genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(GENE_Control)

# 3.3. 选取rpkm值占比最高的五个genus
top_7_genus <- genus_rpkm_sum %>%
  arrange(desc(total_rpkm)) %>%
  head(7)

# 3.3.1. 对每一个genus检查是否在top_7_genus$genus中
genus_rpkm_sum$genus <- ifelse(genus_rpkm_sum$genus %in% top_7_genus$genus, 
                               genus_rpkm_sum$genus, 
                               "others") # 如果不在，就赋值others方便后面统计

# 3.3.2. 按genus分组并计算total_rpkm的总和
genus_rpkm_sum <- genus_rpkm_sum %>%
  group_by(genus) %>%
  summarise(total_rpkm = sum(total_rpkm))

# !!! Revise Genus Name and Change to NON-CODE name
top_7_genus$genus
#genus_rpkm_sum[7,1] <- "Muribaculaceae"

# 4. 开始画 Stacked Bar Plot
# 4.1. 自定义颜色
#carto_pal(n=10, "Earth") # referrence color

## GMPS
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
#genus_rpkm_sum[7,1] <- "Muribaculaceae"


## HPRT 
custom_colors <- c(
  "Alistipes" = "#A16928", 
  "Akkermansia" = "brown4",
  #"Bifidobacterium" = "#C5A06B", 
  "Clostridium" = "#00BBFF",
  "Dysgonomonas" = "#00FFBB",
  "Enterocloster" = "#66a182",
  "Escherichia" = "#DADEB5",
  "Phocaeicola" = "#B5C8B8",
  "Prevotella" = "#20793c",
  "Roseburia" = "#5E9CA8",
  #"Muribaculaceae" = "#2887A1",
  "Ruminococcus" = "darkorange3",
  "others" = "azure4"
)
#genus_rpkm_sum[7,1] <- "Ruminococcus"
#genus_rpkm_sum[7,1] <- "Clostridium"
## IMPDH
custom_colors <- c(
  "Acidaminococcus" = "coral3",
  "Alistipes" = "#A16928", 
  #"Akkermansia" = "brown4",
  "Bacteroides" = "cornflowerblue",
  "Clostridium" = "#00BBFF",
  #"Bifidobacterium" = "#C5A06B", 
  #"Enterocloster" = "#66a182",
  "Escherichia" = "#DADEB5",
  "Phocaeicola" = "#B5C8B8",
  "Prevotella" = "#20793c",
  "Roseburia" = "#5E9CA8",
  #"Muribaculaceae" = "#2887A1",
  "Ruminococcus" = "darkorange3",
  "Sutterella" = "darkorchid4",
  "others" = "azure4"
)
#genus_rpkm_sum[7,1] <- "Ruminococcus"
#genus_rpkm_sum[6,1] <- "Ruminococcus"
#genus_rpkm_sum[7,1] <- "Clostridium"
## TPMT
custom_colors <- c(
  "Bacilli" = "burlywood", 
  "others" = "azure4"
)
#genus_rpkm_sum$genus[1] <-  "Bacilli"

genus_rpkm_sum
## Xdh
custom_colors <- c(
  #"Acidaminococcus" = "coral3",
  #"Alistipes" = "#A16928", 
  #"Akkermansia" = "brown4",
  #"Bacteroides" = "cornflowerblue",
  #"Bifidobacterium" = "#C5A06B", 
  "Coprococcus" = "#CFA127",
  "Cloacibacillus" = "#A5B7FF",
  "Enterocloster" = "#66a182",
  "Escherichia" = "#DADEB5",
  "Firmicutes" = "#03436A",
  "Hungatella" = "#D1E231",
  "Oscillospiraceae" = "#FD8CC1FF",
  "Phascolarctobacterium" = "darkolivegreen4",
  #"Phocaeicola" = "#B5C8B8",
  #"Prevotella" = "#20793c",
  "Roseburia" = "#5E9CA8",
  #"Muribaculaceae" = "#2887A1",
  "Ruminococcus" = "darkorange3",
  #"Sutterella" = "darkorchid4",
  "Scatomorpha" = "#725663FF",
  "Veillonella" = "chocolate4",
  "others" = "azure4"
)
#treatments
genus_rpkm_sum[4,1] <- "Phascolarctobacterium"
genus_rpkm_sum[6,1] <- "Ruminococcus"
genus_rpkm_sum[7,1] <- "Veillonella"
#controls
genus_rpkm_sum[1,1] <- "Oscillospiraceae"
genus_rpkm_sum[2,1] <- "Firmicutes"
genus_rpkm_sum[6,1] <- "Ruminococcus"
genus_rpkm_sum[7,1] <- "Cloacibacillus"
# 4.2. Assign colors to missing Genus
# 从 custom_colors 创建一个数据框，包含所有可能的 genus
all_genus <- data.frame(genus = names(custom_colors))

# 合并数据框
genus_rpkm_sum <- all_genus %>%
  full_join(genus_rpkm_sum, by = "genus") %>%
  replace_na(list(total_rpkm = 0))

# Use forcats package to force order labels
genus_rpkm_sum$genus <- factor(genus_rpkm_sum$genus, levels = names(custom_colors))


# 4.3. Stacked Bar Plot

ggplot(genus_rpkm_sum, aes(x = "", y = total_rpkm, fill = genus)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = custom_colors, drop=FALSE) + 
  #coord_polar(theta = "y") +  # 转换为堆叠的圆环图(if percentage)，移除这行来得到常规的堆叠条形图
  labs(
    title = "Averaged Abundance Contributed by \
    Each Taxa (XDH)",
    x = NULL,
    y = NULL,
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())+  # 隐藏X轴的标签
  #coord_cartesian(ylim = c(0, 0.2)) + # no negative Y value
  theme(text=element_text(size=24), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold"))

