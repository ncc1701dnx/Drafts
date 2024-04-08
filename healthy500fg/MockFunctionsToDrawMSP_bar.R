library(dplyr)
library(tidyr)
library(ggplot2)
library(rcartocolor)
library(forcats) ## Force Lable in pie/bar plot in order

## metadata
meta_healthy500fg <- read.table("Healthy_500FG/Healthy_500FG.metadata.csv", sep = ";", header = T, check.names = F)

## dehydrogenae
GENE_msp <- read.table("Healthy_500FG/dehydrogenae/P_vulgatus_3_oxo_5_alpha_steroid_4_dehydrogenae.faa.Healthy_500FG.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("Healthy_500FG/dehydrogenae/P_vulgatus_3_oxo_5_alpha_steroid_4_dehydrogenae.faa.Healthy_500FG.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)


CheckMergeMSPrpkm <- function(GENE_msp, GENE_rpkm, META, metrics){
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
  
  # 2. 过滤
  # 2.1. filter NA values
  filtered_data_GENE <- combined_data_GENE %>%
    filter(!is.na(genus))
  
  selected_META <- META %>%
    select(sample_id, metrics)
  
  filtered_data_GENE <- filtered_data_GENE %>%
    left_join(selected_META, by = c("sample_id" = "sample_id"))
  
  filtered_data_GENE <- filtered_data_GENE %>%
    filter(!is.na(metrics))
  filtered_data_GENE
}

filtered_data_GENE <- CheckMergeMSPrpkm(GENE_msp, GENE_rpkm, meta_healthy500fg, "gender")

## If with multiple variables -> merge into two variables
#filtered_data_GENE <- filtered_data_GENE %>%
#  mutate(filtered_data_GENE, IM_type = ifelse(IM_type %in% c("6MP", "6MP, MTX", "AZA", "AZA, MTX"), 1, 0))


GENE_M <- filtered_data_GENE %>%
  filter(gender == "M")

GENE_F <- filtered_data_GENE %>%
  filter(gender == "F")

dim(GENE_M)

dim(GENE_F)

## Male
pie_plot_matrix <- GENE_M
input <- "M"
## Female
pie_plot_matrix <- GENE_F
input <- "F"

if(T){
  # 3.1. 计算每一个genus的rpkm值的总和
  genus_rpkm_sum <- pie_plot_matrix %>%
    group_by(genus) %>%
    summarise(total_rpkm = sum(rpkm_value))
  # 3.2 assign value based on input file
  if(input == "M"){
    genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(GENE_M)
  } else if(input == "F"){
    genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(GENE_F)
  }
  
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
}

# !!! Revise Genus Name and Change to NON-CODE name
top_7_genus$genus
print(genus_rpkm_sum)

## dehydrogenae
genus_rpkm_sum$genus[which(genus_rpkm_sum$genus == "CAG-485")] <- "Sangeribacter"
## assign color
custom_colors <- c(
  "Alistipes" = "#A16928", 
  "Bacteroides" = "#C5A06B", 
  #"Enterocloster" = "#66a182",
  #"Escherichia" = "#DADEB5",
  "Phocaeicola" = "#B5C8B8",
  #"Prevotella" = "#20793c",
  "Parabacteroides" = "#5E9CA8",
  "Tidjanibacter" = "#2887A1",
  "Sangeribacter" = "brown4",
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


# 4.3. Stacked Bar Plot

ggplot(genus_rpkm_sum, aes(x = "", y = total_rpkm, fill = genus)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = custom_colors, drop=FALSE) + 
  #coord_polar(theta = "y") +  # 转换为堆叠的圆环图(if percentage)，移除这行来得到常规的堆叠条形图
  labs(
    title = "Averaged Abundance Contributed by \
    Each Taxa (Female)",
    x = NULL,
    y = NULL,
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())+  # 隐藏X轴的标签
  #coord_cartesian(ylim = c(0, 0.2)) + # force Y axis the same when have large differences
  theme(text=element_text(size=24), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold"))
