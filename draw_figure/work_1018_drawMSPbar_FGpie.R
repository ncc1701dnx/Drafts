library(dplyr)
library(ggplot2)
library(RColorBrewer)

## read files
## dehydrogenae
GENE_msp <- read.table("MGX_out/P_vulgatus_Agnsk/P_vulgatus_3_oxo_5_alpha_steroid_4_dehydrogenae.faa.Healthy_500FG.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/P_vulgatus_Agnsk/P_vulgatus_3_oxo_5_alpha_steroid_4_dehydrogenae.faa.Healthy_500FG.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)
GENE_FG <- read.table("MGX_out/P_vulgatus_Agnsk/P_vulgatus_3_oxo_5_alpha_steroid_4_dehydrogenae.faa.Healthy_500FG.70cov_30id.blastp.out.FG_ann", sep = "\t", header = T, check.names = F)
## beta-glucuronidae
GENE_msp <- read.table("MGX_out/P_vulgatus_Agnsk/P_vulgatus_beta-glucuronidae.faa.Healthy_500FG.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/P_vulgatus_Agnsk/P_vulgatus_beta-glucuronidae.faa.Healthy_500FG.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)
GENE_FG <- read.table("MGX_out/P_vulgatus_Agnsk/P_vulgatus_beta-glucuronidae.faa.Healthy_500FG.70cov_30id.blastp.out.FG_ann", sep = "\t", header = T, check.names = F)


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

head(filtered_data_GENE)

pie_plot_matrix <- filtered_data_GENE

# 3.1. 计算每一个genus的rpkm值的总和
genus_rpkm_sum <- pie_plot_matrix %>%
  group_by(genus) %>%
  summarise(total_rpkm = sum(rpkm_value))
genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(filtered_data_GENE)

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
genus_rpkm_sum
## dehydrogenae
custom_colors <- c(
  "Alistipes" = "#A16928", 
  "Bacteroides" = "#2887A1",
  "CAG-485(Prevotella)" = "#C5A06B", 
  "Parabacteroides" = "#66a182",
  "Phocaeicola" = "#DADEB5",
  "Tidjanibacter" = "#20793c",
  #"Phocaeicola" = "#B5C8B8",
  #"Faecalibacterium" = "#5E9CA8",
  #"Ruthenibacterium" = "brown4",
  "others" = "azure4"
)
genus_rpkm_sum[3,1] <- "CAG-485(Prevotella)"


## beta-glucuronidae
custom_colors <- c(
  "Akkermansia" = "#B5C8B8",
  "Alistipes" = "#A16928", 
  "Bacteroides" = "#2887A1",
  #"CAG-485(Prevotella)" = "#C5A06B", 
  #"Parabacteroides" = "#66a182",
  "Phocaeicola" = "#DADEB5",
  #"Tidjanibacter" = "#20793c",
  "Prevotella" = "#5E9CA8",
  "Ruminococcus" = "darkorange3",
  "UBA7173(Porphyromonadaceae)" = "brown4",
  "others" = "azure4"
)
genus_rpkm_sum[6,1] <- "Ruminococcus"
genus_rpkm_sum[7,1] <- "UBA7173(Porphyromonadaceae)"

all_genus <- data.frame(genus = names(custom_colors))

# 合并数据框
genus_rpkm_sum_plot <- all_genus %>%
  full_join(genus_rpkm_sum, by = "genus") %>%
  replace_na(list(total_rpkm = 0))

# Use forcats package to force order labels
genus_rpkm_sum_plot$genus <- factor(genus_rpkm_sum_plot$genus, levels = names(custom_colors))


# 4.3. Stacked Bar Plot

ggplot(genus_rpkm_sum_plot, aes(x = "", y = total_rpkm, fill = genus)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = custom_colors, drop=FALSE) + 
  #coord_polar(theta = "y") +  # 转换为堆叠的圆环图(if percentage)，移除这行来得到常规的堆叠条形图
  labs(
    title = "Averaged Abundance Contributed by \
    Each Taxa (beta-glucuronidae)",
    x = NULL,
    y = NULL,
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())+  # 隐藏X轴的标签
 # coord_cartesian(ylim = c(0, 1.5)) + # no negative Y value
  theme(text=element_text(size=24), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold"))


# 1. 使用table生成频数
tab <- table(GENE_FG$FG)

# 2. 转化为数据框
data <- as.data.frame(tab)

# 3. 使用ggplot2生成堆叠条形图
p <- ggplot(data, aes(x = factor(1), y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1/sum(data$Freq))) +
  coord_polar(theta = "y") +  # 转换为堆叠的圆环图(if percentage)，移除这行来得到常规的堆叠条形图
  theme_minimal() +
  labs(title = "Functional Gene Annotation (beta-glucuronidae)",
    x = "", y = "Percentage", fill = "Category") +
  #coord_flip() +  # 翻转坐标，使之成为水平的堆叠条形图
  theme(axis.title.y=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        panel.grid=element_blank()) + 
  # coord_cartesian(ylim = c(0, 1.5)) + # no negative Y value
  theme(text=element_text(size=24), #调整字号
        legend.text = element_text(size = 14),
        legend.position="bottom",
        strip.text = element_text(face = "bold"))

# 手动添加颜色
n_colors <- length(unique(GENE_FG$FG))
color_palette <- brewer.pal(n_colors, "Set3")
p <- p + scale_fill_manual(values = color_palette)

print(p)

