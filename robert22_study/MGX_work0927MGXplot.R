
### metadata

colnames(metadata_robert)

# GMPS
GENE_msp <- read.table("MGX_out/Robert/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/Robert/GMPS/uniprotkb_GMPS_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)
dim(GENE_msp)
dim(GENE_rpkm)
## HPRT
GENE_msp <- read.table("MGX_out/Robert/HPRT/uniprotkb_HPRT_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/Robert/HPRT/uniprotkb_HPRT_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)
dim(GENE_msp)
dim(GENE_rpkm)
# IMPDH
GENE_msp <- read.table("MGX_out/Robert/IMPDH/uniprotkb_IMPDH_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/Robert/IMPDH/uniprotkb_IMPDH_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)
dim(GENE_msp)
dim(GENE_rpkm)
# TPMT
GENE_msp <- read.table("MGX_out/Robert/TPMT/TPMT_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/Robert/IMPDH/uniprotkb_IMPDH_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)
dim(GENE_msp)
dim(GENE_rpkm)
# XDH
GENE_msp <- read.table("MGX_out/Robert/XDH/Xdh_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.msp_ann", sep = "\t", header = T, check.names = F)
GENE_rpkm <- read.table("MGX_out/Robert/XDH/Xdh_AZA6MP_related.faa.IBD_RobertHM_2022.70cov_30id.blastp.out.rpkm", sep = "\t", header = T, check.names = F)
dim(GENE_msp)
dim(GENE_rpkm)


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

selected_META <- metadata_robert %>%
  select(sample_id, IM_type)

filtered_data_GENE <- filtered_data_GENE %>%
  left_join(selected_META, by = c("sample_id" = "sample_id"))

filtered_data_GENE <- filtered_data_GENE %>%
  filter(!is.na(IM_type))
dim(filtered_data_GENE)

colnames(filtered_data_GENE)
table(filtered_data_GENE$IM_type)

filtered_data_GENE <- filtered_data_GENE %>%
  mutate(filtered_data_GENE, IM_type = ifelse(IM_type %in% c("6MP", "6MP, MTX", "AZA", "AZA, MTX"), 1, 0))

#GENE_A <- filtered_data_GENE %>%
#  filter(IM_type == 1)

#GENE_B <- filtered_data_GENE %>%
#  filter(IM_type == 0)

#dim(GENE_M)

#dim(GENE_F)

GENE_Treatment <- filtered_data_GENE %>%
   filter(IM_type == 1)
  
GENE_Control <- filtered_data_GENE %>%
   filter(IM_type == 0)


t_result <- t.test(GENE_Treatment$rpkm_value, GENE_Control$rpkm_value)
p_result <- t_result$p.value
p_result
combined_data <- rbind(GENE_Treatment, GENE_Control)
combined_data$group <- factor(c(rep("Treatment",nrow(GENE_Treatment)), rep("Control", nrow(GENE_Control))))

## 画散点图和box plot
ggplot(combined_data, aes(x=group, y=rpkm_value)) +
  geom_boxplot(lwd = 1.5) + 
  #eom_jitter(position=position_jitter(0.2), size=2) +
  #geom_jitter(size=4, width = 0.1, height = 0) +
  scale_color_manual(values = my_colors) +
  labs(title="Boxplot and Scatter plot of rpkm_value", y="rpkm_value", x="")+
  theme_minimal()+
  theme(text=element_text(size=22),
        legend.position="right",
        strip.text = element_text(face = "bold"))

t_test_results <- combined_data %>%
  group_by(genus) %>%
  summarise(p_value = t.test(rpkm_value[group == "Treatment"], rpkm_value[group == "Control"])$p.value)

fold_changes <- combined_data %>%
  group_by(genus) %>%
  summarise(fold_change = mean(rpkm_value[group == "Treatment"], na.rm = TRUE) / 
              mean(rpkm_value[group == "Control"], na.rm = TRUE))

final_results <- merge(t_test_results, fold_changes, by="genus")
final_results_sorted <- final_results %>%
  arrange(p_value)

final_results_filtered <- final_results_sorted %>%
  filter(fold_change != 0)

final_results_filtered

# 写入所有结果到CSV文件
write.csv(final_results_filtered, file = "outputs/robert_genetable/XHD_all_results.csv", row.names = FALSE)

# 提取前十个结果
top_10_results <- final_results_filtered[1:10,]

# 写入前十个结果到CSV文件
write.csv(top_10_results, file = "outputs/robert_genetable/XDH_top_10_results.csv", row.names = FALSE)














pie_plot_matrix <- GENE_A
pie_plot_matrix <- GENE_B

# 3.1. 计算每一个genus的rpkm值的总和
genus_rpkm_sum <- pie_plot_matrix %>%
  group_by(genus) %>%
  summarise(total_rpkm = sum(rpkm_value))
genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(GENE_A)
genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(GENE_B)

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
## GMPS
custom_colors <- c(
  "Agathobacter" = "#A16928", 
  "Akkermansia" = "#2887A1",
  "Bacteroides" = "#C5A06B", 
  "Bifidobacterium" = "#66a182",
  "Blautia" = "#DADEB5",
  "Escherichia" = "#20793c",
  "Phocaeicola" = "#B5C8B8",
  "Faecalibacterium" = "#5E9CA8",
  "Ruthenibacterium" = "brown4",
  "others" = "azure4"
)
genus_rpkm_sum[4,1] <- "Blautia"

# HPRT
custom_colors <- c(
  "Agathobacter" = "#A16928", 
  #"Akkermansia" = "#2887A1",
  "Bacteroides" = "#C5A06B", 
  "Bifidobacterium" = "#66a182",
  #"Blautia" = "#DADEB5",
  "Dialister" = "darkorchid4",
  "Escherichia" = "#20793c",
  "Faecalibacterium" = "#5E9CA8",
  "Phocaeicola" = "#B5C8B8",
  #"Ruthenibacterium" = "brown4",
  "Roseburia" = "cornflowerblue",
  "Ruminococcus" = "darkorange3",
  "Streptococcus" = "#725663FF",
  "others" = "azure4"
)
genus_rpkm_sum[6,1] <- "Ruminococcus"

# IMPDH
custom_colors <- c(
  "Agathobacter" = "#A16928", 
  #"Akkermansia" = "#2887A1",
  "Bacteroides" = "#C5A06B", 
  "Bifidobacterium" = "#66a182",
  #"Blautia" = "#DADEB5",
  "Dialister" = "darkorchid4",
  "Escherichia" = "#20793c",
  "Faecalibacterium" = "#5E9CA8",
  "Klebsiella" = "#D1E231",
  "Phocaeicola" = "#B5C8B8",
  "Streptococcus" = "#725663FF",
  "others" = "azure4"
)

# XDH
custom_colors <- c(
  "Agathobacter" = "#A16928", 
  #"Akkermansia" = "#2887A1",
  #"Bacteroides" = "#C5A06B", 
  #"Bifidobacterium" = "#66a182",
  "Blautia" = "#DADEB5",
  #"Dialister" = "darkorchid4",
  "Enterocloster" = "#66a182",
  "Escherichia" = "#20793c",
  "Faecalibacterium" = "#5E9CA8",
  "Flavonifractor" = "darkolivegreen4",
  #"Phocaeicola" = "#B5C8B8",
  #"Ruthenibacterium" = "brown4",
  #"Roseburia" = "cornflowerblue",
  "Ruminococcus" = "darkorange3",
  #"Streptococcus" = "#725663FF",
  "Veillonella" = "chocolate4",
  "others" = "azure4"
)
#treatments
genus_rpkm_sum[2,1] <- "Blautia"
genus_rpkm_sum[7,1] <- "Ruminococcus"
#controls
genus_rpkm_sum[2,1] <- "Blautia"
genus_rpkm_sum[6,1] <- "Ruminococcus"

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
# 4.2. Assign colors to missing Genus
# 从 custom_colors 创建一个数据框，包含所有可能的 genus
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
    Each Taxa (impdh)",
    x = NULL,
    y = NULL,
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())+  # 隐藏X轴的标签
  coord_cartesian(ylim = c(0, 1.5)) + # no negative Y value
  theme(text=element_text(size=24), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold"))
