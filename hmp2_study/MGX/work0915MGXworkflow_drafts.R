GMPS_rpkm_long <- GMPS_rpkm %>%
  pivot_longer(cols = -sample_id, names_to = "gc_id", values_to = "rpkm_value")

# 1. 数据整合
combined_data_GMPS <- GMPS_blastp %>%
  left_join(GMPS_msp, by = c("sseqid" = "gc_id")) %>%
  left_join(GMPS_rpkm_long, by = c("sseqid" = "gc_id"))

# 2. 过滤
filtered_data_GMPS <- combined_data_GMPS %>%
  filter(!is.na(genus))

# 3. 统计
sum(filtered_data_GMPS$rpkm_value)

# 3.1. 计算每一个genus的rpkm值的总和
genus_rpkm_sum <- filtered_data_GMPS %>%
  group_by(genus) %>%
  summarise(total_rpkm = sum(rpkm_value))

# 3.2. 计算总的rpkm值以及每个genus的rpkm值在总rpkm中的占比
total_rpkm <- sum(genus_rpkm_sum$total_rpkm)
genus_rpkm_sum <- genus_rpkm_sum %>%
  mutate(rpkm_ratio = total_rpkm / total_rpkm)

# 3. 选取rpkm值占比最高的五个genus
top_5_genus <- genus_rpkm_sum %>%
  arrange(desc(rpkm_ratio)) %>%
  head(5)

# ggplot2 pie chart
pie_without_unspecified <- ggplot(top_5_genus, aes(x = "", y = total_rpkm, fill = genus)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_carto_d(name = "Genus", palette = "Earth") + #使用好看的颜色 #d to discrete
  labs(x = "", y = "",title = "Pie chart of Genus occurance in GMPS") + 
  theme_bw() +
  theme(axis.ticks = element_blank()) + #去掉上下突出的横线
  theme(axis.text.x = element_blank()) + #去掉周围的文字
  #theme(legend.title = element_blank()) + # 去掉legend的标签
  #scale_fill_discrete(breaks = genus_counts_without_unspecified$genus, labels = genusPersent, col = ) + # 自定义标签（加上百分比）
  theme(panel.grid = element_blank()) + #去掉坐标轴 
  theme(panel.border = element_blank()) + #去掉外围方框
  theme(text=element_text(size=30), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold"))

print(pie_without_unspecified)

length(unique(filtered_data_GMPS$sseqid))
length(unique(meta_data$External.ID))
which(filtered_data_GMPS$sseqid %in% meta_data$External.ID)

which(meta_data$External.ID ==strsplit(filtered_data_GMPS$sample_id[1], split = "_")[[1]][1])

length(which(filtered_data_GMPS$sampleID %in% meta_MGXdata$External.ID))
dim(filtered_data_GMPS)
unique(filtered_data_GMPS$sampleID[which(!filtered_data_GMPS$sampleID %in% meta_data$External.ID)])
unique(filtered_data_GMPS$sampleID[which(filtered_data_GMPS$sampleID %in% meta_MGXdata$External.ID)])
which(meta_data$External.ID == "CSM5FZ3R")
which(meta_data$External.ID == "CSM5FZ3V")
which(meta_data$External.ID == "CSM67UAW")
which(meta_data$External.ID == "CSM79HOT")
unique(filtered_data_GMPS$sample_id)[(unique(filtered_data_GMPS$sample_id) %in% unique(meta_data$External.ID))]

table(meta_MGXdata$any5asa)
table(meta_MGXdata$aza6mp)
meta_MGXdata$External.ID

#### 从这个开始
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

colnames(filtered_data_GMPS)
filtered_data_GMPS <- filtered_data_GMPS[, -c(1:16,21)]

GMPS_Treatment <- filtered_data_GMPS %>%
  filter(aza6mp ==1)

GMPS_Control <- filtered_data_GMPS %>%
  filter(aza6mp ==0)

# 3. 统计
sum(GMPS_Treatment$rpkm_value)
sum(GMPS_Control$rpkm_value)

pie_plot_matrix <- GMPS_Control

# 3.1. 计算每一个genus的rpkm值的总和
genus_rpkm_sum <- pie_plot_matrix %>%
  group_by(genus) %>%
  summarise(total_rpkm = sum(rpkm_value))


# 3.2. 计算总的rpkm值以及每个genus的rpkm值在总rpkm中的占比
total_rpkm <- sum(genus_rpkm_sum$total_rpkm)
genus_rpkm_sum <- genus_rpkm_sum %>%
  mutate(rpkm_ratio = total_rpkm / total_rpkm)

# 3. 选取rpkm值占比最高的五个genus
top_5_genus <- genus_rpkm_sum %>%
  arrange(desc(rpkm_ratio)) %>%
  head(7)

top_5_genus[1,1] <- "Bacteroides"

# ggplot2 pie chart
pie_without_unspecified <- ggplot(top_5_genus, aes(x = "", y = total_rpkm, fill = genus)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_carto_d(name = "Genus", palette = "Earth") + #使用好看的颜色 #d to discrete
  # scale_fill_carto_d(name = "Genus", palette = "Pastel") +
  labs(x = "", y = "",title = "GMPS Genus in AZA/6MP Controls") + 
  theme_bw() +
  theme(axis.ticks = element_blank()) + #去掉上下突出的横线
  theme(axis.text.x = element_blank()) + #去掉周围的文字
  #theme(legend.title = element_blank()) + # 去掉legend的标签
  #scale_fill_discrete(breaks = genus_counts_without_unspecified$genus, labels = genusPersent, col = ) + # 自定义标签（加上百分比）
  theme(panel.grid = element_blank()) + #去掉坐标轴 
  theme(panel.border = element_blank()) + #去掉外围方框
  theme(text=element_text(size=32), #调整字号
        legend.position="right",
        strip.text = element_text(face = "bold"))

print(pie_without_unspecified)
