
install.packages("ggthemes")
install.packages("cowplot")
library(cowplot)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(tidyr)
setwd("../")

## read tables
df1 <- read.csv(header = T, file = "ATP Titer C16 A92 PERK 20240730/tb1.csv")
df2 <- read.csv(header = T, file = "ATP Titer C16 A92 PERK 20240730/tb2.csv")
df3 <- read.csv(header = T, file = "ATP Titer C16 A92 PERK 20240730/tb3.csv")
df4 <- read.csv(header = T, file = "ATP Titer C16 A92 PERK 20240730/tb4.csv")
df5 <- read.csv(header = T, file = "ATP Titer C16 A92 PERK 20240730/tb5.csv")
df6 <- read.csv(header = T, file = "ATP Titer C16 A92 PERK 20240730/tb6.csv")

# 计算每个批次的对照组（0浓度）的平均值
mean_control1 <- mean(df1$X0)
mean_control2 <- mean(df2$X0)
mean_control3 <- mean(df3$X0)
mean_control4 <- mean(df4$X0)
mean_control5 <- mean(df5$X0)
mean_control6 <- mean(df6$X0)

# 调整数据以消除批次效应
df1_adjusted <- sweep(df1, 2, mean_control1, "/")
df2_adjusted <- sweep(df2, 2, mean_control2, "/")
df3_adjusted <- sweep(df3, 2, mean_control3, "/")
df4_adjusted <- sweep(df4, 2, mean_control4, "/")
df5_adjusted <- sweep(df5, 2, mean_control5, "/")
df6_adjusted <- sweep(df6, 2, mean_control6, "/")

# 合并数据集
df_treatment1 <- rbind(df1_adjusted, df2_adjusted)
df_treatment2 <- rbind(df3_adjusted, df4_adjusted)
df_treatment3 <- rbind(df5_adjusted, df6_adjusted)

# 转换数据格式以便于作图和分析
df_treatment1_long <- pivot_longer(df_treatment1, cols = -X0, names_to = "concentration", values_to = "value")
df_treatment2_long <- pivot_longer(df_treatment2, cols = -X0, names_to = "concentration", values_to = "value")
df_treatment3_long <- pivot_longer(df_treatment3, cols = -X0, names_to = "concentration", values_to = "value")

# 定义作图和回归分析函数
plot_regression <- function(df, treatment_name) {
  p <- ggplot(df, aes(x = concentration, y = value)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 1), se = FALSE) +
    labs(title = paste("Linear Regression for Treatment", treatment_name))
  print(p)
  
  # 进行回归分析并返回模型摘要
  model <- lm(value ~ concentration, data = df)
  summary(model)
}

# 对每个处理执行回归分析和作图
plot_regression(df_treatment1_long, "1")
plot_regression(df_treatment2_long, "2")
plot_regression(df_treatment3_long, "3")

# 对每个处理进行ANOVA分析
anova_treatment1 <- aov(value ~ concentration, data = df_treatment1_long)
summary(anova_treatment1)

anova_treatment2 <- aov(value ~ concentration, data = df_treatment2_long)
summary(anova_treatment2)

anova_treatment3 <- aov(value ~ concentration, data = df_treatment3_long)
summary(anova_treatment3)

# 假设df_treatment1_long, df_treatment2_long, df_treatment3_long已经包含了浓度(concentration)和值(value)两列
df_treatment1_long$treatment <- 'A92'
df_treatment2_long$treatment <- 'PERK'
df_treatment3_long$treatment <- 'C16'

# 对每个数据帧的concentration列进行修改，去掉'X'
df_treatment1_long$concentration <- gsub("X10", "10", df_treatment1_long$concentration)
df_treatment2_long$concentration <- gsub("X10", "10", df_treatment2_long$concentration)
df_treatment3_long$concentration <- gsub("X10", "10", df_treatment3_long$concentration)

# 合并数据集
df_combined <- bind_rows(df_treatment1_long, df_treatment2_long, df_treatment3_long)


# ANOVA模型，包括交互作用
anova_combined <- aov(value ~ treatment * concentration, data = df_combined)
summary(anova_combined)

p_combined <- ggplot(df_combined, aes(x=concentration, y=value, color=treatment)) +
  geom_boxplot(aes(fill = treatment, alpha = 0.7)) +
  geom_point() +
  geom_line() +
  facet_wrap(~treatment) +
  labs(title = "Interaction of Treatment and Concentration",
       x = "Concentration",
       y = "Relative Luminiscence") +
  theme_wsj() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 调整X轴标签方向
        axis.title.x = element_text(face="bold", colour="#333333", size=20),  # 自定义X轴标题样式
        axis.title.y = element_text(face="bold", colour="#333333", size=20)) +  # 自定义Y轴标题样式
  scale_x_discrete(labels = function(x) sub("X10\\^", "10^", x)) # 移除X轴标签中的'X'

final_plot <- ggdraw(p_combined) +
  draw_label("TWO WAY ANOVA P Vals\nTreatment: 0.51\nConcentration: 4*10^-7\nTreat*Conce: 0.03",
             x = 0.8, y = 0.9, hjust = 1, fontface = "bold", color = "black", size = 12)

print(final_plot)

# 使用tukey HSD进行精细化分析
tukey_result1 <- TukeyHSD(anova_treatment1)
tukey_result2 <- TukeyHSD(anova_treatment2)
tukey_result3 <- TukeyHSD(anova_treatment3)
tukey_result_combinde <- TukeyHSD(anova_combined)
print(tukey_result1)


# 创建显著性矩阵的数据框
#第一组
sig_matrix <- data.frame(
  concentration1 = c("X10^-8", "X10^-8", "X10^-8", "X10^-7", "X10^-7", "X10^-6"),
  concentration2 = c("X10^-7", "X10^-6", "X10^-5", "X10^-6", "X10^-5", "X10^-5"),
  p_value = c(0.9756107, 0.9965498, 0.0980082, 0.9970214, 0.0379583, 0.0612324)
)
#第二组
sig_matrix <- data.frame(
  concentration1 = c("X10^-6", "X10^-7", "X10^-8", "X10^-7", "X10^-8", "X10^-8"),
  concentration2 = c("X10^-5", "X10^-5", "X10^-5", "X10^-6", "X10^-6", "X10^-7"),
  p_value = c(0.2133742, 0.2854482, 0.0172922, 0.9982133, 0.6780288, 0.5725159)
)
#第三组
sig_matrix <- data.frame(
  concentration1 = c("X10^-6", "X10^-7", "X10^-8", "X10^-7", "X10^-8", "X10^-8"),
  concentration2 = c("X10^-5", "X10^-5", "X10^-5", "X10^-6", "X10^-6", "X10^-7"),
  p_value = c(0.0108447, 0.0002143, 0.0404600, 0.5598169, 0.9558561, 0.2756959)
)

# 转换p值为显著性等级
sig_matrix$significance <- cut(sig_matrix$p_value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                               labels=c("<0.001", "<0.01", "<0.05", "No Signif"))

p1 <- ggplot(sig_matrix, aes(x = concentration1, y = concentration2, fill = significance)) +
  geom_tile(color = "white") +  # 使用白色分隔线
  scale_fill_manual(values = c("<0.001" = "#B95A58", "<0.01" = "#E29957", "<0.05" = "#4292C6", "No Signif" = "#4A5E65")) +
  labs(title = "Tukey HSD Test Significance Levels in A92", x = "Concentration (Mol)", y = "Concentration(Mol)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 调整X轴标签方向
  theme_wsj() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 调整X轴标签方向
        axis.title.x = element_text(face="bold", colour="#333333", size=20),  # 自定义X轴标题样式
        axis.title.y = element_text(face="bold", colour="#333333", size=20)) +  # 自定义Y轴标题样式
  scale_x_discrete(labels = function(x) sub("X10\\^", "10^", x)) +  # 移除X轴标签中的'X'
  scale_y_discrete(labels = function(x) sub("X10\\^", "10^", x))  # 移除Y轴标签中的'X'

p2 <- ggplot(sig_matrix, aes(x = concentration1, y = concentration2, fill = significance)) +
  geom_tile(color = "white") +  # 使用白色分隔线
  scale_fill_manual(values = c("<0.001" = "#B95A58", "<0.01" = "#E29957", "<0.05" = "#4292C6", "No Signif" = "#4A5E65")) +
  labs(title = "Tukey HSD Test Significance Levels in PERK", x = "Concentration (Mol)", y = "Concentration(Mol)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 调整X轴标签方向
  theme_wsj() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 调整X轴标签方向
        axis.title.x = element_text(face="bold", colour="#333333", size=20),  # 自定义X轴标题样式
        axis.title.y = element_text(face="bold", colour="#333333", size=20)) +  # 自定义Y轴标题样式
  scale_x_discrete(labels = function(x) sub("X10\\^", "10^", x)) +  # 移除X轴标签中的'X'
  scale_y_discrete(labels = function(x) sub("X10\\^", "10^", x))  # 移除Y轴标签中的'X'
p3 <- ggplot(sig_matrix, aes(x = concentration1, y = concentration2, fill = significance)) +
  geom_tile(color = "white") +  # 使用白色分隔线
  scale_fill_manual(values = c("<0.001" = "#B95A58", "<0.01" = "#E29957", "<0.05" = "#4292C6", "No Signif" = "#4A5E65")) +
  labs(title = "Tukey HSD Test Significance Levels in C16", x = "Concentration (Mol)", y = "Concentration(Mol)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 调整X轴标签方向
  theme_wsj() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 调整X轴标签方向
        axis.title.x = element_text(face="bold", colour="#333333", size=20),  # 自定义X轴标题样式
        axis.title.y = element_text(face="bold", colour="#333333", size=20)) +  # 自定义Y轴标题样式
  scale_x_discrete(labels = function(x) sub("X10\\^", "10^", x)) +  # 移除X轴标签中的'X'
  scale_y_discrete(labels = function(x) sub("X10\\^", "10^", x))  # 移除Y轴标签中的'X'


# 合并图形
# 组合四个图形
combined_plot <- plot_grid(p1, p2, p3, nrow = 2, ncol = 2, labels = "auto")
print(combined_plot)
