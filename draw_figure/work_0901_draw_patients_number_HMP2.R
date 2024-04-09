#### HMP2 patient distribution（find difference with paper）
# Load ggplot2
library(ggplot2)

# initialiying number of dots and samples
n_rows = 8
n_cols = ceiling(117 / n_rows)

n_shared = 8
n_abundance_based = 5
n_metadata_based = 14
n_other = 117 - n_shared - n_abundance_based - n_metadata_based

# generate the data frame
data <- data.frame(
  x = integer(0),
  y = integer(0),
  group = character(0)
)

# fill the data frame
counter = 1
for (col in 1:n_cols) {
  for (row in 1:n_rows) {
    if (counter <= n_other) {
      group = "other patients"
    } else if (counter <= n_other + n_metadata_based) {
      group = "metadata-based unique patients"
    } else if (counter <= n_other + n_metadata_based + n_shared) {
      group = "shared patients"
    } else {
      group = "abundance-based unique patients"
    }
    
    data <- rbind(data, data.frame(x = col, y = -row, group = group))
    counter = counter + 1
    
    if (counter > 117) {
      break
    }
  }
  if (counter > 117) {
    break
  }
}

my_colors <- brewer.pal(11, "Spectral")
my_colors <- my_colors[c( 9, 10, 11, 7)] # we can manuallz define the clolor settings

# Draw plots with ggplot2
ggplot(data, aes(x = x, y = y, color = group)) +
  geom_point(size = 20, shape = 15) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  labs(title = "Different ways of subset group give different results") +
  scale_color_manual(values = my_colors) +
  theme(text=element_text(size=20, face = "bold"),
        legend.position="right",
        strip.text = element_text(face = "bold"))

### cohort1

29+15
# 
n_rows = 5
n_cols = ceiling(44 / n_rows)

n_AZA_6MP_user_remission = 3
n_AZA_6MP_user = 15 - n_AZA_6MP_user_remission  # exclude those patients who are in remission
#n_UC_patient = 40 - n_AZA_6MP_user - n_AZA_6MP_user_remission
n_UC_patient = 29

# 
data <- data.frame(
  x = integer(0),
  y = integer(0),
  group = character(0)
)

# fill our data frame
counter = 1
for (row in 1:n_rows) {
  for (col in 1:n_cols) {
    if (counter <= n_UC_patient) {
      group = "UC patient"
    } else if (counter <= n_UC_patient + n_AZA_6MP_user) {
      group = "AZA/6MP user"
    } else {
      group = "AZA/6MP user got remission"
    }
    
    data <- rbind(data, data.frame(x = col, y = -row, group = group))
    counter = counter + 1
    
    if (counter > 44) {
      break
    }
  }
  if (counter > 44) {
    break
  }
}

# 从RColorBrewer包中获取颜色
my_colors <- brewer.pal(3, "Set2")

# 使用ggplot2进行绘图
ggplot(data, aes(x = x, y = y, color = group)) +
  geom_point(size = 20, shape = 15) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  labs(title = "Cohort2 Patient Distribution") +
  #scale_color_manual(values = my_colors) +
  coord_fixed(ratio = 1.5) + 
  theme(text=element_text(size=20, face = "bold"),
        legend.position="right",
        strip.text = element_text(face = "bold"))


