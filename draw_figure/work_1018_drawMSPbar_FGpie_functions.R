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

########################################################
## Step 1: visualize Taxonomy anotation based on RPKM ##
########################################################

## Step 1.1 

### This function to integrate MSP file and RPKM file
CheckMergeMSPrpkm <- function(GENE_msp, GENE_rpkm){
  ## eliminate the mspid NA rows
  GENE_msp <- GENE_msp %>%
    filter(!is.na(msp_id))
  
  ## write the code to extract the genus data from msp
  GENE_msp$genus <- sapply(strsplit(GENE_msp$gtdbtk_ann, " "), '[', 1)
  
  ############################################
  #### bind 2 files: MSP and RPKM togather ###
  GENE_rpkm_long <- GENE_rpkm %>%
    pivot_longer(cols = -sample_id, names_to = "gc_id", values_to = "rpkm_value")
  
  # 1. integrate
  combined_data_GENE <- GENE_msp %>%
    #left_join(GENE_msp, by = c("sseqid" = "gc_id")) %>%
    left_join(GENE_rpkm_long, by = c("gc_id" = "gc_id"))
  
  # 2. filter NA values
  filtered_data_GENE <- combined_data_GENE %>%
    filter(!is.na(genus))
  filtered_data_GENE
}

## Examples ##
#filtered_data_GENE <- CheckMergeMSPrpkm(GENE_msp, GENE_rpkm)
#dim(filtered_data_GENE)
#head(filtered_data_GENE)

## Step 1.2

### This function convert filtered_data_GENE to an intermediate file to manually set colors and change names 
pie_plot_matrix <- filtered_data_GENE

ConFilteredGENE <- function(pie_plot_matrix){
  # 3.1 calculate RPKM sum of each genus
  genus_rpkm_sum <- pie_plot_matrix %>%
    group_by(genus) %>%
    summarise(total_rpkm = sum(rpkm_value))
  genus_rpkm_sum$total_rpkm <- genus_rpkm_sum$total_rpkm/nrow(pie_plot_matrix)
  
  # 3.3.2. choose first 7 genus
  top_7_genus <- genus_rpkm_sum %>%
    arrange(desc(total_rpkm)) %>%
    head(7)
  
  # 3.2.1. check if a genus in top_7_genus
  genus_rpkm_sum$genus <- ifelse(genus_rpkm_sum$genus %in% top_7_genus$genus, 
                                 genus_rpkm_sum$genus, 
                                 "others") # if not, count as "others"
  
  # 3.2.2. calculate by genus and sum up the rpkm value of each genus
  genus_rpkm_sum <- genus_rpkm_sum %>%
    group_by(genus) %>%
    summarise(total_rpkm = sum(total_rpkm))
  genus_rpkm_sum
}

## Examples ##
#genus_rpkm_sum <- ConFilteredGENE(pie_plot_matrix = filtered_data_GENE)

## Step 1.3

# !!! Revise Genus Name and Change to NON-CODE name
genus_rpkm_sum

## Examples ##
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

## Step 1.4

# 1.4.1
all_genus <- data.frame(genus = names(custom_colors))

# bind data.frames
genus_rpkm_sum_plot <- all_genus %>%
  full_join(genus_rpkm_sum, by = "genus") %>%
  replace_na(list(total_rpkm = 0))

# 1.4.2
# Use forcats package to force order labels
genus_rpkm_sum_plot$genus <- factor(genus_rpkm_sum_plot$genus, levels = names(custom_colors))


# 1.4.3. Stacked Bar Plot

ggplot(genus_rpkm_sum_plot, aes(x = "", y = total_rpkm, fill = genus)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = custom_colors, drop=FALSE) + 
  #coord_polar(theta = "y") +  # change to pie plot (if percentage)
  labs(
    title = "Averaged Abundance Contributed by \
    Each Taxa (beta-glucuronidae)",
    x = NULL,
    y = NULL,
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())+  # hide x axis label
  # coord_cartesian(ylim = c(0, 1.5)) + # no negative Y value
  theme(text=element_text(size=24), 
        legend.position="right",
        strip.text = element_text(face = "bold"))

##################################################
## Step 2: visualize functional Gene Annotation ##
##################################################

DrawFGpie <- function(targetFGtable){
  # 1. use table to get occurance number
  tab <- table(targetFGtable$FG)
  # 2. convert to data.frame
  data <- as.data.frame(tab)
  
  # 3. generate pie plot using ggplot2
  p <- ggplot(data, aes(x = factor(1), y = Freq, fill = Var1)) + 
    geom_bar(stat = "identity", width = 0.5) + 
    scale_y_continuous(labels = scales::percent_format(scale = 1/sum(data$Freq))) +
    coord_polar(theta = "y") +  # change to pie plot (if percentage)
    theme_minimal() +
    labs(title = "Functional Gene Annotation",
         x = "", y = "Percentage", fill = "Category") +
    #coord_flip() +  # flip axels, get lied down stacked bar plot (conflict with theta = "y") 
    ## eliminate axis and grids
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
  
  # 4. set colors (Categories of FG annotation must not above 12)
  n_colors <- length(unique(GENE_FG$FG))
  color_palette <- brewer.pal(n_colors, "Set3") ## maximum color number: 12
  p <- p + scale_fill_manual(values = color_palette)
  
  # 5. plot!
  print(p)
}

DrawFGpie(targetFGtable = GENE_FG)
