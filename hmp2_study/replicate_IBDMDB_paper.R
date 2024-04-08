newids <- c("PSM7J19B", "PSM7J19J", "PSM6XBSE", "PSM6XBVM", "PSM6XBSK", "PSM6XBUG", "PSM6XBRK", "PSM7J1CU", "MSMA26AZ", "MSMB4LZ4",
            "MSM6J2IG", "MSM6J2Q3", "MSM5LLDI", "MSM5LLDS", "HSMA33OJ", "HSMA33M8", "HSMA33OZ", "HSMA33MI", "HSM5MD73", "HSM6XRS8", "HSM6XRRV", 
            "HSM6XRVO", "CSM79HRG", "CSM7KOO9", "CSM5MCXL", "CSM67UDN")
which(treatments_ID %in% newids)
meta_data$any5asa[meta_data$Participant.ID == "C3004"]
meta_data$visit_num[meta_data$Participant.ID == "C3004"]
which(MBX_results$Metabolite == "154.0502_3.83")
head(MBX_results$Metabolite,100)
names(MBX_results)
MBX_results$RT[which(MBX_results$m.z == 154.0502)]


############################# METABOLOMICS ########################################################


## metabolomics file is edited somewhat in excel because i moved the m/z ratios to the names for the unnamed compounds and got rid of prime values, was too lazy to do this in R , but verified 12/15/20 ####

metabolomics <- read.table("/home/rsm34/5asa/HMP2_metabolomics_edit.txt", header=T, sep="\t", check.names=FALSE)

#5asa 

metab_5asa <- metabolomics %>% select(-"Pooled QC sample CV",-"m/z",-"HMDB (*Representative ID)",-"Method",-"RT",-"Compound") %>%
  rename(SampleID=Metabolite) %>% 
  filter(SampleID == "carboxyibuprofen" | SampleID == "metronidazole" | SampleID == "salicylate" | SampleID == "sulfapyridine" | SampleID == "154.0502_3.83" |
           SampleID == "132.0454_1.66" | SampleID == "194.046_3.93" | SampleID == "154.0502_3.83" | SampleID == "196.0481_7.88" | SampleID == "196.0533_2.89" |
           SampleID == "196.0579_7.11" | SampleID == "196.0602_5.19" | SampleID == "196.0609_2.81" |
           SampleID == "196.0609_5.03" | SampleID == "196.061_1.91" | SampleID == "196.0612_1.69" | SampleID == "196.0615_5.22" |
           SampleID == "196.0631_4.36" | SampleID == "196.0633_4.86" | SampleID == "196.0647_3.68" | SampleID == "196.0721_7.08" |
           SampleID == "196.073_6.08" | SampleID == "196.0798_4.09" | SampleID == "196.081_3.52" | SampleID == "196.0827_4.94" |
           SampleID == "196.0833_4.8" | SampleID == "196.0837_4.57" | SampleID == "196.0837_5.06" | SampleID == "196.0956_7.06" |
           SampleID == "196.0959_7.05" | SampleID == "196.0968_8.17" | SampleID == "196.097_4.83" | SampleID == "196.0972_2.73" |
           SampleID == "196.0972_2.96" | SampleID == "196.0972_3.71" | SampleID == "196.0972_3.78" | SampleID == "196.0972_4.47" |
           SampleID == "196.0972_5.57" |  SampleID == "196.0974_5.86" | SampleID == "196.0974_7.14" | SampleID == "196.0974_7.25")

metab_rotated = setNames(data.frame(t(metab_5asa[,-1])), metab_5asa[,1])
metab_rotated<- cbind (SampleID=row.names(metab_rotated),metab_rotated)
row.names(metab_rotated) <- NULL
mesal_metab_rotated <- metab_rotated
mesal_metab_rotated$SampleID <- as.character(mesal_metab_rotated$SampleID)

write.table(mesal_metab_rotated,"mesal_metab_rotated2.txt",sep="\t",row.names=FALSE,quote=F)

#all metab rotated 

metab_all <- metabolomics %>% select(-"Pooled QC sample CV",-"m/z",-"HMDB (*Representative ID)",-"Method",-"RT",-"Compound") %>%
  rename(SampleID=Metabolite) 

metab_rotated_all = setNames(data.frame(t(metab_all[,-1])), metab_all[,1])
metab_rotated_all<- cbind (SampleID=row.names(metab_rotated_all),metab_rotated_all)
row.names(metab_rotated_all) <- NULL
metab_rotated_all$SampleID <- as.character(metab_rotated_all$SampleID)

write.table(metab_rotated_all,"metab_rotated_all.txt",sep="\t",row.names=FALSE,quote=F)

metab_rotated <- read.table("/home/rsm34/5asa/metab_rotated_all.txt",header=T,sep="\t")
metab_rotated2 <- metab_rotated %>% mutate_all(~replace(., is.na(.), 0)) #make sure NA is = 0. 
write.table(metab_rotated2 , "metab_data_nozero2.txt", sep="\t", quote=FALSE,row.names=FALSE)
            
#now link with metabolomics data, after pre-processing script made it 
metab_rotated2 <- read.table("/home/rsm34/5asabb3/archive/metabolite/metab_data_nozero2.txt",header=T,sep="\t")
metab_IBD<- inner_join(meta_data_IBD,metab_rotated2,by="SampleID") %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% relocate(binary154, .after = IBDusers)

metab_154196<- metab_rotated2 %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% relocate(binary154, .after = SampleID) %>% relocate(X154.0502_3.83, .after=binary154) %>% 
  mutate(binary196 = if_else(X196.0609_2.81 > 10000000,1,0)) %>% relocate(binary196, .after = X154.0502_3.83) %>% relocate(X196.0609_2.81, .after=binary196) %>% select(,1:5)
write.table(metab_154196, "metabs_export.txt", sep="\t", quote=FALSE,row.names=FALSE) ### to be used for other analyses as  short cut 

########################## MERGE ##################################

metab_metag <- left_join(metab_rotated,gene_rotated,by="SampleID")
metab_metag <- na.omit(metab_metag)
write.table(metab_metag,"metab_metag.txt",sep="\t",row.names=FALSE,quote=F)




######################################################################################################
# R programs for creating "metabolite" plots					#
# ============================================					#
# Purposes:  new use boxplots,    						#
#              networks, identificaiton of unannotated compounds, r2 parsing		#
######################################################################################################

#references: 
#https://yardstick.tidymodels.org/reference/roc_curve.html
#https://davetang.org/muse/2017/03/16/matrix-to-adjacency-list-in-r/
#https://www.r-graph-gallery.com/257-input-formats-for-network-charts.html
#http://mr.schochastics.net/netVizR.html check this out 
#https://kateto.net/networks-r-igraph
#https://rdrr.io/cran/KMDA/man/spearman.group.html
##https://data.library.virginia.edu/getting-started-with-factor-analysis/


############################################################### CODE ###################################################

library(dplyr)
library(purrr) 
#library(table1) 
library(tidyr)
library(magrittr)
#library(sm) 
library(ggplot2)
library(corrplot)
#library(caret) 
library(vegan)
library(grid)
library(cowplot)
library(viridis) 
library(nlme)
library(stringr) 
library(pROC)
library(yardstick)
library(reshape2)  
library(ggpubr)
library(viridis)
#install.packages("ggpubr")
##### STEP 1: Identify changes in metabolites before and after adminstration ####

#read in "light file"  
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt",header=T,sep="\t") 

#read in metadata file 
meta_data <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t")
meta_data$SampleID <- as.character(meta_data$SampleID)
cols <- c("oral5asa_lc","diagnosis","biologics","site_name","Antibiotics","bowelsurg","Participant.ID","steroids","pr5asa","bond5asa","sex","dysbiosis")
meta_data[cols] <- lapply(meta_data[cols], factor) #converts these to a factor variable)
meta_data <- meta_data %>% select(c("SampleID","diagnosis","Participant.ID","dysbiosis","site_name","visit_num","oral5asa_lc","sex","consent_age","Antibiotics")) %>% filter(diagnosis == "CD" | diagnosis == "UC")  

mesal_metab<- inner_join(metab_select,meta_data,by="SampleID") %>% drop_na(diagnosis)
#write.table(mesal_metab, "newusersID.txt", sep="\t", quote=FALSE,row.names=FALSE) #this file helps us find new users shown below(based on binary154)


#I then selected participants who were new users of 5-asa based on metabolite, verified 12.17, can see below  
newids <- c("PSM7J19B", "PSM7J19J", "PSM6XBSE", "PSM6XBVM", "PSM6XBSK", "PSM6XBUG", "PSM6XBRK", "PSM7J1CU", "MSMA26AZ", "MSMB4LZ4",
            "MSM6J2IG", "MSM6J2Q3", "MSM5LLDI", "MSM5LLDS", "HSMA33OJ", "HSMA33M8", "HSMA33OZ", "HSMA33MI", "HSM5MD73", "HSM6XRS8", "HSM6XRRV", 
            "HSM6XRVO", "CSM79HRG", "CSM7KOO9", "CSM5MCXL", "CSM67UDN")


#subset

newiddata <- mesal_metab %>% filter(SampleID %in% newids)
newiddata %>% arrange(Participant.ID,binary154) %>% select(binary154,Participant.ID,visit_num,SampleID,X154.0502_3.83) #this shows the list n=13
#   binary154 Participant.ID visit_num SampleID X154.0502_3.83
#1          0          C3004        13 CSM5MCXL         283552
#2          1          C3004        20 CSM67UDN      507071413
#3          0          C3031#        11 CSM79HRG         172164
#4          1          C3031        14 CSM7KOO9      129818420
#5          0          H4014        13 HSM6XRRV         290011
#6          1          H4014        19 HSM6XRVO     1530601385
#7          0          H4015         8 HSM5MD73         168868
#8          1          H4015        13 HSM6XRS8     3121179635
#9          0          H4035        23 HSMA33OZ        2352694
#10         1          H4035        29 HSMA33MI       82747705
#11         0          H4040#        18 HSMA33OJ         390204
#12         1          H4040        29 HSMA33M8     1468731677
#13         0          M2008#         4 MSM5LLDI        1773492
#14         1          M2008         9 MSM5LLDS      643195341
#15         0          M2028        13 MSM6J2IG         557418
#16         1          M2028        19 MSM6J2Q3      109827300
#17         0          M2071#        12 MSMA26AZ        1036302
#18         1          M2071        16 MSMB4LZ4      323254933
#19         0          P6009         4 PSM6XBRK         169904
#20         1          P6009        25 PSM7J1CU      456525486
#21         0          P6010         8 PSM6XBSK         254310
#22         1          P6010        12 PSM6XBUG      879384985
#23         0          P6012#         4 PSM6XBSE        1059704
#24         1          P6012         7 PSM6XBVM       12783515
#25         0          P6016#         6 PSM7J19B         268885
#26         1          P6016        11 PSM7J19J     1171158465


#find out how many weeks on average between pre and post mesalamine (for text)
weeks <- read.csv("/home/rsm34/5asa/hmp2_metadata_load.csv",stringsAsFactors=FALSE) %>% select("External.ID","week_num") %>% rename(SampleID=External.ID)

weeks2 <- inner_join(newiddata,weeks,by="SampleID") %>% distinct(SampleID,.keep_all=TRUE)
weeks_pre <- weeks2 %>% filter(binary154==0) %>% arrange(Participant.ID) 
weeks_post <- weeks2 %>% filter(binary154==1) %>% arrange(Participant.ID) 
mean(weeks_post$week_num-weeks_pre$week_num) #13.1 
sd(weeks_post$week_num-weeks_pre$week_num) #8.7 

#read in FULL METABOLOMICS file after pre-processing script made it 
#########################
metab_rotated2 <- read.table("/home/rsm34/5asa/metab_data_nozero2.txt",header=T,sep="\t")

###### testing ##########

metabs_newids <- inner_join(newiddata,metab_rotated2,by="SampleID") %>% arrange(Participant.ID)%>% arrange(binary154)

#paired wilcoxon for metabolite data  
df = NULL 
df <- metabs_newids[,15:ncol(metabs_newids)]

kmetabs2 = NULL
kmetabs2 = as.data.frame(matrix(-9,ncol(df),2))
colnames(kmetabs2) = c("feature","p")

for(i in c(1:ncol(df)))
{	b<-wilcox.test(df[,i] ~ metabs_newids$binary154, paired=TRUE)
kmetabs2[i,1]= colnames(df)[i]
kmetabs2[i,2]=b$p.value
}
kmetabs2$FDR <- p.adjust(kmetabs2$p,method="fdr")
kmetabs2 <- kmetabs2 %>% arrange((p))
#write.table(kmetabs2 , "newusersMetabs_all.txt", sep="\t", quote=FALSE,row.names=FALSE)
newusersMetabs <- kmetabs2
newusersmetab <- newusersMetabs


## for future plotting
metabs_plotting2 <- metabs_newids %>% select(matches("SampleID|Participant.ID|nicotinuric.acid|X196.0609_2.81.x|binary154|X1.methylnicotinamide|nicotinate|X154.0502_3.83.x|X312.073_6.98|X458.1946_1.76|X189.0771_3.98|X2.aminoadipate|X373.1254_0.71|X242.0458_3"))
names(metabs_plotting2) <- gsub(".x","",names(metabs_plotting2))
newusers_metabs <- metabs_plotting2 %>% arrange(binary154,Participant.ID)
#write.table(newusers_metabs,"newusers_metabs.txt",sep="\t",quote=FALSE,row.names=FALSE)


################################### 1) Fig 2B: PLOTTING FOR NEW USERS ###############################################
###### write a file out for the new users plot and for the uniref analyses ########

mesallabels <- c("pre", "post")


P2 <- ggpaired(newusers_metabs, x = "binary154", y = "nicotinuric.acid",
               color = "binary154", line.color = "gray", line.size = 0.4,
               palette = c("#3C5488FF", "#DC0000FF"),short.panel.labs = FALSE, xlab="5-ASA use",ylab="rel.abund.",title = "nicotinuric.acid",point.size = 1.5,legend="none") +theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(labels=mesallabels)

P3 <- ggpaired(newusers_metabs, x = "binary154", y = "X2.aminoadipate",
               color = "binary154", line.color = "gray", line.size = 0.4,
               palette = c("#3C5488FF", "#DC0000FF"),short.panel.labs = FALSE, xlab="5-ASA use",ylab="rel.abund.",title = "X2.aminoadipate",point.size = 1.5,legend="none") +theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(labels=mesallabels)

pairedplotfxn <- function(metabolite,caption) {
  ggpaired(newusers_metabs, x = "binary154", y = metabolite,
           color = "binary154", line.color = "gray", line.size = 0.4,
           palette = c("#3C5488FF", "#DC0000FF"),short.panel.labs = FALSE, xlab="5-ASA use",ylab="abund.",title = caption,point.size = 1.5,legend="none",font.x=12,font.y=12,font.tickslab=12,font.main=12) +theme(plot.title = element_text(hjust = 0.5))+
    scale_x_discrete(labels=mesallabels)
}

library(gridExtra)
pp1<-pairedplotfxn("X154.0502_3.83","5-ASA")
pp2<-pairedplotfxn("X2.aminoadipate","2-aminoadipate")
pp3<-pairedplotfxn("X196.0609_2.81","N-Ac-5-ASA")

pp4<-pairedplotfxn("X312.073_6.98","312.073\n module:hippurate")
pp5<-pairedplotfxn("X458.1946_1.76","458.1946\n module:Glycerol 3-phosphate")
pp6<-pairedplotfxn("X189.0771_3.98","189.0771\nmodule:Lithocholate")

pp7<-pairedplotfxn("nicotinuric.acid","Nicotinuric Acid")
pp8<-pairedplotfxn("X1.methylnicotinamide","N1-Methyl-nicotinamide")
pp9<-pairedplotfxn("nicotinate","Nicotinate")


fig2b <-grid.arrange(pp1,pp3,pp2,pp7,pp9,pp8,nrow=2) #5.8x6.9 pdf 

ggsave(filename='/home/rsm34/5asabb3/submission/figures/fig2b.png', plot=fig2b, width = 6.9, height = 5.8, dpi = 600) 