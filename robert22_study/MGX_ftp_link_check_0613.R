getwd()

## 读取cohort 1和2的文件进入R
cohort1_ftp <- read.csv("inputs/cohort1_MGX_ftp.csv", header = T)
cohort2_ftp <- read.csv("inputs/cohort2_MGX_ftp.csv", header = T)
dim(cohort1_ftp)
dim(cohort2_ftp)
which(cohort1_ftp$run_accession %in% cohort2_ftp$run_accession)
which(cohort1_ftp$sample_title %in% cohort2_ftp$sample_title)
cohort1_ftp$submitted_ftp
which(cohort1_ftp$submitted_ftp %in% ";")
grep(";", cohort1_ftp$submitted_ftp)
grep(";", cohort2_ftp$submitted_ftp)
ftp_links <- cohort1_ftp$submitted_ftp[grep(";", cohort1_ftp$submitted_ftp)]
ftp_links <- c(ftp_links,cohort2_ftp$submitted_ftp[grep(";", cohort2_ftp$submitted_ftp)])
length(ftp_links)
length(grep(";", cohort1_ftp$submitted_ftp))+length(grep(";", cohort2_ftp$submitted_ftp))

ftp_links ## Those links in this file are all Samples that have both R1 and R2 sequences

ftp_links <- unlist(strsplit(ftp_links, ";"))
writeLines(ftp_links, "outputs/MGX_ftp.txt")
writeLines(ftp_links[1:2], "outputs/test.txt")
cohort1_metadata <- read.csv("Metadata/metadata_cohort1.csv",stringsAsFactors = FALSE, header = T)
cohort2_metadata <- read.csv("Metadata/metadata_cohort2.csv",stringsAsFactors = FALSE, header = T)
cohort2_metadata <- cohort2_metadata[grep("UC",cohort2_metadata$Diagnosis),]
tube_names <- c(cohort2_metadata$tube_id)
length(tube_names)
tube_names

sr <- as.character(unlist(strsplit(ftp_links[100], "/")[[1]][6]))
sr1 <- strsplit(sr, "[.]")[[1]][2]
sr <- ""
ftp_subs <- c()
for(i in 1:length(ftp_links)){
  sr <- as.character(unlist(strsplit(ftp_links[i], "/")[[1]][6]))
  sr1 <- strsplit(sr, "[.]")[[1]][2]
  ftp_subs <- c(ftp_subs, sr1)
}
ftp_subs
ftp_links[which(ftp_subs %in% tube_names)]
length(which(ftp_subs %in% tube_names))
length(which(tube_names %in% ftp_subs))
length(tube_names)
## cohort2 57/63
tube_ids <- c(cohort1_metadata$id)
sr <- ""
cohort1_ids <- c()
for (i in 1: length(cohort1_ftp$sample_title)) {
  sr <- strsplit(cohort1_ftp$sample_title[i],'[.]')[[1]][2]
  cohort1_ids <- c(cohort1_ids, sr)
}
length(which(tube_ids %in% cohort1_ids))
# all cohort 1 have R1 and R2
(which(ftp_subs %in% tube_names))
cohort1_ftp[which(cohort1_ids %in% tube_ids),]

# 读取用ssh脚本下载的所有ftp文件的size文件
ftp_sizes <- read.csv("inputs/MGX_ftp_size.csv", header = F)
colnames(ftp_sizes) <- c('URLs', 'Sizes (in bytes)')
class(ftp_sizes$`Sizes (in bytes)`[1])
which(ftp_sizes$`Sizes (in bytes)`< 100000000)
ftp_size_sub100 <- ftp_sizes[which(ftp_sizes$`Sizes (in bytes)` < 100000000), ]
ftp_size_over100 <- ftp_sizes[which(ftp_sizes$`Sizes (in bytes)` >= 100000000), ]
dim(ftp_size_sub100)
dim(ftp_size_over100)
dim(ftp_sizes)

for(i in 1:length(ftp_links)){
  sr <- as.character(unlist(strsplit(ftp_links[i], "/")[[1]][6]))
  sr1 <- strsplit(sr, "[.]")[[1]][2]
  ftp_subs <- c(ftp_subs, sr1)
}
unlist(strsplit(cohort1_ftp$submitted_ftp, ";"))
length(which(ftp_size_over100$URLs %in% unlist(strsplit(cohort1_ftp$submitted_ftp, ";"))))
length(which(ftp_size_over100$URLs %in% unlist(strsplit(cohort2_ftp$submitted_ftp, ";"))))

ftp_links_cohort1 <- unlist(strsplit(cohort1_ftp$submitted_ftp, ";"))
ftp_links_cohort2 <- unlist(strsplit(cohort2_ftp$submitted_ftp, ";"))
ftp_subs_cohort1 <- c()
for(i in 1:length(ftp_links_cohort1)){
  sr <- as.character(unlist(strsplit(ftp_links_cohort1[i], "/")[[1]][6]))
  sr1 <- strsplit(sr, "[.]")[[1]][2]
  ftp_subs_cohort1 <- c(ftp_subs_cohort1, sr1)
}

ftp_subs_cohort1
length(which(ftp_subs_cohort1 %in% cohort1_metadata$id))
length(which(ftp_subs_cohort1 == "BLANK1"))

ftp_subs_cohort2 <- c()
for(i in 1:length(ftp_links_cohort2)){
  sr <- as.character(unlist(strsplit(ftp_links_cohort2[i], "/")[[1]][6]))
  sr1 <- strsplit(sr, "[.]")[[1]][2]
  ftp_subs_cohort2 <- c(ftp_subs_cohort2, sr1)
}

over100_cohort2 <- ftp_size_over100$URLs[which(ftp_size_over100$URLs %in% unlist(strsplit(cohort2_ftp$submitted_ftp, ";")))]
ftp_subs_cohort2 <- c()
for(i in 1:length(over100_cohort2)){
  sr <- as.character(unlist(strsplit(over100_cohort2[i], "/")[[1]][6]))
  sr1 <- strsplit(sr, "[.]")[[1]][2]
  ftp_subs_cohort2 <- c(ftp_subs_cohort2, sr1)
}
length(which(ftp_subs_cohort2 %in% cohort2_metadata$tube_id))

ftp_subs_cohort2[which(! ftp_subs_cohort2 %in% cohort2_metadata$tube_id)]
ftp_subs_cohort2 <- sub("000", "", ftp_subs_cohort2)
MGX_results_cohort2 <- read.csv("Metadata/MetaG_gOTU_table.csv", header = T)


MGX_results_cohort2[,ncol(MGX_results_cohort2)]

