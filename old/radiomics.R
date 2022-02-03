suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))
suppressMessages(library(matrixStats))


tils <- readRDS("til.rds")
tils$sample_id <- rownames(tils)
cin <- readRDS("chr_instability.rds")
cin <- cin[-which(grepl(cin$overal_cluster, pattern = "Intermediate")), ]
radiomics <- read.csv2("../../Datasets/TRACERx/Radiomics/TracerX_Stats_DNN_Features_radiomics_20190516.csv",colClasses = "character")

# merge everything
merged_df <- merge(cin, tils , by = "sample_id")
merged_df$short_sample_id <- substr(merged_df$sample_id, 0,8)
merged_df <- merge(merged_df, radiomics , by.x = "short_sample_id", by.y = "Patient")
merged_df$unique_id <- seq(1, nrow(merged_df))

columns_weight <- paste0("X.",seq(5, 261)) # make columns names represent weights on NN
for_heatmap <- as.data.frame(merged_df[, c(columns_weight, "unique_id")]) # get all columns representing weights of NN and unique ID
ids <- for_heatmap$unique_id # remember IDs
for_heatmap <- sapply( for_heatmap, as.numeric ) # transform everything to numeric matrix
exclude <- which(is.na(colMeans2(as.matrix(for_heatmap[, -ncol(for_heatmap)]))))  # find colMean
for_heatmap <- as.data.frame(for_heatmap[, -exclude]) # remove columns which have colmeans == NA
for_heatmap <- for_heatmap[, -which(colnames(for_heatmap) == "unique_id")] # remove id from df for heat map
 
ids <- ids[-exclude] # find unique IDs matching the ones from merged data and without NAs
rownames(for_heatmap) <- ids

# filter merged data by those IDS so I can take TIL, wgii, other...
tmp_merged_df <- merged_df %>%
  filter(unique_id %in% ids)

# make annotation for heatmap
df_columns <- data.frame(til = tmp_merged_df$total.til.score.danahaer.x, tmp_merged_df$overal_cluster, as.numeric(tmp_merged_df$CNN.Prediction.Probability), tmp_merged_df$til_cluster_w2)
colnames(df_columns) <- c("til", "Overal Cluster", "CNN pred. Prob.", "TIl Cluster")
rownames(df_columns) <- tmp_merged_df$unique_id


pheatmap(t(for_heatmap),
         cluster_cols = T,
         cluster_rows = T,
         annotation_col = df_columns,
         scale = "row",
         clustering_method = "ward.D2", show_rownames = F, show_colnames = F, 
         cutree_cols = 2,
         cutree_rows = 2)


