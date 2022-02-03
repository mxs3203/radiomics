library(readxl)
library(pheatmap)
library(factoextra)
library(cluster)
library(tidyverse)
library(plotly)
setwd("~/GenomeDK_local/CancerEvolution/phd/Analysis/radiomics")

features_of_interest = c("original_glszm_ZoneEntropy", "original_glszm_SmallAreaEmphasis", "original_glszm_LargeAreaHighGrayLevelEmphasis", "original_glszm_ZonePercentage",
                         "original_glszm_LargeAreaEmphasis","original_firstorder_Variance", "original_firstorder_Entropy","original_firstorder_RobustMeanAbsoluteDeviation",
                         "original_firstorder_Uniformity","original_gldm_SmallDependenceLowGrayLevelEmphasis", "original_gldm_DependenceVariance", "original_gldm_LargeDependenceEmphasis",
                         "original_gldm_SmallDependenceHighGrayLevelEmphasis", "original_glszm_GrayLevelNonUniformityNormalized", "original_gldm_DependenceEntropy", 
                         "original_gldm_GrayLevelVariance", "original_glszm_GrayLevelVariance",
                         "original_glcm_Idn","original_glcm_Imc1", "original_glcm_ClusterProminence",
                         "original_glcm_SumEntropy", "original_glcm_Correlation","original_glcm_InverseVariance", "original_glcm_Contrast", "original_glcm_Idmn","original_glcm_ClusterShade",
                         "original_glcm_JointAverage", "original_glszm_SmallAreaLowGrayLevelEmphasis", "original_glszm_LowGrayLevelZoneEmphasis","original_glszm_HighGrayLevelZoneEmphasis",
                         "original_glszm_LargeAreaLowGrayLevelEmphasis", "original_gldm_LargeDependenceLowGrayLevelEmphasis",
                         "original_glszm_ZoneVariance",  "original_firstorder_Median", "original_firstorder_Skewness",
                         "original_glcm_Autocorrelation", "original_glcm_JointEntropy", "original_glcm_DifferenceEntropy","original_glcm_DifferenceVariance",
                         "original_glcm_ClusterTendency", "original_glcm_JointEnergy", "original_glcm_DifferenceAverage", "original_glszm_SmallAreaHighGrayLevelEmphasis",
                         "original_gldm_SmallDependenceEmphasis", "original_gldm_LargeDependenceHighGrayLevelEmphasis")
features_of_interest[which(duplicated(features_of_interest))]

pyrad = read_excel("../../Datasets/TRACERx/Radiomics/Radiomics20200801/PyRadiomics_20200718/original.xlsx")
colnames(pyrad)
tmp_pyrad = pyrad[,c(features_of_interest, "sampleid")]
rownames(tmp_pyrad) = tmp_pyrad$sampleid
tmp_pyrad = na.omit(tmp_pyrad)
pheatmap(cor(tmp_pyrad[,features_of_interest]), 
         fontsize = 8)



set.seed(27)
pc <- umap::umap(tmp_pyrad[,features_of_interest],
                 n_neighbors = 20,
                 min_dist = 0.9, 
                 spread = 0.96, 
                 n_components = 2)
pc <- as.data.frame(pc$layout)
tmp_pyrad$V1 = pc$V1
tmp_pyrad$V2 = pc$V2

ggplot(pc, aes(x = V1, y = V2)) +
  geom_point() + 
  theme_minimal()

fviz_nbclust(pc, kmeans, method = "wss", k.max = 8) + theme_minimal() + ggtitle("the Elbow Method")
gap_stat <- clusGap(pc, FUN = kmeans, nstart = 30, K.max = 8, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")
fviz_nbclust(pc, kmeans, method = "silhouette", k.max = 8) + theme_minimal() + ggtitle("The Silhouette Plot")



km2 <- kmeans(pc,centers = 3, nstart = 100)
fviz_cluster(km2, data = pc, ellipse.type = "convex") + theme_minimal() + ggtitle("K = 3") 
tmp_pyrad$cluster  = as.factor(km2$cluster)
ggplot(tmp_pyrad, aes(x = V1, y = V2, color = cluster)) +
  geom_point() + 
  theme_minimal()

pheatmap(cor(tmp_pyrad[,features_of_interest], tmp_pyrad[,c("V1", "V2")]))

corr_matrix = cor(tmp_pyrad[,features_of_interest])


# how to shrink the feature space
# version 1: Mean of features
groupLargeArea = c("original_glszm_LargeAreaHighGrayLevelEmphasis", "original_glszm_LargeAreaEmphasis", "original_glszm_ZoneVariance")
groupSmallArea = c("original_gldm_SmallDependenceLowGrayLevelEmphasis", "original_glszm_SmallAreaLowGrayLevelEmphasis", "original_glszm_LowGrayLevelZoneEmphasis")
groupHighGray = c("original_glszm_HighGrayLevelZoneEmphasis", "original_glszm_SmallAreaHighGrayLevelEmphasis", "original_glcm_JointAverage", "original_glcm_Autocorrelation")
groupDifference = c("original_glszm_ZonePercentage", "original_gldm_SmallDependenceEmphasis", "original_glcm_DifferenceEntropy", "original_glcm_DifferenceAverage")
groupEntropy = c("original_glcm_JointEntropy", "original_glcm_SumEntropy", "original_firstorder_Entropy")
groupVariance = c("original_glcm_ClusterTendency", "original_firstorder_Variance", "original_gldm_GrayLevelVariance")
cols_to_remove = c(groupLargeArea,groupSmallArea, groupHighGray, groupDifference, groupEntropy, groupVariance)

tmp_pyrad$groupLargeArea = rowMeans(tmp_pyrad[,groupLargeArea])
tmp_pyrad$groupSmallArea = rowMeans(tmp_pyrad[,groupSmallArea])
tmp_pyrad$groupHighGray = rowMeans(tmp_pyrad[,groupHighGray])
tmp_pyrad$groupDifference = rowMeans(tmp_pyrad[,groupDifference])
tmp_pyrad$groupEntropy = rowMeans(tmp_pyrad[,groupEntropy])
tmp_pyrad$groupVariance = rowMeans(tmp_pyrad[,groupVariance])
# Version 2: Take 1 and remove the rest
groupLargeArea = c("original_glszm_LargeAreaEmphasis", "original_glszm_ZoneVariance")
groupSmallArea = c("original_glszm_SmallAreaLowGrayLevelEmphasis", "original_glszm_LowGrayLevelZoneEmphasis")
groupHighGray = c("original_glszm_SmallAreaHighGrayLevelEmphasis", "original_glcm_JointAverage", "original_glcm_Autocorrelation")
groupDifference = c( "original_gldm_SmallDependenceEmphasis", "original_glcm_DifferenceEntropy", "original_glcm_DifferenceAverage")
groupEntropy = c( "original_glcm_SumEntropy", "original_firstorder_Entropy")
groupVariance = c( "original_firstorder_Variance", "original_gldm_GrayLevelVariance")
cols_to_remove = c(groupLargeArea,groupSmallArea, groupHighGray, groupDifference, groupEntropy, groupVariance)


features_of_interest_new = features_of_interest
features_of_interest_new = features_of_interest_new[-which(features_of_interest_new %in% cols_to_remove)]
tmp_pyrad = tmp_pyrad[, c(features_of_interest_new,"sampleid")]

rownames(tmp_pyrad) = tmp_pyrad$sampleid
pheatmap(cor(tmp_pyrad[,features_of_interest_new]), 
         fontsize = 8)

set.seed(27)
pc <- umap::umap(scale(tmp_pyrad[,features_of_interest_new], center = T, scale = T),
                 n_neighbors = 40,
                 min_dist = 0.1, 
                 spread = 0.24, 
                 n_components = 3)
pc <- as.data.frame(pc$layout)
tmp_pyrad$V1 = pc$V1
tmp_pyrad$V2 = pc$V2
tmp_pyrad$V3 = pc$V3

ggplot(pc, aes(x = V1, y = V2)) +
  geom_point() + 
  theme_minimal()

fviz_nbclust(pc, kmeans, method = "silhouette", k.max = 8) + theme_minimal() + ggtitle("The Silhouette Plot")
km2 <- kmeans(pc,centers = 4, nstart = 100)
fviz_cluster(km2, data = pc, ellipse.type = "convex") + theme_minimal() + ggtitle("K = 4") 
tmp_pyrad$cluster  = as.factor(km2$cluster)
ggplot(tmp_pyrad, aes(x = V1, y = V2, color = cluster)) +
  geom_point() + 
  theme_minimal()

plot_ly(pc, x = ~V1, y = ~V2, z = ~V3, color = ~km2$cluster) %>% 
  layout(scene = list(xaxis = list(title = 'Simple Texture and Uniformity '),
                      yaxis = list(title = 'Complex texture with Variance'),
                      zaxis = list(title = 'Large Dependence and uniformity')))

pheatmap(cor(tmp_pyrad[,features_of_interest_new], tmp_pyrad[,c("V1", "V2", "V3")]))

ann_df = data.frame(cluster = tmp_pyrad$cluster)
rownames(ann_df) = tmp_pyrad$sampleid

for_heatmap = as.data.frame(scale(t(tmp_pyrad[,features_of_interest_new]),center = T, scale = T))
colnames(for_heatmap) = tmp_pyrad$sampleid
order = tmp_pyrad %>% arrange(cluster) %>% dplyr::select("sampleid","cluster")

pheatmap(for_heatmap %>% dplyr::select(order$sampleid), 
         scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         clustering_method = "ward.D2",
         annotation_col = ann_df,
         color = c("#7879FF", "#A3A3FF","#FFF192","#F6BDC0", "#F1959B", "#F07470", "#EA4C46","#DC1C13"),
         breaks = seq(-4,4, by = 2),
         fontsize = 8)

