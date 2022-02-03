library(readxl)
library(pheatmap)
library(factoextra)
library(cluster)
library(tidyverse)
library(plotly)
library(ggmosaic)
library(ggpubr)
library(ggbiplot)
library(NMF)
suppressMessages(library(survminer))
suppressMessages(library(survplot))
source(".././cnsignatures/functions.R")

indx100 = get(load("../../Datasets/TRACERx/Annotations/indx100.RData"))
setwd("~/GenomeDK_local/CancerEvolution/phd/Analysis/radiomics")
ith <- read.delim2("../../Datasets/TRACERx/Annotations/comb_ith_analysis.20170127.txt", sep = "\t")
ith$SNV_py <- as.numeric(as.character(ith$SNV_py))
ith$SCNA <- as.numeric(as.character(ith$SCNA))
ith$sampleID <- rownames(ith)
ith2 <- read.delim2("../../Datasets/TRACERx/Annotations/scna_ith_trx500.20180731.txt")
ith2$FractionTrunk <- as.numeric(as.character(ith2$FractionTrunk))
ith2$Heterogeneous <- as.numeric(as.character(ith2$Heterogeneous))
ith2$Homogeneous <- as.numeric(as.character(ith2$Homogeneous))

clin_ = as.data.frame(get(load(("../../Datasets/TRACERx/Clinical/tracerx_clin.20180701.RData"))))
clin_$ID = rownames(clin_)
clin_$Lesion1SizePath = as.numeric(as.character(clin_$Lesion1SizePath))
clin = read_excel("~/Downloads/TRACERx all patients- July 2020 version 2.xlsx")
clin$ID = paste0(substr(clin$Shorter_ID,0,3), "0",substr(clin$Shorter_ID,4, 8))
pyrad = merge(pyrad, clin, by.y = "ID", by.x = "sampleid", all.x = T )
pyrad = merge(pyrad, clin_ %>% select("ID", "Lesion1SizePath"), by.y = "ID", by.x = "sampleid", all.x = T )

phyl_groups = read.delim2("../../Datasets/TRACERx/Annotations/phyl_groups.txt")
phyl_groups$sampleid = rownames(phyl_groups)

mutations = read.csv2("../../Datasets/TRACERx/haveDriverMuts.csv", sep=",")

oracle= read.csv2("../../Datasets/TRACERx/Annotations/2020-12-08TX421_oracle.csv", sep = ",")
oracle$sampleid = paste0(substr(oracle$PatientID, 3,5), "0", substr(oracle$PatientID, 6,8))
oracle$RiskScore <- as.numeric(as.character(oracle$RiskScore))
oracle = oracle %>% 
  dplyr::group_by(sampleid) %>%
  dplyr::summarize(RiskScore = mean(RiskScore, na.rm = T),
            n = n(),
            lows = sum(bin == "Low"),
            highs = sum(bin == "High"),
            bin = ifelse(highs >= lows, "High", "Low"))

extra_info = read.delim2("../../Datasets/TRACERx/Annotations/NewITH/20210722_20210531_tracerx_sample_table_primary_421_annotated.tsv")
new_ith = read.delim2("../../Datasets/TRACERx/Annotations/NewITH/20210805_Evo_subgroups.tsv")
new_ith$sampleid = paste0(substr(as.character(new_ith$tumour_id), 0, 3),"0", substr(as.character(new_ith$tumour_id), 4, 6))

pyrad = read_excel("../../Datasets/TRACERx/Radiomics/Radiomics20200801/PyRadiomics_20200718/original.xlsx")

pyrad$pathology = case_when(
  pyrad$pathhistr_tracerx_lesion1form == "Invasive adenocarcinoma" ~ "Invasive adenocarcinoma",
  pyrad$pathhistr_tracerx_lesion1form == "Squamous cell carcinoma" ~ "Squamous cell carcinoma",
  TRUE ~ "Other"
)

features = colnames(pyrad)[52:151]
pyrad = pyrad[-is.na(pyrad),]

pyrad = pyrad[complete.cases(pyrad[,features]), ]

pc = prcomp(pyrad[,features], scale. = T, center = T)
pc_res = as.data.frame(pc$x)

ggbiplot(pc, labels=rownames(pc))


minMax <- function(x){
  (x - min(x))/(max(x)-min(x))
}

pyrad[,features] = apply(pyrad[,features], 2, minMax)
features = features[-grepl("firstorder", features )]

kernel = "original"
features_of_interest = c(
  paste0(kernel ,"_glcm_SumEntropy"),
  paste0(kernel ,"_glcm_Idn"),
  paste0(kernel ,"_glcm_JointEntropy"),
  paste0(kernel ,"_glcm_DifferenceEntropy"),
  paste0(kernel ,"_glcm_DifferenceAverage"),
  #paste0(kernel ,"_glcm_DifferenceVariance"),
  #paste0(kernel ,"_glszm_SmallAreaEmphasis"),
  paste0(kernel ,"_glcm_Idmn"),
  #paste0(kernel ,"_glrlm_GrayLevelNonUniformityNormalized"),
  paste0(kernel ,"_glrlm_LongRunEmphasis"),
  paste0(kernel ,"_gldm_SmallDependenceEmphasis"),
  paste0(kernel ,"_gldm_LargeDependenceEmphasis"),
  paste0(kernel ,"_glrlm_RunPercentage")
  #paste0(kernel ,"_glrlm_ShortRunEmphasis"),
  #paste0(kernel ,"_glszm_LargeAreaEmphasis"),
  #paste0(kernel ,"_glcm_Imc1"),
  #paste0(kernel ,"_glcm_Autocorrelation")
)

pyrad[,features_of_interest] = pyrad[,features_of_interest]
pyrad$original_glcm_Idn = -pyrad$original_glcm_Idn
pyrad$original_glcm_Idmn = -pyrad$original_glcm_Idmn
pyrad$original_glrlm_LongRunEmphasis = -pyrad$original_glrlm_LongRunEmphasis
pyrad$original_gldm_LargeDependenceEmphasis = -pyrad$original_gldm_LargeDependenceEmphasis


pc = prcomp(pyrad[,features_of_interest],center = T)
pc_res = as.data.frame(pc$x)
pc_res$dist_pc = as.matrix(dist(pc_res))[nrow(pc_res), ]
pc_res$sampleid = pyrad$sampleid

umap_obj = umap(pyrad[,features_of_interest], n_components = 3)
umap_res = as.data.frame(umap_obj$layout)
umap_res$dist_umap = as.matrix(dist(umap_res))[nrow(umap_res), ]
umap_res$sampleid = pyrad$sampleid

nmf_res = NMF::nmf(pyrad[,features_of_interest], 3) 
nmf_res = as.data.frame(nmf_res@fit@W)
nmf_res$dist_nmf = as.matrix(dist(nmf_res))[nrow(nmf_res), ]
nmf_res$sampleid = pyrad$sampleid

res = merge(nmf_res,umap_res, by = "sampleid" )
res = merge(res,pc_res, by = "sampleid" )

p1 = ggplot(res, aes (x = dist_nmf, dist_umap)) +
   geom_point() +
  theme_minimal()
p2 = ggplot(res, aes (x = dist_nmf, dist_pc)) +
  geom_point() +
  theme_minimal()
p3 = ggplot(res, aes (x = dist_pc, dist_umap)) +
  geom_point() +
  theme_minimal()
ggarrange(p1,p2,p3)

ggbiplot(pc, labels=rownames(pc)) + 
  theme_pubclean()


pheatmap(cor(pyrad[,features_of_interest]), 
         fontsize = 8)


# res = IntNMF::nmf.opt.k(dat = list(as.matrix(pyrad[,features_of_interest])),
#                   n.runs = 30,
#                   n.fold = 3,
#                   k.range = 2:10,
#                   result = TRUE,
#                   make.plot = TRUE,
#                   progress = TRUE)

results = data.frame()
for (k in 2:10){
  rss = c()
  rss_rand = c()
  for (i in 1:1){
    # pyrad data
    nmf_res = NMF::nmf(pyrad[,features_of_interest], k) 
    x <- pyrad[1:nrow(nmf_res@fit@W),features_of_interest]
    rss <- c(rss, rss(nmf_res@fit, x))
    # random data
    #  rmatrix(nrow(pyrad[,features_of_interest]),ncol(pyrad[,features_of_interest]) )
    rand_matrix = rmatrix(nrow(pyrad[,features_of_interest]),ncol(pyrad[,features_of_interest]) )#pyrad[sample(nrow(pyrad[,features_of_interest])),sample(features_of_interest)]
    nmf_res = NMF::nmf(rand_matrix, k) 
    rss_rand <- c(rss_rand, rss(nmf_res@fit, rand_matrix))
  }
  rss_df = data.frame(rss=rss, rss_rand=rss_rand)
  rss_df$k = k
  results = rbind(results, rss_df)
}



results = results %>% 
  dplyr::group_by(k) %>% 
  dplyr::summarise(rss = mean(rss), rss_rand = mean(rss_rand))
ggplot(results ) + 
  geom_point(aes(x = k, y = rss, color ="RSS")) + 
  geom_line(aes(x = k, y = rss, color ="RSS")) + 
  geom_point(aes(x = k, y = rss_rand, color ="Rand RSS")) + 
  geom_line(aes(x = k, y = rss_rand, color ="Rand RSS")) +
  theme_pubclean()







melted_H = reshape2::melt(nmf_res@fit@H)
ggplot(melted_H, aes(x = Var2, y = value, fill = as.factor(Var1))) + 
  geom_bar(stat='identity', position='dodge') + 
  theme_pubclean()+ 
  theme(axis.text.x = element_text(angle = 90))

pyrad$volume = as.numeric(pyrad$lesion1sizepath)
pyrad$volume_from_pyrad = pyrad$original_shape_MeshVolume
pyrad$diameter = pyrad$original_shape_Maximum2DDiameterSlice

pyrad$original_glrlm_LongRunEmphasis = pyrad$original_glrlm_LongRunEmphasis/pyrad$volume_from_pyrad
pyrad$radITH = rowMeans(pyrad[,features_of_interest], na.rm = T)

Q = 3
pyrad$volume_group = gtools::quantcut(pyrad$volume, q=Q, na.rm=TRUE)
pyrad$diameter_group = gtools::quantcut(pyrad$diameter, q=Q, na.rm=TRUE)
pyrad$radITH_group = gtools::quantcut(pyrad$radITH, q=Q, na.rm=TRUE)

predicted_W = read.delim2("~/GenomeDK_local/CancerEvolution/phd/Analysis/radiomics/NMF_Radiomics/tracerx_predicted.csv",sep = ",")
pyrad = merge(pyrad, predicted_W, by = "sampleid")
colnames(pyrad)[(length(colnames(pyrad))-2):length(colnames(pyrad))] = c("S1", "S2","S3")
pyrad$S1 = as.numeric(pyrad$S1)
pyrad$S2 = as.numeric(pyrad$S2)
pyrad$S3 = as.numeric(pyrad$S3)

p1 = ggplot(pyrad, aes(x = S1,  y = volume)) +
  geom_point() + 
  geom_smooth(method = "lm")+
  stat_cor()+
  scale_y_log10()+
  theme_minimal() +
  scale_color_brewer(palette="Set1") 
p2 = ggplot(pyrad, aes(x = S2,  y = volume)) +
  geom_point() + 
  geom_smooth(method = "lm")+
  stat_cor()+
  scale_y_log10()+
  theme_minimal() +
  scale_color_brewer(palette="Set1") 
p3 = ggplot(pyrad, aes(x = S3,  y = volume)) +
  geom_point() + 
  geom_smooth(method = "lm")+
  stat_cor()+
  scale_y_log10()+
  theme_minimal() +
  scale_color_brewer(palette="Set1") 

ggarrange(p1,p2, p3)

pyrad = merge(pyrad, new_ith,by = "sampleid", all.x = T)
p1 = ggplot(pyrad %>% filter(!is.na(evo_major_group)), aes(x = evo_major_group,  y = S1, fill = evo_minor_group )) +
  geom_boxplot() +
  ylab("S1") +
  xlab("ITH Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 
p2 = ggplot(pyrad %>% filter(!is.na(evo_major_group)), aes(x = evo_major_group,  y = S2, fill = evo_minor_group )) +
  geom_boxplot() +
  ylab("S2") +
  xlab("ITH Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 
p3 = ggplot(pyrad %>% filter(!is.na(evo_major_group)), aes(x = evo_major_group,  y = S3, fill = evo_minor_group )) +
  geom_boxplot() +
  ylab("S3") +
  xlab("ITH Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 

ggarrange(p1,p2, p3, ncol=3)


colnames(phyl_groups) = c ("phyl_g", "sampleid")
pyrad = merge(pyrad, phyl_groups,by = "sampleid", all.x = T)

p1=ggplot(pyrad %>% filter(!is.na(phyl_g)), aes(x = phyl_g,  y = S1, fill = phyl_g )) +
  geom_boxplot() +
  ylab("S1") +
  xlab("PHYL Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 
p2=ggplot(pyrad %>% filter(!is.na(phyl_g)), aes(x = phyl_g,  y = S2, fill = phyl_g )) +
  geom_boxplot() +
  ylab("S2") +
  xlab("PHYL Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 
p3=ggplot(pyrad %>% filter(!is.na(phyl_g)), aes(x = phyl_g,  y = S3, fill = phyl_g )) +
  geom_boxplot() +
  ylab("S3") +
  xlab("PHYL Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 
ggarrange(p1,p2)


pyrad = merge(pyrad, oracle,by = "sampleid", all.x = T)

p1= ggplot(pyrad, aes(x = S1,  y = RiskScore)) +
  geom_point() +
  ylab("ORACLE Risk Score") +
  xlab("S1") +
  geom_smooth(method = 'lm') +
  theme_minimal() +
  stat_cor() +
  scale_fill_brewer(palette="Set1") 
p2= ggplot(pyrad, aes(x = S2,  y = RiskScore)) +
  geom_point() +
  ylab("ORACLE Risk Score") +
  xlab("S2") +
  geom_smooth(method = 'lm') +
  theme_minimal() +
  stat_cor() +
  scale_fill_brewer(palette="Set1") 
p3= ggplot(pyrad, aes(x = S3,  y = RiskScore)) +
  geom_point() +
  ylab("ORACLE Risk Score") +
  xlab("S3") +
  geom_smooth(method = 'lm') +
  theme_minimal() +
  stat_cor() +
  scale_fill_brewer(palette="Set1") 

ggarrange(p1,p2, p3)


pyrad$quantcut_group = gtools::quantcut(pyrad$S3, q=2, na.rm=TRUE)
pyrad = merge(pyrad, mutations, by.y = "SampleID", by.x = "sampleid", all.x = T)
ggplot(pyrad %>% filter(pathology != "Other"), aes(x = S3,  
                                                       y = volume_from_pyrad, 
                                                       color = haveEGFR )) +
  geom_point() + 
  scale_y_log10()+
  #scale_x_log10()+
  geom_hline(yintercept = median(pyrad$volume_from_pyrad, na.rm = T) ,linetype ="dotted")+
  geom_vline(xintercept = median(pyrad$S3, na.rm = T),linetype ="dotted" ) +
  ylab("Log10(Volume)") +
  theme_minimal() +
  scale_color_brewer(palette="Set1")+ 
  facet_wrap(~pathology)


pyrad = merge(pyrad, ith2, by.x = "sampleid", by.y = "SampleID", all.x = T)
p1= ggplot(pyrad, aes(x = S1,  y = Homogeneous)) +
  geom_point() + 
  stat_cor()+
  geom_smooth(method = "lm", se= F)+
  ylab("Homogeneous") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") 
p2= ggplot(pyrad, aes(x = S1,  y = Heterogeneous)) +
  geom_point() + 
  stat_cor()+
  geom_smooth(method = "lm",se = F)+
  ylab("Heterogeneious") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") 
p3= ggplot(pyrad, aes(x = S1,  y = FractionTrunk )) +
  geom_point() + 
  geom_smooth(method = "lm", se = F)+
  ylab("FractionTrunk") +
  stat_cor()+
  theme_minimal() +
  scale_color_brewer(palette="Set1") 


ggarrange(p1,p2,p3, nrow = 1,ncol=3)


vega = read.delim2("../../Datasets/sanchez-vega-pws_1026.csv", sep = ";")
muts = get(load("../../Datasets/TRACERx/WES/patientMutTable_20200622.RData"))
muts$SampleID = sub('LTX','LTX0',substr(muts$SampleID, 0,8))
muts = muts %>% 
  filter(DriverMut == TRUE )

results = data.frame()
samples = c()
for (sample in unique(muts$SampleID)){
  samples = c(samples, sample)
  tmp_muts = muts %>% filter(SampleID == sample)
  have_path_mut = c()
  for (geneGroup in unique(vega$Pathway_pretty)) {
    vega_group = vega %>% filter(Pathway_pretty == geneGroup) 
    if (any(vega_group$Gene %in% tmp_muts$Hugo_Symbol == T)) {
      have_path_mut = c(have_path_mut, 1)
    } else {
      have_path_mut = c(have_path_mut, 0)
    }
  }
  results = rbind(results, have_path_mut)
}
colnames(results) = as.character(unique(vega$Pathway))
rownames(results) = samples
results$SampleID = rownames(results)


pyrad = merge(pyrad, results, by.y = "SampleID", by.x = "sampleid", all.x = T)

pyrad$quantcut_group = gtools::quantcut(pyrad$S3, q=3, na.rm=TRUE)
tmp = pyrad %>% filter(pathology == "Invasive adenocarcinoma")
for (col in unique(vega$Pathway)){
  if(length(unique(tmp[,col])) >= 2){
    
    test = fisher.test(table(tmp$radITH_group, tmp[,col])) 
    if (test$p.value <= 0.08){
      print(col)
      print(test)
    }
    
  }
}


tmp = pyrad %>% filter(pathology == "Invasive adenocarcinoma", !is.na(p53))

p1 = ggplot(tmp,
       aes(x = as.factor(p53),y = radITH,fill = as.factor(p53 ))) +
  geom_boxplot() + 
  scale_y_log10()+
  ylab("radITH") +
  stat_compare_means(method = "wilcox") +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  ggtitle("Adeno: p53")
p2 = ggplot(tmp,
            aes(x = as.factor(p53),y = volume,fill = as.factor(p53 ))) +
  geom_boxplot() + 
  scale_y_log10()+
  ylab("Volume") +
  stat_compare_means(method = "wilcox") +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  ggtitle("Adeno:p53")
p3 = ggplot(tmp,
            aes(x = as.factor(p53),y = diameter,fill = as.factor(p53 ))) +
  geom_boxplot() + 
  scale_y_log10()+
  ylab("2D Diameter") +
  stat_compare_means(method = "wilcox") +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  ggtitle("Adeno: p53")

fisher.test(table(tmp$p53, tmp$radITH_group))
p4 = ggplot(data = tmp) + 
  geom_mosaic(aes(x=product(p53, radITH_group), fill = radITH_group)) + 
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  ggtitle("Fisher: 0.45")
fisher.test(table(tmp$p53, tmp$volume_group))
p5 = ggplot(data = tmp) + 
  geom_mosaic(aes(x=product(p53, volume_group), fill = volume_group)) + 
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  ggtitle("Fisher: 0.043")

fisher.test(table(tmp$p53, tmp$diameter_group))
p6 = ggplot(data = tmp) + 
  geom_mosaic(aes(x=product(p53, diameter_group), fill = diameter_group)) + 
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  ggtitle("Fisher: 0.051")

ggarrange(p1,p2,p3,p4,p5,p6, ncol = 3, nrow = 2)


surv = read.delim2("TCGA_survival_data_clean.txt")

  
p3 = ggplot(tmp,
            aes(x = radITH,y = volume_from_pyrad,color = as.factor(pi3k ))) +
  geom_point()+
  scale_y_log10()+
  ylab("Volume") +
  stat_cor()+
  geom_smooth(method = "lm") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") +
  ggtitle("Adeno: PI3K")




ggplot(pyrad %>% filter(pathology == "Squamous cell carcinoma"), aes(x = S1,  
                                                   y = volume, 
                                                   color = as.factor(wnt ))) +
  geom_point() + 
  scale_y_log10()+
  #scale_x_log10()+
  geom_hline(yintercept = median(pyrad$volume, na.rm = T) ,linetype ="dotted")+
  geom_vline(xintercept = median(pyrad$S1, na.rm = T),linetype ="dotted" ) +
  ylab("Log10(Volume)") +
  theme_minimal() +
  scale_color_brewer(palette="Set1")+ 
  facet_wrap(~pathology)

clin = read_excel("~/Downloads/TRACERx all patients- July 2020 version 2.xlsx")
clin$ID = paste0(substr(clin$Shorter_ID,0,3), "0",substr(clin$Shorter_ID,4, 8))
pyrad = merge(pyrad, clin, by.y = "ID", by.x = "sampleid", all.x = T )

pyrad$quantcut_group = gtools::quantcut(pyrad$S3, q=3, na.rm=TRUE)
tmp = pyrad %>% filter(pathology == "Invasive adenocarcinoma")
surv_os <- Surv(as.numeric(as.character(tmp[,'lung_specific_time']))/365, as.numeric(as.character(tmp[,'cens_lung_specific'])))
fit_os <- survfit(survplot::censor(surv_os, 5)~quantcut_group, data = tmp)
makeSurvPlot(fit_os, 
             "Adeno", 
             legen_title = "S3 quantiles", 
             ylab = "Lung spec. surv.",legend_coord = c(0.2,0.5),colors = c()) 



median_vol = median(pyrad$volume)
median_radITH = median(pyrad$S3)

pyrad$volume_group_med = ifelse(pyrad$volume <= median_vol, "Low", "High")
pyrad$ith_group_med = ifelse(pyrad$S3 <= median_radITH, "Low", "High")
pyrad$iith_volume_group = paste0("Volume: ",pyrad$volume_group_med, "- S3: ",pyrad$ith_group_med )


# "Invasive adenocarcinoma"
# "Squamous cell carcinoma"
pyrad$quantcut_group = gtools::quantcut(pyrad$S3, q=3, na.rm=TRUE)
tmp = pyrad %>% filter(pathology =="Invasive adenocarcinoma", !is.na(pi3k))

ggplot(tmp, 
            aes(x = S1,  y = volume, color = as.factor(rtk_kras))) +
  geom_point() + 
  scale_y_log10()+
  #scale_x_log10()+
  geom_hline(yintercept = median(pyrad$volume, na.rm = T), linetype='dotted', col = 'red') + 
  geom_vline(xintercept = median(pyrad$S1, na.rm = T), linetype='dotted', col = 'red') + 
  ylab("Log10(Volume)") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") 

ggplot(data = tmp) + 
  geom_mosaic(aes(x=product(pi3k, quantcut_group), fill = quantcut_group)) + 
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") 

ggarrange(p1,p2, nrow = 1)


colnames(clin)
require(pROC)
library(caret)

pyrad %>% dplyr::group_by(pyrad$haveEGFR) %>% dplyr::summarise(n = n())
params_to_predict = c("rtk_kras", "haveEGFR", "p53", "pi3k", 
                      "cell_cycle","haveKRAS", "haveNFE2L2", "haveSTK11", "nrf2")
res = data.frame()
for (p in params_to_predict){
  print(p)
  acc_test = c()
  acc_train = c()
  for (i in 1:50){
    features = c("S1","S2","S3","volume_from_pyrad", "pathology","sex")
    response_var = p
    
    for_modeling = pyrad[complete.cases(pyrad[,c(features, response_var)]),]
    for_modeling$pathology = as.factor(for_modeling$pathology)
    for_modeling$sex = as.factor(for_modeling$sex)
    
    train.index <- createDataPartition(for_modeling[,response_var], p = .75, list = FALSE)
    train <- for_modeling[ train.index,]
    test  <- for_modeling[-train.index,]
    
    X = train[,features]
    Y = as.factor(train[,response_var])
    model = train(
      x = X,
      y = Y,
      trControl = trainControl(method = "cv", number = 3),
      method = "glm",
      family = "binomial"
    )
    y_hat_prob = predict(model, test[,features],"prob")
    
    rf.roc<-roc(test[,response_var],y_hat_prob[,1])
    acc_test <- c(acc_test, auc(rf.roc))
    
    model$importance
  }
  t = data.frame(acc_test)
  colnames(t) = c("AUC")
  t$param = p
  res = rbind(res,t)
}


    


ggplot(res, aes(x = param, y = AUC)) +
       geom_boxplot() + ylab("AUC")+ 
      theme_minimal()



bins = as.data.frame(apply(X = pyrad[,features_of_interest], MARGIN = 2, function(x){
  gtools::quantcut(x, q=2, na.rm=TRUE)
}))



write.csv2(bins,"~/PycharmProjects/PythonLDA/pyrad_bins.csv")
