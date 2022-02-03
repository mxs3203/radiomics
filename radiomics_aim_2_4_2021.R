library(readxl)
library(pheatmap)
library(factoextra)
library(cluster)
library(tidyverse)
library(plotly)
library(ggmosaic)
library(ggpubr)
source(".././cnsignatures/functions.R")
#####################################################################################
# This functions takes gene symbols and returns a vector of entrez or ensembl ID's with symbols as names
#####################################################################################
map2entrez <- function(gene_symbols, translate_into='entrez'){
  require(org.Hs.eg.db)
  entrez.alias <- as.list(org.Hs.egALIAS2EG)
  entrez.symbol <- as.list(org.Hs.egSYMBOL2EG)
  
  entrez <- gene_symbols
  for(i in 1:length(entrez)){
    tmp <- entrez.symbol[entrez[i]][[1]][1]
    if(is.null(tmp)){
      tmp <- entrez.alias[entrez[i]][[1]][1]
    }
    if(is.null(tmp)){
      tmp <- NA
    }
    names(entrez)[i] <- tmp
    out <- entrez
  }
  if(translate_into=='ensembl'){
    entrez.ensembl <- as.list(org.Hs.egENSEMBL)
    ensembl <- entrez
    for(i in 1:length(ensembl)){
      tmp <- entrez.ensembl[names(ensembl[i])][[1]][1]
      if(is.null(tmp)){
        tmp <- NA
      }
      names(ensembl)[i] <- tmp
    }
    out <- ensembl
  }
  return(out)
}
# necessary gene definition
immune_genes <- read.table("../Immune/20190102/selected_markers.txt", sep = '\t', header = TRUE, as.is = TRUE)
immune_gene_groups <- c('CD8 T cells', 'B-cells', 'CD45', 'Cytotoxic cells', 'DC', 'Exhausted CD8', 'Macrophages', 'Mast cells', 'Neutrophils', 'NK CD56dim cells', 'NK cells', 'T-cells','Th1 cells', 'Treg', 'CD4')

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

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

phyl_groups = read.delim2("../../Datasets/TRACERx/Annotations/phyl_groups.txt")
phyl_groups$sampleid = rownames(phyl_groups)

colnames(pyrad)

features_of_interest = c(
                         "original_glcm_SumEntropy", 
                         "original_glcm_Idn",
                         "original_glcm_JointEntropy", 
                         "original_glcm_DifferenceEntropy",
                         #"original_glcm_DifferenceVariance",
                          "original_glcm_DifferenceAverage",
                         # "original_glszm_SmallAreaEmphasis",
                         "original_glcm_Idmn",
                         #"original_glrlm_GrayLevelNonUniformityNormalized",
                         "original_glrlm_LongRunEmphasis",
                         #"original_glrlm_ShortRunEmphasis",
                         "original_glrlm_RunPercentage",
                         "original_gldm_SmallDependenceEmphasis",
                         "original_gldm_LargeDependenceEmphasis")
                         

mutations = read.csv2("../../Datasets/TRACERx/haveDriverMuts.csv", sep=",")

oracle= read.csv2("../../Datasets/TRACERx/Annotations/2020-12-08TX421_oracle.csv", sep = ",")
oracle$sampleid = paste0(substr(oracle$PatientID, 3,5), "0", substr(oracle$PatientID, 6,8))
oracle$RiskScore <- as.numeric(as.character(oracle$RiskScore))
oracle = oracle %>% 
  group_by(sampleid) %>%
  summarize(RiskScore = mean(RiskScore, na.rm = ),
            n = n(),
            lows = sum(bin == "Low"),
            highs = sum(bin == "High"),
            bin = ifelse(highs >= lows, "High", "Low"))

extra_info = read.delim2("../../Datasets/TRACERx/Annotations/NewITH/20210722_20210531_tracerx_sample_table_primary_421_annotated.tsv")
new_ith = read.delim2("../../Datasets/TRACERx/Annotations/NewITH/20210805_Evo_subgroups.tsv")
new_ith$sampleid = paste0(substr(as.character(new_ith$tumour_id), 0, 3),"0", substr(as.character(new_ith$tumour_id), 4, 6))

pyrad = read_excel("../../Datasets/TRACERx/Radiomics/Radiomics20200801/PyRadiomics_20200718/original.xlsx")
colnames(pyrad)
pyrad$pathology = case_when(
  pyrad$pathhistr_tracerx_lesion1form == "Invasive adenocarcinoma" ~ "Invasive adenocarcinoma",
  pyrad$pathhistr_tracerx_lesion1form == "Squamous cell carcinoma" ~ "Squamous cell carcinoma",
  TRUE ~ "Other"
)

tmp_pyrad = pyrad[,c(features_of_interest, "sampleid", "original_shape_MeshVolume", "pathology", "original_shape_Maximum2DDiameterSlice")]
rownames(tmp_pyrad) = tmp_pyrad$sampleid
tmp_pyrad = na.omit(tmp_pyrad)
tmp_pyrad$original_glcm_Idn = -tmp_pyrad$original_glcm_Idn
tmp_pyrad$original_glcm_Idmn = -tmp_pyrad$original_glcm_Idmn
tmp_pyrad$original_glrlm_LongRunEmphasis = -tmp_pyrad$original_glrlm_LongRunEmphasis
tmp_pyrad$original_gldm_LargeDependenceEmphasis = -tmp_pyrad$original_gldm_LargeDependenceEmphasis
tmp_pyrad$original_glrlm_GrayLevelNonUniformityNormalized = -tmp_pyrad$original_glrlm_GrayLevelNonUniformityNormalized

pheatmap(cor(tmp_pyrad[,features_of_interest]), 
         fontsize = 8)


clin = read_excel("~/Downloads/TRACERx all patients- July 2020 version 2.xlsx")
clin$ID = paste0(substr(clin$Shorter_ID,0,3), "0",substr(clin$Shorter_ID,4, 8))
tmp_pyrad = merge(tmp_pyrad, clin, by.y = "ID", by.x = "sampleid", all.x = T )


tmp_pyrad$volume = tmp_pyrad$original_shape_MeshVolume/1000.0
tmp_pyrad$diameter = tmp_pyrad$original_shape_Maximum2DDiameterSlice/1000.0
corrplot::corrplot(cor(tmp_pyrad[,features_of_interest], tmp_pyrad[,"volume"]), method = "number")

tmp_pyrad[,features_of_interest] = apply(tmp_pyrad[,features_of_interest],2, FUN=normalize)

# Normalizing by volume
tmp_pyrad$original_glrlm_LongRunEmphasis = tmp_pyrad$original_glrlm_LongRunEmphasis/tmp_pyrad$volume

# row mean
tmp_pyrad$rad_ith_mean = rowMeans(tmp_pyrad[,features_of_interest])

#tmp_pyrad[,"volume"] = apply(tmp_pyrad[,"volume"],2, FUN=normalize)
hist(tmp_pyrad$volume, breaks = 200)
quantile(tmp_pyrad$volume)
Q = 3
tmp_pyrad$volume_group = gtools::quantcut(tmp_pyrad$volume, q=Q, na.rm=TRUE)
tmp_pyrad$quantcut_group = gtools::quantcut(tmp_pyrad$rad_ith_mean, q=Q, na.rm=TRUE)
tmp_pyrad$diameter_group = gtools::quantcut(tmp_pyrad$Lesion1_size_pathology, q=Q, na.rm=TRUE)


ggplot(tmp_pyrad, aes(x = rad_ith_mean,  y = volume, color = diameter_group )) +
  geom_point() + 
  scale_y_log10()+
  geom_smooth()+
  #scale_x_log10()+
  theme_minimal() +
  scale_color_brewer(palette="Set1") 

tmp_pyrad = merge(tmp_pyrad, ith2,by.x = "sampleid", by.y = "SampleID", all.x = T)
tmp_pyrad = merge(tmp_pyrad, new_ith,by = "sampleid", all.x = T)
p1= ggplot(tmp_pyrad, aes(x = diameter,  y = Homogeneous)) +
  geom_point() + 
  stat_cor()+
  scale_x_log10()+
  geom_smooth(method = "lm", se= F)+
  ylab("Homogeneous") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") 
p2= ggplot(tmp_pyrad, aes(x = diameter,  y = Heterogeneous)) +
  geom_point() + 
  stat_cor()+
  scale_x_log10()+
  geom_smooth(method = "lm",se = F)+
  ylab("Heterogeneious") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") 
p3= ggplot(tmp_pyrad, aes(x = diameter,  y = FractionTrunk )) +
  geom_point() + 
  scale_x_log10()+
  geom_smooth(method = "lm", se = F)+
  ylab("FractionTrunk") +
  stat_cor()+
  theme_minimal() +
  scale_color_brewer(palette="Set1") 


ggarrange(p1,p2,p3, nrow = 1,ncol=3)

p1 = ggplot(tmp_pyrad %>% filter(!is.na(evo_major_group)), aes(x = evo_major_group,  y = rad_ith_mean, fill = evo_minor_group )) +
  geom_boxplot() +
  ylab("iITH(mean)") +
  xlab("ITH Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 

p2 = ggplot(tmp_pyrad %>% filter(!is.na(evo_major_group)), aes(x = evo_major_group,  y = volume, fill = evo_minor_group )) +
  geom_boxplot() +
  ylab("Volume") +
  xlab("ITH Groups") +
  scale_y_log10()+
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 

p3 = ggplot(tmp_pyrad %>% filter(!is.na(evo_major_group)), aes(x = evo_major_group,  y = rad_ith_mean, fill = evo_major_group )) +
  geom_boxplot() +
  ylab("iITH(mean)") +
  xlab("ITH Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 
p4 = ggplot(tmp_pyrad %>% filter(!is.na(evo_major_group)), aes(x = evo_major_group,  y = volume, fill = evo_major_group )) +
  geom_boxplot() +
  ylab("Volume") +
  xlab("ITH Groups") +
  scale_y_log10()+
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 
ggarrange(p1,p3,p2,p4)

features_of_interest

p1=ggplot(tmp_pyrad, aes(x = original_glcm_SumEntropy, y = rad_ith_mean, color = volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p2=ggplot(tmp_pyrad, aes(x = original_glcm_Idn, y = rad_ith_mean, color = volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p3=ggplot(tmp_pyrad, aes(x = original_glcm_JointEntropy, y = rad_ith_mean, color =volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p4=ggplot(tmp_pyrad, aes(x = original_glcm_DifferenceEntropy, y = rad_ith_mean, color = volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p5=ggplot(tmp_pyrad, aes(x = original_glrlm_LongRunEmphasis, y = rad_ith_mean, color = volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p6=ggplot(tmp_pyrad, aes(x = original_glcm_DifferenceAverage, y = rad_ith_mean, color = volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p7=ggplot(tmp_pyrad, aes(x = original_glrlm_RunPercentage, y = rad_ith_mean, color = volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p8=ggplot(tmp_pyrad, aes(x = original_glcm_Idmn, y = rad_ith_mean, color =volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p9=ggplot(tmp_pyrad, aes(x = original_gldm_LargeDependenceEmphasis, y = rad_ith_mean, color = volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")
p10=ggplot(tmp_pyrad, aes(x = original_gldm_SmallDependenceEmphasis, y = rad_ith_mean, color = volume_group)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_pubclean()  + ylab("iITH(mean)")


ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow = 3,ncol = 4)


colnames(phyl_groups) = c ("phyl_g", "sampleid")
tmp_pyrad = merge(tmp_pyrad, phyl_groups,by = "sampleid", all.x = T)
ggplot(tmp_pyrad) + 
  geom_mosaic(aes(x = product(quantcut_group, phyl_g), fill = phyl_g), na.rm=TRUE) +
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  ylab("iITH Groups") + xlab("Phyl Groups")

ggplot(tmp_pyrad %>% filter(!is.na(phyl_g)), aes(x = phyl_g,  y = diameter, fill = phyl_g )) +
  geom_boxplot() +
  ylab("Diameter") +
  xlab("PHYL Groups") +
  scale_y_log10()+
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 


tmp_pyrad = merge(tmp_pyrad, oracle,by = "sampleid", all.x = T)
ggplot(tmp_pyrad) + 
  geom_mosaic(aes(x = product(volume_group, bin), fill = bin), na.rm=TRUE) +
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  ylab("iITH Groups") + xlab("Oracle")

ggplot(tmp_pyrad, aes(x = quantcut_group,  y = RiskScore, fill = quantcut_group )) +
  geom_boxplot() +
  ylab("ORACLE Risk Score") +
  xlab("radITH Groups") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 

ggplot(tmp_pyrad %>% filter(rad_ith_mean > 0.01), aes(x = rad_ith_mean,  y = RiskScore)) +
  geom_point() +
  ylab("ORACLE Risk Score") +
  xlab("radITH") +
  theme_minimal() +
  scale_x_log10()+
  stat_cor() + geom_smooth(method = "lm") + 
  scale_fill_brewer(palette="Set1") 



quantile(tmp_pyrad$rad_ith_mean)
tmp_pyrad$rad_ith_mean[which(tmp_pyrad$rad_ith_mean > 1)] = 1
tmp_pyrad$quantcut_group = gtools::quantcut(tmp_pyrad$rad_ith_mean, q=3, na.rm=TRUE)
tmp_pyrad_mut = merge(tmp_pyrad, mutations, by.y = "SampleID", by.x = "sampleid", all.x = T)

median(tmp_pyrad$rad_ith_mean)

ggplot(tmp_pyrad %>% filter(pathology != "Other"), aes(x = rad_ith_mean,  
                                                       y = volume, 
                                                       color = evo_major_group)) +
  geom_point() + 
  scale_y_log10()+
  #scale_x_log10()+
  xlab("iITH(mean)") +
  geom_hline(yintercept = 15.57925 ,linetype ="dotted")+
  geom_vline(xintercept =0.4720226,linetype ="dotted" ) +
  ylab("Log10(Volume)") +
  theme_minimal() +
  scale_color_brewer(palette="Set1")+ 
  facet_wrap(~pathology)

tmp = tmp_pyrad %>% filter(pathology != "Other")
plot(table(tmp$pathology, tmp$volume_group), main = "")
ggplot(tmp) + 
  geom_mosaic(aes(x = product(volume_group, pathology), fill = pathology), na.rm=TRUE) +
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  ylab("Vol Groups") + 
  geom_text(data = layer_data(ggplot2::last_plot(), 1) %>% filter(.wt > 0),
            aes(x = (xmin + xmax) / 2,
                y = (ymin + ymax) / 2,
                label = .wt)) 

p1 = ggplot(tmp %>% filter(pathology == "Invasive adenocarcinoma")) + 
  geom_mosaic(aes(x = product(volume_group, evo_major_group), fill = evo_major_group), na.rm=TRUE) +
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  ylab("iITH Groups")
p2 = ggplot(tmp %>% filter(pathology == "Squamous cell carcinoma")) + 
  geom_mosaic(aes(x = product(volume_group, evo_major_group), fill = evo_major_group), na.rm=TRUE) +
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  ylab("iITH Groups")
ggarrange(p1,p2)

tmp = tmp_pyrad_clin %>% filter(histology_group == "Squamous cell carcinoma")
tmp = tmp #%>% filter(FractionTrunk_Group == "High")
surv_os <- Surv(as.numeric(as.character(tmp[,'lung_specific_time']))/365, as.numeric(as.character(tmp[,'cens_lung_specific'])))
fit_os <- survfit(survplot::censor(surv_os, 5)~diameter_group, data = tmp)
makeSurvPlot(fit_os, 
                  paste0(''), 
                  legen_title = "Diameter Group", 
                  ylab = "Lung spec. surv.",legend_coord = c(0.2,0.5),colors = c()) 




tmp = tmp_pyrad_clin# %>% filter(histology_group == "Adenocarcinoma")
tmp$iith_group = ifelse(tmp$rad_ith_mean < median(tmp$rad_ith_mean), "Low", "High")
tmp$FractionTrunk_Group = ifelse(tmp$FractionTrunk < median(tmp$FractionTrunk, na.rm = T), "Low", "High")

tmp_tmp = tmp %>% filter(FractionTrunk_Group == "High")
surv_os <- Surv(as.numeric(as.character(tmp_tmp[,'lung_specific_time']))/365, as.numeric(as.character(tmp_tmp[,'cens_lung_specific'])))
fit_os <- survfit(survplot::censor(surv_os, 5)~iith_group, data = tmp_tmp)
p1 = makeSurvPlot(fit_os, 
             paste0('ITH High'), 
             legen_title = "iITH Group", 
             ylab = "Lung spec. surv.",legend_coord = c(0.2,0.5),colors = c()) 

tmp_tmp = tmp %>% filter(FractionTrunk_Group == "Low")
surv_os <- Surv(as.numeric(as.character(tmp_tmp[,'lung_specific_time']))/365, as.numeric(as.character(tmp_tmp[,'cens_lung_specific'])))
fit_os <- survfit(survplot::censor(surv_os, 5)~iith_group, data = tmp_tmp)
p2 = makeSurvPlot(fit_os, 
                  paste0('ITH Low'), 
                  legen_title = "iITH Group", 
                  ylab = "Lung spec. surv.",legend_coord = c(0.2,0.5),colors = c()) 

plots = list()
plots[[1]] = p1
plots[[2]] = p2
pdf(file = "~/Desktop/iITH_Surv.pdf",  width = 6, height = 10)
arrange_ggsurvplots(plots, ncol = 1, nrow = 2)
dev.off()

saveRDS(object = plots, file = "~/Desktop/list_of_two_surv_plots.rds")




# Sanchez Vega mutations

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

tmp_pyrad_mut = merge(tmp_pyrad, results, by.y = "SampleID", by.x = "sampleid", all.x = T)


tmp = tmp_pyrad_mut %>% filter(pathology == "Squamous cell carcinoma")
for (col in unique(vega$Pathway)){
  if(length(unique(tmp[,col])) >= 2){
    
    test = fisher.test(table(tmp$quantcut_group, tmp[,col])) 
    if (test$p.value <= 0.08){
      print(col)
      print(test)
    }
    
  }
}



median_vol = median(tmp_pyrad_mut$volume)
median_radITH = median(tmp_pyrad_mut$rad_ith_mean, na.rm = T)
tmp_pyrad_mut$volume_group_med = ifelse(tmp_pyrad_mut$volume <= median_vol, "Low", "High")
tmp_pyrad_mut$ith_group_med = ifelse(tmp_pyrad_mut$rad_ith_mean <= median_radITH, "Low", "High")
tmp_pyrad_mut$iith_volume_group = paste0("Volume: ",tmp_pyrad_mut$volume_group_med, "- iITH: ",tmp_pyrad_mut$ith_group_med )



tmp = tmp_pyrad_mut %>% filter(pathology == "Squamous cell carcinoma", !is.na(hippo))
p1 = ggplot(tmp, 
       aes(x = rad_ith_mean,  y = volume, color = as.factor(hippo))) +
  geom_point() + 
  scale_y_log10()+
  #scale_x_log10()+
  geom_hline(yintercept = median_vol, linetype='dotted', col = 'red') + 
  geom_vline(xintercept = median_radITH, linetype='dotted', col = 'red') + 
  xlab("iITH(sum)") +
  ylab("Log10(Volume)") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") 
p2 = ggplot(data = tmp) + 
  geom_mosaic(aes(x=product(hippo, volume_group), fill = volume_group)) + 
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") 

ggarrange(p1,p2, nrow = 1)


tmp_pyrad = merge(tmp_pyrad, muts %>% 
                    group_by(SampleID) %>% 
                    summarise(tmb = n()),by.x="sampleid",by.y="SampleID",all.x = T)
ggplot(tmp_pyrad , aes(x = volume_group,  y = tmb, fill = volume_group )) +
  geom_boxplot() +
  ylab("log10(# of driver muts)") +
  xlab("Vol") +
  scale_y_log10() +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 
ggplot(tmp_pyrad , aes(x = volume,  y = tmb)) +
  geom_point() +
  ylab("log10(# of driver muts)") +
  xlab("Log10(Volume)") +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal() +
  stat_cor() +
  geom_smooth(method = "lm") +
  scale_fill_brewer(palette="Set1") 





suppressMessages(library(survminer))
suppressMessages(library(survplot))

clin = read_excel("~/Downloads/TRACERx all patients- July 2020 version 2.xlsx")
clin$ID = paste0(substr(clin$Shorter_ID,0,3), "0",substr(clin$Shorter_ID,4, 8))
tmp_pyrad_clin = merge(tmp_pyrad, clin, by.y = "ID", by.x = "sampleid", all.x = T )
tmp_pyrad_clin$vol_iith_group = paste0(tmp_pyrad_clin$diameter_group, "-", tmp_pyrad_clin$volume_group)

tmp_pyrad_clin %>% group_by(vol_iith_group) %>% summarise(n = n())
tmp = tmp_pyrad_clin %>% filter(histology_group == "Adenocarcinoma")
surv_os <- Surv(as.numeric(as.character(tmp[,'lung_specific_time']))/365, as.numeric(as.character(tmp[,'cens_lung_specific'])))
fit_os <- survfit(survplot::censor(surv_os, 5)~vol_iith_group, data = tmp)
makeSurvPlot(fit_os, 
             paste0('Adeno'), 
             legen_title = "Vol-radITH groups", 
             ylab = "Lung spec. surv.",legend_coord = c(0.2,0.5),colors = c()) 

surv_fit = coxph(surv_os ~  rad_ith_mean + smoking_status + histology_group + NSCLCstage + age
                 , data = tmp, robust = T, model = T)
surv_forest_plot(surv_fit, title = "", n = nrow(tmp))


tmp_pyrad_mut = merge(tmp_pyrad, results, by.y = "SampleID", by.x = "sampleid")


# Heatmap of features
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
hist = "Invasive adenocarcinoma"
# Squamous cell carcinoma


tmp = tmp_pyrad_mut %>% arrange(rad_ith_mean) %>% filter(pathology == hist)
ann_df <- data.frame(Volume = log10(as.numeric(tmp$original_shape_MeshVolume)),
                     iITH_Group = as.factor(tmp$quantcut_group),
                     pi3k = as.factor(tmp$pi3k),
                     cell_cycle = as.factor(tmp$cell_cycle),
                     nrf2 = as.factor(tmp$nrf2))
rownames(ann_df) <- tmp$sampleid
for_heatmap <- as.data.frame(t(tmp[,features_of_interest]))
colnames(for_heatmap) <- tmp$sampleid
pheatmap(for_heatmap,
         scale = "row",
         annotation_col = ann_df,
         cluster_rows = F,
         cluster_cols = F,
         main = hist,
         col = coul,show_colnames = F
)


tracerx_exp <- readRDS("../../Datasets/TRACERx/RNAseq/log_tpm_summarized_by_sample_exp.rds")
colnames(tracerx_exp) = gsub(substring(colnames(tracerx_exp),3,8),pattern = "LTX", replacement = "LTX0")


hallmarks = hallmarkEnrichment(tracerx_exp,method = "ssgsea") # “gsva”, “ssgsea”, “zscore”, “plage” 
hallmarks_m = reshape2::melt(hallmarks) 
hallmarks_m = merge(hallmarks_m, tmp_pyrad_mut %>% select("sampleid","diameter","rad_ith_mean", "pathology", "volume"), by.y = "sampleid",by.x = "SampleID", all.x = T )
hallmarks_m$volume_group_med = ifelse(hallmarks_m$volume <= median(hallmarks_m$volume, na.rm = T), "Low", "High")
hallmarks_m$ith_group_med = ifelse(hallmarks_m$rad_ith_mean <= median(hallmarks_m$rad_ith_mean, na.rm = T), "Low", "High")
hallmarks_m$iith_volume_group = paste0("Volume: ",hallmarks_m$volume_group_med, "- iITH: ",hallmarks_m$ith_group_med )

first_half = colnames(hallmarks)[1:25]
second_half = colnames(hallmarks)[25:50]


ggplot(hallmarks_m %>% filter( pathology == "Invasive adenocarcinoma", variable %in% second_half ),
       aes(x = variable, y = value, fill = iith_volume_group)) + 
  geom_boxplot()+ 
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  stat_compare_means(hide.ns = T, label = "p.signif") + 
  theme(axis.text.x = element_text(angle = 90))

ggplot(hallmarks_m %>% filter( pathology == "Squamous cell carcinoma", variable %in% second_half ),
       aes(x = diameter, y = value)) + 
  geom_point()+ 
  theme_pubclean() +  
  scale_fill_brewer(palette="Set1") +
  stat_cor() +
  geom_smooth(method = "lm") +
  facet_wrap(~variable)




hist = "Invasive adenocarcinoma"
tmp = tmp_pyrad_mut %>% filter(pathology == hist)
ann_df <- data.frame(Volume = log10(as.numeric(tmp$volume)),
                     iITH_Group = as.factor(tmp$quantcut_group),
                     iITH = as.numeric(tmp$rad_ith_mean))
rownames(ann_df) <- tmp$sampleid
for_heatmap = tmp_pyrad_mut %>% 
  filter(pathology == hist) %>% 
  select(c(vega$Pathway, "rad_ith_mean")) %>% 
  select("rad_ith_mean","pi3k","cell_cycle","nrf2", everything()) %>% 
  arrange(desc(nrf2), desc(cell_cycle), desc(pi3k)) %>%
  select(c(vega$Pathway)) %>% 
  select("pi3k","cell_cycle","nrf2", everything()) %>% 
  t()
colnames(for_heatmap) <- tmp$sampleid
pheatmap(for_heatmap,
         scale = "row",
         annotation_col = ann_df,
         cluster_rows = F,
         cluster_cols = F,
         main = hist,
         cellheight = 8,
         cellwidth = 4,
         fontsize = 7,
         color = c("#FEE38C", "#B10026"),
         breaks = c(0, 0.5, 1),
         show_colnames = F
)


# Expression Analysis 
# GEP

tracerx_exp_df = as.data.frame(t(tracerx_exp))
"HLA-DQA1" %in% colnames(tracerx_exp_df)
colnames(tracerx_exp_df)[startsWith(prefix = "HLA",x = colnames(tracerx_exp_df))]

# GEP signature
# This is the list of genes and their coefficient. They should be multiplied with normalised expression and summed up to get the GEP score
gep_genes <- unlist(list(CCL5=0.008346, CD27=0.072293, CD274=0.042853, CD276=-0.0239, CD8A=0.031021, CMKLR1=0.151253, CXCL9=0.074135, CXCR6=0.004313, `HLA-DQA1`=0.020091, `HLA-DRB1`=0.058806, `HLA-E`=0.07175, IDO1=0.060679, LAG3=0.123895, NKG7=0.075524, PDCD1LG2=0.003734, PSMB10=0.032999, STAT1=0.250229, TIGIT=0.084767))

score = 0
for (g in names(gep_genes)){
  print(g)
  # Get the gene "g" by column name from DF and multiply by coef. from list
  score = score + tracerx_exp_df[,g] * gep_genes[g][[1]]
}
tracerx_exp_df$gep_score = score
tracerx_exp_df$sampleid = rownames(tracerx_exp_df)

tracerx_exp_df[,c("sampleid","gep_score")]
tmp_pyrad = merge(tmp_pyrad,tracerx_exp_df %>% select("sampleid","gep_score"), by = "sampleid",  all.x = T)
ggplot(tmp_pyrad, 
       aes(x = volume,  y = gep_score)) +
  geom_point() + 
  scale_y_log10()+
  scale_x_log10()+
  xlab("Vol") +
  ylab("GEP") +
  stat_cor() + 
  geom_smooth(method="lm") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") 

ggplot(tmp_pyrad , aes(x = volume_group,  y = gep_score, fill = quantcut_group )) +
  geom_boxplot() +
  ylab("GEP") +
  xlab("iITH") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 


# TILs

tracerx_exp_df = as.data.frame(t(tracerx_exp))
TIL <- calculateTIL(tracerx_exp,immune_genes = immune_genes, dataType = "TPM")
TIL$sampleid = rownames(TIL)
tmp_pyrad = merge(tmp_pyrad,TIL, by = "sampleid", all.x = T)
ggplot(tmp_pyrad, 
       aes(x = diameter,  y = total.til.score.danahaer)) +
  geom_point() + 
 scale_x_log10()+
  ylab("TIL") +
  stat_cor() + 
  geom_smooth(method="lm") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") 

ggplot(tmp_pyrad , aes(x = quantcut_group,  y = total.til.score.danahaer, fill = quantcut_group )) +
  geom_boxplot() +
  ylab("TIL") +
  xlab("iITH") +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 


# Extra info
extra_info$sampleid = gsub(as.character(extra_info$patient_name), pattern = "LTX", replacement = "LTX0")
tmp_pyrad = merge(tmp_pyrad, extra_info,by="sampleid",all.x = T)
tmp_pyrad$purity = as.numeric(as.character(tmp_pyrad$purity))
tmp_pyrad$mean_tumour_pliody = as.numeric(as.character(tmp_pyrad$mean_tumour_pliody))
tmp_pyrad$mean_clonal_vaf = as.numeric(as.character(tmp_pyrad$mean_clonal_vaf))
ggplot(tmp_pyrad , aes(x = volume,  y = purity  )) +
  geom_point() +
  ylab("Purity") +
  xlab("Vol") +
  geom_smooth(method="lm")+
  scale_x_log10() +
  stat_cor()+  theme_minimal() +
  scale_fill_brewer(palette="Set1") 
ggplot(tmp_pyrad , aes(x = quantcut_group,  y = mean_clonal_vaf, fill = quantcut_group )) +
  geom_boxplot() +
  ylab("log10(Mean Clonal VAF)") +
  xlab("iITH") +
  scale_y_log10() +
  theme_minimal() +
  stat_compare_means() +
  scale_fill_brewer(palette="Set1") 

ggplot(tmp_pyrad , aes(x = volume,  y = mean_clonal_vaf)) +
  geom_point() +
  ylab("(Mean Clonal VAF)") +
  xlab("Volume") +
  scale_x_log10() +
  theme_minimal() +
  geom_smooth(method = "lm") + 
  stat_cor() +
  scale_fill_brewer(palette="Set1") 





mean_of_exp = mean(as.matrix(tracerx_exp))
sds_means = data.frame(sds = GMCM:::rowSds(as.matrix(tracerx_exp)), means = rowMeans(tracerx_exp))
sds_means_filtered = sds_means %>% filter(sds > 1, means > mean_of_exp)




tracerx_exp_df = as.data.frame(t(tracerx_exp[which(rownames(tracerx_exp) %in% rownames(sds_means_filtered)),]))
tracerx_exp_df$sampleid = rownames(tracerx_exp_df)



genes = colnames(tracerx_exp_df)[2:length(colnames(tracerx_exp_df))-1]
tracerx_exp_df = merge(tracerx_exp_df,tmp_pyrad_mut %>% dplyr::select("sampleid","rad_ith_mean"), by = "sampleid")

adeno_ids = tmp_pyrad_mut %>% filter(pathology == "Invasive adenocarcinoma") %>% select("sampleid") %>% c()
squamous_ids = tmp_pyrad_mut %>% filter(pathology == "Squamous cell carcinoma") %>% select("sampleid") %>% c()

tmp = tracerx_exp_df #%>% filter(sampleid %in% squamous_ids)
cor_res = data.frame()
for (gene in genes){
  cor = cor.test(tmp[,"rad_ith_mean"],tmp[,gene])
  cor_res = rbind(cor_res, data.frame(p = cor$p.value, cor = cor$estimate[[1]], name = gene))
}

cor_res = cor_res %>% 
  mutate(p_adj = p.adjust(p, method = "bonferroni")) %>% 
  filter(p_adj <= 0.05)


tracerx_exp_df_heatmap = tracerx_exp_df[,cor_res$name]
genes = colnames(tracerx_exp_df_heatmap)[2:length(colnames(tracerx_exp_df_heatmap))]
rownames(tracerx_exp_df_heatmap) = tracerx_exp_df$sampleid
tracerx_exp_df_heatmap$sampleid =  tracerx_exp_df$sampleid


hist = "Invasive adenocarcinoma"
tmp = tmp_pyrad_mut %>% dplyr::filter(pathology == hist)
tmp = merge(tmp, tracerx_exp_df_heatmap, by  = "sampleid")
ann_df <- data.frame(Volume = log10(as.numeric(tmp$original_shape_MeshVolume)),
                     iITH_Group = as.factor(tmp$quantcut_group),
                     iITH = as.numeric(tmp$rad_ith_mean),
                     Histology = tmp$pathology)
rownames(ann_df) <- tmp$sampleid
for_heatmap = tmp %>% dplyr::arrange(rad_ith_mean) %>% dplyr::select(genes) %>% t()
colnames(for_heatmap) <- tmp$sampleid
pheatmap(for_heatmap,
         scale = "row",
         annotation_col = ann_df,
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F
)




### FGSVA

library(org.Hs.eg.db)

tmp <- cor(tracerx_exp_df[,genes],tracerx_exp_df[,'rad_ith_mean'],method='spearman')[,1]
a<- map2entrez(names(tmp))
tmp <- setNames(as.numeric(VAT_SD_Mean[,2]),(VAT_SD_Mean[,1]))
pathways <- reactomePathways(VAT_SD_Mean[,1])
fgseaRes <- fgsea(pathways, tmp)
head(fgseaRes)

plotEnrichment(pathways[["Immune System"]],
               tmp) + labs(title="Immune System")



topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
fig <- plotGseaTable(pathways[topPathways], tmp, fgseaRes, gseaParam=0.5,render=T)















