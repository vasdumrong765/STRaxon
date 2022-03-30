
library(tidyverse)
library(readxl)
library(scales)
library(eulerr)
library(ComplexHeatmap)
library(circlize)
setwd('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_v2')

################################################
### Part 1: Plotting expression by clusters
################################################

# plotting clusters for STR axon full proteins
protein_cluster <- read_tsv("part4_cluster_analysis/rawd/20220201_STRaxon_full_proteome.timeCourse_ALL.tsv", col_names = TRUE) %>%
  mutate(mean_neonate = (plex1.neonate1 + plex1.neonate2 + plex2.neonate3 + plex2.neonate4 + plex2.neonate5)/5)

protein_cluster_long <- protein_cluster %>%
  filter(!is.na(cluster)) %>%
  #filter(padjust < 0.01) %>%
  select(-c('padjust')) %>%
  pivot_longer(cols = -c('Protein', 'cluster', 'Gene Symbol', 'mean_neonate'),
               names_to = 'Sample',
               values_to = 'Abundance') %>%
  separate(col = Sample,
           c(NA, 'BioReplicate')) %>%
  mutate(Condition = str_extract(BioReplicate, pattern = '(.*)(?=(\\d))')) %>%
  mutate(Abundance_norm = Abundance - mean_neonate)
  
# reorder factors
protein_cluster_long$Condition <- factor(protein_cluster_long$Condition,
                                       levels=c('neonate', 'earlypostnatal', 'preweanling', 'adult'))

protein_cluster_long$BioReplicate <- factor(protein_cluster_long$BioReplicate,
                                          levels = c('neonate1', 'neonate2', 'neonate3', 'neonate4', 'neonate5',
                                                     'earlypostnatal1', 'earlypostnatal2', 'earlypostnatal3', 
                                                     'earlypostnatal4', 'earlypostnatal5',
                                                     'preweanling1', 'preweanling2', 'preweanling3', 'preweanling4', 
                                                     'preweanling5',
                                                     'adult1', 'adult2', 'adult3', 'adult4', 'adult5'))

# collapse time points version 2
protein_cluster_long %>% 
  ggplot(aes(x = factor(Condition), y = Abundance_norm, 
             group = factor(Protein))) +
  stat_summary(fun = mean, aes(group=factor(Protein)), geom = 'line', size=0.1, colour = 'gray', alpha=0.8) +
  facet_wrap(~cluster) +
  stat_summary(fun = mean, aes(group=1), geom = 'line', size=2, colour = 'coral2') +
  geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  coord_cartesian(ylim = c(-2,2)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "none") 


# for plotting cluster heatmap
protein_cluster_long_wide <- protein_cluster_long %>%
  select(-c('Abundance', 'mean_neonate', 'Condition')) %>%
  pivot_wider(names_from = 'BioReplicate', 
              values_from = Abundance_norm) %>%
  arrange(cluster)

col_fun <- colorRamp2(c(-2.5, 0, 2.5), c('#003049', 'white', '#9D0208'))
col_fun(seq(-3,3))

Heatmap(protein_cluster_long_wide[,-c(1:3)], cluster_rows = TRUE, cluster_columns = FALSE, 
        split = protein_cluster_long_wide$cluster, col=col_fun)


##################################################
## Part 2 Qualtiy control genes
##################################################

QC_genes <- c('Dcx', 'Nefm', 'Syp', 'Rab3a', 'Slc17a7')

protein_cluster_long %>% 
  dplyr::filter(`Gene Symbol` %in% QC_genes) %>%
  ggplot(aes(x = factor(Condition), y = Abundance_norm, 
             group = factor(Protein))) +
  stat_summary(fun = mean, aes(group=factor(Protein), color = factor(`Gene Symbol`)), 
               geom = 'line', size=0.5) +
  geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  #coord_cartesian(ylim = c(-2,2.5)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "right") 

##################################################
## Part 3 Phospho clustering
##################################################

# plotting clusters for STR axon full proteins
phospho_cluster <- read_tsv("part4_cluster_analysis/rawd/20220201_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro.tsv", col_names = TRUE) %>%
  mutate(mean_neonate = (Plex1.neonate1 + Plex1.neonate2 + Plex2.neonate3 + Plex2.neonate4 + Plex2.neonate5)/5)

phospho_cluster_long <- phospho_cluster %>%
  filter(!is.na(cluster)) %>%
  #filter(padjust < 0.01) %>%
  select(-c('padjust')) %>%
  pivot_longer(cols = -c('pep', 'cluster', 'prot', 'feature', 'mean_neonate'),
               names_to = 'Sample',
               values_to = 'Abundance') %>%
  separate(col = Sample,
           c(NA, 'BioReplicate')) %>%
  mutate(Condition = str_extract(BioReplicate, pattern = '(.*)(?=(\\d))')) %>%
  mutate(Abundance_norm = Abundance - mean_neonate)

# reorder factors
phospho_cluster_long$Condition <- factor(phospho_cluster_long$Condition,
                                         levels=c('neonate', 'earlypostnatal', 'preweanling', 'adult'))

phospho_cluster_long$BioReplicate <- factor(phospho_cluster_long$BioReplicate,
                                            levels = c('neonate1', 'neonate2', 'neonate3', 'neonate4', 'neonate5',
                                                       'earlypostnatal1', 'earlypostnatal2', 'earlypostnatal3', 
                                                       'earlypostnatal4', 'earlypostnatal5',
                                                       'preweanling1', 'preweanling2', 'preweanling3', 
                                                       'preweanling4', 'preweanling5',
                                                       'adult1', 'adult2', 'adult3', 'adult4', 'adult5'))

# collapse time points version 2
phospho_cluster_long %>% 
  ggplot(aes(x = factor(Condition), y = Abundance_norm, 
             group = factor(pep))) +
  stat_summary(fun = mean, aes(group=factor(pep)), geom = 'line', size=0.1, colour = 'gray', alpha=0.8) +
  facet_wrap(~cluster) +
  stat_summary(fun = mean, aes(group=1), geom = 'line', size=2, colour = 'coral2') +
  geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  coord_cartesian(ylim = c(-1.5,1.5)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "none") 

# for plotting cluster heatmap
phospho_cluster_long_wide <- phospho_cluster_long %>%
  select(-c('Abundance', 'mean_neonate', 'Condition')) %>%
  pivot_wider(names_from = 'BioReplicate', 
              values_from = Abundance_norm) %>%
  arrange(cluster)

col_fun <- colorRamp2(c(-2.5, 0, 2.5), c('#003049', 'white', '#9D0208'))
col_fun(seq(-3,3))

Heatmap(phospho_cluster_long_wide[,-c(1:4)], cluster_rows = TRUE, cluster_columns = FALSE, 
        split = phospho_cluster_long_wide$cluster, col=col_fun)

##################################################
## Part 4 Disease risk gene analysis
##################################################

###########################################
# Webgestal ORA analysis for customed disease database
# 8-cluster build
###########################################

# import results
#raw_disease_webgestal<- read_excel('part4_cluster_analysis/webgestal_results/webgestal_riskgenes_20220202_2274bg.xlsx')
raw_disease_webgestal<- read_excel('part4_cluster_analysis/webgestal_results/20220217_combined_disease_webgestal.xlsx')


raw_disease_webgestal$FDR <- as.numeric(raw_disease_webgestal$FDR)
raw_disease_webgestal$pValue <- as.numeric(raw_disease_webgestal$pValue)
raw_disease_webgestal$symbol <- factor(raw_disease_webgestal$symbol,
                                       levels = c('Scz', 'PD', 'MS', 'MDD', 'Gli', 'EP', 'DD', 'BP', 'ASD', 'AD'))
#levels = c('AD', 'ASD', 'BP', 'DD', 'Ep', 'Gli', 'MDD', 'MS', 'PD', 'Scz'))

raw_disease_webgestal2 <- raw_disease_webgestal %>% 
  filter(!(symbol == 'Gli')) %>%
  filter(FDR < 0.3)

ggplot(raw_disease_webgestal2, aes(x = symbol, y = factor(cluster))) +
  geom_point(aes(x = symbol, y = factor(cluster), 
                 #size = as.numeric(enrichmentRatio),
                 size = as.numeric(overlap),
                 #colour = FDR)) +
                 #colour = (-1)*log10(FDR))) +
                 colour = (-1)*log10(pValue))) +  
  scale_color_gradient(low = "gray80", high = "coral3") +
  #scale_color_gradient(low = "gray80", high = "coral3", limits=c(0, NA)) +
  scale_size_continuous(range = c(3,8)) +
  scale_y_discrete(limits = as.factor(c(8:1))) +
  #scale_y_reverse() +
  theme_bw() +
  ylab("Gene cluster") +
  xlab("Disease symbol") +
  labs(color="-log10(pvalue)", size="Number\nof overlap genes")
#theme(legend.position = 'none')


###########################################
# Webgestal ORA analysis for Reactome pathways
# 8-cluster build
###########################################

# import results
raw_reactome_webgestal <- read_excel('part4_cluster_analysis/webgestal_results/webgestal_reactome_20220202_2274bg.xlsx')

raw_reactome_webgestal$order <- factor(seq_along(raw_reactome_webgestal$geneSet), ordered = TRUE)
#fct_rev(raw_reactome_webgestal$order)

# only pick top 5
toppick_vector <- c(c(1:5), c(11:15), c(21:25), c(31:35), c(41:45), c(51:55), c(61:65), c(71:75))
raw_reactome_webgestal_select <- raw_reactome_webgestal[toppick_vector,]
raw_reactome_webgestal_select$order <- droplevels(raw_reactome_webgestal_select$order)


raw_reactome_webgestal_select %>% 
  ggplot(aes(x = order, y = (-1)*log10(pValue))) +
  geom_bar(aes(fill = pValue, alpha = (-1)*pValue), stat = 'identity') +
  scale_x_discrete(limits = rev(levels(raw_reactome_webgestal_select$order))) +
  geom_text(label = raw_reactome_webgestal[toppick_vector,]$description, y = 0.1, hjust = 'left', size = 2) +
  scale_fill_gradient(low = "gray70", high = "gray80") +
  coord_flip() +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     #axis.text.x = element_blank(),
                     axis.text.x = element_text(angle = 0, hjust = 1),
                     legend.position = "none") 

###########################################
# Risk gene cluster analysis
###########################################

# ASD risk genes in cluster 2
ASD_risk_cluster2 <- raw_disease_webgestal2$userId[raw_disease_webgestal2$cluster == 2] %>% 
  str_split(pattern = ';') %>% unlist()

# Protein associated with Reactome pathways
# Pick Sodium/Calcium Exchangers, Calcineurin activates NFAT, IP3, IP4, PKA
cluster2_selected_reactome <- dplyr::filter(raw_reactome_webgestal, cluster == 2) %>%
  dplyr::filter(description %in% c('Sodium/Calcium exchangers', 'Calcineurin activates NFAT', 'Synthesis of IP3 and IP4 in the cytosol', 'Reduction of cytosolic Ca++ levels', 'FCERI mediated Ca+2 mobilization', 'PKA-mediated phosphorylation of CREB', 'PKA activation')) %>% select(userId) 
  
  str_split(pattern = ';') %>% unlist()
cluster2_selected_reactome

c('P0DP27', 'Q6P069', 'P48453', 'O08586', 'O88444')



