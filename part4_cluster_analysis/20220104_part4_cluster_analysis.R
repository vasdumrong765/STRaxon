
library(tidyverse)
library(readxl)
library(scales)
library(eulerr)
library(ComplexHeatmap)
library(circlize)
setwd('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_code')

################################################
### Part 1: Plotting expression by clusters
################################################

# plotting clusters for STR axon full proteins
protein_cluster <- read_tsv("part4_cluster_analysis/rawd/20210811_STRaxon_full_proteome.timeCourse_8clusterMembership_maSigPro.tsv", col_names = TRUE) %>%
  mutate(mean_neonate = (plex1.neonate1 + plex1.neonate2 + plex2.neonate3 + plex2.neonate4 + plex2.neonate5)/5)

protein_cluster_long <- protein_cluster %>%
  filter(!is.na(cluster)) %>%
  #filter(padjust < 0.01) %>%
  select(-c('padjust')) %>%
  pivot_longer(cols = -c('Protein', 'cluster', 'Gene.Symbol', 'mean_neonate'),
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
  coord_cartesian(ylim = c(-2,2.5)) +
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
  dplyr::filter(Gene.Symbol %in% QC_genes) %>%
  ggplot(aes(x = factor(Condition), y = Abundance_norm, 
             group = factor(Protein))) +
  stat_summary(fun = mean, aes(group=factor(Protein), color = factor(Gene.Symbol)), 
               geom = 'line', size=0.5) +
  geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  #coord_cartesian(ylim = c(-2,2.5)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "right") 

