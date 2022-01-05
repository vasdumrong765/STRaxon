library(tidyverse)
library(ggpubr) #correlation
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(rmotifx)
#library(ggseqlogo)

##########################################
# Import and filter peptides for protein with complete cases
##########################################
# Anastasia used 20210811_phosphopeptide_STRaxon_full.txt as input for before unrolling 1108 phosphoproteins already intersected with proteins

# 2686 proteins with missing values
# 1864 phosphoproteins
# 1108 phosphoproteins in proteins
#phosphopeptide_STRaxon_full <- raw.phospeptide.pd_norm %>%
#  dplyr::filter(`Master Protein Accessions` %in% STRaxon_full_proteome$Protein)

#phosphopeptide_raw <- read.table('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_heatmap/20210811_phosphopeptide_STRaxon_full.txt', stringsAsFactors = FALSE, header = TRUE)
#length(unique(phosphopeptide_raw$Master.Protein.Accessions))


# Unrolled phosphopeptides from Anastasia
# Later used in Corr analysis without removing NAs
phosphopeptide_unNorm_reformat <- read_tsv('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_heatmap/20210811_phosphopeptide_STRaxon_full.reformatted.tsv', col_names = TRUE) %>%
  mutate(feature = paste(modPep, `Master Protein Accessions`, sep = '.'))

#1035 phosphoproteins used in Corr analysis
#4433 modPep.prot or phosphopeptides used in Corr analysis
length(unique(phosphopeptide_unNorm_reformat$`Master Protein Accessions`))
length(unique(phosphopeptide_unNorm_reformat$feature))

#phosphopeptide_unNorm_reformat <- read_tsv('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_heatmap/20210811_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro_3sept21.tsv', col_names = TRUE)
#length(unique(phosphopeptide_unNorm_reformat$prot))

#make sure that all 1035 phosphoproteins are in fact in STRaxon_full_proteome
#sum(unique(phosphopeptide_unNorm_reformat$`Master Protein Accessions`) %in% STRaxon_full_proteome$Protein)

# plotting clusters for STR axon full proteins
protein_cluster <- read_tsv("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/clusters_analysis/rawd/20210811_STRaxon_full_proteome.timeCourse_8clusterMembership_maSigPro.tsv", col_names = TRUE) %>%
  mutate(mean_neonate = (plex1.neonate1 + plex1.neonate2 + plex2.neonate3 + plex2.neonate4 + plex2.neonate5)/5)

#only consider phosphoproteins in complete 2274
#only consider phosphopeptides with complete cases
phosphopeptide_unNorm_reformat_in2274 <- phosphopeptide_unNorm_reformat %>%
  filter(`Master Protein Accessions` %in% protein_cluster$Protein)
vars <- colnames(phosphopeptide_unNorm_reformat_in2274)[6:25]
phosphopeptide_unNorm_reformat_in2274 <- phosphopeptide_unNorm_reformat_in2274 %>% drop_na(any_of(vars))

# 707 phosphoproteins
# 2374 unique phosphopeptides
length(unique(phosphopeptide_unNorm_reformat_in2274$`Master Protein Accessions`))
length(unique(phosphopeptide_unNorm_reformat_in2274$feature))


##########################################
# Import and filter peptides for protein with complete cases
##########################################

# phosphopeptide-vs-protein correlation input
# directly from Anastasia
# use all peptides with missing values from out/20210811_phosphopeptide_STRaxon_full.reformatted.tsv
phospho_corr_spearman <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_heatmap/20210811_phosphopeptide_STRaxon_full.spearmanCorr.txt", header = TRUE)

# how many unique proteins: 1006
length(unique(phospho_corr_spearman$protein))

# only consider peptides without missing values in proteins without missing values
# 707 phosphoproteins
phospho_corr_spearman_mapped <- phospho_corr_spearman %>%
  filter(protein %in% phosphopeptide_unNorm_reformat_in2274$`Master Protein Accessions`) %>%
  mutate(modPep.prot = paste(pep, protein, sep='.')) #%>%
  #filter(modPep.prot %in% phosphopeptide_unique_residue$modPep.prot)

length(unique(phospho_corr_spearman_mapped$modPep.prot))
#length(unique(phospho_corr_spearman_mapped$protein))


##########################################
# Phosphopeptide tally and correlation heatmap
##########################################

# find the number of maximum peptides for protein in 707 phosphoproteins
phosphopep_tally <- phospho_corr_spearman_mapped %>%
  group_by(protein) %>% tally() %>% arrange(desc(n))

max(phosphopep_tally$n)
median(phosphopep_tally$n)

#phosphopep_tally %>%
#  ggplot(aes(n)) +
#  geom_histogram(binwidth = 1) +
#  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_blank())

phosphopep_tally %>%
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(0,50)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())


# create a matrix for heatmap upto 98 peptide per protein
# matrix (m x n)
heatmap_matrix <- matrix(0, nrow = 707, ncol=98)

# get a vector of correlation values for each protein
rank_corr <- function(i){
  # for protein i and peptide j
  # sort correlation
  corr_vector <- phospho_corr_spearman_mapped[phospho_corr_spearman_mapped == phosphopep_tally$protein[i],'corr'] %>%
    sort(decreasing = FALSE)
  
  # assign values to matrix m x n for protein i and peptide j
  for (j in c(1:phosphopep_tally$n[i])) {
    heatmap_matrix[i,j] <- corr_vector[j]
  }
  return(heatmap_matrix)
}

# loop through each protein i 
for (i in seq_along(phosphopep_tally$protein)){
  heatmap_matrix <- rank_corr(i)
}
rownames(heatmap_matrix) <- phosphopep_tally$protein
colnames(heatmap_matrix) <- seq(1:dim(heatmap_matrix)[2])


####### heatmap with ggplot2
# long format for heatmap data
#heatmap_matrix_long <- as.data.frame(heatmap_matrix) %>% rownames_to_column(var = "X") %>%
#  pivot_longer(cols = -c("X"),
#               names_to = "Y",
#               values_to = "Z")

# make sure to rank X by max number of peptides
#heatmap_matrix_long$X <- factor(heatmap_matrix_long$X,
#                                   levels=row.names(heatmap_matrix))

#heatmap_matrix_long %>%
#  ggplot(aes(x=X, y=as.numeric(Y), fill = Z)) +
#  geom_tile(color = "white") +
#  scale_fill_gradient2(low = 'red', high = 'black', mid = 'white') +
#  #coord_cartesian(ylim = c(0,15)) +
#  xlab('Protein') +
#  ylab('Phosphopeptide') +
#  labs(fill = "Spearman Corr") +
#  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_blank(),
#                     axis.text.x = element_blank())


####### heatmap with Complex Heatmap
col_fun <- colorRamp2(c(1, 0, -1), c('black', 'white', 'red'))
col_fun(seq(-3,3))

# upto 25 peptides
Heatmap(t(heatmap_matrix[,1:25]), cluster_rows = TREU, cluster_columns = FALSE, 
        split = 
        show_column_names = FALSE,show_row_names = TRUE, use_raster = TRUE, col=col_fun)
Heatmap(protein_cluster_long_wide[,-c(1:3)], cluster_rows = TRUE, cluster_columns = FALSE, 
        split = protein_cluster_long_wide$cluster, col=col_fun)
# full heatmap
# Heatmap(t(heatmap_matrix), cluster_rows = FALSE, cluster_columns = FALSE, 
#        show_column_names = FALSE,show_row_names = FALSE, use_raster = TRUE, col=col_fun)



##########################################
# Selecting decorrelated peptides for manual investigations related to ASD
##########################################

# annotate which residues are being modified
# for (i in seq_along(phosphopeptide_unique_residue$modPep)){
#  if (phosphopeptide_unique_residue$modPep.prot[i] %in% phospho_corr_spearman_mapped$modPep.prot){
#    location <- which(phosphopeptide_unique_residue$modPep.prot[i] == phospho_corr_spearman_mapped$modPep.prot)
#    phospho_corr_spearman_mapped$modPep.prot.res[location] <- phosphopeptide_unique_residue$modPep.prot.res[i]
#    rm(location)
#  } else {
#    print('peptide not found')
#    i = i+1
#  }
#}

# filter for decorrelation with FDR 0.1
decorrelated_pep <- phospho_corr_spearman_mapped %>%
  filter(corr < 0) %>%
  #filter(corr >= 0) %>%
  filter(fdr < 0.1)


# top decorr peptides
decorrelated_pep_top <- decorrelated_pep %>%
  arrange(fdr)




# import 172 ASD proteins
#protein_cluster_ASD_tally_Uniprot <- read_excel('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/clusters_analysis/ASD_tally_mapped.xlsx')

# are any of the decorrelated peptides in ASD proteins
# 7 phosphopeptides
#decorrelated_pep_ASD <- decorrelated_pep %>%
#  filter(protein %in% protein_cluster_ASD_tally_Uniprot$Entry)

#phospho_protein_decor_ASD <- unique(decorrelated_pep_ASD$protein)

#write.csv(decorrelated_pep_ASD, '/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_heatmap/decorrelated_pep_ASD.csv')



##########################################
# Motif analysis
##########################################
#############################
# foreground data
#############################
# Phosphopep sequence from MS/MS was prealigned using a python script with 's/t' as a centered residue of 15aa
# For rmotifx convert all to uppercase
phosphopeptide_unNorm_reformat_motif <- read_tsv('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_analysis_motif/20210811_phosphopeptide_STRaxon_full.reformatted_withmotif.tsv', col_names = TRUE)
phosphopeptide_unNorm_reformat_motif <- phosphopeptide_unNorm_reformat_motif[,-38] 
phosphopeptide_unNorm_reformat_motif <- phosphopeptide_unNorm_reformat_motif %>%
  mutate(phosY = ifelse(grepl(pattern = '\\[Y', modStr), 1, 0)) %>%
  mutate(modPep.prot = paste(modPep, `Master Protein Accessions`, sep = '.'))

fg.seqs <- toupper(phosphopeptide_unNorm_reformat_motif$modPep_motif_Y)

#############################
# Background data 1
# PTMsig.db for mouse v1.9.0.gmt
#############################
bg.database <-  read_excel('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_analysis_motif/db_ptm.sig.db.all.v1.9.0/data_PTMsigDB_all_sites_v1.9.0.xlsx')

# 15aa background database all sequences
bg.seqs1 <- bg.database$site.flanking

# Find overrepresented motifs
mot1 = motifx(fg.seqs, bg.seqs1, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
head(mot1)

#ggseqlogo(data = mot1$motif, seq_type = 'aa')

#############################
# Background data 2
# PhosphoSitePlus phosphorylation site database
#############################

#use readLines to inspect file format
#bg.database.psp <- readLines('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_analysis_motif/Phosphorylation_site_dataset.gz')
bg.database.psp <- read.table('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_analysis_motif/Phosphorylation_site_dataset.gz', skip = 3, sep = '\t', header = TRUE)

bg.seqs2 <- toupper(bg.database.psp$SITE_...7_AA)

# Find overrepresented motifs
mot2 = motifx(fg.seqs, bg.seqs2, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
head(mot2)

#ggseqlogo(data = mot2$motif, seq_type = 'aa')

#############################
# Motifx for decorrelated pep
#############################

phosphopeptide_unNorm_reformat_motif_decorr <- phosphopeptide_unNorm_reformat_motif %>%
  filter(modPep.prot %in% decorrelated_pep$modPep.prot)

fg.seqs_decorr <- toupper(phosphopeptide_unNorm_reformat_motif_decorr$modPep_motif_Y)

# Find overrepresented motifs
mot1_decorr = motifx(fg.seqs_decorr, bg.seqs1, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
head(mot1_decorr)
mot2_decorr = motifx(fg.seqs_decorr, bg.seqs2, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
head(mot2_decorr)


#write.csv(phosphopeptide_unNorm_reformat_motif_decorr, '/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_analysis_motif/decorrelated_pep_10pcFDR.csv')


#############################
# Cluster Plots for Phosphopeptides
#############################

# plotting clusters for STR axon full proteins
phospho_cluster <- read_tsv("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/clusters_analysis/rawd/20210811_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro_3sept21.tsv", col_names = TRUE) %>%
  mutate(mean_neonate = (Plex1.neonate1 + Plex1.neonate2 + Plex2.neonate3 + Plex2.neonate4 + Plex2.neonate5)/5)

#4 clusters
# phospho_cluster4 <- read_tsv("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/clusters_analysis/rawd/20210811_phosphopeptide_STRaxon_full_ALL_4clusters_maSigPro_3sept21.tsv", col_names = TRUE)
# phospho_cluster <- merge(phospho_cluster[,-c(2)], phospho_cluster4[,1:2], by = 'feature')

# keep all peptides
phospho_cluster_long <- phospho_cluster %>%
  #filter(!is.na(cluster)) %>%
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


#############################
# Profile plots for decorr-peptide-protein pairs
#############################

phosphopeptide_unNorm_reformat_motif_decorr <- phosphopeptide_unNorm_reformat_motif %>%
  #filter(modPep.prot %in% decorrelated_pep$modPep.prot)
  filter(modPep.prot %in% decorrelated_pep_top$modPep.prot)
phosphopeptide_unNorm_reformat_motif_decorr_long <- phosphopeptide_unNorm_reformat_motif_decorr %>%
  mutate(mean_neonate = (Plex1.neonate1 + Plex1.neonate2 + Plex2.neonate3 + Plex2.neonate4 + Plex2.neonate5)/5) %>%
  select(c(3,6:25, 35,40,41)) %>%
  pivot_longer(cols = -c('Master Protein Accessions', 'modPep', 'modPep.prot', 'mean_neonate'),
               names_to = 'Sample',
               values_to = 'Abundance') %>%
  separate(col = Sample,
           c(NA, 'BioReplicate')) %>%
  mutate(Condition = str_extract(BioReplicate, pattern = '(.*)(?=(\\d))')) %>%
  mutate(Abundance_norm = Abundance - mean_neonate)
  

# phosphopeptides decor
phospho_cluster_long_decorr <- phosphopeptide_unNorm_reformat_motif_decorr_long
phospho_cluster_long_decorr$Condition <- factor(phospho_cluster_long_decorr$Condition,
                                         levels=c('neonate', 'earlypostnatal', 'preweanling', 'adult'))

phospho_cluster_long_decorr$BioReplicate <- factor(phospho_cluster_long_decorr$BioReplicate,
                                            levels = c('neonate1', 'neonate2', 'neonate3', 'neonate4', 'neonate5',
                                                       'earlypostnatal1', 'earlypostnatal2', 'earlypostnatal3', 
                                                       'earlypostnatal4', 'earlypostnatal5',
                                                       'preweanling1', 'preweanling2', 'preweanling3', 
                                                       'preweanling4', 'preweanling5',
                                                       'adult1', 'adult2', 'adult3', 'adult4', 'adult5'))

phospho_cluster_long_decorr %>% 
  ggplot(aes(x = factor(Condition), y = Abundance_norm, 
             group = factor(modPep))) +
  stat_summary(fun = mean, aes(group=factor(modPep)), geom = 'line', size=0.5) +
  #facet_wrap(~cluster) +
  #stat_summary(fun = mean, aes(group=1), geom = 'line', size=2, colour = 'coral2') +
  geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "none") 


#sum(unique(phospho_cluster_long_decorr$`Master Protein Accessions`) %in% unique(protein_cluster_long$Protein))



#############################
# PatternPlots for decorr-peptide-protein pairs
#############################

## patternPlots3 function
# protein maybe not be detected at all in the phosphodata. So we only use the intersection
phospho_prot_decorr_unique <- intersect(unique(phospho_cluster_long_decorr$`Master Protein Accessions`),
                                    unique(protein_cluster_long$Protein))

protein_cluster_long_decorr <- protein_cluster_long %>%
  filter(Protein %in% phospho_prot_decorr_unique)


patternPlots_decorr <- function(i){
  data1 <- protein_cluster_long_decorr %>% 
    filter(protein_cluster_long_decorr$Protein %in% phospho_prot_decorr_unique[i])
  data2 <- phospho_cluster_long_decorr %>% filter(phospho_cluster_long_decorr$`Master Protein Accessions` %in% phospho_prot_decorr_unique[i])
  
  p1 <- data1 %>% ggplot(aes(x = factor(Condition), y = Abundance_norm)) +
    stat_summary(fun = mean, aes(group=1), geom = 'line', size=2, colour = 'coral2') +
    geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
    stat_summary(data = data2, fun = mean, aes(group=modPep.prot, colour=modPep.prot), 
                 geom = 'line', size=0.5, alpha=0.6) +
    ggtitle(paste(phospho_prot_decorr_unique[i], unique(data1$Gene.Symbol), sep=', ')) +
    theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_blank(),
                       axis.text.x = element_text(angle = 45, hjust = 1))
  return(p1)
}

#setwd('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/clusters_analysis/patterns/')
#dev.new()
#pdf(file = "/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_heatmap/phospeptidePatterns_decorr_modPep.pdf", width = 6, height = 6, onefile = TRUE)
#for (i in seq_along(phospho_prot_decorr_unique)) {
#  plot(patternPlots_decorr(i))
#}
#dev.off()

length(unique(phospho_corr_spearman_mapped$protein))
length(unique(decorrelated_pep$protein))




#############################
# Cluster Plots for Phosphopeptides
#############################

# plotting clusters for STR axon full proteins
phospho_cluster <- read_tsv("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/clusters_analysis/rawd/20210811_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro_3sept21.tsv", col_names = TRUE) %>%
  mutate(mean_neonate = (Plex1.neonate1 + Plex1.neonate2 + Plex2.neonate3 + Plex2.neonate4 + Plex2.neonate5)/5)

#4 clusters
# phospho_cluster4 <- read_tsv("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/clusters_analysis/rawd/20210811_phosphopeptide_STRaxon_full_ALL_4clusters_maSigPro_3sept21.tsv", col_names = TRUE)
# phospho_cluster <- merge(phospho_cluster[,-c(2)], phospho_cluster4[,1:2], by = 'feature')


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

#write.csv(protein_cluster_long, '/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/enrichment_cluster_analysis/protein_4cluster_wide.csv')

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