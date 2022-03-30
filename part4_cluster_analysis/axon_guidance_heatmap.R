library(tidyverse)
library(readxl)
library(limma)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
setwd('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_v2')

# plotting clusters for STR axon full proteins
protein_cluster <- read_tsv("part4_cluster_analysis/rawd/20220201_STRaxon_full_proteome.timeCourse_ALL.tsv", col_names = TRUE) %>%
  mutate(mean_neonate = (plex1.neonate1 + plex1.neonate2 + plex2.neonate3 + plex2.neonate4 + plex2.neonate5)/5) %>%
  mutate(mean_earlypostnatal = (plex1.earlypostnatal1 + plex1.earlypostnatal2 + 
                                  plex1.earlypostnatal3 + plex2.earlypostnatal4 + 
                                  plex2.earlypostnatal5)/5) %>%
  mutate(mean_preweanling = (plex1.preweanling1 + plex1.preweanling2 + 
                               plex2.preweanling3 + plex2.preweanling4 + plex2.preweanling5)/5) %>%
  mutate(mean_adult = (plex1.adult1 + plex1.adult2 + plex1.adult3 + 
                         plex2.adult4 + plex2.adult5)/5)

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


######################################################
# maSigPro protein-level + LIMMA analysis
######################################################

targets_axon <- data.frame(c('neonate_1', 'neonate_2', 'neonate_3', 'neonate_4', 'neonate_5',
                             'earlypostnatal_1', 'earlypostnatal_2', 'earlypostnatal_3', 
                             'earlypostnatal_4', 'earlypostnatal_5',
                             'preweanling_1', 'preweanling_2', 'preweanling_3', 
                             'preweanling_4', 'preweanling_5',
                             'adult_1', 'adult_2', 'adult_3', 'adult_4', 'adult_5'))
targets_axon$Condition <- str_extract(targets_axon[,1], pattern = '^[^_]+(?=_)')
colnames(targets_axon) <- c('sample', 'condition')

# design matrix
lev <- c("neonate","earlypostnatal","preweanling","adult")
f <- factor(targets_axon$condition, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- lev
rownames(design) <- targets_axon$sample
design

# make the contrast
contrast <- makeContrasts(EarlypostnatalvsNeonate = earlypostnatal-neonate, 
                          PreweanlingvsNeonate = preweanling-neonate, 
                          AdultvsNeonate = adult-neonate, 
                          levels = design)
contrast

# do the linear model fitting
data_limma <- protein_cluster[,c(5:24)]
colnames(data_limma) <- targets_axon$sample
rownames(data_limma) <- protein_cluster$Protein

fit <- lmFit(data_limma, design)

# get the fit for the contrast of interest
fit2 <- contrasts.fit(fit, contrast)

# do the empirical Bayes moderation of the test statistic (with trended variance)
fit2 <- eBayes(fit2, trend = TRUE)

# let's see how many up and down candidates, and the top tags
tmp1 <- topTable(fit2, coef='EarlypostnatalvsNeonate',
                 sort.by = 'none', number = Inf)
tmp2 <- topTable(fit2, coef='PreweanlingvsNeonate',
                 sort.by = 'none', number = Inf)
tmp3 <- topTable(fit2, coef='AdultvsNeonate',
                 sort.by = 'none', number = Inf)

tmp1$Protein <- rownames(data_limma) #earlypostnatal
tmp2$Protein <- rownames(data_limma) #preweanling
tmp3$Protein <- rownames(data_limma) #adult


######################################################
# Kinases found in this dataset
######################################################

# Import Uniprot kinase family database
Uniprot_KinaseFam <- read_excel('part4_cluster_analysis/rawd/pkinfam_20210916.xlsx')

# kinases found in the dataset
pkinfam_in_protein_wide <- Uniprot_KinaseFam %>%
  merge(protein_cluster, by.x = 'Uniprot_accession', by.y = 'Protein')

# axon guidance pathway node IDs
mmu04360_ID <- read_excel("part4_cluster_analysis/rawd/mmu04360_keggID.xlsx", col_names = FALSE)
colnames(mmu04360_ID) <- c('KEGG', 'Gene')

mmu04360_ID$KEGG <- paste0('mmu:', mmu04360_ID$KEGG)

mmu04360_ID$Uniprot <- NA
for (i in seq_along(mmu04360_ID$Gene)){
  if(sum(grepl(mmu04360_ID$Gene[i], protein_cluster$Gene.Symbol)) == 0){next}
  tmp <- protein_cluster$Protein[grepl(mmu04360_ID$Gene[i], protein_cluster$Gene.Symbol)]
  tmp <- paste(tmp, sep = '', collapse = ';')
  mmu04360_ID$Uniprot[i] <- tmp
  rm(tmp)
}
mmu04360_ID$Detect <- !is.na(mmu04360_ID$Uniprot)


# kinases that are found in our dataset that also in the axon guidance pathway
pkinfam_in_axon <- pkinfam_in_protein_wide %>%
  dplyr::filter(`Gene Symbol` %in% mmu04360_ID$Gene)

pkinfam_in_axon_long <- pkinfam_in_axon[,c(1:28)] %>%
  pivot_longer(cols = c(8:27), 
               names_to = 'BioReplicate',
               values_to = 'Abundance') %>%
  mutate(Abundance_norm = Abundance - mean_neonate) %>%
  #mutate(Condition = str_extract(BioReplicate, pattern = '(.*)(?=(\\d))'))
  mutate(Condition = str_extract(BioReplicate, pattern = 'neonate|earlypostnatal|preweanling|adult'))

pkinfam_in_axon_norm_mean <- pkinfam_in_axon[,c(1:7,28:31)] %>%
  mutate(norm_neonate = mean_neonate - mean_neonate) %>%
  mutate(norm_earlypostnatal = mean_earlypostnatal - mean_neonate) %>%
  mutate(norm_preweanling = mean_preweanling - mean_neonate) %>%
  mutate(norm_adult = mean_adult - mean_neonate)


# complex heatmap
# with normalized abundance to neonate
col_fun <- colorRamp2(c(-2.5, 0, 2.5), c('#00296B', 'gray90', '#9D0208'))
col_fun(seq(-3,3))

Heatmap(as.matrix(pkinfam_in_axon_norm_mean[,c(12:15)]),
        cluster_rows = TRUE, cluster_columns = FALSE, 
        split = pkinfam_in_axon$kinase_family,
        row_labels = pkinfam_in_axon_norm_mean$`Gene Symbol`,
        col = col_fun)

# complex heatmap
# with mean abundance
col_fun <- colorRamp2(c(2, 8, 12), c('#00296B', 'gray90', '#9D0208'))
col_fun(seq(-3,3))

Heatmap(as.matrix(pkinfam_in_axon_norm_mean[,c(8:11)]),
        cluster_rows = TRUE, cluster_columns = FALSE, 
        split = pkinfam_in_axon_norm_mean$kinase_family,
        row_labels = pkinfam_in_axon_norm_mean$`Gene Symbol`,
        col = col_fun)

######################################################
# Phosphatases found in this dataset
######################################################

# axon guidance pathway node IDs
mmu_DEPOD_uniprot <- read_excel("part4_cluster_analysis/rawd/PPases_in_DEPOD_201906_mmu_mapped.xlsx")

DEPOD_protein_wide <- protein_cluster %>%
  dplyr::filter(Protein %in% mmu_DEPOD_uniprot$Entry) %>%
  dplyr::filter(!(`Gene Symbol` %in% pkinfam_in_axon$`Gene Symbol`))

# phosphatases found in the dataset
DEPOD_protein_in_axon_wide <- DEPOD_protein_wide %>%
  dplyr::filter(`Gene Symbol` %in% mmu04360_ID$Gene)

# with mean abundance
col_fun <- colorRamp2(c(2, 8, 12), c('#00296B', 'gray90', '#9D0208'))
col_fun(seq(-3,3))

Heatmap(as.matrix(DEPOD_protein_in_axon_wide[,c(25:28)]),
        cluster_rows = TRUE, cluster_columns = FALSE, 
        row_labels = DEPOD_protein_in_axon_wide$`Gene Symbol`,
        col = col_fun)


######################################################
# Kinases/Phosphatases found in this dataset
######################################################

kinases_axon <- pkinfam_in_axon[,c(1,2,5:7,28:31)]
colnames(kinases_axon)[1] <- 'Protein'

phosphatases_axon <- DEPOD_protein_in_axon_wide[,c(1:4,25:28)]
phosphatases_axon$kinase_family <- 'Phosphatase'

KPN_axon <- rbind(kinases_axon, phosphatases_axon) %>%
  mutate(norm_neonate = mean_neonate - mean_neonate) %>%
  mutate(norm_earlypostnatal = mean_earlypostnatal - mean_neonate) %>%
  mutate(norm_preweanling = mean_preweanling - mean_neonate) %>%
  mutate(norm_adult = mean_adult - mean_neonate)

# with mean abundance
col_fun <- colorRamp2(c(2, 8, 12), c('#00296B', 'gray90', '#9D0208'))
col_fun(seq(-3,3))

Heatmap(as.matrix(KPN_axon[,c(6:9)]),
        cluster_rows = TRUE, cluster_columns = FALSE, 
        row_labels = KPN_axon$`Gene Symbol`,
        split = KPN_axon$kinase_family,
        col = col_fun)

col_fun <- colorRamp2(c(-2.5, 0, 2.5), c('#00296B', 'gray90', '#9D0208'))
col_fun(seq(-3,3))

Heatmap(as.matrix(KPN_axon[,c(10:13)]),
        cluster_rows = TRUE, cluster_columns = FALSE, 
        row_labels = KPN_axon$`Gene Symbol`,
        split = KPN_axon$kinase_family,
        col = col_fun)

# non-KPN_axon
non_KPN_axon <- protein_cluster %>%
  dplyr::filter(`Gene Symbol` %in% mmu04360_ID$Gene) %>%
  dplyr::filter(!(Protein %in% KPN_axon$Protein)) %>%
  mutate(norm_neonate = mean_neonate - mean_neonate) %>%
  mutate(norm_earlypostnatal = mean_earlypostnatal - mean_neonate) %>%
  mutate(norm_preweanling = mean_preweanling - mean_neonate) %>%
  mutate(norm_adult = mean_adult - mean_neonate)

col_fun <- colorRamp2(c(-2.5, 0, 2.5), c('#00296B', 'gray90', '#9D0208'))
col_fun(seq(-3,3))

Heatmap(as.matrix(non_KPN_axon[,c(29:32)]),
        cluster_rows = TRUE, cluster_columns = FALSE, 
        row_labels = non_KPN_axon$`Gene Symbol`,
        col = col_fun)


tmp <- rbind(KPN_axon[,c(4,11:13)], non_KPN_axon[,c(3,30:32)])
#write.csv(tmp, "axon_guidance_log2FC.csv")


######################################################
# Pathview
######################################################

# compute the average log2 of each condition
protein_expression_norm <- read.csv("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/KEGG_Pathway_mapping/protein_cluster_wide_pathview_input_normtoneonate.csv") %>%
  mutate(mean_neonate = (neonate1 + neonate2 + neonate3 + neonate4 + neonate5)/5) %>%
  mutate(mean_earlypostnatal = 
           (earlypostnatal1 + earlypostnatal2 + earlypostnatal3 + earlypostnatal4 + earlypostnatal5)/5) %>%
  mutate(mean_preweanling = 
           (preweanling1 + preweanling2 + preweanling3 + preweanling4 + preweanling5)/5) %>%
  mutate(mean_adult = (adult1 + adult2 + adult3 + adult4 + adult5)/5)

#colnames(protein_expression_norm) <- 
#  c('Protein',
#    'neonate_1', 'neonate_2', 'neonate_3', 'neonate_4', 'neonate_5',
#    'earlypostnatal_1', 'earlypostnatal_2', 'earlypostnatal_3', 'earlypostnatal_4', 'earlypostnatal_5',
#    'preweanling_1', 'preweanling_2', 'preweanling_3', 'preweanling_4', 'preweanling_5',
#    'adult_1', 'adult_2', 'adult_3', 'adult_4', 'adult_5')

protein_expression_KEGGid <- read.csv("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/KEGG_Pathway_mapping/protein_cluster_wide_pathview_input_KEGGidconvert.csv")

protein_expression_norm_input <- merge(protein_expression_norm[,c(1,22:25)], 
                                       protein_expression_KEGGid[,c(1:2)], by = 'Protein')
rownames(protein_expression_norm_input) <- protein_expression_norm_input$KEGGID
protein_expression_norm_input <- protein_expression_norm_input %>% 
  dplyr::select(-c('KEGGID', 'Protein')) %>% as.matrix()
colnames(protein_expression_norm_input) <- c('time_1', 'time_2', 'time_3', 'time_4')


setwd("/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/KEGG_Pathway_mapping")
pv.out <- pathview(gene.data = protein_expression_norm_input,
                   gene.idtype = 'kegg',
                   species = 'mmu',
                   pathway.id = '04810',
                   match.data = F, multi.state = T, same.layer = F, #kegg.native = F,
                   low = 'gray90', mid = 'gray50', high = 'coral3')









#####





#
pkinfam_in_axon_long %>%
  ggplot(aes(Condition, Abundance)) +
  #ggplot(aes(Condition, Abundance_norm)) +
  geom_boxplot() +
  geom_point(size=1, shape=19, position = position_jitter()) +
  #geom_boxplot(aes(color=factor(Protein))) +
  facet_wrap(~Gene.Symbol, ncol = 7) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1))

pkinfam_in_axon_long$Condition <- factor(pkinfam_in_axon_long$Condition,
                                         levels = c('neonate', 'earlypostnatal', 'preweanling', 'adult'))

pkinfam_in_axon_long %>% 
  ggplot(aes(x = factor(Condition), y = Abundance, 
             group = factor(Gene.Symbol))) +
  stat_summary(fun = mean, aes(group=factor(Gene.Symbol)), 
               geom = 'line', size=0.5, colour = 'gray') + #, alpha=0.8) +
  facet_wrap(~kinase_family) +
  #stat_summary(fun = mean, aes(group=1), geom = 'line', size=2, colour = 'coral2') +
  #geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  #coord_cartesian(ylim = c(-2,2.5)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "right") 
                     #legend.position = "none") 


######################################################
# Up and down regulated kinases
######################################################

FCcutoff <- 2.5
sig_pkinfam_earlypostnatal <- tmp1 %>%
  dplyr::filter(Protein %in% pkinfam_in_protein_wide$Uniprot_accession) %>%
  dplyr::filter(adj.P.Val < 0.01) %>%
  dplyr::filter(abs(logFC) > log2(FCcutoff))

sig_pkinfam_preweanling <- tmp2 %>%
  dplyr::filter(Protein %in% pkinfam_in_protein_wide$Uniprot_accession) %>%
  dplyr::filter(adj.P.Val < 0.01) %>%
  dplyr::filter(abs(logFC) > log2(FCcutoff))

sig_pkinfam_adult <- tmp3 %>%
  dplyr::filter(Protein %in% pkinfam_in_protein_wide$Uniprot_accession) %>%
  dplyr::filter(adj.P.Val < 0.01) %>%
  dplyr::filter(abs(logFC) > log2(FCcutoff))


#venn.diagram(
#  x = list(sig_pkinfam_earlypostnatal$Protein, 
#           sig_pkinfam_preweanling$Protein, 
#           sig_pkinfam_adult$Protein),
#  category.names = c("Earlypostnatal" , "Preweanling " , "Adult"),
#  filename = 'sig_pkinfam_venn.png',
#  output = FALSE)


######################################################
# Overlap sig kinases in PKIN interactions?
######################################################

overlap_pkinfam_sig <- intersect(sig_pkinfam_earlypostnatal$Protein, sig_pkinfam_preweanling$Protein)
overlap_pkinfam_sig <- intersect(overlap_pkinfam_sig, sig_pkinfam_adult$Protein)

load('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/Network_PHONEMES/PHONEMES_STRaxon/KS_mmu_db.RData')

bn_overlap <- data.frame(bn) %>% 
  dplyr::filter(K.AC %in% overlap_pkinfam_sig) %>% 
  group_by(K.AC) %>% tally()
  
# write.table(overlap_pkinfam_sig, file = '/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/Network_PHONEMES/PHONEMES_STRaxon/overlap_pkinfam_sig.txt')





######################################################
# Plotting trends
######################################################

plot_prot_list <- c('P39688', 'P70206', 'Q3UH93', 'Q80UG2','Q9QUR8')
plot_prot_list <- c('O35495', 'P49615', 'Q04735', 'Q3UTQ8', 'Q8K0D0')
plot_prot_list <- c('Q2NL51', 'Q9WV60')

plot_long_df <- protein_cluster_long %>%
  dplyr::filter(Protein %in% plot_prot_list)

plot_long_df %>%
  ggplot(aes(Condition, Abundance)) +
  #ggplot(aes(Condition, Abundance_norm)) +
  geom_boxplot() +
  geom_point(size=1, shape=19, position = position_jitter()) +
  #geom_boxplot(aes(color=factor(Protein))) +
  facet_wrap(~Gene.Symbol, ncol = 7) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1))


#















####

#############

pkinfam_in_protein_clust %>% 
  ggplot(aes(x = factor(Condition), y = factor(Protein), fill = Abundance_norm)) +
  geom_tile() + 
  scale_fill_gradient2(low="blue", mid='gray90', high="red") +
  #stat_summary(fun = mean, aes(group=factor(Protein)), geom = 'tile') +
  #stat_summary(fun = mean, aes(group=factor(Protein)), geom = 'tile', size=0.1, colour = 'gray', alpha=0.8) +
  facet_wrap(~kinase_family) +
  #stat_summary(fun = mean, aes(group=1), geom = 'line', size=2, colour = 'coral2') +
  #geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  #coord_cartesian(ylim = c(-2,2.5)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "none") 



pkinfam_in_protein_clust %>% 
  ggplot(aes(x = factor(Condition), y = Abundance_norm, 
             group = factor(Protein))) +
  stat_summary(fun = mean, aes(group=factor(Protein)), geom = 'line', size=0.1, colour = 'gray', alpha=0.8) +
  #facet_wrap(~cluster) +
  #stat_summary(fun = mean, aes(group=1), geom = 'line', size=2, colour = 'coral2') +
  geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  #coord_cartesian(ylim = c(-2,2.5)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "none") 

pkinfam_in_protein_clust


## for top kinases in KEA that detected in our dataset
topKEA <- c('AKT1', 'RPS6KA1', 'MAPK1', 'MAP3K11', 'DAPK1', 'FYN', 'MAP4K5', 'STYK1', 'MAP3K5', 'MET',
            'CDK2', 'AURKA', 'PRKCI', 'STK11', 'DYRK2', 'FYN', 'CDK5', 'PRKCA', 'PRKACA')

Uniprot_KinaseFam_topKEA <- Uniprot_KinaseFam %>%
  filter(Entry_name %in% unique(topKEA))

pkinfam_in_topKEA_wide <- Uniprot_KinaseFam_topKEA %>%
  merge(protein_cluster, by.x = 'Uniprot_accession', by.y = 'Protein')


pkinfam_in_topKEA_clust <- protein_cluster_long %>%
  filter(Protein %in% Uniprot_KinaseFam_topKEA$`Uniprot_accession`) %>%
  merge(Uniprot_KinaseFam_topKEA, by.x = 'Protein', by.y = 'Uniprot_accession')


pkinfam_in_topKEA_clust %>% 
  ggplot(aes(x = factor(Condition), y = Abundance_norm, 
             group = factor(Protein))) +
  stat_summary(fun = mean, aes(group=factor(Protein), colour=factor(Gene.Symbol)), geom = 'line', size=0.75) +
  #facet_wrap(~cluster) +
  #stat_summary(fun = mean, aes(group=1), geom = 'line', size=2, colour = 'coral2') +
  geom_hline(yintercept = 0, colour = 'dodgerblue3', linetype=3, size=0.5) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1))

