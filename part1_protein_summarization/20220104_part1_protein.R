## Dumrongprechachan et al 2022
## Part 1 Protein Summarization

## start up library
library(tidyverse)
library(MSstatsTMT)
library(data.table)
library(ggpubr)
library(ggrepel)
library(gplots)
library(edgeR)
library(psych)
library(reader)
library(robustbase)
library(expss)
library(limma)
library(readxl)
library(RColorBrewer)
source("/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Rfunctions/med_SL_norm.R")
source("/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive//Rfunctions/plot_label.R")
library(ComplexHeatmap)
library(mclust)

setwd('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_code')

#########################
## Load data - STR axon flow through SwissProt only
#########################

# Read PSMs PD output into R
# S/N export from PD
raw.pd_frac <- read.table("part1_protein_summarization/04162021_STRaxon_flowthrough_first_injection_PSMs.txt", sep="\t", header=TRUE) %>% 
  dplyr::filter(Isolation.Interference.... < 70) %>%
  # remove entries with quant info excluded by method
  dplyr::filter(Quan.Info == '') %>%
  # 1% peptide FDR
  dplyr::filter(Percolator.q.Value < 0.01)


# Read PD protein tabs into R
# For gene name mapping
raw.protein.pd_frac <- read_excel(path="part1_protein_summarization/04162021_STRaxon_flowthrough_first_injection_Proteins.xlsx")

# MSstatsTMT annotation file into R
annotation.pd_frac <- read.table("part1_protein_summarization/PD_STRaxon_annotation.txt", sep="\t", header=TRUE)


#########################
## MSstatsTMT PSM selection
#########################

# Converting PD output to MSstats format
# Remove Proteins with 1 Feature
# Remember that MSstatsTMT only use unique peptides for protein summarization and quantifications
# Shared peptides are not currently implemented.
input.pd <- PDtoMSstatsTMTFormat(raw.pd_frac,
                                 annotation.pd_frac,
                                 which.proteinid = 'Master.Protein.Accessions',
                                 useNumProteinsColumn = FALSE,
                                 useUniquePeptide = TRUE,
                                 rmPSM_withMissing_withinRun = FALSE,
                                 rmPSM_withfewMea_withinRun = TRUE, #features with 1-2 measurements
                                 rmProtein_with1Feature = TRUE,
                                 summaryforMultipleRows = sum)
# Removed features/PSMs that are shared between multiple proteins (PSMs with multiple Uniprot entries)
input.pd <- input.pd %>% filter(!grepl(";", ProteinName))



#########################
## Data summarized with the 'LogSum' MSstatsTMT workflow without additional peptide filtering
## No global median peptide normalization prior to non-specific cutoff analysis
#########################

# protein abundance are aggregated by summing the PSM intensities then log2 transformed.
quant.msstats0 <- proteinSummarization(input.pd,
                                       method="LogSum",
                                       global_norm=FALSE,
                                       reference_norm=TRUE,
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE,                  
                                       MBimpute = FALSE,
                                       maxQuantileforCensored = NULL)

# remove negative log2 abundance after normalization
quant.msstats0 <- quant.msstats0 %>% dplyr::filter(Abundance > 0)


##################################################
## Filter 1 - remove proteins that are non-specifically enriched during pull-down based on P18 data
##################################################
# Remove proteins that are not positively enriched over the IP control background (Cre-negative samples)
# This list focuses only on the STR P18 samples

levels(droplevels(quant.msstats0$Condition))
#create a comparison matrix
comparison1 <- matrix(c(-1, 1,  0, 0, 0, 0, 0,
                         0, 0, -1, 0, 0, 1, 0),
                      nrow=2, ncol=7, byrow=TRUE)
colnames(comparison1) <- c("CTX_negative","CTX_H2B", 
                           "STR_negative", "P5", "P11", "P18", "P50")
row.names(comparison1) <- c("CTX_H2B-CTX_negative", 
                            "P18-STR_negative")
comparison1

#perform MSstats group comparison as indicated above
test.contrast1 <- groupComparisonTMT(data = quant.msstats0,
                                     contrast.matrix = comparison1,
                                     moderated = TRUE,
                                     adj.method = 'BH',
                                     remove_norm_channel = TRUE,
                                     remove_empty_channel = TRUE)

# volcano plot for filter 1
volcano_df1 <- ggplot(data = test.contrast1,
                      aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "red")) +
  geom_hline(yintercept = (-1)*log10(0.005), linetype = "dotted", size = 0.6) +
  geom_vline(xintercept = log2(2.5), linetype = "dotted", size = 0.6) +
  facet_wrap(~Label, scale = "free", ncol = 6)

volcano_df1

# impose some cutoff in P18-STR_negative
# proteins must have 2.5 fold change over the Cre-negative control in 
# must have q-value < 0.01

non_specific_enrich <- test.contrast1 %>%
  dplyr::filter(Label == 'P18-STR_negative') %>%
  dplyr::filter(adj.pvalue >= 0.005 | (log2FC < log2(2.5)))


# remarks - for protein with One Condition Missing
# is that condition missing in the Cre-negative control samples?
protein_missing_OneCondition <- test.contrast1 %>%
  dplyr::filter(Label == 'P18-STR_negative') %>%
  dplyr::filter(issue == 'oneConditionMissing')

protein_missing_OneCondition_df <- quant.msstats0 %>%
  dplyr::filter(Condition %in% c('STR_negative', 'P18')) %>%
  dplyr::filter(Protein %in% protein_missing_OneCondition$Protein)

# it looks like that One Condition was missed only in the Cre_negative control
# makes sense since the Cre_negative control should have none to background level signal in theory
unique(protein_missing_OneCondition_df$Condition)

# histogram of log2FC over control
test.contrast1 %>%  ggplot(aes(x = log2FC)) +
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~Label) +
  #coord_cartesian(xlim = c(0,50)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())



##################################################
## Global median normalization at protein-level + apply Filter 1
## Limma-based correction ofb atch effects
##################################################
# Remove highly non-specific enrichment
# Compare CTX_H2B with P18 STR to define STR enriched proteins

# protein abundance are aggregated by summing the PSM intensities then log2 transformed.
quant.msstats1 <- proteinSummarization(input.pd,
                                       method="LogSum",
                                       global_norm=TRUE,  #median normalize peptide signal
                                       reference_norm=TRUE,
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE,                  
                                       MBimpute = FALSE,
                                       maxQuantileforCensored = NULL)

tmp <- quant.msstats1 
# remove negative log2 abundance after normalization and non-specific enrichment
quant.msstats1 <- tmp %>% 
  dplyr::filter(Abundance > 0) %>%
  dplyr::filter(!(Protein %in% non_specific_enrich$Protein)) %>%
  dplyr::filter(!(Protein %in% protein_missing_OneCondition$Protein))
rm(tmp)

# remove Norm, CTX_negative, STR_negative samples
quant.msstats1_wide <- quant.msstats1 %>%
  dplyr::filter(!(Condition %in% c('Norm', 'CTX_negative', 'STR_negative'))) %>%
  dplyr::select(Protein, Abundance, BioReplicate) %>%
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)

# re-arrange columns
col_order <- c('P5_1', 'P5_2', 'P5_3', 'P5_4', 'P5_5',
               'P11_1', 'P11_2', 'P11_3', 'P11_4', 'P11_5',
               'P18_1', 'P18_2', 'P18_3', 'P18_4', 'P18_5',
               'P50_1', 'P50_2', 'P50_3', 'P50_4', 'P50_5',
               'CTX_H2B_1', 'CTX_H2B_2', 'CTX_H2B_3')

quant.msstats1_wide <- cbind(quant.msstats1_wide[,1],
                             quant.msstats1_wide[, col_order])

# look at boxplot before protein level median normalization
# color by condition
colors = c(rep('gray75',5), 
           rep('gray55',5),
           rep('gray35',5),
           rep('gray15',5),
           rep('steelblue4',3))
boxplot(quant.msstats1_wide[,-1], col = colors)

# protein-level median normalization and Limma-based batch effect correction
# MSstatsTMT Abundance already in log2

# batch information
batch_colors <- c(1, 1, 2, 2, 2,
                  1, 1, 1, 2, 2,
                  1, 1, 2, 2, 2, 
                  1, 1, 1, 2, 2,
                  1, 1, 2)

quant.msstats1_wide_med <- removeBatchEffect(x = quant.msstats1_wide[,-1],
                                             batch = batch_colors) %>% med_norm() %>% data.frame()

quant.msstats1_wide_med <- cbind(quant.msstats1_wide[,1], quant.msstats1_wide_med)
colnames(quant.msstats1_wide_med)[1] <- 'Protein'

# check boxplot after med_norm
colors2 <-  c(rep('#DD6E42',5), 
              rep('#E8DAB2',5),
              rep('#4F6D7A',5),
              rep('#C0D6DF',5),
              rep('#EAEAEA',3))
boxlabels <- c('neonate1', 'neonate2', 'neonate3', 'neonate4', 'neonate5',
               'earlypostnatal1', 'earlypostnatal2', 'earlypostnatal3', 
               'earlypostnatal4', 'earlypostnatal5',
               'preweanling1', 'preweanling2', 'preweanling3', 
               'preweanling4', 'preweanling5',
               'adult1', 'adult2', 'adult3', 'adult4', 'adult5',
               'ctx_preweanling1', 'ctx_preweanling2', 'ctx_preweanling3')
boxplot(quant.msstats1_wide_med[,-1], col = colors2, las = 2, cex.axis = 0.7)
boxplot(quant.msstats1_wide_med[,-1], col = batch_colors)


# check heatmap
# no NA for heatmap
# Looks like samples of the sample condition clustered together!
heatmap_quant.msstats1_med <- quant.msstats1_wide_med %>% drop_na()
Heatmap(matrix = as.matrix(heatmap_quant.msstats1_med[,-1]), show_row_names = FALSE, split = 8)


# ready for MSstatsTMT comparison
# convert data back to long format
quant.msstats1_wide_med_long <- pivot_longer(quant.msstats1_wide_med,
                                               cols = -Protein,
                                               names_to = 'BioReplicate',
                                               values_to = 'Abundance')
quant.msstats1_wide_med_long <- merge(quant.msstats1_wide_med_long,
                                        quant.msstats1[,-3], # remove original Abundance column
                                        by=c('Protein','BioReplicate'))

##################################################
## Filter 2 - Define axon-enriched proteins
##################################################

levels(droplevels(quant.msstats1_wide_med_long$Condition))

#create a comparison matrix
comparison2 <- matrix(c(-1, 0, 0, 1, 0),
                      nrow=1, ncol=5, byrow=TRUE)
colnames(comparison2) <- c("CTX_H2B", "P5", "P11", "P18", "P50")
row.names(comparison2) <- c("P18-CTX_H2B")
comparison2

test.contrast2 <- groupComparisonTMT(data = quant.msstats1_wide_med_long,
                                     contrast.matrix = comparison2,
                                     moderated = TRUE,
                                     adj.method = 'BH',
                                     remove_norm_channel = FALSE,
                                     remove_empty_channel = TRUE)

# volcano plot for filter 2
volcano_df2 <- ggplot(data = test.contrast2,
                      aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "red")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6)

volcano_df2

# axon enriched proteins are those with log2FC > 0 & adj.pvalue < 0.05
axon_enrich_protein <- test.contrast2 %>%
  dplyr::filter(log2FC > 0) %>%
  dplyr::filter(adj.pvalue < 0.05)


##################################################
## Normalized full dataset for Anastasia
##################################################

# rename to include plex information
# plex1_ or plex2_
batch_order <- c('plex1', 'plex1', 'plex2', 'plex2', 'plex2',
                  'plex1', 'plex1', 'plex1', 'plex2', 'plex2',
                  'plex1', 'plex1', 'plex2', 'plex2', 'plex2',
                  'plex1', 'plex1', 'plex1', 'plex2', 'plex2')

bioreplicates <- c('neonate1', 'neonate2', 'neonate3', 'neonate4', 'neonate5',
                   'earlypostnatal1', 'earlypostnatal2', 'earlypostnatal3', 'earlypostnatal4', 'earlypostnatal5',
                   'preweanling1', 'preweanling2', 'preweanling3', 'preweanling4', 'preweanling5',
                   'adult1', 'adult2', 'adult3', 'adult4', 'adult5')

STRaxon_full_proteome <- quant.msstats1_wide_med %>%
  select(-c('CTX_H2B_1', 'CTX_H2B_2', 'CTX_H2B_3'))

columnnames <- paste(batch_order, bioreplicates, sep = '.')
colnames(STRaxon_full_proteome)[-1] <- columnnames

# mapped ID information to the dataframe
mapped_ID <-  raw.protein.pd_frac %>%
  dplyr::select(Accession,
                Description,
                `# PSMs`, # PSMs associated with this protein; only unique PSMs are used in MSstats summarization
                `# Unique Peptides`,
                `Biological Process`,
                `Cellular Component`,
                `Molecular Function`,
                `Gene Symbol`)

STRaxon_full_proteome <- merge(x=STRaxon_full_proteome,
                               y=mapped_ID, by.x="Protein", by.y = "Accession")


##################################################
## CV calculations
##################################################

CV_values <- quant.msstats1_wide_med_long %>% group_by(Protein, Condition) %>%
  summarise(cond_avg = mean(Abundance, na.rm=TRUE), cond_sd = sd(Abundance, na.rm=TRUE)) %>%
  mutate(cond_cv = cond_sd*100/cond_avg)

CV_values$Condition <- factor(CV_values$Condition,
                              levels = c('P5', 'P11', 'P18', 'P50','CTX_H2B'))

my_colors <- c('#DD6E42','#E8DAB2','#4F6D7A','#C0D6DF','#EAEAEA')

# traditional boxplots with outliers removed
# only for STR axon samples
CV_values %>% 
  ggplot(aes(x = Condition, y = cond_cv)) +
  geom_boxplot(aes(fill=Condition),notch = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  ylim(c(0,40)) + ylab('% CV') +
  ggtitle("CV distributions") +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

# CV density plots
CV_values %>% ggplot(aes(x = cond_cv, color = Condition)) +
  geom_density() +
  scale_colour_manual(values = my_colors) +
  xlab('% CV') +
  coord_cartesian(xlim = c(0, 50)) +
  ggtitle("CV distributions") +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())


##################################################
## Plot DAVID GO analysis results
##################################################

# DAVID GO results
# calculate -log10(FDR)
# sort by Fold.Enrichment

GOCC_axon <- read.table("part1_protein_summarization/DAVID_20210801/enriched_CC.txt", sep = '\t', header = TRUE) %>% 
  mutate(logFDR = -log10(FDR)) %>%
  arrange(desc(logFDR)) %>%
  separate(Term, c('GO_accession', 'Term'), sep='~')

GOCC_nonaxon <-  read.table("part1_protein_summarization/DAVID_20210801/unenriched_CC.txt", sep = '\t', header = TRUE) %>% 
  mutate(logFDR = -log10(FDR)) %>%
  arrange(desc(logFDR)) %>%
  separate(Term, c('GO_accession', 'Term'), sep='~')


# top13 terms
GOCC_axon_topterm <- GOCC_axon[1:13,] %>% arrange(Fold.Enrichment)
GOCC_nonaxon_topterm <- GOCC_nonaxon[1:13,] %>% arrange(Fold.Enrichment)

# code adapted from https://bio-protocol.org/bio101/e3429
# top 13 terms
ggplot(GOCC_axon_topterm, aes(x = Term, y = Fold.Enrichment)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "azure4", size=.5) +
    geom_point(aes(x = Term, y = Fold.Enrichment,size = Count, colour = logFDR), alpha=.7) +
  scale_x_discrete(limits = GOCC_axon_topterm$Term) +
  scale_color_gradient(low = "dodgerblue2", high = "coral2", limits=c(0, NA)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5) +
  xlab("GO cellular component") +
  ylab("Fold enrichment") +
  labs(color="-log10(FDR)", size="Number\nof genes")


ggplot(GOCC_nonaxon_topterm, aes(x = Term, y = Fold.Enrichment)) +
  geom_hline(yintercept = 1, linetype="dashed", color = "azure4", size=.5) +
  geom_point(aes(x = Term, y = Fold.Enrichment,size = Count, colour = logFDR), alpha=.7) +
  scale_x_discrete(limits = GOCC_nonaxon_topterm$Term) +
  scale_color_gradient(low = "dodgerblue2", high = "coral2", limits=c(0, NA)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5) +
  xlab("GO cellular component") +
  ylab("Fold enrichment") +
  labs(color="-log10(FDR)", size="Number\nof genes")


##################################################
## Tally axonal proteins
##################################################

uniprot_mapped_axon_enriched <- read_excel("part1_protein_summarization/uniprot_20210826_STRaxon_all_tally_for_axon.xlsx") %>%
  mutate(quick_axon = grepl(pattern = 'axon', `Gene ontology (cellular component)`)) %>%
  mutate(quick_presy = grepl(pattern = 'presynap', `Gene ontology (cellular component)`)) %>%
  mutate(quick_gcone = grepl(pattern = 'growth cone', `Gene ontology (cellular component)`)) %>%
  mutate(quick_neuproj = grepl(pattern = 'neuron projection', `Gene ontology (cellular component)`)) %>%
  mutate(quick_axolemma = grepl(pattern = 'axolemma', `Gene ontology (cellular component)`)) %>%
  mutate(quick_synapv = grepl(pattern = 'synaptic vesicle', `Gene ontology (cellular component)`)) %>%
  
  mutate(all_axon = (quick_axon | quick_presy | quick_gcone | quick_neuproj | quick_axolemma |
                       quick_synapv))

# sum(uniprot_mapped_axon_enriched$all_axon)


##################################################
## volcano plot for filter 2
##################################################

test.contrast2_mapped <- merge(test.contrast2, mapped_ID,
                               by.x="Protein", by.y="Accession") %>%
  filter(!is.na(adj.pvalue)) %>%
  plot_label(FCvalue = 5,
             FDRsig = 0.05,
             FDRplot = 1e-10,
             direction = 'pos')

for (i in seq_along(test.contrast2_mapped$Protein)) {
  if (test.contrast2_mapped$Protein[i] %in% uniprot_mapped_axon_enriched$Entry){
    test.contrast2_mapped$axon_tally[i] <- 
      uniprot_mapped_axon_enriched[uniprot_mapped_axon_enriched$Entry == test.contrast2_mapped$Protein[i],]$all_axon}
  else {
    test.contrast2_mapped$axon_tally[i] <- FALSE
  }}

# sum(test.contrast2_mapped$axon_tally)

test.contrast2_mapped <- test.contrast2_mapped %>%
  mutate(axon_label_select = ifelse(axon_tally == TRUE, Protein_label_select, NA)) %>%
  mutate(colour_code = ifelse(axon_tally == TRUE, 1, ifelse(sig == 1, 2, 3))) %>%
  mutate(alpha_code = ifelse(axon_tally == TRUE, 2, ifelse(sig == 1, 1, 1)))


# simple version
test.contrast2_mapped %>% ggplot(aes(x = log2FC, y = (-1)*log10(pvalue))) +
  geom_point(aes(fill = factor(sig), colour = factor(colour_code), alpha = alpha_code), shape =21) +
  scale_fill_manual(values = c('gray','#66B2FF')) +
  scale_colour_manual(values = c('#DD6E42', '#66B2FF', 'gray')) + 
  scale_alpha(range = c(0.5,1)) +
  scale_shape_identity() +
  geom_hline(yintercept = (-1)*log10(0.05), 
             linetype = "dotted", size = 0.6) +
  ylab('q-value') +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     legend.position = 'none')


# complicated version
volcano_df2_label <- test.contrast2_mapped %>% 
  ggplot(aes(x = log2FC, y = (-1)*log10(pvalue))) +
  geom_point(aes(colour = factor(sig), fill = factor(axon_tally), shape =21)) +
  scale_colour_manual(values = c('gray', 'black')) +
  scale_fill_manual(values = c('gray', 'coral2')) +
  scale_shape_identity() +
  geom_hline(yintercept = (-1)*log10(0.05), 
             linetype = "dotted", size = 0.6) +
  ylab('q-value') +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_blank()) +
  geom_text_repel(aes(label = axon_label_select), 
                  max.overlaps = Inf,
                  force = 5, #direction = "y",
                  nudge_x = 0.3, nudge_y = 0.3,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.5, #segment line thickness
                  size=4, fontface = 'italic')
volcano_df2_label


##################################################
## Rank plots
##################################################

# rank plot
# rank based on log2FC
test.contrast2_mapped_rank <- test.contrast2_mapped %>% arrange(log2FC)
test.contrast2_mapped_rank$Rank <- seq_along(test.contrast2_mapped_rank$Protein)


# histogram tally of GOaxon
test.contrast2_mapped_rank$csum_Goaxon <- 0
test.contrast2_mapped_rank$csum_notGoaxon <- 0
test.contrast2_mapped_rank <- test.contrast2_mapped_rank %>% arrange(desc(log2FC))
test.contrast2_mapped_rank$Rank_rev <- seq_along(test.contrast2_mapped_rank$Protein)
for (i in seq_along(test.contrast2_mapped_rank$Rank)){
  # counting
  test.contrast2_mapped_rank$csum_Goaxon[i] <- sum(test.contrast2_mapped_rank$axon_tally[1:i])
  test.contrast2_mapped_rank$csum_notGoaxon[i] <- sum(!(test.contrast2_mapped_rank$axon_tally[1:i]))
  # compute 'TPR' and 'FPR'
  test.contrast2_mapped_rank$TPR[i] <- test.contrast2_mapped_rank$csum_Goaxon[i] / 605
  test.contrast2_mapped_rank$FPR[i] <- test.contrast2_mapped_rank$csum_notGoaxon[i] / 2078
  # compute TPR - FPR
  test.contrast2_mapped_rank$csum_diff[i] <- 
    test.contrast2_mapped_rank$csum_Goaxon[i] - test.contrast2_mapped_rank$csum_notGoaxon[i]
  test.contrast2_mapped_rank$rate_diff[i] <- test.contrast2_mapped_rank$TPR[i] - test.contrast2_mapped_rank$FPR[i]
}


################
# rank plot
test.contrast2_mapped_rank <- test.contrast2_mapped_rank %>% arrange(log2FC)
test.contrast2_mapped_rank %>% ggplot(mapping = aes(x = Rank_rev, y = log2FC)) + 
  geom_point(aes(colour=factor(axon_tally)), size = 1.5, show.legend = FALSE) + 
  #scale_shape_manual(values = c(1, 21)) +
  #scale_fill_manual(values = c("gray", "coral2")) +
  #scale_colour_continuous(low = "gray", high = "red") +
  scale_colour_manual(values = c("gray", "coral2")) +
  geom_hline(yintercept = 0,linetype = "dotted", size = 0.6) +
  #geom_text_repel(aes(label = Protein_label_select2), 
  #                max.overlaps = Inf,
  #                force = 3, #direction = "y",
  #                #nudge_x = 0.2, nudge_y = 0.2,
  #                min.segment.length = unit(0, 'lines'), #draw every segment
  #                segment.size = 0.2, #segment line thickness
  #                size=4,
  #                fontface = 'italic') + #font size
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())


################
# cumulative sum GOaxon and notGOaxon plots

# 4x4 figure
test.contrast2_mapped_rank %>% 
  ggplot(mapping = aes(x = Rank_rev, y = csum_Goaxon)) + 
  geom_point(colour='coral2', size=1) +
  geom_point(aes(x = Rank_rev, y = csum_notGoaxon), colour='gray', size=1) +
  #scale_y_reverse() +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

# rate of annotation 4x3
test.contrast2_mapped_rank %>% 
  ggplot(mapping = aes(x = Rank_rev, y = TPR)) + 
  geom_point(colour='coral2', size=1) +
  geom_point(aes(x = Rank_rev, y = FPR), colour='gray', size=1) +
  scale_x_reverse() +
  coord_flip() +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

test.contrast2_mapped_rank %>% 
  ggplot(mapping = aes(x = Rank_rev, y = rate_diff)) + 
  geom_point(colour='black', size=1) +
  scale_x_reverse() +
  coord_flip() +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())


STRaxon_full_proteome_noNA <- STRaxon_full_proteome %>%
  drop_na(plex1.neonate1:plex2.adult5)

#write_csv(STRaxon_full_proteome_noNA, file = "/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/FT_SwissProt_only/20211020_STRaxon_full_proteome_noNA_2276.csv")
