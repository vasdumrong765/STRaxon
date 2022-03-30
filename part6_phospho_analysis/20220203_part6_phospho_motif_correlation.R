## part of this code was originally written by A2IDEA, modified by V.D.

library(tidyverse)
library(ggpubr) #correlation
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(rmotifx)
setwd('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_v2')


##########################################
## Motif analysis
##########################################


# Unrolled phosphopeptides from Anastasia
# Later used in Corr analysis without removing NAs
phosphopeptide_unNorm_reformat <- read_tsv('part4_cluster_analysis/rawd/20220201_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro.tsv', col_names = TRUE)

#no NA df
#774 proteins
#2462 modPep.prot
length(unique(phosphopeptide_unNorm_reformat$prot))
length(unique(phosphopeptide_unNorm_reformat$feature))

# plotting clusters for STR axon full proteins
# 2274 proteins
protein_cluster <- read_tsv("part4_cluster_analysis/rawd/20220201_STRaxon_full_proteome.timeCourse_ALL.tsv", 
                            col_names = TRUE) %>%
  mutate(mean_neonate = (plex1.neonate1 + plex1.neonate2 + plex2.neonate3 + plex2.neonate4 + plex2.neonate5)/5) %>%
  mutate(mean_adult = (plex1.adult1 + plex1.adult2 + plex1.adult3 + plex2.adult4 + plex2.adult5)/5) %>%
  mutate(adult.neonate = mean_adult - mean_neonate)



#only consider phosphoproteins in complete 2274
#only consider phosphopeptides with complete cases
phosphopeptide_unNorm_reformat_in2274 <- phosphopeptide_unNorm_reformat %>%
  filter(prot %in% protein_cluster$Protein)
vars <- colnames(phosphopeptide_unNorm_reformat_in2274)[6:25]
phosphopeptide_unNorm_reformat_in2274 <- phosphopeptide_unNorm_reformat_in2274 %>% drop_na(any_of(vars))

# 708 phosphoproteins
# 2353 unique phosphopeptides
length(unique(phosphopeptide_unNorm_reformat_in2274$prot))
length(unique(phosphopeptide_unNorm_reformat_in2274$feature))



# find the number of maximum features for protein in 707 phosphoproteins
phosphopep_tally <- phosphopeptide_unNorm_reformat_in2274 %>%
  group_by(prot) %>% tally() %>% arrange(desc(n))

max(phosphopep_tally$n)
median(phosphopep_tally$n)
phosphopep_tally %>%
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(0,30)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())


#############################
# foreground data
#############################
# Phosphopep sequence from MS/MS was prealigned using a python script with 's/t' as a centered residue of 15aa
# For rmotifx convert all to uppercase
phosphopeptide_unNorm_reformat_motif <- read_excel('part6_phospho_analysis/20220201_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro_withmotif.xlsx')

#phosphopeptide_unNorm_reformat_motif <- phosphopeptide_unNorm_reformat_motif[,-38] 
phosphopeptide_unNorm_reformat_motif <- phosphopeptide_unNorm_reformat_motif %>%
  mutate(phosY = ifelse(grepl(pattern = '\\[Y', modStr), 1, 0))

fg.seqs <- toupper(phosphopeptide_unNorm_reformat_motif$centered_motif)


#############################
# Background data 2
# PhosphoSitePlus phosphorylation site database
#############################

#use readLines to inspect file format
#bg.database.psp <- readLines('/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/phospho_analysis_motif/Phosphorylation_site_dataset.gz')
bg.database.psp <- read.table('part6_phospho_analysis/Phosphorylation_site_dataset.gz', skip = 3, sep = '\t', header = TRUE)

bg.seqs2 <- toupper(bg.database.psp$SITE_...7_AA)

# Find overrepresented motifs
mot2 = motifx(fg.seqs, bg.seqs2, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
head(mot2)

#write.csv(mot2, 'part6_phospho_analysis/motif_analysis_results.csv', row.names = F)

# plot motif analysis
mot2 %>% 
  ggplot(aes(x = motif, y = fold.increase)) +
  geom_bar(aes(fill = fold.increase, alpha = score), stat = 'identity') +
  #scale_x_discrete(limits = rev(levels(raw_reactome_webgestal_select$order))) +
  geom_text(label = mot2$motif, y = 0.1, hjust = 'left', size = 4) +
  scale_fill_gradient(low = "gray70", high = "gray80") +
  coord_flip() +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
#axis.text.x = element_text(angle = 0, hjust = 1),
#legend.position = "none") 




##########################################
## Correletion analysis 
##########################################

## Read in peptide data
pep <- read_tsv("part4_cluster_analysis/rawd/20220201_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro.tsv", col_names = TRUE) %>% dplyr::select(c('prot', 'pep', 6:25))
colnames(pep)[1] <- "Protein"
colnames(pep)[2] <- "modPep"

## Read in protein data
prot <- read_tsv("part4_cluster_analysis/rawd/20220201_STRaxon_full_proteome.timeCourse_ALL.tsv", col_names = TRUE)
prot <- prot[,c(1, 5:24)]

#pepI <- 1
#METHOD <- 'spearman'

getCorr <- function(pepI, METHOD="pearson") {
  
  protID <- pep$Protein[pepI] ## get out the parent protein for this peptide
  
  pepValues <- suppressMessages(reshape2::melt(pep[pepI, -c(1,2)])) ## get out the quant data for this peptide
  colnames(pepValues) <- c("timePoint", "pep")
  pepValues[,1] <- tolower(pepValues[,1])
  
  ## remove any NA cases from the pepValue data.frame
  if(any(is.na(pepValues[,2]))) {
    j <- which(is.na(pepValues[,2]))
    pepValues <- pepValues[-j,]
    rm(j)
  }
  
  protValues <- suppressMessages(reshape2::melt(prot[ prot$Protein == protID, -1])) ## get out quant data for parent protein
  colnames(protValues) <- c("timePoint", "prot")
  protValues[,1] <- tolower(protValues[,1])
  
  ## remove any NA cases from the protValue data.frame
  if(any(is.na(protValues[,2]))) {
    j <- which(is.na(protValues[,2]))
    protValues <- protValues[-j,]
    rm(j)
  }
  
  ## identify timepoints in common between the protein and the peptide
  common_timepoints <- intersect(pepValues[,1], protValues[,1])
  
  ## filter both data sets to only contain these common timepoints
  pepValues <- filter(pepValues, timePoint %in% common_timepoints)
  protValues <- filter(protValues, timePoint %in% common_timepoints)
  
  d <- inner_join(x=protValues, y=pepValues, by="timePoint")
  
  ## create output data.frame 
  #ret <- data.frame(protein=protID, 
  #				  pep=pep$modPep[pepI],
  #				  neonate=NA,
  #				  earlypostnatal=NA,
  #				  preweanling=NA,
  #				  adult=NA,
  #				  stringsAsFactors=F)
  
  ret <- data.frame(protein = protID,
                    pep = pep$modPep[pepI],
                    corr = NA,
                    pvalue = NA,
                    stringsAsFactors = F)
  
  if(nrow(d) > 0) {
    ## compute correlations
    ## for(tp in c("neonate", "earlypostnatal", "preweanling", "adult")) {
    #	g <- grep(tp, d$timePoint)
    #	d2 <- d[g, ]
    
    #	g <- grep(colnames(ret))
    #  ret[1, g] <- cor(x=d2$prot, y=d2$pep, method=METHOD)
    #ret$corr <- cor(x=d$prot, y=d$pep, method=METHOD)
    tmp <- cor.test(x=d$prot, y=d$pep, method=METHOD)
    ret$corr <- tmp$estimate
    ret$pvalue <- tmp$p.value
    rm(tmp)
  }
  return(ret)
}

# Iterate over each phospho-peptide (Spearman)
out <- data.frame()
for(i in 1:nrow(pep)) {
  tmp <- getCorr(i, METHOD="spearman")
  if(!is.na(tmp)){ out <- rbind(out, tmp)}
  #if( all( !is.na(tmp[1,c(3:6)]) ) ) { out <- rbind(out, tmp) }
}
rm(tmp, i)

p_adjust <- p.adjust(out$pvalue, method = 'BH')
out <- cbind(out, p_adjust)
out$pep.prot <- paste(out$pep, out$protein, sep = '.')

# include more info on mod residues
straxon_data <- read_tsv('part5_PHONEMES_pc/rawd/20220201_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro_v3.tsv')
tmp <- merge(out, straxon_data[,c(1,29,31)], by.x = 'pep.prot', by.y = 'feature')
tmp2 <- merge(tmp, protein_cluster[,c(1,27)], by.x = 'protein', by.y = 'Protein')

out <- tmp2


## Write to disk
## write.table(out, file="20220203_phosphopeptide_STRaxon_full.spearmanCorr.txt", sep="\t", row.names=F)
# rm(out)



##########################################
# Define decorrelated peptides
##########################################

# correlation histogram
out %>% ggplot(aes(x=corr)) +
  geom_histogram() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())


# filter out NA correlation
# filter for decorrelation with FDR 0.1
decorrelated_pep <- out %>%
  filter(!is.na(corr)) %>%   # 2109 correlation
  
  filter(corr < 0) %>%       # 309 negative corr
  filter(p_adjust < 0.1) %>% # 59 padjusted negative corr
  filter(adult.neonate < 0)  # 42/59 protein up, phospho down // 17/59 is the opposite


  #filter(corr >= 0) %>%      # 1800 positive corr
  #filter(p_adjust < 0.1) %>% # 1344 padjusted negative corr
  #filter(adult.neonate > 0)  # 737/1344 protein up, phospho down // 607 is the opposite

# remember it's true that all decorrelated peptides belong to 707 proteins
# but not all decorrelated peptides belong to phosphopeptide_unNorm_reformat
# some decorrelated pep comes from entry with missing values
# thus are automatically excluded from phosphopeptide_unNorm_reformat_in2274 which is a list of all phosphopeptides that without missing values
length(unique(decorrelated_pep$modPep.prot))
length(intersect(decorrelated_pep$modPep.prot, phosphopeptide_unNorm_reformat_in2274$feature))


##########################################
# Export with more details
##########################################

# more details phosphopeptide df
phosphopeptide_MORE <- read_excel('part5_PHONEMES_pc/rawd/20220201_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro_withmotif.xlsx')

phosphopeptide_MORE <- merge(phosphopeptide_MORE, protein_cluster[,c(1,3)], by.x = 'prot', by.y = 'Protein')
#write.csv(phosphopeptide_MORE, '20220217_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro_withmotif_genes.csv',row.names=F)

corr_analysis_result <- merge(out, protein_cluster[,c(1,3)], by.x = 'protein', by.y = 'Protein')
#write.csv(out, file="20220203_phosphopeptide_STRaxon_full.spearmanCorr.csv",row.names=F)


##########################################
# Disease and Decorrelated Peptides
##########################################

risk_genes_MLM <- read_excel('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_v2/part6_phospho_analysis/Data_File_S5.xlsx',sheet = 2)


tmp <- data.frame(colnames(risk_genes_MLM), c(''))
colnames(tmp) <- c('disease', 'no_of_decor')
tmp$protein <- c('')

tmp$no_of_decor[1] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$`Developmental Delay`))
tmp$no_of_decor[2] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$Epilepsy))
tmp$no_of_decor[3] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$ADHD))
tmp$no_of_decor[4] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$`Major Depressive Disorder`))
tmp$no_of_decor[5] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$`Multiple Sclerosis`))
tmp$no_of_decor[6] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$`Parkinsons Disease`))
tmp$no_of_decor[7] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$`Alzheimer's`))
tmp$no_of_decor[8] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$Glioma))
tmp$no_of_decor[9] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$`Bipolar Disorder`))
tmp$no_of_decor[10] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$`Autism Spectrum Disorders`))
tmp$no_of_decor[11] <- length(intersect(decorrelated_pep$protein, risk_genes_MLM$Schizophrenia))

tmp$protein[1] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$`Developmental Delay`), collapse=',')
tmp$protein[2] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$Epilepsy), collapse=',')
tmp$protein[3] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$ADHD), collapse=',')
tmp$protein[4] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$`Major Depressive Disorder`), collapse=',')
tmp$protein[5] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$`Multiple Sclerosis`), collapse=',')
tmp$protein[6] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$`Parkinsons Disease`), collapse=',')
tmp$protein[7] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$`Alzheimer's`), collapse=',')
tmp$protein[8] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$Glioma), collapse=',')
tmp$protein[9] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$`Bipolar Disorder`), collapse=',')
tmp$protein[10] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$`Autism Spectrum Disorders`), collapse=',')
tmp$protein[11] <- paste0(intersect(decorrelated_pep$protein, risk_genes_MLM$Schizophrenia), collapse=',')

#write.csv(tmp, 'decor_pep_prot_disease_protein_down.csv', row.names = FALSE)
all_decorrelated_pep <- out %>%
  filter(!is.na(corr)) %>%   # 2109 correlation
  filter(corr < 0) %>%       
  filter(p_adjust < 0.1)

decor_pep_in_disease <- read_excel('part6_phospho_analysis/decor_pep_prot_disease_UniprotMapped.xlsx')

all_decorrelated_pep_disease <- all_decorrelated_pep %>%
  filter(protein %in% decor_pep_in_disease$Entry)

all_decorrelated_pep_disease <- merge(all_decorrelated_pep_disease,
                                      decor_pep_in_disease,
                                      by.x = 'protein', by.y = 'Entry')
#write.csv(all_decorrelated_pep_disease, 'all_decorrelated_pep_disease.csv', row.names =F)
