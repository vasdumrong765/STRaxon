## Dumrongprechachan et al 2022
## Part 2 Phosphopeptide Summarization

## start up library
library(tidyverse)
library(readxl)
library(robustbase)
library(limma)
library(ComplexHeatmap)
source("/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Rfunctions/med_SL_norm.R")
setwd('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_v2')


#########################
## Load data - STR axon Phospho MS2 SwissProt only
#########################

# Read PD peptide groups tabs into R
# Only retain high confidence peptides
# Each plex was run twice (survey and with exclusion list). Files were treated as fractions during PD analysis.

raw.phospeptide.pd <- read_excel("part2_phos_summarization/STR_Axon_1_and_2_Phos_two_injections_Peptide_Groups_Correct.xlsx") %>%
  dplyr::filter(Confidence == 'High')

# Read annotation file
annotation.pd_phos <- read.table("part2_phos_summarization/PD_STRaxon_annotation.txt", sep="\t", header=TRUE)

# Replace the abundance column names with bioreplicates names
sample_names <- annotation.pd_phos %>%
  dplyr::filter(Fraction == 1) %>%
  dplyr::mutate(columns = paste(paste0('Plex',Mixture), BioReplicate, sep = '.'))

sample_names$columns <- gsub('P50_', 'adult', sample_names$columns)
sample_names$columns <- gsub('P18_', 'preweanling', sample_names$columns)
sample_names$columns <- gsub('P11_', 'earlypostnatal', sample_names$columns)
sample_names$columns <- gsub('P5_', 'neonate', sample_names$columns)
sample_names$columns <- gsub('CTX_H2B_', 'CTXH2B', sample_names$columns)
sample_names$columns <- gsub('CTX_negative_', 'CTXnegative', sample_names$columns)
sample_names$columns <- gsub('STR_negative_', 'STRnegative', sample_names$columns)
unique(sample_names$columns)

colnames(raw.phospeptide.pd)[16:47] <- sample_names$columns


#########################
## Load data - STR axon enriched and unenrich proteins
#########################

# Only retain peptides associated with filtered proteins
# STRaxon_full_proteome <- read.table(file = 'part2_phos_summarization/20210811_STRaxon_full_proteome.txt', sep = '\t', header = TRUE) %>% data.frame()

STRaxon_full_proteome <- read_csv(file = 'part2_phos_summarization/20220130_STRaxon_full_proteome.csv') %>% data.frame()



#########################
## Perform peptide and batch normalization
#########################

# select useful columns
# remove Cre_negative columns
raw.phospeptide.pd_norm <- raw.phospeptide.pd[,c(1,3,4,11,12,13,16:47)]
raw.phospeptide.pd_norm <- raw.phospeptide.pd_norm[,!grepl('negative',colnames(raw.phospeptide.pd_norm))]

names(raw.phospeptide.pd_norm)


# reference channel normalization adapted from P. Wilmart (IRS normalization)
raw.phospeptide.pd_norm$norm_average <- 
  apply(raw.phospeptide.pd_norm[,c(19,31)], 1, function(x) exp(mean(log(x), na.rm = TRUE)))

# calculate norm factor for plex1 and plex2
# compute the scaling factor for each peptide within each TMT set
raw.phospeptide.pd_norm$factor_plex1 <- 
   raw.phospeptide.pd_norm$norm_average / raw.phospeptide.pd_norm$Plex1.Norm_1
raw.phospeptide.pd_norm$factor_plex2 <- 
  raw.phospeptide.pd_norm$norm_average / raw.phospeptide.pd_norm$Plex2.Norm_2

# multiply scaling factor to each protein within plex
raw.phospeptide.pd_norm[,c(7:19)] <- raw.phospeptide.pd_norm[,c(7:19)] * raw.phospeptide.pd_norm$factor_plex1
raw.phospeptide.pd_norm[,c(20:31)] <- raw.phospeptide.pd_norm[,c(20:31)] * raw.phospeptide.pd_norm$factor_plex2

# remove reference channel related columns
raw.phospeptide.pd_norm <- raw.phospeptide.pd_norm[,-c(19,31:34)]


# perform limma-based batch effect correction and global median normalization
# batch information
# Intensity must be in log2 before batch correction and med_norm()
batch_colors <- c(rep(1,12), rep(2,11))

#phospep_norm_temp <- removeBatchEffect(x = log2(raw.phospeptide.pd_norm[,c(7:29)]),
#                                       batch = batch_colors) %>% med_norm()

phospep_norm_temp <- removeBatchEffect(x = log2(raw.phospeptide.pd_norm[,c(7:29)]),
                                       batch = batch_colors) %>% med_subtraction()

# assign normalized corrected log2 intensity to the dataframe
raw.phospeptide.pd_norm[,c(7:29)] <- phospep_norm_temp

# replace non finite elements with NA
# replace NaN with NA
# remove all negative signals with NA
raw.phospeptide.pd_norm[,c(7:29)][sapply(raw.phospeptide.pd_norm[,c(7:29)], is.infinite)] <- NA
raw.phospeptide.pd_norm[,c(7:29)][sapply(raw.phospeptide.pd_norm[,c(7:29)], is.nan)] <- NA
raw.phospeptide.pd_norm[,c(7:29)][raw.phospeptide.pd_norm[,c(7:29)] <= 0] <- NA

# re-arrange columns by names
raw.phospeptide.pd_norm <- raw.phospeptide.pd_norm[,
      c("Annotated Sequence", "Modifications", "Master Protein Accessions",
        "Positions in Master Proteins", "Modifications in Master Proteins",
        "Plex1.neonate1", "Plex1.neonate2", "Plex2.neonate3", "Plex2.neonate4", "Plex2.neonate5",
        "Plex1.earlypostnatal1", "Plex1.earlypostnatal2", "Plex1.earlypostnatal3", "Plex2.earlypostnatal4", 
        "Plex2.earlypostnatal5",
        "Plex1.preweanling1", "Plex1.preweanling2", "Plex2.preweanling3", "Plex2.preweanling4",
        "Plex2.preweanling5",
        "Plex1.adult1", "Plex1.adult2", "Plex1.adult3", "Plex2.adult4", "Plex2.adult5")]

# condition colors
colors = c(rep('gray75',5),
           rep('gray55',5),
           rep('gray35',5),
           rep('gray15',5))
boxplot(raw.phospeptide.pd_norm[,c(6:25)], col = colors)

# batch information
batch_colors <- c(1, 1, 2, 2, 2,
                  1, 1, 1, 2, 2,
                  1, 1, 2, 2, 2, 
                  1, 1, 1, 2, 2,
                  1, 1, 2)
boxplot(raw.phospeptide.pd_norm[,c(6:25)], col = batch_colors)

colors2 <-  c(rep('#DD6E42',5), 
              rep('#E8DAB2',5),
              rep('#4F6D7A',5),
              rep('#C0D6DF',5),
              rep('#EAEAEA',3))
boxplot(raw.phospeptide.pd_norm[,c(6:25)], col = colors2, las = 2, cex.axis = 0.7)


#########################
## Filter phosphopeptides for only proteins that matter for the analysis of time points
#########################

raw.phospeptide.pd_norm <- cbind(raw.phospeptide.pd_norm, raw.phospeptide.pd[,c(5:10, 14:15)])


phosphopeptide_STRaxon_full <- raw.phospeptide.pd_norm %>%
  dplyr::filter(`Master Protein Accessions` %in% STRaxon_full_proteome$Protein)


# clustering using all phosphopeptides
# it seems like there is a plex effect
heatmap_phospho_STRaxon <- raw.phospeptide.pd_norm[,c(6:25)] %>% drop_na()
colnames(heatmap_phospho_STRaxon) 
Heatmap(matrix = as.matrix(heatmap_phospho_STRaxon), show_row_names = FALSE)

# 2709 proteins
# 1864 phosphoproteins
# 1124 phosphoproteins in proteins
length(unique(STRaxon_full_proteome$Protein))
length(unique(raw.phospeptide.pd_norm$`Master Protein Accessions`))
length(unique(phosphopeptide_STRaxon_full$`Master Protein Accessions`))


# writing ALL phosphopeptides
# write.table(phosphopeptide_STRaxon_full,
#            file = '/Users/vasdumrong/Box/VasD_projects/Project_STRaxon/analysis/Phospho_MS2/20210811_phosphopeptide_STRaxon_full.txt',
#            sep = '\t', row.names = FALSE)

#write.csv(phosphopeptide_STRaxon_full, file = 'part2_phos_summarization/20220131_phosphopeptide_STRaxon_full.csv')








