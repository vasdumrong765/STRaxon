library(tidyverse)
library(limma)
library(OmnipathR)
library(readxl)
setwd('E:/Vas/Data_analysis/PHONEMES_straxon/straxon_GSK3b_Fyn_Cdk5/')

#########################
## Prerequisite 1: Use OmpipathR to build kinase-substrate Network for the mouse proteome
#########################

# OmnipathR
# Obtaining the PTMS relations online from OmnipathR package and retaining only phosphorylation/dephosphorylation modifications
# For mouse organism = 10090
ptms_mmu <- get_signed_ptms(
  enzsub = import_omnipath_enzsub(organism = 10090),
  interactions = import_omnipath_interactions(organism = 10090))


ptms_mmu <- ptms_mmu[which(ptms_mmu$modification%in%c("phosphorylation", "dephosphorylation")), ]


# re-map to UNIPROT
# unique_KS <- data.frame(unique(c(ptms_mmu$enzyme, ptms_mmu$substrate)))
# write.csv(unique_KS, file = 'C:/Users/vasdu/Desktop/PHONEMES_STRaxon/unique_KS.csv' ,row.names = FALSE)
# rm(unique_KS)
unique_KS_maped <- read_excel('./rawd/unique_KS_mapped.xlsx')

# only include KS interactions that mapped to Swissprot
ptms_mmu_mapped <- merge(ptms_mmu, unique_KS_maped[,c(2,3)], 
                         by.x = 'enzyme', by.y = 'Entry')
ptms_mmu_mapped <- merge(ptms_mmu_mapped, unique_KS_maped[,c(2,3)], 
                         by.x = 'substrate', by.y = 'Entry')

# Initializing an empty PKN for PHONEMEeS
bn <- matrix(data = , nrow = nrow(ptms_mmu_mapped), ncol = 8)
colnames(bn) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

# Filling the initialized PKN with information from OmniPath
bn[, 1] <- ptms_mmu_mapped$substrate
bn[, 2] <- ptms_mmu_mapped$`Entry name.y` #substrate Entry name
bn[, 3] <- ptms_mmu_mapped$enzyme
bn[, 4] <- ptms_mmu_mapped$`Entry name.x` #enzyme Entry name
bn[, 5] <- ptms_mmu_mapped$residue_type
bn[, 6] <- ptms_mmu_mapped$residue_offset
bn[, 7] <- paste0("e", 1:nrow(bn))
bn[, 8] <- paste0(ptms_mmu_mapped$`Entry name.y`, "_", 
                  ptms_mmu_mapped$residue_type, ptms_mmu_mapped$residue_offset)

# Making the PKN as a data-frame object
allD <- as.data.frame(bn)
allD$S.AC <- as.character(allD$S.AC)
allD$S.ID <- as.character(allD$S.ID)
allD$K.AC <- as.character(allD$K.AC)
allD$K.ID <- as.character(allD$K.ID)
allD$res <- as.character(allD$res)
allD$pos <- as.character(allD$pos)
allD$SID <- paste0("e", 1:nrow(allD))
allD$S.cc <- as.character(allD$S.cc)


# Still need to convert genes to UniprotID
# First export
# write.csv(allD, file = 'allD_mouse.csv' ,row.names = FALSE)





#########################
## Prerequisite 2: Build GMMres object class from data
#########################

# In this step, we will reformat and analyze data to get two objects
# "GMM.ID", "GMM","GMM.wFC"

# 2.1
# Prepare data.IDmap
# Load data without FC calculation or q-values
straxon_data <- read_csv(file = './rawd/20211114_unpacked_modResidue.csv' , col_names = TRUE)

straxon_data$FirstmodResidue <- str_extract(straxon_data$modResidue, pattern = '([^;]+)')
straxon_data$First_site <- substr(straxon_data$FirstmodResidue, start = 1, stop = 1)
straxon_data$First_ResNum <- str_extract(straxon_data$FirstmodResidue, pattern = '[0-9]+')
straxon_data$countPhos <- str_count(straxon_data$modResidue, ';') + 1

targets_axon <- data.frame(c('neonate_1', 'neonate_2', 'neonate_3', 'neonate_4', 'neonate_5',
    'earlypostnatal_1', 'earlypostnatal_2', 'earlypostnatal_3', 
    'earlypostnatal_4', 'earlypostnatal_5',
    'preweanling_1', 'preweanling_2', 'preweanling_3', 
    'preweanling_4', 'preweanling_5',
    'adult_1', 'adult_2', 'adult_3', 'adult_4', 'adult_5'))
targets_axon$Condition <- str_extract(targets_axon[,1], pattern = '^[^_]+(?=_)')
colnames(targets_axon) <- c('sample', 'condition')

# get the uniprot mapping for straxon_data
straxon_uniprotmapping <- read_excel('./rawd/20211114_modResidue_uniprotmapping.xlsx')
straxon_data <- merge(straxon_data, straxon_uniprotmapping[, c(2,3,6)], 
                      by.x = 'Master Protein Accessions', by.y = 'Entry')

# reformat to create data.IDmap
straxon_data_reformat <- straxon_data[, c(1, 65:70, 7:25)]
straxon_data_reformat$phosID <- paste0(straxon_data_reformat$`Entry name`,'_',
                                       straxon_data_reformat$FirstmodResidue,'_',
                                       straxon_data_reformat$countPhos)
straxon_data_reformat$Prot_site <- paste0(straxon_data_reformat$`Entry name`,'_',
                                          straxon_data_reformat$FirstmodResidue)

data.IDmap <- straxon_data_reformat %>%
  dplyr::select(c(phosID, FirstmodResidue, Prot_site))
colnames(data.IDmap) <- c("dataID", "UPID", "S.cc")


# Sanity Check Are there any overlap between S.cc mmu database and experimental data
length(intersect(allD$S.cc, data.IDmap$S.cc))


# 2.2
# Preparing values for GMM object

# 2.2.1 Use LIMMA, lmFit takes log2-transformed input
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
data_limma <- straxon_data[,c(7:26)]
colnames(data_limma) <- targets_axon$sample

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

tmp1$PhosID <- straxon_data_reformat$phosID
tmp2$PhosID <- straxon_data_reformat$phosID
tmp3$PhosID <- straxon_data_reformat$phosID


# log2FC
matrix_data <- merge(tmp1[,c(1,7,5)], tmp2[,c(1,7,5)], by = 'PhosID')
matrix_data <- merge(matrix_data, tmp3[,c(1,7,5)], by = 'PhosID')
matrix_data <- matrix_data[,c(1,2,4,6,3,5,7)]
colnames(matrix_data) <- c('PhosID','earlypostnatal', 'preweanling', 'adult',
                           'adjPval_earlypostnatal', 'adjPval_preweanling', 'adjPval_adult')

##### placeholder
# for duplicate entries, pick the first one
matrix_data$PhosID <- as.factor(matrix_data$PhosID)
matrix_data <- matrix_data %>%
  dplyr::filter(!duplicated(PhosID))

rownames(matrix_data) <- matrix_data$PhosID


matrix_qval <- matrix_data[,c(5,6,7)] #adjpvalue
matrix_data <- matrix_data[,c(2,3,4)] #log2FC
rm(tmp1, tmp2, tmp3)

# 2.2.2 
threshQ <- 0.01 # defining threshold q-value

## details
##
# assigning corresponding q-value to each measurement
# setting missing q-values to 1
tp1.p <- as.numeric(as.matrix(matrix_qval)[, 1])
tp1.p[which(is.na(tp1.p))] <- 1     

# scoring vector tp1.lo, status vector tp1.c
# initializing scoring vector to be assigned to each measurement
# initializing status vector to be assigned to each measurement (Control - C/Perturbed P)
tp1.lo <- rep(0, nrow(matrix_data))                 
names(tp1.lo) <- rownames(matrix_data)
tp1.c <- rep("C", length(tp1.p)) 

# assigning the Perturbed - P status to measurements with q-value < threshQ
# assigning scores to each measurement (negative score for P nodes and positive scores for C nodes)
tp1.c[which(tp1.p < threshQ)] <- "P"
tp1.lo <- log(x = tp1.p/threshQ, base = 2)   
table(tp1.c)


# initializing regulation status tp1.s
tp1.s <- rep("OK", length(tp1.p))            
names(tp1.s) <- rownames(matrix_data)

# giving FP status to measurements with lower fc than 1.5 and OK status to measurement with higher fc than 1.5
# Not implemented in this analysis
# tp1.s[which(min > -0.58)] <- "FP"            
# tp1.s[which(max > 0.58)] <- "OK"
# tp1.s[which(abs(numData.Avg[,1]) < 1)] <- "FP"
for(i in 1:length(tp1.p)){
  if(tp1.c[i]=="C"){
    tp1.s[i] <- "OK"
  }
}
table(tp1.s)


## repeat above procedure to all other conditions
tp2.p <- as.numeric(as.matrix(matrix_qval)[, 2])
tp2.p[which(is.na(tp2.p))] <- 1     
tp2.lo <- rep(0, nrow(matrix_data))                 
names(tp2.lo) <- rownames(matrix_data)
tp2.c <- rep("C", length(tp2.p)) 
tp2.c[which(tp2.p < threshQ)] <- "P"
tp2.lo <- log(x = tp2.p/threshQ, base = 2)   
table(tp2.c)
tp2.s <- rep("OK", length(tp2.p))            
names(tp2.s) <- rownames(matrix_data)
for(i in 1:length(tp2.p)){
  if(tp2.c[i]=="C"){
    tp2.s[i] <- "OK"
  }
}
table(tp2.s)

##
tp3.p <- as.numeric(as.matrix(matrix_qval)[, 3])
tp3.p[which(is.na(tp3.p))] <- 1     
tp3.lo <- rep(0, nrow(matrix_data))                 
names(tp3.lo) <- rownames(matrix_data)
tp3.c <- rep("C", length(tp3.p)) 
tp3.c[which(tp3.p < threshQ)] <- "P"
tp3.lo <- log(x = tp3.p/threshQ, base = 2)   
table(tp3.c)
tp3.s <- rep("OK", length(tp3.p))            
names(tp3.s) <- rownames(matrix_data)
for(i in 1:length(tp3.p)){
  if(tp3.c[i]=="C"){
    tp3.s[i] <- "OK"
  }
}
table(tp3.s)




# 2.2.3 assembling GMM objects

###
# Preparing GMM objects
GMM<-vector("list", dim(matrix_data)[1])
names(GMM)<-rownames(matrix_data)
for(i in 1:length(GMM)){
  GMM[[i]]<-rbind(
    c(as.character(tp1.lo[i]), as.character(tp1.c[i]), as.character(tp1.p[i]), as.character(tp1.s[i])),
    c(as.character(tp2.lo[i]), as.character(tp2.c[i]), as.character(tp2.p[i]), as.character(tp2.s[i])),
    c(as.character(tp3.lo[i]), as.character(tp3.c[i]), as.character(tp3.p[i]), as.character(tp3.s[i]))
  )
  colnames(GMM[[i]])<-c("Indiv", "clus","FCvCaPval","status")
  rownames(GMM[[i]])<-c("tp_1earlypostnatal", "tp_2preweanling", "tp_3adult")
}

# with fold change... but in this case just log2 intensities
GMM.wFC<-vector("list", dim(matrix_data)[1])
names(GMM.wFC)<-rownames(matrix_data)
for(i in 1:length(GMM.wFC)){
  GMM.wFC[[i]]<-rbind(
    c(as.character(tp1.lo[i]), as.character(tp1.c[i]), as.character(tp1.p[i]), as.character(tp1.s[i]), as.character(matrix_data[i, 1])),
    c(as.character(tp2.lo[i]), as.character(tp2.c[i]), as.character(tp2.p[i]), as.character(tp2.s[i]), as.character(matrix_data[i, 2])),
    c(as.character(tp3.lo[i]), as.character(tp3.c[i]), as.character(tp3.p[i]), as.character(tp3.s[i]), as.character(matrix_data[i, 3]))
  )
  colnames(GMM.wFC[[i]])<-c("Indiv", "clus", "FCvCaPval", "status", "FCvC")
  rownames(GMM.wFC[[i]])<-c("tp_1earlypostnatal", "tp_2preweanling", "tp_3adult")
}

###
# Saving GMM object as a list

GMM.ID <- data.IDmap
colnames(GMM.ID) <- c("dataID", "UPID", "S.cc")
GMM.ID <- as.data.frame(GMM.ID)

save(list=c("GMM.ID", "GMM","GMM.wFC"), file="dataGMM_STRaxon.RData")
save(list=c("allD", "ptms_mmu_mapped", "bn"), file="KS_mmu_db.RData")

