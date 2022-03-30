# Time-Point Analysis for the Modeling the STR axon dataset
# version 2 - using kinases that are significantly changing with time?

# Load packages
library(PHONEMeS)
library(BioNet)
library(igraph)
library(hash)
library(dplyr)
library(readxl)
library(readr)
library(XML)
# cplex 
path_to_executable_solver <- 'C:/Program Files/IBM/ILOG/CPLEX_Studio201/cplex/bin/x64_win64/cplex.exe'
setwd('E:/Vas/Data_analysis/20220203_PHONEMES_update/')

# Loading database and data-object
load(file = "dataGMM_STRaxon.RData")
load(file = "KS_mmu_db.RData")
rm(bn)


# Preparing background network as a PHONEMeS input
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataInput<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)

# Choose the conditions for each time-point
conditions <- list(c("tp_1earlypostnatal"), c("tp_2preweanling"), c("tp_3adult"))
names(conditions) <- c("tp_1earlypostnatal", "tp_2preweanling", "tp_3adult")

# Choose the targets for each time-point
# Identifying top most regulated kinases overlapped across 3 time points using LIMMA with abs(log2FC) > log2(1.5)
# 29 kinases were differentially regulated
# overlap_pkinfam_sig <- read.table('overlap_pkinfam_sig.txt')
# overlap_pkinfam_sig <- c('O35099', 'P11798', 'P28028', 'P54754', 'P63318', 'Q6PHZ2', 'Q8VDU5', 'Q9JM52')

# Picking Fyn, Gsk3b, and Cdk5
overlap_pkinfam_sig <- c('P39688','Q9WV60','P49615')

# how many allD interactions for overlap_pkinfam_sig?
allD_for_pkinfam_sig <- allD %>%
  #dplyr::filter(K.AC %in% overlap_pkinfam_sig$x)
  dplyr::filter(K.AC %in% overlap_pkinfam_sig)

# check substrate-kinase interactions detected in the dataset
overlap_substrate <- intersect(allD_for_pkinfam_sig$S.cc, GMM.ID$S.cc)
overlap_allD <- allD %>%
  dplyr::filter(S.cc %in% overlap_substrate)
overlap_allD_count <- overlap_allD %>% group_by(K.ID) %>% tally() %>% arrange(desc(n))


#targetKinases <- overlap_allD_count[1:10, 'K.ID']$K.ID
#targetKinases <- overlap_allD_count$K.ID
#targets.P <- list(tp_1earlypostnatal=targetKinases, 
#                  tp_2preweanling=targetKinases, 
#                  tp_3adult=targetKinases)

targets.P <- list(tp_1earlypostnatal=c("FYN_MOUSE", "GSK3B_MOUSE", "CDK5_MOUSE"), 
                  tp_2preweanling=c("FYN_MOUSE", "GSK3B_MOUSE", "CDK5_MOUSE"), 
                  tp_3adult=c("FYN_MOUSE", "GSK3B_MOUSE", "CDK5_MOUSE"))

targets.P <- list(tp_1earlypostnatal=c("GSK3B_MOUSE"), 
                  tp_2preweanling=c("GSK3B_MOUSE"), 
                  tp_3adult=c("GSK3B_MOUSE"))

targets.P <- list(tp_1earlypostnatal=c("CDK5_MOUSE"), 
                  tp_2preweanling=c("CDK5_MOUSE"), 
                  tp_3adult=c("CDK5_MOUSE"))

targets.P <- list(tp_1earlypostnatal=c("FYN_MOUSE"), 
                  tp_2preweanling=c("FYN_MOUSE"), 
                  tp_3adult=c("FYN_MOUSE"))

# Next we assign the experimental conditions from which we get the measurements at each specific time-point
# In this case we have the same one experimental condition for each time-point
experiments <- list(tp1=c(1), tp2=c(2), tp3=c(3))



# Running multiple time-point variant of PHONEMeS and retain only those interactions which have a weight higher than 10/appear at least 10 times in the network
# separate solutions we have obtained out of the 100 runs we have set to perform (nIter=100).
# separate solutions we have obtained out of the 50 runs we have set to perform (nIter=50).
set.seed(383789)
resultsMulti = runPHONEMeS_mult(targets.P = targets.P, conditions = conditions, inputObj = dataInput, experiments = experiments, bg = bg, nIter = 50, solverPath = path_to_executable_solver)


nodeAttribudes <- assignAttributes(sif = resultsMulti[, c(1, 2, 4)], inputObj = dataInput, targets = targets.P, writeAttr = FALSE)

write.table(x = resultsMulti, file = "straxon_network_FYN.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = resultsMulti[which(as.numeric(resultsMulti[, 2])>=10), ], file = "straxon_network_FYN_obs10.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = nodeAttribudes, file = "straxon_attributes_FYN.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
