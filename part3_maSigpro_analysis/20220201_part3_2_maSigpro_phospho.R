## Dumrongprechachan et al 2022
## Part 3 maSigpro phospho
## original code by A2IDEA, modified by V.D.

library(tidyverse)
library(maSigPro)
setwd('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_v2/part3_maSigpro_analysis')

## True means that redundant peptides will be collapsed
## False means redundant peptides will be removed
tryHard <- FALSE

getClusterMembership <- function(group, K=8) {
  ## group: one of the values from names(sigs$sig.genes)
  ## K: number of clusters (default is 8 but this is arbitrary)
  
  x <- see.genes(sigs$sig.genes[[group]], 
                 time.col=2,
                 repl.col=1,
                 group.cols=c(3:6),
                 show.fit=F, 
                 k=K,
                 dis=design$dis)
  
  clusterMembers <- data.frame(feature=names(x$cut), 
                               peptide=NA,
                               prot=NA,
                               cluster=x$cut, 
                               stringsAsFactors=F,
                               row.names=NULL)
  clusterMembers$peptide <- gsub("(.+)\\..+$", "\\1", clusterMembers$feature)
  clusterMembers$prot <- gsub(".+\\.(.+)$", "\\1", clusterMembers$feature)
  
  
  clusterMembers <- arrange(clusterMembers, cluster, prot)    
  
  ## This produces a matrix with protein IDs as rows and cluster labels as columns
  ## numbers in each cell tell you how many peptides that protein has in each cluster
  ## PROBLEM: how do you interpret a protein having different peptides in different clusters?
  peptide_usage <- table(clusterMembers$prot, clusterMembers$cluster)
  
  ret <- list(membership=clusterMembers, usage=peptide_usage)
  return(ret)
}

########################
# Load data
########################
# normalized phosphopeptide intensities (batch effect corrected + med_norm)
ms2 <- read_csv('rawd/20220131_phosphopeptide_STRaxon_full.csv')
phos_reformat <- read.delim('rawd/20220201_phosphopeptide_reformatted.txt', sep = '\t')

tmp <- merge(ms2, phos_reformat,
             by.x = c("Annotated Sequence", "Master Protein Accessions", "Positions in Master Proteins"),
             by.y = c("seq", "Master.Protein.Accessions", "Positions.in.Master.Proteins"))
ms2 <- tmp
rm(tmp)


# you need to rearrange the columns
g <- grep("modPep", colnames(ms2))
j <- grep("Master Protein", colnames(ms2)) ## start point
k <- grep("Qvality PEP", colnames(ms2)) - 1## stop point

tmp <- ms2[,c(g,j:k)]
ms2 <- select(tmp, -c(`Positions in Master Proteins`, `Modifications in Master Proteins`))
colnames(ms2)[c(1,2)] <- c("pep", "prots")
rm(tmp, g, j, k)



########################
# Unroll peptides
########################
# function to unroll the protein-peptide listings
unroll <- function(i) {
  all_prots <- strsplit(ms2$prots[i], split=";", fixed=T)[[1]]
  all_prots <- sapply(all_prots, function(x) { trimws(gsub("-\\d+$", "", x)) } )
  all_prots <- unique(all_prots)
  
  ## all isoforms of the same protein, collapse to single value
  if(length(all_prots) == 1) {
    ret <- ms2[i,]
    ret$prots <- all_prots
  } else {
    ## multiple distinct protein IDs, unroll them into new rows
    ret <- data.frame()
    for(cur_prot in all_prots) {
      tmp <- ms2[i,]
      tmp$prots <- cur_prot
      ret <- rbind(ret, tmp)
    }
  }
  
  return(ret)
}

# unroll the master protein accessions so that you have only 1 peptide + protein pairing
d <- data.frame()
for(i in 1:nrow(ms2)) {
  d <- rbind(d, unroll(i))
}
ms2 <- d
rm(d, i, unroll)


########################
# Remove duplicates
########################
# Each protein + peptide must only be represented once.
# For protein + peptide pairings that occur more than once, take the row instance with the largest `rowSum`.
# You can use the `tryHard = TRUE` argument to use this filtering approach.

#tmp <- data.frame(ms2, rowSum=apply(ms2[,-c(1,2)], 1, sum), stringsAsFactors=F)
tmp <- data.frame(ms2, rowSum=apply(ms2[,-c(1,2,3)], 1, sum), stringsAsFactors=F)
tmp$k <- paste0(tmp$pep,".",tmp$prots)
t1 <- select(tmp, k, rowSum)

t2 <- group_by(t1, k) %>% summarize(freq=n(), maxSum=max(rowSum))

if(tryHard) { 
  ## collapse peptides that are indistinguishable 
  ## We will select the row for the peptide that has the largest rowSum value
  t3 <- inner_join(x=t1, y=t2, by=c("k"="k", "rowSum"="maxSum"))
  
  d <- inner_join(x=tmp, y=t3, by=c("k", "rowSum"))
  
  ms2 <- select(d, -rowSum, -k, -freq)
  rm(t3)
  
} else {
  ## remove peptides that are indistinguishable
  keep_peps <- filter(t2, freq == 1) %>% as.data.frame()
  keep_peps <- keep_peps[,1]
  keep_peps <- unique(keep_peps)
  
  d <- filter(tmp, k %in% keep_peps)
  
  ms2 <- select(d, -rowSum, -k) 
}

rm(tmp, t1, t2, d)

# Keep complete cases
tmp <- as.matrix(ms2[,-c(1,2)])
j <- complete.cases(tmp)

ms2 <- ms2[j,]
rm(tmp, j)

ms2 <- ms2[,-3]


########################
# Study design
########################
pheno <- data.frame(lab=colnames(ms2)[-c(1,2)], 
                    plex=0, state=NA, rep=0,
                    stringsAsFactors=F, 
                    row.names=colnames(ms2)[-c(1,2)])

pheno$state <- gsub(".*\\.(.+)\\d+$", "\\1", pheno$lab)
pheno$plex <- as.integer(gsub("Plex(\\d+)\\..*", "\\1", pheno$lab))
pheno$rep <- as.integer(gsub("Plex\\d+\\.\\w+(\\d+)$", "\\1", pheno$lab))
pheno <- pheno[,-1]

pheno$state <- factor(pheno$state, levels=c("neonate", "earlypostnatal", "preweanling", "adult"))

## maSigPro code starts here\
## time points:
## "neonate"=1, "earlypostnatal"=2, "preweanling"=3, "adult"=4

d <- dplyr::select(pheno, rep, plex)
d$Time <- as.numeric(pheno$state)
d$Neonate=0
d$Earlypostnatal=0
d$Preweanling=0
d$Adult=0

## You need to annotate the comparison columns
d$Neonate[ grep("neonate", rownames(d)) ] <- 1
d$Earlypostnatal[ grep("early", rownames(d)) ] <- 1
d$Preweanling[ grep("pre", rownames(d)) ] <- 1
d$Adult[ grep("adult", rownames(d)) ] <- 1

d <- dplyr::select(d, -plex)

design <- make.design.matrix(edesign=d, 
                             degree=2,
                             time.col=2, ## column for time
                             repl.col=1, ## column for replicates
                             group.cols=c(3,4,5,6))

rm(d)


########################
# maSigPro regression model
########################
d <- as.matrix(select(ms2, -pep, -prots))
rownames(d) <- paste0(ms2$pep,".",ms2$prots)


## Perform regression fit for each protein
fit <- p.vector(data=d, 
                design=design, 
                min.obs=1,
                #Q = 0.05,
                Q = 1, 
                MT.adjust="BH",
                counts=F, 
                family=gaussian())


tstep <- T.fit(fit, 
               step.method = "backward", 
               min.obs=1, 
               alfa = 0.05,
               family=gassian())


## Get out all peptides+proteins that are significant
#sigs <- get.siggenes(tstep, rsq = 0.6, vars = "all")
sigs <- get.siggenes(tstep, rsq = 0, vars = "all")

# 8 clusters
sigs_all <- see.genes(sigs$sig.genes, show.fit = T, dis =design$dis,
                      cluster.method="hclust" ,cluster.data = 1, k = 8)

map <- select(ms2, pep, prots) %>%
  mutate(feature = paste(pep, prots, sep = ".")) %>%
  rename(prot = prots)

tmp <- as.data.frame(sigs_all$cut) %>%
  tibble::rownames_to_column("feature")  %>%
  dplyr::rename(cluster = "sigs_all$cut") %>%
  left_join(map, by = "feature")


########################
# Export results
########################
# extract coefficient, p-value, R2 statistics
tmp2 <- merge(sigs[['sig.genes']]$coefficients[,2:4], sigs[['sig.genes']]$sig.pvalues[,c(4:6, 1:2)],
              by = 0, all = TRUE)
colnames(tmp2)[1] <- 'pep.Prot'

tmp3 <- merge(tmp, tmp2, by.x = 'feature', by.y = 'pep.Prot')

res <- cbind(padjust=fit$p.adjusted[ fit$p.adjusted < 1 ], fit$SELEC)
res <- data.frame(feature=rownames(res),
                  pep=gsub("(.+)\\.\\w+", "\\1", rownames(res)),
                  prot=gsub(".+\\.(\\w+)", "\\1", rownames(res)),
                  res,
                  row.names=NULL,
                  stringsAsFactors=F)
res1 <- right_join(x=tmp, y=res,by= c("feature", "pep", "prot"))

tmp10 <- phos_reformat[,-6]
colnames(tmp10)[c(3,6)] <- c('prot', 'pep')
res2 <- left_join(res1, tmp10, by = c('pep', 'prot'))

res3 <- res2
res3$FirstmodResidue <- str_extract(res3$Modifications.in.Master.Proteins, pattern = '([^;]+)')
res3$FirstmodResidue <- str_extract(res3$FirstmodResidue, pattern = '\\[(.*?)\\]')
res3$FirstmodResidue <- str_replace(res3$FirstmodResidue, pattern = '\\[', replacement = '')
res3$FirstmodResidue <- str_replace(res3$FirstmodResidue, pattern = '\\]', replacement = '')
res3$First_site <- substr(res3$FirstmodResidue, start = 1, stop = 1)
res3$First_ResNum <- str_extract(res3$FirstmodResidue, pattern = '[0-9]+')
#res3$countPhos <- str_count(res3$modResidue, ';') + 1

#x <- inner_join(x=cluster_members, y=res, by="feature")
#x <- arrange(x, prot, pep, group, cluster)

#write.table(x, file=outF, sep="\t", row.names=F)
#write.table(res2, file= "out/20220201_phosphopeptide_STRaxon_full_ALL_8clusters_maSigPro_v2.tsv", sep="\t", row.names=F)

rm(d)


########################
# note - hinge analysis
########################

set.seed(1234)
k.max <- 25

wss <- sapply(1:k.max, function(k) { kmeans(d, k, nstart=50, iter.max=15)$tot.withinss })

plot(1:k.max, y=wss, type="b", pch=19, xaxt="n", xlab="# Clusters K", ylab="Total within-cluster sum of squares")
axis(1, at=seq(1,k.max), las=1)

rm(k.max, wss)
