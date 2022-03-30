## Dumrongprechachan et al 2022
## Part 3 maSigpro proteins
## original code by A2IDEA, modified by V.D.

suppressMessages(library(dplyr))
options(dplyr.summarise.inform = FALSE) ## stop summarize warning in dplyr
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(limma)) ## for batch correction
suppressMessages(library(maSigPro)) ## package we are using for differential expression
suppressMessages(library(ggpubr))
rm(list=ls());  #source("~/.Rprofile");

setwd('/Volumes/GoogleDrive/My Drive/Vas_Klab_GoogleDrive/Project_STRaxon/analysis_v2/part3_maSigpro_analysis')
#setwd("~/Documents/A2IDEAS/UPMC/SOW_2021-7-03-004/PhosphoTimeCourse") # location dependent


#rawd <- read.delim("rawd/20210811_STRaxon_full_proteome.txt", as.is=T)
rawd <- read_csv(file = 'rawd/20220130_STRaxon_full_proteome.csv')

## final output file name
# outF <-  "out/20210811_STRaxon_full_proteome.timeCourse.tsv"
# out_clusters = gsub(".tsv", ".fullProteome_clusters.pdf", outF)

## you need to rearrange the columns
ms2 <- dplyr::select(rawd, Protein, starts_with("plex"))
colnames(ms2)[1] <- "prots"

## keep complete cases
tmp <- as.matrix(ms2[,-1])
j <- complete.cases(tmp)

ms2 <- ms2[j,]
rm(tmp, j)


## check data with PCA
pheno <- data.frame(lab=colnames(ms2)[-1], 
                    plex=0, state=NA, rep=0,
                    stringsAsFactors=F, 
                    row.names=colnames(ms2)[-1])

pheno$state <- gsub("plex\\d+\\.(.+)\\d+$", "\\1", pheno$lab)
pheno$plex <- as.integer(gsub("plex(\\d+)\\..*", "\\1", pheno$lab))
pheno$rep <- as.integer(gsub("plex\\d+\\.\\w+(\\d+)$", "\\1", pheno$lab))
pheno <- pheno[,-1]

pheno$state <- factor(pheno$state, levels=c("neonate", "earlypostnatal", "preweanling", "adult"))

mat <- as.matrix(ms2[,-1])
rownames(mat) <- ms2[,1]$prots


df_pca <- prcomp( t(mat), scale=T)

m <- data.frame(PC1=df_pca$x[,1], PC2=df_pca$x[,2], stringsAsFactors=F)

m <- cbind(m, pheno) 
m$plex <- as.factor(m$plex)
m$rep <- as.factor(m$rep)

x <- summary(df_pca)$importance
pc1 <- round((x[2,1] * 100), 2)
PC1_pct <- paste0("PC1 (", pc1, "%)")
pc2 <- round((x[2,2] * 100), 2)
PC2_pct <- paste0("PC2 (", pc2, "%)")
rm(x, pc1, pc2, df_pca)


## constant plotting variables to be used for all plots
theme_obj <-  theme(plot.title = element_text(size = 24, face = "bold"), 
                    axis.text = element_text(size = 8), 
                    legend.position = "right", 
                    legend.title = element_text(face = "bold", size = 10), 
                    legend.text = element_text(size = 8))


#======================== make 3 PC plots =================================
p1 <- ggplot(m) +
  geom_point(aes(x=PC1, y=PC2, color=plex), alpha = 0.75, size = 3) + 
  ggtitle("Plex") + 
  xlab(PC1_pct) + ylab(PC2_pct) +
  theme_bw() + theme_obj


p2 <- ggplot(m) +
  geom_point(aes(x=PC1, y=PC2, color=state), alpha=0.75, size=3) + 
  ggtitle("State") + 
  xlab(PC1_pct) + ylab(PC2_pct) +
  theme_bw() + theme_obj

p3 <- ggplot(m) +
  geom_point(aes(x=PC1, y=PC2, color=rep), alpha=0.75, size=3) + 
  ggtitle("Replicate") + 
  xlab(PC1_pct) + ylab(PC2_pct) +
  theme_bw() + theme_obj

p <- cowplot::plot_grid(p1, p2, p3, nrow=1, ncol=3)
ggsave(plot = p, filename = "figures/full_protein_timeCourse_PCAall.pdf", units = "in", height = 7, width = 11)

rm(p1, p2, p3, mat, m, PC1_pct, PC2_pct, theme_obj, p)


############################
## maSigPro code starts here
############################

## time points:
## "neonate"=1, "earlypostnatal"=2, "preweanling"=3, "adult"=4

d <- select(pheno, rep, plex)
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

d <- select(d, -plex)

design <- make.design.matrix(edesign=d, 
                             degree=1,
                             time.col=2, ## column for time
                             repl.col=1, ## column for replicates
                             group.cols=c(3,4,5,6))

rm(d)

d <- as.matrix(ms2[,-1])
rownames(d) <-ms2[,1]$prots


## Perform regression fit for each protein
fit <- p.vector(data=d, 
                design=design, 
                min.obs=1, 
                counts=F, 
                Q = 0.99, # to print out different genes with significance 
                family=gaussian())


tstep <- T.fit(fit, 
               step.method = "backward", 
               min.obs=1, 
               alfa = 0.05, 
               family=gaussian())

## Get out all proteins that are significant
## 
## 

map <- dplyr::select(rawd, Protein, `Gene Symbol`)

sigs <- get.siggenes(tstep, rsq = 0, vars = "all")

#suma2Venn(sigs$summary)

sigs_all <- see.genes(sigs$sig.genes, show.fit = T, dis =design$dis,
                      cluster.method="hclust" ,cluster.data = 1, k = 8)

tmp <- as.data.frame(sigs_all$cut) %>%
  tibble::rownames_to_column("Protein")  %>%
  dplyr::rename(cluster = "sigs_all$cut") %>%
  left_join(map, by = "Protein")


#write.table(tmp, file= "out/20220201_STRaxon_full_proteome.timeCourse_8clusterMembership_maSigPro_ALL.tsv", sep="\t", row.names=F)


#see.genes(sigs$sig.genes, show.fit = T, dis =design$dis,
#          cluster.method="hclust" ,cluster.data = 1, k = 8)



## record clusters
res <- cbind(padjust=fit$p.adjusted[ fit$p.adjusted < 1 ], fit$SELEC)
res <- data.frame(prot=rownames(res), 
                  res,
                  row.names=NULL,
                  stringsAsFactors=F)


map <- select(rawd, Protein, Gene.Symbol)
res <- inner_join(x=map, y=res, by=c("Protein"="prot"))
res1 <- right_join(x=tmp, y=res,by=c("Protein", "Gene Symbol"))
#res <- arrange(res, group, cluster)


#write.table(res, file=outF, sep="\t", row.names=F)
write.table(res1, file= "out/20220201_STRaxon_full_proteome.timeCourse_ALL.tsv", sep="\t", row.names=F)

rm(d)

