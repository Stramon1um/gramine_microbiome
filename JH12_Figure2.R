#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Figure 2 used in Maver et al., manuscript
#  Revision 11/20 
#  mauro.maver@unibz.it
#  d.bulgarelli@dundee.ac.uk  
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#required packages 
library ("phyloseq")
library("DESeq2")
library ("ggplot2")
source("tern_e.R")
library("grid")

#R session info
sessionInfo()

#import the Deseq file
#G0
JH12_G0_cds <- readRDS("JH12_G0_genus_cds_2.rds")
JH12_G0_cds

#execute the differential count analysis with the function DESeq 
JH12_cds_test <- DESeq(JH12_G0_cds, fitType="local", betaPrior = FALSE) 

#define the OTUs differentially enriched in the rhizosphere samples
Morex_rhizo_G0 <- results(JH12_cds_test, contrast = c("Description", "QUARRY", "MOREXQUARRY")) 
Barke_rhizo_G0 <- results(JH12_cds_test, contrast = c("Description", "QUARRY", "BARKEQUARRY")) 

#inspect the result files
Morex_rhizo_G0
Barke_rhizo_G0

#extract  OTUs whose adjusted p.value in a given comparison is below 0.05 
Morex_rhizo_FDR005 <- Morex_rhizo_G0[(rownames(Morex_rhizo_G0)[which(Morex_rhizo_G0$padj <0.05)]), ]
Barke_rhizo_FDR005 <- Barke_rhizo_G0[(rownames(Barke_rhizo_G0)[which(Barke_rhizo_G0$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, rhizosphere is at the second term of comparison)
Morex_rhizo_enriched <-  Morex_rhizo_G0[(rownames(Morex_rhizo_G0)[which(Morex_rhizo_G0$log2FoldChange < 0)]), ]
Barke_rhizo_enriched <-  Barke_rhizo_G0[(rownames(Barke_rhizo_G0)[which(Barke_rhizo_G0$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Morex_rhizo_enriched_FDR005_G0 <- intersect(rownames(Morex_rhizo_FDR005), rownames(Morex_rhizo_enriched))
Barke_rhizo_enriched_FDR005_G0 <- intersect(rownames(Barke_rhizo_FDR005), rownames(Barke_rhizo_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Morex_rhizo_enriched_FDR005_G0)
length(Barke_rhizo_enriched_FDR005_G0)

#Identify the genotype effect: using filtered data
G0_diagnostic <- counts(JH12_G0_cds)[unique(union(Morex_rhizo_enriched_FDR005_G0,Barke_rhizo_enriched_FDR005_G0)), ]
JH12_G0_cds_diagnostic <- JH12_G0_cds[rownames(G0_diagnostic), ]
JH12_G0_cds_diagnostic

#execute the differential count analysis with the function DESeq 
JH12_cds_diagnostic_test <- DESeq(JH12_G0_cds_diagnostic, fitType="local", betaPrior = FALSE) 

#define the OTUs differentially enriched in the rhizosphere samples
Morex_Barke_G0 <- results(JH12_cds_diagnostic_test, contrast = c("Description", "MOREXQUARRY", "BARKEQUARRY")) 

#inspect the result files
Morex_Barke_G0

#extract  OTUs whose adjusted p.value in a given comparison is below 0.05 
Morex_Barke_FDR005 <- Morex_Barke_G0[(rownames(Morex_Barke_G0)[which(Morex_Barke_G0$padj <0.05)]), ]

#enriched in rhizosphere of either genotypes
Morex_enriched <-  Morex_Barke_G0[(rownames(Morex_Barke_G0)[which(Morex_Barke_G0$log2FoldChange > 0)]), ]
Barke_enriched <-  Morex_Barke_G0[(rownames(Morex_Barke_G0)[which(Morex_Barke_G0$log2FoldChange < 0)]), ]

#differentially enriched between genotypes
Morex_enriched_G0_005 <- intersect(rownames(Morex_enriched), rownames(Morex_Barke_FDR005))
Barke_enriched_G0_005 <- intersect(rownames(Barke_enriched), rownames(Morex_Barke_FDR005))

#cumulative effect
sum(length(Morex_enriched_G0_005), length(Barke_enriched_G0_005))

#cumulative abundances
G0_effect <- sum(Morex_Barke_G0[union(Morex_enriched_G0_005,Barke_enriched_G0_005), ]$baseMean)
#proportion of reads
G0_proportion <- G0_effect/(sum(Morex_Barke_G0$baseMean)) * 100
G0_proportion

#Ternary plot
#extract the means of the genotypes for plotting
G0_base_mean <- sapply(levels(JH12_cds_diagnostic_test$Description), function(lvl) rowMeans(counts(JH12_cds_diagnostic_test,normalized=TRUE)[,JH12_cds_diagnostic_test$Description == lvl] ) )
#Build a ternary matrix for ternary plots
mean_Bulk <- as.data.frame(G0_base_mean[, 3])
mean_Morex<- as.data.frame(G0_base_mean[, 2])
mean_Barke <- as.data.frame(G0_base_mean[, 1])
temat_1A <- cbind(mean_Morex, mean_Barke, mean_Bulk)
colnames(temat_1A) <- c("Morex","Barke","Bulk")
dim(temat_1A)
#colors and plotting
dev.off()
fig_colors <- ifelse(rownames(temat_1A) %in% Morex_enriched_G0_005, "magenta","darkgrey")
names(fig_colors) <- rownames(temat_1A)
fig_colors[Barke_enriched_G0_005] <- "blue"

#FIGURE 2A
tern_e(temat_1A, scale = 1, prop=T, col=fig_colors, grid_color="black", labels_color="black", pch=19, main="G0 effect")

#save this plot for figure generation
#G24
JH12_G24_cds <- readRDS("JH12_G24_genus_cds_2.rds")
JH12_G24_cds

#execute the differential count analysis with the function DESeq 
JH12_cds_test <- DESeq(JH12_G24_cds, fitType="local", betaPrior = FALSE) 

#define the OTUs differentially enriched in the rhizosphere samples
Morex_rhizo_G24 <- results(JH12_cds_test, contrast = c("Description", "QUARRY", "MOREXQUARRY")) 
Barke_rhizo_G24 <- results(JH12_cds_test, contrast = c("Description", "QUARRY", "BARKEQUARRY")) 

#inspect the result files
Morex_rhizo_G24
Barke_rhizo_G24

#extract  OTUs whose adjusted p.value in a given comparison is below 0.05 
Morex_rhizo_FDR005 <- Morex_rhizo_G24[(rownames(Morex_rhizo_G24)[which(Morex_rhizo_G24$padj <0.05)]), ]
Barke_rhizo_FDR005 <- Barke_rhizo_G24[(rownames(Barke_rhizo_G24)[which(Barke_rhizo_G24$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, rhizosphere is at the second term of comparison)
Morex_rhizo_enriched <-  Morex_rhizo_G24[(rownames(Morex_rhizo_G24)[which(Morex_rhizo_G24$log2FoldChange < 0)]), ]
Barke_rhizo_enriched <-  Barke_rhizo_G24[(rownames(Barke_rhizo_G24)[which(Barke_rhizo_G24$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Morex_rhizo_enriched_FDR005_G24 <- intersect(rownames(Morex_rhizo_FDR005), rownames(Morex_rhizo_enriched))
Barke_rhizo_enriched_FDR005_G24 <- intersect(rownames(Barke_rhizo_FDR005), rownames(Barke_rhizo_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Morex_rhizo_enriched_FDR005_G24)
length(Barke_rhizo_enriched_FDR005_G24)

#Identify the genotype effect: using filtered data
G24_diagnostic <- counts(JH12_G24_cds)[unique(union(Morex_rhizo_enriched_FDR005_G24,Barke_rhizo_enriched_FDR005_G24)), ]
JH12_G24_cds_diagnostic <- JH12_G24_cds[rownames(G24_diagnostic), ]
JH12_G24_cds_diagnostic

#execute the differential count analysis with the function DESeq 
JH12_cds_diagnostic_test <- DESeq(JH12_G24_cds_diagnostic, fitType="local", betaPrior = FALSE) 

#define the OTUs differentially enriched in the rhizosphere samples
Morex_Barke_G24 <- results(JH12_cds_diagnostic_test, contrast = c("Description", "MOREXQUARRY", "BARKEQUARRY")) 

#inspect the result files
Morex_Barke_G24

#extract  OTUs whose adjusted p.value in a given comparison is below 0.05 
Morex_Barke_FDR005 <- Morex_Barke_G24[(rownames(Morex_Barke_G24)[which(Morex_Barke_G24$padj <0.05)]), ]

#enriched in rhizosphere of either genotypes
Morex_enriched <-  Morex_Barke_G24[(rownames(Morex_Barke_G24)[which(Morex_Barke_G24$log2FoldChange > 0)]), ]
Barke_enriched <-  Morex_Barke_G24[(rownames(Morex_Barke_G24)[which(Morex_Barke_G24$log2FoldChange < 0)]), ]

#differentially enriched between genotypes
Morex_enriched_G24_005 <- intersect(rownames(Morex_enriched), rownames(Morex_Barke_FDR005))
Barke_enriched_G24_005 <- intersect(rownames(Barke_enriched), rownames(Morex_Barke_FDR005))

#cumulative effect
sum(length(Morex_enriched_G24_005), length(Barke_enriched_G24_005))

#Ternary plot
#extract the means of the genotypes for plotting
G24_base_mean <- sapply(levels(JH12_cds_diagnostic_test$Description), function(lvl) rowMeans(counts(JH12_cds_diagnostic_test,normalized=TRUE)[,JH12_cds_diagnostic_test$Description == lvl] ) )
#Build a ternary matrix for ternary plots
mean_Bulk <- as.data.frame(G24_base_mean[, 3])
mean_Morex<- as.data.frame(G24_base_mean[, 2])
mean_Barke <- as.data.frame(G24_base_mean[, 1])
temat_1A <- cbind(mean_Morex, mean_Barke, mean_Bulk)
colnames(temat_1A) <- c("Morex","Barke","Bulk")
dim(temat_1A)
#colors and plotting
dev.off()
fig_colors <- ifelse(rownames(temat_1A) %in% Morex_enriched_G24_005, "magenta","darkgrey")
names(fig_colors) <- rownames(temat_1A)
fig_colors[Barke_enriched_G24_005] <- "blue"

#FIGURE 2B
tern_e(temat_1A, scale = 1, prop=T, col=fig_colors, grid_color="black", labels_color="black", pch=19, main="G24 effect")
#save this plot for figure generation

#G46
JH12_G46_cds <- readRDS("JH12_G46_genus_cds_2.rds")
JH12_G46_cds

#execute the differential count analysis with the function DESeq 
JH12_cds_test <- DESeq(JH12_G46_cds, fitType="local", betaPrior = FALSE) 

#define the OTUs differentially enriched in the rhizosphere samples
Morex_rhizo_G46 <- results(JH12_cds_test, contrast = c("Description", "QUARRY", "MOREXQUARRY")) 
Barke_rhizo_G46 <- results(JH12_cds_test, contrast = c("Description", "QUARRY", "BARKEQUARRY")) 

#inspect the result files
Morex_rhizo_G46
Barke_rhizo_G46

#extract  OTUs whose adjusted p.value in a given comparison is below 0.05 
Morex_rhizo_FDR005 <- Morex_rhizo_G46[(rownames(Morex_rhizo_G46)[which(Morex_rhizo_G46$padj <0.05)]), ]
Barke_rhizo_FDR005 <- Barke_rhizo_G46[(rownames(Barke_rhizo_G46)[which(Barke_rhizo_G46$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, rhizosphere is at the second term of comparison)
Morex_rhizo_enriched <-  Morex_rhizo_G46[(rownames(Morex_rhizo_G46)[which(Morex_rhizo_G46$log2FoldChange < 0)]), ]
Barke_rhizo_enriched <-  Barke_rhizo_G46[(rownames(Barke_rhizo_G46)[which(Barke_rhizo_G46$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Morex_rhizo_enriched_FDR005_G46 <- intersect(rownames(Morex_rhizo_FDR005), rownames(Morex_rhizo_enriched))
Barke_rhizo_enriched_FDR005_G46 <- intersect(rownames(Barke_rhizo_FDR005), rownames(Barke_rhizo_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Morex_rhizo_enriched_FDR005_G46)
length(Barke_rhizo_enriched_FDR005_G46)

#Identify the genotype effect: using filtered data
G46_diagnostic <- counts(JH12_G46_cds)[unique(union(Morex_rhizo_enriched_FDR005_G46,Barke_rhizo_enriched_FDR005_G46)), ]
JH12_G46_cds_diagnostic <- JH12_G46_cds[rownames(G46_diagnostic), ]
JH12_G46_cds_diagnostic

#execute the differential count analysis with the function DESeq 
JH12_cds_diagnostic_test <- DESeq(JH12_G46_cds_diagnostic, fitType="local", betaPrior = FALSE) 

#define the OTUs differentially enriched in the rhizosphere samples
Morex_Barke_G46 <- results(JH12_cds_diagnostic_test, contrast = c("Description", "MOREXQUARRY", "BARKEQUARRY")) 

#inspect the result files
Morex_Barke_G46

#extract  OTUs whose adjusted p.value in a given comparison is below 0.05 
Morex_Barke_FDR005 <- Morex_Barke_G46[(rownames(Morex_Barke_G46)[which(Morex_Barke_G46$padj <0.05)]), ]

#enriched in rhizosphere of either genotypes
Morex_enriched <-  Morex_Barke_G46[(rownames(Morex_Barke_G46)[which(Morex_Barke_G46$log2FoldChange > 0)]), ]
Barke_enriched <-  Morex_Barke_G46[(rownames(Morex_Barke_G46)[which(Morex_Barke_G46$log2FoldChange < 0)]), ]

#differentially enriched between genotypes
Morex_enriched_G46_005 <- intersect(rownames(Morex_enriched), rownames(Morex_Barke_FDR005))
Barke_enriched_G46_005 <- intersect(rownames(Barke_enriched), rownames(Morex_Barke_FDR005))

#cumulative effect
sum(length(Morex_enriched_G46_005), length(Barke_enriched_G46_005))

#color codiding
G46_genera <- union(Morex_enriched_G46_005, Barke_enriched_G46_005)


#Ternary plot
#extract the means of the genotypes for plotting
G46_base_mean <- sapply(levels(JH12_cds_diagnostic_test$Description), function(lvl) rowMeans(counts(JH12_cds_diagnostic_test,normalized=TRUE)[,JH12_cds_diagnostic_test$Description == lvl] ) )
#Build a ternary matrix for ternary plots
mean_Bulk <- as.data.frame(G46_base_mean[, 3])
mean_Morex<- as.data.frame(G46_base_mean[, 2])
mean_Barke <- as.data.frame(G46_base_mean[, 1])
temat_1A <- cbind(mean_Morex, mean_Barke, mean_Bulk)
colnames(temat_1A) <- c("Morex","Barke","Bulk")
dim(temat_1A)
#colors and plotting
dev.off()
fig_colors <- ifelse(rownames(temat_1A) %in% Morex_enriched_G46_005, "magenta","darkgrey")
names(fig_colors) <- rownames(temat_1A)
fig_colors[Barke_enriched_G46_005] <- "blue"

#FIGURE 2C
tern_e(temat_1A, scale = 1, prop=T, col=fig_colors, grid_color="black", labels_color="black", pch=19, main="G46 effect")
#save this plot for figure generation

#import taxonomic annotations for pie chart
dat_tax <- read.delim("JH12_dat_tax_noPlants_ordered.txt", row.names = 1)
#G0 - FIGURE 2D
dat_tax_G0 <- dat_tax[union(Morex_enriched_G0_005, Barke_enriched_G0_005), ]
ggplot(dat_tax_G0, aes(x=factor(1), fill=Class))+
  geom_bar(width = 1)+
  coord_polar("y")

#G24 - FIGURE 2E
dat_tax_G24 <- dat_tax[union(Morex_enriched_G24_005, Barke_enriched_G24_005), ]
ggplot(dat_tax_G24, aes(x=factor(1), fill=Class))+
  geom_bar(width = 1)+
  coord_polar("y")

#G46 - FIGURE 2F
dat_tax_G46 <- dat_tax[union(Morex_enriched_G46_005, Barke_enriched_G46_005), ]
ggplot(dat_tax_G46, aes(x=factor(1), fill=Class))+
  geom_bar(width = 1)+
  coord_polar("y")

#end
