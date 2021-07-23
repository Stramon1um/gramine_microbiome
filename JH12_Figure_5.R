#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Figure 5 used in Maver et al., manuscript
#  Revision 07/21 
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

colorblind_Palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#R session info
sessionInfo()

JH12_data_phyloseq_genus_rare <- readRDS(file = "JH12_data_phyloseq_genus_0221_Silva138.rds")


#import the Deseq file
#G0
JH12_G0_cds <- readRDS("JH12_G0_genus_cds_0221.rds")
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


Morex_prova <- as.data.frame(Morex_enriched_G0_005)
Barke_prova <- as.data.frame(Barke_enriched_G0_005)

Morex_prova$phylum <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_prova$Morex_enriched_G0_005,"Phylum"])
Barke_prova$phylum <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_prova$Barke_enriched_G0_005,"Phylum"])

morex_prova_proteo <- subset(Morex_prova, phylum == "Proteobacteria")
nrow(morex_prova_proteo)
morex_prova_actino <- subset(Morex_prova, phylum == "Actinobacteriota")
nrow(morex_prova_actino)
morex_prova_bacte <- subset(Morex_prova, phylum == "Bacteroidota")
nrow(morex_prova_bacte)
morex_prova_others <- subset(Morex_prova, phylum != "Proteobacteria" & phylum != "Actinobacteriota" & phylum != "Bacteroidota")
nrow(morex_prova_others)

barke_prova_proteo <- subset(Barke_prova, phylum == "Proteobacteria")
nrow(barke_prova_proteo)
barke_prova_actino <- subset(Barke_prova, phylum == "Actinobacteriota")
nrow(barke_prova_actino)
barke_prova_bacte <- subset(Barke_prova, phylum == "Bacteroidota")
nrow(barke_prova_bacte)
barke_prova_others <- subset(Barke_prova, phylum != "Proteobacteria" & phylum != "Actinobacteriota" & phylum != "Bacteroidota")
nrow(barke_prova_others)

Morex_proteo <- as.character(morex_prova_proteo$Morex_enriched_G0_005)
Morex_actino <- as.character(morex_prova_actino$Morex_enriched_G0_005)
Morex_bacte <- as.character(morex_prova_bacte$Morex_enriched_G0_005)

Barke_proteo <- as.character(barke_prova_proteo$Barke_enriched_G0_005)
Barke_actino <- as.character(barke_prova_actino$Barke_enriched_G0_005)
Barke_bacte <- as.character(barke_prova_bacte$Barke_enriched_G0_005)


fig_colors <- ifelse(rownames(temat_1A) %in% Morex_enriched_G0_005 | rownames(temat_1A) %in% Barke_enriched_G0_005, "black","darkgrey")
names(fig_colors) <- rownames(temat_1A)
fig_colors[Morex_proteo] <- "#009E73"
fig_colors[Morex_actino] <- "red"
fig_colors[Morex_bacte] <- "#0072B2"

fig_colors[Barke_proteo] <- "#009E73"
fig_colors[Barke_actino] <- "red"
fig_colors[Barke_bacte] <- "#0072B2"


tern_e(temat_1A, scale = 1, prop=T, col=fig_colors, grid_color="black", labels_color="black", pch=19, main="G0 effect phyla v3")
#save this plot for figure generation

Morex_proteo_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_proteo, "Genus"])
Morex_proteo_taxa
Morex_actino_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_actino, "Genus"])
Morex_actino_taxa
Morex_bacte_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_bacte, "Genus"])
Morex_bacte_taxa

Barke_proteo_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_proteo,"Genus"])
Barke_proteo_taxa
Barke_actino_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_actino,"Genus"])
Barke_actino_taxa
Barke_bacte_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_bacte,"Genus"])
Barke_bacte_taxa


#G24
JH12_G24_cds <- readRDS("JH12_G24_genus_cds_0221.rds")
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


Morex_prova_g24 <- as.data.frame(Morex_enriched_G24_005)
Barke_prova_g24 <- as.data.frame(Barke_enriched_G24_005)

Morex_prova_g24$phylum <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_prova_g24$Morex_enriched_G24_005,"Phylum"])
Barke_prova_g24$phylum <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_prova_g24$Barke_enriched_G24_005,"Phylum"])

morex_prova_g24_proteo <- subset(Morex_prova_g24, phylum == "Proteobacteria")
nrow(morex_prova_g24_proteo)
morex_prova_g24_actino <- subset(Morex_prova_g24, phylum == "Actinobacteriota")
nrow(morex_prova_g24_actino)
morex_prova_g24_bacte <- subset(Morex_prova_g24, phylum == "Bacteroidota")
nrow(morex_prova_g24_bacte)
morex_prova_g24_others <- subset(Morex_prova_g24, phylum != "Proteobacteria" & phylum != "Actinobacteriota" & phylum != "Bacteroidota")
nrow(morex_prova_g24_others)

barke_prova_g24_proteo <- subset(Barke_prova_g24, phylum == "Proteobacteria")
nrow(barke_prova_g24_proteo)
barke_prova_g24_actino <- subset(Barke_prova_g24, phylum == "Actinobacteriota")
nrow(barke_prova_g24_actino)
barke_prova_g24_bacte <- subset(Barke_prova_g24, phylum == "Bacteroidota")
nrow(barke_prova_g24_bacte)
barke_prova_g24_others <- subset(Barke_prova_g24, phylum != "Proteobacteria" & phylum != "Actinobacteriota" & phylum != "Bacteroidota")
nrow(barke_prova_g24_others)

Morex_proteo_g24 <- as.character(morex_prova_g24_proteo$Morex_enriched_G24_005)
Morex_actino_g24 <- as.character(morex_prova_g24_actino$Morex_enriched_G24_005)
Morex_bacte_g24 <- as.character(morex_prova_g24_bacte$Morex_enriched_G24_005)

Barke_proteo_g24 <- as.character(barke_prova_g24_proteo$Barke_enriched_G24_005)
Barke_actino_g24 <- as.character(barke_prova_g24_actino$Barke_enriched_G24_005)
Barke_bacte_g24 <- as.character(barke_prova_g24_bacte$Barke_enriched_G24_005)


fig_colors_g24 <- ifelse(rownames(temat_1A) %in% Morex_enriched_G24_005 | rownames(temat_1A) %in% Barke_enriched_G24_005, "black","darkgrey")
names(fig_colors_g24) <- rownames(temat_1A)
fig_colors_g24[Morex_proteo_g24] <- "#009E73"
fig_colors_g24[Morex_actino_g24] <- "red"
fig_colors_g24[Morex_bacte_g24] <- "#0072B2"

fig_colors_g24[Barke_proteo_g24] <- "#009E73"
fig_colors_g24[Barke_actino_g24] <- "red"
fig_colors_g24[Barke_bacte_g24] <- "#0072B2"


tern_e(temat_1A, scale = 1, prop=T, col=fig_colors_g24, grid_color="black", labels_color="black", pch=19, main="G24 effect v3")
#save this plot for figure generation

Morex_proteo_G24_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_proteo_g24, "Genus"])
Morex_proteo_G24_taxa
Morex_actino_G24_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_actino_g24, "Genus"])
Morex_actino_G24_taxa
Morex_bacte_G24_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_bacte_g24, "Genus"])
Morex_bacte_G24_taxa

Barke_proteo_G24_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_proteo_g24,"Genus"])
Barke_proteo_G24_taxa
Barke_actino_G24_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_actino_g24,"Genus"])
Barke_actino_G24_taxa
Barke_bacte_G24_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_bacte_g24,"Genus"])
Barke_bacte_G24_taxa


#G46
JH12_G46_cds <- readRDS("JH12_G46_genus_cds_0221.rds")
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


Morex_prova_g46 <- as.data.frame(Morex_enriched_G46_005)
Barke_prova_g46 <- as.data.frame(Barke_enriched_G46_005)

Morex_prova_g46$phylum <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_prova_g46$Morex_enriched_G46_005,"Phylum"])
Barke_prova_g46$phylum <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_prova_g46$Barke_enriched_G46_005,"Phylum"])

morex_prova_g46_proteo <- subset(Morex_prova_g46, phylum == "Proteobacteria")
nrow(morex_prova_g46_proteo)
morex_prova_g46_actino <- subset(Morex_prova_g46, phylum == "Actinobacteriota")
nrow(morex_prova_g46_actino)
morex_prova_g46_bacte <- subset(Morex_prova_g46, phylum == "Bacteroidota")
nrow(morex_prova_g46_bacte)
morex_prova_g46_others <- subset(Morex_prova_g46, phylum != "Proteobacteria" & phylum != "Actinobacteriota" & phylum != "Bacteroidota")
nrow(morex_prova_g46_others)

barke_prova_g46_proteo <- subset(Barke_prova_g46, phylum == "Proteobacteria")
nrow(barke_prova_g46_proteo)
barke_prova_g46_actino <- subset(Barke_prova_g46, phylum == "Actinobacteriota")
nrow(barke_prova_g46_actino)
barke_prova_g46_bacte <- subset(Barke_prova_g46, phylum == "Bacteroidota")
nrow(barke_prova_g46_bacte)
barke_prova_g46_others <- subset(Barke_prova_g46, phylum != "Proteobacteria" & phylum != "Actinobacteriota" & phylum != "Bacteroidota")
nrow(barke_prova_g46_others)

Morex_proteo_g46 <- as.character(morex_prova_g46_proteo$Morex_enriched_G46_005)
Morex_actino_g46 <- as.character(morex_prova_g46_actino$Morex_enriched_G46_005)
Morex_bacte_g46 <- as.character(morex_prova_g46_bacte$Morex_enriched_G46_005)

Barke_proteo_g46 <- as.character(barke_prova_g46_proteo$Barke_enriched_G46_005)
Barke_actino_g46 <- as.character(barke_prova_g46_actino$Barke_enriched_G46_005)
Barke_bacte_g46 <- as.character(barke_prova_g46_bacte$Barke_enriched_G46_005)


fig_colors_g46 <- ifelse(rownames(temat_1A) %in% Morex_enriched_G46_005 | rownames(temat_1A) %in% Barke_enriched_G46_005, "black","darkgrey")
names(fig_colors_g46) <- rownames(temat_1A)
fig_colors_g46[Morex_proteo_g46] <- "#009E73"
fig_colors_g46[Morex_actino_g46] <- "red"
fig_colors_g46[Morex_bacte_g46] <- "#0072B2"

fig_colors_g46[Barke_proteo_g46] <- "#009E73"
fig_colors_g46[Barke_actino_g46] <- "red"
fig_colors_g46[Barke_bacte_g46] <- "#0072B2"


tern_e(temat_1A, scale = 1, prop=T, col=fig_colors_g46, grid_color="black", labels_color="black", pch=19, main="G46 effect v3")
#save this plot for figure generation


Morex_proteo_G46_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_proteo_g46, "Genus"])
Morex_proteo_G46_taxa
Morex_actino_G46_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_actino_g46, "Genus"])
Morex_actino_G46_taxa
Morex_bacte_G46_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Morex_bacte_g46, "Genus"])
Morex_bacte_G46_taxa

Barke_proteo_G46_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_proteo_g46,"Genus"])
Barke_proteo_G46_taxa
Barke_actino_G46_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_actino_g46,"Genus"])
Barke_actino_G46_taxa
Barke_bacte_G46_taxa <- as.character(tax_table(JH12_data_phyloseq_genus_rare)[Barke_bacte_g46,"Genus"])
Barke_bacte_G46_taxa
