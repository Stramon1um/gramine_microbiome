#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Figure 3 and Figure 4 used in Maver et al., manuscript
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
library ("UpSetR")
library("ggplot2")


#
##
### JH12 UPSetR plots --> Figure 3

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

#G24
JH12_G24_cds <- readRDS(file = "JH12_G24_genus_cds_0221.rds")
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

#G46
JH12_G46_cds <- readRDS(file = "JH12_G46_genus_cds_0221.rds")
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

#Prepare the data for UpSetR visualisation
#Morex
Morex_rhizo_enriched_FDR005_G0_counts <- as.data.frame(Morex_rhizo_G0[Morex_rhizo_enriched_FDR005_G0, 1])
dim(Morex_rhizo_enriched_FDR005_G0_counts)
rownames(Morex_rhizo_enriched_FDR005_G0_counts) <- Morex_rhizo_enriched_FDR005_G0
colnames(Morex_rhizo_enriched_FDR005_G0_counts) <- c("counts_Morex_G0")
Morex_rhizo_enriched_FDR005_G0_counts[Morex_rhizo_enriched_FDR005_G0_counts > 1] <- 1
Morex_rhizo_enriched_FDR005_G0_counts
dim(Morex_rhizo_enriched_FDR005_G0_counts)
Morex_rhizo_enriched_FDR005_G24_counts <- as.data.frame(Morex_rhizo_G24[Morex_rhizo_enriched_FDR005_G24, 1])
dim(Morex_rhizo_enriched_FDR005_G24_counts)
rownames(Morex_rhizo_enriched_FDR005_G24_counts) <- Morex_rhizo_enriched_FDR005_G24
colnames(Morex_rhizo_enriched_FDR005_G24_counts) <- c("counts_Morex_G24")
Morex_rhizo_enriched_FDR005_G24_counts[Morex_rhizo_enriched_FDR005_G24_counts > 1] <- 1
Morex_rhizo_enriched_FDR005_G24_counts
dim(Morex_rhizo_enriched_FDR005_G24_counts)
Morex_rhizo_enriched_FDR005_G46_counts <- as.data.frame(Morex_rhizo_G46[Morex_rhizo_enriched_FDR005_G46, 1])
dim(Morex_rhizo_enriched_FDR005_G46_counts)
rownames(Morex_rhizo_enriched_FDR005_G46_counts) <- Morex_rhizo_enriched_FDR005_G46
colnames(Morex_rhizo_enriched_FDR005_G46_counts) <- c("counts_Morex_G46")
Morex_rhizo_enriched_FDR005_G46_counts[Morex_rhizo_enriched_FDR005_G46_counts > 1] <- 1
Morex_rhizo_enriched_FDR005_G46_counts
dim(Morex_rhizo_enriched_FDR005_G46_counts)
#combine the datasets: note they have unequal values
#define a list of unique OTUs
OTU_list <- unique(c(rownames(Morex_rhizo_enriched_FDR005_G0_counts), rownames(Morex_rhizo_enriched_FDR005_G24_counts)))
OTU_list <- unique(c(rownames(Morex_rhizo_enriched_FDR005_G46_counts), OTU_list))
length(OTU_list)
#G0
G0_eriched_merging <- as.data.frame(Morex_rhizo_enriched_FDR005_G0_counts[OTU_list, ])
colnames(G0_eriched_merging) <- c("Genera_G0")
row.names(G0_eriched_merging) <- as.vector(OTU_list)
#G24
G24_eriched_merging <- as.data.frame(Morex_rhizo_enriched_FDR005_G24_counts[OTU_list, ])
colnames(G24_eriched_merging) <- c("Genera_G24")
row.names(G24_eriched_merging) <- as.vector(OTU_list)
#G46
G46_eriched_merging <- as.data.frame(Morex_rhizo_enriched_FDR005_G46_counts[OTU_list, ])
colnames(G46_eriched_merging) <- c("Genera_G46")
row.names(G46_eriched_merging) <- as.vector(OTU_list)
#Merge the dataset
Morex_OTUs <- cbind(G0_eriched_merging, G24_eriched_merging)
Morex_OTUs <- cbind(Morex_OTUs, G46_eriched_merging)
#set NA to 0
Morex_OTUs[is.na(Morex_OTUs)] <- 0
dim(Morex_OTUs)
#visualisation
plot_UPSET_Morex<-upset(Morex_OTUs, sets = c("Genera_G0", "Genera_G24", "Genera_G46"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on",text.scale = 2, point.size = 4,  line.size = 1, sets.x.label = "Enriched genera", mainbar.y.label = "Genera simultaneously enriched",)

plot_UPSET_Morex  # Figure 3A

#Barke
Barke_rhizo_enriched_FDR005_G0_counts <- as.data.frame(Barke_rhizo_G0[Barke_rhizo_enriched_FDR005_G0, 1])
dim(Barke_rhizo_enriched_FDR005_G0_counts)
rownames(Barke_rhizo_enriched_FDR005_G0_counts) <- Barke_rhizo_enriched_FDR005_G0
colnames(Barke_rhizo_enriched_FDR005_G0_counts) <- c("counts_Barke_G0")
Barke_rhizo_enriched_FDR005_G0_counts[Barke_rhizo_enriched_FDR005_G0_counts > 1] <- 1
Barke_rhizo_enriched_FDR005_G0_counts
dim(Barke_rhizo_enriched_FDR005_G0_counts)
Barke_rhizo_enriched_FDR005_G24_counts <- as.data.frame(Barke_rhizo_G24[Barke_rhizo_enriched_FDR005_G24, 1])
dim(Barke_rhizo_enriched_FDR005_G24_counts)
rownames(Barke_rhizo_enriched_FDR005_G24_counts) <- Barke_rhizo_enriched_FDR005_G24
colnames(Barke_rhizo_enriched_FDR005_G24_counts) <- c("counts_Barke_G24")
Barke_rhizo_enriched_FDR005_G24_counts[Barke_rhizo_enriched_FDR005_G24_counts > 1] <- 1
Barke_rhizo_enriched_FDR005_G24_counts
dim(Barke_rhizo_enriched_FDR005_G24_counts)
Barke_rhizo_enriched_FDR005_G46_counts <- as.data.frame(Barke_rhizo_G46[Barke_rhizo_enriched_FDR005_G46, 1])
dim(Barke_rhizo_enriched_FDR005_G46_counts)
rownames(Barke_rhizo_enriched_FDR005_G46_counts) <- Barke_rhizo_enriched_FDR005_G46
colnames(Barke_rhizo_enriched_FDR005_G46_counts) <- c("counts_Barke_G46")
Barke_rhizo_enriched_FDR005_G46_counts[Barke_rhizo_enriched_FDR005_G46_counts > 1] <- 1
Barke_rhizo_enriched_FDR005_G46_counts
dim(Barke_rhizo_enriched_FDR005_G46_counts)
#combine the datasets: note they have unequal values
#define a list of unique OTUs
OTU_list <- unique(c(rownames(Barke_rhizo_enriched_FDR005_G0_counts), rownames(Barke_rhizo_enriched_FDR005_G24_counts)))
OTU_list <- unique(c(rownames(Barke_rhizo_enriched_FDR005_G46_counts), OTU_list))
length(OTU_list)
#G0
G0_eriched_merging <- as.data.frame(Barke_rhizo_enriched_FDR005_G0_counts[OTU_list, ])
colnames(G0_eriched_merging) <- c("Genera_G0")
row.names(G0_eriched_merging) <- as.vector(OTU_list)
#G24
G24_eriched_merging <- as.data.frame(Barke_rhizo_enriched_FDR005_G24_counts[OTU_list, ])
colnames(G24_eriched_merging) <- c("Genera_G24")
row.names(G24_eriched_merging) <- as.vector(OTU_list)
#G46
G46_eriched_merging <- as.data.frame(Barke_rhizo_enriched_FDR005_G46_counts[OTU_list, ])
colnames(G46_eriched_merging) <- c("Genera_G46")
row.names(G46_eriched_merging) <- as.vector(OTU_list)
#Merge the dataset
Barke_OTUs <- cbind(G0_eriched_merging, G24_eriched_merging)
Barke_OTUs <- cbind(Barke_OTUs, G46_eriched_merging)
#set NA to 0
Barke_OTUs[is.na(Barke_OTUs)] <- 0
#visualisation
plot_UPSET_Barke<-upset(Barke_OTUs, sets = c("Genera_G0", "Genera_G24", "Genera_G46"), sets.bar.color = "#56B4E9",
                        order.by = "freq", empty.intersections = "on",text.scale = 2, point.size = 4,  line.size = 1, sets.x.label = "Enriched genera", mainbar.y.label = "Genera simultaneously enriched")

plot_UPSET_Barke #Figure 3B

#
##
### JH12 HEATMAP --> Figure 4

#Identify the "core OTUs" that are a) enriched in G24 and G46 versus G0 and b) conserved across the two genotypes
#Morex
Morex_G_treatment_enriched <- intersect(rownames(Morex_rhizo_enriched_FDR005_G24_counts), rownames(Morex_rhizo_enriched_FDR005_G46_counts))
#remove the OTUs from the former list that are enriched in G0 (final number should be 209)
Morex_G_treatment_enriched_2 <- setdiff(Morex_G_treatment_enriched, rownames(Morex_rhizo_enriched_FDR005_G0_counts))
length(Morex_G_treatment_enriched_2)
#Barke
Barke_G_treatment_enriched <- intersect(rownames(Barke_rhizo_enriched_FDR005_G24_counts), rownames(Barke_rhizo_enriched_FDR005_G46_counts))
#remove the OTUs from the former list that are enriched in G0 (final number should be 202)
Barke_G_treatment_enriched_2 <- setdiff(Barke_G_treatment_enriched, rownames(Barke_rhizo_enriched_FDR005_G0_counts))
length(Barke_G_treatment_enriched_2)
#Common across the two genotype
Elite_G_treatment_enriched <- intersect(Morex_G_treatment_enriched_2, Barke_G_treatment_enriched_2)
length(Elite_G_treatment_enriched)
#Identify the "core OTUs" that are a) depleted in G24 and G46 versus G0 and b) conserved across the two genotypes
#Morex
Morex_G_treatment_depleted <- setdiff(rownames(Morex_rhizo_enriched_FDR005_G0_counts), rownames(Morex_rhizo_enriched_FDR005_G24_counts)) 
Morex_G_treatment_depleted_2 <- setdiff(rownames(Morex_rhizo_enriched_FDR005_G0_counts), rownames(Morex_rhizo_enriched_FDR005_G46_counts)) 
Morex_G_treatment_depleted_3 <- intersect(Morex_G_treatment_depleted, Morex_G_treatment_depleted_2)
#inspect the dataset, it should be 158
length(Morex_G_treatment_depleted_3)
#Barke
Barke_G_treatment_depleted <- setdiff(rownames(Barke_rhizo_enriched_FDR005_G0_counts), rownames(Barke_rhizo_enriched_FDR005_G24_counts)) 
Barke_G_treatment_depleted_2 <- setdiff(rownames(Barke_rhizo_enriched_FDR005_G0_counts), rownames(Barke_rhizo_enriched_FDR005_G46_counts)) 
Barke_G_treatment_depleted_3 <- intersect(Barke_G_treatment_depleted, Barke_G_treatment_depleted_2)
#inspect the dataset, it should be 286
length(Barke_G_treatment_depleted_3)
#combine the dataset
Elite_G_treatment_depleted <- intersect(Morex_G_treatment_depleted_3, Barke_G_treatment_depleted_3)
length(Elite_G_treatment_depleted)

#Define the genotype effect: estblish a set of pair-wise comparison and intersect with the Elite/Depleted dataset
#Morex
JH12_Morex_cds <- readRDS(file = "JH12_Morex_cds_0221.rds")
JH12_Morex_cds
#execute the differential count analysis with the function DESeq 
JH12_cds_test <- DESeq(JH12_Morex_cds, fitType="local", betaPrior = FALSE) 

#define the OTUs differentially enriched in the rhizosphere samples
Morex_G24 <- results(JH12_cds_test, contrast = c("Treatments", "G0", "G24")) 
Morex_G46 <- results(JH12_cds_test, contrast = c("Treatments", "G0", "G46")) 

#inspect the result files
Morex_G24
Morex_G46

#extract  OTUs whose adjusted p.value in a given comparison is below 0.05 
Morex_G24_FDR005 <- Morex_G24[(rownames(Morex_G24)[which(Morex_G24$padj <0.05)]), ]
Morex_G46_FDR005 <- Morex_G46[(rownames(Morex_G46)[which(Morex_G46$padj <0.05)]), ]
#Inspect the files
Morex_G24_FDR005
Morex_G46_FDR005
#no OTUs differentially enriched in rhizosphere samples at G24 and 1 at G46

#Barke
JH12_Barke_cds <- readRDS(file = "JH12_Barke_cds_0221.rds")
JH12_Barke_cds
#execute the differential count analysis with the function DESeq 
JH12_cds_test <- DESeq(JH12_Barke_cds, fitType="local", betaPrior = FALSE) 

#define the OTUs differentially enriched in the rhizosphere samples
Barke_G24 <- results(JH12_cds_test, contrast = c("Treatments", "G0", "G24")) 
Barke_G46 <- results(JH12_cds_test, contrast = c("Treatments", "G0", "G46")) 

#inspect the result files
Barke_G24
Barke_G46

#extract  OTUs whose adjusted p.value in a given comparison is below 0.05 
Barke_G24_FDR005 <- Barke_G24[(rownames(Barke_G24)[which(Barke_G24$padj <0.05)]), ]
Barke_G46_FDR005 <- Barke_G46[(rownames(Barke_G46)[which(Barke_G46$padj <0.05)]), ]

#inspect the files
Barke_G24_FDR005
Barke_G46_FDR005

#How many of these OTUs are enriched from soil?

#intersect enriched
Barke_G24_FDR005_rhizo_enriched <- intersect(rownames(Barke_G24_FDR005), rownames(Barke_rhizo_enriched_FDR005_G24_counts))
Barke_G46_FDR005_rhizo_enriched <- intersect(rownames(Barke_G46_FDR005), rownames(Barke_rhizo_enriched_FDR005_G46_counts))

#Proportion of reads
#G24
Barke_G24_FDR005_rhizo_enriched_counts <- as.data.frame(Barke_G24[Barke_G24_FDR005_rhizo_enriched, 1])
dim(Barke_G24_FDR005_rhizo_enriched_counts)
rownames(Barke_G24_FDR005_rhizo_enriched_counts) <- Barke_G24_FDR005_rhizo_enriched
Barke_G24_FDR005_rhizo_enriched_counts_total <- as.data.frame(Barke_G24[rownames(Barke_rhizo_enriched_FDR005_G24_counts), 1])
dim(Barke_G24_FDR005_rhizo_enriched_counts_total)
rownames(Barke_G24_FDR005_rhizo_enriched_counts_total) <- rownames(Barke_rhizo_enriched_FDR005_G24_counts)
ratio_Barke_reads_G24 <- colSums(Barke_G24_FDR005_rhizo_enriched_counts)/colSums(Barke_G24_FDR005_rhizo_enriched_counts_total) * 100
ratio_Barke_reads_G24
#G46
Barke_G46_FDR005_rhizo_enriched_counts <- as.data.frame(Barke_G46[Barke_G46_FDR005_rhizo_enriched, 1])
dim(Barke_G46_FDR005_rhizo_enriched_counts)
rownames(Barke_G46_FDR005_rhizo_enriched_counts) <- Barke_G46_FDR005_rhizo_enriched
Barke_G46_FDR005_rhizo_enriched_counts_total <- as.data.frame(Barke_G46[rownames(Barke_rhizo_enriched_FDR005_G46_counts), 1])
dim(Barke_G46_FDR005_rhizo_enriched_counts_total)
rownames(Barke_G46_FDR005_rhizo_enriched_counts_total) <- rownames(Barke_rhizo_enriched_FDR005_G46_counts)
ratio_Barke_reads_G46 <- colSums(Barke_G46_FDR005_rhizo_enriched_counts)/colSums(Barke_G46_FDR005_rhizo_enriched_counts_total) * 100
ratio_Barke_reads_G46


JH12_data_phyloseq_genus_rare <- readRDS(file = "JH12_data_phyloseq_genus_0221_Silva138.rds")

Elite_G_treatment_enriched_heatmap <- prune_taxa(Elite_G_treatment_enriched,JH12_data_phyloseq_genus_rare)

plot_heatmap<-plot_heatmap(Elite_G_treatment_enriched_heatmap, method = "PCoA", taxa.label = "Order", sample.label ="Treatments", low="#66CCFF", high="#000033", na.value = "white")+
  facet_wrap(Treatments~Microhabitat~Description, nrow = 1, scales = "free_x", strip.position = "bottom")+
  theme_bw()+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_blank(),
        strip.placement = "outside",
        strip.background =  element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size=9),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        legend.position="bottom",
        legend.title = element_text(size = 16),
        axis.ticks.x = element_blank())

plot_heatmap[["labels"]]$fill <- "Reads abundance"

plot_heatmap

##
###
#### % cumulative relative abundance genera 


Morex <- subset_samples(JH12_data_phyloseq_genus_rare, Description =="MOREXQUARRY")
Barke <- subset_samples(JH12_data_phyloseq_genus_rare, Description  =="BARKEQUARRY")

Morex_G <- prune_taxa(Elite_G_treatment_enriched,Morex)
Morex_G
Barke_G <- prune_taxa(Elite_G_treatment_enriched,Barke)
Barke_G

Morex_pc = sum(sample_sums(Morex_G))/sum(sample_sums(Morex))*100
Morex_pc
Barke_pc = sum(sample_sums(Barke_G))/sum(sample_sums(Barke))*100
Barke_pc


#### inspect 18 most abundant genera

Elite_G_treatment_enriched_heatmap@tax_table

write.table(Elite_G_treatment_enriched_heatmap@tax_table, "heatmap_table.csv", sep = ",")
