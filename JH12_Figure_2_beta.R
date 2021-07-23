#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Beta diversity Fig. 2 and Fig. S3 used in Maver et al., manuscript
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
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library("agricolae")

colorblind_Palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#import the Phyloseq object 2 file
JH12_beta_div <- readRDS(file = "JH12_data_phyloseq_genus_0221_Silva138.rds")
JH12_beta_div

design <- read.delim("Map_JH12_2.txt", sep = "\t", header=TRUE, row.names=1)
design <- design[sample_names(JH12_beta_div), ]
design


#
##
###
##### PCoA bray --> Figure Suppplementary S3

#PCoA bray distance
JH12_beta_div_bray <- ordinate(JH12_beta_div, "PCoA", "bray")
plot_ordination(JH12_beta_div, JH12_beta_div_bray , color = "Treatments")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH12_beta_div, JH12_beta_div_bray , shape ="Description", color = "Treatments")
p = p + geom_point(size = 4, alpha = 1)
p + scale_colour_manual(name = "Concentration", labels = c("G0", "G24", "G46"), values = colorblind_Palette)+
  scale_shape_discrete(name = "Sample type", labels = c("Barke", "Morex", "Bulk"))+
  ggtitle("PCoA 16S data, Bray distance")+
  theme_bw()+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_text(color="Black", size=16),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position="right",
        legend.title = element_text(size = 16))



#BC distance adonis
BC <- phyloseq::distance(JH12_beta_div, "bray")
#Microhabitat effect
adonis(BC ~ Microhabitat * Treatments, data= design, permutations = 5000)
#Description effect
adonis(BC ~ Description * Treatments, data= design, permutations = 5000)



#
##
###
#### CAP Bray-Curtis --> Figure 2

#constrained ordination
JH12_CAP <- ordinate(JH12_beta_div, "CAP", "bray", ~ Treatments * Description)
plot_ordination(JH12_beta_div, JH12_CAP, color = "Treatments")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH12_beta_div, JH12_CAP, shape ="Description", color = "Treatments")
p = p + geom_point(size = 5, alpha = 1)
p + scale_colour_manual(name = "Concentration", labels = c("G0", "G24", "G46"), values = colorblind_Palette)+
  scale_shape_discrete(name = "Sample type", labels = c("Barke", "Morex", "Bulk"))+
  theme_bw()+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_text(color="Black", size=16),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position="right",
        legend.title = element_text(size = 16))


anova(JH12_CAP, permutations = how(nperm=5000))