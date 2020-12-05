#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Beta diversity figures used in Maver et al., manuscript
#  Revision 12/20 
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

#import the Phyloseq object 2 file
JH12_beta_div <- readRDS(file = "JH12_data_phyloseq_genus_rare.rds")
JH12_beta_div

design <- read.delim("Map_JH12_2.txt", sep = "\t", header=TRUE, row.names=1)
design <- design[sample_names(JH12_beta_div), ]
design


#
##
###
##### PCoA bray

#PCoA bray distance
JH12_beta_div_bray <- ordinate(JH12_beta_div, "PCoA", "bray")
plot_ordination(JH12_beta_div, JH12_beta_div_bray , color = "Treatments")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH12_beta_div, JH12_beta_div_bray , shape ="Description", color = "Treatments")
p = p + geom_point(size = 4, alpha = 1)
p + scale_colour_discrete(name = "Concentration", labels = c("G0", "G24", "G46"))+
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
##### PCoA weighted unifrac

#PCoA weighted unifrac distance
#info Unifrac: https://en.wikipedia.org/wiki/UniFrac
JH12_beta_div_WU <- ordinate(JH12_beta_div, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH12_beta_div, JH12_beta_div_WU , color = "Treatments")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH12_beta_div, JH12_beta_div_WU , shape ="Description", color = "Treatments")
p = p + geom_point(size = 4, alpha = 1)
p + scale_colour_discrete(name = "Concentration", labels = c("G0", "G24", "G46"))+
  scale_shape_discrete(name = "Sample type", labels = c("Barke", "Morex", "Bulk"))+
  ggtitle("PCoA 16S data, weighted Unifrac distance")+
  theme_bw()+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_text(color="Black", size=16),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position="right",
        legend.title = element_text(size = 16))


#WU distance adonis
WU <- phyloseq::distance(JH12_beta_div, "unifrac", weighted= TRUE)
#Microhabita effect
adonis(WU ~ Microhabitat * Treatments, data= design, permutations = 5000)
#Description effect
adonis(WU ~ Treatments * Description, data= design, permutations = 5000)


#
##
###
#### CAP Bray-Curtis

#constrained ordination
JH12_CAP <- ordinate(JH12_beta_div, "CAP", "bray", ~ Treatments * Description)
plot_ordination(JH12_beta_div, JH12_CAP, color = "Treatments")

#png("CAP_16S_Bray.png", units="px", width=4960, height=3508, res=300) ##A4 paper size

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH12_beta_div, JH12_CAP, shape ="Description", color = "Treatments")
p = p + geom_point(size = 5, alpha = 1)
p + scale_colour_discrete(name = "Concentration", labels = c("G0", "G24", "G46"))+
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

#dev.off()

anova(JH12_CAP, permutations = how(nperm=5000))

#anova(JH12_CAP ~ Treatments * Description, permutations = how(nperm=5000))
