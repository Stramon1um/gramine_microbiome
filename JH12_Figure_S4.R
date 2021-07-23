#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Figure Supplementary S4 used in Maver et al., manuscript
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
library("grid")
library("gridExtra")

colorblind_Palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#import the Phyloseq object 2 file
JH12_alfa_div <- readRDS(file = "JH12_data_phyloseq_ASV_0221_Silva138.rds")
JH12_alfa_div


###
#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_rare <- read.delim("JH12_data_phyloseq_rare_table_counts_140219.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
###


#Index calculations
JH12_alpha_rare <-  estimate_richness(JH12_alfa_div, measures = c("Observed", "Shannon", "Chao1"))
JH12_alpha_rare

#design file
design <- read.delim("Map_JH12_2.txt", sep = "\t", header=TRUE, row.names=1)
design

design_treatment <- as.data.frame(design[, 5])
rownames(design_treatment) <- rownames(design)
colnames(design_treatment) <- c("treatment")
design_treatment

#description
design_description <- as.data.frame(design[, 6])
rownames(design_description) <- rownames(design)
colnames(design_description) <- c("description")
design_description 

#data frame Genotype_Description
design_TD <- cbind(design_treatment, design_description)

#remove nursery sample
design_TD <- design_TD[colnames(dat_count_rare), ] 

#
##
###
#### OBSERVED ASVs

#Observed ASVs
JH12_alpha_rare_Observed <- as.data.frame(JH12_alpha_rare[ ,1])
rownames(JH12_alpha_rare_Observed) <- rownames(JH12_alpha_rare)
colnames(JH12_alpha_rare_Observed) <- c("Observed")

#Combine the dataset sample description and Observed ASVs
JH12_alpha_rare_Observed_TD <- cbind(design_TD, JH12_alpha_rare_Observed)
JH12_alpha_rare_Observed_TD <- as.data.frame(JH12_alpha_rare_Observed_TD)
JH12_alpha_rare_Observed_TD$treatment
JH12_alpha_rare_Observed_TD$description

#Order the levels according to a defined order
JH12_alpha_rare_Observed_TD$treatment <- ordered(JH12_alpha_rare_Observed_TD$treatment, levels=c("G0", "G24", "G46")) 
JH12_alpha_rare_Observed_TD$description <- ordered(JH12_alpha_rare_Observed_TD$description, levels=c("QUARRY", "BARKEQUARRY", "MOREXQUARRY")) 


description_names <- c(
  `QUARRY` = "Bulk soil",
  `BARKEQUARRY` = "Barke",
  `MOREXQUARRY` = "Morex")

p1<-ggplot(JH12_alpha_rare_Observed_TD, aes(x=treatment, y=Observed, fill=treatment))+
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 4.5)+
  #geom_jitter(size = 3, shape = 16)+
  geom_point(size=3, shape = 16)+
  scale_fill_manual(name = "Concentration", labels = c("G0", "G24", "G46"), values = colorblind_Palette)+
  facet_wrap(~ description, labeller = as_labeller(description_names))+
  theme_bw()+
  ylab("Observed ASVs")+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.position="none",
        legend.title = element_text(size = 16))


#
##
###
#### CHAO1

#Chao1 ASVs
JH12_alpha_rare_Chao1 <- as.data.frame(JH12_alpha_rare[ ,2])
rownames(JH12_alpha_rare_Chao1) <- rownames(JH12_alpha_rare)
colnames(JH12_alpha_rare_Chao1) <- c("Chao1")

#Combine the dataset sample description and Chao1 ASVs
JH12_alpha_rare_Chao1_TD <- cbind(design_TD, JH12_alpha_rare_Chao1)
JH12_alpha_rare_Chao1_TD <- as.data.frame(JH12_alpha_rare_Chao1_TD)
JH12_alpha_rare_Chao1_TD$treatment
JH12_alpha_rare_Chao1_TD$description

#Order the levels according to a defined order
JH12_alpha_rare_Chao1_TD$treatment <- ordered(JH12_alpha_rare_Chao1_TD$treatment, levels=c("G0", "G24", "G46")) 
JH12_alpha_rare_Chao1_TD$description <- ordered(JH12_alpha_rare_Chao1_TD$description, levels=c("QUARRY", "BARKEQUARRY", "MOREXQUARRY")) 


p2<-ggplot(JH12_alpha_rare_Chao1_TD, aes(x=treatment, y=Chao1, fill=treatment))+
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 4.5)+
  #geom_jitter(size = 3, shape = 16)+
  geom_point(size=3, shape = 16)+
  scale_fill_manual(name = "Concentration", labels = c("G0", "G24", "G46"), values = colorblind_Palette)+
  facet_wrap(~ description, labeller = as_labeller(description_names))+
  theme_bw()+
  ylab("Chao1")+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.position="none",
        legend.title = element_text(size = 16))


#
##
###
#### SHANNON

#Shannon ASVs
JH12_alpha_rare_Shannon <- as.data.frame(JH12_alpha_rare[ ,4])
rownames(JH12_alpha_rare_Shannon) <- rownames(JH12_alpha_rare)
colnames(JH12_alpha_rare_Shannon) <- c("Shannon")

#Combine the dataset sample description and Shannon ASVs
JH12_alpha_rare_Shannon_TD <- cbind(design_TD, JH12_alpha_rare_Shannon)
JH12_alpha_rare_Shannon_TD <- as.data.frame(JH12_alpha_rare_Shannon_TD)
JH12_alpha_rare_Shannon_TD$treatment
JH12_alpha_rare_Shannon_TD$description

#Order the levels according to a defined order
JH12_alpha_rare_Shannon_TD$treatment <- ordered(JH12_alpha_rare_Shannon_TD$treatment, levels=c("G0", "G24", "G46")) 
JH12_alpha_rare_Shannon_TD$description <- ordered(JH12_alpha_rare_Shannon_TD$description, levels=c("QUARRY", "BARKEQUARRY", "MOREXQUARRY")) 


p3<-ggplot(JH12_alpha_rare_Shannon_TD, aes(x=treatment, y=Shannon, fill=treatment))+
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 4.5)+
  #geom_jitter(size = 3, shape = 16)+
  geom_point(size=3, shape = 16)+
  scale_fill_manual(name = "Concentration", labels = c("G0", "G24", "G46"), values = colorblind_Palette)+
  facet_wrap(~ description, labeller = as_labeller(description_names))+
  theme_bw()+
  ylab("Shannon")+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.position="bottom",
        legend.title = element_text(size = 16))

# Merging the 3 plots

lay5 <- rbind(c(1),
              c(2),
              c(3))

grid.arrange(p1, p2, p3, layout_matrix = lay5,
             top = textGrob("",gp=gpar(fontsize=14,font=2)))