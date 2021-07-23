#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Figure Supplementary S2 used in Maver et al., manuscript
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

library("ggplot2")
library("agricolae")


colorblind_Palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

dry_w<-read.csv('dry_weight.csv',header = TRUE, sep=",")

dry_w_subset = subset(dry_w, genotype!="Bulk")


# Figure S2
ggplot(dry_w_subset, aes(x=treatment, y=dw_g))+
  #geom_bar(aes(fill=treatment),width = 0.5, stat = "count")+
  geom_boxplot(aes(fill=treatment))+
  facet_wrap(~genotype, strip.position = "top")+
  #scale_y_continuous(limits = c(0,600), expand = c(0, 0))+
  #coord_flip() +
  #scale_fill_brewer(palette="Paired")+
  #ggtitle("Genotypes response to Gramine application")+
  xlab("genotype/treatments")+
  ylab("Plant shoot dry weight (g)")+
  scale_fill_manual(name = "Concentration", labels = c("G0", "G24", "G46"), values = colorblind_Palette)+
  theme_bw()+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.position="bottom",
        legend.title = element_text(size = 16))

#ANOVA
dry_w_stats <- aov(dw_g ~ treatment, data = dry_w_subset)
summary(dry_w_stats)

#ANOVA
dry_w_stats <- aov(dw_g ~ treatment, data = subset(dry_w_subset, genotype=="Barke"))
summary(dry_w_stats)

#tukey test
TukeyHSD(dry_w_stats, ordered = TRUE)
tukey.test2 <- HSD.test(dry_w_stats, trt = ('treatment'))
tukey.test2

#ANOVA
dry_w_stats <- aov(dw_g ~ treatment, data = subset(dry_w_subset, genotype=="Morex"))
summary(dry_w_stats)

#tukey test
TukeyHSD(dry_w_stats, ordered = TRUE)
tukey.test2 <- HSD.test(dry_w_stats, trt = ('treatment'))
tukey.test2
