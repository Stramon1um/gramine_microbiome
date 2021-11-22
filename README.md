# Gramine and Microbiome

This repository contains scripts to reproduce analysis and figures presented in the manuscript [Maver et al. 2021](https://peerj.com/articles/12498/), describing the effects of an external addition of the alkaloid gramine on the rhizosphere microbiota recruitment in modern barley plants. 

List of files:

- ```gramine_microbiome.Rproj``` -  Rproject file of this repository

- ```dada2_asv_scripts``` -  folder containing files, scripts and html export of notebook to reproduce the ASVs picking strategy described in the manuscript

- ```sessioninfo.text``` -  packages used for calculations

- ```Map_JH12_2.txt``` - Mapping file of the experiment

- ```JH12_object_construction_v0621_ASVs.R``` and ```JH12_object_construction_v0221_genus.R``` - scripts for quality filtering and phyloseq object construction

Main figures:

- ```JH12_Figure_1_alpha.R``` - script for Figure 1 alpha-diversity calculations 

- ```JH12_Figure_2_beta.R``` - script for Figure 2 beta-diversity calculations + Supp. Figure S3 beta-diversity

- ```JH12_Figure_3_4.R``` - script for Figures 3 DeSeq, UpSetR plot and Figure 4 Heatmap

- ```JH12_Figure_5.R``` - script for Figure 5 Ternary plots

Supplementary figures:

- ```JH12_Figure_S1.R``` - script for Supp. Figure S1 Adsorption models

- ```JH12_Figure_S2.R``` - script for Supp. Figure S2 Dry Weight

- ```JH12_Figure_S4.R``` - script for Supp. Figure S4 alpha-diversity ASVs calculations 

- ```JH12_Figure_S5.R``` - script for Supp. Figure S5 beta-diversity calculations ASVs calculations 