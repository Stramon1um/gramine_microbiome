#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Phyloseq objects and DESeq objects used in Maver et al., manuscript
#  Revision 08/20 
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
library ("ape")

#retrieve R and package versions and compare to the uploaded file in gitHub for the reproducibility of the code
sessionInfo()


#############################################################
#import the count matrix and the desing file
#############################################################


#OTU table generated using QIIME 1.9.0. 
dat_info <- read.delim("JH12_otu_table_SILVA132_97_nc2_noCMC.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the file 
dim(dat_info)
colnames(dat_info)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the MM identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_info[, 1:45]))
OTU_97_reads

#the following information will be reported as general descriptor of the sequencing effort
max(OTU_97_reads)
min(OTU_97_reads)
median(OTU_97_reads)

OTU_97_reads_sum <- sum(OTU_97_reads)
OTU_97_reads_sum 

#total reads
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:45]))
OTU_97_reads_sum

#generate a OTU table with count information only
dat_count_noPlants <- as.data.frame(dat_info[, 1:45])

#design file
design <- read.delim("Map_JH12_2.txt", sep = "\t", header=TRUE, row.names=1)
design

#Create a new taxa table
dat_tax_noPlants <- as.data.frame(dat_info[rownames(dat_info), 46])
rownames(dat_tax_noPlants) <- rownames(dat_info)
colnames(dat_tax_noPlants) <- c("ConsensusLineage")
colnames(dat_tax_noPlants)
#save the above file and in excel we will create a new tax table where each column represents a taxonomic rank
#write.table(dat_tax_noPlants, file="JH12_dat_tax_noPlants.txt", sep="\t")

#############################################################
#Genererate the phyloseq object
#Data required: dat_count_noplants; design, JH07_noPlants_ordered.txt, and 97_otus.tree.gz
#############################################################

#The OTU Table counts
JH12_OTU <- otu_table(dat_count_noPlants, taxa_are_rows=TRUE)

#The taxonomy information
#Note that the file JH16_dat_tax_noPlants_ordered.txt has been generated from the output of lines 69-75  
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
JH12_taxa_ordered <- read.delim ("JH12_dat_tax_noPlants_ordered.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
JH12_taxa <- tax_table(as.matrix(JH12_taxa_ordered))
dim(JH12_taxa)

#The mapping file 
JH12_map <- sample_data(design)

#The phylogenetic tree: from SILVA132
JH12_tree <- read_tree("97_otus.tre")

#merge the files and create the phyloseq object
JH12_data_phyloseq <- merge_phyloseq(JH12_OTU, JH12_taxa, JH12_map,  JH12_tree)

#inspect the generated data
JH12_data_phyloseq
sum(colSums(otu_table(JH12_data_phyloseq)))
dim(dat_count_noPlants)
sum(colSums(dat_count_noPlants))

#abundance filtering: remove OTUs tallyig less than 10 reads in at least 11% of the samples
JH12_data_phyloseq_2 = filter_taxa(JH12_data_phyloseq, function(x) sum(x > 10) > (0.11*length(x)), TRUE)
JH12_data_phyloseq_2
sort(sample_sums(JH12_data_phyloseq_2))

##ratio filtered reads/total reads
ratio <- sum(sample_sums(JH12_data_phyloseq_2))/sum(sample_sums(JH12_data_phyloseq))*100
ratio

#aggregate samples at genus level
JH12_data_phyloseq_genus <- tax_glom(JH12_data_phyloseq_2, taxrank= "Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
#compare the two objects
#OTUs
JH12_data_phyloseq_2
sort(sample_sums(JH12_data_phyloseq_2))
hist(sample_sums(JH12_data_phyloseq_2))
#Genera
JH12_data_phyloseq_genus
sort(sample_sums(JH12_data_phyloseq_genus))
hist(sample_sums(JH12_data_phyloseq_genus))

#rarefy at even sequencing depth (25482) as differences between MM12 (least sequences sample) and MM23 (most sequenced sample) is >10X
JH12_data_phyloseq_genus_rare <- rarefy_even_depth(JH12_data_phyloseq_genus, rngseed=TRUE)
#JH12_data_phyloseq_genus_rare
#sort(sample_sums(JH12_data_phyloseq_genus_rare))
#hist(sample_sums(JH12_data_phyloseq_genus_rare))

#For reproducibility of the code save the above object for phyloseq calculation
#saveRDS(JH12_data_phyloseq_genus_rare, file = "JH12_data_phyloseq_genus_rare.rds")


JH12_data_phyloseq_genus_rare <- readRDS(file = "JH12_data_phyloseq_genus_rare.rds")


#Deseq objects
#subset for Treatment levels
JH12_data_phyloseq_G0 <- subset_samples(JH12_data_phyloseq_genus_rare, Treatments =="G0")
JH12_data_phyloseq_G24 <- subset_samples(JH12_data_phyloseq_genus_rare, Treatments =="G24")
JH12_data_phyloseq_G46 <- subset_samples(JH12_data_phyloseq_genus_rare, Treatments =="G46")

#G0
#extract count data 
JH12_counts_integer <- otu_table(JH12_data_phyloseq_G0)
countData = as.data.frame(JH12_counts_integer)
colnames(JH12_counts_integer)

#the design file containing sample information
colData = design[colnames(JH12_counts_integer), ]
rownames(colData)

#construct a DESeq dataset combining count data and sample information
JH12_G0_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)
#saveRDS(JH12_G0_cds, file = "JH12_G0_genus_cds_2.rds")

#G24
#extract count data 
JH12_counts_integer <- otu_table(JH12_data_phyloseq_G24)
countData = as.data.frame(JH12_counts_integer)
colnames(JH12_counts_integer)

#the design file containing sample information
colData = design[colnames(JH12_counts_integer), ]
rownames(colData)

#construct a DESeq dataset combining count data and sample information
JH12_G24_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)
#saveRDS(JH12_G24_cds, file = "JH12_G24_genus_cds_2.rds")

#G46
#extract count data 
JH12_counts_integer <- otu_table(JH12_data_phyloseq_G46)
countData = as.data.frame(JH12_counts_integer)
colnames(JH12_counts_integer)

#the design file containing sample information
colData = design[colnames(JH12_counts_integer), ]
rownames(colData)

#construct a DESeq dataset combining count data and sample information
JH12_G46_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)
#saveRDS(JH12_G46_cds, file = "JH12_G46_genus_cds_2.rds")

#end



#subset for Genotype
JH12_data_phyloseq_Morex <- subset_samples(JH12_data_phyloseq_genus_rare, Description =="MOREXQUARRY")
JH12_data_phyloseq_Barke <- subset_samples(JH12_data_phyloseq_genus_rare, Description  =="BARKEQUARRY")

#Morex
#extract count data 
JH12_counts_integer <- otu_table(JH12_data_phyloseq_Morex)
countData = as.data.frame(JH12_counts_integer)
colnames(JH12_counts_integer)

#the design file containing sample information
colData = design[colnames(JH12_counts_integer), ]
rownames(colData)

#construct a DESeq dataset combining count data and sample information
JH12_Morex_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Treatments)
#saveRDS(JH12_Morex_cds, file = "JH12_Morex_cds_0520.rds")

#Barke
#extract count data 
JH12_counts_integer <- otu_table(JH12_data_phyloseq_Barke)
countData = as.data.frame(JH12_counts_integer)
colnames(JH12_counts_integer)

#the design file containing sample information
colData = design[colnames(JH12_counts_integer), ]
rownames(colData)

#construct a DESeq dataset combining count data and sample information
JH12_Barke_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Treatments)
#saveRDS(JH12_Barke_cds, file = "JH12_Barke_cds_0520.rds")
