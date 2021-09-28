#!/usr/bin/env Rscript

#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 24
#$ -mods l_hard mfree 64G

# Final stage of dada2 pipeline: Merging paired reads, 
# chimera removal and taxonomic assignment

num_threads=24

check_len<-function(x){
	if(nchar(x)<=50) {
		cat(paste('Short sequence removed: ',x,"\n"))
		return(FALSE)
	} else {
		return(TRUE)
	}
}

getN <- function(x) sum(getUniques(x))

summarise_read_counts<-function(read_counts) {

	# creates summary output of read counts at different points
	# of pipeline, and also plot...

	read_counts<-read_counts %>% 
	  merge(dadaF_counts,on='sample',all.x=TRUE) %>%
	  merge(dadaR_counts,on='sample',all.x=TRUE) %>%
	  merge(m,on='sample',all.x=TRUE) %>%
	  merge(seq_counts,on='sample',all.x=TRUE)
	
	read_counts[is.na(read_counts)] <- 0
	rownames(read_counts)<-read_counts$sample
	 read_counts<-read_counts %>% select (c(-sample))
	
	colnames(read_counts)<-c('Input', 'Filtered', 'denoisedF', 'denoisedR','Merged','Non-chimeric')
	write.table(read_counts,file='read_counts.txt',sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

	read_counts$sample<-row.names(read_counts)
	sample_count=length(read_counts$sample)

	read_counts<-read_counts %>% 
		separate(sample,into=c('prefix','sample_id'),sep="\\.",convert = TRUE,remove=FALSE) %>% 
		arrange(sample_id) %>% 
		select (-c(sample_id, prefix))
	read_counts$sample=factor(read_counts$sample,levels=read_counts$sample)

	read_counts <- read_counts %>% 
		pivot_longer( -sample, names_to='variable', values_to='value',
					  values_transform = list(value = as.integer)) %>%
					  mutate(value=value/1000) 

	read_counts$variable=factor(read_counts$variable,
		levels=c('Input','Filtered','denoisedF','denoisedR','Merged','Non-chimeric'))
	

	plot<-ggplot(read_counts,aes(x=variable,y=value,group=sample))+
	  geom_line()+
	  facet_wrap(~sample,ncol=ceiling(sqrt(sample_count)))+
	  theme_bw()+
	  ylab('Reads x1000')+
	  theme(axis.title.x=element_blank(),
      axis.text.x=element_text(angle=90,hjust=1))

	ggsave(plot,filename='read_counts.pdf',device='pdf',units='cm',height=30,width=30)
}

plot_errs<-function() {

	# Creates plot of observed/predicted errors using dada2::plotErrors()
	cat('Plotting errors\n')
	errF<-readRDS('cache/err_F.rds')
	errR<-readRDS('cache/err_R.rds')

	F_errplot<-plotErrors(errF,nominalQ=TRUE)+ggtitle('Forward Reads')
	R_errplot<-plotErrors(errR,nominalQ=TRUE)+ggtitle('Reverse Reads')

	error_plots<-grid.arrange(F_errplot,R_errplot,ncol=2)
	ggsave('plots/error_rate_plots.png',plot=error_plots,width=40,height=20,units='cm')

}

suppressPackageStartupMessages(library("getopt"))

optspec=matrix(c(
	'name',				'n', 1, 'character',	'Name of job',
	'taxonomy', 		't', 1, 'character',	paste0('Path to taxonomy database to use for classification'),
	'merge_overlap',	'm', 2, 'character',	'Minimum overlap for merging reads (default: 12)',
	'help',     		'h', 0, 'logical',		'Display help'
),byrow=TRUE,ncol=5)

opt=getopt(optspec)
if (!is.null(opt$help)) {
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$name)) {
	cat("Error: no name argument provided\n")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$taxonomy)) {
	cat("Error: no taxonomy database provided")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$merge_overlap)) {opt$merge_overlap = '12'}

suppressPackageStartupMessages({
	library("tidyverse")
	library("ShortRead")
	library("ggplot2")
	library("gridExtra")
	library("dada2")
	library("phyloseq")
	library("dplyr")
})

writeLines(capture.output(sessionInfo()), "02_taxassign.sessionInfo.txt")

set.seed(100)

read_path=opt$reads
filt_path="filtered"
ref_fasta=opt$taxonomy

if(!file.exists('cache/samdf.rds')) {
	cat('samdf.rds does not exist. Please run dada2_00.R to generate this file')
	stop()
} else {
	samdf<-readRDS('cache/samdf.rds')
}

if(!file.exists('cache/read_lengths.txt')){
	cat('read_lengths.txt does not exist. Please run dada2_00.R to generate this file')
	stop()
} else {
	read_lengths<-readRDS('cache/read_lengths.txt')
}

if(!file.exists('cache/sample_names.rds')){
	cat('sample_names.rds does not exist. Please run dada2_00.R to generate this file')
	stop()
} else {
	sample.names<-readRDS('cache/sample_names.rds')
}

if(!file.exists('cache/trimming_results.rds')){
	cat('trimming_results.rds does not exist. Please run dada2_00.R to generate this file')
	stop()
} else {
	read_counts<-readRDS('cache/trimming_results.rds')
	sample<-sapply(strsplit(rownames(read_counts),"_"),'[',1)
	read_counts<-data.frame(cbind(read_counts,sample))
}

## minboot defined bootstrap threshold - recommended for 50 if reads <250 bp, otherwise 80
if(read_lengths[1]<250) {
	min_boot=50
} else {
	min_boot=80
}

cat("\ndada2 analysis\n==============\n")
cat(paste("Taxonomy:",opt$taxonomy,"\n"))
cat(paste("Merge overlap:",opt$merge_overlap,"\n"))
cat(paste("Minimum bootstrap confidence:", min_boot,"\n\n"))

if(!file.exists('cache/derep_F.rds')) {
	cat('derep_F.rds does not exist. Please run dada2_01.R to generate this file')
	stop()
} else {
	cat('Reloading cached forward deplicated reads\n')
	derepFs<-readRDS('cache/derep_F.rds')
}

if(!file.exists('cache/derep_R.rds')) {
	cat('derep_R.rds does not exist. Please run dada2_01.R to generate this file')
	stop()
} else {
	cat('Reloading cached reverse deplicated reads\n')
	derepRs<-readRDS('cache/derep_R.rds')
}

if(!file.exists('cache/dada_F.rds')) {
	cat('dada_F.rds does not exist. Please run dada2_01.R to generate this file')
	stop()
} else {
	cat('Reloading cached dada2 forward object\n')
	dadaF<-readRDS('cache/dada_F.rds')
	dadaF_counts<-data.frame(sapply(dadaF,getN))
	dadaF_counts$sample<-rownames(dadaF_counts)
}

if(!file.exists('cache/dada_R.rds')) {	
	cat('dada_R.rds does not exist. Please run dada2_01.R to generate this file')
	stop()
} else {
	cat('Reloading cached dada2 reverse object\n')
	dadaR<-readRDS('cache/dada_R.rds')
	dadaR_counts<-data.frame(sapply(dadaR,getN))
	dadaR_counts$sample<-rownames(dadaR_counts)
}

if(!file.exists('cache/seqtab.rds')) {
	message('\nMerging reads\n=============')
	mergers<-mergePairs(dadaF, derepFs, dadaR, derepRs, minOverlap=opt$merge_overlap, trimOverhang=TRUE, verbose=TRUE)
	seqtab.all<-makeSequenceTable(mergers)
	m<-data.frame(sapply(mergers,getN))
	m$sample<-rownames(m)
	saveRDS(mergers,file='cache/mergers.rds')

	message("\nRemoving chimeras\n=================")
	seqtab<-removeBimeraDenovo(seqtab.all,multithread=num_threads)
	removed<-round(100-(sum(seqtab)/sum(seqtab.all)*100),digits=2)
	cat(paste('Proportion of reads removed: ',as.character(removed),"%\n",sep=""))

	message("\nChecking for short sequences\n============================")
	short_seqs<-sapply(colnames(seqtab),check_len)
	seqtab<-seqtab[,short_seqs]
	
	seq_counts<-data.frame(rowSums(seqtab))
	seq_counts$sample<-rownames(seq_counts)
	
	saveRDS(seqtab,'cache/seqtab.rds')
	
} else {
	message('Reloading cached seqtab')
	seqtab<-readRDS('cache/seqtab.rds')
	mergers<-readRDS('cache/mergers.rds')
	m<-data.frame(sapply(mergers,getN))
	m$sample<-rownames(m)
	seq_counts<-data.frame(rowSums(seqtab))
	seq_counts$sample<-rownames(seq_counts)
}

plot_errs()
summarise_read_counts(read_counts)

if(!file.exists('cache/taxtab.rds')) {
	message('\nAssigning taxonomy\n===================')
	taxtab<-assignTaxonomy(seqtab, refFasta = ref_fasta, minBoot=min_boot, multithread=num_threads,verbose=TRUE)
	if (length(colnames(taxtab))==6) {
		colnames(taxtab)<-c("Kingdom","Phylum","Class","Order","Family","Genus")
	} else if (length(colnames(taxtab))==7) {
		colnames(taxtab)<-c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
	} else {
		message("Error: Unrecognised taxtab columns: ", colnames(taxtab))
		quit()
	}
	saveRDS(taxtab,file='cache/taxtab.rds')
} else {
	message('Reloading cached taxtab')
	taxtab<-readRDS('cache/taxtab.rds')
}

if(!file.exists(paste(opt$name,'_dada2.rds',sep=''))) {
	message('\nCreating phyloseq object\n========================')
	ps<-phyloseq(
		tax_table(taxtab),
		sample_data(samdf),
		otu_table(t(seqtab),taxa_are_rows=TRUE)
	)
	saveRDS(ps,file=paste(opt$name,'_dada2.rds',sep=''))
}

