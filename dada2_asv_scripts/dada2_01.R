#!/usr/bin/env Rscript

#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 48
#$ -jc 4week
#$ -mods l_hard mfree 64G
#$ -adds l_hard h_vmem 64G
#$ -t 1-2

# Script to generate phyloseq objects for 16S amplicon sequencing using dada2,

SGE_TASK_ID<-Sys.getenv('SGE_TASK_ID')
if (SGE_TASK_ID == 1){
    dir='F'
    long_dir='forward'
} else {
    dir='R'
    long_dir='reverse'
}

library("getopt")

num_threads=48

optspec=matrix(c(
	'metadata', 'm', 1, 'character', 'Path to metadata file',
	'help',     'h', 0, 'logical',   'Display help'
),byrow=TRUE,ncol=5)

opt=getopt(optspec)

if (!is.null(opt$help)) {
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$metadata)) {
	cat("Error: no metadata argument provided\n")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

library("ShortRead")
library("ggplot2")
library("gridExtra")
library("dada2")

writeLines(capture.output(sessionInfo()), paste("01_dada2.",dir,".sessionInfo.txt",sep=''))

set.seed(100)

filt_path="filtered"
metadata_file=opt$metadata

if(!file.exists('cache')) {
	cat('Error: cache directory does not exist. Have you run 00_dada2_QC.R?')
    q(status=1)
}

message("dada2 analysis\n==============")
message(paste("Metadata:",metadata_file))
message(paste("Direction:",long_dir,"\n"))

fns<-sort(list.files(filt_path, full.names=TRUE))
if (dir=='F') {
    filts<-fns[grepl("_R1",fns)]
} else {
    filts<-fns[grepl("_R2",fns)]
}

samples<-read.csv(metadata_file,header = TRUE,sep="\t")
# Qiime mapping files start with '#', but this screws up the column naming, so reset 1st column name correctly...
colnames(samples)[1]<-'SampleID'
rownames(samples)<-samples$SampleID

# Assume that sample name is first element of '_' delimited filename...
sample.names<-sapply(strsplit(basename(fns),"_"),'[',1)
sample.names<-unique(sample.names)

# Determine median read length of forward and reverse reads based on 1st fastq file
reads<-readFastq(filts[[1]])
read_length<-median(width(reads))
names(filts)<-sample.names

if(!file.exists(paste("cache/err_",dir,".rds"))) {
	message('\nLearning errors...\n===============')
	err<-learnErrors(filts,multithread=num_threads)
	errplot<-plotErrors(err,nominalQ=TRUE)
	saveRDS(err,paste('cache/err_',dir,'.rds',sep=""))
} else {
	message('\nReloading cached error profile')
	err<-readRDS(paste('cache/err_',dir,'.rds'))
}

if(!file.exists(paste('cache/derep_',dir,'.rds'))) {
	message('\nDereplicating reads\n==================')
	dereps<-derepFastq(filts)
	names(dereps) <- sample.names
	saveRDS(dereps,paste('cache/derep_',dir,'.rds',sep=""))
} else {
	message('\nReloading cached deplicated reads')
	dereps<-readRDS(paste('cache/derep_',dir,'.rds',sep=""))
}

if(!file.exists(paste('cache/dada_',dir,'.rds',sep=""))) {
	message('\nInferring sample composition\n===========================')
	# Use pooling which is somewhat slower, but increases ability to identify
	# rare ASVs
	dada_out <- dada(dereps, err=err, pool=TRUE, selfConsist=TRUE, MAX_CONSIST=20, multithread=num_threads, verbose=TRUE)
	saveRDS(dada_out,paste('cache/dada_',dir,'.rds',sep=""))
} 
message("\nConvergence on error model...")
dada2:::checkConvergence(dada_out[[1]])
