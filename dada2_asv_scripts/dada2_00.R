#!/usr/bin/env Rscript

#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 24
#$ -mods l_hard mfree 200G

# 1st stage of dada2 workflow: QC and trimming

suppressPackageStartupMessages(library("getopt"))
num_threads=24

optspec=matrix(c(
	'reads',    'r', 1, 'character', 'Path to directory of fastq files',
	'metadata', 'm', 1, 'character', 'Path to metadata file',
	'truncF',   'F', 2, 'integer',   'Length to truncate forward reads to (default: no truncation)',
	'truncR',   'R', 2, 'integer',   'Length to truncate forward reads to (default: no truncation)',
	'maxeeF',   'q', 2, 'integer',   'Maximum expected errors for forward reads (default: 2)',
	'maxeeR',   'w', 2, 'integer',   'Maximum expected errors for reverse reads (default: 2)',
	'help',     'h', 0, 'logical',   'Display help'
),byrow=TRUE,ncol=5)

opt=getopt(optspec)

if (!is.null(opt$help)) {
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$reads)) {
	cat("Error: no reads argument provided\n")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$metadata)) {
	cat("Error: no metadata argument provided\n")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$maxeeF)) {
	maxeeF=2
} else {
	maxeeF=opt$maxeeF
}

if (is.null(opt$maxeeR)) {
	maxeeR=2
} else {
	maxeeR=opt$maxeeF
}

suppressPackageStartupMessages({
	library("ShortRead")
	library("ggplot2")
	library("gridExtra")
	library("dada2")
})

writeLines(capture.output(sessionInfo()), "00_trimming.sessionInfo.txt")

set.seed(100)

read_path=opt$reads
filt_path="filtered"
metadata_file=opt$metadata

if(!file.exists('cache')) {
	dir.create('cache')
}

message("dada2 read trimming\n==============")
message(paste("Read directory: ",read_path))
message(paste("Metadata:",metadata_file))

fns<-sort(list.files(read_path, full.names=TRUE))
fnFs<-fns[grepl("_R1",fns)]
fnRs<-fns[grepl("_R2",fns)]

# Determine median read length of forward and reverse reads based on 1st fastq file
reads<-readFastq(fnFs[[1]])
forward_length<-median(width(reads))
reads<-readFastq(fnRs[[1]])
reverse_length<-median(width(reads))
read_lengths=c(forward_length,reverse_length)
saveRDS(read_lengths,file='cache/read_lengths.txt')

if (length(opt$truncF)) {
	truncFlen<-opt$truncF
} else {
	truncFlen<-0
}

if (length(opt$truncR)) {
	truncRlen<-opt$truncR
} else {
	truncRlen<-0
}

message(paste("Truncate forward:",truncFlen))
message(paste("Truncate reverse:",truncRlen))
message(paste("Max expected errors forward:",maxeeF))
message(paste("Max expected errors reverse:",maxeeR))

message("\nRead lengths\n============")
message(paste0("Forward: ",forward_length))
message(paste0("Reverse: ",reverse_length))
samples<-read.csv(metadata_file,header = TRUE,sep="\t")

# Qiime mapping files start with '#', but this screws up the column naming, so reset 1st column name correctly...
colnames(samples)[1]<-'SampleID'
rownames(samples)<-samples$SampleID
saveRDS(samples,file='cache/samdf.rds')
# Assume that sample name is first element of '_' delimited filename...
sample.names<-sapply(strsplit(basename(fnFs),"_"),'[',1)
saveRDS(sample.names,file='cache/sample_names.rds')

filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
names(filtFs)<-sample.names
names(filtRs)<-sample.names

if(!file.exists('trimming_results.rds')){
	message("\nFiltering reads:\n================")
	message(paste0("Forward length: ",truncFlen))
	message(paste0("Reverse length: ",truncRlen))
	if (truncFlen>0) {
		trim_results <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncFlen,truncRlen),
			maxN=0, maxEE=c(maxeeF,maxeeR), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=num_threads)
	} else {
		trim_results <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
			maxN=0, maxEE=c(maxeeF,maxeeR), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=num_threads)
	}
	saveRDS(trim_results,'cache/trimming_results.rds')
} else {
	message('\nRead filtering has already been completed.\n')
	message('Remove trimming_results.rds to rerun\n')
}
