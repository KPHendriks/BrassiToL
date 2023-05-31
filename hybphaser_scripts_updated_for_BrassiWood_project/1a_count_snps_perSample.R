##################
### SNPs count ###
##################

### THIS IS AN UPDATED VERSION OF THE SCRIPT TO BE HANDLED PER SAMPLE BY GNU PARALLEL IN AN HPC ENVIRONMENT
### Kasper Hendriks, 13 January 2022
### Please note that this script was written to be used with the input from GNU Parallel in a bash script as follows:
### parallel -j 67 Rscript ../hybphaser_scripts_updated_for_BrassiWood_project/1a_count_snps_perSample.R config_inclusive.txt :::: ../namelist.txt


####Added by Kasper Hendriks, January 2022.
#When running on a by sample mode on a HPC using GNU Parallel, load sample name from namelist.txt file as follows.
args = commandArgs(trailingOnly=TRUE)
config_file<-args[1]
sample<-args[2]

# load config
#if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)
library(seqinr)
library(stringr)

# generate folders

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)
dir.create(output_Robjects, showWarnings = F, recursive = T)

####KasperH: additionally add folder to same per sample results first
output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset,"tab_sample_temp")
dir.create(output_Robjects, showWarnings = F, recursive = T)


#####################################################################################
### Counting polymorphic sites (masked as ambiguity codes in consensus sequences) ###
#####################################################################################

# getting targets and sample names
if(targets_file_format == "AA"){
  targets <- read.fasta(fasta_file_with_targets, seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
} else if(targets_file_format == "DNA"){
  targets <- read.fasta(fasta_file_with_targets, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
} else {
  print("Warning! Target file type not set properly. Should be 'DNA' or 'AA'!")
}
targets_name <- unique(gsub(".*-","",labels(targets)))

samples <- readLines(path_to_namelist)

if(intronerated_contig=="yes"){
  intronerated_name <- "intronerated" 
  intronerated_underscore <- "_"
} else {
  intronerated_name <- ""
  intronerated_underscore <- ""
}

# function for counting ambiguities
seq_stats <- function(file){
  fasta <- read.fasta(file, as.string=TRUE, set.attributes = FALSE )
  seq <- gsub("N|[?]|-","",fasta[[1]])
  c(round(str_length(seq),0),(str_count(seq,"Y|K|R|S|M|y|k|r|s|m|w") + str_count(seq,"W|D|H|B|V|w|d|h|b|v")*2) )
}


# fill tables with information on snps and sequence length for each sample and locus 
start_time <- Sys.time()

print(paste(sample," ", round(difftime(Sys.time(),start_time, units="secs"),0),"s", sep="") )

tab_sample <- data.frame(targets=targets_name, seq_length=NA, ambis=NA, ambi_prop=NA)

consensus_files <- list.files(file.path(path_to_output_folder,"01_data",sample, paste(intronerated_name,intronerated_underscore,"consensus", sep="")),full.names = T)

for(consensus_file in consensus_files) {
  
  gene <- gsub("(_intronerated|).fasta","",gsub(".*/","",consensus_file))
  
  file.path(path_to_output_folder,"01_data",sample, "consensus", paste(gene,intronerated_underscore,intronerated_name,".fasta",sep=""))
  
  if(file.info(consensus_file)$size !=0){
    stats <- seq_stats(consensus_file)
    tab_sample$seq_length[grep(gene,tab_sample$targets)] <- stats[1]
    tab_sample$ambis[grep(gene,tab_sample$targets)] <- stats[2]  
  } else {
    stats <- c(NA,NA)
    print(paste("WARNING: consensus file ",consensus_file," not found."))
    next
  }  

}

tab_sample$ambi_prop <- tab_sample$ambis/tab_sample$seq_length  

### Generate output table for this sample as an Robject

saveRDS(tab_sample,file=file.path(output_Robjects,paste0(sample,".Rds")), version = 2)
