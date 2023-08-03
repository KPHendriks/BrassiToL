#################
### SNPs count ###
##################

### THIS IS AN UPDATED VERSION OF THE SCRIPT TO BE HANDLED PER SAMPLE BY GNU PARALLEL IN AN HPC ENVIRONMENT
### Kasper Hendriks, 13 January 2022

####Added by Kasper Hendriks, January 2022.
#When running on a by sample mode on a HPC using GNU Parallel, load sample name from namelist.txt file as follows.
args = commandArgs(trailingOnly=TRUE)
config_file<-args[1]

# load config
#if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)
library(seqinr)
library(stringr)

# define the path to the output folder
output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)

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

# generate empty tables for SNPS and sequence length
tab_snps <- data.frame(loci=targets_name)
tab_length <- data.frame(loci=targets_name)

#Now loop through the list of samples and collect results into main dataframes

for (sample in samples){
  
  # collect the data
  if(file.exists(file.path(output_Robjects,paste0("tab_sample_temp/",sample,".Rds")))){
    tab_sample <- readRDS(file=file.path(output_Robjects,paste0("tab_sample_temp/",sample,".Rds")))
  } else if (file.exists(file.path(path_to_output_folder,paste0("tab_sample_temp/",sample,".Rds")))) {
    tab_sample <- readRDS(file=file.path(path_to_output_folder,paste0("tab_sample_temp/",sample,".Rds")))
  } else {
    print(paste("WARNING: Rds data from script 1a for sample ",sample," not found."))
    next
  }
  
  # bind the data to the general dataframes
  tab_snps <- cbind(tab_snps,tab_sample$ambi_prop)
  tab_length <- cbind(tab_length, tab_sample$seq_length)
  
  # update column name to correspond to this sample
  colnames(tab_snps)[ncol(tab_snps)] <- sample
  colnames(tab_length)[ncol(tab_snps)] <- sample
  
  # rm the tab_sample object
  rm(tab_sample)
  
}

# set the row names for the dataframes
rownames(tab_snps) <- targets_name
rownames(tab_length) <- targets_name

# remove first columns from dataframes for proper use in script 1b
tab_snps <- tab_snps[,-1]
tab_length <- tab_length[,-1]

### Generate output tables and save data in Robjects

saveRDS(tab_snps,file=file.path(output_Robjects,"Table_SNPs.Rds"), version = 2)
saveRDS(tab_length,file=file.path(output_Robjects,"Table_consensus_length.Rds"), version = 2)

